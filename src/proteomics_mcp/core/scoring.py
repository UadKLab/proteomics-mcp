from __future__ import annotations

import dataclasses
import math

import numpy as np

from proteomics_mcp.config import DEFAULT_METRIC_PROFILE, DEFAULT_MIN_PREDICTED_INTENSITY, METRIC_PROFILES
from proteomics_mcp.core.acquisition import acquisition_metadata_summary, precursor_target_mz
from proteomics_mcp.core.annotation import aa_sequence_from_proforma, annotate_observed_ions
from proteomics_mcp.core.metrics import cosine_similarity, hyperscore, mean, pearson_correlation, spectral_angle, xcorr_like
from proteomics_mcp.models.schemas import CandidatePSM, MatchedPeak, ObservedSpectrum, PredictedFragment, PredictedSpectrum, SpectrumMetricResult

NEUTRON_MASS = 1.00335484


def ppm_error(observed_mz: float, predicted_mz: float) -> float:
    return ((observed_mz - predicted_mz) / predicted_mz) * 1e6


def match_predicted_to_observed(
    observed: ObservedSpectrum,
    predicted_fragments: list[PredictedFragment],
    *,
    tolerance: float = 0.05,
    tolerance_unit: str = "Da",
    one_to_one: bool = True,
) -> list[MatchedPeak]:
    """Match predicted fragments to observed peaks.

    Every (fragment, peak) pair within tolerance is a candidate, ordered closest first; ties go to the more intense
    peak, then the lower m/z, then the more intense predicted fragment.

    one_to_one (metric profile v2): pairs are assigned in that order, so each peak explains at most one fragment and
    each fragment at most one peak (#74). Without competing fragments, every fragment gets its nearest peak. When two
    fragments compete for a peak (e.g. isobaric b3/y2), the closer one gets it, and the other takes its next-closest
    unclaimed peak within tolerance, if any.

    Otherwise (legacy_v1): every fragment takes its nearest peak, even one that another fragment also takes, so a
    peak's intensity can be counted more than once.
    """
    observed_mz = np.asarray(observed.mz, dtype=float)
    observed_intensity = np.asarray(observed.intensity, dtype=float)
    order = np.argsort(observed_mz, kind="stable")
    sorted_mz = observed_mz[order]
    candidates: list[tuple[float, float, float, float, int, int]] = []
    for pred_idx, fragment in enumerate(predicted_fragments):
        window = tolerance if tolerance_unit == "Da" else abs(fragment.mz) * tolerance / 1e6
        left = int(np.searchsorted(sorted_mz, fragment.mz - window, side="left"))
        right = int(np.searchsorted(sorted_mz, fragment.mz + window, side="right"))
        for obs_idx in order[left:right]:
            obs_idx = int(obs_idx)
            obs_mz_value = float(observed_mz[obs_idx])
            err = ppm_error(obs_mz_value, fragment.mz) if tolerance_unit == "ppm" else obs_mz_value - fragment.mz
            if abs(err) <= tolerance:
                candidates.append((abs(err), -float(observed_intensity[obs_idx]), obs_mz_value, -float(fragment.intensity), pred_idx, obs_idx))
    assigned: dict[int, int] = {}
    claimed_peaks: set[int] = set()
    for *_, pred_idx, obs_idx in sorted(candidates):
        if pred_idx in assigned or (one_to_one and obs_idx in claimed_peaks):
            continue
        assigned[pred_idx] = obs_idx
        claimed_peaks.add(obs_idx)
    matches: list[MatchedPeak] = []
    for pred_idx, best_idx in sorted(assigned.items()):
        fragment = predicted_fragments[pred_idx]
        obs_mz_value = float(observed_mz[best_idx])
        err = ppm_error(obs_mz_value, fragment.mz) if tolerance_unit == "ppm" else obs_mz_value - fragment.mz
        matches.append(
            MatchedPeak(
                annotation=fragment.annotation,
                predicted_mz=fragment.mz,
                observed_mz=obs_mz_value,
                mz_error=err,
                mz_error_unit=tolerance_unit,
                predicted_intensity=fragment.intensity,
                observed_intensity=float(observed_intensity[best_idx]),
                ion_type=fragment.ion_type,
                ordinal=fragment.ordinal,
                fragment_charge=fragment.charge,
                observed_index=best_idx,
                predicted_index=pred_idx,
            )
        )
    return matches


def _precursor_signal_fraction(observed: ObservedSpectrum) -> float | None:
    if observed.precursor_mz is None or not observed.intensity:
        return None
    total = sum(observed.intensity)
    if total <= 0:
        return None
    precursor_sum = sum(i for mz, i in zip(observed.mz, observed.intensity) if abs(mz - float(observed.precursor_mz)) <= 1.5)
    return precursor_sum / total


def _isotope_intensity(observed: ObservedSpectrum, matches: list[MatchedPeak], tolerance: float) -> float:
    matched_indices = {m.observed_index for m in matches if m.observed_index is not None}
    observed_mz = np.asarray(observed.mz, dtype=float)
    order = np.argsort(observed_mz, kind="stable")
    sorted_mz = observed_mz[order]
    claimed: set[int] = set()
    extra = 0.0
    for match in matches:
        target_mz = match.predicted_mz + NEUTRON_MASS / max(match.fragment_charge, 1)
        best_idx = None
        best_delta = float("inf")
        left = int(np.searchsorted(sorted_mz, target_mz - tolerance, side="left"))
        right = int(np.searchsorted(sorted_mz, target_mz + tolerance, side="right"))
        for idx in order[left:right]:
            idx = int(idx)
            mz = float(observed_mz[idx])
            if idx in matched_indices or idx in claimed:
                continue
            delta = abs(mz - target_mz)
            if delta <= tolerance and delta < best_delta:
                best_idx = idx
                best_delta = delta
        if best_idx is not None:
            claimed.add(best_idx)
            extra += observed.intensity[best_idx]
    return extra


def similarity_vectors(
    matches: list[MatchedPeak], fragments: list[PredictedFragment], metric_profile: str = DEFAULT_METRIC_PROFILE
) -> tuple[list[float], list[float]]:
    """Return the (observed, predicted) intensity vectors that spectral angle, cosine and Pearson compare.

    legacy_v1 uses the matched predicted ions only, so a candidate matching a single peak scores 1.0 (#72).
    v2 uses every predicted ion, with zero observed intensity where no peak matched.
    """
    if metric_profile == "legacy_v1":
        return [m.observed_intensity for m in matches], [m.predicted_intensity for m in matches]
    if metric_profile == "v2":
        observed = [0.0] * len(fragments)
        for match in matches:
            observed[match.predicted_index] = match.observed_intensity
        return observed, [f.intensity for f in fragments]
    raise ValueError(f"Unknown metric_profile {metric_profile!r}; expected one of {', '.join(METRIC_PROFILES)}.")


def _similarity_metrics(
    matches: list[MatchedPeak], fragments: list[PredictedFragment], metric_profile: str
) -> tuple[float | None, float | None, float | None]:
    observed, predicted = similarity_vectors(matches, fragments, metric_profile)
    angle, cosine, pearson = spectral_angle(observed, predicted), cosine_similarity(observed, predicted), pearson_correlation(observed, predicted)
    if metric_profile == "v2" and not matches and any(value > 0 for value in predicted):
        # Nothing matched: the observed vector is all zeros, so the spectra share no signal.
        angle, cosine = 0.0, 0.0
    return angle, cosine, pearson


def scored_fragment_indices(fragments: list[PredictedFragment], metric_profile: str, min_predicted_intensity: float) -> list[int]:
    """Indices of the predicted fragments that are scored.

    Predictors give every possible ion an intensity, and ions predicted to be absent get a tiny positive one (Prosit
    down to ~1e-8, MS2PIP 2**x - 0.001) or 0. Under v2, fragments without a finite positive intensity, and those below
    min_predicted_intensity times the most intense one, are left out (#83): they would count as predicted ions, and
    could be matched to noise peaks. legacy_v1 scores every fragment, as 0.1.10 did, for every predictor.
    """
    if metric_profile == "legacy_v1":
        return list(range(len(fragments)))
    with_signal = [index for index, fragment in enumerate(fragments) if math.isfinite(fragment.intensity) and fragment.intensity > 0]
    if not with_signal or min_predicted_intensity <= 0:
        return with_signal
    threshold = min_predicted_intensity * max(fragments[index].intensity for index in with_signal)
    return [index for index in with_signal if fragments[index].intensity >= threshold]


def precursor_mass_errors(peptidoform: str, charge: int | None, observed_precursor_mz: float | None) -> tuple[float | None, float | None]:
    if charge is None or charge <= 0 or observed_precursor_mz is None or not math.isfinite(observed_precursor_mz):
        return None, None
    seq = aa_sequence_from_proforma(peptidoform)
    if seq is None:
        return None, None
    theoretical_mz = float(seq.getMZ(int(charge)))
    if theoretical_mz == 0:
        return None, None
    error_da = observed_precursor_mz - theoretical_mz
    error_ppm = error_da / theoretical_mz * 1e6
    return error_ppm, error_da


def score_candidate_spectrum(
    *,
    candidate: CandidatePSM,
    observed: ObservedSpectrum,
    predicted: PredictedSpectrum,
    tolerance: float = 0.05,
    tolerance_unit: str = "Da",
    metric_profile: str = DEFAULT_METRIC_PROFILE,
    min_predicted_intensity: float = DEFAULT_MIN_PREDICTED_INTENSITY,
    annotation_ion_types: str = "by",
    annotation_max_isotope: int = 1,
    annotation_max_ion_charge: int | None = None,
    annotation_neutral_losses: bool = False,
    use_observed_annotation: bool = True,
) -> SpectrumMetricResult:
    scored = scored_fragment_indices(predicted.fragments, metric_profile, min_predicted_intensity)
    fragments = [predicted.fragments[index] for index in scored]
    matches = match_predicted_to_observed(
        observed, fragments, tolerance=tolerance, tolerance_unit=tolerance_unit, one_to_one=metric_profile != "legacy_v1"
    )
    angle, cosine, pearson = _similarity_metrics(matches, fragments, metric_profile)
    matched_observed = [m.observed_intensity for m in matches]
    matched_predicted = [m.predicted_intensity for m in matches]
    predicted_mz = [f.mz for f in fragments]
    predicted_int = [f.intensity for f in fragments]
    xcorr = xcorr_like(observed.mz, observed.intensity, predicted_mz, predicted_int)
    log10_xcorr = math.log10(xcorr) if xcorr is not None and xcorr > 0 else None
    log2_xcorr = math.log2(xcorr) if xcorr is not None and xcorr > 0 else None
    total_intensity = sum(observed.intensity)
    matched_intensity = sum(matched_observed)
    isotope_intensity = _isotope_intensity(observed, matches, tolerance) if tolerance_unit == "Da" else 0.0
    explained = (matched_intensity + isotope_intensity) / total_intensity if total_intensity > 0 else None
    ppm_errors = [ppm_error(m.observed_mz, m.predicted_mz) for m in matches if m.predicted_mz != 0]
    da_errors = [m.observed_mz - m.predicted_mz for m in matches]
    total_predicted = len(fragments)
    charge = candidate.precursor_charge or observed.precursor_charge or predicted.precursor_charge
    observed_precursor_mz = candidate.precursor_mz or precursor_target_mz(observed)
    precursor_error_ppm, precursor_error_da = precursor_mass_errors(candidate.peptidoform, charge, observed_precursor_mz)
    hyperscore_value = hyperscore(len(matches), matched_observed, matched_predicted)
    log2_hyperscore = math.log2(hyperscore_value) if hyperscore_value is not None and hyperscore_value > 0 else None
    annotation_summary: dict[str, object] = {"engine": "internal_predicted_fragment_matcher", "matched_predicted_fragment_count": len(matches)}
    if use_observed_annotation:
        try:
            ion_annotation = annotate_observed_ions(
                observed,
                peptidoform=candidate.peptidoform,
                tolerance=tolerance,
                tolerance_unit=tolerance_unit,
                ion_types=annotation_ion_types,
                max_isotope=annotation_max_isotope,
                max_ion_charge=annotation_max_ion_charge,
                neutral_losses=annotation_neutral_losses,
            )
            annotation_summary.update(ion_annotation)
        except Exception as exc:
            annotation_summary.update({"engine": "proteomics_mcp_terminal_ion_annotator", "annotation_error": str(exc)})
    return SpectrumMetricResult(
        candidate=candidate,
        spectral_angle=angle,
        cosine_similarity=cosine,
        pearson_correlation=pearson,
        xcorr=xcorr,
        log10_xcorr=log10_xcorr,
        log2_xcorr=log2_xcorr,
        hyperscore=hyperscore_value,
        log2_hyperscore=log2_hyperscore,
        matched_ion_count=len(matches),
        total_predicted_ions=total_predicted,
        matched_ion_fraction=(len(matches) / total_predicted) if total_predicted else None,
        explained_intensity_fraction=explained,
        precursor_signal_fraction=_precursor_signal_fraction(observed),
        fragment_mass_error_mean_ppm=mean(ppm_errors),
        fragment_mass_error_abs_mean_ppm=mean(abs(v) for v in ppm_errors),
        fragment_mass_error_mean_da=mean(da_errors),
        fragment_mass_error_abs_mean_da=mean(abs(v) for v in da_errors),
        precursor_mass_error_ppm=precursor_error_ppm,
        precursor_mass_error_da=precursor_error_da,
        # predicted_index refers to predicted.fragments, which the renderers draw, not to the scored subset
        matched_peaks=[dataclasses.replace(match, predicted_index=scored[match.predicted_index]) for match in matches],
        metric_profile=metric_profile,
        predictor_metadata={"predictor_name": predicted.predictor_name, "predictor_version": predicted.predictor_version, "model": predicted.model},
        annotation_settings={
            "tolerance": tolerance,
            "tolerance_unit": tolerance_unit,
            "min_predicted_intensity": 0.0 if metric_profile == "legacy_v1" else min_predicted_intensity,
            "annotation_ion_types": annotation_ion_types,
            "annotation_max_isotope": annotation_max_isotope,
            "annotation_max_ion_charge": annotation_max_ion_charge,
            "annotation_neutral_losses": annotation_neutral_losses,
        },
        annotation_summary=annotation_summary,
        acquisition_provenance={"source": observed.source, **observed.metadata, "summary": acquisition_metadata_summary(observed, predicted)},
    )
