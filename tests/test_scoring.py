from __future__ import annotations

import math

import numpy as np
import pytest

from proteomics_mcp.config import load_config
from proteomics_mcp.core.acquisition import acquisition_metadata_summary
from proteomics_mcp.core.scoring import match_predicted_to_observed, precursor_mass_errors, score_candidate_spectrum
from proteomics_mcp.models.schemas import CandidatePSM, ObservedSpectrum, PredictedFragment, PredictedSpectrum


def test_match_uses_closest_peak_within_da_tolerance():
    observed = ObservedSpectrum(mz=[99.99, 100.01, 100.015], intensity=[10, 50, 100])
    fragments = [PredictedFragment(mz=100.0, intensity=1.0, annotation="b1+", ion_type="b", ordinal=1)]

    [match] = match_predicted_to_observed(observed, fragments, tolerance=0.02, tolerance_unit="Da")

    assert match.observed_mz == 100.01
    assert match.observed_intensity == 50


def test_score_candidate_spectrum_returns_broad_metrics():
    candidate = CandidatePSM(candidate_id="c1", source_engine="test", peptidoform="PEPTIDE", precursor_charge=2)
    observed = ObservedSpectrum(mz=[100.0, 200.0, 300.0], intensity=[10.0, 20.0, 30.0], precursor_mz=300.0, precursor_charge=2)
    predicted = PredictedSpectrum(
        peptidoform="PEPTIDE",
        precursor_charge=2,
        predictor_name="mock",
        fragments=[
            PredictedFragment(mz=100.0, intensity=0.5, annotation="b1+", ion_type="b", ordinal=1),
            PredictedFragment(mz=200.0, intensity=1.0, annotation="y1+", ion_type="y", ordinal=1),
        ],
    )

    result = score_candidate_spectrum(candidate=candidate, observed=observed, predicted=predicted)

    assert result.metric_profile == "v2"
    assert result.matched_ion_count == 2
    assert result.total_predicted_ions == 2
    assert result.matched_ion_fraction == 1.0
    assert result.spectral_angle is not None
    assert result.cosine_similarity is not None
    assert result.xcorr is not None
    assert result.log2_xcorr is not None
    assert result.hyperscore is not None
    assert result.log2_hyperscore is not None
    assert result.precursor_signal_fraction == 0.5
    assert result.precursor_mass_error_ppm is not None
    assert result.annotation_summary["engine"] == "proteomics_mcp_terminal_ion_annotator"
    assert "ion_type_counts" in result.annotation_summary


def test_precursor_mass_error_is_zero_for_theoretical_openms_mz():
    import pyopenms  # type: ignore

    theoretical_mz = float(pyopenms.AASequence.fromString("PEPTIDE").getMZ(2))

    ppm, da = precursor_mass_errors("PEPTIDE", 2, theoretical_mz)

    assert ppm is not None
    assert da is not None
    assert abs(ppm) < 1e-9
    assert abs(da) < 1e-12


def test_score_candidate_uses_precursor_target_metadata_for_precursor_error():
    import pyopenms  # type: ignore

    theoretical_mz = float(pyopenms.AASequence.fromString("PEPTIDE").getMZ(2))
    candidate = CandidatePSM(candidate_id="c1", source_engine="test", peptidoform="PEPTIDE", precursor_charge=2)
    observed = ObservedSpectrum(
        mz=[100.0],
        intensity=[10.0],
        precursor_charge=2,
        metadata={"attributes": [{"accession": "MS:1000827", "name": "isolation window target m/z", "value": theoretical_mz}]},
    )
    predicted = PredictedSpectrum(
        peptidoform="PEPTIDE",
        precursor_charge=2,
        predictor_name="mock",
        fragments=[PredictedFragment(mz=100.0, intensity=1.0, annotation="b1+", ion_type="b", ordinal=1)],
    )

    result = score_candidate_spectrum(candidate=candidate, observed=observed, predicted=predicted)

    assert result.precursor_mass_error_ppm is not None
    assert result.precursor_mass_error_da is not None
    assert abs(result.precursor_mass_error_ppm) < 1e-9
    assert abs(result.precursor_mass_error_da) < 1e-12


def test_acquisition_separates_instrument_model_from_detector_type():
    observed = ObservedSpectrum(
        mz=[100.0],
        intensity=[10.0],
        metadata={
            "attributes": [
                {"accession": "MS:1000463", "name": "instrument", "value": "Orbitrap Fusion Lumos"},
                {"accession": "MS:1000484", "name": "orbitrap"},
            ]
        },
    )

    summary = acquisition_metadata_summary(observed)

    assert summary["instrument_model"] == "Orbitrap Fusion Lumos"
    assert summary["detector_type"] == "orbitrap"


def test_acquisition_does_not_treat_instrument_model_as_detector_type():
    observed = ObservedSpectrum(
        mz=[100.0],
        intensity=[10.0],
        metadata={"attributes": [{"accession": "MS:1002732", "name": "Orbitrap Fusion Lumos"}]},
    )

    summary = acquisition_metadata_summary(observed)

    assert summary["instrument_model"] == "Orbitrap Fusion Lumos"
    assert summary["detector_type"] is None


def _ten_y_ions() -> PredictedSpectrum:
    fragments = [
        PredictedFragment(mz=100.0 + 100.0 * i, intensity=1.0, annotation=f"y{i}+", ion_type="y", ordinal=i) for i in range(1, 11)
    ]
    return PredictedSpectrum(peptidoform="PEPTIDEK", precursor_charge=2, predictor_name="mock", fragments=fragments)


def _score(observed: ObservedSpectrum, predicted: PredictedSpectrum, metric_profile: str):
    candidate = CandidatePSM(candidate_id="c1", source_engine="test", peptidoform="PEPTIDEK", precursor_charge=2)
    return score_candidate_spectrum(
        candidate=candidate, observed=observed, predicted=predicted, metric_profile=metric_profile, use_observed_annotation=False
    )


def test_similarity_counts_unmatched_predicted_ions_as_zero():
    # One of ten equally intense predicted ions matches a peak (#72).
    observed = ObservedSpectrum(mz=[200.0], intensity=[50.0])

    result = _score(observed, _ten_y_ions(), "v2")

    assert result.matched_ion_count == 1
    assert result.cosine_similarity == pytest.approx(1 / math.sqrt(10))
    assert result.spectral_angle == pytest.approx(1 - 2 * math.acos(1 / math.sqrt(10)) / math.pi)
    assert result.spectral_angle < 0.25


def test_pearson_includes_unmatched_predicted_ions():
    predicted = _ten_y_ions()
    for i, fragment in enumerate(predicted.fragments, start=1):
        fragment.intensity = float(i)
    observed = ObservedSpectrum(mz=[1000.0, 1100.0], intensity=[40.0, 50.0])

    result = _score(observed, predicted, "v2")

    expected = np.corrcoef([0.0] * 8 + [40.0, 50.0], [float(i) for i in range(1, 11)])[0, 1]
    assert result.pearson_correlation == pytest.approx(expected)
    assert _score(observed, predicted, "legacy_v1").pearson_correlation == pytest.approx(1.0)


def test_legacy_v1_profile_compares_matched_ions_only():
    observed = ObservedSpectrum(mz=[200.0], intensity=[50.0])

    result = _score(observed, _ten_y_ions(), "legacy_v1")

    assert result.metric_profile == "legacy_v1"
    assert result.spectral_angle == pytest.approx(1.0)
    assert result.cosine_similarity == pytest.approx(1.0)


def test_similarity_is_one_when_every_predicted_ion_matches_in_proportion():
    observed = ObservedSpectrum(mz=[100.0 + 100.0 * i for i in range(1, 11)], intensity=[7.0] * 10)

    for profile in ("v2", "legacy_v1"):
        result = _score(observed, _ten_y_ions(), profile)
        assert result.spectral_angle == pytest.approx(1.0)
        assert result.cosine_similarity == pytest.approx(1.0)


def test_similarity_without_matches_is_zero():
    observed = ObservedSpectrum(mz=[55.0], intensity=[50.0])

    result = _score(observed, _ten_y_ions(), "v2")

    assert result.matched_ion_count == 0
    assert result.spectral_angle == 0.0
    assert result.cosine_similarity == 0.0
    assert result.pearson_correlation is None


def test_unknown_metric_profile_is_rejected():
    with pytest.raises(ValueError, match="metric_profile"):
        load_config(overrides={"metric_profile": "v3"})
    with pytest.raises(ValueError, match="metric_profile"):
        _score(ObservedSpectrum(mz=[200.0], intensity=[50.0]), _ten_y_ions(), "v3")


def test_default_metric_profile_is_v2():
    assert load_config().metric_profile == "v2"


def _two_fragments(b_mz: float, y_mz: float) -> PredictedSpectrum:
    return PredictedSpectrum(
        peptidoform="PEPTIDEK",
        precursor_charge=2,
        predictor_name="mock",
        fragments=[
            PredictedFragment(mz=b_mz, intensity=1.0, annotation="b4+", ion_type="b", ordinal=4),
            PredictedFragment(mz=y_mz, intensity=0.5, annotation="y8++", ion_type="y", ordinal=8, charge=2),
        ],
    )


def test_one_peak_explains_at_most_one_predicted_ion():
    # Two predicted ions within tolerance of a single peak (#74).
    observed = ObservedSpectrum(mz=[500.01], intensity=[100.0])

    result = _score(observed, _two_fragments(500.00, 500.02), "v2")

    assert result.matched_ion_count == 1
    assert result.explained_intensity_fraction == pytest.approx(1.0)
    assert len({m.observed_index for m in result.matched_peaks}) == len(result.matched_peaks)


def test_legacy_v1_matches_a_peak_to_every_ion_in_tolerance():
    observed = ObservedSpectrum(mz=[500.01], intensity=[100.0])

    result = _score(observed, _two_fragments(500.00, 500.02), "legacy_v1")

    assert result.matched_ion_count == 2
    assert result.explained_intensity_fraction == pytest.approx(2.0)


def test_competing_ion_takes_its_next_closest_unclaimed_peak():
    observed = ObservedSpectrum(mz=[500.01, 500.04], intensity=[100.0, 30.0])

    matches = match_predicted_to_observed(observed, _two_fragments(500.00, 500.012).fragments, tolerance=0.05)

    by_annotation = {m.annotation: m.observed_mz for m in matches}
    assert by_annotation == {"y8++": 500.01, "b4+": 500.04}


def test_one_to_one_matching_keeps_nearest_peaks_without_competition():
    observed = ObservedSpectrum(mz=[300.0, 300.02, 600.0], intensity=[10.0, 50.0, 20.0])
    fragments = _two_fragments(300.01, 600.01).fragments

    one_to_one = match_predicted_to_observed(observed, fragments, tolerance=0.05)
    nearest = match_predicted_to_observed(observed, fragments, tolerance=0.05, one_to_one=False)

    assert [(m.annotation, m.observed_index) for m in one_to_one] == [(m.annotation, m.observed_index) for m in nearest]


def _with_near_absent_ion() -> PredictedSpectrum:
    # y2 is predicted to be absent (1e-5 of the base peak), as predictors report for most possible ions.
    return PredictedSpectrum(
        peptidoform="PEPTIDEK",
        precursor_charge=2,
        predictor_name="mock",
        fragments=[
            PredictedFragment(mz=200.0, intensity=1.0, annotation="y1+", ion_type="y", ordinal=1),
            PredictedFragment(mz=300.0, intensity=1e-5, annotation="y2+", ion_type="y", ordinal=2),
            PredictedFragment(mz=400.0, intensity=0.5, annotation="y3+", ion_type="y", ordinal=3),
        ],
    )


def test_near_absent_predicted_ions_are_not_scored():
    # The noise peak at 300.0 would otherwise count as a matched ion and as explained intensity (#83).
    observed = ObservedSpectrum(mz=[200.0, 300.0, 400.0], intensity=[100.0, 40.0, 50.0])

    result = _score(observed, _with_near_absent_ion(), "v2")

    assert result.total_predicted_ions == 2
    assert result.matched_ion_count == 2
    assert result.explained_intensity_fraction == pytest.approx(150.0 / 190.0)
    # matched peaks still index the full predicted spectrum, which renderers draw
    assert sorted(m.predicted_index for m in result.matched_peaks) == [0, 2]
    assert result.annotation_settings["min_predicted_intensity"] == 0.001


def test_legacy_v1_scores_every_predicted_ion():
    observed = ObservedSpectrum(mz=[200.0, 300.0, 400.0], intensity=[100.0, 40.0, 50.0])

    result = _score(observed, _with_near_absent_ion(), "legacy_v1")

    assert result.total_predicted_ions == 3
    assert result.matched_ion_count == 3
    assert result.explained_intensity_fraction == pytest.approx(1.0)


def test_min_predicted_intensity_zero_scores_every_ion():
    observed = ObservedSpectrum(mz=[200.0, 300.0, 400.0], intensity=[100.0, 40.0, 50.0])
    candidate = CandidatePSM(candidate_id="c1", source_engine="test", peptidoform="PEPTIDEK", precursor_charge=2)

    result = score_candidate_spectrum(
        candidate=candidate, observed=observed, predicted=_with_near_absent_ion(), min_predicted_intensity=0.0, use_observed_annotation=False
    )

    assert result.total_predicted_ions == 3


@pytest.mark.parametrize("value", [-0.1, 1.0, 2.0])
def test_min_predicted_intensity_must_be_a_fraction(value):
    with pytest.raises(ValueError, match="min_predicted_intensity"):
        load_config(overrides={"min_predicted_intensity": value})
