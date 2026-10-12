from __future__ import annotations

import sys
from copy import deepcopy
from functools import lru_cache
from typing import Any

from proteomics_mcp.core.ion_annotation import parse_ion_annotation, parse_ms2pip_ion_type
from proteomics_mcp.models.schemas import PredictedFragment, PredictedSpectrum

class MissingOpenMPError(RuntimeError):
    """MS2PIP is installed but XGBoost cannot load the OpenMP runtime (libomp)."""


LIBOMP_HINT = (
    "MS2PIP could not load XGBoost because the OpenMP runtime (libomp) is missing. "
    "On macOS, install it with `brew install libomp`, then restart proteomics-mcp."
)


def _libomp_error(exc: Exception) -> MissingOpenMPError | None:
    """Return an actionable error if xgboost failed to load libomp on macOS, otherwise None."""
    if sys.platform == "darwin" and "libomp" in str(exc):
        return MissingOpenMPError(LIBOMP_HINT)
    return None


def _resolve_ms2pip_model(model: str, precursor_charge: int, max_fragment_charge: int | None) -> str:
    if max_fragment_charge is None:
        max_fragment_charge = max(1, precursor_charge - 1)
    if max_fragment_charge < 2:
        return model
    charge2_models = {"HCD": "HCDch2", "HCD2019": "HCDch2", "CID": "CIDch2"}
    return charge2_models.get(model, model)


def _predicted_spectrum_from_result(
    result: Any,
    peptidoform: str,
    precursor_charge: int,
    *,
    model: str = "HCD",
    max_fragment_charge: int | None = None,
) -> PredictedSpectrum:
    requested_max_fragment_charge = max_fragment_charge
    if max_fragment_charge is None:
        max_fragment_charge = max(1, precursor_charge - 1)
    resolved_model = _resolve_ms2pip_model(model, precursor_charge, max_fragment_charge)
    if result.theoretical_mz is None or result.predicted_intensity is None:
        raise RuntimeError(f"MS2PIP returned no fragment prediction for {peptidoform}/{precursor_charge}.")
    fragments: list[PredictedFragment] = []
    for ion_type in sorted(result.theoretical_mz.keys()):
        mzs = result.theoretical_mz[ion_type]
        intensities = result.predicted_intensity[ion_type]
        base_type, charge = parse_ms2pip_ion_type(ion_type)
        suffix = "+" * charge
        for idx, (mz, log_intensity) in enumerate(zip(mzs, intensities), start=1):
            annotation = f"{base_type}{idx}{suffix}"
            ion, ordinal, parsed_charge = parse_ion_annotation(annotation)
            if parsed_charge > max_fragment_charge:
                continue
            linear_intensity = max(0.0, (2 ** float(log_intensity)) - 0.001)
            fragments.append(
                PredictedFragment(
                    mz=float(mz),
                    intensity=linear_intensity,
                    annotation=annotation,
                    ion_type=ion,
                    ordinal=ordinal,
                    charge=parsed_charge,
                )
            )
    if fragments:
        max_intensity = max(f.intensity for f in fragments) or 1.0
        for fragment in fragments:
            fragment.intensity = fragment.intensity / max_intensity
    fragments.sort(key=lambda f: f.mz)
    return PredictedSpectrum(
        peptidoform=peptidoform,
        precursor_charge=precursor_charge,
        fragments=fragments,
        predictor_name="ms2pip",
        predictor_version=None,
        model=resolved_model,
        metadata={
            "requested_model": model,
            "resolved_model": resolved_model,
            "requested_max_fragment_charge": requested_max_fragment_charge,
            "max_fragment_charge": max_fragment_charge,
            "fragment_charge_rule": "charge <= precursor_charge - 1" if requested_max_fragment_charge is None else "explicit",
            "available_fragment_charges": sorted({fragment.charge for fragment in fragments}),
        },
    )


@lru_cache(maxsize=1024)
def _predict_fragments_ms2pip_cached(
    peptidoform: str,
    precursor_charge: int,
    model: str,
    max_fragment_charge: int | None,
) -> PredictedSpectrum:
    try:
        from ms2pip.core import predict_single  # type: ignore
    except ImportError as exc:
        raise RuntimeError("MS2PIP prediction requires ms2pip to be installed.") from exc
    except Exception as exc:  # xgboost raises XGBoostError, not ImportError, when libomp is missing
        actionable = _libomp_error(exc)
        if actionable is None:
            raise
        raise actionable from exc
    requested_max_fragment_charge = max_fragment_charge
    resolved_model = _resolve_ms2pip_model(model, precursor_charge, max_fragment_charge)
    result = predict_single(f"{peptidoform}/{precursor_charge}", model=resolved_model)
    return _predicted_spectrum_from_result(
        result,
        peptidoform,
        precursor_charge,
        model=model,
        max_fragment_charge=requested_max_fragment_charge,
    )


def predict_fragments_ms2pip(
    peptidoform: str,
    precursor_charge: int,
    *,
    model: str = "HCD",
    max_fragment_charge: int | None = None,
) -> PredictedSpectrum:
    """Predict one fragment spectrum with local MS2PIP and a bounded process cache."""
    return deepcopy(_predict_fragments_ms2pip_cached(peptidoform, precursor_charge, model, max_fragment_charge))


def predict_fragments_ms2pip_batch(
    requests: list[tuple[str, int]],
    *,
    model: str = "HCD",
    max_fragment_charge: int | None = None,
    processes: int | None = None,
) -> list[PredictedSpectrum]:
    """Predict a batch with MS2PIP while retaining input order and per-charge models."""
    if not requests:
        return []
    try:
        from ms2pip.core import predict_batch  # type: ignore
        from psm_utils import PSM, PSMList, Peptidoform  # type: ignore
    except ImportError as exc:
        raise RuntimeError("MS2PIP batch prediction requires ms2pip and psm-utils to be installed.") from exc
    except Exception as exc:  # xgboost raises XGBoostError, not ImportError, when libomp is missing
        actionable = _libomp_error(exc)
        if actionable is None:
            raise
        raise actionable from exc

    grouped: dict[str, list[tuple[int, str, int]]] = {}
    for index, (peptidoform, precursor_charge) in enumerate(requests):
        resolved_model = _resolve_ms2pip_model(model, precursor_charge, max_fragment_charge)
        grouped.setdefault(resolved_model, []).append((index, peptidoform, precursor_charge))

    predictions: list[PredictedSpectrum | None] = [None] * len(requests)
    for resolved_model, entries in grouped.items():
        psms = PSMList(
            psm_list=[
                PSM(peptidoform=Peptidoform(f"{peptidoform}/{precursor_charge}"), spectrum_id=str(index))
                for index, peptidoform, precursor_charge in entries
            ]
        )
        results = predict_batch(psms, model=resolved_model, processes=processes)
        if len(results) != len(entries):
            raise RuntimeError(f"MS2PIP returned {len(results)} results for {len(entries)} requests.")
        for (index, peptidoform, precursor_charge), result in zip(entries, results, strict=True):
            predictions[index] = _predicted_spectrum_from_result(
                result,
                peptidoform,
                precursor_charge,
                model=model,
                max_fragment_charge=max_fragment_charge,
            )
    if any(prediction is None for prediction in predictions):
        raise RuntimeError("MS2PIP did not return a prediction for every request.")
    return [prediction for prediction in predictions if prediction is not None]
