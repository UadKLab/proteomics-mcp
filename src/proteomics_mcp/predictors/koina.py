from __future__ import annotations

import json
import math
from functools import lru_cache
from pathlib import Path
import re
import tempfile
from typing import Any
from urllib import parse, request

from proteomics_mcp.adapters.modifications import normalize_peptidoform_to_unimod_style
from proteomics_mcp.core.ion_annotation import parse_ion_annotation
from proteomics_mcp.models.schemas import PredictedFragment, PredictedSpectrum

_DEFAULT_BATCH_SIZE = 64


def _urlopen(req: request.Request, timeout: float) -> Any:
    return request.urlopen(req, timeout=timeout)


def _koina_cache_dir(cache_dir: str | None) -> Path:
    if cache_dir:
        path = Path(cache_dir).expanduser()
    else:
        path = Path(tempfile.gettempdir()) / "proteomics_mcp_koina_cache"
    path.mkdir(parents=True, exist_ok=True)
    return path


def _cache_key(
    *,
    peptidoform: str,
    precursor_charge: int,
    model: str,
    collision_energy: float,
    fragmentation_type: str,
    instrument_type: str,
    max_fragment_charge: int | None,
) -> str:
    payload = {
        "schema_version": 2,
        "peptidoform": peptidoform,
        "precursor_charge": precursor_charge,
        "model": model,
        "collision_energy": collision_energy,
        "fragmentation_type": fragmentation_type,
        "instrument_type": instrument_type,
        "max_fragment_charge": max_fragment_charge,
    }
    return json.dumps(payload, sort_keys=True, separators=(",", ":"))


def _cache_path(cache_dir: str | None, key: str) -> Path:
    import hashlib

    digest = hashlib.sha1(key.encode("utf-8")).hexdigest()
    return _koina_cache_dir(cache_dir) / f"{digest}.json"


def _count_residues(peptidoform: str) -> int:
    return sum(1 for char in peptidoform if char.isalpha() and char.isupper())


def _validate_prosit_constraints(peptidoform: str, precursor_charge: int, model: str) -> None:
    if not model.startswith("Prosit_"):
        return
    residue_count = _count_residues(peptidoform)
    if residue_count > 30:
        raise ValueError(f"{model} currently expects peptide length <= 30 residues; got {residue_count}.")
    if precursor_charge > 6:
        raise ValueError(f"{model} currently expects precursor charge <= 6; got {precursor_charge}.")


def _normalize_peptidoform_for_koina(peptidoform: str) -> str:
    normalized = normalize_peptidoform_to_unimod_style(peptidoform) or peptidoform
    return normalized.replace("--", "-")


def _infer_triton_shape(data: Any) -> list[int]:
    if isinstance(data, list):
        if not data:
            return [0]
        if isinstance(data[0], list):
            return [len(data)] + _infer_triton_shape(data[0])
        return [len(data)]
    return [1]


def _ensure_row_lists(values: Any, batch_size: int) -> list[list[Any]]:
    if not isinstance(values, list):
        return [[values] for _ in range(batch_size)]
    if batch_size == 1 and (not values or not isinstance(values[0], list)):
        return [values]
    if len(values) == batch_size and all(isinstance(entry, list) for entry in values):
        return values
    if len(values) == batch_size and not any(isinstance(entry, list) for entry in values):
        return [[entry] for entry in values]
    return [values[i] if i < len(values) and isinstance(values[i], list) else [values[i]] for i in range(min(batch_size, len(values)))]


def _reshape_flat_values(values: list[Any], shape: list[int]) -> Any:
    if not shape:
        return values[0] if values else None
    if len(shape) == 1:
        return values[: shape[0]]
    chunk = 1
    for dim in shape[1:]:
        chunk *= dim
    rows = []
    for index in range(shape[0]):
        start = index * chunk
        rows.append(_reshape_flat_values(values[start : start + chunk], shape[1:]))
    return rows


def _decode_outputs(payload: dict[str, Any]) -> dict[str, Any]:
    outputs: dict[str, Any] = {}
    for output in payload.get("outputs", []):
        if not isinstance(output, dict):
            continue
        name = str(output.get("name") or "")
        data = output.get("data")
        shape = output.get("shape")
        if name and isinstance(data, list) and isinstance(shape, list):
            expected = 1
            for dim in shape:
                expected *= int(dim)
            if expected == len(data):
                outputs[name] = _reshape_flat_values(data, [int(dim) for dim in shape])
                continue
        if name:
            outputs[name] = data
    return outputs


def _first_present(mapping: dict[str, Any], names: list[str]) -> Any:
    for name in names:
        if name in mapping and mapping[name] is not None:
            return mapping[name]
    return None


def _build_annotation(ion_type: Any, ordinal: Any, charge: Any) -> str:
    ion = str(ion_type)
    number = int(ordinal)
    state = int(charge) if charge is not None else 1
    return f"{ion}{number}{'+' * max(state, 1)}"


def _has_predicted_signal(intensity: Any) -> bool:
    """Prosit-style models predict an intensity for every annotated ion; 0 (or less) means the ion is not expected.

    Such fragments are left out (#83): they add nothing to the predicted spectrum, but would count as predicted ions
    and could be matched to noise peaks.
    """
    value = float(intensity)
    return math.isfinite(value) and value > 0


def _predicted_spectrum_from_rows(
    *,
    peptidoform: str,
    precursor_charge: int,
    model: str,
    max_fragment_charge: int | None,
    mz_row: list[Any],
    intensity_row: list[Any],
    annotation_row: list[Any] | None,
    ion_type_row: list[Any] | None,
    ordinal_row: list[Any] | None,
    charge_row: list[Any] | None,
    collision_energy: float,
    fragmentation_type: str,
    instrument_type: str,
    server_url: str,
    cache_hit: bool,
    model_version: str | None = None,
) -> PredictedSpectrum:
    if annotation_row is None:
        if ion_type_row is None or ordinal_row is None:
            raise ValueError("Koina response did not include fragment annotations or enough fields to reconstruct them.")
        if charge_row is None:
            charge_row = [1] * len(mz_row)
        annotation_row = [_build_annotation(ion, ordinal, fragment_charge) for ion, ordinal, fragment_charge in zip(ion_type_row, ordinal_row, charge_row)]
    if len(mz_row) != len(intensity_row) or len(mz_row) != len(annotation_row):
        raise ValueError("Koina fragment m/z, intensity, and annotation arrays must have equal lengths.")
    fragments: list[PredictedFragment] = []
    dropped_without_signal = 0
    allowed_charge = max_fragment_charge if max_fragment_charge is not None else max(precursor_charge - 1, 1)
    for mz, intensity, annotation in zip(mz_row, intensity_row, annotation_row):
        if float(mz) <= 0:
            # Prosit-style outputs use a negative m/z sentinel for unavailable ions.
            continue
        ion_type, ordinal, fragment_charge = parse_ion_annotation(str(annotation))
        if fragment_charge > allowed_charge:
            continue
        if not _has_predicted_signal(intensity):
            dropped_without_signal += 1
            continue
        fragments.append(
            PredictedFragment(
                mz=float(mz),
                intensity=float(intensity),
                annotation=str(annotation),
                ion_type=ion_type,
                ordinal=ordinal,
                charge=fragment_charge,
            )
        )
    if fragments:
        max_intensity = max(fragment.intensity for fragment in fragments) or 1.0
        for fragment in fragments:
            fragment.intensity = fragment.intensity / max_intensity
    fragments.sort(key=lambda fragment: fragment.mz)
    return PredictedSpectrum(
        peptidoform=peptidoform,
        precursor_charge=precursor_charge,
        fragments=fragments,
        predictor_name="koina",
        # Koina serves versioned Triton models; record the version that answered, e.g. "Prosit_2020_intensity_HCD@1"
        predictor_version=f"{model}@{model_version}" if model_version else None,
        model=model,
        metadata={
            "server_url": server_url,
            "model_version": model_version,
            "dropped_zero_intensity_fragments": dropped_without_signal,
            "collision_energy": collision_energy,
            "fragmentation_type": fragmentation_type,
            "instrument_type": instrument_type,
            "max_fragment_charge": allowed_charge,
            "fragment_charge_rule": "charge <= precursor_charge - 1" if max_fragment_charge is None else "explicit max_fragment_charge",
            "cache_hit": cache_hit,
            "available_fragment_charges": sorted({fragment.charge for fragment in fragments}),
        },
    )


def _format_input_values(values: list[Any], shape: list[Any]) -> tuple[list[int], list[Any]]:
    """Format batch values for the one-column tensor convention used by Koina."""
    dimensions = [int(dimension) for dimension in shape if int(dimension) != -1]
    if len(shape) >= 2:
        width = dimensions[-1] if dimensions else 1
        if width != 1:
            raise ValueError(f"Only scalar per-spectrum Koina inputs are supported; received shape {shape}.")
        return [len(values), 1], [[value] for value in values]
    return [len(values)], values


def _koina_request_payload(
    items: list[tuple[str, int]],
    *,
    model_inputs: list[dict[str, Any]],
    collision_energy: float,
    fragmentation_type: str,
    instrument_type: str,
) -> dict[str, Any]:
    sequences = [sequence for sequence, _ in items]
    charges = [charge for _, charge in items]
    values_by_name = {
        "peptide_sequences": sequences,
        "precursor_charges": charges,
        "collision_energies": [float(collision_energy)] * len(items),
        "fragmentation_types": [fragmentation_type] * len(items),
        "instrument_types": [instrument_type] * len(items),
    }
    inputs: list[dict[str, Any]] = []
    for model_input in model_inputs:
        name = str(model_input.get("name") or "")
        if name not in values_by_name:
            raise ValueError(f"Koina model requires unsupported input '{name}'.")
        shape, data = _format_input_values(values_by_name[name], list(model_input.get("shape") or [-1]))
        inputs.append(
            {
                "name": name,
                "shape": shape,
                "datatype": str(model_input.get("datatype") or "BYTES"),
                "data": data,
            }
        )
    if not inputs:
        raise ValueError("Koina model metadata did not advertise any inputs.")
    return {"inputs": inputs}


@lru_cache(maxsize=64)
def _describe_koina_model_cached(model: str, server_url: str, timeout_seconds: float) -> dict[str, Any]:
    endpoint = server_url.rstrip("/") + f"/v2/models/{parse.quote(model, safe='')}"
    req = request.Request(endpoint, headers={"Accept": "application/json"})
    with _urlopen(req, timeout_seconds) as response:
        payload = json.loads(response.read().decode("utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Koina returned invalid metadata for model '{model}'.")
    return payload


def describe_koina_model(
    *,
    model: str,
    server_url: str = "https://koina.wilhelmlab.org",
    timeout_seconds: float = 30.0,
) -> dict[str, Any]:
    """Return live Triton metadata for one Koina model without running inference."""
    return _describe_koina_model_cached(model, server_url.rstrip("/"), float(timeout_seconds))


def _request_koina_predictions(
    items: list[tuple[str, int]],
    *,
    model: str,
    server_url: str,
    collision_energy: float,
    fragmentation_type: str,
    instrument_type: str,
    timeout_seconds: float,
) -> dict[str, Any]:
    model_metadata = describe_koina_model(
        model=model,
        server_url=server_url,
        timeout_seconds=timeout_seconds,
    )
    payload = _koina_request_payload(
        items,
        model_inputs=model_metadata.get("inputs", []),
        collision_energy=collision_energy,
        fragmentation_type=fragmentation_type,
        instrument_type=instrument_type,
    )
    endpoint = server_url.rstrip("/") + f"/v2/models/{parse.quote(model, safe='')}/infer"
    encoded = json.dumps(payload).encode("utf-8")
    req = request.Request(endpoint, data=encoded, headers={"Content-Type": "application/json"})
    with _urlopen(req, timeout_seconds) as response:
        return json.loads(response.read().decode("utf-8"))


def predict_fragments_koina_batch(
    items: list[tuple[str, int]],
    *,
    model: str,
    max_fragment_charge: int | None = None,
    server_url: str = "https://koina.wilhelmlab.org",
    collision_energy: float = 30.0,
    fragmentation_type: str = "HCD",
    instrument_type: str = "QE",
    timeout_seconds: float = 30.0,
    batch_size: int = _DEFAULT_BATCH_SIZE,
    cache_dir: str | None = None,
) -> list[PredictedSpectrum]:
    normalized_items = [(_normalize_peptidoform_for_koina(peptidoform), int(charge)) for peptidoform, charge in items]
    for peptidoform, precursor_charge in normalized_items:
        _validate_prosit_constraints(peptidoform, precursor_charge, model)

    cache_hits: dict[int, PredictedSpectrum] = {}
    pending: list[tuple[int, str, int, Path]] = []
    for index, (peptidoform, precursor_charge) in enumerate(normalized_items):
        key = _cache_key(
            peptidoform=peptidoform,
            precursor_charge=precursor_charge,
            model=model,
            collision_energy=float(collision_energy),
            fragmentation_type=fragmentation_type,
            instrument_type=instrument_type,
            max_fragment_charge=max_fragment_charge,
        )
        path = _cache_path(cache_dir, key)
        payload = json.loads(path.read_text(encoding="utf-8")) if path.exists() else None
        # entries cached before model versions were recorded (<= 0.1.9) are predicted again
        if payload is not None and "model_version" in (payload.get("metadata") or {}):
            cached = [PredictedFragment(**fragment) for fragment in payload.get("fragments", [])]
            # entries cached by <= 0.1.10 still contain zero-intensity fragments
            fragments = [fragment for fragment in cached if _has_predicted_signal(fragment.intensity)]
            metadata = dict(payload.get("metadata") or {})
            metadata["dropped_zero_intensity_fragments"] = metadata.get("dropped_zero_intensity_fragments", 0) + len(cached) - len(fragments)
            cache_hits[index] = PredictedSpectrum(
                peptidoform=payload["peptidoform"],
                precursor_charge=int(payload["precursor_charge"]),
                fragments=fragments,
                predictor_name=payload.get("predictor_name", "koina"),
                predictor_version=payload.get("predictor_version"),
                model=payload.get("model"),
                metadata={**metadata, "cache_hit": True},
            )
        else:
            pending.append((index, peptidoform, precursor_charge, path))

    predictions: list[PredictedSpectrum | None] = [None] * len(items)
    for index, predicted in cache_hits.items():
        predictions[index] = predicted

    for offset in range(0, len(pending), max(int(batch_size), 1)):
        chunk = pending[offset : offset + max(int(batch_size), 1)]
        response = _request_koina_predictions(
            [(peptidoform, precursor_charge) for _, peptidoform, precursor_charge, _ in chunk],
            model=model,
            server_url=server_url,
            collision_energy=collision_energy,
            fragmentation_type=fragmentation_type,
            instrument_type=instrument_type,
            timeout_seconds=timeout_seconds,
        )
        outputs = _decode_outputs(response)
        mz_values = _first_present(outputs, ["mz", "mz_values", "mzs"])
        intensity_values = _first_present(outputs, ["intensities", "intensity", "intensity_values"])
        if mz_values is None or intensity_values is None:
            raise ValueError("Koina response must include both fragment m/z and intensity outputs.")
        mz_rows = _ensure_row_lists(mz_values, len(chunk))
        intensity_rows = _ensure_row_lists(intensity_values, len(chunk))
        annotation_rows = _first_present(outputs, ["annotation", "annotations", "ion_annotations"])
        ion_type_rows = _first_present(outputs, ["ion_types", "fragment_types"])
        ordinal_rows = _first_present(outputs, ["ordinals", "fragment_numbers", "fragment_ordinals"])
        charge_rows = _first_present(outputs, ["charges", "fragment_charges"])
        if annotation_rows is not None:
            annotation_rows = _ensure_row_lists(annotation_rows, len(chunk))
        if ion_type_rows is not None:
            ion_type_rows = _ensure_row_lists(ion_type_rows, len(chunk))
        if ordinal_rows is not None:
            ordinal_rows = _ensure_row_lists(ordinal_rows, len(chunk))
        if charge_rows is not None:
            charge_rows = _ensure_row_lists(charge_rows, len(chunk))
        if len(mz_rows) != len(chunk) or len(intensity_rows) != len(chunk):
            raise ValueError("Koina response batch size does not match the submitted items.")

        for row_index, (prediction_index, peptidoform, precursor_charge, cache_path) in enumerate(chunk):
            predicted = _predicted_spectrum_from_rows(
                peptidoform=peptidoform,
                precursor_charge=precursor_charge,
                model=model,
                max_fragment_charge=max_fragment_charge,
                mz_row=mz_rows[row_index],
                intensity_row=intensity_rows[row_index],
                annotation_row=annotation_rows[row_index] if annotation_rows is not None else None,
                ion_type_row=ion_type_rows[row_index] if ion_type_rows is not None else None,
                ordinal_row=ordinal_rows[row_index] if ordinal_rows is not None else None,
                charge_row=charge_rows[row_index] if charge_rows is not None else None,
                collision_energy=collision_energy,
                fragmentation_type=fragmentation_type,
                instrument_type=instrument_type,
                server_url=server_url,
                cache_hit=False,
                model_version=str(response["model_version"]) if response.get("model_version") else None,
            )
            predictions[prediction_index] = predicted
            cache_path.write_text(json.dumps(predicted.model_dump(), indent=2), encoding="utf-8")

    if any(prediction is None for prediction in predictions):
        raise ValueError("Koina returned no prediction for one or more submitted items.")
    return [prediction for prediction in predictions if prediction is not None]


def predict_fragments_koina(
    peptidoform: str,
    precursor_charge: int,
    *,
    model: str = "Prosit_2020_intensity_HCD",
    max_fragment_charge: int | None = None,
    server_url: str = "https://koina.wilhelmlab.org",
    collision_energy: float = 30.0,
    fragmentation_type: str = "HCD",
    instrument_type: str = "QE",
    timeout_seconds: float = 30.0,
    cache_dir: str | None = None,
) -> PredictedSpectrum:
    return predict_fragments_koina_batch(
        [(peptidoform, precursor_charge)],
        model=model,
        max_fragment_charge=max_fragment_charge,
        server_url=server_url,
        collision_energy=collision_energy,
        fragmentation_type=fragmentation_type,
        instrument_type=instrument_type,
        timeout_seconds=timeout_seconds,
        batch_size=1,
        cache_dir=cache_dir,
    )[0]
