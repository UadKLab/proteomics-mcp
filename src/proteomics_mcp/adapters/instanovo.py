from __future__ import annotations

import ast
import math
import re
from pathlib import Path
from typing import Any

import pandas as pd

from proteomics_mcp.adapters.modifications import normalize_peptidoform_to_unimod_style
from proteomics_mcp.core.spectrum_loader import parse_scan_number
from proteomics_mcp.models.schemas import CandidatePSM


def _is_missing(value: Any) -> bool:
    return value is None or (isinstance(value, float) and math.isnan(value))


def _as_float(value: Any) -> float | None:
    if _is_missing(value) or value == "":
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _as_int(value: Any) -> int | None:
    if _is_missing(value) or value == "":
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return None


def _parse_float_list(value: Any) -> list[float] | None:
    if _is_missing(value) or value == "":
        return None
    if isinstance(value, list):
        return [_safe_float(v) for v in value if _safe_float(v) is not None]
    if isinstance(value, str):
        try:
            parsed = ast.literal_eval(value)
        except (SyntaxError, ValueError):
            return None
        if isinstance(parsed, list):
            return [_safe_float(v) for v in parsed if _safe_float(v) is not None]
    return None


def _safe_float(value: Any) -> float | None:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if math.isfinite(out) else None


def _probability_from_score(value: float | None) -> float | None:
    if value is None:
        return None
    if 0.0 <= value <= 1.0:
        return value
    if value <= 0.0:
        return max(0.0, min(1.0, math.exp(value)))
    return None


def normalize_peptidoform_to_instanovo_style(peptidoform: str | None) -> str | None:
    return normalize_peptidoform_to_unimod_style(peptidoform)


def _split_tokens(peptidoform: str, tokenised: Any = None) -> list[str]:
    if isinstance(tokenised, str) and tokenised.strip():
        try:
            parsed = ast.literal_eval(tokenised)
            if isinstance(parsed, list):
                return [str(x) for x in parsed]
        except (SyntaxError, ValueError):
            pass
        if "," in tokenised:
            return [part.strip() for part in tokenised.split(",") if part.strip()]
    tokens: list[str] = []
    i = 0
    while i < len(peptidoform):
        ch = peptidoform[i]
        if ch == "[":
            end = peptidoform.find("]", i + 1)
            if end == -1:
                i += 1
            else:
                if tokens:
                    tokens[-1] += peptidoform[i : end + 1]
                else:
                    tokens.append(peptidoform[i : end + 1])
                i = end + 1
            continue
        if ch.isalpha():
            tokens.append(ch)
        i += 1
    return tokens


def _token_summary(peptidoform: str, token_log_probs: list[float] | None, tokenised: Any = None) -> dict[str, Any]:
    if not token_log_probs:
        return {
            "token_log_probabilities": None,
            "token_log_probability_mean": None,
            "token_log_probability_min": None,
            "token_confidences": None,
            "token_confidence_mean": None,
            "token_confidence_min": None,
            "token_confidence_min_position": None,
            "token_confidence_min_token": None,
        }
    confidences = [math.exp(v) for v in token_log_probs]
    tokens = _split_tokens(peptidoform, tokenised)
    min_idx = min(range(len(confidences)), key=lambda idx: confidences[idx])
    return {
        "token_log_probabilities": token_log_probs,
        "token_log_probability_mean": sum(token_log_probs) / len(token_log_probs),
        "token_log_probability_min": min(token_log_probs),
        "token_confidences": confidences,
        "token_confidence_mean": sum(confidences) / len(confidences),
        "token_confidence_min": confidences[min_idx],
        "token_confidence_min_position": min_idx + 1,
        "token_confidence_min_token": tokens[min_idx] if min_idx < len(tokens) else None,
    }


def _detect_schema(columns: set[str]) -> str:
    if {"preds", "log_probs"}.issubset(columns):
        return "legacy_simple"
    if any(col.startswith("diffusion_") for col in columns) or any(col.startswith("instanovo_") for col in columns):
        return "documented_combined"
    if "predictions" in columns and any(col.startswith("predictions_beam_") for col in columns):
        return "instanovo_beam"
    if "predictions" in columns and "log_probs" in columns:
        return "instanovo_predictions"
    if "prediction" in columns and "log_probability" in columns:
        return "instanovo_predictions"
    raise ValueError("Could not detect an InstaNovo prediction schema from table columns.")


def _default_prediction_column(row: pd.Series, schema: str) -> tuple[str | None, str | None, Any]:
    if schema == "legacy_simple":
        return _get_str(row, "preds"), "log_probs", row.get("preds_tokenised")
    if "predictions" in row.index and not _is_missing(row.get("predictions")):
        return _get_str(row, "predictions"), "log_probs", row.get("predictions_tokenised")
    if "diffusion_predictions" in row.index and not _is_missing(row.get("diffusion_predictions")):
        return _get_str(row, "diffusion_predictions"), "diffusion_log_probabilities", row.get("diffusion_predictions_tokenised")
    if "instanovo_predictions" in row.index and not _is_missing(row.get("instanovo_predictions")):
        return _get_str(row, "instanovo_predictions"), "instanovo_log_probabilities", row.get("instanovo_predictions_tokenised")
    if "prediction" in row.index and not _is_missing(row.get("prediction")):
        return _get_str(row, "prediction"), "log_probability", row.get("prediction_tokenised")
    return None, None, None


def _get_str(row: pd.Series, column: str) -> str | None:
    if column not in row.index or _is_missing(row[column]):
        return None
    value = str(row[column]).strip()
    return value or None


def _first_float(row: pd.Series, columns: tuple[str, ...]) -> float | None:
    for column in columns:
        if column in row.index:
            value = _as_float(row.get(column))
            if value is not None:
                return value
    return None


def _beam_log_probability(row: pd.Series, beam_index: int) -> float | None:
    candidates = (
        f"predictions_log_probability_beam_{beam_index}",
        f"log_probs_beam_{beam_index}",
        f"log_probability_beam_{beam_index}",
        f"beam_{beam_index}_log_probability",
    )
    return _first_float(row, candidates)


def _candidate_id(row: pd.Series, path: str | Path, row_idx: int) -> str:
    for column in ("prediction_id", "global_index", "id", "spectrum_id"):
        if column in row.index and not _is_missing(row.get(column)):
            return str(row.get(column))
    return f"{Path(path).stem}:{row_idx}"


def _scan_number(row: pd.Series) -> int | None:
    for column in ("scan_number", "scan", "first_scan", "scan_id"):
        if column in row.index:
            value = _as_int(row.get(column))
            if value is not None:
                return value
    if "spectrum_index" in row.index:
        value = _as_int(row.get("spectrum_index"))
        if value is not None:
            # InstaNovo's labelled `*_PSM` exports retain the zero-based mzML
            # spectrum-list index, whereas mzML scan identifiers are one-based.
            sample = _get_str(row, "sample")
            return value + 1 if sample and sample.endswith("_PSM") else value
    for column in ("spectrum_id", "id", "usi"):
        if column in row.index:
            value = parse_scan_number(row.get(column))
            if value is not None:
                return value
    return None


def _spectrum_id(row: pd.Series) -> str | None:
    for column in ("spectrum_id", "id", "scan_id"):
        value = _get_str(row, column)
        if value:
            return value
    scan_number = _scan_number(row)
    if scan_number is not None:
        return str(scan_number)
    return None


def _experiment_name(row: pd.Series, default_experiment_name: str | None = None) -> str | None:
    for column in ("experiment_name", "group", "run_name", "source_name", "sample", "source_file"):
        value = _get_str(row, column)
        if value:
            if value.endswith("_PSM"):
                value = value[: -len("_PSM")]
            if value.lower().endswith((".raw", ".mzml", ".mgf")):
                value = Path(value).stem
            return value
    spectrum_id = _spectrum_id(row)
    if spectrum_id and ".raw" in spectrum_id.lower():
        return Path(spectrum_id).stem
    return default_experiment_name


def load_instanovo_predictions(
    path: str | Path,
    limit: int | None = None,
    default_experiment_name: str | None = None,
) -> list[CandidatePSM]:
    """Load legacy and current InstaNovo CSV outputs into canonical candidates.

    ``default_experiment_name`` supplies raw-file provenance for minimal exports
    that only contain spectrum indices and prediction fields.
    """
    frame = pd.read_csv(path, nrows=limit)
    schema = _detect_schema(set(frame.columns))
    candidates: list[CandidatePSM] = []
    for row_idx, row in frame.iterrows():
        peptidoform, score_col, tokenised = _default_prediction_column(row, schema)
        if not peptidoform:
            continue
        peptidoform = normalize_peptidoform_to_instanovo_style(peptidoform)
        raw_score = _as_float(row.get(score_col)) if score_col else None
        next_beam_score = _first_float(row, ("next_beam", "nextbeam", "nextscore", "next_score", "next_log_probability", "next_log_prob", "next_best_score"))
        if next_beam_score is None and schema == "instanovo_beam":
            beam0 = _beam_log_probability(row, 0)
            beam1 = _beam_log_probability(row, 1)
            if beam1 is not None:
                next_beam_score = beam1
            elif raw_score is not None and beam0 is not None and not math.isclose(raw_score, beam0):
                next_beam_score = beam0
        token_col = "token_log_probs"
        if token_col not in row.index:
            token_col = f"{score_col.replace('log_probabilities', 'token_log_probabilities')}" if score_col else "token_log_probs"
        token_log_probs = _parse_float_list(row.get(token_col)) if token_col in row.index else None
        token_data = _token_summary(peptidoform, token_log_probs, tokenised)
        confidence = _probability_from_score(raw_score)
        next_beam = _probability_from_score(next_beam_score)
        next_beam_margin = (confidence - next_beam) if confidence is not None and next_beam is not None else None
        scan_number = _scan_number(row)
        candidate_id = _candidate_id(row, path, row_idx)
        candidates.append(
            CandidatePSM(
                candidate_id=candidate_id,
                source_engine="instanovo",
                source_engine_version=None,
                experiment_name=_experiment_name(row, default_experiment_name),
                scan_number=scan_number,
                spectrum_id=_spectrum_id(row),
                usi=_get_str(row, "usi"),
                peptidoform=peptidoform,
                precursor_mz=_as_float(row.get("precursor_mz")),
                precursor_charge=_as_int(row.get("precursor_charge")),
                rank=0,
                raw_score=raw_score,
                raw_score_name=score_col,
                next_beam=next_beam,
                next_beam_margin=next_beam_margin,
                confidence=confidence,
                confidence_type="sequence_probability_from_log_probability" if confidence is not None else None,
                source_columns={"schema": schema, **{str(k): row[k] for k in row.index}},
                **token_data,
            )
        )
    return candidates
