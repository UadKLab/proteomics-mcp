from __future__ import annotations

import ast
import logging
import math
import re
from pathlib import Path
from typing import Any

import pandas as pd

from proteomics_mcp.adapters.modifications import normalize_peptidoform_to_unimod_style
from proteomics_mcp.adapters.table_io import describe_columns, read_candidate_table
from proteomics_mcp.core.spectrum_loader import parse_scan_number
from proteomics_mcp.models.schemas import CandidatePSM

logger = logging.getLogger("proteomics_mcp")


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


def _probability_from_log_probability(value: float | None) -> float | None:
    """Sequence probability from an InstaNovo log-probability.

    InstaNovo scores are natural-log probabilities, so a value just above 0 (numerical noise in InstaNovo+ output, e.g.
    1.9e-08) means a probability of about 1, not 1.9e-08.
    """
    if value is None:
        return None
    return math.exp(min(value, 0.0))


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
    # Positions are only meaningful when there is one value per token (not so for some degenerate beam predictions).
    aligned = len(tokens) == len(confidences)
    return {
        "token_log_probabilities": token_log_probs,
        "token_log_probability_mean": sum(token_log_probs) / len(token_log_probs),
        "token_log_probability_min": min(token_log_probs),
        "token_confidences": confidences,
        "token_confidence_mean": sum(confidences) / len(confidences),
        "token_confidence_min": confidences[min_idx],
        "token_confidence_min_position": min_idx + 1 if aligned else None,
        "token_confidence_min_token": tokens[min_idx] if aligned else None,
    }


# Prediction columns of the InstaNovo output formats, in order of preference, each with its log-probability, tokenised
# and token log-probability columns (first present one wins):
#   final_prediction         InstaNovo+ 1.1.x refined from transformer predictions; the model InstaNovo selected
#   predictions              1.1.0 transformer (with log_probabilities), 1.2.x (with log_probs; refined in combined runs)
#   preds                    1.0.0 and 1.1.x transformer
#   diffusion_predictions    InstaNovo+ 1.1.x; empty where refinement produced nothing
#   transformer_predictions  the transformer prediction that InstaNovo+ 1.1.x refined
#   instanovo_predictions    the transformer prediction in 1.2.x combined runs
#   prediction               minimal exports
_PREDICTION_SOURCES: tuple[tuple[str, tuple[str, ...], str, tuple[str, ...]], ...] = (
    ("final_prediction", ("final_log_probabilities",), "final_prediction_tokenised", ()),
    ("predictions", ("log_probs", "log_probabilities"), "predictions_tokenised", ("token_log_probs", "token_log_probabilities")),
    ("preds", ("log_probs",), "preds_tokenised", ("token_log_probs",)),
    ("diffusion_predictions", ("diffusion_log_probabilities",), "diffusion_predictions_tokenised", ()),
    (
        "transformer_predictions",
        ("transformer_log_probabilities",),
        "transformer_predictions_tokenised",
        ("transformer_token_log_probabilities",),
    ),
    (
        "instanovo_predictions",
        ("instanovo_prediction_log_probability", "instanovo_log_probabilities"),
        "instanovo_predictions_tokenised",
        ("instanovo_prediction_token_log_probabilities", "instanovo_token_log_probabilities"),
    ),
    ("prediction", ("log_probability",), "prediction_tokenised", ("token_log_probability", "token_log_probabilities")),
)


def _align_token_log_probs(
    token_log_probs: list[float] | None, n_tokens: int, *, legacy_beam_layout: bool, reversed_order: bool
) -> list[float] | None:
    """Token log-probabilities in sequence order, one per token.

    Beam-search files with the legacy layout (preds_beam_k / token_log_probs_k columns, InstaNovo 1.0.x and 1.1.x) store
    one extra, trailing value per sequence (end of sequence), and 1.1.x stores the rest in reverse order. Checked against
    greedy runs of the same spectra: after this alignment, the arrays of identical predictions agree for ~95-98% of
    them (Zenodo 10.5281/zenodo.20756892).
    """
    if not token_log_probs or not legacy_beam_layout or len(token_log_probs) != n_tokens + 1:
        return token_log_probs
    aligned = token_log_probs[:-1]
    return aligned[::-1] if reversed_order else aligned


def _detect_schema(columns: set[str]) -> str:
    if {"preds", "log_probs"}.issubset(columns):
        return "legacy_simple"
    if any(col.startswith(("diffusion_", "instanovo_", "transformer_", "final_")) for col in columns):
        return "documented_combined"
    if "predictions" in columns and any(col.startswith("predictions_beam_") for col in columns):
        return "instanovo_beam"
    if "predictions" in columns and ({"log_probs", "log_probabilities"} & columns):
        return "instanovo_predictions"
    if "prediction" in columns and "log_probability" in columns:
        return "instanovo_predictions"
    raise ValueError("Could not detect an InstaNovo prediction schema from table columns.")


def _default_prediction_column(row: pd.Series, schema: str) -> tuple[str | None, str | None, Any, str | None, str | None]:
    """The row's prediction: (peptide, score column, tokenised value, token log-probability column, prediction column)."""
    for prediction_column, score_columns, tokenised_column, token_columns in _PREDICTION_SOURCES:
        peptide = _get_str(row, prediction_column)
        if not peptide:
            continue
        score_column = next((column for column in score_columns if column in row.index), None)
        token_column = next((column for column in token_columns if column in row.index), None)
        return peptide, score_column, row.get(tokenised_column), token_column, prediction_column
    return None, None, None, None, None


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
        f"instanovo_predictions_log_probability_beam_{beam_index}",
        f"log_probs_beam_{beam_index}",
        f"log_probability_beam_{beam_index}",
        f"beam_{beam_index}_log_probability",
        # InstaNovo+ 1.1.x with verbose prediction columns: one column per diffusion prediction
        f"log_probability_{beam_index}",
    )
    return _first_float(row, candidates)


def _beam_prediction(row: pd.Series, beam_index: int) -> str | None:
    for column in (f"predictions_beam_{beam_index}", f"instanovo_predictions_beam_{beam_index}", f"preds_beam_{beam_index}", f"prediction_{beam_index}"):
        value = _get_str(row, column)
        if value:
            return value
    return None


def _next_beam_log_probability(row: pd.Series, peptide: str, raw_score: float | None, schema: str) -> float | None:
    """Log-probability of the runner-up beam, when the beams rank alternatives to the chosen prediction."""
    explicit = _first_float(row, ("next_beam", "nextbeam", "nextscore", "next_score", "next_log_probability", "next_log_prob", "next_best_score"))
    if explicit is not None:
        return explicit
    beam0_peptide = _beam_prediction(row, 0)
    if beam0_peptide is not None:
        # In combined runs the chosen prediction can be the refined one, and the beams the transformer's alternatives.
        return _beam_log_probability(row, 1) if beam0_peptide == peptide else None
    if schema == "instanovo_beam":
        beam0 = _beam_log_probability(row, 0)
        if raw_score is not None and beam0 is not None and not math.isclose(raw_score, beam0):
            return beam0
    return None


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
    frame = read_candidate_table(path, limit=limit)
    try:
        schema = _detect_schema(set(frame.columns))
    except ValueError as exc:
        raise ValueError(f"{exc} Columns found: {describe_columns(frame)}.") from exc
    candidates: list[CandidatePSM] = []
    skipped_without_prediction = 0
    # For MGF input (and mzML input in 1.0.0), InstaNovo writes the spectrum's 0-based position in the file as
    # scan_number, not a scan number. Scan numbers start at 1, so a 0 in the column marks positions. A subset of such
    # a file that lacks the first spectrum cannot be recognised this way.
    legacy_beam_layout = any(column.startswith(("preds_beam_", "token_log_probs_")) for column in frame.columns)
    # InstaNovo 1.1.x adds delta_mass_ppm; 1.0.x, which keeps beam token values in sequence order, does not.
    reversed_beam_tokens = legacy_beam_layout and "delta_mass_ppm" in frame.columns
    scan_values = pd.to_numeric(frame["scan_number"], errors="coerce") if "scan_number" in frame.columns else None
    scan_number_is_position = scan_values is not None and bool((scan_values == 0).any())
    for row_idx, row in frame.iterrows():
        raw_peptide, score_col, tokenised, token_col, prediction_col = _default_prediction_column(row, schema)
        if not raw_peptide:
            skipped_without_prediction += 1
            continue
        peptidoform = normalize_peptidoform_to_instanovo_style(raw_peptide)
        raw_score = _as_float(row.get(score_col)) if score_col else None
        next_beam_score = _next_beam_log_probability(row, raw_peptide, raw_score, schema)
        token_log_probs = _parse_float_list(row.get(token_col)) if token_col else None
        token_log_probs = _align_token_log_probs(
            token_log_probs,
            len(_split_tokens(peptidoform, tokenised)),
            legacy_beam_layout=legacy_beam_layout,
            reversed_order=reversed_beam_tokens,
        )
        token_data = _token_summary(peptidoform, token_log_probs, tokenised)
        confidence = _probability_from_log_probability(raw_score)
        next_beam = _probability_from_log_probability(next_beam_score)
        next_beam_margin = (confidence - next_beam) if confidence is not None and next_beam is not None else None
        scan_number = _scan_number(row)
        spectrum_index = None
        if scan_number_is_position:
            spectrum_index, scan_number = scan_number, None
        candidate_id = _candidate_id(row, path, row_idx)
        candidates.append(
            CandidatePSM(
                candidate_id=candidate_id,
                source_engine="instanovo",
                source_engine_version=None,
                experiment_name=_experiment_name(row, default_experiment_name),
                scan_number=scan_number,
                spectrum_id=_spectrum_id(row),
                spectrum_index=spectrum_index,
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
                source_columns={"schema": schema, "prediction_column": prediction_col, **{str(k): row[k] for k in row.index}},
                **token_data,
            )
        )
    if skipped_without_prediction:
        # InstaNovo writes an empty prediction where decoding found no valid sequence (e.g. every beam failed).
        logger.info("%s: %d of %d rows have no prediction and were skipped", path, skipped_without_prediction, len(frame))
    return candidates
