from __future__ import annotations

from pathlib import Path

import pandas as pd

from proteomics_mcp.adapters.instanovo import _as_float, _as_int, _is_missing, _probability_from_score, normalize_peptidoform_to_instanovo_style
from proteomics_mcp.adapters.table_io import describe_columns, read_candidate_table
from proteomics_mcp.models.schemas import CandidatePSM


def _get_str(row: pd.Series, column: str) -> str | None:
    if column not in row.index or _is_missing(row[column]):
        return None
    value = str(row[column]).strip()
    return value or None


def load_casanovo_predictions(path: str | Path, *, limit: int | None = None) -> list[CandidatePSM]:
    """Load Casanovo prediction exports into canonical candidates."""
    frame = read_candidate_table(path, limit=limit)
    if not {"predictions", "prediction"} & set(frame.columns):
        # Without it every row is skipped and the table would load as empty, silently (#140).
        raise ValueError(f"Missing required Casanovo column 'predictions'. Columns found: {describe_columns(frame)}.")
    candidates: list[CandidatePSM] = []
    for row_idx, row in frame.iterrows():
        peptidoform = _get_str(row, "predictions") or _get_str(row, "prediction")
        if not peptidoform:
            continue
        peptidoform = normalize_peptidoform_to_instanovo_style(peptidoform)
        raw_score = _as_float(row.get("log_probs"))
        confidence = _probability_from_score(raw_score)
        candidate_id = str(row.get("prediction_id")) if "prediction_id" in row.index and not _is_missing(row.get("prediction_id")) else f"{Path(path).stem}:{row_idx}"
        candidates.append(
            CandidatePSM(
                candidate_id=candidate_id,
                source_engine="casanovo",
                source_engine_version=None,
                experiment_name=_get_str(row, "group") or _get_str(row, "experiment_name"),
                scan_number=_as_int(row.get("scan_number")),
                spectrum_id=_get_str(row, "spectrum_id"),
                usi=_get_str(row, "usi"),
                peptidoform=peptidoform,
                precursor_mz=_as_float(row.get("precursor_mz")),
                precursor_charge=_as_int(row.get("precursor_charge")),
                rank=0,
                raw_score=raw_score,
                raw_score_name="log_probs" if raw_score is not None else None,
                next_beam=None,
                next_beam_margin=None,
                confidence=confidence,
                confidence_type="sequence_probability_from_log_probability" if confidence is not None else None,
                source_columns={"schema": "casanovo_predictions", **{str(k): row[k] for k in row.index}},
            )
        )
    return candidates
