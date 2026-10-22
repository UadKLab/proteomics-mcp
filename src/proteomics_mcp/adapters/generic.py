from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import pandas as pd

from proteomics_mcp.adapters.instanovo import _as_float, _as_int, _parse_float_list, _token_summary
from proteomics_mcp.adapters.table_io import read_candidate_table
from proteomics_mcp.models.schemas import CandidatePSM


def _confidence_from_score(value: float | None, score_type: str) -> float | None:
    if value is None:
        return None
    if score_type == "confidence_0_1":
        return max(0.0, min(1.0, value))
    if score_type == "posterior_error_probability":
        return max(0.0, min(1.0, 1.0 - value))
    if score_type in {"evalue", "score_unscaled", "unknown"}:
        return None
    raise ValueError(f"Unsupported score_type: {score_type}")


def _text_or_none(value: Any) -> str | None:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    text = str(value).strip()
    return text or None


def load_generic_candidates(
    path: str | Path,
    *,
    peptide_column: str,
    scan_column: str | None = None,
    score_column: str | None = None,
    next_beam_column: str | None = None,
    next_beam_margin_column: str | None = None,
    next_score_column: str | None = None,
    next_score_margin_column: str | None = None,
    score_type: str = "unknown",
    charge_column: str | None = None,
    precursor_mz_column: str | None = None,
    spectrum_id_column: str | None = None,
    run_column: str | None = None,
    experiment_name_column: str | None = None,
    usi_column: str | None = None,
    token_log_probabilities_column: str | None = None,
    source_engine: str = "generic",
    limit: int | None = None,
) -> list[CandidatePSM]:
    frame = read_candidate_table(path, limit=limit)
    if run_column and experiment_name_column and run_column != experiment_name_column:
        raise ValueError("run_column and experiment_name_column name the same field; give one of them.")
    run_column = run_column or experiment_name_column
    missing = [
        col
        for col in [
            peptide_column,
            scan_column,
            score_column,
            next_beam_column or next_score_column,
            next_beam_margin_column or next_score_margin_column,
            charge_column,
            precursor_mz_column,
            spectrum_id_column,
            run_column,
            usi_column,
            token_log_probabilities_column,
        ]
        if col and col not in frame.columns
    ]
    if missing:
        raise ValueError(f"Missing required mapped columns: {', '.join(missing)}")
    candidates: list[CandidatePSM] = []
    for row_idx, row in frame.iterrows():
        peptide = str(row[peptide_column]).strip()
        if not peptide or peptide.lower() == "nan":
            continue
        raw_score = _as_float(row.get(score_column)) if score_column else None
        next_beam_col = next_beam_column or next_score_column
        next_beam_margin_col = next_beam_margin_column or next_score_margin_column
        next_beam = _as_float(row.get(next_beam_col)) if next_beam_col else None
        next_beam_margin = _as_float(row.get(next_beam_margin_col)) if next_beam_margin_col else None
        if next_beam_margin is None and raw_score is not None and next_beam is not None:
            next_beam_margin = raw_score - next_beam
        token_log_probs = _parse_float_list(row.get(token_log_probabilities_column)) if token_log_probabilities_column else None
        token_data = _token_summary(peptide, token_log_probs)
        candidates.append(
            CandidatePSM(
                candidate_id=f"{Path(path).stem}:{row_idx}",
                source_engine=source_engine,
                peptidoform=peptide,
                scan_number=_as_int(row.get(scan_column)) if scan_column else None,
                spectrum_id=str(row.get(spectrum_id_column)) if spectrum_id_column else None,
                # The PSM's run, so it finds its spectrum in a multi-run table (#142)
                experiment_name=_text_or_none(row.get(run_column)) if run_column else None,
                usi=str(row.get(usi_column)) if usi_column else None,
                precursor_mz=_as_float(row.get(precursor_mz_column)) if precursor_mz_column else None,
                precursor_charge=_as_int(row.get(charge_column)) if charge_column else None,
                rank=0,
                raw_score=raw_score,
                raw_score_name=score_column,
                next_beam=next_beam,
                next_beam_margin=next_beam_margin,
                confidence=_confidence_from_score(raw_score, score_type),
                confidence_type=score_type,
                source_columns={str(k): row[k] for k in row.index},
                **token_data,
            )
        )
    return candidates
