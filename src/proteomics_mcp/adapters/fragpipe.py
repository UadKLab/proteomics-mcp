from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

import pandas as pd

from proteomics_mcp.adapters.instanovo import _as_float, _as_int, normalize_peptidoform_to_instanovo_style
from proteomics_mcp.adapters.modifications import build_peptidoform, parse_fragpipe_assigned_modifications
from proteomics_mcp.models.schemas import CandidatePSM


def _is_missing(value: Any) -> bool:
    return value is None or (isinstance(value, float) and math.isnan(value))


def _clean_text(value: Any) -> str | None:
    if _is_missing(value):
        return None
    text = str(value).strip()
    return text or None


def _spectrum_parts(value: Any) -> tuple[str | None, int | None, int | None]:
    text = _clean_text(value)
    if not text:
        return None, None, None
    match = re.match(r"^(?P<run>.+?)\.(?P<scan>\d+)\.(?P<scan_repeat>\d+)\.(?P<charge>\d+)$", text)
    if match:
        return match.group("run"), int(match.group("scan")), int(match.group("charge"))
    fallback = re.search(r"^(?P<run>.+?)[.:](?P<scan>\d+)$", text)
    if fallback:
        return fallback.group("run"), int(fallback.group("scan")), None
    return text, None, None


def _peptidoform_from_row(row: pd.Series) -> str | None:
    peptide = _clean_text(row.get("Peptide"))
    assigned_mods = parse_fragpipe_assigned_modifications(row.get("Assigned Modifications"), peptide)
    if peptide and assigned_mods:
        return build_peptidoform(peptide, residue_mods=assigned_mods)
    for column in ("Modified Peptide", "Peptide"):
        value = _clean_text(row.get(column))
        if value and value.lower() != "nan":
            return normalize_peptidoform_to_instanovo_style(value)
    return None


def load_fragpipe_psms(path: str | Path, *, limit: int | None = None) -> list[CandidatePSM]:
    frame = pd.read_csv(path, sep="\t", nrows=limit)
    required = {"Spectrum", "Peptide"}
    missing = sorted(col for col in required if col not in frame.columns)
    if missing:
        raise ValueError(f"Missing required FragPipe columns: {', '.join(missing)}")

    candidates: list[CandidatePSM] = []
    for row_idx, row in frame.iterrows():
        peptidoform = _peptidoform_from_row(row)
        if not peptidoform:
            continue
        run_name, scan_number, spectrum_charge = _spectrum_parts(row.get("Spectrum"))
        precursor_charge = _as_int(row.get("Charge")) or spectrum_charge
        precursor_mz = _as_float(row.get("Calibrated Observed M/Z"))
        if precursor_mz is None:
            precursor_mz = _as_float(row.get("Observed M/Z"))
        raw_score = _as_float(row.get("Hyperscore"))
        probability = _as_float(row.get("Probability"))

        candidates.append(
            CandidatePSM(
                candidate_id=f"{Path(path).stem}:{row_idx}",
                source_engine="fragpipe",
                source_engine_version=None,
                experiment_name=run_name,
                scan_number=scan_number,
                spectrum_id=_clean_text(row.get("Spectrum")),
                usi=None,
                peptidoform=peptidoform,
                precursor_mz=precursor_mz,
                precursor_charge=precursor_charge,
                rank=0,
                raw_score=raw_score,
                raw_score_name="Hyperscore" if raw_score is not None else None,
                # FragPipe Nextscore is preserved upstream but is not a de novo beam confidence.
                next_beam=None,
                next_beam_margin=None,
                confidence=max(0.0, min(1.0, probability)) if probability is not None else None,
                confidence_type="confidence_0_1" if probability is not None else None,
                source_columns={str(k): row[k] for k in row.index},
            )
        )
    return candidates
