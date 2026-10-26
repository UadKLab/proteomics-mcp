from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

import pandas as pd

from proteomics_mcp.adapters.table_io import describe_columns, read_candidate_table
from proteomics_mcp.adapters.instanovo import _as_float, _as_int, normalize_peptidoform_to_instanovo_style
from proteomics_mcp.adapters.modifications import (
    MASS_MATCH_TOLERANCE,
    build_peptidoform,
    nearest_builtin_unimod,
    parse_fragpipe_modifications,
    resolve_mod_token,
)
from proteomics_mcp.core.annotation import AA_RESIDUE_MASS, H2O_MASS, HYDROGEN_MASS
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


# FragPipe's "Modified Peptide" writes total masses: of the residue ("M[147]"), of the N-terminal group, i.e. H plus the
# modification ("n[43]" for acetyl), or of the C-terminal group, OH plus the modification. Some tables carry deltas
# instead ("M[15.994915]"), so both readings are tried.
_TERMINAL_BASE_MASS = {"n": HYDROGEN_MASS, "c": H2O_MASS - HYDROGEN_MASS}
# Integer-rounded masses carry up to 0.5 Da of rounding; they are snapped to the nearest built-in modification.
_ROUNDED_MASS_TOLERANCE = 0.5


def _delta_token(value: str, base_mass: float) -> str:
    """Modification token for a bracketed mass that is either a total (base + modification) or a delta."""
    mass = float(value)
    rounded = re.fullmatch(r"\d+", value.strip()) is not None
    # A mass-adding modification makes the total at least the base mass; a smaller value is more likely a delta.
    readings = [mass - base_mass, mass] if mass >= base_mass - _ROUNDED_MASS_TOLERANCE else [mass, mass - base_mass]
    # Common modifications first, for both readings: the full Unimod catalogue has an entry near almost any mass.
    for delta in readings:
        known = nearest_builtin_unimod(delta, _ROUNDED_MASS_TOLERANCE if rounded else MASS_MATCH_TOLERANCE)
        if known:
            return known
    if not rounded:
        for delta in readings:
            token = resolve_mod_token(f"{delta:+.6f}")
            if token and token.startswith("UNIMOD:"):
                return token
    return f"{readings[0]:+.6f}"


def _peptidoform_from_modified_peptide(value: str) -> str | None:
    """Peptidoform from a FragPipe "Modified Peptide" such as "n[43]PEPTM[147]K", whose masses are totals, not deltas."""
    match = re.fullmatch(r"(?:n\[(?P<n>[\d.]+)\])?(?P<body>.*?)(?:c\[(?P<c>[\d.]+)\])?", value.strip())
    if not match:
        return normalize_peptidoform_to_instanovo_style(value)
    residues: list[str] = []
    residue_mods: dict[int, list[str]] = {}
    for residue, total in re.findall(r"([A-Z])(?:\[([\d.]+)\])?", match.group("body")):
        residues.append(residue)
        if total:
            residue_mods[len(residues)] = [_delta_token(total, AA_RESIDUE_MASS.get(residue, 0.0))]
    if "".join(residues) != re.sub(r"\[[\d.]+\]", "", match.group("body")):
        return normalize_peptidoform_to_instanovo_style(value)  # not the notation above; keep every token visible
    nterm = [_delta_token(match.group("n"), _TERMINAL_BASE_MASS["n"])] if match.group("n") else []
    cterm = [_delta_token(match.group("c"), _TERMINAL_BASE_MASS["c"])] if match.group("c") else []
    return build_peptidoform("".join(residues), residue_mods=residue_mods, nterm_mods=nterm, cterm_mods=cterm)


def _peptidoform_from_row(row: pd.Series) -> str | None:
    """Peptidoform from "Assigned Modifications" (deltas, including N-term/C-term), else from "Modified Peptide"."""
    peptide = _clean_text(row.get("Peptide"))
    residue_mods, nterm_mods, cterm_mods, unparsed = parse_fragpipe_modifications(row.get("Assigned Modifications"), peptide)
    modified_peptide = _clean_text(row.get("Modified Peptide"))
    if modified_peptide and modified_peptide.lower() == "nan":
        modified_peptide = None
    if peptide and unparsed and modified_peptide:
        return _peptidoform_from_modified_peptide(modified_peptide)
    if peptide and (residue_mods or nterm_mods or cterm_mods or unparsed):
        return build_peptidoform(
            peptide, residue_mods=residue_mods, nterm_mods=nterm_mods, cterm_mods=cterm_mods, unlocalized_mods=unparsed
        )
    if modified_peptide:
        return _peptidoform_from_modified_peptide(modified_peptide)
    return normalize_peptidoform_to_instanovo_style(peptide) if peptide else None


def load_fragpipe_psms(path: str | Path, *, limit: int | None = None) -> list[CandidatePSM]:
    frame = read_candidate_table(path, limit=limit)
    required = {"Spectrum", "Peptide"}
    missing = sorted(col for col in required if col not in frame.columns)
    if missing:
        raise ValueError(f"Missing required FragPipe columns: {', '.join(missing)}. Columns found: {describe_columns(frame)}.")

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
