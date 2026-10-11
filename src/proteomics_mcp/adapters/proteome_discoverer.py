from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

from proteomics_mcp.adapters.instanovo import _as_float, _as_int, normalize_peptidoform_to_instanovo_style
from proteomics_mcp.adapters.modifications import build_peptidoform
from proteomics_mcp.adapters.table_io import read_candidate_table
from proteomics_mcp.models.schemas import CandidatePSM


def _is_missing(value: Any) -> bool:
    return value is None or (isinstance(value, float) and math.isnan(value))


def _clean_text(value: Any) -> str | None:
    if _is_missing(value):
        return None
    text = str(value).strip()
    return text or None


def _run_name(value: Any) -> str | None:
    text = _clean_text(value)
    if not text:
        return None
    path_name = Path(text)
    stem = path_name.stem
    return stem or text


def _peptide_core(value: Any) -> str | None:
    text = _clean_text(value)
    if not text:
        return None
    parts = text.split(".")
    return parts[1] if len(parts) >= 3 else text


def _terminal_mod_name(mod_text: str) -> str | None:
    names = re.findall(r"\(([^()]+)\)", mod_text)
    for name in reversed(names):
        clean = str(name).strip()
        if clean.lower() not in {"prot", "pep", "protein"}:
            return clean
    return None


def _peptidoform_from_annotated_sequence(annotated_sequence: Any, modifications: Any = None) -> str | None:
    core = _peptide_core(annotated_sequence)
    if not core:
        return None
    mods_text = _clean_text(modifications) or ""
    if "met-loss" in mods_text.lower() and core and core[0].lower() == "m":
        core = core[1:]
    residues = [ch.upper() for ch in core if ch.isalpha()]
    if not residues:
        return None

    residue_mods: dict[int, list[str]] = {}
    nterm_mods: list[str] = []
    cterm_mods: list[str] = []
    for part in [segment.strip() for segment in mods_text.split(";") if segment.strip()]:
        residue_match = re.match(r"^(?P<residue>[A-Za-z])(?P<position>\d+)\((?P<mod>.+)\)$", part)
        if residue_match:
            position = int(residue_match.group("position"))
            mod_name = residue_match.group("mod").strip()
            if 1 <= position <= len(residues):
                residue_mods.setdefault(position, []).append(mod_name)
            continue
        if part.startswith("N-Term"):
            mod_name = _terminal_mod_name(part)
            if mod_name and mod_name.lower() != "met-loss":
                nterm_mods.append(mod_name)
            continue
        if part.startswith("C-Term"):
            mod_name = _terminal_mod_name(part)
            if mod_name:
                cterm_mods.append(mod_name)

    return build_peptidoform("".join(residues), residue_mods=residue_mods, nterm_mods=nterm_mods, cterm_mods=cterm_mods)


def load_proteome_discoverer_psms(path: str | Path, *, limit: int | None = None, sheet_name: str | int = 0) -> list[CandidatePSM]:
    frame = read_candidate_table(path, limit=limit, sheet_name=sheet_name)
    required = {"Annotated Sequence", "Spectrum File", "First Scan", "Charge", "m/z [Da]"}
    missing = sorted(col for col in required if col not in frame.columns)
    if missing:
        raise ValueError(f"Missing required Proteome Discoverer columns: {', '.join(missing)}")

    candidates: list[CandidatePSM] = []
    for row_idx, row in frame.iterrows():
        peptidoform = _peptidoform_from_annotated_sequence(row.get("Annotated Sequence"), row.get("Modifications"))
        if not peptidoform:
            continue
        experiment_name = _run_name(row.get("Spectrum File"))
        scan_number = _as_int(row.get("First Scan"))
        charge = _as_int(row.get("Charge"))
        precursor_mz = _as_float(row.get("m/z [Da]"))
        raw_score = _as_float(row.get("XCorr"))
        pep = _as_float(row.get("Percolator PEP"))
        confidence = max(0.0, min(1.0, 1.0 - pep)) if pep is not None else None

        spectrum_id = None
        if experiment_name and scan_number is not None:
            spectrum_id = f"{experiment_name}:{scan_number}"

        candidates.append(
            CandidatePSM(
                candidate_id=f"{Path(path).stem}:{row_idx}",
                source_engine="proteome_discoverer",
                source_engine_version=None,
                experiment_name=experiment_name,
                scan_number=scan_number,
                spectrum_id=spectrum_id,
                usi=None,
                peptidoform=peptidoform,
                precursor_mz=precursor_mz,
                precursor_charge=charge,
                rank=_as_int(row.get("Rank")),
                raw_score=raw_score,
                raw_score_name="XCorr" if raw_score is not None else None,
                # DeltaScore is a search-engine score difference, not a de novo beam confidence.
                next_beam=None,
                next_beam_margin=None,
                confidence=confidence,
                confidence_type="posterior_error_probability" if pep is not None else None,
                source_columns={str(k): row[k] for k in row.index},
            )
        )
    return candidates
