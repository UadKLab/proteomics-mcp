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


# PSM and peptide tables write "<count>x<name> [<sites>]", e.g. "1xTMT6plex [N-Term]; 2xCarbamidomethyl [C5; C12]";
# sites may carry a localisation probability ("S5(99.2)") or be unlocalised ("S/T").
_COUNTED_MOD_RE = re.compile(r"(?P<count>\d+)x(?P<name>[^\[;]+?)\s*\[(?P<sites>[^\]]*)\]")


def _counted_modifications(
    mods_text: str, residues: list[str]
) -> tuple[dict[int, list[str]], list[str], list[str], list[str], str]:
    """Parse the "NxName [sites]" notation; return (residue, N-term, C-term, unlocalised modifications, remaining text)."""
    residue_mods: dict[int, list[str]] = {}
    nterm: list[str] = []
    cterm: list[str] = []
    unlocalized: list[str] = []
    for match in _COUNTED_MOD_RE.finditer(mods_text):
        name = match.group("name").strip()
        placed = 0
        for site in [s.strip() for s in match.group("sites").split(";") if s.strip()]:
            site = re.sub(r"\([\d.]+\)$", "", site).strip()
            residue_site = re.fullmatch(r"(?P<residue>[A-Za-z])(?P<position>\d+)", site)
            position = int(residue_site.group("position")) if residue_site else 0
            if residue_site and 1 <= position <= len(residues) and residues[position - 1] == residue_site.group("residue").upper():
                residue_mods.setdefault(int(residue_site.group("position")), []).append(name)
                placed += 1
            elif site.lower().startswith("n-term"):
                nterm.append(name)
                placed += 1
            elif site.lower().startswith("c-term"):
                cterm.append(name)
                placed += 1
        unlocalized.extend([name] * max(0, int(match.group("count")) - placed))
    return residue_mods, nterm, cterm, unlocalized, _COUNTED_MOD_RE.sub("", mods_text)


def _peptidoform_from_annotated_sequence(annotated_sequence: Any, modifications: Any = None) -> str | None:
    core = _peptide_core(annotated_sequence)
    if not core:
        return None
    mods_text = _clean_text(modifications) or ""
    met_removed = "met-loss" in mods_text.lower() and core and core[0].lower() == "m"
    if met_removed:
        core = core[1:]
    residues = [ch.upper() for ch in core if ch.isalpha()]
    if not residues:
        return None

    residue_mods, nterm_mods, cterm_mods, unlocalized_mods, mods_text = _counted_modifications(mods_text, residues)
    if met_removed:
        # The methionine is already removed from the sequence; keep only what follows the loss ("Met-loss+Acetyl" -> "Acetyl").
        nterm_mods = [re.sub(r"^met-loss\+?", "", name, flags=re.IGNORECASE) for name in nterm_mods]
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
            if mod_name and met_removed:
                mod_name = re.sub(r"^met-loss\+?", "", mod_name, flags=re.IGNORECASE)
            if mod_name and mod_name.lower() != "met-loss":
                nterm_mods.append(mod_name)
            continue
        if part.startswith("C-Term"):
            mod_name = _terminal_mod_name(part)
            if mod_name:
                cterm_mods.append(mod_name)
            continue
        unlocalized_mods.append(part)

    return build_peptidoform(
        "".join(residues),
        residue_mods=residue_mods,
        nterm_mods=[name for name in nterm_mods if name],
        cterm_mods=cterm_mods,
        unlocalized_mods=unlocalized_mods,
    )


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
