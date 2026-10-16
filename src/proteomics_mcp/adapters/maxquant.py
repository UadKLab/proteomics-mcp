from __future__ import annotations

import logging
import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

import pandas as pd

from proteomics_mcp.adapters.instanovo import _as_float, _as_int, normalize_peptidoform_to_instanovo_style
from proteomics_mcp.adapters.modifications import build_peptidoform, matching_parenthesis
from proteomics_mcp.core.annotation import peptidoform_neutral_mass
from proteomics_mcp.models.schemas import CandidatePSM

logger = logging.getLogger("proteomics_mcp")

# MaxQuant writes "Mass" with four decimals; a missing or wrong modification is off by far more (e.g. 57.02 Da).
MASS_TOLERANCE_DA = 0.01

NO_FIXED_MODIFICATIONS_WARNING = (
    "MaxQuant fixed modifications unknown: no parameters.txt or mqpar.xml next to the msms table and no "
    "fixed_modifications given. \"Modified sequence\" lists only variable modifications, so fixed ones such as "
    "carbamidomethylation of cysteine are missing from the peptidoform."
)


def _is_missing(value: Any) -> bool:
    return value is None or (isinstance(value, float) and math.isnan(value))


def _clean_text(value: Any) -> str | None:
    if _is_missing(value):
        return None
    text = str(value).strip()
    return text or None


def parse_fixed_modifications(value: str | list[str] | tuple[str, ...] | None) -> list[str]:
    """Fixed modifications in MaxQuant notation, e.g. "Carbamidomethyl (C)", from a list or a ";"-separated string."""
    if value is None:
        return []
    items = value.split(";") if isinstance(value, str) else list(value)
    return [str(item).strip() for item in items if str(item).strip()]


def find_maxquant_parameters(msms_path: str | Path) -> Path | None:
    """The parameters.txt or mqpar.xml of a MaxQuant run, next to its msms table.

    Prefixed copies count too, as in renamed deposits: "maxquant_msms.txt" goes with "maxquant_parameters.txt".
    """
    path = Path(msms_path)
    stem = path.stem
    prefix = stem[: -len("msms")] if stem.lower().endswith("msms") else ""
    names = [f"{prefix}parameters.txt", "parameters.txt", f"{prefix}mqpar.xml", "mqpar.xml"]
    for name in dict.fromkeys(names):
        candidate = path.with_name(name)
        if candidate.is_file():
            return candidate
    return None


def read_fixed_modifications(parameters_path: str | Path) -> list[str] | None:
    """Fixed modifications from a MaxQuant parameters.txt ("Fixed modifications" row) or mqpar.xml.

    Returns None when the file does not state them, and an empty list when it states that there are none.
    """
    path = Path(parameters_path)
    if path.suffix.lower() == ".xml":
        root = ET.parse(path).getroot()
        groups = root.findall(".//fixedModifications")
        if not groups:
            return None
        names = [(item.text or "").strip() for group in groups for item in group.findall("string")]
        return list(dict.fromkeys(name for name in names if name))
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        key, _, value = line.partition("\t")
        if key.strip().lower() == "fixed modifications":
            return parse_fixed_modifications(value)
    return None


def _fixed_modification_sites(modification: str) -> tuple[str, str | None, bool]:
    """(residue letters, terminus "N"/"C"/None, protein-terminal) from a MaxQuant name such as "Carbamidomethyl (C)"."""
    match = re.search(r"\(([^()]*)\)\s*$", modification)
    sites = match.group(1) if match else ""
    lowered = sites.lower()
    terminus = "N" if "n-term" in lowered else "C" if "c-term" in lowered else None
    words = re.sub(r"(?i)\b(any|protein)\b|[nc]-term", " ", sites).split()
    residues = "".join(word for word in words if word.isalpha() and word.isupper())
    return residues, terminus, "protein" in lowered


def _apply_fixed_modifications(
    residues: list[str],
    residue_mods: dict[int, list[str]],
    nterm_mods: list[str],
    cterm_mods: list[str],
    fixed_modifications: list[str],
) -> list[str]:
    """Add fixed modifications to every matching site that carries no variable modification; return warnings."""
    warnings: list[str] = []
    for modification in fixed_modifications:
        letters, terminus, protein_terminal = _fixed_modification_sites(modification)
        if protein_terminal:
            warnings.append(f"Fixed modification {modification!r} applies to protein termini, which msms.txt does not show; not applied.")
            continue
        if terminus == "N":
            if not nterm_mods and (not letters or residues[0] in letters):
                nterm_mods.append(modification)
            continue
        if terminus == "C":
            if not cterm_mods and (not letters or residues[-1] in letters):
                cterm_mods.append(modification)
            continue
        if not letters:
            warnings.append(f"Fixed modification {modification!r} names no residue; not applied.")
            continue
        for position, residue in enumerate(residues, start=1):
            if residue in letters and position not in residue_mods:
                residue_mods[position] = [modification]
    return warnings


def _peptidoform_from_modified_sequence(
    modified_sequence: Any, fallback_sequence: Any = None, fixed_modifications: list[str] | None = None
) -> str | None:
    """Peptidoform from a MaxQuant "Modified sequence", e.g. "_(ac)PEPM(ox)K_" or "_(Acetyl (Protein N-term))PEPM(Oxidation (M))K_".

    A modification before the first residue is N-terminal; one whose name mentions "C-term" is C-terminal. Names may
    contain nested parentheses. Text that is not understood stays in the peptidoform as an explicit token.
    "Modified sequence" lists only variable modifications; fixed_modifications (MaxQuant names, e.g.
    "Carbamidomethyl (C)") are added to every matching site that has no variable modification (#101).
    """
    return _peptidoform_and_warnings(modified_sequence, fallback_sequence, fixed_modifications)[0]


def _peptidoform_and_warnings(
    modified_sequence: Any, fallback_sequence: Any = None, fixed_modifications: list[str] | None = None
) -> tuple[str | None, list[str]]:
    text = _clean_text(modified_sequence)
    if text:
        text = text.replace("_", "")
        residues: list[str] = []
        residue_mods: dict[int, list[str]] = {}
        nterm_mods: list[str] = []
        cterm_mods: list[str] = []
        unlocalized_mods: list[str] = []
        i = 0
        while i < len(text):
            ch = text[i]
            if ch == "(":
                end = matching_parenthesis(text, i)
                if end is None:
                    unlocalized_mods.append(text[i:])
                    break
                name = text[i + 1 : end].strip()
                if not residues:
                    nterm_mods.append(name)
                elif "c-term" in name.lower():
                    cterm_mods.append(name)
                else:
                    residue_mods.setdefault(len(residues), []).append(name)
                i = end + 1
                continue
            if ch.isalpha():
                residues.append(ch.upper())
            i += 1
        if residues:
            warnings = _apply_fixed_modifications(residues, residue_mods, nterm_mods, cterm_mods, fixed_modifications or [])
            peptidoform = build_peptidoform(
                "".join(residues),
                residue_mods=residue_mods,
                nterm_mods=nterm_mods,
                cterm_mods=cterm_mods,
                unlocalized_mods=unlocalized_mods,
            )
            return peptidoform, warnings
    fallback = _clean_text(fallback_sequence)
    if not fallback:
        return None, []
    if fixed_modifications:
        return _peptidoform_and_warnings(f"_{fallback}_", None, fixed_modifications)
    return normalize_peptidoform_to_instanovo_style(fallback), []


def _mass_check_warning(peptidoform: str, reported_mass: float | None) -> str | None:
    """Compare the peptidoform's monoisotopic mass with MaxQuant's "Mass" column, an independent reference."""
    if reported_mass is None:
        return None
    computed = peptidoform_neutral_mass(peptidoform)
    if computed is None:
        return f"Mass of {peptidoform} cannot be computed to compare with MaxQuant's Mass ({reported_mass:.4f} Da)."
    difference = computed - reported_mass
    if abs(difference) > MASS_TOLERANCE_DA:
        return (f"Mass of {peptidoform} is {computed:.4f} Da but MaxQuant reports {reported_mass:.4f} Da ({difference:+.4f} Da): "
                "a modification is missing or wrong.")
    return None


def load_maxquant_psms(
    path: str | Path, *, limit: int | None = None, fixed_modifications: str | list[str] | None = None
) -> list[CandidatePSM]:
    """Read a MaxQuant msms.txt.

    fixed_modifications: MaxQuant names, e.g. "Carbamidomethyl (C)" (a list, or ";"-separated); "" means none. When not
    given, they are read from the parameters.txt or mqpar.xml next to the table. If neither is available, every
    candidate carries a warning instead of a silently guessed default.

    Each peptidoform's mass is checked against MaxQuant's "Mass" column (when present); a difference above
    MASS_TOLERANCE_DA adds a warning to that candidate.
    """
    frame = pd.read_csv(path, sep="\t", nrows=limit)
    required = {"Raw file", "Scan number", "Modified sequence", "Charge", "m/z"}
    missing = sorted(col for col in required if col not in frame.columns)
    if missing:
        raise ValueError(f"Missing required MaxQuant columns: {', '.join(missing)}")

    table_warnings: list[str] = []
    if fixed_modifications is not None:
        fixed = parse_fixed_modifications(fixed_modifications)
    else:
        parameters_path = find_maxquant_parameters(path)
        fixed = read_fixed_modifications(parameters_path) if parameters_path is not None else None
        if fixed is None:
            fixed = []
            table_warnings.append(NO_FIXED_MODIFICATIONS_WARNING)
            logger.warning("%s: %s", path, NO_FIXED_MODIFICATIONS_WARNING)
        else:
            logger.info("%s: fixed modifications from %s: %s", path, parameters_path.name, ", ".join(fixed) or "none")

    candidates: list[CandidatePSM] = []
    mass_mismatches = 0
    for row_idx, row in frame.iterrows():
        peptidoform, row_warnings = _peptidoform_and_warnings(row.get("Modified sequence"), row.get("Sequence"), fixed)
        if not peptidoform:
            continue
        mass_warning = _mass_check_warning(peptidoform, _as_float(row.get("Mass")))
        if mass_warning:
            row_warnings.append(mass_warning)
            mass_mismatches += 1

        experiment_name = _clean_text(row.get("Raw file"))
        scan_number = _as_int(row.get("Scan number"))
        raw_score = _as_float(row.get("Score"))
        pep = _as_float(row.get("PEP"))
        confidence = max(0.0, min(1.0, 1.0 - pep)) if pep is not None else None
        spectrum_id = f"{experiment_name}:{scan_number}" if experiment_name and scan_number is not None else None

        candidates.append(
            CandidatePSM(
                candidate_id=f"{Path(path).stem}:{row_idx}",
                source_engine="maxquant",
                source_engine_version=None,
                experiment_name=experiment_name,
                scan_number=scan_number,
                spectrum_id=spectrum_id,
                usi=None,
                peptidoform=peptidoform,
                precursor_mz=_as_float(row.get("m/z")),
                precursor_charge=_as_int(row.get("Charge")),
                rank=None,
                raw_score=raw_score,
                raw_score_name="Score" if raw_score is not None else None,
                # MaxQuant score deltas are not de novo beam confidences.
                next_beam=None,
                next_beam_margin=None,
                confidence=confidence,
                confidence_type="posterior_error_probability" if pep is not None else None,
                source_columns={str(k): row[k] for k in row.index},
                warnings=[*table_warnings, *row_warnings],
            )
        )
    if mass_mismatches:
        logger.warning("%s: %d of %d peptidoforms differ from MaxQuant's Mass by more than %.2f Da",
                       path, mass_mismatches, len(candidates), MASS_TOLERANCE_DA)
    return candidates
