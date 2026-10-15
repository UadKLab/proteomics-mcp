from __future__ import annotations

import math
import re
from pathlib import Path
from typing import Any

import pandas as pd

from proteomics_mcp.adapters.instanovo import _as_float, _as_int, normalize_peptidoform_to_instanovo_style
from proteomics_mcp.adapters.modifications import build_peptidoform, matching_parenthesis
from proteomics_mcp.models.schemas import CandidatePSM


def _is_missing(value: Any) -> bool:
    return value is None or (isinstance(value, float) and math.isnan(value))


def _clean_text(value: Any) -> str | None:
    if _is_missing(value):
        return None
    text = str(value).strip()
    return text or None


def _peptidoform_from_modified_sequence(modified_sequence: Any, fallback_sequence: Any = None) -> str | None:
    """Peptidoform from a MaxQuant "Modified sequence", e.g. "_(ac)PEPM(ox)K_" or "_(Acetyl (Protein N-term))PEPM(Oxidation (M))K_".

    A modification before the first residue is N-terminal; one whose name mentions "C-term" is C-terminal. Names may
    contain nested parentheses. Text that is not understood stays in the peptidoform as an explicit token.
    """
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
            return build_peptidoform(
                "".join(residues),
                residue_mods=residue_mods,
                nterm_mods=nterm_mods,
                cterm_mods=cterm_mods,
                unlocalized_mods=unlocalized_mods,
            )
    fallback = _clean_text(fallback_sequence)
    return normalize_peptidoform_to_instanovo_style(fallback) if fallback else None


def load_maxquant_psms(path: str | Path, *, limit: int | None = None) -> list[CandidatePSM]:
    frame = pd.read_csv(path, sep="\t", nrows=limit)
    required = {"Raw file", "Scan number", "Modified sequence", "Charge", "m/z"}
    missing = sorted(col for col in required if col not in frame.columns)
    if missing:
        raise ValueError(f"Missing required MaxQuant columns: {', '.join(missing)}")

    candidates: list[CandidatePSM] = []
    for row_idx, row in frame.iterrows():
        peptidoform = _peptidoform_from_modified_sequence(row.get("Modified sequence"), row.get("Sequence"))
        if not peptidoform:
            continue

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
            )
        )
    return candidates
