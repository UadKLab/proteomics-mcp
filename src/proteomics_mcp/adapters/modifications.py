from __future__ import annotations

import math
import re
from functools import lru_cache
from typing import Any

NAME_TO_UNIMOD = {
    "acetyl": "UNIMOD:1",
    "carbamidomethyl": "UNIMOD:4",
    "carbamyl": "UNIMOD:5",
    "deamidated": "UNIMOD:7",
    "phospho": "UNIMOD:21",
    "gln->pyro-glu": "UNIMOD:28",
    "glu->pyro-glu": "UNIMOD:27",
    "pyro-glu": "UNIMOD:28",
    "oxidation": "UNIMOD:35",
    "methyl": "UNIMOD:34",
    "formyl": "UNIMOD:122",
    "ox": "UNIMOD:35",
    "ph": "UNIMOD:21",
    "de": "UNIMOD:7",
    "ac": "UNIMOD:1",
    "cam": "UNIMOD:4",
    "cm": "UNIMOD:4",
    "pgq": "UNIMOD:28",
    "pge": "UNIMOD:27",
}

UNIMOD_TO_MASS = {
    "UNIMOD:1": 42.01056468403,
    "UNIMOD:4": 57.021463734999996,
    "UNIMOD:5": 43.005813,
    "UNIMOD:7": 0.9840155848,
    "UNIMOD:21": 79.9663308946,
    "UNIMOD:27": -18.01056468403,
    "UNIMOD:28": -17.02654910112,
    "UNIMOD:35": 15.99491461957,
    "UNIMOD:34": 14.01565006446,
    "UNIMOD:122": 27.99491461957,
}

UNIMOD_TO_DISPLAY_NAME = {
    "UNIMOD:1": "acetyl",
    "UNIMOD:4": "cam",
    "UNIMOD:5": "carbamyl",
    "UNIMOD:7": "deamidated",
    "UNIMOD:21": "ph",
    "UNIMOD:27": "pyroglu",
    "UNIMOD:28": "pyroglu",
    "UNIMOD:35": "ox",
    "UNIMOD:34": "methyl",
    "UNIMOD:122": "formyl",
}

MASS_MATCH_TOLERANCE = 0.01
_UNIMOD_CATALOG: Any | None = None


def _load_unimod_catalog() -> Any | None:
    global _UNIMOD_CATALOG
    if _UNIMOD_CATALOG is not None:
        return _UNIMOD_CATALOG
    try:
        from pyteomics.mass import Unimod  # type: ignore
    except ImportError:
        return None
    try:
        _UNIMOD_CATALOG = Unimod()
    except Exception:
        _UNIMOD_CATALOG = None
    return _UNIMOD_CATALOG


def _find_unimod_entry(token: str | None) -> Any | None:
    if token is None:
        return None
    catalog = _load_unimod_catalog()
    if catalog is None:
        return None
    clean = str(token).strip()
    if not clean:
        return None
    if _parse_mass_token(clean) is not None:
        return None
    if clean.upper().startswith("UNIMOD:"):
        record_id = clean.split(":", 1)[1].strip()
        try:
            return next((entry for entry in catalog.mods if str(entry.get("record_id")) == record_id), None)
        except Exception:
            return None
    try:
        mod = catalog.by_title(clean)
    except Exception:
        mod = None
    if mod is None:
        try:
            mod = catalog.by_name(clean)
        except Exception:
            mod = None
    if isinstance(mod, list):
        mod = next(
            (
                entry
                for entry in mod
                if str(entry.get("title", "")).lower() == clean.lower()
                or str(entry.get("full_name", "")).lower() == clean.lower()
            ),
            mod[0] if mod else None,
        )
    return mod


def resolve_unimod_id(token: str | None) -> str | None:
    if token is None:
        return None
    clean = str(token).strip()
    if not clean:
        return None
    if clean.upper().startswith("UNIMOD:"):
        suffix = clean.split(":", 1)[1].strip()
        return f"UNIMOD:{suffix}" if suffix else None
    mapped = NAME_TO_UNIMOD.get(clean.lower())
    if mapped:
        return mapped
    mod = _find_unimod_entry(clean)
    if mod is not None:
        return f"UNIMOD:{int(mod['record_id'])}"
    mass = _parse_mass_token(clean)
    if mass is None:
        return None
    return _resolve_mass_to_unimod(mass)


@lru_cache(maxsize=1024)
def unimod_mass(token: str | None) -> float | None:
    resolved = resolve_unimod_id(token)
    if resolved is None:
        return None
    if resolved in UNIMOD_TO_MASS:
        return UNIMOD_TO_MASS[resolved]
    mod = _find_unimod_entry(resolved)
    if mod is None:
        return None
    try:
        return float(mod["mono_mass"])
    except Exception:
        return None


def format_peptidoform_for_display(peptidoform: str | None) -> str | None:
    """Render normalized modification tokens in compact InstaNovo-style notation."""
    if peptidoform is None:
        return None
    text = str(peptidoform).strip()
    if not text:
        return None

    def replace_mod(match: Any) -> str:
        token = str(match.group(1)).strip()
        resolved = resolve_unimod_id(token) or token.upper()
        display = UNIMOD_TO_DISPLAY_NAME.get(resolved, resolved)
        return f"({display})"

    return re.sub(r"\[([^\]]+)\]", replace_mod, text)


def resolve_mod_token(token: str | None) -> str | None:
    if token is None:
        return None
    clean = str(token).strip()
    if not clean:
        return None
    resolved = resolve_unimod_id(clean)
    if resolved is not None:
        return resolved
    mass = _parse_mass_token(clean)
    if mass is None:
        return clean
    match = _resolve_mass_to_unimod(mass)
    return match or clean


def normalize_peptidoform_to_unimod_style(peptidoform: str | None) -> str | None:
    if peptidoform is None:
        return None
    text = str(peptidoform).strip()
    if not text:
        return None

    def replace_mod(match: Any) -> str:
        token = resolve_mod_token(str(match.group(1)))
        return f"[{token}]" if token else ""

    text = text.replace(" ", "")
    # Current InstaNovo exports can encode residue modifications as M(ox).
    # Convert only parentheses attached to an amino-acid residue so unrelated
    # free-text punctuation is not interpreted as a modification.
    text = re.sub(
        r"(?<=[A-Z])\(([^()]+)\)",
        lambda match: f"[{resolve_mod_token(str(match.group(1))) or str(match.group(1)).strip()}]",
        text,
    )
    text = re.sub(r"\[([^\]]+)\]", replace_mod, text)
    text = re.sub(r"^(\[[^\]]+\])(?=[A-Z])", r"\1-", text)
    return text


def build_peptidoform(
    sequence: str | None,
    *,
    residue_mods: dict[int, list[str]] | None = None,
    nterm_mods: list[str] | None = None,
    cterm_mods: list[str] | None = None,
) -> str | None:
    if sequence is None:
        return None
    residues = [ch.upper() for ch in str(sequence) if ch.isalpha()]
    if not residues:
        return None
    residue_mods = residue_mods or {}
    nterm_mods = nterm_mods or []
    cterm_mods = cterm_mods or []
    peptide = "".join(
        f"{residue}{''.join(f'[{resolve_mod_token(mod) or mod}]' for mod in residue_mods.get(idx, []))}"
        for idx, residue in enumerate(residues, start=1)
    )
    if nterm_mods:
        peptide = f"{''.join(f'[{resolve_mod_token(mod) or mod}]' for mod in nterm_mods)}-{peptide}"
    if cterm_mods:
        peptide = f"{peptide}-{''.join(f'[{resolve_mod_token(mod) or mod}]' for mod in cterm_mods)}"
    return normalize_peptidoform_to_unimod_style(peptide)


def parse_fragpipe_assigned_modifications(value: Any, sequence: str | None) -> dict[int, list[str]]:
    text = str(value).strip() if value is not None and not (isinstance(value, float) and math.isnan(value)) else ""
    if not text or not sequence:
        return {}
    residues = [ch.upper() for ch in str(sequence) if ch.isalpha()]
    mods: dict[int, list[str]] = {}
    for part in [segment.strip() for segment in text.split(",") if segment.strip()]:
        match = re.match(r"(?P<position>\d+)(?P<residue>[A-Za-z])\((?P<mod>[^()]+)\)$", part)
        if not match:
            continue
        position = int(match.group("position"))
        residue = match.group("residue").upper()
        raw_mod = match.group("mod").strip()
        if 1 <= position <= len(residues) and residues[position - 1] == residue:
            mods.setdefault(position, []).append(raw_mod)
    return mods


def _parse_mass_token(token: str) -> float | None:
    match = re.fullmatch(r"[+-]?\d+(?:\.\d+)?", token.strip())
    if not match:
        return None
    try:
        return float(token)
    except ValueError:
        return None


def _resolve_mass_to_unimod(mass: float) -> str | None:
    best_match: str | None = None
    best_delta: float | None = None
    mass_table = dict(UNIMOD_TO_MASS)
    catalog = _load_unimod_catalog()
    if catalog is not None:
        try:
            for entry in catalog.mods:
                if not entry.get("approved", True):
                    continue
                try:
                    mass_table[f"UNIMOD:{int(entry['record_id'])}"] = float(entry["mono_mass"])
                except Exception:
                    continue
        except Exception:
            pass
    for unimod, ref_mass in mass_table.items():
        delta = abs(mass - ref_mass)
        if delta <= MASS_MATCH_TOLERANCE and (best_delta is None or delta < best_delta):
            best_match = unimod
            best_delta = delta
    return best_match
