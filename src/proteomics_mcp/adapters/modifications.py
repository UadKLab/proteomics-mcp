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
    "deamidation": "UNIMOD:7",
    "amidated": "UNIMOD:2",
    "amidation": "UNIMOD:2",
    "dimethyl": "UNIMOD:36",
    "gg": "UNIMOD:121",
    "glygly": "UNIMOD:121",
    "itraq4plex": "UNIMOD:214",
    "ammonia-loss": "UNIMOD:385",
    "tmt6plex": "UNIMOD:737",
    "tmt10plex": "UNIMOD:737",
    "tmt11plex": "UNIMOD:737",
    "tmtpro": "UNIMOD:2016",
}

UNIMOD_TO_MASS = {
    "UNIMOD:1": 42.01056468403,
    "UNIMOD:4": 57.021463734999996,
    "UNIMOD:5": 43.005813,
    "UNIMOD:7": 0.9840155848,
    "UNIMOD:21": 79.9663308946,
    "UNIMOD:27": -18.01056468403,
    # Same mass as UNIMOD:28; listed first so a bare -17.03 resolves to the generic ammonia loss, not the Gln-specific one.
    "UNIMOD:385": -17.02654910112,
    "UNIMOD:28": -17.02654910112,
    "UNIMOD:35": 15.99491461957,
    "UNIMOD:34": 14.01565006446,
    "UNIMOD:122": 27.99491461957,
    "UNIMOD:2": -0.984016,
    "UNIMOD:36": 28.0313,
    "UNIMOD:121": 114.042927,
    "UNIMOD:214": 144.102063,
    "UNIMOD:737": 229.162932,
    "UNIMOD:2016": 304.207146,
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
    "UNIMOD:2": "amidated",
    "UNIMOD:36": "dimethyl",
    "UNIMOD:121": "gg",
    "UNIMOD:214": "itraq4plex",
    "UNIMOD:385": "nh3loss",
    "UNIMOD:737": "tmt",
    "UNIMOD:2016": "tmtpro",
}

MASS_MATCH_TOLERANCE = 0.01
# A signed or unsigned mass, including forms like "+.98" (InstaNovo) and "+42.011" (Casanovo).
_MASS_TOKEN = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)"
_BARE_MASS_RE = re.compile(r"[+-](?:\d+(?:\.\d*)?|\.\d+)")
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


def _without_site_suffix(name: str) -> str:
    match = re.fullmatch(r"(.+?)\s+\([^()]*\)", name)
    return match.group(1).strip() if match else name


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
    # MaxQuant writes "Title (Sites)", e.g. "Oxidation (M)" or "Acetyl (Protein N-term)".
    title = _without_site_suffix(clean)
    if title != clean and NAME_TO_UNIMOD.get(title.lower()):
        return NAME_TO_UNIMOD[title.lower()]
    mod = _find_unimod_entry(clean)
    if mod is None and title != clean:
        mod = _find_unimod_entry(title)
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

    text = bracket_modification_tokens(text)
    text = re.sub(r"\[([^\]]+)\]", replace_mod, text)
    # Leading tokens are N-terminal: "[UNIMOD:1]-PEPTIDE". Tokens followed by "?" have an unknown position.
    text = re.sub(r"^((?:\[[^\]]+\])+)(?=[A-Z])", r"\1-", text)
    return text


def unresolved_modification_tokens(peptidoform: str | None) -> list[str]:
    """Modification tokens of a normalized peptidoform that are neither a UNIMOD accession nor a mass shift.

    normalize_peptidoform_to_unimod_style keeps such tokens as written ("[Foo]"); a mass that matches no UNIMOD entry
    stays a mass ("[+123.4567]") and still scores, so it is not listed.
    """
    if not peptidoform:
        return []
    return [token for token in re.findall(r"\[([^\]]+)\]", str(peptidoform)) if not token.strip().upper().startswith("UNIMOD:") and _parse_mass_token(token) is None]


def matching_parenthesis(text: str, start: int) -> int | None:
    """Index of the ")" that closes the "(" at text[start], allowing nesting such as "(Oxidation (M))"."""
    depth = 0
    for index in range(start, len(text)):
        if text[index] == "(":
            depth += 1
        elif text[index] == ")":
            depth -= 1
            if depth == 0:
                return index
    return None


def _outside_brackets(text: str, transform: Any) -> str:
    parts = re.split(r"(\[[^\]]*\])", text)
    return "".join(part if part.startswith("[") and part.endswith("]") else transform(part) for part in parts)


def bracket_modification_tokens(peptidoform: str) -> str:
    """Put every modification token of a peptidoform string in square brackets, without resolving it.

    Handles parenthesised tokens anywhere, including nested ones and N- or C-terminal ones ("M(ox)", "(+42.01)PEPTIDE",
    "M(Oxidation (M))", "PEPTIDE-(+.98)"), and unbracketed signed masses ("+42.011PEPM+15.995K"). Spaces outside
    brackets are removed.
    """
    out: list[str] = []
    i = 0
    while i < len(peptidoform):
        ch = peptidoform[i]
        if ch == "[":
            end = peptidoform.find("]", i + 1)
            if end == -1:
                out.append(peptidoform[i:])
                break
            out.append(peptidoform[i : end + 1])
            i = end + 1
            continue
        if ch == "(":
            end = matching_parenthesis(peptidoform, i)
            if end is None:
                out.append(peptidoform[i:])
                break
            out.append(f"[{peptidoform[i + 1 : end].strip()}]")
            i = end + 1
            continue
        out.append(ch)
        i += 1
    return _outside_brackets(
        "".join(out), lambda part: _BARE_MASS_RE.sub(lambda match: f"[{match.group(0)}]", part.replace(" ", ""))
    )


def build_peptidoform(
    sequence: str | None,
    *,
    residue_mods: dict[int, list[str]] | None = None,
    nterm_mods: list[str] | None = None,
    cterm_mods: list[str] | None = None,
    unlocalized_mods: list[str] | None = None,
) -> str | None:
    """Build a peptidoform from a sequence and modifications by 1-based residue position.

    unlocalized_mods are written as ProForma unknown-position tokens ("[Phospho]?PEPTIDE"); they keep the modification
    visible, and parse_peptidoform flags it as unsupported.
    """
    if sequence is None:
        return None
    residues = [ch.upper() for ch in str(sequence) if ch.isalpha()]
    if not residues:
        return None
    residue_mods = residue_mods or {}
    nterm_mods = nterm_mods or []
    cterm_mods = cterm_mods or []

    def token(mod: str) -> str:
        # Square brackets delimit tokens, so they cannot appear inside one.
        return f"[{(resolve_mod_token(mod) or mod).replace('[', '(').replace(']', ')')}]"

    peptide = "".join(
        f"{residue}{''.join(token(mod) for mod in residue_mods.get(idx, []))}" for idx, residue in enumerate(residues, start=1)
    )
    if nterm_mods:
        peptide = f"{''.join(token(mod) for mod in nterm_mods)}-{peptide}"
    if cterm_mods:
        peptide = f"{peptide}-{''.join(token(mod) for mod in cterm_mods)}"
    if unlocalized_mods:
        peptide = f"{''.join(f'{token(mod)}?' for mod in unlocalized_mods)}{peptide}"
    return normalize_peptidoform_to_unimod_style(peptide)


def parse_fragpipe_modifications(
    value: Any, sequence: str | None
) -> tuple[dict[int, list[str]], list[str], list[str], list[str]]:
    """Parse a FragPipe "Assigned Modifications" value such as "N-term(42.0106), 5M(15.9949)".

    Returns (residue modifications by 1-based position, N-terminal, C-terminal, parts not understood).
    """
    text = str(value).strip() if value is not None and not (isinstance(value, float) and math.isnan(value)) else ""
    if not text or not sequence:
        return {}, [], [], []
    residues = [ch.upper() for ch in str(sequence) if ch.isalpha()]
    mods: dict[int, list[str]] = {}
    nterm: list[str] = []
    cterm: list[str] = []
    unparsed: list[str] = []
    for part in [segment.strip() for segment in text.split(",") if segment.strip()]:
        match = re.fullmatch(r"(?P<position>\d+)(?P<residue>[A-Za-z])\((?P<mod>[^()]+)\)", part)
        if match:
            position = int(match.group("position"))
            if 1 <= position <= len(residues) and residues[position - 1] == match.group("residue").upper():
                mods.setdefault(position, []).append(match.group("mod").strip())
            else:
                unparsed.append(part)
            continue
        terminal = re.fullmatch(r"(?P<term>[NnCc])-?term\((?P<mod>[^()]+)\)", part)
        if terminal:
            (nterm if terminal.group("term").upper() == "N" else cterm).append(terminal.group("mod").strip())
            continue
        unparsed.append(part)
    return mods, nterm, cterm, unparsed


def parse_fragpipe_assigned_modifications(value: Any, sequence: str | None) -> dict[int, list[str]]:
    """Residue modifications only; see parse_fragpipe_modifications for terminal ones."""
    return parse_fragpipe_modifications(value, sequence)[0]


def nearest_builtin_unimod(mass: float, tolerance: float) -> str | None:
    """The built-in modification closest to a mass delta, within tolerance (no catalogue lookup)."""
    best = min(UNIMOD_TO_MASS.items(), key=lambda item: abs(item[1] - mass))
    return best[0] if abs(best[1] - mass) <= tolerance else None


def _parse_mass_token(token: str) -> float | None:
    match = re.fullmatch(_MASS_TOKEN, token.strip())
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
