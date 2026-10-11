from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import Any

import numpy as np

from proteomics_mcp.adapters.modifications import unimod_mass
from proteomics_mcp.models.schemas import ObservedSpectrum

PROTON_MASS = 1.007276466812
HYDROGEN_MASS = 1.00782503223
H2O_MASS = 18.01056468403
NH3_MASS = 17.02654910112
CO_MASS = 27.99491461957
NEUTRON_MASS = 1.00335483507

AA_RESIDUE_MASS = {
    "A": 71.037113805,
    "R": 156.101111050,
    "N": 114.042927470,
    "D": 115.026943065,
    "C": 103.009184505,
    "E": 129.042593135,
    "Q": 128.058577540,
    "G": 57.021463735,
    "H": 137.058911875,
    "I": 113.084064015,
    "L": 113.084064015,
    "K": 128.094963050,
    "M": 131.040484645,
    "F": 147.068413945,
    "P": 97.052763875,
    "S": 87.032028435,
    "T": 101.047678505,
    "W": 186.079312980,
    "Y": 163.063328575,
    "V": 99.068413945,
}

SIDECHAIN_ION_LOSSES = {
    "C": {"primary": ("SH", 32.979896), "alt": None},
    "D": {"primary": ("CO2H", 44.997654), "alt": None},
    "E": {"primary": ("C2H3O2", 59.013304), "alt": None},
    "I": {"primary": ("C2H5", 29.039125), "alt": ("CH3", 15.023475)},
    "K": {"primary": ("C3H8N", 58.065674), "alt": None},
    "L": {"primary": ("C3H7", 43.054775), "alt": None},
    "M": {"primary": ("C2H5S", 61.011196), "alt": None},
    "N": {"primary": ("CONH2", 44.013639), "alt": None},
    "Q": {"primary": ("C2H4NO", 58.029289), "alt": None},
    "V": {"primary": ("CH3", 15.023475), "alt": None},
}

UNIMOD_TO_OPENMS = {
    "1": "Acetyl (N-term)",
    "4": "Carbamidomethyl (C)",
    "7": {"N": "Deamidated (N)", "Q": "Deamidated (Q)", "R": "Deamidated (R)"},
    "21": {"S": "Phospho (S)", "T": "Phospho (T)", "Y": "Phospho (Y)"},
    "35": "Oxidation (M)",
}

UNIMOD_MASS_DELTA = {
    "1": 42.01056468403,
    "4": 57.021463734999996,
    "7": 0.9840155848,
    "21": 79.9663308946,
    "35": 15.99491461957,
}

MOD_NAME_MASS_DELTA = {
    "acetyl": 42.01056468403,
    "carbamidomethyl": 57.021463734999996,
    "deamidated": 0.9840155848,
    "phospho": 79.9663308946,
    "oxidation": 15.99491461957,
}
RADICAL_DOT = "\u2022"

REPORTER_ION_MZ = {
    "TMT126": 126.127726,
    "TMT127N": 127.124761,
    "TMT127C": 127.131081,
    "TMT128N": 128.128116,
    "TMT128C": 128.134436,
    "TMT129N": 129.131471,
    "TMT129C": 129.137790,
    "TMT130N": 130.134825,
    "TMT130C": 130.141145,
    "TMT131N": 131.138180,
    "TMT131C": 131.144500,
    "TMT132N": 132.141535,
    "TMT132C": 132.147855,
    "TMT133N": 133.144890,
    "TMT133C": 133.151210,
    "TMT134N": 134.148245,
    "TMT134C": 134.154565,
    "TMT135N": 135.151600,
    "iTRAQ114": 114.110679,
    "iTRAQ115": 115.107714,
    "iTRAQ116": 116.111069,
    "iTRAQ117": 117.114424,
    "iTRAQ113": 113.107324,
    "iTRAQ118": 118.111459,
    "iTRAQ119": 119.114814,
    "iTRAQ121": 121.121524,
}

IMMONIUM_ION_MZ = {
    "A": 44.0494757,
    "R": 129.1140232,
    "N": 87.0552899,
    "D": 88.0393046,
    "C": 76.0221029,
    "E": 102.0549547,
    "Q": 101.0709400,
    "G": 30.0348256,
    "H": 110.0712691,
    "I": 86.0964256,
    "L": 86.0964256,
    "K": 101.1073246,
    "M": 104.0528467,
    "F": 120.0807755,
    "P": 70.0651257,
    "S": 60.0443903,
    "T": 74.0600403,
    "W": 159.0916735,
    "Y": 136.0756901,
    "V": 72.0807755,
}

STANDARD_NEUTRAL_LOSSES = {"H2O": H2O_MASS, "NH3": NH3_MASS}
PHOSPHO_NEUTRAL_LOSSES = {"H3PO4": 97.976895573}


@dataclass(slots=True)
class ParsedPeptidoform:
    residues: list[str]
    residue_masses: list[float]
    residue_mod_names: list[list[str]]
    residue_mod_deltas: list[list[float]]
    nterm_delta: float = 0.0
    cterm_delta: float = 0.0
    unsupported_modification: bool = False


def ppm_error(observed_mz: float, theoretical_mz: float) -> float:
    return ((observed_mz - theoretical_mz) / theoretical_mz) * 1e6


def _clean_peptidoform(peptidoform: str) -> str:
    text = peptidoform.strip()
    if "/" in text:
        text = text.rsplit("/", 1)[0]
    flanked = re.match(r"^(?:\[[^\]]*\]|[A-Za-z-]+)\.(.+)\.(?:\[[^\]]*\]|[A-Za-z-]+)$", text)
    if flanked:
        text = flanked.group(1)
    return text


def _mod_delta(token: str, residue: str | None = None) -> float | None:
    normalized = token.strip().strip("[]")
    if not normalized:
        return None
    unimod_delta = unimod_mass(normalized)
    if unimod_delta is not None:
        return unimod_delta
    try:
        if normalized[0] in "+-":
            return float(normalized)
    except (IndexError, ValueError):
        return None
    try:
        numeric = float(normalized)
    except ValueError:
        numeric = None
    if numeric is not None and math.isfinite(numeric):
        return numeric
    lowered = normalized.lower()
    for name, delta in MOD_NAME_MASS_DELTA.items():
        if lowered.startswith(name.lower()):
            return delta
    return None


def _mod_name(token: str) -> str:
    normalized = token.strip().strip("[]")
    if normalized.upper().startswith("UNIMOD:"):
        return normalized.upper().replace("UNIMOD:", "UNIMOD:")
    return normalized


def parse_peptidoform(peptidoform: str) -> ParsedPeptidoform:
    text = _clean_peptidoform(peptidoform)
    residues: list[str] = []
    residue_masses: list[float] = []
    residue_mod_names: list[list[str]] = []
    residue_mod_deltas: list[list[float]] = []
    nterm_delta = 0.0
    cterm_delta = 0.0
    unsupported_modification = False
    i = 0
    pending_nterm = True
    while i < len(text):
        ch = text[i]
        if ch == "[":
            end = text.find("]", i + 1)
            if end == -1:
                unsupported_modification = True
                break
            token = text[i + 1 : end]
            delta = _mod_delta(token, residues[-1] if residues else None)
            if residues:
                residue_mod_names[-1].append(_mod_name(token))
                if delta is not None:
                    residue_mod_deltas[-1].append(delta)
                    residue_masses[-1] += delta
                else:
                    unsupported_modification = True
            elif pending_nterm:
                if delta is not None:
                    nterm_delta += delta
                else:
                    unsupported_modification = True
            else:
                if delta is not None:
                    cterm_delta += delta
                else:
                    unsupported_modification = True
            i = end + 1
            if i < len(text) and text[i] == "-":
                i += 1
            continue
        if "A" <= ch <= "Z" and ch in AA_RESIDUE_MASS:
            residues.append(ch)
            residue_masses.append(AA_RESIDUE_MASS[ch])
            residue_mod_names.append([])
            residue_mod_deltas.append([])
            pending_nterm = False
        i += 1
    return ParsedPeptidoform(
        residues=residues,
        residue_masses=residue_masses,
        residue_mod_names=residue_mod_names,
        residue_mod_deltas=residue_mod_deltas,
        nterm_delta=nterm_delta,
        cterm_delta=cterm_delta,
        unsupported_modification=unsupported_modification,
    )


def residues_from_peptidoform(peptidoform: str) -> list[str]:
    """Return residue letters from a ProForma-like peptide string.

    The internal annotator currently ignores modification mass shifts. Modified
    precursor masses still use pyOpenMS through aa_sequence_from_proforma where
    available; fragment annotation remains intentionally simple and transparent.
    """
    return parse_peptidoform(peptidoform).residues


def aa_sequence_from_proforma(peptidoform: str):
    try:
        import pyopenms  # type: ignore
    except ImportError:
        return None

    text = _clean_peptidoform(peptidoform)
    residues: list[str] = []
    mods: list[str | None] = []
    nterm_mod: str | None = None
    i = 0
    pending_nterm = True
    while i < len(text):
        if text.startswith("[UNIMOD:", i):
            end = text.find("]", i)
            if end == -1:
                return None
            unimod = text[i + 8 : end]
            mapped = UNIMOD_TO_OPENMS.get(unimod)
            if pending_nterm and isinstance(mapped, str):
                nterm_mod = mapped
            i = end + 1
            if i < len(text) and text[i] == "-":
                i += 1
            continue
        ch = text[i]
        if "A" <= ch <= "Z":
            residues.append(ch)
            mods.append(None)
            pending_nterm = False
            i += 1
            if text.startswith("[UNIMOD:", i):
                end = text.find("]", i)
                if end == -1:
                    return None
                unimod = text[i + 8 : end]
                mapped = UNIMOD_TO_OPENMS.get(unimod)
                if isinstance(mapped, dict):
                    mods[-1] = mapped.get(ch)
                elif isinstance(mapped, str):
                    mods[-1] = mapped
                i = end + 1
            continue
        i += 1
    if not residues:
        return None
    try:
        seq = pyopenms.AASequence.fromString("".join(residues))
        if nterm_mod:
            seq.setNTerminalModification(nterm_mod)
        for idx, mod in enumerate(mods):
            if mod:
                seq.setModification(idx, mod)
        return seq
    except Exception:
        return None


def _charge_label(charge: int) -> str:
    return "+" * max(charge, 1)


def _charged_mz(neutral_mass: float, charge: int) -> float:
    return (neutral_mass + charge * PROTON_MASS) / charge


BASE_ION_TYPES = {"a", "b", "c", "x", "y", "z"}
EXOTIC_ION_TYPES = {"a+1", "z-dot", "d", "d-alt", "w", "w-alt"}
AUXILIARY_ION_TYPES = {"immonium", "reporter", "precursor"}
ALL_ION_TYPES = BASE_ION_TYPES | EXOTIC_ION_TYPES | AUXILIARY_ION_TYPES


def _parse_ion_types(ion_types: str) -> list[str]:
    text = str(ion_types or "").strip().lower()
    if text in {"all", "exotic"}:
        return ["a", "a+1", "b", "c", "d", "d-alt", "x", "y", "z", "z-dot", "w", "w-alt", "immonium", "reporter", "precursor"] if text == "all" else ["a+1", "d", "d-alt", "z-dot", "w", "w-alt"]
    alias_map = {
        "i": "immonium",
        "immonium": "immonium",
        "r": "reporter",
        "reporter": "reporter",
        "p": "precursor",
        "precursor": "precursor",
    }
    tokens: list[str] = []
    normalized = text.replace("z_dot", "z-dot").replace("zdot", "z-dot").replace("z*", "z-dot").replace("z.", "z-dot")
    for word, alias in alias_map.items():
        if word in normalized.split(","):
            tokens.append(alias)
            normalized = normalized.replace(word, "")
    for multi in ("d-alt", "w-alt", "a+1", "z-dot"):
        if multi in normalized:
            tokens.append(multi)
            normalized = normalized.replace(multi, "")
    if "," in normalized or " " in normalized:
        for token in re.split(r"[\s,]+", normalized):
            if not token:
                continue
            if token in alias_map:
                tokens.append(alias_map[token])
            elif token in BASE_ION_TYPES or token in {"d", "w"}:
                tokens.append(token)
    else:
        for ion in normalized:
            if ion in alias_map:
                tokens.append(alias_map[ion])
            elif ion in BASE_ION_TYPES or ion in {"d", "w"}:
                tokens.append(ion)
    return list(dict.fromkeys(tokens))


def _neutral_loss_options(parsed: ParsedPeptidoform, neutral_losses: dict[str, float] | bool | None) -> dict[str, float]:
    if isinstance(neutral_losses, dict):
        return {str(k): float(v) for k, v in neutral_losses.items()}
    if not neutral_losses:
        return {}
    losses = dict(STANDARD_NEUTRAL_LOSSES)
    mod_names = {name.lower() for names in parsed.residue_mod_names for name in names}
    if any("phospho" in name or name == "unimod:21" for name in mod_names):
        losses.update(PHOSPHO_NEUTRAL_LOSSES)
    return losses


def _ion_neutral_mass(ion_type: str, prefix_mass: float, suffix_mass: float, residue: str | None = None) -> tuple[float, str | None] | None:
    if ion_type == "a":
        return prefix_mass - CO_MASS, None
    if ion_type == "a+1":
        return prefix_mass - CO_MASS + HYDROGEN_MASS, None
    if ion_type == "b":
        return prefix_mass, None
    if ion_type == "c":
        return prefix_mass + NH3_MASS, None
    if ion_type in {"d", "d-alt"}:
        loss_key = "alt" if ion_type == "d-alt" else "primary"
        loss = SIDECHAIN_ION_LOSSES.get(str(residue), {}).get(loss_key) if residue else None
        if loss is None:
            return None
        return prefix_mass - CO_MASS - float(loss[1]) + HYDROGEN_MASS, str(loss[0])
    if ion_type == "x":
        return suffix_mass + H2O_MASS + CO_MASS - 2 * HYDROGEN_MASS, None
    if ion_type == "y":
        return suffix_mass + H2O_MASS, None
    if ion_type == "z":
        return suffix_mass + H2O_MASS - NH3_MASS, None
    if ion_type == "z-dot":
        return suffix_mass + H2O_MASS - NH3_MASS + HYDROGEN_MASS, None
    if ion_type in {"w", "w-alt"}:
        loss_key = "alt" if ion_type == "w-alt" else "primary"
        loss = SIDECHAIN_ION_LOSSES.get(str(residue), {}).get(loss_key) if residue else None
        if loss is None:
            return None
        return suffix_mass + H2O_MASS - NH3_MASS - float(loss[1]) + HYDROGEN_MASS, str(loss[0])
    return None


def generate_terminal_fragments(
    peptidoform: str,
    *,
    ion_types: str = "by",
    max_ion_charge: int = 1,
    max_isotope: int = 0,
    neutral_losses: dict[str, float] | bool | None = None,
) -> list[dict[str, Any]]:
    parsed = parse_peptidoform(peptidoform)
    residues = parsed.residues
    if len(residues) < 2:
        return []
    residue_masses = parsed.residue_masses
    prefix_masses: list[float] = []
    running = float(parsed.nterm_delta)
    for mass in residue_masses[:-1]:
        running += mass
        prefix_masses.append(running)
    suffix_masses: list[float] = []
    running = float(parsed.cterm_delta)
    for mass in reversed(residue_masses[1:]):
        running += mass
        suffix_masses.append(running)
    suffix_masses.reverse()

    fragments: list[dict[str, Any]] = []
    allowed = [ion for ion in _parse_ion_types(ion_types) if ion not in AUXILIARY_ION_TYPES]
    loss_options = _neutral_loss_options(parsed, neutral_losses)
    peptide_length = len(residues)
    for ordinal, (prefix_mass, suffix_mass) in enumerate(zip(prefix_masses, suffix_masses), start=1):
        for ion_type in allowed:
            n_term_ion = ion_type in {"a", "a+1", "b", "c", "d", "d-alt"}
            residue = residues[ordinal - 1] if n_term_ion else residues[peptide_length - ordinal]
            ion_ordinal = ordinal if n_term_ion else peptide_length - ordinal
            ion = _ion_neutral_mass(ion_type, prefix_mass, suffix_mass, residue=residue)
            if ion is None:
                continue
            neutral_mass, sidechain_loss = ion
            if neutral_mass <= 0:
                continue
            loss_variants: list[tuple[str | None, float]] = [(None, neutral_mass)]
            for loss_name, loss_mass in loss_options.items():
                candidate_mass = neutral_mass - loss_mass
                if candidate_mass > 0:
                    loss_variants.append((loss_name, candidate_mass))
            for neutral_loss_name, variant_mass in loss_variants:
                for charge in range(1, max(1, int(max_ion_charge)) + 1):
                    base_mz = _charged_mz(variant_mass, charge)
                    for isotope in range(0, max(0, int(max_isotope)) + 1):
                        mz = base_mz + isotope * NEUTRON_MASS / charge
                        isotope_label = f"[M+{isotope}]" if isotope else ""
                        if ion_type == "a+1":
                            annotation = f"a{ion_ordinal}{RADICAL_DOT}{_charge_label(charge)}"
                        elif ion_type == "z-dot":
                            annotation = f"z{ion_ordinal}{RADICAL_DOT}{_charge_label(charge)}"
                        else:
                            annotation = f"{ion_type}{ion_ordinal}{_charge_label(charge)}"
                        if neutral_loss_name:
                            annotation += f"-{neutral_loss_name}"
                        annotation += isotope_label
                        fragments.append(
                            {
                                "annotation": annotation,
                                "mz": mz,
                                "ion_type": ion_type,
                                "ordinal": ion_ordinal,
                                "charge": charge,
                                "isotope": isotope,
                                "neutral_loss": neutral_loss_name,
                                "residue": residue,
                                "sidechain_loss": sidechain_loss,
                            }
                        )
    return fragments


def generate_immonium_ions(peptidoform: str) -> list[dict[str, Any]]:
    parsed = parse_peptidoform(peptidoform)
    ions: list[dict[str, Any]] = []
    seen: set[tuple[str, int]] = set()
    for ordinal, residue in enumerate(parsed.residues, start=1):
        mz = IMMONIUM_ION_MZ.get(residue)
        if mz is None or (residue, ordinal) in seen:
            continue
        seen.add((residue, ordinal))
        ions.append(
            {
                "annotation": f"Immonium({residue})",
                "mz": mz,
                "ion_type": "immonium",
                "ordinal": ordinal,
                "charge": 1,
                "isotope": 0,
                "residue": residue,
            }
        )
    return ions


def generate_reporter_ions() -> list[dict[str, Any]]:
    return [
        {
            "annotation": name,
            "mz": mz,
            "ion_type": "reporter",
            "ordinal": None,
            "charge": 1,
            "isotope": 0,
        }
        for name, mz in REPORTER_ION_MZ.items()
    ]


def generate_precursor_related_ions(
    peptidoform: str,
    *,
    precursor_charge: int | None,
    max_isotope: int = 0,
    neutral_losses: dict[str, float] | bool | None = None,
) -> list[dict[str, Any]]:
    seq = aa_sequence_from_proforma(peptidoform)
    if seq is None:
        return []
    charge_limit = max(1, int(precursor_charge or 1))
    parsed = parse_peptidoform(peptidoform)
    loss_options = _neutral_loss_options(parsed, neutral_losses)
    ions: list[dict[str, Any]] = []
    for charge in range(1, charge_limit + 1):
        neutral_mass = float(seq.getMZ(1)) - PROTON_MASS
        if neutral_mass is None or neutral_mass <= 0:
            continue
        loss_variants: list[tuple[str | None, float]] = [(None, neutral_mass)]
        for loss_name, loss_mass in loss_options.items():
            candidate_mass = neutral_mass - loss_mass
            if candidate_mass > 0:
                loss_variants.append((loss_name, candidate_mass))
        for neutral_loss_name, variant_mass in loss_variants:
            for isotope in range(0, max(0, int(max_isotope)) + 1):
                mz = _charged_mz(variant_mass, charge) + isotope * NEUTRON_MASS / charge
                annotation = f"p{_charge_label(charge)}"
                if neutral_loss_name:
                    annotation += f"-{neutral_loss_name}"
                if isotope:
                    annotation += f"[M+{isotope}]"
                ions.append(
                    {
                        "annotation": annotation,
                        "mz": mz,
                        "ion_type": "precursor",
                        "ordinal": None,
                        "charge": charge,
                        "isotope": isotope,
                        "neutral_loss": neutral_loss_name,
                    }
                )
    return ions


def generate_theoretical_ions(
    peptidoform: str,
    *,
    ion_types: str = "by",
    precursor_charge: int | None = None,
    max_ion_charge: int = 1,
    max_isotope: int = 0,
    neutral_losses: dict[str, float] | bool | None = None,
) -> list[dict[str, Any]]:
    allowed = set(_parse_ion_types(ion_types))
    fragments = generate_terminal_fragments(
        peptidoform,
        ion_types=",".join(sorted(allowed - AUXILIARY_ION_TYPES)),
        max_ion_charge=max_ion_charge,
        max_isotope=max_isotope,
        neutral_losses=neutral_losses,
    )
    if "immonium" in allowed:
        fragments.extend(generate_immonium_ions(peptidoform))
    if "reporter" in allowed:
        fragments.extend(generate_reporter_ions())
    if "precursor" in allowed:
        fragments.extend(
            generate_precursor_related_ions(
                peptidoform,
                precursor_charge=precursor_charge,
                max_isotope=max_isotope,
                neutral_losses=neutral_losses,
            )
        )
    return fragments


def _tolerance_error(observed_mz: float, theoretical_mz: float, tolerance_unit: str) -> float:
    return ppm_error(observed_mz, theoretical_mz) if tolerance_unit == "ppm" else observed_mz - theoretical_mz


def annotate_observed_ions(
    observed: ObservedSpectrum,
    *,
    peptidoform: str,
    tolerance: float = 0.05,
    tolerance_unit: str = "Da",
    ion_types: str = "by",
    max_isotope: int = 0,
    max_ion_charge: int | None = None,
    neutral_losses: dict[str, float] | bool | None = None,
) -> dict[str, Any]:
    """Annotate an observed spectrum against internally generated terminal ions."""
    charge_limit = max_ion_charge or max((observed.precursor_charge or 2) - 1, 1)
    theoretical = generate_theoretical_ions(
        peptidoform,
        ion_types=ion_types,
        precursor_charge=observed.precursor_charge,
        max_ion_charge=charge_limit,
        max_isotope=max_isotope,
        neutral_losses=neutral_losses,
    )
    labels = [""] * len(observed.mz)
    matches: list[dict[str, Any]] = []
    counts: dict[str, int] = {}
    intensity_by_ion_type: dict[str, float] = {}
    total_intensity = sum(float(i) for i in observed.intensity)
    annotated_intensity = 0.0
    indexed_theoretical = sorted(enumerate(theoretical), key=lambda item: float(item[1]["mz"]))
    theoretical_mz = np.asarray([float(fragment["mz"]) for _, fragment in indexed_theoretical], dtype=float)
    for obs_idx, obs_mz in enumerate(observed.mz):
        best: dict[str, Any] | None = None
        best_abs_error = math.inf
        best_original_index = math.inf
        window = tolerance if tolerance_unit == "Da" else abs(float(obs_mz)) * tolerance / 1e6
        left = int(np.searchsorted(theoretical_mz, float(obs_mz) - window, side="left"))
        right = int(np.searchsorted(theoretical_mz, float(obs_mz) + window, side="right"))
        for original_index, fragment in indexed_theoretical[left:right]:
            err = _tolerance_error(float(obs_mz), float(fragment["mz"]), tolerance_unit)
            abs_err = abs(err)
            if abs_err <= tolerance and (abs_err < best_abs_error or (abs_err == best_abs_error and original_index < best_original_index)):
                best = fragment
                best_abs_error = abs_err
                best_original_index = original_index
        if best is None:
            continue
        label = str(best["annotation"])
        labels[obs_idx] = label
        ion_type = str(best["ion_type"])
        counts[ion_type] = counts.get(ion_type, 0) + 1
        observed_intensity = float(observed.intensity[obs_idx])
        intensity_by_ion_type[ion_type] = intensity_by_ion_type.get(ion_type, 0.0) + observed_intensity
        annotated_intensity += observed_intensity
        matched_theoretical_mz = float(best["mz"])
        matches.append(
            {
                "observed_index": obs_idx,
                "observed_mz": float(obs_mz),
                "observed_intensity": observed_intensity,
                "annotation": label,
                "theoretical_mz": matched_theoretical_mz,
                "mz_error": _tolerance_error(float(obs_mz), matched_theoretical_mz, tolerance_unit),
                "mz_error_ppm": ppm_error(float(obs_mz), matched_theoretical_mz),
                "mz_error_unit": tolerance_unit,
                "ion_type": ion_type,
                "ordinal": best["ordinal"],
                "charge": best["charge"],
                "isotope": best["isotope"],
                "residue": best.get("residue"),
                "sidechain_loss": best.get("sidechain_loss"),
            }
        )
    parsed_ion_types = _parse_ion_types(ion_types)
    return {
        "engine": "proteomics_mcp_terminal_ion_annotator",
        "spectrum_id": observed.spectrum_id,
        "scan_number": observed.scan_number,
        "peptidoform": peptidoform,
        "tolerance": tolerance,
        "tolerance_unit": tolerance_unit,
        "ion_types": ",".join(parsed_ion_types),
        "max_isotope": max_isotope,
        "max_ion_charge": charge_limit,
        "neutral_losses": neutral_losses,
        "annotations": labels,
        "matches": matches,
        "annotated_peak_count": len(matches),
        "ion_type_counts": counts,
        "annotated_intensity": annotated_intensity,
        "total_observed_intensity": total_intensity,
        "annotated_intensity_fraction": annotated_intensity / total_intensity if total_intensity > 0 else None,
        "ion_type_intensity": intensity_by_ion_type,
        "total_theoretical_ions": len(theoretical),
        "unsupported_modification_warning": parse_peptidoform(peptidoform).unsupported_modification,
    }
