"""Comet XCorr and the X!Tandem hyperscore, reimplemented from the engines' source code (#51).

Both scores compare the observed spectrum with the candidate's *theoretical* fragment ions, not with a predicted
spectrum, so they do not depend on the fragment-intensity predictor. They follow, step by step:

- Comet 2026.01.1 (github.com/UWPR/Comet, tag v2026.01.1): CometPreprocess.cpp (LoadIons, MakeCorrData,
  Preprocess) and CometSearch.cpp (XcorrScore), with the masses of CometMassSpecUtils.cpp.
- X! Tandem 2015.12.15.2 (the Bioconda xtandem source): mspectrumcondition.cpp (condition), mscore.cpp
  (add_B, add_Y, score) and mscore_tandem.cpp (blur, dot, hfactor, hconvert), with the masses of masscalc.cpp.

validation/engine_scores/ checks them against the engines themselves, run from pinned BioContainers images.
Each function cites the engine code it follows; where an engine works in single precision, so does this module.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from proteomics_mcp.core.annotation import parse_peptidoform
from proteomics_mcp.models.schemas import ObservedSpectrum

COMET_VERSION = "2026.01.1"
XTANDEM_VERSION = "2015.12.15.2"

# --- residue masses ------------------------------------------------------------------------------------------------

# Elemental formulas of the residues, as both engines build them (Comet: CometMassSpecUtils::AssignMass;
# X! Tandem: msequtilities::set_aa). Elements: C, H, N, O, S, Se.
_RESIDUE_FORMULAS = {
    "G": (2, 3, 1, 1, 0, 0),
    "A": (3, 5, 1, 1, 0, 0),
    "S": (3, 5, 1, 2, 0, 0),
    "P": (5, 7, 1, 1, 0, 0),
    "V": (5, 9, 1, 1, 0, 0),
    "T": (4, 7, 1, 2, 0, 0),
    "C": (3, 5, 1, 1, 1, 0),
    "U": (3, 5, 1, 1, 0, 1),
    "L": (6, 11, 1, 1, 0, 0),
    "I": (6, 11, 1, 1, 0, 0),
    "N": (4, 6, 2, 2, 0, 0),
    "D": (4, 5, 1, 3, 0, 0),
    "Q": (5, 8, 2, 2, 0, 0),
    "K": (6, 12, 2, 1, 0, 0),
    "E": (5, 7, 1, 3, 0, 0),
    "M": (5, 9, 1, 1, 1, 0),
    "H": (6, 7, 3, 1, 0, 0),
    "F": (9, 9, 1, 1, 0, 0),
    "R": (6, 12, 4, 1, 0, 0),
    "Y": (9, 9, 1, 2, 0, 0),
    "W": (11, 10, 2, 1, 0, 0),
    "O": (12, 19, 3, 2, 0, 0),
}


def _residue_masses(c: float, h: float, n: float, o: float, s: float, se: float) -> dict[str, float]:
    return {aa: f[0] * c + f[1] * h + f[2] * n + f[3] * o + f[4] * s + f[5] * se for aa, f in _RESIDUE_FORMULAS.items()}


# Comet: CometMassSpecUtils.h / CometMassSpecUtils::AssignMass; PROTON_MASS from CometDataInternal.h
_COMET_H, _COMET_O, _COMET_C, _COMET_N, _COMET_S, _COMET_SE = 1.007825035, 15.99491463, 12.0, 14.0030740, 31.9720707, 79.9165196
COMET_PROTON = 1.00727646688
COMET_RESIDUE_MASS = _residue_masses(_COMET_C, _COMET_H, _COMET_N, _COMET_O, _COMET_S, _COMET_SE)
_COMET_H2O = _COMET_H + _COMET_H + _COMET_O
_COMET_CO = _COMET_C + _COMET_O
_COMET_NH3 = _COMET_N + 3 * _COMET_H
_COMET_NH2 = _COMET_N + 2 * _COMET_H

# X! Tandem: masscalc::masscalc; m_dProton from msequtilities.cpp
_XT_H, _XT_O, _XT_N, _XT_SE, _XT_C, _XT_S = 1.007825035, 15.99491463, 14.003074, 79.916520, 12.0, 31.9720707
XTANDEM_PROTON = 1.007276
XTANDEM_RESIDUE_MASS = _residue_masses(_XT_C, _XT_H, _XT_N, _XT_O, _XT_S, _XT_SE)
_XT_H2O = 2 * _XT_H + _XT_O


@dataclass(frozen=True, slots=True)
class _Peptide:
    residues: list[str]
    residue_masses: list[float]  # engine residue mass plus the residue's modification deltas
    nterm_delta: float
    cterm_delta: float


def _peptide(peptidoform: str, residue_table: dict[str, float]) -> _Peptide | None:
    """Residues with the engine's masses and the peptidoform's modification deltas, or None if not understood."""
    try:
        parsed = parse_peptidoform(peptidoform)
    except Exception:
        return None
    if parsed.unsupported_modification or not parsed.residues:
        return None
    masses: list[float] = []
    for residue, deltas in zip(parsed.residues, parsed.residue_mod_deltas):
        base = residue_table.get(residue.upper())
        if base is None:
            return None
        masses.append(base + sum(deltas))
    return _Peptide([r.upper() for r in parsed.residues], masses, parsed.nterm_delta, parsed.cterm_delta)


# --- Comet XCorr ---------------------------------------------------------------------------------------------------

_COMET_ION_SERIES = ("a", "b", "c", "x", "y", "z", "z1")  # Comet's ION_SERIES_A..ION_SERIES_Z1 order


@dataclass(frozen=True, slots=True)
class CometXcorrParams:
    """The comet.params settings XCorr depends on; defaults are those of the file `comet -p` writes."""

    fragment_bin_tol: float = 0.02
    fragment_bin_offset: float = 0.0
    theoretical_fragment_ions: int = 0  # 0 = use flanking peaks, 1 = M peak only
    ion_series: tuple[str, ...] = ("b", "y")
    max_fragment_charge: int = 3
    remove_precursor_peak: int = 0  # 0 = no, 1 = precursor, 2 = precursor and charge-reduced precursors
    remove_precursor_tolerance: float = 1.5
    minimum_intensity: float = 0.0
    percentage_base_peak: float = 0.0
    peptide_mass_tolerance: float = 20.0  # peptide_mass_tolerance_upper; sets the binned array's size
    peptide_mass_units: int = 2  # 0 = amu, 1 = mmu, 2 = ppm
    precursor_tolerance_type: int = 1  # 0 = MH+, 1 = precursor m/z (only used for amu/mmu)

    def __post_init__(self) -> None:
        unknown = set(self.ion_series) - set(_COMET_ION_SERIES)
        if unknown:
            raise ValueError(f"Unknown Comet ion series {sorted(unknown)}; expected some of {', '.join(_COMET_ION_SERIES)}.")
        if self.remove_precursor_peak not in (0, 1, 2):
            raise ValueError("remove_precursor_peak must be 0, 1 or 2.")


# The two settings the comet.params template recommends.
COMET_PRESETS = {
    "high_res": CometXcorrParams(),
    "low_res": CometXcorrParams(fragment_bin_tol=1.0005, fragment_bin_offset=0.4, theoretical_fragment_ions=1),
}
_FLOAT_ZERO = 1e-6  # CometDataInternal.h; sparse-matrix entries within +-FLOAT_ZERO are not stored
_XCORR_OFFSET = 75  # iXcorrProcessingOffset


def _comet_bin(mass: float, params: CometXcorrParams) -> int:
    # Common.h: #define BIN(dMass) (int)((dMass)*dInverseBinWidth + dOneMinusBinOffset)
    return int(mass * (1.0 / params.fragment_bin_tol) + (1.0 - params.fragment_bin_offset))


def _comet_mass_cushion(exp_mass: float, params: CometXcorrParams) -> float:
    # CometPreprocess::GetMassCushion
    if params.peptide_mass_units == 0:
        cushion = params.peptide_mass_tolerance * (8 if params.precursor_tolerance_type == 1 else 1)
    elif params.peptide_mass_units == 1:
        cushion = params.peptide_mass_tolerance * 0.001 * (8 if params.precursor_tolerance_type == 1 else 1)
    else:
        cushion = params.peptide_mass_tolerance * exp_mass / 1e6
    return cushion + 5.0


def comet_fast_xcorr_spectrum(
    mz: list[float] | np.ndarray, intensity: list[float] | np.ndarray, precursor_mz: float, charge: int, params: CometXcorrParams
) -> tuple[np.ndarray, int]:
    """The preprocessed spectrum XCorr reads from (Comet's pfFastXcorrData after the sparse matrix), and its size."""
    exp_mass = charge * precursor_mz - (charge - 1.0) * COMET_PROTON  # dExpPepMass, the experimental [M+H]+
    array_size = int((exp_mass + _comet_mass_cushion(exp_mass, params)) * (1.0 / params.fragment_bin_tol))
    raw = np.zeros(array_size, dtype=np.float64)

    # LoadIons
    mz_values = np.asarray(mz, dtype=np.float64)
    intensities = np.asarray(intensity, dtype=np.float64)
    cutoff = params.minimum_intensity
    if 0.0 < params.percentage_base_peak <= 1.0 and intensities.size:
        cutoff = max(params.percentage_base_peak * float(intensities.max()), params.minimum_intensity)
    highest_ion = 0
    highest_intensity = 0.0
    precursors = [(exp_mass + (j - 1.0) * COMET_PROTON) / j for j in range(1, charge + 1)]
    for ion, value in zip(mz_values[::-1], intensities[::-1]):  # Comet reads peaks in reverse order
        if not (value >= cutoff and value > 0.0) or ion >= exp_mass + 50.0:
            continue
        bin_ = _comet_bin(float(ion), params)
        value = math.sqrt(float(value))
        highest_ion = max(highest_ion, bin_)
        if bin_ >= array_size or value <= raw[bin_]:
            continue
        if params.remove_precursor_peak == 1 and not abs(ion - precursors[-1]) > params.remove_precursor_tolerance:
            continue
        if params.remove_precursor_peak == 2 and any(abs(ion - p) < params.remove_precursor_tolerance for p in precursors):
            continue
        raw[bin_] = value
        highest_intensity = max(highest_intensity, value)

    # MakeCorrData: 10 windows, each normalised to 50; bins at or below 5% of the highest intensity are dropped
    corr = np.zeros(array_size, dtype=np.float64)
    window = highest_ion // 10 + 1
    for w in range(10):
        start, stop = w * window, min((w + 1) * window, highest_ion + 1, array_size)
        if start >= stop:
            continue
        segment = raw[start:stop]
        window_max = float(segment.max()) if segment.size else 0.0
        if window_max > 0.0:
            keep = segment > 0.05 * highest_intensity
            corr[start:stop] = np.where(keep, segment * (50.0 / window_max), 0.0)

    # Preprocess: y'_i = y_i - (sum of y over i-75..i+75, without y_i) / 150, truncated at the array's ends
    cumulative = np.concatenate(([0.0], np.cumsum(corr)))
    index = np.arange(array_size)
    lo = np.maximum(index - _XCORR_OFFSET, 0)
    hi = np.minimum(index + _XCORR_OFFSET, array_size - 1)
    fast = (cumulative[hi + 1] - cumulative[lo] - corr) * (1.0 / (2 * _XCORR_OFFSET))
    transformed = corr - fast
    data = transformed.astype(np.float32)
    data[0] = np.float32(0.0)
    if params.theoretical_fragment_ions == 0:
        # flanking peaks: pfFastXcorrData[i] += (float)(0.5*y'_{i-1}) + (float)(0.5*y'_{i+1}), in single precision
        left = np.zeros(array_size, dtype=np.float32)
        right = np.zeros(array_size, dtype=np.float32)
        left[1:] = (transformed[:-1] * 0.5).astype(np.float32)
        right[1:-1] = (transformed[2:] * 0.5).astype(np.float32)
        data[1:] = (data[1:] + left[1:]) + right[1:]
    data[np.abs(data) <= _FLOAT_ZERO] = 0.0
    return data, array_size


def comet_fragment_mz(peptide: _Peptide, ion_series: str, length: int, charge: int) -> float:
    """m/z of a Comet fragment ion of `length` residues (CometMassSpecUtils::GetFragmentIonMass)."""
    if ion_series in ("a", "b", "c"):
        mass = COMET_PROTON + peptide.nterm_delta + sum(peptide.residue_masses[:length])
        mass += {"a": -_COMET_CO, "b": 0.0, "c": _COMET_NH3}[ion_series]
    else:
        mass = COMET_PROTON + _COMET_H2O + peptide.cterm_delta + sum(peptide.residue_masses[len(peptide.residue_masses) - length :])
        mass += {"x": _COMET_CO - 2 * _COMET_H, "y": 0.0, "z": -_COMET_NH2, "z1": -_COMET_NH2 + _COMET_H}[ion_series]
    return (mass + (charge - 1.0) * COMET_PROTON) / charge


def comet_xcorr(
    observed: ObservedSpectrum,
    peptidoform: str,
    *,
    precursor_mz: float | None = None,
    charge: int | None = None,
    params: CometXcorrParams = COMET_PRESETS["high_res"],
) -> float | None:
    """Comet XCorr of the peptidoform against the observed spectrum, rounded to 3 decimals as Comet reports it.

    Returns None when the peptidoform cannot be expressed in Comet's terms (an unknown residue or modification) or
    the precursor m/z or charge is missing.
    """
    precursor_mz = precursor_mz if precursor_mz is not None else observed.precursor_mz
    charge = charge if charge is not None else observed.precursor_charge
    peptide = _peptide(peptidoform, COMET_RESIDUE_MASS)
    if peptide is None or precursor_mz is None or not charge or len(peptide.residues) < 2:
        return None
    data, array_size = comet_fast_xcorr_spectrum(observed.mz, observed.intensity, float(precursor_mz), int(charge), params)
    max_charge = 1 if charge == 1 else min(int(charge) - 1, params.max_fragment_charge)
    claimed: set[int] = set()
    xcorr = 0.0
    series = [s for s in _COMET_ION_SERIES if s in params.ion_series]
    for fragment_charge in range(1, max_charge + 1):
        for ion_series in series:
            for length in range(1, len(peptide.residues)):
                bin_ = _comet_bin(comet_fragment_mz(peptide, ion_series, length, fragment_charge), params)
                # Each bin counts once (pbDuplFragment); bins outside the spectrum's array score 0
                if bin_ <= 0 or bin_ in claimed:
                    continue
                claimed.add(bin_)
                if bin_ < array_size:
                    xcorr += float(data[bin_])
    xcorr *= 0.005  # "Scale intensities to 50 and divide score by 1E4."
    return round(xcorr * 1000.0) / 1000.0


# --- X! Tandem hyperscore -----------------------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class XTandemHyperscoreParams:
    """The X! Tandem input parameters the hyperscore depends on; defaults are those of default_input.xml."""

    fragment_monoisotopic_mass_error: float = 0.4
    fragment_monoisotopic_mass_error_units: str = "Daltons"  # or "ppm"
    dynamic_range: float = 100.0
    total_peaks: int = 50
    minimum_fragment_mz: float = 150.0
    minimum_peaks: int = 15
    # "refine, spectrum synthesis" (default yes) weights ions by the residues at the cleavage, in every pass
    spectrum_synthesis: bool = True
    # X! Tandem's "residue, modification mass": residues carrying only these keep their letter. Any other modification
    # is scored as X! Tandem scores a potential modification, on a lower-case residue that takes no synthesis weight.
    fixed_modifications: tuple[tuple[str, float], ...] = (("C", 57.021464),)

    def __post_init__(self) -> None:
        if self.fragment_monoisotopic_mass_error_units not in ("Daltons", "ppm"):
            raise ValueError("fragment_monoisotopic_mass_error_units must be 'Daltons' or 'ppm'.")


XTANDEM_DEFAULTS = XTandemHyperscoreParams()
_F32 = np.float32


def _xtandem_factorials() -> list[float]:
    # mscore_tandem::mscore_tandem builds 0!..63! by repeated double multiplication; hfactor caps at 63
    table, value = [1.0], 1.0
    for a in range(1, 64):
        value *= float(a)
        table.append(value)
    return table


_XTANDEM_FACTORIAL = _xtandem_factorials()


def xtandem_condition_spectrum(
    mz: list[float] | np.ndarray, intensity: list[float] | np.ndarray, parent_mh: float, charge: int, params: XTandemHyperscoreParams
) -> list[tuple[np.float32, np.float32]] | None:
    """X! Tandem's spectrum conditioning (mspectrumcondition::condition), on single-precision peaks.

    Returns the conditioned (m/z, intensity) peaks sorted by m/z, or None if X! Tandem rejects the spectrum.
    """
    if charge > 4:  # m_fMaxZ
        return None
    order = np.argsort(np.asarray(mz, dtype=np.float64), kind="stable")
    peaks = [[_F32(mz[i]), _F32(intensity[i])] for i in order]

    def merge_close(peaks: list[list[np.float32]], width: float) -> list[list[np.float32]]:
        # remove_isotopes / clean_isotopes: within `width` of the peak being written, keep the most intense; the
        # reference m/z is not updated when a later peak replaces it
        if len(peaks) < 2:
            return peaks
        kept: list[list[np.float32]] = []
        write = list(peaks[0])
        reference = write[0]
        for read in peaks[1:]:
            if (read[0] - reference) >= _F32(width) or read[0] < _F32(200.0):
                kept.append(write)
                write = list(read)
                reference = write[0]
            elif read[1] > write[1]:
                write = list(read)
        kept.append(write)
        return kept

    peaks = merge_close(peaks, 0.95)
    # remove_parent: within 2/z below or above the parent m/z (X! Tandem uses 1.00727 here)
    parent_mz = _F32(1.00727 + (parent_mh - 1.00727) / charge)
    window = _F32(_F32(2.0) / _F32(charge))
    peaks = [p for p in peaks if not ((parent_mz - p[0] >= 0.0 and parent_mz - p[0] < window) or (p[0] - parent_mz > 0.0 and p[0] - parent_mz < window))]
    # remove_low_masses: drop the leading peaks at or below the minimum fragment m/z
    first = next((i for i, p in enumerate(peaks) if p[0] > _F32(params.minimum_fragment_mz)), len(peaks))
    peaks = peaks[first:]
    # dynamic_range: scale so the base peak is dynamic_range, and drop peaks below 1
    if peaks:
        factor = _F32(max(p[1] for p in peaks) / _F32(params.dynamic_range))
        scaled = [[p[0], _F32(p[1] / factor)] for p in peaks]
        peaks = [p for p in scaled if p[1] >= _F32(1.0)]
    if len(peaks) < params.minimum_peaks:
        return None
    peaks = merge_close(peaks, 1.5)
    # remove_small: the total_peaks most intense (std::sort by intensity, descending), then back to m/z order
    if len(peaks) > params.total_peaks:
        peaks = sorted(peaks, key=lambda p: -float(p[1]))[: params.total_peaks]
    peaks.sort(key=lambda p: float(p[0]))
    return [(p[0], p[1]) for p in peaks]


def _xtandem_error(params: XTandemHyperscoreParams) -> float:
    # mscore::set_fragment_error: in ppm, the error is taken at m/z 200
    if params.fragment_monoisotopic_mass_error_units == "ppm":
        return 200.0 * params.fragment_monoisotopic_mass_error / 1e6
    return params.fragment_monoisotopic_mass_error


def _xtandem_blur(peaks: list[tuple[np.float32, np.float32]], params: XTandemHyperscoreParams) -> list[tuple[int, np.float32]]:
    """mscore_tandem::blur: each peak's intensity at its integer bin and its neighbours (more for ppm errors)."""
    width = 1.0
    convert = _xtandem_error(params) / width
    factor = 200.0 / convert
    w = int(0.1 + width)
    out: list[tuple[int, np.float32]] = []
    last_bin, last_intensity = 0, _F32(0.0)
    for mz, intensity in peaks:
        if not intensity > 0.5:
            continue
        value = int(float(mz) / convert)
        a = -w
        if params.fragment_monoisotopic_mass_error_units == "ppm":
            a = int(a * value / factor - 0.5)
            if a > -w:
                a = -w
        upper = -a
        while a <= upper:
            if last_bin == value + a:
                if last_intensity < intensity:
                    out[-1] = (out[-1][0], intensity)
            else:
                last_bin, last_intensity = value + a, intensity
                out.append((last_bin, intensity))
            a += 1
    return out


def _xtandem_dot(spectrum: list[tuple[int, np.float32]], ions: list[tuple[int, float]]) -> tuple[np.float32, int]:
    """mscore_tandem::dot: the merge walk over the blurred spectrum and the theoretical ions, as X! Tandem does it.

    Returns the convolution score (single precision) and the number of ions that landed on a spectrum entry.
    """
    bins = [b for b, _ in spectrum]
    seq = [b for b, _ in ions] + [0]
    size, count = len(bins), len(ions)
    step = max(int(0.5 + size / count), 1) if count else 1
    t = a = matched = 0
    score = _F32(0.0)
    lseq = seq[0]
    while lseq != 0 and t != size:
        advance_seq = False
        if bins[t] < lseq:
            advance_seq = True
            while t + step < size and bins[t + step] < lseq:
                t += step
            t += 1
            while t < size and bins[t] < lseq:
                t += 1
        elif bins[t] > lseq:
            a += 1
            lseq = seq[a]
            while bins[t] > lseq and lseq != 0:
                a += 1
                lseq = seq[a]
        if lseq == 0 or t == size:
            break
        if bins[t] == lseq:
            score = _F32(score + _F32(spectrum[t][1] * _F32(ions[a][1])))
            matched += 1
        if advance_seq:
            a += 1
            lseq = seq[a]
        else:
            t += 1
    return score, matched


# msequtilities::synthesis(true): b- and y-ion score factors by residue; every other residue, and every residue when
# synthesis is off, has factor 1
_XTANDEM_SYNTHESIS_B = {"D": 5.0, "N": 2.0, "V": 3.0, "E": 3.0, "Q": 2.0, "I": 3.0, "L": 3.0}
_XTANDEM_SYNTHESIS_Y = {"P": 5.0}


def _xtandem_letters(peptide: _Peptide, params: XTandemHyperscoreParams) -> list[str]:
    """The residues as X! Tandem holds them: lower case where a modification other than a fixed one is present."""
    fixed = {residue: mass for residue, mass in params.fixed_modifications}
    letters = []
    for residue, mass in zip(peptide.residues, peptide.residue_masses):
        delta = mass - XTANDEM_RESIDUE_MASS[residue]
        unmodified = abs(delta) < 1e-6 or (residue in fixed and abs(delta - fixed[residue]) < 1e-6)
        letters.append(residue if unmodified else residue.lower())
    return letters


def _xtandem_factor(table: dict[str, float], letter: str, params: XTandemHyperscoreParams) -> float:
    return table.get(letter, 1.0) if params.spectrum_synthesis else 1.0


def _xtandem_b_ions(peptide: _Peptide, letters: list[str], charge: int, we: float, params: XTandemHyperscoreParams) -> list[tuple[int, float]]:
    # mscore::add_B: weight = B factor of the residue before the cleavage x Y factor of the one after; b2 x3, or x10
    # when the second residue is an (unmodified) P
    value = peptide.nterm_delta
    ions = []
    for a in range(len(letters) - 1):
        value += peptide.residue_masses[a]
        weight = float(_F32(_xtandem_factor(_XTANDEM_SYNTHESIS_B, letters[a], params)) * _F32(_xtandem_factor(_XTANDEM_SYNTHESIS_Y, letters[a + 1], params)))
        if a == 1:
            weight = float(_F32(weight) * _F32(10 if letters[1] == "P" else 3))
        ions.append((int((XTANDEM_PROTON + value / charge) * we), weight))
    return ions


def _xtandem_y_ions(
    peptide: _Peptide, letters: list[str], charge: int, we: float, params: XTandemHyperscoreParams, *, largest_only: bool = False
) -> list[tuple[int, float]]:
    # mscore::add_Y: weight = Y factor of the residue after the cleavage x B factor of the one before; y(n-2) x3, or
    # x10 when the second residue is an (unmodified) P. largest_only is add_Y(0, 2): only a < 5, without that factor
    value = _XT_H2O + peptide.cterm_delta
    ions = []
    for a in range(len(letters) - 1, 0, -1):
        value += peptide.residue_masses[a]
        mz_bin = int((XTANDEM_PROTON + value / charge) * we)
        weight = float(_F32(_xtandem_factor(_XTANDEM_SYNTHESIS_Y, letters[a], params)) * _F32(_xtandem_factor(_XTANDEM_SYNTHESIS_B, letters[a - 1], params)))
        if largest_only:
            if a < 5:
                ions.append((mz_bin, weight))
            continue
        if a == 2:
            weight = float(_F32(weight) * _F32(10 if letters[1] == "P" else 3))
        ions.append((mz_bin, weight))
    return ions


def xtandem_hyperscore(
    observed: ObservedSpectrum,
    peptidoform: str,
    *,
    precursor_mz: float | None = None,
    charge: int | None = None,
    params: XTandemHyperscoreParams = XTANDEM_DEFAULTS,
) -> float | None:
    """The X! Tandem hyperscore as X! Tandem reports it: 4·log10(H) to one decimal (mscore_tandem::hconvert).

    H = (sum over b and y ions of observed intensity x ion weight) x the product of the two largest per-charge
    factors N_b! N_y!, scored on X! Tandem's conditioned spectrum (mscore::score). Returns None when X! Tandem would
    reject the spectrum, or the peptidoform or precursor cannot be used.
    """
    raw = xtandem_hyperscore_raw(observed, peptidoform, precursor_mz=precursor_mz, charge=charge, params=params)
    if raw is None:
        return None
    reported = _F32(0.0) if raw <= 0.0 else _F32(_F32(4.0) * _F32(math.log10(float(raw))))
    return float(f"{float(reported):.1f}")


def xtandem_hyperscore_raw(
    observed: ObservedSpectrum,
    peptidoform: str,
    *,
    precursor_mz: float | None = None,
    charge: int | None = None,
    params: XTandemHyperscoreParams = XTANDEM_DEFAULTS,
) -> np.float32 | None:
    """H before X! Tandem's 4·log10 conversion (m_fHyper)."""
    precursor_mz = precursor_mz if precursor_mz is not None else observed.precursor_mz
    charge = charge if charge is not None else observed.precursor_charge
    peptide = _peptide(peptidoform, XTANDEM_RESIDUE_MASS)
    if peptide is None or precursor_mz is None or not charge or len(peptide.residues) < 2:
        return None
    charge = int(charge)
    parent_mh = (float(precursor_mz) - XTANDEM_PROTON) * charge + XTANDEM_PROTON  # loadmspectrum (MGF)
    peaks = xtandem_condition_spectrum(observed.mz, observed.intensity, parent_mh, charge, params)
    if peaks is None:
        return None
    spectrum = _xtandem_blur(peaks, params)
    letters = _xtandem_letters(peptide, params)
    we = 1.0 / _xtandem_error(params)
    charge_limit = 2 if charge == 1 else charge
    factors = [1.0] * (charge + 1)
    total = 0.0
    for ion_type in ("y", "b"):  # mscore::score walks T_Y (0x01) before T_B (0x02)
        value = 0.0
        for fragment_charge in range(1, charge_limit):
            if ion_type == "y":
                ions = _xtandem_y_ions(peptide, letters, fragment_charge, we, params)
            else:
                ions = _xtandem_b_ions(peptide, letters, fragment_charge, we, params)
            score, matched = _xtandem_dot(spectrum, ions)
            value += float(score)
            if fragment_charge == 1 and ion_type == "y" and charge == 2:
                extra_score, extra_matched = _xtandem_dot(spectrum, _xtandem_y_ions(peptide, letters, 2, we, params, largest_only=True))
                value += float(extra_score)
                matched += extra_matched
            factors[fragment_charge] *= _XTANDEM_FACTORIAL[min(matched, 63)]
        total += value
    factors.sort()
    best = factors[-1] * factors[-2] if len(factors) > 2 else factors[-1]
    hyper = best * total
    return _F32(min(hyper, float(np.finfo(np.float32).max)))
