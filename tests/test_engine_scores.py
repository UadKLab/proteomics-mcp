"""Comet XCorr and X! Tandem hyperscore (#51): each preprocessing step, and the engines' own scores on the examples."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from proteomics_mcp.core import engine_scores as es
from proteomics_mcp.models.schemas import ObservedSpectrum

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = json.loads((ROOT / "tests" / "fixtures" / "engine_reference_examples.json").read_text(encoding="utf-8"))


def _example_spectrum(name: str) -> ObservedSpectrum:
    """The payload's spectrum exactly as the engines read it from the MGF: m/z to 6, intensity to 4 decimals."""
    payload = json.loads((ROOT / "examples" / f"{name}.json").read_text(encoding="utf-8"))
    observed, candidate = payload["observed"], payload["metrics"]["candidate"]
    pairs = [(float(f"{mz:.6f}"), float(f"{value:.4f}")) for mz, value in zip(observed["mz"], observed["intensity"]) if value > 0]
    precursor_mz = candidate.get("precursor_mz") or observed.get("precursor_mz")
    return ObservedSpectrum(mz=[mz for mz, _ in pairs], intensity=[value for _, value in pairs], precursor_mz=float(precursor_mz))


# --- engines' own scores ---------------------------------------------------------------------------------------------


@pytest.mark.parametrize("row", REFERENCE["comet"], ids=lambda row: f"{row['payload'].split('/')[-1][:30]}-{row['rank']}")
def test_comet_xcorr_matches_comet(row):
    params = es.CometXcorrParams(**REFERENCE["comet_params"])
    observed = _example_spectrum(row["payload"])

    assert es.comet_xcorr(observed, row["peptidoform"], charge=row["charge"], params=params) == pytest.approx(row["xcorr"], abs=1e-9)


@pytest.mark.parametrize("row", REFERENCE["xtandem"], ids=lambda row: row["payload"].split("/")[-1][:40])
def test_xtandem_hyperscore_matches_xtandem(row):
    observed = _example_spectrum(row["payload"])

    assert es.xtandem_hyperscore(observed, row["peptidoform"], charge=row["charge"]) == pytest.approx(row["hyperscore"], abs=1e-9)


def test_reference_covers_decoys_and_lower_ranks():
    # Agreement on wrong candidates matters as much as on the best one
    assert sum(row["decoy"] for row in REFERENCE["comet"]) > 50
    assert {row["rank"] for row in REFERENCE["comet"]} == set(range(1, 11))


# --- Comet steps -----------------------------------------------------------------------------------------------------


def test_comet_bin_follows_the_bin_macro():
    high, low = es.COMET_PRESETS["high_res"], es.COMET_PRESETS["low_res"]
    assert es._comet_bin(100.0, high) == 5001  # (int)(100/0.02 + 1.0)
    assert es._comet_bin(100.0, low) == int(100.0 / 1.0005 + 0.6)


def test_fast_xcorr_subtracts_the_mean_of_the_150_neighbouring_bins_and_adds_flanking_peaks():
    params = es.CometXcorrParams(peptide_mass_tolerance=3.0, peptide_mass_units=0, precursor_tolerance_type=0)
    data, size = es.comet_fast_xcorr_spectrum([500.0], [400.0], 600.0, 1, params)
    peak = es._comet_bin(500.0, params)

    # One peak: sqrt(400) = 20 is the window maximum, normalised to 50. y' = 50 at the peak and -50/150 within 75 bins.
    assert size == int((600.0 + 3.0 + 5.0) / 0.02)
    assert data[peak] == pytest.approx(50 - 50 / 150, abs=1e-5)  # 50 + 0.5 * (-1/3) on both sides
    assert data[peak + 1] == pytest.approx(-50 / 150 + 0.5 * 50 - 0.5 * 50 / 150, abs=1e-5)
    assert data[peak + 40] == pytest.approx(-50 / 150 * 2, abs=1e-5)
    assert data[peak + 80] == 0.0
    assert data.dtype == np.float32


def test_comet_drops_peaks_at_or_below_five_percent_of_the_base_peak():
    params = es.CometXcorrParams(theoretical_fragment_ions=1)
    data, _ = es.comet_fast_xcorr_spectrum([300.0, 900.0], [10000.0, 24.0], 600.0, 2, params)  # sqrt: 100 and 4.9
    assert data[es._comet_bin(900.0, params)] == 0.0  # dropped, and no other peak within 75 bins


def test_comet_xcorr_sums_each_theoretical_bin_once():
    # Charge 1 and 2 fragments of a 3+ precursor; a bin claimed by an earlier ion (charge first, then a/b/c/x/y/z,
    # then length) is not counted again (CometSearch.cpp, pbDuplFragment)
    params = es.COMET_PRESETS["low_res"]  # wide bins, so b and y ions do share bins
    peptide = es._peptide("GGGGGGGGGR", es.COMET_RESIDUE_MASS)
    observed = _example_spectrum(REFERENCE["comet"][0]["payload"])
    data, size = es.comet_fast_xcorr_spectrum(observed.mz, observed.intensity, observed.precursor_mz, 3, params)
    bins = [es._comet_bin(es.comet_fragment_mz(peptide, series, length, z), params) for z in (1, 2) for series in ("b", "y") for length in range(1, 10)]
    assert len(set(bins)) < len(bins)
    expected = round(0.005 * sum(float(data[b]) for b in dict.fromkeys(bins) if 0 < b < size) * 1000) / 1000
    assert es.comet_xcorr(observed, "GGGGGGGGGR", charge=3, params=params) == pytest.approx(expected, abs=1e-12)


def test_comet_xcorr_is_none_for_residues_comet_does_not_know():
    observed = ObservedSpectrum(mz=[200.0], intensity=[1.0], precursor_mz=400.0, precursor_charge=2)
    assert es.comet_xcorr(observed, "PEPTIBE") is None
    without_charge = ObservedSpectrum(mz=[200.0], intensity=[1.0], precursor_mz=400.0)
    assert es.comet_xcorr(without_charge, "PEPTIDE") is None


# --- X! Tandem steps -------------------------------------------------------------------------------------------------


def _conditioned(mz, intensity, parent_mh=2000.0, charge=2, **params):
    return es.xtandem_condition_spectrum(mz, intensity, parent_mh, charge, es.XTandemHyperscoreParams(minimum_peaks=0, **params))


def test_xtandem_keeps_the_most_intense_peak_within_each_isotope_window():
    # remove_isotopes (0.95 Da from 400.0) keeps 400.5; clean_isotopes (1.5 Da) keeps 402.0 apart from it
    peaks = _conditioned([400.0, 400.5, 400.9, 402.0], [10.0, 50.0, 20.0, 30.0])
    assert [round(float(mz), 1) for mz, _ in peaks] == [400.5, 402.0]
    # within 1.5 Da the second pass keeps only the more intense of the two survivors
    assert [round(float(mz), 1) for mz, _ in _conditioned([400.0, 400.5, 401.2], [10.0, 50.0, 30.0])] == [400.5]


def test_xtandem_removes_the_parent_region_and_low_masses():
    parent_mz = 1.00727 + (2000.0 - 1.00727) / 2
    peaks = _conditioned([140.0, 300.0, parent_mz - 0.5, parent_mz + 0.5, parent_mz + 1.5], [5.0, 5.0, 5.0, 5.0, 5.0])
    assert [round(float(mz), 2) for mz, _ in peaks] == [300.0, round(parent_mz + 1.5, 2)]


def test_xtandem_dynamic_range_drops_peaks_below_one_percent_and_keeps_the_top_50():
    mz = [200.0 + 2 * i for i in range(80)]
    intensity = [1000.0] + [5.0] + [100.0 + i for i in range(78)]
    peaks = _conditioned(mz, intensity)
    assert len(peaks) == 50
    assert all(float(value) >= 1.0 for _, value in peaks)
    assert float(max(value for _, value in peaks)) == pytest.approx(100.0)


def test_xtandem_rejects_spectra_with_too_few_peaks():
    params = es.XTandemHyperscoreParams()  # minimum peaks 15
    assert es.xtandem_condition_spectrum([300.0, 400.0], [1.0, 1.0], 1000.0, 2, params) is None


def test_synthesis_weights_follow_the_residues_at_the_cleavage():
    params = es.XTANDEM_DEFAULTS
    peptide = es._peptide("ADPLK", es.XTANDEM_RESIDUE_MASS)
    letters = es._xtandem_letters(peptide, params)
    b_weights = [weight for _, weight in es._xtandem_b_ions(peptide, letters, 1, 1 / 0.4, params)]
    y_weights = [weight for _, weight in es._xtandem_y_ions(peptide, letters, 1, 1 / 0.4, params)]
    # b1 A|D: 1; b2 D|P: 5 x 5 x 3 (b2); b3 P|L: 1; b4 L|K: 3
    assert b_weights == [1.0, 75.0, 1.0, 3.0]
    # y1 K: L|K 3; y2 LK: P|L 1; y3 PLK: D|P 25 x 3 (y(n-2)); y4: A|D 1
    assert y_weights == [3.0, 1.0, 75.0, 1.0]
    no_synthesis = es.XTandemHyperscoreParams(spectrum_synthesis=False)
    assert [w for _, w in es._xtandem_b_ions(peptide, es._xtandem_letters(peptide, no_synthesis), 1, 2.5, no_synthesis)] == [1.0, 3.0, 1.0, 1.0]


def test_a_potential_modification_removes_the_residue_weight_but_a_fixed_one_does_not():
    params = es.XTANDEM_DEFAULTS
    assert es._xtandem_letters(es._peptide("DC[+57.021464]K", es.XTANDEM_RESIDUE_MASS), params) == ["D", "C", "K"]
    assert es._xtandem_letters(es._peptide("D[+0.984016]CK", es.XTANDEM_RESIDUE_MASS), params) == ["d", "C", "K"]


def test_doubly_charged_precursors_add_the_four_largest_y_ions_at_charge_2():
    peptide = es._peptide("PEPTIDEK", es.XTANDEM_RESIDUE_MASS)
    letters = es._xtandem_letters(peptide, es.XTANDEM_DEFAULTS)
    extra = es._xtandem_y_ions(peptide, letters, 2, 2.5, es.XTANDEM_DEFAULTS, largest_only=True)
    full = es._xtandem_y_ions(peptide, letters, 2, 2.5, es.XTANDEM_DEFAULTS)
    assert [b for b, _ in extra] == [b for b, _ in full][-4:]


def test_hyperscore_is_reported_as_4_log10_to_one_decimal():
    observed = _example_spectrum(REFERENCE["xtandem"][0]["payload"])
    row = REFERENCE["xtandem"][0]
    raw = es.xtandem_hyperscore_raw(observed, row["peptidoform"], charge=row["charge"])
    assert es.xtandem_hyperscore(observed, row["peptidoform"], charge=row["charge"]) == round(4 * np.log10(float(raw)), 1)
