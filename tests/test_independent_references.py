"""Checks against references the package does not produce itself (#179, from the #33 checklist).

- fragment and peptide masses against pyteomics, with modifications given as elemental formulas;
- peptidoform masses against the masses MaxQuant, FragPipe and Proteome Discoverer report for the reference PSMs, and
  precursor mass errors against the spectra's own precursor m/z and the engines' own errors;
- invariants of scoring and annotation on random spectra: explained intensity at most 1, each observed peak explains
  at most one predicted ion and each predicted ion at most one peak.

The whole supplement tables are checked the same way in tests/e2e/test_batch1_tables.py.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from proteomics_mcp.adapters.fragpipe import load_fragpipe_psms
from proteomics_mcp.adapters.maxquant import load_maxquant_psms
from proteomics_mcp.adapters.proteome_discoverer import load_proteome_discoverer_psms
from proteomics_mcp.core.annotation import PROTON_MASS, annotate_observed_ions, generate_terminal_fragments, peptidoform_neutral_mass
from proteomics_mcp.core.scoring import precursor_isotope_corrected_errors, score_candidate_spectrum
from proteomics_mcp.models.schemas import CandidatePSM, ObservedSpectrum, PredictedFragment, PredictedSpectrum

mass = pytest.importorskip("pyteomics.mass")

EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
NEUTRON = 1.00335483507
# The built-in modification masses are Unimod's, printed to 6 decimals: within 0.5 µDa of the formulas' masses
MASS_TOLERANCE = 1e-5

# Modifications as elemental formulas: their masses come from pyteomics' element table, not from the package's
FORMULAS = {
    "UNIMOD:1": "H2C2O",  # Acetyl
    "UNIMOD:2": "HNO-1",  # Amidated
    "UNIMOD:4": "H3C2NO",  # Carbamidomethyl
    "UNIMOD:7": "H-1N-1O",  # Deamidated
    "UNIMOD:21": "HPO3",  # Phospho
    "UNIMOD:35": "O",  # Oxidation
    "UNIMOD:121": "H6C4N2O2",  # GG
    "UNIMOD:737": "H20C8C[13]4NN[15]O2",  # TMT6plex
}

# (sequence, {1-based position: accession}, N-terminal accession, C-terminal accession)
PEPTIDOFORMS = [
    ("PEPTIDEK", {}, None, None),
    ("ACDEFGHIK", {2: "UNIMOD:4"}, "UNIMOD:1", None),
    ("PEPMSTIDEK", {4: "UNIMOD:35", 5: "UNIMOD:21"}, None, None),
    ("SAMPLEKR", {7: "UNIMOD:737"}, "UNIMOD:737", None),
    ("NQGHTYLLK", {2: "UNIMOD:7", 5: "UNIMOD:21", 6: "UNIMOD:21"}, None, None),
    ("LGGKESTAR", {4: "UNIMOD:121"}, None, None),
    ("PEPTIDEK", {}, None, "UNIMOD:2"),
]


def _proforma(sequence: str, mods: dict[int, str], nterm: str | None, cterm: str | None) -> str:
    body = "".join(f"{residue}[{mods[i]}]" if i in mods else residue for i, residue in enumerate(sequence, start=1))
    return f"{f'[{nterm}]-' if nterm else ''}{body}{f'-[{cterm}]' if cterm else ''}"


def _delta(accession: str | None) -> float:
    return mass.calculate_mass(formula=FORMULAS[accession]) if accession else 0.0


@pytest.mark.parametrize(("sequence", "mods", "nterm", "cterm"), PEPTIDOFORMS, ids=[_proforma(*p) for p in PEPTIDOFORMS])
def test_fragment_masses_match_pyteomics(sequence, mods, nterm, cterm):
    fragments = generate_terminal_fragments(_proforma(sequence, mods, nterm, cterm), ion_types="a,b,c,x,y,z,z-dot", max_ion_charge=3)
    length = len(sequence)

    assert len(fragments) == (length - 1) * 7 * 3
    for fragment in fragments:
        ordinal, charge, ion_type = fragment["ordinal"], fragment["charge"], fragment["ion_type"]
        if ion_type in {"a", "b", "c"}:
            residues, delta = sequence[:ordinal], _delta(nterm) + sum(_delta(mods.get(i)) for i in range(1, ordinal + 1))
        else:
            residues, delta = sequence[-ordinal:], _delta(cterm) + sum(_delta(mods.get(i)) for i in range(length - ordinal + 1, length + 1))
        expected = mass.fast_mass(residues, ion_type=ion_type, charge=charge) + delta / charge
        assert fragment["mz"] == pytest.approx(expected, abs=MASS_TOLERANCE), fragment["annotation"]


@pytest.mark.parametrize(("sequence", "mods", "nterm", "cterm"), PEPTIDOFORMS, ids=[_proforma(*p) for p in PEPTIDOFORMS])
def test_peptide_masses_match_pyteomics(sequence, mods, nterm, cterm):
    expected = mass.fast_mass(sequence) + _delta(nterm) + _delta(cterm) + sum(_delta(accession) for accession in mods.values())

    assert peptidoform_neutral_mass(_proforma(sequence, mods, nterm, cterm)) == pytest.approx(expected, abs=MASS_TOLERANCE)


# The search-engine reference PSMs: each payload keeps its row of the engine's table, which is read again here with
# today's adapters. (adapter, file name, options, engine mass column, offset to neutral mass, mass tolerance in Da: the
# engines print 4, 3 and 5 decimals)
ENGINES = {
    "17_maxquant_top5": (load_maxquant_psms, "msms.txt", {"fixed_modifications": "Carbamidomethyl (C)"}, "Mass", 0.0, 1e-4),
    "18_fragpipe_top5": (load_fragpipe_psms, "psm.tsv", {}, "Calculated Peptide Mass", 0.0, 6e-4),
    "19_proteome_discoverer_top5": (load_proteome_discoverer_psms, "psms.tsv", {}, "Theo. MH+ [Da]", -PROTON_MASS, 1e-4),
}
ENGINE_PAYLOADS = sorted(path for folder in ENGINES for path in EXAMPLES.glob(f"{folder}/*_payload.json"))


def _reread(payload_path: Path, tmp_path: Path) -> tuple[CandidatePSM, dict]:
    loader, name, options, *_ = ENGINES[payload_path.parent.name]
    payload = json.loads(payload_path.read_text(encoding="utf-8"))
    table = tmp_path / name
    pd.DataFrame([payload["metrics"]["candidate"]["source_columns"]]).to_csv(table, sep="\t", index=False)
    [candidate] = loader(table, **options)
    return candidate, payload


@pytest.mark.skipif(not ENGINE_PAYLOADS, reason="examples/ (reference payloads) is only available in a repository checkout")
@pytest.mark.parametrize("payload_path", ENGINE_PAYLOADS, ids=[f"{p.parent.name[:2]}:{p.name.split('_')[3]}" for p in ENGINE_PAYLOADS])
def test_peptidoform_mass_matches_the_engine(payload_path, tmp_path):
    candidate, payload = _reread(payload_path, tmp_path)
    *_, mass_column, offset, tolerance = ENGINES[payload_path.parent.name]

    reported = float(payload["metrics"]["candidate"]["source_columns"][mass_column]) + offset
    assert peptidoform_neutral_mass(candidate.peptidoform) == pytest.approx(reported, abs=tolerance)
    assert candidate.warnings == []


@pytest.mark.skipif(not ENGINE_PAYLOADS, reason="examples/ (reference payloads) is only available in a repository checkout")
@pytest.mark.parametrize("payload_path", ENGINE_PAYLOADS, ids=[f"{p.parent.name[:2]}:{p.name.split('_')[3]}" for p in ENGINE_PAYLOADS])
def test_precursor_error_of_accepted_psms_is_a_few_ppm(payload_path, tmp_path):
    candidate, payload = _reread(payload_path, tmp_path)
    row = payload["metrics"]["candidate"]["source_columns"]

    # against the precursor m/z in the spectrum file, which the package did not compute
    isotope, ppm, _ = precursor_isotope_corrected_errors(candidate.peptidoform, candidate.precursor_charge, payload["observed"]["precursor_mz"])
    assert isotope == 0 and abs(ppm) < 5

    # against the engine's own error for the precursor m/z it reports
    _, ours, _ = precursor_isotope_corrected_errors(candidate.peptidoform, candidate.precursor_charge, candidate.precursor_mz)
    if payload_path.parent.name.startswith("17_"):
        engine = float(row["Mass Error [ppm]"])
    elif payload_path.parent.name.startswith("18_"):
        calculated, delta = float(row["Calculated Peptide Mass"]), float(row["Delta Mass"])
        engine = (float(row["Calibrated Observed Mass"]) - calculated - round(delta / NEUTRON) * NEUTRON) / calculated * 1e6
    else:
        engine = float(row["DeltaM [ppm]"])
    # MaxQuant's m/z and FragPipe's masses are printed with 5 and 4 decimals
    assert ours == pytest.approx(engine, abs=0.1)


# Invariants on random spectra: random peptides, predictions with weak and zero ions, observed peaks that are noise,
# near-duplicates within the tolerance of one ion, and peaks between two close ions


def _random_case(seed: int, tolerance: float, unit: str):
    rng = np.random.default_rng(seed)
    sequence = "".join(rng.choice(list("ACDEFGHIKLMNPQRSTVWY"), size=int(rng.integers(6, 22))))
    peptidoform = "".join(f"{r}[UNIMOD:35]" if r == "M" and rng.random() < 0.5 else f"{r}[UNIMOD:21]" if r in "STY" and rng.random() < 0.2 else r for r in sequence)
    ions = generate_terminal_fragments(peptidoform, ion_types="b,y", max_ion_charge=2)
    # some ions absent (0) or near-absent (below the 0.001 cut-off), as predictors give them
    intensities = rng.random(len(ions))
    intensities[rng.random(len(ions)) < 0.15] = 0.0
    intensities[rng.random(len(ions)) < 0.15] = 1e-6
    fragments = [
        PredictedFragment(mz=ion["mz"], intensity=float(value), annotation=ion["annotation"], ion_type=ion["ion_type"], ordinal=ion["ordinal"], charge=ion["charge"])
        for ion, value in zip(ions, intensities)
    ]
    width = tolerance if unit == "Da" else tolerance * 1e-6 * 1000.0
    picked = [ion["mz"] for ion in ions if rng.random() < 0.6]
    observed_mz = (
        [mz + rng.uniform(-0.9, 0.9) * width for mz in picked]
        + [mz + rng.uniform(-0.9, 0.9) * width for mz in picked if rng.random() < 0.3]  # two peaks within one ion's tolerance
        + [(a + b) / 2 for a, b in zip(sorted(picked), sorted(picked)[1:]) if b - a < width]  # a peak between two ions
        + list(rng.uniform(100, 2000, size=int(rng.integers(0, 60))))  # noise
    )
    observed_mz = sorted(observed_mz)
    observed = ObservedSpectrum(mz=observed_mz, intensity=list(rng.uniform(1, 1e6, size=len(observed_mz))), precursor_charge=3)
    candidate = CandidatePSM(candidate_id=str(seed), source_engine="test", peptidoform=peptidoform, precursor_charge=3)
    predicted = PredictedSpectrum(peptidoform=peptidoform, precursor_charge=3, predictor_name="random", fragments=fragments)
    return candidate, observed, predicted


@pytest.mark.parametrize(("tolerance", "unit"), [(0.05, "Da"), (20.0, "ppm"), (0.5, "Da")])
@pytest.mark.parametrize("seed", range(40))
def test_scoring_invariants(seed, tolerance, unit):
    candidate, observed, predicted = _random_case(seed, tolerance, unit)

    result = score_candidate_spectrum(candidate=candidate, observed=observed, predicted=predicted, tolerance=tolerance, tolerance_unit=unit,
                                      metric_profile="v2", use_observed_annotation=False)

    assert 0.0 <= result.explained_intensity_fraction <= 1.0 + 1e-12
    assert result.matched_ion_count <= min(result.total_predicted_ions, len(observed.mz))
    observed_indices = [match.observed_index for match in result.matched_peaks]
    predicted_indices = [match.predicted_index for match in result.matched_peaks]
    assert len(set(observed_indices)) == len(observed_indices), "an observed peak explains more than one predicted ion"
    assert len(set(predicted_indices)) == len(predicted_indices), "a predicted ion is matched to more than one peak"
    assert all(abs(match.mz_error) <= tolerance + 1e-9 for match in result.matched_peaks)
    for metric in (result.spectral_angle, result.cosine_similarity):
        assert metric is None or 0.0 <= metric <= 1.0 + 1e-12
    assert result.pearson_correlation is None or -1.0 - 1e-12 <= result.pearson_correlation <= 1.0 + 1e-12


@pytest.mark.parametrize("seed", range(40))
def test_annotation_invariants(seed):
    candidate, observed, _ = _random_case(seed, 0.05, "Da")

    annotation = annotate_observed_ions(observed, peptidoform=candidate.peptidoform, tolerance=0.05, ion_types="all", max_isotope=2,
                                        neutral_losses=True)

    annotated = {match["observed_index"] for match in annotation["matches"]}
    assert annotation["annotated_peak_count"] == len(annotated) <= len(observed.mz)
    assert annotation["annotated_intensity"] == pytest.approx(sum(observed.intensity[i] for i in annotated))
    assert 0.0 <= annotation["annotated_intensity_fraction"] <= 1.0 + 1e-12
    assert all(math.isfinite(match["theoretical_mz"]) and abs(match["mz_error"]) <= 0.05 + 1e-9 for match in annotation["matches"])
