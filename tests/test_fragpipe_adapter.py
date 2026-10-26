from __future__ import annotations

import pandas as pd
import pytest

from proteomics_mcp.adapters import fragpipe
from proteomics_mcp.adapters.fragpipe import load_fragpipe_psms
from proteomics_mcp.core.annotation import AA_RESIDUE_MASS


def test_load_fragpipe_psms_extracts_join_fields(tmp_path):
    path = tmp_path / "psm.tsv"
    path.write_text(
        "\t".join(
            [
                "Spectrum",
                "Peptide",
                "Modified Peptide",
                "Charge",
                "Observed M/Z",
                "Calibrated Observed M/Z",
                "Hyperscore",
                "Nextscore",
                "Probability",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "runA.00582.00582.2",
                "KATEVSK",
                "",
                "2",
                "381.7211",
                "381.7226",
                "11.278",
                "9.500",
                "0.9840",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    [candidate] = load_fragpipe_psms(path)

    assert candidate.source_engine == "fragpipe"
    assert candidate.experiment_name == "runA"
    assert candidate.scan_number == 582
    assert candidate.spectrum_id == "runA.00582.00582.2"
    assert candidate.peptidoform == "KATEVSK"
    assert candidate.precursor_charge == 2
    assert candidate.precursor_mz == 381.7226
    assert candidate.raw_score == 11.278
    assert candidate.raw_score_name == "Hyperscore"
    assert candidate.next_beam is None
    assert candidate.next_beam_margin is None
    assert candidate.source_columns["Nextscore"] == 9.5
    assert candidate.confidence == 0.984


def test_load_fragpipe_psms_prefers_modified_peptide(tmp_path):
    path = tmp_path / "psm.tsv"
    path.write_text(
        "\t".join(["Spectrum", "Peptide", "Modified Peptide", "Charge", "Observed M/Z"])
        + "\n"
        + "\t".join(["runA.00042.00042.2", "MPEPTIDE", "M[15.994915]PEPTIDE", "2", "500.0"])
        + "\n",
        encoding="utf-8",
    )

    [candidate] = load_fragpipe_psms(path)

    assert candidate.peptidoform == "M[UNIMOD:35]PEPTIDE"


def test_load_fragpipe_psms_prefers_assigned_modifications_for_numeric_placeholders(tmp_path):
    path = tmp_path / "psm.tsv"
    path.write_text(
        "\t".join(["Spectrum", "Peptide", "Modified Peptide", "Assigned Modifications", "Charge", "Observed M/Z"])
        + "\n"
        + "\t".join(["runA.00042.00042.2", "TGMTMTEK", "TGM[147]TM[147]TEK", "3M(15.9949),5M(15.9949)", "2", "500.0"])
        + "\n",
        encoding="utf-8",
    )

    [candidate] = load_fragpipe_psms(path)

    assert candidate.peptidoform == "TGM[UNIMOD:35]TM[UNIMOD:35]TEK"


# Spectrum names, "Modified Peptide" masses and the fallbacks between FragPipe's peptide columns (#166)


@pytest.mark.parametrize(
    ("spectrum", "expected"),
    [
        ("run_01.01234.01234.3", ("run_01", 1234, 3)),
        ("run.with.dots.00042.00042.2", ("run.with.dots", 42, 2)),
        ("run_01:57", ("run_01", 57, None)),
        ("run_01.57", ("run_01", 57, None)),
        ("not a spectrum name", ("not a spectrum name", None, None)),
        (None, (None, None, None)),
        ("  ", (None, None, None)),
    ],
)
def test_spectrum_parts(spectrum, expected):
    assert fragpipe._spectrum_parts(spectrum) == expected


@pytest.mark.parametrize(
    ("value", "residue", "expected"),
    [
        ("147", "M", "UNIMOD:35"),  # rounded total: M + oxidation
        ("147.0354", "M", "UNIMOD:35"),  # exact total
        ("15.994915", "M", "UNIMOD:35"),  # a delta, smaller than the residue
        ("160", "C", "UNIMOD:4"),  # rounded total: C + carbamidomethyl
    ],
)
def test_delta_token_reads_totals_and_deltas(value, residue, expected):
    assert fragpipe._delta_token(value, AA_RESIDUE_MASS[residue]) == expected


def test_delta_token_looks_up_unknown_exact_masses_in_the_catalogue(monkeypatch):
    looked_up = []
    monkeypatch.setattr(fragpipe, "resolve_mod_token", lambda token: looked_up.append(token) or "UNIMOD:9999")

    assert fragpipe._delta_token("100.123456", AA_RESIDUE_MASS["G"]) == "UNIMOD:9999"
    assert looked_up == [f"{100.123456 - AA_RESIDUE_MASS['G']:+.6f}"]


def test_delta_token_keeps_a_mass_nothing_matches(monkeypatch):
    monkeypatch.setattr(fragpipe, "resolve_mod_token", lambda token: None)

    assert fragpipe._delta_token("100.123456", AA_RESIDUE_MASS["G"]) == f"{100.123456 - AA_RESIDUE_MASS['G']:+.6f}"
    # a rounded mass is never looked up in the catalogue, which has an entry near almost any mass
    monkeypatch.setattr(fragpipe, "resolve_mod_token", lambda token: pytest.fail("rounded masses are not looked up"))
    assert fragpipe._delta_token("160", AA_RESIDUE_MASS["G"]) == f"{160 - AA_RESIDUE_MASS['G']:+.6f}"


def test_modified_peptide_in_another_notation_is_kept_visible():
    # signed masses are not FragPipe's total-mass notation; every token stays in the peptidoform
    assert fragpipe._peptidoform_from_modified_peptide("PEPM[+15.9949]K") == "PEPM[UNIMOD:35]K"


def test_modified_peptide_is_used_when_assigned_modifications_are_not_understood():
    row = pd.Series({"Peptide": "PEPMK", "Assigned Modifications": "4M(15.9949), 9X(1.0)", "Modified Peptide": "PEPM[147]K"})

    assert fragpipe._peptidoform_from_row(row) == "PEPM[UNIMOD:35]K"


def test_a_nan_modified_peptide_is_ignored():
    row = pd.Series({"Peptide": "PEPTIDEK", "Assigned Modifications": "", "Modified Peptide": "nan"})

    assert fragpipe._peptidoform_from_row(row) == "PEPTIDEK"


def test_rows_without_a_peptide_are_skipped_and_observed_mz_is_the_fallback(tmp_path):
    path = tmp_path / "psm.tsv"
    pd.DataFrame(
        [
            {"Spectrum": "run.00010.00010.2", "Peptide": "", "Charge": 2, "Observed M/Z": 500.0},
            {"Spectrum": "run.00011.00011.2", "Peptide": "PEPTIDEK", "Charge": "", "Observed M/Z": 464.7},
        ]
    ).to_csv(path, sep="\t", index=False)

    [candidate] = load_fragpipe_psms(path)

    assert (candidate.scan_number, candidate.precursor_charge, candidate.precursor_mz) == (11, 2, 464.7)
    assert candidate.confidence is None and candidate.raw_score is None
