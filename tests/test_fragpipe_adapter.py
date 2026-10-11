from __future__ import annotations

from proteomics_mcp.adapters.fragpipe import load_fragpipe_psms


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
