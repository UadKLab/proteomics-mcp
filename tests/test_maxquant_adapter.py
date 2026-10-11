from __future__ import annotations

import pytest

from proteomics_mcp.adapters.maxquant import load_maxquant_psms


def test_load_maxquant_psms_extracts_join_and_score_fields(tmp_path):
    path = tmp_path / "msms.txt"
    path.write_text(
        "\t".join(
            [
                "Raw file",
                "Scan number",
                "Sequence",
                "Modified sequence",
                "Charge",
                "m/z",
                "Score",
                "Delta score",
                "PEP",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "01974c_BD1-TUM_missing_first_4_01_01-3xHCD-1h-R4",
                "14345",
                "AAANGLMQTSK",
                "_AAANGLMQTSK_",
                "2",
                "546.27932",
                "131.120",
                "129.790",
                "2.2652e-13",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    [candidate] = load_maxquant_psms(path)

    assert candidate.source_engine == "maxquant"
    assert candidate.experiment_name == "01974c_BD1-TUM_missing_first_4_01_01-3xHCD-1h-R4"
    assert candidate.scan_number == 14345
    assert candidate.spectrum_id == "01974c_BD1-TUM_missing_first_4_01_01-3xHCD-1h-R4:14345"
    assert candidate.peptidoform == "AAANGLMQTSK"
    assert candidate.precursor_charge == 2
    assert candidate.precursor_mz == 546.27932
    assert candidate.raw_score == 131.12
    assert candidate.raw_score_name == "Score"
    assert candidate.next_beam is None
    assert candidate.next_beam_margin is None
    assert candidate.source_columns["Delta score"] == pytest.approx(129.79)
    assert candidate.confidence is not None and candidate.confidence > 0.999999999999


def test_load_maxquant_psms_normalizes_inline_modifications(tmp_path):
    path = tmp_path / "msms.txt"
    path.write_text(
        "\t".join(["Raw file", "Scan number", "Sequence", "Modified sequence", "Charge", "m/z", "Modifications"])
        + "\n"
        + "\t".join(["runA", "244", "AETLMQARKR", "_AETLM(ox)QARKR_", "2", "617.845", "Oxidation (M)"])
        + "\n"
        + "\t".join(["runA", "725", "ALMSPAGMLR", "_ALM(ox)SPAGM(ox)LR_", "2", "531.242", "2 Oxidation (M)"])
        + "\n",
        encoding="utf-8",
    )

    candidates = load_maxquant_psms(path)

    assert candidates[0].peptidoform == "AETLM[UNIMOD:35]QARKR"
    assert candidates[1].peptidoform == "ALM[UNIMOD:35]SPAGM[UNIMOD:35]LR"
