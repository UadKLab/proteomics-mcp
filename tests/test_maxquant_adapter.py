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
    # msms.txt's m/z is theoretical; without MaxQuant's mass errors there is no observed precursor m/z to give
    assert candidate.precursor_mz is None
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


# The observed precursor m/z, rebuilt from MaxQuant's mass errors (#179): the Figure 3B PSM
# AQNSC[UNIMOD:4]PIC[UNIMOD:4]QK, scan 6637, whose mzML precursor is 603.274597 (1.47 ppm from the theoretical 603.27371)


def _one_row(tmp_path, **errors):
    columns = ["Raw file", "Scan number", "Sequence", "Modified sequence", "Charge", "m/z", "Mass", *errors]
    values = ["run", "6637", "AQNSCPICQK", "_AQNSCPICQK_", "2", "603.27371", "1204.5329", *[str(v) for v in errors.values()]]
    path = tmp_path / "msms.txt"
    path.write_text("\t".join(c.replace("_", " ").replace("Mass Error ppm", "Mass Error [ppm]") for c in columns) + "\n" + "\t".join(values) + "\n",
                    encoding="utf-8")
    return path


def test_precursor_mz_is_the_recalibrated_observation(tmp_path):
    path = _one_row(tmp_path, Mass_Error_ppm=-0.020891, Simple_Mass_Error_ppm=1.4698311)

    [candidate] = load_maxquant_psms(path, fixed_modifications="Carbamidomethyl (C)")

    assert candidate.precursor_mz == pytest.approx(603.27371 * (1 - 0.020891e-6), abs=1e-9)
    assert candidate.precursor_mz != 603.27371  # not the theoretical m/z


def test_without_recalibration_the_simple_error_gives_the_instrument_precursor(tmp_path):
    # MaxQuant leaves "Mass Error [ppm]" empty for type MSMS (no isotope pattern)
    path = _one_row(tmp_path, Mass_Error_ppm="", Simple_Mass_Error_ppm=1.4698311)

    [candidate] = load_maxquant_psms(path, fixed_modifications="Carbamidomethyl (C)")

    assert candidate.precursor_mz == pytest.approx(603.274597, abs=2e-6)  # the mzML precursor of scan 6637
