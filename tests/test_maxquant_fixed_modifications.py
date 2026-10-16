"""MaxQuant fixed modifications (#101): "Modified sequence" lists only variable modifications.

The rows are from the manuscript supplement (Zenodo 23295902, processed/proteome_tools/maxquant_msms.txt), whose
maxquant_parameters.txt sets "Fixed modifications: Carbamidomethyl (C)". Expected masses are MaxQuant's own "Mass"
column, checked against pyteomics residue masses plus the Unimod carbamidomethyl and oxidation deltas.
"""

from __future__ import annotations

import re

import pytest
from pyteomics import mass

from proteomics_mcp.adapters.maxquant import (
    MASS_TOLERANCE_DA,
    NO_FIXED_MODIFICATIONS_WARNING,
    find_maxquant_parameters,
    load_maxquant_psms,
    read_fixed_modifications,
)

COLUMNS = ["Raw file", "Scan number", "Sequence", "Modifications", "Modified sequence", "Charge", "m/z", "Mass", "Score", "PEP"]
RAW = "01974c_BD1-TUM_missing_first_4_01_01-3xHCD-1h-R4"
ROWS = [
    [RAW, "6637", "AQNSCPICQK", "Unmodified", "_AQNSCPICQK_", "2", "603.27371", "1204.5329", "313.89", "6.079E-67"],
    [RAW, "38257", "HTYDDEFFTFK", "Unmodified", "_HTYDDEFFTFK_", "2", "725.3197", "1448.6248", "289.22", "4.0077E-60"],
    [RAW, "5597", "KAVEATQVCR", "Unmodified", "_KAVEATQVCR_", "2", "581.30587", "1160.5972", "281.41", "2.3148E-55"],
    [RAW, "21967", "CMDPRLEFEK", "Oxidation (M)", "_CM(ox)DPRLEFEK_", "2", "670.80231", "1339.5901", "82.417", "2.6621E-06"],
]
EXPECTED = {
    6637: "AQNSC[UNIMOD:4]PIC[UNIMOD:4]QK",
    38257: "HTYDDEFFTFK",
    5597: "KAVEATQVC[UNIMOD:4]R",
    21967: "C[UNIMOD:4]M[UNIMOD:35]DPRLEFEK",
}
CARBAMIDOMETHYL, OXIDATION = 57.021464, 15.994915


def _write_msms(directory, name="maxquant_msms.txt"):
    path = directory / name
    path.write_text("\n".join("\t".join(row) for row in [COLUMNS, *ROWS]) + "\n", encoding="utf-8")
    return path


def _write_parameters(directory, name="maxquant_parameters.txt", fixed="Carbamidomethyl (C)"):
    (directory / name).write_text(f"Parameter\tValue\nVersion\t1.6.0.16\nFixed modifications\t{fixed}\n", encoding="utf-8")


def test_reference_masses_need_carbamidomethyl_cysteine():
    """The independent reference: MaxQuant's Mass equals the sequence mass plus 57.021 Da per cysteine (and oxidation)."""
    for row in ROWS:
        sequence, maxquant_mass = row[2], float(row[7])
        expected = mass.fast_mass(sequence, ion_type="M", charge=0) + sequence.count("C") * CARBAMIDOMETHYL
        expected += OXIDATION if "Oxidation" in row[3] else 0.0
        assert maxquant_mass == pytest.approx(expected, abs=1e-3)


def test_fixed_modifications_are_read_from_the_parameters_file_next_to_msms(tmp_path):
    path = _write_msms(tmp_path)
    _write_parameters(tmp_path)
    assert find_maxquant_parameters(path) == tmp_path / "maxquant_parameters.txt"

    candidates = load_maxquant_psms(path)

    assert {c.scan_number: c.peptidoform for c in candidates} == EXPECTED
    assert all(c.warnings == [] for c in candidates)


def test_fixed_modifications_are_read_from_mqpar_xml(tmp_path):
    path = _write_msms(tmp_path, "msms.txt")
    (tmp_path / "mqpar.xml").write_text(
        "<MaxQuantParams><parameterGroups><parameterGroup><fixedModifications>"
        "<string>Carbamidomethyl (C)</string></fixedModifications></parameterGroup></parameterGroups></MaxQuantParams>",
        encoding="utf-8",
    )
    assert read_fixed_modifications(tmp_path / "mqpar.xml") == ["Carbamidomethyl (C)"]
    assert {c.scan_number: c.peptidoform for c in load_maxquant_psms(path)} == EXPECTED


def test_explicit_fixed_modifications_take_precedence(tmp_path):
    path = _write_msms(tmp_path)
    _write_parameters(tmp_path, fixed="")  # a parameters file that says "none" is overridden by the explicit option

    candidates = load_maxquant_psms(path, fixed_modifications="Carbamidomethyl (C)")

    assert {c.scan_number: c.peptidoform for c in candidates} == EXPECTED


def test_without_any_source_the_peptidoforms_are_unchanged_and_carry_a_warning(tmp_path):
    path = _write_msms(tmp_path)

    candidates = load_maxquant_psms(path)

    assert candidates[0].peptidoform == "AQNSCPICQK"
    assert all(NO_FIXED_MODIFICATIONS_WARNING in c.warnings for c in candidates)


def test_an_explicit_empty_list_means_no_fixed_modifications(tmp_path):
    path = _write_msms(tmp_path)

    candidates = load_maxquant_psms(path, fixed_modifications="")

    assert candidates[0].peptidoform == "AQNSCPICQK"
    assert all(NO_FIXED_MODIFICATIONS_WARNING not in c.warnings for c in candidates)


def test_a_variable_modification_on_the_same_residue_is_kept(tmp_path):
    path = tmp_path / "msms.txt"
    path.write_text(
        "\t".join(COLUMNS) + "\n" + "\t".join([RAW, "1", "PEPCK", "Oxidation (C)", "_PEPC(ox)K_", "2", "1", "1", "1", "0.01"]) + "\n",
        encoding="utf-8",
    )

    [candidate] = load_maxquant_psms(path, fixed_modifications=["Carbamidomethyl (C)"])

    assert candidate.peptidoform == "PEPC[UNIMOD:35]K"


def test_terminal_and_protein_terminal_fixed_modifications(tmp_path):
    path = tmp_path / "msms.txt"
    path.write_text(
        "\t".join(COLUMNS) + "\n" + "\t".join([RAW, "1", "PEPTIDEK", "Unmodified", "_PEPTIDEK_", "2", "1", "1", "1", "0.01"]) + "\n",
        encoding="utf-8",
    )

    [tmt] = load_maxquant_psms(path, fixed_modifications="TMT6plex (K);TMT6plex (Any N-term)")
    [protein] = load_maxquant_psms(path, fixed_modifications="Acetyl (Protein N-term)")

    assert tmt.peptidoform == "[UNIMOD:737]-PEPTIDEK[UNIMOD:737]"
    assert protein.peptidoform == "PEPTIDEK"
    assert any("protein termini" in warning for warning in protein.warnings)


def _mass_warnings(candidate):
    return [warning for warning in candidate.warnings if "MaxQuant reports" in warning]


def test_peptidoform_masses_agree_with_maxquant_when_fixed_modifications_are_applied(tmp_path):
    path = _write_msms(tmp_path)
    _write_parameters(tmp_path)

    candidates = load_maxquant_psms(path)

    assert all(_mass_warnings(c) == [] for c in candidates)


def test_missing_fixed_modifications_are_caught_by_the_mass_check(tmp_path):
    """Without the parameters file, each cysteine-containing PSM is off by 57.021 Da per cysteine."""
    path = _write_msms(tmp_path)

    by_scan = {c.scan_number: c for c in load_maxquant_psms(path)}

    assert _mass_warnings(by_scan[38257]) == []  # no cysteine
    for scan, cysteines in ((6637, 2), (5597, 1), (21967, 1)):
        [warning] = _mass_warnings(by_scan[scan])
        difference = float(re.search(r"\(([+-]\d+\.\d+) Da\)", warning).group(1))
        assert difference == pytest.approx(-cysteines * CARBAMIDOMETHYL, abs=1e-3), warning  # Mass has 4 decimals


def test_mass_check_tolerance_and_tables_without_mass(tmp_path):
    columns = [column for column in COLUMNS if column != "Mass"]
    rows = [[value for column, value in zip(COLUMNS, row) if column != "Mass"] for row in ROWS]
    path = tmp_path / "msms.txt"
    path.write_text("\n".join("\t".join(row) for row in [columns, *rows]) + "\n", encoding="utf-8")

    candidates = load_maxquant_psms(path, fixed_modifications="")

    assert MASS_TOLERANCE_DA == 0.01
    assert all(_mass_warnings(c) == [] for c in candidates)  # nothing to compare with


def test_load_candidates_passes_fixed_modifications_over_mcp_and_cli(tmp_path, capsys):
    import asyncio
    import json

    from mcp import Client

    from proteomics_mcp.cli import main
    from proteomics_mcp.server import mcp

    path = _write_msms(tmp_path)

    async def call(arguments):
        async with Client(mcp) as client:
            return await client.call_tool("load_candidates", arguments)

    result = asyncio.run(call({"path": str(path), "adapter": "maxquant", "fixed_modifications": "Carbamidomethyl (C)"}))
    assert not result.is_error, result.content
    rows = result.structured_content["result"]
    assert {row["scan_number"]: row["peptidoform"] for row in rows} == EXPECTED
    assert all(row["warnings"] == [] for row in rows)

    main(["load-candidates", "--path", str(path), "--adapter", "maxquant", "--fixed-modifications", "Carbamidomethyl (C)"])
    assert json.loads(capsys.readouterr().out) == rows

    error = asyncio.run(call({"path": str(path), "adapter": "fragpipe", "fixed_modifications": "Carbamidomethyl (C)"}))
    assert error.is_error and "applies only to adapter='maxquant'" in error.content[0].text
