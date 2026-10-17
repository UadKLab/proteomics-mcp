from __future__ import annotations

import glob

import pandas as pd
import pytest

from proteomics_mcp.adapters.casanovo import load_casanovo_predictions
from proteomics_mcp.adapters.fragpipe import load_fragpipe_psms
from proteomics_mcp.adapters.instanovo import load_instanovo_predictions
from proteomics_mcp.adapters.maxquant import load_maxquant_psms
from proteomics_mcp.adapters.table_io import read_candidate_table

TABLE = pd.DataFrame({"scan_number": [1, 2], "predictions": ["PEPTIDEK", "PEPTLDEK"], "token_log_probs": ["[-0.1, -0.2]", "[-0.3, -0.4]"]})


def _write(frame: pd.DataFrame, path, sep: str) -> None:
    frame.to_csv(path, sep=sep, index=False)


@pytest.mark.parametrize(
    ("name", "sep"),
    [
        ("table.csv", ","),
        ("table.tsv", "\t"),
        ("table.txt", "\t"),
        ("table.tab", "\t"),
        ("table.csv.gz", ","),
        ("table.tsv.gz", "\t"),
        ("table.tsv.bz2", "\t"),
        ("table", "\t"),  # no extension: delimiter from the header line
        ("table.dat", ";"),
        # content that does not match the extension (#140)
        ("tabs.csv", "\t"),
        ("commas.tsv", ","),
        ("commas.txt", ","),
    ],
)
def test_read_candidate_table_formats(tmp_path, name, sep):
    path = tmp_path / name
    _write(TABLE, path, sep)

    frame = read_candidate_table(path)

    assert list(frame.columns) == list(TABLE.columns)
    assert frame["predictions"].tolist() == ["PEPTIDEK", "PEPTLDEK"]
    assert frame["token_log_probs"].tolist() == ["[-0.1, -0.2]", "[-0.3, -0.4]"]


def test_read_candidate_table_respects_limit_and_excel(tmp_path):
    _write(TABLE, tmp_path / "table.tsv", "\t")
    TABLE.to_excel(tmp_path / "table.xlsx", index=False)

    assert len(read_candidate_table(tmp_path / "table.tsv", limit=1)) == 1
    assert read_candidate_table(tmp_path / "table.xlsx")["predictions"].tolist() == ["PEPTIDEK", "PEPTLDEK"]


def test_instanovo_reads_tab_separated_predictions(tmp_path):
    fixture = sorted(glob.glob("tests/fixtures/instanovo_formats/*.csv"))[0]
    tsv = tmp_path / "predictions.tsv"
    pd.read_csv(fixture).to_csv(tsv, sep="\t", index=False)

    from_csv = load_instanovo_predictions(fixture)
    from_tsv = load_instanovo_predictions(tsv)

    assert [c.peptidoform for c in from_tsv] == [c.peptidoform for c in from_csv]
    assert [c.confidence for c in from_tsv] == [c.confidence for c in from_csv]


def test_instanovo_schema_error_lists_the_columns(tmp_path):
    path = tmp_path / "not_instanovo.csv"
    _write(pd.DataFrame({"peptide": ["PEPTIDEK"], "score": [0.5]}), path, ",")

    with pytest.raises(ValueError, match="Columns found: peptide, score"):
        load_instanovo_predictions(path)


@pytest.mark.parametrize("name", ["casanovo.csv", "casanovo.tsv", "casanovo.tsv.gz"])
def test_casanovo_reads_every_separator(tmp_path, name):
    frame = pd.DataFrame({"prediction_id": [0, 1], "predictions": ["PEPTIDEK", "M[Oxidation]PEPTIDEK"], "log_probs": [-0.1, -0.2], "precursor_mz": [466.7, 547.3], "precursor_charge": [2, 2]})
    path = tmp_path / name
    _write(frame, path, "," if ".csv" in name else "\t")

    candidates = load_casanovo_predictions(path)

    assert [c.peptidoform for c in candidates] == ["PEPTIDEK", "M[UNIMOD:35]PEPTIDEK"]


def test_casanovo_without_predictions_column_is_an_error_not_an_empty_table(tmp_path):
    path = tmp_path / "casanovo.csv"
    _write(pd.DataFrame({"prediction_id": [0], "sequence": ["PEPTIDEK"]}), path, ",")

    with pytest.raises(ValueError, match="Missing required Casanovo column 'predictions'. Columns found: prediction_id, sequence"):
        load_casanovo_predictions(path)


def test_maxquant_and_fragpipe_accept_comma_separated_exports(tmp_path):
    maxquant = pd.DataFrame({"Raw file": ["run"], "Scan number": [7], "Modified sequence": ["_PEPM(ox)K_"], "Sequence": ["PEPMK"], "Charge": [2], "m/z": [310.6]})
    _write(maxquant, tmp_path / "msms.csv", ",")
    _write(maxquant, tmp_path / "msms.txt", "\t")
    fragpipe = pd.DataFrame({"Spectrum": ["run.00007.00007.2"], "Peptide": ["PEPMK"], "Assigned Modifications": ["4M(15.9949)"], "Charge": [2]})
    _write(fragpipe, tmp_path / "psm.csv", ",")
    _write(fragpipe, tmp_path / "psm.tsv", "\t")

    for path in ("msms.csv", "msms.txt"):
        [candidate] = load_maxquant_psms(tmp_path / path, fixed_modifications="")
        assert candidate.peptidoform == "PEPM[UNIMOD:35]K"
    for path in ("psm.csv", "psm.tsv"):
        [candidate] = load_fragpipe_psms(tmp_path / path)
        assert candidate.peptidoform == "PEPM[UNIMOD:35]K"
