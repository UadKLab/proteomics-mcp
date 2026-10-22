"""One run-matching rule for joins, lookups and scoring (#142), gzipped file names (#126), and a metadata join that
reports what it matched and the retention-time unit it used (#146)."""

from __future__ import annotations

import gzip
import json
import shutil
from pathlib import Path

import pandas as pd
import pytest

import proteomics_mcp.tools.workflows as workflows
from proteomics_mcp.core import spectrum_loader
from proteomics_mcp.core.runs import candidate_run, run_from_spectrum_id, run_name
from proteomics_mcp.core.spectrum_loader import load_spectrum_from_table
from proteomics_mcp.tools.workflows import build_prediction_batch_tool, enrich_frame_with_spectra_metadata, score_candidates_tool

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("/data/SF_rep1.mzML.gz", "SF_rep1"),
        ("SF_rep1.mzML", "SF_rep1"),
        ("C:\\runs\\SF_rep1.raw", "SF_rep1"),
        ("run.d/", "run"),
        ("A.mgf.bz2", "A"),
        ("QE-HF-01.5000", "QE-HF-01.5000"),  # a dot that is not a file suffix stays
        ("SF_rep1", "SF_rep1"),
        ("", None),
        (float("nan"), None),
    ],
)
def test_run_name(value, expected):
    assert run_name(value) == expected


def test_candidate_run_prefers_the_experiment_then_the_spectrum_id():
    assert candidate_run("A.mzML", "B:12") == "A"
    assert candidate_run(None, "B:12") == "B"
    assert candidate_run(None, "mzspec:PXD1:C:scan:7") == "C"
    assert candidate_run(None, "controllerType=0 controllerNumber=1 scan=12") is None
    assert run_from_spectrum_id("B:12") == ("B", 12)


# --- #126: gzipped inputs keep the run name of the uncompressed file ------------------------------------------------


def test_gzipped_mzml_gets_the_run_name_of_the_uncompressed_file(tmp_path):
    pytest.importorskip("pyopenms")
    gz = tmp_path / "SF_rep1.mzML.gz"
    with (FIXTURES / "mzml_isolation_injection.mzML").open("rb") as source, gzip.open(gz, "wb") as target:
        shutil.copyfileobj(source, target)
    spectrum_loader.mzml_to_feather(gz, tmp_path / "out.feather")

    assert set(pd.read_feather(tmp_path / "out.feather")["source_name"]) == {"SF_rep1"}


def test_mgf_run_name_drops_compression_and_mgf_suffixes(tmp_path):
    mgf = tmp_path / "SF_rep1.mgf"
    mgf.write_text("BEGIN IONS\nTITLE=SF_rep1.42.42.2\nCHARGE=2+\nPEPMASS=500.2\n100.0 10.0\nEND IONS\n", encoding="utf-8")
    spectrum_loader.mgf_to_feather(mgf, tmp_path / "out.feather")

    assert set(pd.read_feather(tmp_path / "out.feather")["source_name"]) == {"SF_rep1"}


# --- #142: one run rule for every layout and every adapter --------------------------------------------------------


def _multi_run_table(path: Path, run_column: str, *, suffix: str = "") -> Path:
    """Two runs that share scan numbers 5 and 7; each spectrum's single peak identifies it."""
    rows = []
    for run, base in (("runA", 100.0), ("runB", 200.0)):
        for scan in (5, 7):
            rows.append({run_column: f"{run}{suffix}", "scan_number": scan, "precursor_mz": 500.0, "precursor_charge": 2,
                         "mz_array": [base + scan], "intensity_array": [1.0]})
    pd.DataFrame(rows).to_feather(path)
    return path


@pytest.mark.parametrize(("run_column", "suffix"), [("source_file", ".mzML"), ("source_name", ""), ("experiment_name", ""), ("source_file", ".mzML.gz")])
def test_lookup_finds_the_spectrum_of_the_named_run_in_every_layout(tmp_path, run_column, suffix):
    table = _multi_run_table(tmp_path / "s.feather", run_column, suffix=suffix)

    assert load_spectrum_from_table(table, scan_number=7, experiment_name="runB").mz == [207.0]
    assert load_spectrum_from_table(table, scan_number=5, experiment_name="runA.mzML").mz == [105.0]
    with pytest.raises(LookupError, match="different spectra match scan 7 in runs runA, runB"):
        load_spectrum_from_table(table, scan_number=7)


def test_a_spectrum_id_of_run_and_scan_is_looked_up_by_run_and_scan(tmp_path):
    table = _multi_run_table(tmp_path / "s.feather", "experiment_name")

    assert load_spectrum_from_table(table, spectrum_id="runB:5").mz == [205.0]


def _generic_table(tmp_path: Path, rows: list[dict], mapping: dict) -> tuple[Path, Path]:
    table = tmp_path / "psms.csv"
    pd.DataFrame(rows).to_csv(table, index=False)
    mapping_path = tmp_path / "mapping.json"
    mapping_path.write_text(json.dumps({"peptide_column": "peptide", "charge_column": "charge", **mapping}), encoding="utf-8")
    return table, mapping_path


@pytest.mark.parametrize("mapping", [{"spectrum_id_column": "spectrum"}, {"run_column": "run", "scan_column": "scan"}])
def test_generic_candidates_score_against_their_own_run(tmp_path, monkeypatch, mapping):
    table = _multi_run_table(tmp_path / "s.feather", "experiment_name")
    rows = [{"peptide": "PEPTIDEK", "charge": 2, "spectrum": f"{run}:{scan}", "run": run, "scan": scan} for run in ("runA", "runB") for scan in (5, 7)]
    candidates, mapping_path = _generic_table(tmp_path, rows, mapping)
    seen = []
    real_score = workflows.score_candidate_spectrum

    def remember(**kwargs):
        seen.append((kwargs["candidate"].spectrum_id, kwargs["candidate"].experiment_name, kwargs["observed"].mz[0]))
        return real_score(**kwargs)

    monkeypatch.setattr(workflows, "score_candidate_spectrum", remember)
    summary = score_candidates_tool(candidates_path=str(candidates), spectra_path=str(table), adapter="generic",
                                    column_mapping_path=str(mapping_path), predictor_backend="mock")

    assert summary["rows_failed"] == 0, summary
    assert sorted(mz for *_, mz in seen) == [105.0, 107.0, 205.0, 207.0]


def test_generic_mapping_rejects_two_different_run_columns(tmp_path):
    candidates, mapping_path = _generic_table(tmp_path, [{"peptide": "PEPTIDEK", "charge": 2, "a": "x", "b": "y"}],
                                              {"run_column": "a", "experiment_name_column": "b"})
    with pytest.raises(ValueError, match="same field"):
        workflows.load_candidates_tool(str(candidates), adapter="generic", column_mapping_path=str(mapping_path))


# --- #146: the metadata join reports what it matched ----------------------------------------------------------------


def _instanovo_layout(path: Path, *, minutes: bool = True) -> Path:
    """The InstaNovo spectra layout: experiment_name, retention_time (here in minutes), no source_file."""
    rows = [{"experiment_name": "runA", "scan_number": scan, "precursor_mz": 500.0, "precursor_charge": 2,
             "retention_time": scan / 2 if minutes else scan * 30.0, "mz_array": [100.0], "intensity_array": [1.0]}
            for scan in range(1, 41)]
    pd.DataFrame(rows).to_feather(path)
    return path


def _frame(rows: list[dict]) -> pd.DataFrame:
    return pd.DataFrame([{"candidate_id": f"c{i}", "peptidoform": "PEPTIDEK", **row} for i, row in enumerate(rows)])


def test_join_keeps_every_row_and_counts_the_unmatched_by_reason(tmp_path):
    spectra = _instanovo_layout(tmp_path / "s.feather")
    frame = _frame([
        {"run_key": "runA", "scan_number": 4},
        {"run_key": "runA", "scan_number": 4000},
        {"run_key": "runZ", "scan_number": 4},
        {"run_key": "runA", "scan_number": None},
        {"run_key": None, "scan_number": 4},
    ])

    joined, summary = enrich_frame_with_spectra_metadata(frame, spectra_path=str(spectra), rt_unit="minutes")

    assert len(joined) == 5
    assert (summary["matched_rows"], summary["unmatched_rows"]) == (1, 4)
    assert summary["unmatched_reasons"] == {"no_run": 1, "no_scan_number": 1, "no_spectrum_for_run_and_scan": 2}
    assert summary["runs_not_in_spectra"] == {"runZ": 1}
    assert (summary["spectra_run_column"], summary["rt_source_field"], summary["rt_unit"]) == ("experiment_name", "retention_time", "minutes")
    assert joined.loc[0, "retention_time_seconds"] == pytest.approx(120.0)  # scan 4: 2.0 min
    assert joined.loc[0, "retention_time_minutes"] == pytest.approx(2.0)


def test_join_warns_when_seconds_look_like_minutes(tmp_path):
    spectra = _instanovo_layout(tmp_path / "s.feather", minutes=True)
    frame = _frame([{"run_key": "runA", "scan_number": scan} for scan in range(1, 41)])

    _, summary = enrich_frame_with_spectra_metadata(frame, spectra_path=str(spectra))

    assert summary["rt_unit"] == "seconds"
    assert any("rt_unit='minutes'" in warning for warning in summary["warnings"])
    _, in_minutes = enrich_frame_with_spectra_metadata(frame, spectra_path=str(spectra), rt_unit="minutes")
    assert in_minutes["warnings"] == []


def test_an_unknown_rt_column_or_unit_is_an_error(tmp_path):
    spectra = _instanovo_layout(tmp_path / "s.feather")
    frame = _frame([{"run_key": "runA", "scan_number": 4}])
    with pytest.raises(ValueError, match="not a column"):
        enrich_frame_with_spectra_metadata(frame, spectra_path=str(spectra), rt_source_field="rt")
    with pytest.raises(ValueError, match="rt_unit"):
        enrich_frame_with_spectra_metadata(frame, spectra_path=str(spectra), rt_unit="hours")


def test_prediction_batch_keeps_candidates_without_scans_and_explains_a_failed_join(tmp_path):
    spectra = _instanovo_layout(tmp_path / "s.feather")
    # A mapping without run or scan columns: every row has a charge, none can be joined
    candidates, mapping_path = _generic_table(tmp_path, [{"peptide": "PEPTIDEK", "charge": 2}, {"peptide": "SAMPLER", "charge": 3}], {})

    manifest = build_prediction_batch_tool(path=str(candidates), output_prefix=str(tmp_path / "batch"), spectra_path=str(spectra),
                                           adapter="generic", column_mapping_path=str(mapping_path), predictor_backend="mock")

    assert manifest["predictions_succeeded"] == 2
    assert manifest["metadata_join"]["matched_rows"] == 0
    assert manifest["metadata_join"]["unmatched_reasons"]["no_run"] == 2


def test_a_failed_join_is_named_when_no_candidate_has_a_charge(tmp_path):
    spectra = _instanovo_layout(tmp_path / "s.feather")
    candidates, mapping_path = _generic_table(tmp_path, [{"peptide": "PEPTIDEK", "charge": None}], {})

    with pytest.raises(ValueError, match="No candidate matched a spectrum: 1 without a run, 0 without a scan number"):
        build_prediction_batch_tool(path=str(candidates), output_prefix=str(tmp_path / "batch"), spectra_path=str(spectra),
                                    adapter="generic", column_mapping_path=str(mapping_path), predictor_backend="mock")
