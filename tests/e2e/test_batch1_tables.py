"""Batch 1, items 1 and 3 of #30: adapters, tables and SDRF on the Zenodo supplement."""

from __future__ import annotations

import csv
import json
from importlib.metadata import version

import pytest
from packaging.version import Version

pytestmark = pytest.mark.e2e

INSTALLED = Version(version("proteomics-mcp"))

# (adapter, table inside the Zenodo supplement). The generic adapter is covered separately.
ADAPTER_TABLES = [
    ("instanovo", "processed/hela/instanovo_labelled_kpreds.csv"),
    ("instanovo", "processed/hela/instanovo_raw_kpreds.csv"),
    ("casanovo", "processed/sbrodae/casanovo_predictions.csv"),
    ("proteome_discoverer", "processed/sbrodae/proteome_discoverer_PSMs.xlsx"),
    ("fragpipe", "processed/arabidopsis/fragpipe_psm.tsv"),
    ("maxquant", "processed/proteome_tools/maxquant_msms.txt"),
]
IDS = [f"{adapter}:{table.rsplit('/', 1)[-1]}" for adapter, table in ADAPTER_TABLES]


@pytest.mark.parametrize(("adapter", "table"), ADAPTER_TABLES, ids=IDS)
def test_load_candidates(zenodo, tool, report, adapter, table):
    is_error, rows = tool("load_candidates", {"path": str(zenodo / table), "adapter": adapter, "limit": 50})
    assert not is_error, rows
    assert rows, "no candidates loaded"
    first = rows[0]
    locator = {key: first.get(key) for key in ("experiment_name", "scan_number", "spectrum_id", "usi")}
    report("tables", f"load_candidates {adapter} {table}", "pass", {
        "rows": len(rows), "peptidoform": first.get("peptidoform"), "precursor_charge": first.get("precursor_charge"),
        "confidence": first.get("confidence"), **locator,
    })
    assert all(row.get("peptidoform") for row in rows)
    assert any(value not in (None, "") for value in locator.values()), "no spectrum locator on the first candidate"


@pytest.mark.parametrize(("adapter", "table"), ADAPTER_TABLES, ids=IDS)
def test_inspect_tables_mcp_and_cli(zenodo, tool, cli, report, adapter, table):
    path = str(zenodo / table)
    is_error, psm = tool("inspect_psm_table", {"path": path, "adapter": adapter})
    assert not is_error, psm
    is_error, peptide = tool("inspect_peptide_table", {"path": path, "adapter": adapter})
    assert not is_error, peptide
    meta = psm["dataset_metadata"]
    assert meta["total_psm_count"] > 0 and meta["distinct_peptide_count"] > 0

    completed = cli("inspect-psm-table", "--path", path, "--adapter", adapter)
    assert completed.returncode == 0, completed.stderr[-2000:]
    assert json.loads(completed.stdout)["dataset_metadata"] == meta, "CLI and MCP disagree"
    report("tables", f"inspect_psm_table / inspect_peptide_table {adapter} {table}", "pass",
           {"dataset_metadata": meta, "cli_parity": True})


@pytest.mark.parametrize(("adapter", "table"), ADAPTER_TABLES, ids=IDS)
def test_render_tables(zenodo, tool, report, tmp_path, adapter, table):
    path = str(zenodo / table)
    outputs = {}
    for name in ("render_psm_table", "render_peptide_table"):
        target = tmp_path / f"{name}.html"
        is_error, result = tool(name, {"path": path, "adapter": adapter, "output_path": str(target)})
        assert not is_error, result
        assert target.is_file() and target.stat().st_size > 1000
        outputs[name] = target.stat().st_size
    report("tables", f"render_psm_table / render_peptide_table {adapter} {table}", "pass", {"html_bytes": outputs})


@pytest.mark.xfail(INSTALLED < Version("0.1.10"), strict=True,
                   reason="#38 (fixed in 0.1.10): MCP and CLI expose no column-mapping arguments for the generic adapter")
def test_generic_adapter_with_column_mapping(zenodo, tool, cli, report, tmp_path):
    """The generic adapter over MCP and the CLI, with a JSON column mapping for the Zenodo scored table."""
    table = zenodo / "processed/arabidopsis/supplementary_scored_psms.csv"
    mapping = tmp_path / "mapping.json"
    mapping.write_text(json.dumps({"peptide_column": "peptidoform", "scan_column": "scan_number", "charge_column": "precursor_charge",
                                   "precursor_mz_column": "precursor_mz", "spectrum_id_column": "spectrum_id",
                                   "score_column": "confidence", "score_type": "confidence_0_1"}))
    is_error, result = tool("load_candidates", {"path": str(table), "adapter": "generic", "column_mapping_path": str(mapping), "limit": 5})
    report("tables", "load_candidates generic with column mapping (MCP)", "fail (#38)" if is_error else "pass", {"result": str(result)[:300]})
    assert not is_error, result
    with table.open(newline="") as handle:
        expected = [(row["peptidoform"], int(row["scan_number"])) for row, _ in zip(csv.DictReader(handle), range(5))]
    assert [(row["peptidoform"], row["scan_number"]) for row in result] == expected
    completed = cli("load-candidates", "--path", str(table), "--adapter", "generic", "--column-mapping", str(mapping), "--limit", "5")
    assert completed.returncode == 0, completed.stderr[-2000:]
    fields = ("peptidoform", "scan_number", "precursor_charge", "precursor_mz", "confidence", "spectrum_id")
    assert [tuple(row[f] for f in fields) for row in json.loads(completed.stdout)] == [tuple(row[f] for f in fields) for row in result]


def test_scored_tables(zenodo, tool, report):
    scored = str(zenodo / "processed/arabidopsis/supplementary_scored_psms.csv")
    sdrf = str(zenodo / "processed/arabidopsis/sdrf.tsv")
    is_error, groups = tool("list_dataset_groups", {"scored_path": scored, "sdrf_path": sdrf})
    assert not is_error, groups
    is_error, summary = tool("summarize_dataset", {"scored_path": scored, "sdrf_path": sdrf})
    assert not is_error, summary
    first = json.loads(json.dumps(summary, default=str))
    peptide = "SHDGEGK"
    is_error, peptide_summary = tool("get_peptide_summary", {"scored_path": scored, "peptide": peptide})
    assert not is_error, peptide_summary
    is_error, peptidoform_summary = tool("get_peptidoform_summary", {"scored_path": scored, "peptidoform": peptide})
    assert not is_error, peptidoform_summary
    report("dataset", "list_dataset_groups / summarize_dataset / get_peptide_summary / get_peptidoform_summary", "pass",
           {"summary_keys": sorted(first)[:12], "groups": str(groups)[:200]})


def test_load_sdrf_local(zenodo, tool, report):
    is_error, local = tool("load_sdrf", {"path": str(zenodo / "processed/arabidopsis/sdrf.tsv")})
    assert not is_error, local
    report("sdrf", "load_sdrf arabidopsis/sdrf.tsv", "pass", {"result": str(local)[:200]})


@pytest.mark.network
@pytest.mark.xfail(INSTALLED < Version("0.1.10"), strict=True, reason="#40 (fixed in 0.1.10): only the first 100 files of a PRIDE project are searched")
@pytest.mark.parametrize("accession", ["PXD014877", "PXD010595"])
def test_retrieve_pride_sdrf(tool, report, accession):
    is_error, remote = tool("retrieve_pride_sdrf", {"project_accession": accession, "timeout_seconds": 60})
    report("sdrf", f"retrieve_pride_sdrf {accession}", "fail (#40)" if is_error else "pass", {"result": str(remote)[:300]})
    assert not is_error, remote


# The whole search-engine tables of the supplement against the engines' own numbers (#179): (adapter, table, the
# engine's mass column, offset to the neutral mass, tolerance in Da as the engine prints it, precursor tolerance in ppm)
ENGINE_TABLES = [
    # MaxQuant's default main-search tolerance; parameters.txt does not record it
    ("maxquant", "processed/proteome_tools/maxquant_msms.txt", "Mass", 0.0, 1e-4, 4.5),
    # fragger.params: precursor_true_tolerance = 20 ppm, isotope_error = 0/1/2/3
    ("fragpipe", "processed/arabidopsis/fragpipe_psm.tsv", "Calculated Peptide Mass", 0.0, 6e-4, 20.0),
    # Proteome Discoverer's default precursor tolerance
    ("proteome_discoverer", "processed/sbrodae/proteome_discoverer_PSMs.xlsx", "Theo. MH+ [Da]", -1.007276466812, 1e-4, 10.0),
]


def _engine_candidates(zenodo, adapter, table):
    from proteomics_mcp.tools.workflows import _load_candidate_models

    return _load_candidate_models(str(zenodo / table), adapter)


@pytest.mark.parametrize(("adapter", "table", "mass_column", "offset", "tolerance", "_"), ENGINE_TABLES, ids=[t[0] for t in ENGINE_TABLES])
def test_peptidoform_masses_match_the_engines(zenodo, report, adapter, table, mass_column, offset, tolerance, _):
    annotation = pytest.importorskip("proteomics_mcp.core.annotation")
    if not hasattr(annotation, "peptidoform_neutral_mass"):
        pytest.skip("the installed package predates peptidoform_neutral_mass (#73)")
    candidates = _engine_candidates(zenodo, adapter, table)

    differences = [annotation.peptidoform_neutral_mass(c.peptidoform) or float("nan") for c in candidates]
    differences = [abs(ours - float(c.source_columns[mass_column]) - offset) for ours, c in zip(differences, candidates)]
    off = [c.peptidoform for c, difference in zip(candidates, differences) if not difference <= tolerance]
    report("tables", f"peptidoform mass vs {adapter} {mass_column}", "pass" if not off else "fail",
           {"psms": len(candidates), "outside": len(off), "max_abs_difference_da": max(differences), "examples": off[:5]})
    assert not off, f"{len(off)} of {len(candidates)} peptidoforms differ from {mass_column} by more than {tolerance} Da: {off[:5]}"


@pytest.mark.parametrize(("adapter", "table", "_m", "_o", "_t", "tolerance_ppm"), ENGINE_TABLES, ids=[t[0] for t in ENGINE_TABLES])
def test_precursor_errors_are_within_the_search_tolerance(zenodo, report, adapter, table, _m, _o, _t, tolerance_ppm):
    """The precursor m/z each engine reports (engine_precursor_mass_error_ppm when scored), not the spectrum's."""
    scoring = pytest.importorskip("proteomics_mcp.core.scoring")
    if 3 not in getattr(scoring, "PRECURSOR_ISOTOPE_ERRORS", ()):
        pytest.skip("the installed package predates the +3 isotope correction and MaxQuant's observed precursor m/z (#179)")
    candidates = _engine_candidates(zenodo, adapter, table)
    if adapter == "maxquant":
        # type MSMS has no isotope pattern, so MaxQuant did not recalibrate it and its precursor may be another peak
        candidates = [c for c in candidates if c.source_columns.get("Type") != "MSMS"]

    errors = [scoring.precursor_isotope_corrected_errors(c.peptidoform, c.precursor_charge, c.precursor_mz)[1] for c in candidates]
    outside = [(c.peptidoform, error) for c, error in zip(candidates, errors) if error is None or abs(error) > tolerance_ppm]
    largest = max(abs(error) for error in errors if error is not None)
    report("tables", f"engine precursor error {adapter} within {tolerance_ppm} ppm", "pass" if not outside else "fail",
           {"psms": len(candidates), "outside": len(outside), "max_abs_ppm": largest, "examples": outside[:5]})
    assert not outside, f"{len(outside)} of {len(candidates)} precursor errors outside {tolerance_ppm} ppm: {outside[:5]}"
