"""Batch 1, items 1 and 3 of #30: adapters, tables and SDRF on the Zenodo supplement."""

from __future__ import annotations

import json

import pytest

pytestmark = pytest.mark.e2e

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


@pytest.mark.xfail(strict=True, reason="#38: MCP and CLI expose no column-mapping arguments for the generic adapter")
def test_generic_adapter_needs_column_mapping_over_mcp(zenodo, tool, report):
    """The generic adapter is documented, but MCP and CLI expose no column-mapping arguments."""
    is_error, result = tool("load_candidates", {"path": str(zenodo / "processed/arabidopsis/supplementary_scored_psms.csv"),
                                                "adapter": "generic", "limit": 5})
    report("tables", "load_candidates generic (MCP)", "fail (#38)" if is_error else "pass", {"result": str(result)[:300]})
    assert not is_error, result


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
@pytest.mark.xfail(strict=True, reason="#40: only the first 100 files of a PRIDE project are searched")
@pytest.mark.parametrize("accession", ["PXD014877", "PXD010595"])
def test_retrieve_pride_sdrf(tool, report, accession):
    is_error, remote = tool("retrieve_pride_sdrf", {"project_accession": accession, "timeout_seconds": 60})
    report("sdrf", f"retrieve_pride_sdrf {accession}", "fail (#40)" if is_error else "pass", {"result": str(remote)[:300]})
    assert not is_error, remote
