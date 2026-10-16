"""limit in the inspect tools sets the preview size only; row caps elsewhere are reported, not silent (#82).

Expected counts come from the input table itself, read with pandas.
"""

from __future__ import annotations

import asyncio
import json

import pandas as pd
import pytest

from proteomics_mcp.tools.workflows import inspect_peptide_table_tool, inspect_psm_table_tool, summarize_numeric_field_tool

COLUMNS = ["Spectrum", "Peptide", "Modified Peptide", "Charge", "Calibrated Observed M/Z", "Hyperscore", "Probability", "Retention", "Protein"]
PEPTIDES = ["PEPTIDEK", "SAMPLER", "LESLIEK", "VATVSLPR", "AQNSCPICQK", "HTYDDEFFTFK"]


@pytest.fixture
def table(tmp_path):
    """30 PSMs of 6 peptides, in an order where the first 5 rows contain only 5 distinct peptides."""
    rows = []
    for i in range(30):
        peptide = PEPTIDES[i % len(PEPTIDES)]
        rows.append([f"runA.{1000 + i:05d}.{1000 + i:05d}.2", peptide, "", 2, 400.0 + i, 10.0 + i, 0.5 + i / 100, 600.0 + i, "P1"])
    path = tmp_path / "psm.tsv"
    pd.DataFrame(rows, columns=COLUMNS).to_csv(path, sep="\t", index=False)
    return path


def _reference(path):
    frame = pd.read_csv(path, sep="\t")
    return len(frame), frame["Peptide"].nunique()


def test_inspect_psm_table_limit_changes_only_the_preview(table):
    psms, peptides = _reference(table)

    full = inspect_psm_table_tool(path=str(table), adapter="fragpipe")
    limited = inspect_psm_table_tool(path=str(table), adapter="fragpipe", limit=5)

    for result in (full, limited):
        assert result["dataset_metadata"]["total_psm_count"] == psms
        assert result["dataset_metadata"]["distinct_peptide_count"] == peptides
        assert "row_limit" not in result and "warnings" not in result
    assert {key: value for key, value in limited.items() if not key.endswith("_preview")} == {
        key: value for key, value in full.items() if not key.endswith("_preview")}
    assert len(full["psm_preview"]) == 10 and len(limited["psm_preview"]) == 5
    assert len(limited["peptide_preview"]) == 5
    assert limited["psm_preview"] == full["psm_preview"][:5]


def test_inspect_peptide_table_limit_changes_only_the_preview(table):
    psms, peptides = _reference(table)

    limited = inspect_peptide_table_tool(path=str(table), adapter="fragpipe", limit=2)

    assert limited["dataset_metadata"]["total_psm_count"] == psms
    assert limited["dataset_metadata"]["distinct_peptide_count"] == peptides
    assert len(limited["peptide_preview"]) == 2


def test_max_rows_reads_only_the_first_rows_and_says_so(table):
    result = inspect_psm_table_tool(path=str(table), adapter="fragpipe", max_rows=5)

    assert result["dataset_metadata"]["total_psm_count"] == 5
    assert result["row_limit"] == {"argument": "max_rows", "row_limit": 5, "rows_analyzed": 5}
    assert len(result["warnings"]) == 1 and "max_rows=5: only the first 5 rows" in result["warnings"][0]


def test_row_caps_of_the_summary_tools_are_reported(table):
    psms, _ = _reference(table)

    full = summarize_numeric_field_tool(path=str(table), adapter="fragpipe", field="confidence")
    capped = summarize_numeric_field_tool(path=str(table), adapter="fragpipe", field="confidence", limit=5)

    assert full["count"] == psms and "warnings" not in full
    assert capped["count"] == 5
    assert capped["row_limit"] == {"argument": "limit", "row_limit": 5, "rows_analyzed": 5}
    assert "limit=5: only the first 5 rows" in capped["warnings"][0]


def test_negative_preview_size_is_rejected(table):
    with pytest.raises(ValueError, match="limit must be 0 or more"):
        inspect_psm_table_tool(path=str(table), adapter="fragpipe", limit=-1)


def test_inspect_psm_table_limit_over_mcp_and_cli(table, capsys):
    from mcp import Client

    from proteomics_mcp.cli import main
    from proteomics_mcp.server import mcp

    psms, _ = _reference(table)

    async def call(arguments):
        async with Client(mcp) as client:
            return await client.call_tool("inspect_psm_table", arguments)

    result = asyncio.run(call({"path": str(table), "adapter": "fragpipe", "limit": 3}))
    assert not result.is_error, result.content
    payload = result.structured_content
    assert payload["dataset_metadata"]["total_psm_count"] == psms and len(payload["psm_preview"]) == 3

    capped = asyncio.run(call({"path": str(table), "adapter": "fragpipe", "max_rows": 4}))
    assert capped.structured_content["row_limit"]["row_limit"] == 4

    main(["inspect-psm-table", "--path", str(table), "--adapter", "fragpipe", "--limit", "3"])
    cli_result = json.loads(capsys.readouterr().out)
    assert cli_result["dataset_metadata"] == payload["dataset_metadata"] and len(cli_result["psm_preview"]) == 3
