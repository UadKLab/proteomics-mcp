"""The generic adapter over MCP and the CLI: a JSON column mapping on every tool and command that takes an adapter (#38)."""

from __future__ import annotations

import argparse
import asyncio
import json

import pandas as pd
import pytest

import proteomics_mcp.tools.workflows as workflows
from proteomics_mcp.cli import build_parser, main
from proteomics_mcp.server import mcp
from proteomics_mcp.tools.workflows import GENERIC_MAPPING_KEYS, build_prediction_batch_tool, load_candidates_tool, score_from_prediction_batch_tool

MAPPING = {"peptide_column": "Sequence", "scan_column": "Scan", "charge_column": "Charge", "score_column": "Score",
           "score_type": "confidence_0_1", "spectrum_id_column": "Run", "source_engine": "mapped"}


@pytest.fixture
def generic_table(tmp_path):
    table = tmp_path / "mapped.csv"
    pd.DataFrame([
        {"Sequence": "PEPTIDE", "Scan": 10, "Charge": 2, "Score": 0.9, "Run": "run"},
        {"Sequence": "PEPTIDEK", "Scan": 11, "Charge": 2, "Score": 0.7, "Run": "run"},
    ]).to_csv(table, index=False)
    mapping = tmp_path / "mapping.json"
    mapping.write_text(json.dumps(MAPPING), encoding="utf-8")
    return table, mapping


def _call_tool(name: str, arguments: dict):
    from mcp import Client

    async def call():
        async with Client(mcp) as client:
            return await client.call_tool(name, arguments)

    return asyncio.run(call())


def test_mapping_keys_are_the_generic_adapter_arguments():
    assert "peptide_column" in GENERIC_MAPPING_KEYS and "score_type" in GENERIC_MAPPING_KEYS
    assert "limit" not in GENERIC_MAPPING_KEYS


def test_load_candidates_with_a_mapping_over_mcp_and_cli(generic_table, capsys):
    table, mapping = generic_table
    result = _call_tool("load_candidates", {"path": str(table), "adapter": "generic", "column_mapping_path": str(mapping)})
    assert not result.is_error, result.content
    rows = result.structured_content["result"]
    assert [(row["peptidoform"], row["scan_number"], row["precursor_charge"], row["confidence"]) for row in rows] == [
        ("PEPTIDE", 10, 2, 0.9), ("PEPTIDEK", 11, 2, 0.7)]
    assert {row["source_engine"] for row in rows} == {"mapped"}

    main(["load-candidates", "--path", str(table), "--adapter", "generic", "--column-mapping", str(mapping)])
    assert json.loads(capsys.readouterr().out) == rows


@pytest.mark.parametrize(
    ("arguments", "message"),
    [
        ({"adapter": "generic"}, "needs a column mapping"),
        ({"adapter": "fragpipe", "mapping": MAPPING}, "applies only to adapter='generic'"),
        ({"adapter": "generic", "mapping": {**MAPPING, "peptide": "Sequence"}}, "Unknown column mapping keys: peptide"),
        ({"adapter": "generic", "mapping": {**MAPPING, "peptide_column": "Missing"}}, "Missing required mapped columns: Missing"),
        ({"adapter": "generic", "mapping": ["Sequence"]}, "must be a JSON object"),
    ],
    ids=["no-mapping", "not-generic", "unknown-key", "unknown-column", "not-an-object"],
)
def test_mapping_errors_are_clear(generic_table, tmp_path, arguments, message):
    table, _ = generic_table
    kwargs = {}
    if "mapping" in arguments:
        path = tmp_path / "bad_mapping.json"
        path.write_text(json.dumps(arguments["mapping"]), encoding="utf-8")
        kwargs["column_mapping_path"] = str(path)
    with pytest.raises(ValueError, match=message):
        load_candidates_tool(str(table), arguments["adapter"], **kwargs)


def _dummy(schema: dict):
    kind = schema.get("type") or next((option.get("type") for option in schema.get("anyOf", []) if option.get("type") != "null"), "string")
    return {"string": "x", "integer": 1, "number": 1.0, "boolean": False, "array": [], "object": {}}.get(kind, "x")


def _adapter_tools() -> dict:
    tools = {tool.name: tool for tool in mcp._tool_manager.list_tools()}
    return {name: tool for name, tool in tools.items() if "adapter" in tool.parameters.get("properties", {})}


def test_every_mcp_tool_with_an_adapter_accepts_a_column_mapping(generic_table, monkeypatch):
    """Each tool forwards column_mapping_path to the generic adapter (what happens after loading does not matter here)."""
    table, mapping = generic_table
    tools = _adapter_tools()
    assert len(tools) >= 20
    seen: dict[str, dict] = {}

    def recording_loader(path, **kwargs):
        seen[current] = kwargs
        raise ValueError("stop after loading")

    monkeypatch.setattr(workflows, "load_generic_candidates", recording_loader)
    for current, tool in tools.items():
        properties = tool.parameters["properties"]
        assert "column_mapping_path" in properties, current
        arguments = {name: _dummy(properties[name]) for name in tool.parameters.get("required", [])}
        arguments.update({"adapter": "generic", "column_mapping_path": str(mapping)})
        for path_argument in ("path", "candidates_path"):
            if path_argument in properties:
                arguments[path_argument] = str(table)
        _call_tool(current, arguments)
    missing = sorted(name for name in tools if seen.get(name, {}).get("peptide_column") != "Sequence")
    assert missing == []


def _commands_with_adapter() -> dict[str, argparse.ArgumentParser]:
    parser = build_parser()
    subparsers = next(action for action in parser._actions if isinstance(action, argparse._SubParsersAction))
    return {name: sub for name, sub in subparsers.choices.items()
            if any("--adapter" in action.option_strings for action in sub._actions)}


def test_every_cli_command_with_an_adapter_accepts_a_column_mapping(generic_table, monkeypatch, capsys):
    table, mapping = generic_table
    commands = _commands_with_adapter()
    assert len(commands) >= 20
    seen: dict[str, dict] = {}

    def recording_loader(path, **kwargs):
        seen[current] = kwargs
        raise ValueError("stop after loading")

    monkeypatch.setattr(workflows, "load_generic_candidates", recording_loader)
    for current, sub in commands.items():
        argv = [current, "--adapter", "generic", "--column-mapping", str(mapping)]
        for action in sub._actions:
            if action.required and action.option_strings:
                option = action.option_strings[0]
                argv += [option, str(table) if "path" in option else "1" if action.type in (int, float) else "x"]
        with pytest.raises((ValueError, SystemExit)):
            main(argv)
    capsys.readouterr()
    missing = sorted(name for name in commands if seen.get(name, {}).get("peptide_column") != "Sequence")
    assert missing == []


def test_generic_prediction_batch_is_scored_with_its_mapping(generic_table, tmp_path, monkeypatch):
    """score_from_prediction_batch reloads the candidates with the mapping recorded in the prediction manifest."""
    table, mapping = generic_table
    monkeypatch.setattr(workflows, "load_spectrum_from_table", lambda *_, **__: workflows.ObservedSpectrum(
        mz=[200.0, 300.0, 400.0], intensity=[10.0, 20.0, 5.0], precursor_mz=400.0, precursor_charge=2))
    manifest = build_prediction_batch_tool(path=str(table), output_prefix=str(tmp_path / "batch"), adapter="generic",
                                           column_mapping_path=str(mapping), predictor_backend="mock")
    assert manifest["column_mapping_path"] == str(mapping)
    assert manifest["predictions_succeeded"] == 2
    spectra = tmp_path / "spectra.feather"
    pd.DataFrame([{"scan_number": scan, "spectrum_id": f"run:{scan}", "mz_array": [200.0, 300.0, 400.0],
                   "intensity_array": [10.0, 20.0, 5.0], "precursor_mz": 400.0, "precursor_charge": 2} for scan in (10, 11)]).to_feather(spectra)
    scored = score_from_prediction_batch_tool(manifest_path=manifest["output_files"]["manifest_json"], spectra_path=str(spectra),
                                              output_prefix=str(tmp_path / "scored"))
    assert scored["candidate_rows_requested"] == 2
