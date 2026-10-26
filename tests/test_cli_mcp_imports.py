from __future__ import annotations

import asyncio
import importlib.util
import inspect
import json
import sys
from pathlib import Path

import pytest

import proteomics_mcp.server as server_module
from proteomics_mcp import __version__
from proteomics_mcp.cli import build_parser
from proteomics_mcp.server import mcp
from proteomics_mcp.resources.docs import SCHEMAS

MANIFEST_PATH = Path(__file__).resolve().parents[1] / "manifest.json"


def test_cli_has_expected_commands():
    parser = build_parser()
    help_text = parser.format_help()
    assert "proteomics-mcp" in help_text
    assert "load-sdrf" in help_text
    assert "retrieve-pride-sdrf" in help_text
    assert "inspect-psm-table" in help_text
    assert "inspect-peptide-table" in help_text
    assert "list-distribution-fields" in help_text
    assert "summarize-numeric-field" in help_text
    assert "summarize-categorical-field" in help_text
    assert "summarize-missing-values" in help_text
    assert "summarize-cv-by-category" in help_text
    assert "bin-numeric-field" in help_text
    assert "bin-2d-numeric-fields" in help_text
    assert "summarize-dashboard" in help_text
    assert "summarize-overlap" in help_text
    assert "summarize-modifications" in help_text
    assert "summarize-quant-quality" in help_text
    assert "summarize-file-coverage" in help_text
    assert "render-summary-html" in help_text
    assert "render-dashboard-html" in help_text
    assert "build-prediction-batch" in help_text
    assert "inspect-prediction-batch" in help_text
    assert "inspect-prediction-batch-status" in help_text
    assert "score-from-prediction-batch" in help_text
    assert "inspect-scored-batch" in help_text
    assert "inspect-scored-batch-status" in help_text
    assert "render-psm-table" in help_text
    assert "render-peptide-table" in help_text
    assert "summarize-dataset" in help_text
    assert "render-dataset-report" in help_text
    assert "render-spectra-batch" in help_text
    assert "mgf-to-feather" in help_text
    assert "bruker-d-to-feather" in help_text
    assert "export-to-instanovo" in help_text


def test_server_module_builds_an_mcp_server_when_mcp_is_installed():
    from mcp.server import MCPServer

    assert isinstance(mcp, MCPServer)


def test_server_module_imports_without_mcp_dependency(monkeypatch):
    # Simulate an environment without the mcp package: a None entry in sys.modules makes its import fail.
    for name in [name for name in sys.modules if name == "mcp" or name.startswith("mcp.")]:
        monkeypatch.setitem(sys.modules, name, None)
    monkeypatch.setitem(sys.modules, "mcp", None)
    # Load a separate copy, so the server module the other tests use is left as it is.
    spec = importlib.util.spec_from_file_location("proteomics_mcp_server_without_mcp", server_module.__file__)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    assert module.MCPServer is None
    assert module.mcp is None
    with pytest.raises(RuntimeError, match="optional 'mcp' package"):
        module.run()


def test_server_exposes_typed_mcp_surface_without_options_json():
    server_source = inspect.getsource(server_module)
    assert "options_json" not in server_source
    assert "method://tooling_surface_v1" in SCHEMAS



def test_mcp_initialize_reports_package_version():
    from mcp import Client

    async def initialize():
        async with Client(mcp) as client:
            return client.server_info

    server_info = asyncio.run(initialize())
    assert server_info.name == "proteomics-mcp"
    assert server_info.version == __version__


def test_actionable_tool_errors_reach_the_client(monkeypatch):
    from mcp import Client

    from proteomics_mcp.predictors.ms2pip import LIBOMP_HINT, MissingOpenMPError

    def missing_libomp(*args, **kwargs):
        raise MissingOpenMPError(LIBOMP_HINT)

    monkeypatch.setattr(server_module, "predict_fragments_tool", missing_libomp)

    async def call():
        async with Client(mcp) as client:
            return await client.call_tool("predict_fragments", {"peptidoform": "PEPTIDE", "charge": 2})

    result = asyncio.run(call())
    assert result.is_error
    assert "brew install libomp" in result.content[0].text


def _call_tool(name: str, arguments: dict):
    from mcp import Client

    async def call():
        async with Client(mcp) as client:
            return await client.call_tool(name, arguments)

    return asyncio.run(call())


@pytest.mark.parametrize(
    ("name", "arguments", "expected"),
    [
        ("load_candidates", {"path": "does/not/exist.csv", "adapter": "fragpipe"}, "FileNotFoundError"),
        ("load_candidates", {"path": "does/not/exist.csv", "adapter": "not_an_adapter"}, "ValueError: Unsupported candidate adapter: not_an_adapter"),
        ("score_spectrum", {"candidate": {"peptidoform": "PEPTIDE"}, "observed": {"mz": [100.0], "intensity": [1.0]}},
         "ValueError: Invalid 'candidate'"),
        ("annotate_observed_spectrum", {"observed": {"mz": "not a list"}, "peptidoform": "PEPTIDE"}, "ValueError: Invalid 'observed'"),
    ],
    ids=["missing-file", "unknown-adapter", "incomplete-candidate", "invalid-spectrum"],
)
def test_expected_tool_errors_reach_the_client(name, arguments, expected):
    """Errors the agent can correct are shown with their message, not only "Error executing tool" (#27)."""
    result = _call_tool(name, arguments)
    assert result.is_error
    assert expected in result.content[0].text


def test_lookup_errors_reach_the_client_without_repr_quotes(monkeypatch):
    def unknown_field(**_):
        raise KeyError("Unknown field 'spectral_angel'.")

    monkeypatch.setattr(server_module, "summarize_numeric_field_tool", unknown_field)
    result = _call_tool("summarize_numeric_field", {"path": "x.csv", "field": "spectral_angel"})
    assert result.is_error
    assert result.content[0].text.endswith("KeyError: Unknown field 'spectral_angel'.")


def test_unexpected_tool_errors_stay_masked(monkeypatch):
    def crash(*args, **kwargs):
        raise RuntimeError("internal state /home/someone/secret")

    monkeypatch.setattr(server_module, "predict_fragments_tool", crash)
    result = _call_tool("predict_fragments", {"peptidoform": "PEPTIDE", "charge": 2})
    assert result.is_error
    # mcp 2.0.x passed the message through; the package requires 2.1, the first release that hides it (#164)
    assert "secret" not in result.content[0].text
    assert "Error executing tool predict_fragments" in result.content[0].text


def test_mcp_tools_have_directory_annotations():
    tools = {tool.name: tool for tool in mcp._tool_manager.list_tools()}
    assert len(tools) >= 50
    for tool in tools.values():
        assert tool.annotations is not None
        assert tool.annotations.title
        assert tool.annotations.read_only_hint is not None
        assert tool.annotations.destructive_hint is not None
    assert tools["inspect_psm_table"].annotations.read_only_hint is True
    assert tools["mgf_to_feather"].annotations.destructive_hint is True
    assert tools["retrieve_usi_spectrum"].annotations.open_world_hint is True

    if not MANIFEST_PATH.is_file():
        pytest.skip("manifest.json is only available in a source tree or sdist")
    manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    assert {tool["name"] for tool in manifest["tools"]} == set(tools)
