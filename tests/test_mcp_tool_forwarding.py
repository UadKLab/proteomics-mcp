"""Every MCP tool passes its arguments on, and every resource is served (#166).

Most behaviour is tested through the *_tool functions; these tests check the MCP layer in front of them: a wrapper that
dropped or renamed an argument would make the agent's setting silently ineffective.
"""

from __future__ import annotations

import asyncio
import hashlib
import inspect
from pathlib import Path
from typing import Any

import pytest

import proteomics_mcp.server as server_module
from proteomics_mcp.resources.docs import SCHEMAS
from proteomics_mcp.server import mcp

pytestmark = pytest.mark.skipif(mcp is None, reason="the mcp package is not installed")

TOOLS = mcp._tool_manager.list_tools() if mcp is not None else []
TOOL_FUNCTIONS = [name for name in dir(server_module) if name.endswith("_tool") and callable(getattr(server_module, name))]


def _value(name: str, schema: dict[str, Any]) -> Any:
    """A value of the schema's type that names its argument, so it can be found again among the forwarded ones."""
    options = [option for option in schema.get("anyOf", [schema]) if option.get("type") != "null"]
    option = options[0] if options else schema
    if "enum" in option:
        return option["enum"][0]
    kind = option.get("type")
    if kind == "array":
        return [_value(name, option.get("items", {}))]
    number = int(hashlib.sha1(name.encode()).hexdigest()[:4], 16) + 1  # distinct per argument
    return {"string": f"value of {name}", "integer": number, "number": number + 0.5, "boolean": True, "object": {"argument": name}}.get(
        kind, f"value of {name}"
    )


def _call(name: str, arguments: dict[str, Any]):
    from mcp import Client

    async def call():
        async with Client(mcp) as client:
            return await client.call_tool(name, arguments)

    return asyncio.run(call())


@pytest.mark.parametrize("tool", TOOLS, ids=[tool.name for tool in TOOLS])
def test_every_argument_reaches_the_tool_function(monkeypatch, tool):
    calls: list[tuple[str, tuple, dict]] = []
    returns_list = "list" in str(inspect.signature(tool.fn).return_annotation)

    def recorder(function_name):
        def record(*args, **kwargs):
            calls.append((function_name, args, kwargs))
            return [] if returns_list else {}

        return record

    for function_name in TOOL_FUNCTIONS:
        monkeypatch.setattr(server_module, function_name, recorder(function_name))
    arguments = {name: _value(name, schema) for name, schema in tool.parameters["properties"].items()}

    result = _call(tool.name, arguments)

    assert not result.is_error, result.content[0].text
    assert len(calls) == 1, f"{tool.name} called {[name for name, _, _ in calls]}"
    _, args, kwargs = calls[0]
    forwarded = [*args, *kwargs.values()]
    # some wrappers combine arguments, e.g. range_min and range_max into value_range
    forwarded += [item for value in forwarded if isinstance(value, tuple) for item in value]
    assert [name for name, value in arguments.items() if value not in forwarded] == []


def test_relative_output_paths_resolve_against_the_output_directory(monkeypatch, tmp_path):
    monkeypatch.setattr(server_module, "_OUTPUT_BASE", str(tmp_path))

    assert server_module._resolve("plots/figure.html") == str(tmp_path / "plots" / "figure.html")
    assert server_module._resolve(str(tmp_path / "absolute.html")) == str(tmp_path / "absolute.html")
    assert server_module._resolve(None) is None
    monkeypatch.setattr(server_module, "_OUTPUT_BASE", None)
    assert server_module._resolve("plots/figure.html") == "plots/figure.html"


def test_every_resource_is_served():
    from mcp import Client

    async def read_all():
        async with Client(mcp) as client:
            listed = await client.list_resources()
            return {str(resource.uri): (await client.read_resource(resource.uri)).contents[0].text for resource in listed.resources}

    served = asyncio.run(read_all())

    assert served == {uri: SCHEMAS[uri] for uri in served}
    assert set(served) == {uri for uri in SCHEMAS if uri.split("://")[0] in {"method", "schema", "config"}}


def test_the_server_refuses_to_run_without_mcp(monkeypatch):
    monkeypatch.setattr(server_module, "mcp", None)

    with pytest.raises(RuntimeError, match="requires the optional 'mcp' package"):
        server_module.run()


def test_output_paths_are_expanded(monkeypatch):
    monkeypatch.setattr(server_module, "_OUTPUT_BASE", "~/pmcp-out")

    assert server_module._resolve("x.html") == str(Path("~/pmcp-out").expanduser() / "x.html")
