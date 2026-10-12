#!/usr/bin/env python3
"""Smoke-test the installed MCP server over stdio.

Launches ``proteomics-mcp serve`` with the ``mcp`` client SDK and checks that:

- ``initialize`` succeeds and reports the expected server version;
- the core tools are listed (and, with ``--manifest``, that the tool names match
  a Claude Desktop ``manifest.json`` exactly);
- one offline tool call works: ``annotate_observed_spectrum`` on the bundled
  fixture ``tests/fixtures/smoke_observed_spectrum.json``.

No network access and no real data are needed. Run it with the Python of the
environment where the package is installed, e.g.

    .venv/bin/python scripts/smoke_mcp.py --manifest manifest.json

``--server-command`` launches another command instead, e.g. the Claude Desktop
extension entry point from manifest.json:

    python scripts/smoke_mcp.py --expect-version 0.1.8 \
        --server-command "uv run --directory . --no-dev extension_server.py"
"""

from __future__ import annotations

import argparse
import asyncio
import json
import shlex
import shutil
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_FIXTURE = ROOT / "tests" / "fixtures" / "smoke_observed_spectrum.json"
CORE_TOOLS = {
    "inspect_psm_table",
    "render_peptide_table",
    "score_spectrum",
    "build_prediction_batch",
    "describe_predictor_model",
    "annotate_observed_spectrum",
}
TIMEOUT_SECONDS = 600  # allows a first-run `uv run` environment sync with --server-command


def _server_command() -> str:
    """Prefer the console script next to the running interpreter (same environment)."""
    bindir = Path(sys.executable).parent
    for name in ("proteomics-mcp", "proteomics-mcp.exe"):
        candidate = bindir / name
        if candidate.is_file():
            return str(candidate)
    found = shutil.which("proteomics-mcp")
    if found is None:
        raise SystemExit("proteomics-mcp console script not found next to the interpreter or on PATH")
    return found


async def _run(args: argparse.Namespace) -> list[str]:
    problems: list[str] = []
    fixture = json.loads(args.fixture.read_text(encoding="utf-8"))
    if args.server_command:
        command, *command_args = shlex.split(args.server_command)
    else:
        command, command_args = _server_command(), ["serve"]
    print(f"server command: {' '.join([command, *command_args])}")
    params = StdioServerParameters(command=command, args=command_args)

    async with stdio_client(params) as (read, write):
        async with ClientSession(read, write) as session:
            init = await session.initialize()
            print(f"initialize: {init.server_info.name} {init.server_info.version}")
            if init.server_info.name != "proteomics-mcp":
                problems.append(f"server name is {init.server_info.name!r}, expected 'proteomics-mcp'")
            if not args.skip_version_check and init.server_info.version != args.expect_version:
                problems.append(f"server version is {init.server_info.version!r}, expected {args.expect_version!r}")

            tools = {tool.name for tool in (await session.list_tools()).tools}
            print(f"tools: {len(tools)}")
            missing = sorted(CORE_TOOLS - tools)
            if missing:
                problems.append(f"core tools missing: {', '.join(missing)}")
            if args.manifest is not None:
                manifest_tools = {tool["name"] for tool in json.loads(args.manifest.read_text(encoding="utf-8"))["tools"]}
                if manifest_tools != tools:
                    problems.append(
                        f"tool names differ from {args.manifest.name}: "
                        f"server-only={sorted(tools - manifest_tools)} manifest-only={sorted(manifest_tools - tools)}"
                    )

            result = await session.call_tool(
                "annotate_observed_spectrum",
                {"observed": fixture["observed"], "peptidoform": fixture["peptidoform"]},
            )
            content = result.structured_content or {}
            content = content.get("result", content)
            expected = fixture["expected"]
            print(
                f"annotate_observed_spectrum: is_error={result.is_error} "
                f"annotated_peak_count={content.get('annotated_peak_count')} ion_type_counts={content.get('ion_type_counts')}"
            )
            if result.is_error:
                problems.append("annotate_observed_spectrum returned an error")
            for key, value in expected.items():
                if content.get(key) != value:
                    problems.append(f"annotate_observed_spectrum {key} is {content.get(key)!r}, expected {value!r}")
    return problems


def main(argv: list[str] | None = None) -> int:
    try:
        installed = version("proteomics-mcp")
    except PackageNotFoundError:
        installed = None
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--expect-version", default=installed, help="expected server version (default: installed version)")
    parser.add_argument("--skip-version-check", action="store_true", help="do not compare the server version")
    parser.add_argument("--manifest", type=Path, help="manifest.json whose tool names must match the server exactly")
    parser.add_argument("--server-command", help="command that starts the MCP server (default: proteomics-mcp serve)")
    parser.add_argument("--fixture", type=Path, default=DEFAULT_FIXTURE, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if not args.skip_version_check and args.expect_version is None:
        parser.error("proteomics-mcp is not installed; pass --expect-version or --skip-version-check")

    problems = asyncio.run(asyncio.wait_for(_run(args), TIMEOUT_SECONDS))
    if problems:
        print("MCP smoke test FAILED:")
        for problem in problems:
            print(f"  - {problem}")
        return 1
    print("MCP smoke test OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
