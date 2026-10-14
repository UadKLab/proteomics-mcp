"""Shared fixtures for end-to-end checks of the installed package (#30).

Run with the interpreter of an environment where proteomics-mcp was installed from
PyPI (never from this source tree):

    python -m pytest -o pythonpath= -m e2e tests/e2e

Environment variables:
    PROTEOMICS_MCP_E2E_CACHE   download cache (default: ~/.cache/proteomics-mcp-e2e)
    PROTEOMICS_MCP_E2E_REPORT  write a JSON report of all recorded checks to this path
"""

from __future__ import annotations

import asyncio
import hashlib
import json
import os
import shutil
import subprocess
import sys
import urllib.request
import zipfile
from importlib.metadata import version
from pathlib import Path
from typing import Any

import pytest

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_DIR = REPO_ROOT / "examples"

ZENODO_RECORD = "23295902"
ZENODO_URL = f"https://zenodo.org/records/{ZENODO_RECORD}/files/supplement.zip?download=1"
ZENODO_MD5 = "1445e7af348b12c703f477d50fb01130"

_RECORDS: list[dict[str, Any]] = []


def _package_location() -> Path:
    import proteomics_mcp

    return Path(proteomics_mcp.__file__).resolve().parent


@pytest.fixture(scope="session", autouse=True)
def installed_package() -> dict[str, str]:
    """End-to-end checks must exercise the installed package, not the source tree."""
    location = _package_location()
    if "site-packages" not in location.parts:
        pytest.exit(f"e2e tests require an installed proteomics-mcp; imported from {location}", returncode=4)
    info = {"proteomics-mcp": version("proteomics-mcp"), "mcp": version("mcp"), "ms2pip": version("ms2pip"),
            "python": sys.version.split()[0], "platform": sys.platform, "location": str(location)}
    record("environment", "installed package", "info", info)
    return info


def record(section: str, item: str, status: str, details: dict[str, Any] | None = None) -> None:
    """Collect a result line for the per-batch report posted on #30."""
    _RECORDS.append({"section": section, "item": item, "status": status, "details": details or {}})


def pytest_sessionfinish(session, exitstatus):
    report = os.environ.get("PROTEOMICS_MCP_E2E_REPORT")
    if report and _RECORDS:
        Path(report).write_text(json.dumps(_RECORDS, indent=1, default=str), encoding="utf-8")


def _cache_dir() -> Path:
    root = Path(os.environ.get("PROTEOMICS_MCP_E2E_CACHE", Path.home() / ".cache" / "proteomics-mcp-e2e"))
    root.mkdir(parents=True, exist_ok=True)
    return root


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


@pytest.fixture(scope="session")
def zenodo() -> Path:
    """The manuscript's Zenodo supplement, downloaded once, md5-checked and extracted."""
    archive = _cache_dir() / f"zenodo-{ZENODO_RECORD}-supplement.zip"
    if not archive.is_file() or _md5(archive) != ZENODO_MD5:
        partial = archive.with_suffix(".part")
        with urllib.request.urlopen(ZENODO_URL, timeout=300) as response, partial.open("wb") as handle:
            shutil.copyfileobj(response, handle)
        partial.replace(archive)
    checksum = _md5(archive)
    assert checksum == ZENODO_MD5, f"Zenodo supplement md5 {checksum} != {ZENODO_MD5}"
    target = _cache_dir() / f"zenodo-{ZENODO_RECORD}"
    if not (target / "supplement" / "raw_file_inventory.tsv").is_file():
        with zipfile.ZipFile(archive) as bundle:
            members = [m for m in bundle.namelist() if not m.startswith("__MACOSX/") and not m.endswith(".DS_Store")]
            bundle.extractall(target, members=members)
    record("inputs", "Zenodo supplement.zip", "info", {"record": ZENODO_RECORD, "md5": checksum,
                                                         "sha256": hashlib.sha256(archive.read_bytes()).hexdigest()})
    return target / "supplement"


@pytest.fixture(scope="session")
def examples_dir() -> Path:
    if not EXAMPLES_DIR.is_dir():
        pytest.skip("examples/ (reference payloads) is only available in a repository checkout")
    return EXAMPLES_DIR


def _structured(result: Any) -> Any:
    content = result.structured_content
    if isinstance(content, dict) and set(content) == {"result"}:
        return content["result"]
    return content


def call_tool(name: str, arguments: dict[str, Any]) -> tuple[bool, Any]:
    """Call an MCP tool of the installed server through the in-memory MCP client.

    Returns (is_error, structured result or error text).
    """
    from mcp import Client

    from proteomics_mcp.server import mcp

    async def _call():
        async with Client(mcp) as client:
            return await client.call_tool(name, arguments)

    result = asyncio.run(_call())
    if result.is_error:
        return True, " ".join(getattr(block, "text", "") for block in result.content)
    return False, _structured(result)


def run_cli(*args: str, timeout: int = 600) -> subprocess.CompletedProcess[str]:
    """Run the installed console script that sits next to the current interpreter."""
    bindir = Path(sys.executable).parent
    script = next((bindir / name for name in ("proteomics-mcp", "proteomics-mcp.exe") if (bindir / name).is_file()), None)
    assert script is not None, f"proteomics-mcp console script not found in {bindir}"
    return subprocess.run([str(script), *args], capture_output=True, text=True, timeout=timeout)


@pytest.fixture(scope="session")
def tool():
    return call_tool


@pytest.fixture(scope="session")
def cli():
    return run_cli


@pytest.fixture(scope="session")
def report():
    return record
