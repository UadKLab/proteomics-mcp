from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

# Set PROTEOMICS_MCP_TEST_INSTALLED=1 to require that the tests import the installed
# package (site-packages) rather than the source tree, e.g. when testing a built wheel.
REQUIRE_INSTALLED = os.environ.get("PROTEOMICS_MCP_TEST_INSTALLED") == "1"


def _package_location() -> Path:
    import proteomics_mcp

    return Path(proteomics_mcp.__file__).resolve().parent


def pytest_report_header(config):
    return f"proteomics_mcp imported from: {_package_location()}"


def pytest_sessionstart(session):
    location = _package_location()
    if REQUIRE_INSTALLED and "site-packages" not in location.parts:
        raise pytest.UsageError(
            f"PROTEOMICS_MCP_TEST_INSTALLED=1 but proteomics_mcp was imported from {location}, not site-packages"
        )


def pytest_collection_modifyitems(config, items):
    """Deselect end-to-end tests unless the -m expression asks for them (e.g. -m e2e)."""
    if "e2e" in (config.getoption("markexpr") or ""):
        return
    selected = [item for item in items if item.get_closest_marker("e2e") is None]
    deselected = [item for item in items if item.get_closest_marker("e2e") is not None]
    if deselected:
        config.hook.pytest_deselected(items=deselected)
        items[:] = selected


@pytest.fixture
def ms2pip_runtime():
    """Skip on macOS when ms2pip cannot load because xgboost cannot find libomp.

    xgboost raises XGBoostError (not ImportError) when libomp.dylib is missing.
    Any other import failure is re-raised so real breakage is not hidden.
    """
    try:
        import ms2pip  # noqa: F401
    except ImportError:
        return
    except Exception as exc:
        if sys.platform == "darwin" and "libomp" in str(exc):
            pytest.skip("ms2pip/xgboost cannot load libomp on macOS; install it with `brew install libomp`")
        raise
