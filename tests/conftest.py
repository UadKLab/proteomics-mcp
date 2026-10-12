from __future__ import annotations

import sys

import pytest


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
