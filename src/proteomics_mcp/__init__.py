"""Proteomics MCP package."""

from importlib.metadata import PackageNotFoundError, version

__all__ = ["__version__"]

# pyproject.toml is the single source of truth; read the installed distribution's version.
try:
    __version__ = version("proteomics-mcp")
except PackageNotFoundError:  # imported from a source tree without being installed
    __version__ = "0+unknown"
