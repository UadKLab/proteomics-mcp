import json
import re
import tomllib
from pathlib import Path

import pytest

from proteomics_mcp import __version__


ROOT = Path(__file__).resolve().parents[1]

pytestmark = pytest.mark.skipif(
    not (ROOT / "pyproject.toml").is_file(),
    reason="release metadata files are only available in a source tree or sdist",
)


def test_package_and_extension_versions_match():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    extension = json.loads((ROOT / "manifest.json").read_text(encoding="utf-8"))

    assert project["version"] == __version__ == extension["version"]
    assert project["name"] == extension["name"]
    assert (ROOT / "docs" / "releases" / f"v{__version__}.md").is_file()


def test_citation_version_matches_package():
    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    match = re.search(r"""^version:\s*["']?([^"'\s#]+)""", citation, flags=re.MULTILINE)

    assert match is not None, "CITATION.cff has no top-level version"
    assert match.group(1) == __version__
