import json
import tomllib
from pathlib import Path

from proteomics_mcp import __version__


ROOT = Path(__file__).resolve().parents[1]


def test_package_and_extension_versions_match():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    extension = json.loads((ROOT / "manifest.json").read_text(encoding="utf-8"))

    assert project["version"] == __version__ == extension["version"]
    assert project["name"] == extension["name"]
    assert (ROOT / "docs" / "releases" / f"v{__version__}.md").is_file()
