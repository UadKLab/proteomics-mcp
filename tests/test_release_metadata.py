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


def test_pyproject_authors_match_citation():
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    citation = (ROOT / "CITATION.cff").read_text(encoding="utf-8")
    cited = [
        f"{given} {family}"
        for given, family in re.findall(r'- given-names: "([^"]+)"\n\s+family-names: "([^"]+)"', citation)
    ]

    assert cited, "CITATION.cff lists no authors"
    assert [author["name"] for author in project["authors"]] == cited


def test_readme_links_work_on_pypi():
    # The README is also the PyPI project description, where relative links break (#98).
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    targets = re.findall(r"\]\(([^)\s]+)\)|(?:href|src)=\"([^\"]+)\"", readme)
    links = [markdown or html for markdown, html in targets]
    assert links
    assert [link for link in links if not link.startswith("https://")] == []
    # The documentation pages are linked on the documentation site, not as Markdown on GitHub.
    assert [link for link in links if "/blob/main/docs/" in link] == []
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    assert project["urls"]["Documentation"] == "https://uadklab.github.io/proteomics-mcp/"
