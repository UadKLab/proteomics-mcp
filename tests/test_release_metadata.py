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


def _check_version_module():
    import importlib.util

    spec = importlib.util.spec_from_file_location("check_version", ROOT / "scripts" / "check_version.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _release_tree(root: Path, heading: str) -> Path:
    (root / "docs" / "releases").mkdir(parents=True)
    (root / "pyproject.toml").write_text('[project]\nname = "proteomics-mcp"\nversion = "9.9.9"\n', encoding="utf-8")
    (root / "manifest.json").write_text('{"version": "9.9.9"}', encoding="utf-8")
    (root / "CITATION.cff").write_text('version: "9.9.9"\n', encoding="utf-8")
    (root / "docs" / "releases" / "v9.9.9.md").write_text("notes\n", encoding="utf-8")
    (root / "CHANGELOG.md").write_text(f"# Changelog\n\n## [Unreleased]\n\n{heading}\n\n### Fixed\n", encoding="utf-8")
    return root


@pytest.mark.parametrize(
    ("heading", "tag", "today", "problem"),
    [
        ("## [9.9.9] - 2026-10-26", None, None, None),
        ("## [9.9.9] - 2026-10-26", "v9.9.9", "2026-10-26", None),
        ("## [9.9.9] - 2026-10-26", "v9.9.9", "2026-10-27", None),  # a time zone apart
        ("## [9.9.9] - 2026-10-26", "v9.9.9", "2026-10-30", "dates 9.9.9 2026-10-26, but it is tagged on 2026-10-30"),
        ("## [9.9.9]", None, None, "no '## [9.9.9] - YYYY-MM-DD' heading"),
        ("## [9.9.8] - 2026-10-26", None, None, "no '## [9.9.9] - YYYY-MM-DD' heading"),
    ],
)
def test_the_changelog_dates_the_release_on_the_tag_day(tmp_path, heading, tag, today, problem):
    import datetime as dt

    check_version = _check_version_module()
    root = _release_tree(tmp_path, heading)

    _, _, problems = check_version.check(root, tag=tag, today=dt.date.fromisoformat(today) if today else None)

    assert [p for p in problems if problem and problem in p] if problem else problems == []
    if problem:
        assert len(problems) == 1


def test_the_repository_changelog_dates_the_current_version():
    check_version = _check_version_module()

    assert check_version._changelog_date(ROOT, __version__) is not None
