#!/usr/bin/env python3
"""Check that every copy of the release version agrees with pyproject.toml.

pyproject.toml ``[project].version`` is the single source of truth. This script
fails if any of the following disagree with it:

- ``manifest.json`` (Claude Desktop extension) ``version``;
- ``CITATION.cff`` ``version``;
- ``install.json``, if it carries a ``version`` anywhere;
- the git tag, when given with ``--tag`` (expected ``v<version>``);
- the installed distribution, when ``--installed`` is given;

and it fails if ``docs/releases/v<version>.md`` (GitHub release notes) is missing, or if ``CHANGELOG.md`` has no
``## [<version>] - YYYY-MM-DD`` heading. With ``--tag`` that date must also be the day the tag is checked (UTC, give or
take a day for time zones): the release prep sets it, and it has to be updated if the tag comes later.

Usage:
    python scripts/check_version.py
    python scripts/check_version.py --tag v0.1.9
    python scripts/check_version.py --installed
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import re
import sys
import tomllib
from importlib.metadata import PackageNotFoundError, version as installed_version
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
DISTRIBUTION = "proteomics-mcp"


def _pyproject_version(root: Path) -> str:
    project = tomllib.loads((root / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    if "version" in project.get("dynamic", []) or "version" not in project:
        raise SystemExit("pyproject.toml must declare a static [project].version")
    return project["version"]


def _citation_version(root: Path) -> str | None:
    text = (root / "CITATION.cff").read_text(encoding="utf-8")
    match = re.search(r"""^version:\s*["']?([^"'\s#]+)["']?\s*(?:#.*)?$""", text, flags=re.MULTILINE)
    return match.group(1) if match else None


def _find_versions(value: Any, path: str = "") -> list[tuple[str, str]]:
    found: list[tuple[str, str]] = []
    if isinstance(value, dict):
        for key, item in value.items():
            child = f"{path}.{key}" if path else key
            if key == "version" and isinstance(item, str):
                found.append((child, item))
            else:
                found.extend(_find_versions(item, child))
    elif isinstance(value, list):
        for index, item in enumerate(value):
            found.extend(_find_versions(item, f"{path}[{index}]"))
    return found


def _changelog_date(root: Path, version: str) -> dt.date | None:
    """The date of the version's CHANGELOG.md heading, "## [X.Y.Z] - YYYY-MM-DD"."""
    text = (root / "CHANGELOG.md").read_text(encoding="utf-8")
    match = re.search(rf"^## \[{re.escape(version)}\] - (\d{{4}}-\d{{2}}-\d{{2}})\s*$", text, flags=re.MULTILINE)
    return dt.date.fromisoformat(match.group(1)) if match else None


def check(
    root: Path, tag: str | None = None, installed: bool = False, today: dt.date | None = None
) -> tuple[str, list[str], list[str]]:
    expected = _pyproject_version(root)
    problems: list[str] = []
    checked: list[str] = [f"pyproject.toml: {expected}"]

    manifest_version = json.loads((root / "manifest.json").read_text(encoding="utf-8")).get("version")
    checked.append(f"manifest.json: {manifest_version}")
    if manifest_version != expected:
        problems.append(f"manifest.json version is {manifest_version!r}, expected {expected!r}")

    citation_version = _citation_version(root)
    checked.append(f"CITATION.cff: {citation_version}")
    if citation_version != expected:
        problems.append(f"CITATION.cff version is {citation_version!r}, expected {expected!r}")

    install_json = root / "install.json"
    if install_json.is_file():
        install_versions = _find_versions(json.loads(install_json.read_text(encoding="utf-8")))
        checked.append(f"install.json: {', '.join(f'{k}={v}' for k, v in install_versions) or 'not versioned'}")
        for key, value in install_versions:
            if value != expected:
                problems.append(f"install.json {key} is {value!r}, expected {expected!r}")

    if tag is not None:
        checked.append(f"git tag: {tag}")
        if tag != f"v{expected}":
            problems.append(f"tag is {tag!r}, expected 'v{expected}'")

    if installed:
        try:
            installed_value = installed_version(DISTRIBUTION)
        except PackageNotFoundError:
            installed_value = None
        checked.append(f"installed {DISTRIBUTION}: {installed_value}")
        if installed_value != expected:
            problems.append(f"installed {DISTRIBUTION} is {installed_value!r}, expected {expected!r}")

    released = _changelog_date(root, expected)
    checked.append(f"CHANGELOG.md: {released or 'no dated heading'}")
    if released is None:
        problems.append(f"CHANGELOG.md has no '## [{expected}] - YYYY-MM-DD' heading")
    elif tag is not None:
        today = today or dt.datetime.now(dt.timezone.utc).date()
        # A day either way for time zones: the tag is pushed from a local clock, the check runs in UTC
        if abs((today - released).days) > 1:
            problems.append(f"CHANGELOG.md dates {expected} {released}, but it is tagged on {today}: set the release date before tagging")

    notes = root / "docs" / "releases" / f"v{expected}.md"
    checked.append(f"docs/releases/v{expected}.md: {'present' if notes.is_file() else 'missing'}")
    if not notes.is_file():
        problems.append(f"release notes missing: docs/releases/v{expected}.md")

    return expected, checked, problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--tag", help="git tag to check, e.g. v0.1.9")
    parser.add_argument("--installed", action="store_true", help=f"also check the installed {DISTRIBUTION} version")
    parser.add_argument("--root", type=Path, default=ROOT, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)

    expected, checked, problems = check(args.root, tag=args.tag, installed=args.installed)
    for line in checked:
        print(f"  {line}")
    if problems:
        print(f"Version check FAILED for {expected}:")
        for problem in problems:
            print(f"  - {problem}")
        return 1
    print(f"Version check OK: {expected}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
