#!/usr/bin/env python3
"""Check that a built wheel contains every package module and declared resource.

The package is built with ``include-package-data = false``, so a runtime file
under ``src/proteomics_mcp`` reaches the wheel only if it is listed in
``[tool.setuptools.package-data]``. This script fails if:

- a ``.py`` module under ``src/proteomics_mcp`` is missing from the wheel;
- a file matched by ``[tool.setuptools.package-data]`` is missing from the wheel;
- a non-``.py`` file under ``src/proteomics_mcp`` is not declared as package data
  (it would be silently left out of the wheel);
- the ``proteomics-mcp`` console script is not declared in the wheel metadata.

Usage: python scripts/check_wheel.py dist/proteomics_mcp-<version>-py3-none-any.whl
"""

from __future__ import annotations

import argparse
import sys
import tomllib
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
PACKAGE = "proteomics_mcp"
CONSOLE_SCRIPT = "proteomics-mcp = proteomics_mcp.cli:main"


def _source_files() -> list[Path]:
    package_root = SRC / PACKAGE
    return sorted(
        path
        for path in package_root.rglob("*")
        if path.is_file() and "__pycache__" not in path.parts and path.suffix not in {".pyc", ".pyo"}
    )


def _declared_package_data() -> set[str]:
    pyproject = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    package_data = pyproject.get("tool", {}).get("setuptools", {}).get("package-data", {})
    declared: set[str] = set()
    for package, patterns in package_data.items():
        package_root = SRC / Path(*package.split("."))
        for pattern in patterns:
            for path in package_root.glob(pattern):
                if path.is_file():
                    declared.add(path.relative_to(SRC).as_posix())
    return declared


def check(wheel_path: Path) -> list[str]:
    with zipfile.ZipFile(wheel_path) as wheel:
        names = set(wheel.namelist())
        entry_points = next((name for name in names if name.endswith(".dist-info/entry_points.txt")), None)
        entry_points_text = wheel.read(entry_points).decode("utf-8") if entry_points else ""

    declared = _declared_package_data()
    problems: list[str] = []
    for path in _source_files():
        relative = path.relative_to(SRC).as_posix()
        if path.suffix == ".py":
            if relative not in names:
                problems.append(f"module missing from wheel: {relative}")
        elif relative not in declared:
            problems.append(f"resource not declared in [tool.setuptools.package-data]: {relative}")
    for relative in sorted(declared):
        if relative not in names:
            problems.append(f"declared resource missing from wheel: {relative}")
    if CONSOLE_SCRIPT not in {line.strip() for line in entry_points_text.splitlines()}:
        problems.append(f"console script not declared in wheel metadata: {CONSOLE_SCRIPT}")
    return problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("wheel", type=Path, help="path to the built .whl file")
    args = parser.parse_args(argv)

    problems = check(args.wheel)
    if problems:
        print(f"{args.wheel.name}: {len(problems)} problem(s)")
        for problem in problems:
            print(f"  - {problem}")
        return 1
    modules = sum(1 for path in _source_files() if path.suffix == ".py")
    print(f"{args.wheel.name}: OK ({modules} modules, {len(_declared_package_data())} declared resources)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
