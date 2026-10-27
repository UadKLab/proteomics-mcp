#!/usr/bin/env python3
"""Write a shields.io endpoint badge for the total coverage in a ``coverage json`` report (#166).

The coverage test leg writes it on every push to main, the docs workflow publishes it with the documentation site,
and the README badge reads it through https://img.shields.io/endpoint.

Usage:
    python scripts/coverage_badge.py coverage.json coverage-badge.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any

# (lowest percentage, colour), best first; below the last one the badge is red
COLOURS = ((90.0, "brightgreen"), (80.0, "green"), (70.0, "yellowgreen"), (60.0, "yellow"), (50.0, "orange"))


def badge(report: dict[str, Any]) -> dict[str, Any]:
    """The badge for a report's combined line and branch coverage, the figure fail_under is checked against."""
    shown = round(float(report["totals"]["percent_covered"]), 1)  # the colour follows the figure the badge shows
    colour = next((name for lowest, name in COLOURS if shown >= lowest), "red")
    return {"schemaVersion": 1, "label": "coverage", "message": f"{shown:.1f}%", "color": colour}


def main(argv: list[str] | None = None) -> int:
    args = sys.argv[1:] if argv is None else argv
    if len(args) != 2:
        print(__doc__.strip().splitlines()[-1].strip(), file=sys.stderr)
        return 2
    report = json.loads(Path(args[0]).read_text(encoding="utf-8"))
    Path(args[1]).write_text(json.dumps(badge(report)) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main())
