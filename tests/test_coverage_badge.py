"""The README coverage badge file (#166)."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "coverage_badge.py"

pytestmark = pytest.mark.skipif(not SCRIPT.is_file(), reason="scripts/ is only available in a source tree or sdist")


def _module():
    spec = importlib.util.spec_from_file_location("coverage_badge", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize(
    ("percent", "message", "colour"),
    [(88.65, "88.7%", "green"), (90.0, "90.0%", "brightgreen"), (79.99, "80.0%", "green"), (79.94, "79.9%", "yellowgreen"), (42.0, "42.0%", "red")],
)
def test_badge_shows_the_combined_coverage(percent, message, colour):
    badge = _module().badge({"totals": {"percent_covered": percent}})

    assert badge == {"schemaVersion": 1, "label": "coverage", "message": message, "color": colour}


def test_the_script_writes_a_shields_endpoint_file(tmp_path):
    report = tmp_path / "coverage.json"
    report.write_text(json.dumps({"totals": {"percent_covered": 91.2, "covered_lines": 1}}), encoding="utf-8")
    output = tmp_path / "coverage-badge.json"

    assert _module().main([str(report), str(output)]) == 0
    assert json.loads(output.read_text(encoding="utf-8"))["message"] == "91.2%"
    assert _module().main([str(report)]) == 2
