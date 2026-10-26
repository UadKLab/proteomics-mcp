"""Every PSM Comet and X! Tandem reported in validation/engine_scores/golden/ is reproduced (#51).

The golden file holds the engines' own scores, from the pinned containers; validate.py --run-engines regenerates it,
and .github/workflows/engine-validation.yml checks that it still does. This test needs neither Docker nor network.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

VALIDATION = Path(__file__).resolve().parents[1] / "validation" / "engine_scores"
if not (VALIDATION / "golden" / "psms.tsv").is_file() or not (VALIDATION.parents[1] / "examples").is_dir():
    pytest.skip("validation/ and examples/ are only available in a source checkout", allow_module_level=True)

spec = importlib.util.spec_from_file_location("engine_validation", VALIDATION / "validate.py")
validate = importlib.util.module_from_spec(spec)
sys.modules["engine_validation"] = validate
spec.loader.exec_module(validate)


@pytest.fixture(scope="module")
def rescored(tmp_path_factory):
    return validate.rescore(validate.read_golden(), tmp_path_factory.mktemp("engine_validation"))


def test_every_engine_psm_is_reproduced(rescored):
    failed = [(row["run"], row["title"], row["peptidoform"], row["engine_score"], row["ours"]) for row in rescored if not row["ok"]]
    assert not failed, failed[:10]


def test_the_validation_covers_both_engines_presets_and_wrong_candidates(rescored):
    summary = validate.summarize(rescored)
    assert set(summary) == {run.name for run in validate.RUNS}
    assert all(item["psms"] > 50 for item in summary.values())
    assert summary["comet_high_res"]["decoy_psms"] > 100 and summary["comet_low_res"]["decoy_psms"] > 100
    assert summary["xtandem_reversed"]["decoy_psms"] == summary["xtandem_reversed"]["psms"]


def test_committed_inputs_are_the_ones_the_golden_file_was_made_from():
    import hashlib
    import json

    recorded = json.loads((VALIDATION / "golden" / "versions.json").read_text(encoding="utf-8"))["inputs_sha256"]
    for name, digest in recorded.items():
        path = next(p for p in (VALIDATION / "data" / name, VALIDATION / "params" / name) if p.is_file())
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, name
