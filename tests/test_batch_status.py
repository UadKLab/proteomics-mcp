"""Batch status while a batch runs and after it was interrupted (#46)."""

from __future__ import annotations

import json
import threading
import time
from pathlib import Path

import pytest

import proteomics_mcp.tools.workflows as workflows
from proteomics_mcp.predictors.mock import mock_predicted_spectrum
from proteomics_mcp.tools.workflows import (
    build_prediction_batch_tool,
    inspect_prediction_batch_status_tool,
    score_from_prediction_batch_tool,
)


@pytest.fixture
def candidates(tmp_path):
    path = tmp_path / "preds.csv"
    path.write_text(
        "experiment_name,scan_number,precursor_charge,preds,log_probs\n"
        "run,10,2,PEPTIDE,-0.1\nrun,11,2,PEPTIDEK,-0.2\nrun,12,2,SAMPLER,-0.3\n",
        encoding="utf-8",
    )
    return path


def _manifest_path(prefix: Path) -> Path:
    return prefix.with_name(prefix.name + "_prediction_manifest.json")


def test_status_can_be_inspected_while_the_batch_runs(candidates, tmp_path, monkeypatch):
    first_done, release = threading.Event(), threading.Event()
    calls = []

    def slow_predictor(config, peptidoform, charge):
        calls.append(peptidoform)
        if len(calls) == 2:
            first_done.set()
            release.wait(10)
        return mock_predicted_spectrum(peptidoform, charge)

    monkeypatch.setattr(workflows, "_predict_fragment_spectrum", slow_predictor)
    prefix = tmp_path / "batch"
    worker = threading.Thread(target=build_prediction_batch_tool,
                              kwargs={"path": str(candidates), "output_prefix": str(prefix), "predictor_backend": "mock"})
    worker.start()
    try:
        assert first_done.wait(10)
        running = inspect_prediction_batch_status_tool(str(_manifest_path(prefix)))
        assert json.loads(_manifest_path(prefix).read_text())["phase"] == "running"
    finally:
        release.set()
        worker.join(10)

    assert running["phase"] == "predicting" and running["completed"] == 1 and running["total"] == 3
    assert running["seconds_since_update"] >= 0
    finished = inspect_prediction_batch_status_tool(str(_manifest_path(prefix)))
    assert finished["phase"] == "completed" and finished["succeeded"] == 3
    assert "seconds_since_update" not in finished
    assert json.loads(_manifest_path(prefix).read_text())["phase"] == "completed"


def test_interrupted_batch_reports_its_last_status_and_cannot_be_scored(candidates, tmp_path, monkeypatch):
    def interrupted(config, peptidoform, charge):
        if peptidoform == "PEPTIDEK":
            raise KeyboardInterrupt  # not an Exception: stops the batch like a killed process would
        return mock_predicted_spectrum(peptidoform, charge)

    monkeypatch.setattr(workflows, "_predict_fragment_spectrum", interrupted)
    prefix = tmp_path / "batch"
    with pytest.raises(KeyboardInterrupt):
        build_prediction_batch_tool(path=str(candidates), output_prefix=str(prefix), predictor_backend="mock")
    time.sleep(0.05)

    status = inspect_prediction_batch_status_tool(str(_manifest_path(prefix)))
    assert status["phase"] == "predicting" and status["completed"] == 1 and status["seconds_since_update"] > 0
    with pytest.raises(ValueError, match="has not finished"):
        score_from_prediction_batch_tool(manifest_path=str(_manifest_path(prefix)), spectra_path=str(tmp_path / "x.feather"),
                                         output_prefix=str(tmp_path / "scored"))


def test_finished_manifest_keeps_its_fields_and_adds_the_phase(candidates, tmp_path):
    manifest = build_prediction_batch_tool(path=str(candidates), output_prefix=str(tmp_path / "batch"), predictor_backend="mock")
    assert manifest["phase"] == "completed"
    for key in ("batch_id", "input_path", "adapter", "predictor_backend", "unique_prediction_requests", "predictions_succeeded",
                "predictions_failed", "batch_seconds", "output_files"):
        assert key in manifest
    assert set(manifest["output_files"]) == {"predictions_csv", "fragments_csv", "prediction_map_csv", "failures_csv", "log_jsonl",
                                             "status_json", "manifest_json"}
