from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

import proteomics_mcp.tools.workflows as workflows
from proteomics_mcp.tools.workflows import score_candidates_tool, score_spectrum_tool, summarize_score_metrics_tool


def test_score_candidates_tool_scores_table(monkeypatch, tmp_path):
    candidates_path = tmp_path / "preds.csv"
    candidates_path.write_text(
        "experiment_name,scan_number,preds,log_probs\nrun,10,PEPTIDE,-0.1\n",
        encoding="utf-8",
    )
    spectra_path = tmp_path / "spectra.feather"
    frame = pd.DataFrame(
        [
            {
                "scan_number": 10,
                "spectrum_id": "run:10",
                "mz_array": [200.0, 300.0],
                "intensity_array": [10.0, 20.0],
                "precursor_mz": 500.0,
                "precursor_charge": 2,
            }
        ]
    )
    monkeypatch.setattr(workflows, "load_spectrum_from_table", lambda *_, **__: workflows.ObservedSpectrum(mz=[200.0, 300.0], intensity=[10.0, 20.0], precursor_mz=500.0, precursor_charge=2))
    output_path = tmp_path / "scores.csv"
    status_path = tmp_path / "scores_status.json"

    summary = score_candidates_tool(
        candidates_path=str(candidates_path),
        spectra_path=str(spectra_path),
        output_path=str(output_path),
        predictor_backend="mock",
        status_path=str(status_path),
    )

    assert summary["rows_scored"] == 1
    assert summary["rows_failed"] == 0
    assert summary["rows_preview"][0]["spectral_angle"] is not None
    assert summary["rows_preview"][0]["experiment_name"] == "run"
    assert "match_quality" in summary["metric_groups"]
    assert "ion_evidence" in summary["metric_groups"]
    assert "mass_error" in summary["metric_groups"]
    assert "acquisition" in summary["metric_groups"]
    assert "isolation_target_mz" in summary["rows_preview"][0]
    assert "precursor_charge" in summary["rows_preview"][0]
    assert output_path.exists()
    assert summary["status_path"] == str(status_path)
    status = json.loads(Path(status_path).read_text(encoding="utf-8"))
    assert status["phase"] == "completed"
    assert status["completed"] == 1
    assert status["succeeded"] == 1


def test_score_candidates_tool_supports_fragpipe_adapter(monkeypatch, tmp_path):
    candidates_path = tmp_path / "psm.tsv"
    candidates_path.write_text(
        "Spectrum\tPeptide\tCharge\tObserved M/Z\tHyperscore\tNextscore\tProbability\n"
        "runA.00582.00582.2\tKATEVSK\t2\t381.7211\t11.278\t9.500\t0.9840\n",
        encoding="utf-8",
    )
    seen: dict[str, object] = {}

    def fake_loader(path, *, scan_number=None, spectrum_id=None, experiment_name=None):
        seen["path"] = path
        seen["scan_number"] = scan_number
        seen["spectrum_id"] = spectrum_id
        seen["experiment_name"] = experiment_name
        return workflows.ObservedSpectrum(mz=[200.0, 300.0], intensity=[10.0, 20.0], precursor_mz=500.0, precursor_charge=2)

    monkeypatch.setattr(workflows, "load_spectrum_from_table", fake_loader)
    output_path = tmp_path / "scores.csv"

    summary = score_candidates_tool(
        candidates_path=str(candidates_path),
        spectra_path=str(tmp_path / "spectra.feather"),
        output_path=str(output_path),
        adapter="fragpipe",
        predictor_backend="mock",
    )

    assert summary["rows_scored"] == 1
    assert seen["scan_number"] == 582
    assert seen["experiment_name"] == "runA"
    assert seen["spectrum_id"] == "runA.00582.00582.2"


def test_score_candidates_tool_supports_proteome_discoverer_adapter(monkeypatch, tmp_path):
    candidates_path = tmp_path / "pd.xlsx"
    pd.DataFrame(
        [
            {
                "Annotated Sequence": "[R].KTQTVcNFTDGALVQHQEWDGKESTITR.[K]",
                "Spectrum File": "20231127_sample_run.raw",
                "First Scan": 11770,
                "Charge": 4,
                "m/z [Da]": 813.15192,
                "XCorr": 9.55,
                "DeltaScore": 0.7717,
                "Percolator PEP": 3.739192e-28,
            }
        ]
    ).to_excel(candidates_path, index=False)
    seen: dict[str, object] = {}

    def fake_loader(path, *, scan_number=None, spectrum_id=None, experiment_name=None):
        seen["path"] = path
        seen["scan_number"] = scan_number
        seen["spectrum_id"] = spectrum_id
        seen["experiment_name"] = experiment_name
        return workflows.ObservedSpectrum(mz=[200.0, 300.0], intensity=[10.0, 20.0], precursor_mz=500.0, precursor_charge=2)

    monkeypatch.setattr(workflows, "load_spectrum_from_table", fake_loader)
    output_path = tmp_path / "scores.csv"

    summary = score_candidates_tool(
        candidates_path=str(candidates_path),
        spectra_path=str(tmp_path / "spectra.feather"),
        output_path=str(output_path),
        adapter="proteome_discoverer",
        predictor_backend="mock",
    )

    assert summary["rows_scored"] == 1
    assert seen["scan_number"] == 11770
    assert seen["experiment_name"] == "20231127_sample_run"
    assert seen["spectrum_id"] == "20231127_sample_run:11770"


def test_score_candidates_tool_supports_maxquant_adapter(monkeypatch, tmp_path):
    candidates_path = tmp_path / "msms.txt"
    candidates_path.write_text(
        "Raw file\tScan number\tSequence\tModified sequence\tCharge\tm/z\tScore\tDelta score\tPEP\n"
        "runA\t14345\tAAANGLMQTSK\t_AAANGLMQTSK_\t2\t546.27932\t131.120\t129.790\t2.2652e-13\n",
        encoding="utf-8",
    )
    seen: dict[str, object] = {}

    def fake_loader(path, *, scan_number=None, spectrum_id=None, experiment_name=None):
        seen["path"] = path
        seen["scan_number"] = scan_number
        seen["spectrum_id"] = spectrum_id
        seen["experiment_name"] = experiment_name
        return workflows.ObservedSpectrum(mz=[200.0, 300.0], intensity=[10.0, 20.0], precursor_mz=500.0, precursor_charge=2)

    monkeypatch.setattr(workflows, "load_spectrum_from_table", fake_loader)
    output_path = tmp_path / "scores.csv"

    summary = score_candidates_tool(
        candidates_path=str(candidates_path),
        spectra_path=str(tmp_path / "spectra.feather"),
        output_path=str(output_path),
        adapter="maxquant",
        predictor_backend="mock",
    )

    assert summary["rows_scored"] == 1
    assert seen["scan_number"] == 14345
    assert seen["experiment_name"] == "runA"
    assert seen["spectrum_id"] == "runA:14345"


def test_score_candidates_tool_supports_casanovo_adapter(monkeypatch, tmp_path):
    candidates_path = tmp_path / "casanovo.csv"
    candidates_path.write_text(
        "prediction_id,predictions,log_probs,precursor_mz,precursor_charge,group,scan_number,spectrum_id\n"
        "0,DM[Oxidation]HGVTSHLSNQELQDLVEFMK,-0.24043669456534025,858.73639,3,runA,321,runA:321\n",
        encoding="utf-8",
    )
    seen: dict[str, object] = {}

    def fake_loader(path, *, scan_number=None, spectrum_id=None, experiment_name=None):
        seen["path"] = path
        seen["scan_number"] = scan_number
        seen["spectrum_id"] = spectrum_id
        seen["experiment_name"] = experiment_name
        return workflows.ObservedSpectrum(mz=[200.0, 300.0], intensity=[10.0, 20.0], precursor_mz=500.0, precursor_charge=2)

    monkeypatch.setattr(workflows, "load_spectrum_from_table", fake_loader)
    output_path = tmp_path / "scores.csv"

    summary = score_candidates_tool(
        candidates_path=str(candidates_path),
        spectra_path=str(tmp_path / "spectra.feather"),
        output_path=str(output_path),
        adapter="casanovo",
        predictor_backend="mock",
    )

    assert summary["rows_scored"] == 1
    assert seen["scan_number"] == 321
    assert seen["experiment_name"] == "runA"
    assert seen["spectrum_id"] == "runA:321"


def test_score_spectrum_tool_returns_agent_metric_summary():
    result = score_spectrum_tool(
        candidate={"candidate_id": "c1", "source_engine": "test", "peptidoform": "PEPTIDE", "precursor_charge": 2, "confidence": 0.9},
        observed={
            "mz": [100.0],
            "intensity": [10.0],
            "precursor_charge": 2,
            "metadata": {"attributes": [{"accession": "MS:1000422", "name": "higher energy beam-type collision-induced dissociation"}]},
        },
        predicted={
            "peptidoform": "PEPTIDE",
            "precursor_charge": 2,
            "predictor_name": "mock",
            "model": "HCD",
            "fragments": [{"mz": 100.0, "intensity": 1.0, "annotation": "b1+", "ion_type": "b", "ordinal": 1}],
        },
        annotation_ion_types="abyczximpr",
    )

    summary = result["metric_summary"]

    assert summary == summarize_score_metrics_tool(result)
    assert summary["match_quality"]["spectral_angle"] is not None
    assert summary["ion_evidence"]["matched_ion_count"] == 1
    assert "precursor_mass_error_ppm" in summary["mass_error"]
    assert summary["candidate"]["confidence_percent"] == 90.0
    assert summary["acquisition"]["fragmentation_method"] == "HCD"
    assert "instrument_model" in summary["acquisition"]
    assert result["annotation_settings"]["annotation_ion_types"] == "abyczximpr"
