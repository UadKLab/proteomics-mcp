from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from proteomics_mcp.core.dataset import (
    build_dashboard_summary,
    bin_numeric_distribution,
    bin_2d_numeric_distribution,
    candidate_models_to_frame,
    distribution_field_catalog,
    get_peptide_summary_row,
    get_peptidoform_summary_row,
    inspection_peptide_table,
    inspection_psm_table,
    render_dashboard_bundle_html,
    render_peptide_table_html,
    render_prediction_accuracy_html,
    render_psm_table_html,
    render_summary_html,
    render_table_browse_html,
    render_dataset_report_html,
    summarize_candidate_inspection_frame,
    summarize_categorical_distribution,
    summarize_cv_by_category,
    summarize_file_coverage,
    summarize_dataset_path,
    summarize_missing_values,
    summarize_modifications,
    summarize_numeric_distribution,
    summarize_overlap,
    summarize_prediction_accuracy,
    summarize_quant_quality,
    write_dataset_outputs,
)
from proteomics_mcp.tools.workflows import (
    bin_numeric_field_tool,
    bin_2d_numeric_fields_tool,
    build_id_dashboard_tool,
    build_prediction_batch_tool,
    inspect_prediction_batch_status_tool,
    inspect_prediction_batch_tool,
    score_from_prediction_batch_tool,
    inspect_scored_batch_status_tool,
    inspect_scored_batch_tool,
    inspect_peptide_table_tool,
    inspect_psm_table_tool,
    list_distribution_fields_tool,
    list_dataset_groups_tool,
    render_peptide_table_tool,
    render_psm_table_tool,
    render_dashboard_html_tool,
    render_summary_html_tool,
    render_dataset_report_tool,
    summarize_categorical_field_tool,
    summarize_cv_by_category_tool,
    summarize_dashboard_tool,
    summarize_file_coverage_tool,
    summarize_missing_values_tool,
    summarize_modifications_tool,
    summarize_numeric_field_tool,
    summarize_overlap_tool,
    summarize_prediction_accuracy_tool,
    render_prediction_accuracy_html_tool,
    summarize_quant_quality_tool,
    summarize_dataset_tool,
)
from proteomics_mcp.models.schemas import CandidatePSM
from proteomics_mcp.predictors.mock import mock_predicted_spectrum


def _scored_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "candidate_id": "c1",
                "source_engine": "instanovo",
                "experiment_name": "runA",
                "scan_number": 100,
                "spectrum_id": "runA:100",
                "peptidoform": "PEPTIDE",
                "precursor_charge": 2,
                "confidence": 0.95,
                "raw_score": -0.5,
                "next_beam": -1.0,
                "next_beam_margin": 0.5,
                "spectral_angle": 0.80,
                "cosine_similarity": 0.81,
                "pearson_correlation": 0.78,
                "xcorr": 4.0,
                "hyperscore": 5.0,
                "matched_ion_count": 12,
                "matched_ion_fraction": 0.60,
                "explained_intensity_fraction": 0.55,
                "annotated_peak_count": 22,
                "annotated_intensity_fraction": 0.66,
                "fragmentation_method": "HCD",
                "instrument_model": "Orbitrap Exploris",
                "collision_energy": 28.0,
                "ion_type_counts": {"b": 10, "y": 12},
            },
            {
                "candidate_id": "c2",
                "source_engine": "instanovo",
                "experiment_name": "runA",
                "scan_number": 101,
                "spectrum_id": "runA:101",
                "peptidoform": "PEP[UNIMOD:21]TIDE",
                "precursor_charge": 3,
                "confidence": 0.90,
                "raw_score": -0.4,
                "next_beam": -0.9,
                "next_beam_margin": 0.5,
                "spectral_angle": 0.70,
                "cosine_similarity": 0.71,
                "pearson_correlation": 0.69,
                "xcorr": 3.6,
                "hyperscore": 4.7,
                "matched_ion_count": 11,
                "matched_ion_fraction": 0.55,
                "explained_intensity_fraction": 0.50,
                "annotated_peak_count": 20,
                "annotated_intensity_fraction": 0.61,
                "fragmentation_method": "HCD",
                "instrument_model": "Orbitrap Exploris",
                "collision_energy": 29.0,
                "ion_type_counts": {"b": 8, "y": 10},
            },
            {
                "candidate_id": "c3",
                "source_engine": "fragpipe",
                "experiment_name": "runB",
                "scan_number": 200,
                "spectrum_id": "runB:200",
                "peptidoform": "PEPTIDE",
                "precursor_charge": 2,
                "confidence": None,
                "raw_score": 42.0,
                "next_beam": 38.0,
                "next_beam_margin": 4.0,
                "spectral_angle": 0.77,
                "cosine_similarity": 0.79,
                "pearson_correlation": 0.74,
                "xcorr": 5.2,
                "hyperscore": 6.5,
                "matched_ion_count": 13,
                "matched_ion_fraction": 0.64,
                "explained_intensity_fraction": 0.58,
                "annotated_peak_count": 24,
                "annotated_intensity_fraction": 0.68,
                "fragmentation_method": "HCD",
                "instrument_model": "Orbitrap Exploris",
                "collision_energy": 30.0,
                "ion_type_counts": {"b": 11, "y": 13},
            },
            {
                "candidate_id": "c4",
                "source_engine": "proteome_discoverer",
                "experiment_name": "runC",
                "scan_number": 300,
                "spectrum_id": "runC:300",
                "peptidoform": "M[UNIMOD:35]PEPTIDE",
                "precursor_charge": 4,
                "confidence": 0.99,
                "raw_score": 10.0,
                "next_beam": 7.0,
                "next_beam_margin": 3.0,
                "spectral_angle": 0.92,
                "cosine_similarity": 0.90,
                "pearson_correlation": 0.89,
                "xcorr": 6.0,
                "hyperscore": 7.0,
                "matched_ion_count": 18,
                "matched_ion_fraction": 0.75,
                "explained_intensity_fraction": 0.71,
                "annotated_peak_count": 30,
                "annotated_intensity_fraction": 0.80,
                "fragmentation_method": "HCD",
                "instrument_model": "Orbitrap Exploris",
                "collision_energy": 27.0,
                "ion_type_counts": {"b": 14, "y": 16},
            },
        ]
    )


def test_dataset_summary_collapses_peptides_and_preserves_peptidoforms(tmp_path):
    scored_path = tmp_path / "scores.csv"
    _scored_frame().to_csv(scored_path, index=False)

    aggregation = summarize_dataset_path(scored_path)

    assert aggregation.dataset_summary["dataset_metadata"]["total_psm_count"] == 4
    assert aggregation.dataset_summary["dataset_metadata"]["distinct_peptide_count"] == 2
    assert aggregation.dataset_summary["dataset_metadata"]["distinct_peptidoform_count"] == 3

    peptide = get_peptide_summary_row(aggregation, "PEPTIDE")
    assert peptide is not None
    assert peptide["total_psm_count"] == 3
    assert peptide["distinct_peptidoform_count"] == 2
    assert peptide["modified_psm_count"] == 1
    assert aggregation.dataset_summary["dataset_metadata"]["distinct_peptide_precursor_count"] == 3
    assert aggregation.dataset_summary["dataset_metadata"]["distinct_peptidoform_precursor_count"] == 3
    assert peptide["best_candidate_id"] == "c1"
    assert peptide["per_file_psm_counts"] == {"runA": 2, "runB": 1}
    assert peptide["per_charge_counts"] == {"2": 2, "3": 1}

    peptidoform = get_peptidoform_summary_row(aggregation, "PEP[UNIMOD:21]TIDE")
    assert peptidoform is not None
    assert peptidoform["total_psm_count"] == 1
    assert peptidoform["best_candidate_id"] == "c2"


def test_best_psm_selection_prefers_confidence_then_raw_score(tmp_path):
    frame = pd.DataFrame(
        [
            {"candidate_id": "high_raw", "peptidoform": "PEPTIDE", "experiment_name": "runA", "spectrum_id": "runA:1", "confidence": None, "raw_score": 100.0, "spectral_angle": 0.99, "matched_ion_count": 20},
            {"candidate_id": "has_conf", "peptidoform": "PEPTIDE", "experiment_name": "runA", "spectrum_id": "runA:2", "confidence": 0.80, "raw_score": 1.0, "spectral_angle": 0.10, "matched_ion_count": 1},
            {"candidate_id": "conf_tie_better_raw", "peptidoform": "KPEPTIDE", "experiment_name": "runA", "spectrum_id": "runA:3", "confidence": 0.50, "raw_score": 2.0, "spectral_angle": 0.20, "matched_ion_count": 2},
            {"candidate_id": "conf_tie_lower_raw", "peptidoform": "KPEPTIDE", "experiment_name": "runA", "spectrum_id": "runA:4", "confidence": 0.50, "raw_score": 1.0, "spectral_angle": 0.90, "matched_ion_count": 9},
        ]
    )
    scored_path = tmp_path / "scores.csv"
    frame.to_csv(scored_path, index=False)

    aggregation = summarize_dataset_path(scored_path)

    peptide = get_peptide_summary_row(aggregation, "PEPTIDE")
    other = get_peptide_summary_row(aggregation, "KPEPTIDE")

    assert peptide is not None and peptide["best_candidate_id"] == "has_conf"
    assert other is not None and other["best_candidate_id"] == "conf_tie_better_raw"


def test_dataset_summary_supports_sdrf_and_writes_outputs(tmp_path):
    scored_path = tmp_path / "scores.csv"
    _scored_frame().to_csv(scored_path, index=False)
    sdrf_path = tmp_path / "meta.sdrf.tsv"
    sdrf_path.write_text(
        "assay name\tcomment[data file]\tsource name\tcharacteristics[organism]\n"
        "assayA\trunA.raw\tsampleA\tHomo sapiens\n"
        "assayB\trunB.raw\tsampleB\tHomo sapiens\n",
        encoding="utf-8",
    )

    summary = summarize_dataset_tool(scored_path=str(scored_path), output_prefix=str(tmp_path / "dataset"), sdrf_path=str(sdrf_path))

    assert summary["sdrf_summary"]["available"] is True
    assert "sample_summary_csv" in summary["output_files"]
    assert Path(summary["output_files"]["dataset_summary_json"]).exists()


def test_dataset_report_and_group_tools(tmp_path):
    scored_path = tmp_path / "scores.csv"
    _scored_frame().to_csv(scored_path, index=False)
    aggregation = summarize_dataset_path(scored_path)
    report_path = tmp_path / "report.html"

    html = render_dataset_report_html(aggregation, str(report_path))
    assert report_path.exists()
    assert "Dataset report" in html
    assert "Top confident peptides" in html
    assert "Confidence vs spectrum-evidence outliers" in html

    tool_result = render_dataset_report_tool(scored_path=str(scored_path), output_path=str(tmp_path / "tool_report.html"))
    groups = list_dataset_groups_tool(scored_path=str(scored_path))

    assert tool_result["format"] == "html"
    assert groups["dataset"]["total_psm_count"] == 4
    assert "runA" in groups["runs"]


def test_candidate_inspection_builds_psm_and_peptide_tables(tmp_path):
    candidate_path = tmp_path / "preds.csv"
    candidate_path.write_text(
        "experiment_name,scan_number,preds,log_probs,token_log_probs,precursor_charge,predictions_log_probability_beam_1\n"
        "runA,10,PEPTIDE,-0.1,\"[-0.01, -0.02]\",2,-0.2\n"
        "runB,11,PEP[UNIMOD:21]TIDE,-0.2,\"[-0.03, -0.04]\",3,-0.4\n",
        encoding="utf-8",
    )

    result = inspect_psm_table_tool(path=str(candidate_path), adapter="instanovo")

    assert result["dataset_metadata"]["total_psm_count"] == 2
    assert result["dataset_metadata"]["distinct_peptide_count"] == 1
    assert result["dataset_metadata"]["distinct_peptide_precursor_count"] == 2
    assert result["dataset_metadata"]["distinct_peptidoform_precursor_count"] == 2
    assert "output_files" not in result
    assert result["peptide_preview"][0]["distinct_peptidoform_count"] == 2
    assert result["peptide_preview"][0]["distinct_charge_count"] == 2
    assert "best_spectral_angle" not in result["peptide_preview"][0]
    assert "source_columns" not in result["psm_preview"][0]


def test_candidate_inspection_tables_are_id_only():
    candidates = [
        CandidatePSM(
            candidate_id="c1",
            source_engine="instanovo",
            peptidoform="PEPTIDE",
            experiment_name="runA",
            scan_number=1,
            confidence=0.9,
            precursor_charge=2,
            token_confidence_mean=0.95,
            token_confidence_min=0.91,
            token_confidence_min_position=3,
            token_confidence_min_token="P",
        ),
    ]
    inspection = summarize_candidate_inspection_frame(candidate_models_to_frame(candidates))
    psm = inspection_psm_table(inspection)
    peptide = inspection_peptide_table(inspection)

    assert "source_columns" not in psm.columns
    assert "token_log_probabilities" not in psm.columns
    assert "best_spectral_angle" not in peptide.columns
    assert "fragmentation_methods" not in peptide.columns
    assert "best_token_confidence_mean" in peptide.columns


def test_candidate_inspection_returns_single_peptide_summary(tmp_path):
    candidates = [
        CandidatePSM(candidate_id="c1", source_engine="instanovo", peptidoform="PEPTIDE", experiment_name="runA", scan_number=1, confidence=0.9, precursor_charge=2),
        CandidatePSM(candidate_id="c2", source_engine="instanovo", peptidoform="PEP[UNIMOD:21]TIDE", experiment_name="runB", scan_number=2, confidence=0.8, precursor_charge=3),
    ]
    inspection = summarize_candidate_inspection_frame(candidate_models_to_frame(candidates))
    assert inspection.dataset_summary["dataset_metadata"]["distinct_peptide_count"] == 1
    candidate_path = tmp_path / "preds.csv"
    candidate_path.write_text(
        "experiment_name,scan_number,preds,log_probs,precursor_charge,predictions_log_probability_beam_1\n"
        "runA,1,PEPTIDE,-0.1,2,-0.2\n"
        "runB,2,PEP[UNIMOD:21]TIDE,-0.2,3,-0.4\n",
        encoding="utf-8",
    )

    summary = inspect_peptide_table_tool(path=str(candidate_path), adapter="instanovo", peptide="PEPTIDE")
    assert summary["distinct_peptidoform_count"] == 2
    assert summary["best_peptidoform"] == "PEPTIDE"


def test_render_table_tools_write_html_only(tmp_path):
    candidate_path = tmp_path / "preds.csv"
    candidate_path.write_text(
        "experiment_name,scan_number,preds,log_probs,precursor_charge,predictions_log_probability_beam_1\n"
        "runA,1,PEPTIDE,-0.1,2,-0.2\n"
        "runB,2,PEP[UNIMOD:21]TIDE,-0.2,3,-0.4\n",
        encoding="utf-8",
    )

    psm_result = render_psm_table_tool(path=str(candidate_path), adapter="instanovo", output_path=str(tmp_path / "psm.html"))
    peptide_result = render_peptide_table_tool(path=str(candidate_path), adapter="instanovo", output_path=str(tmp_path / "peptide.html"))

    assert Path(psm_result["output_path"]).exists()
    assert Path(peptide_result["output_path"]).exists()
    assert psm_result["format"] == "html"
    assert peptide_result["format"] == "html"


def test_render_table_helpers_write_html_only(tmp_path):
    candidates = [
        CandidatePSM(candidate_id="c1", source_engine="instanovo", peptidoform="PEPTIDE", experiment_name="runA", scan_number=1, confidence=0.9, precursor_charge=2),
        CandidatePSM(candidate_id="c2", source_engine="instanovo", peptidoform="PEP[UNIMOD:21]TIDE", experiment_name="runB", scan_number=2, confidence=0.8, precursor_charge=3),
    ]
    inspection = summarize_candidate_inspection_frame(candidate_models_to_frame(candidates))

    psm_outputs = render_psm_table_html(inspection, tmp_path / "psm_only.html")
    peptide_outputs = render_peptide_table_html(inspection, tmp_path / "peptide_only.html")

    assert Path(psm_outputs["output_path"]).exists()
    assert Path(peptide_outputs["output_path"]).exists()


def test_render_table_browse_html_writes_simple_html(tmp_path):
    frame = pd.DataFrame(
        [
            {"sequence": "PEPTIDE", "peptidoforms": ["PEPTIDE", "PEP[UNIMOD:21]TIDE"], "confidence": 0.99},
            {"sequence": "KPEPTIDE", "peptidoforms": '["KPEPTIDE"]', "confidence": 0.95},
        ]
    )
    output = tmp_path / "table.html"

    html = render_table_browse_html(frame, title="PSM table", output_path=output)

    assert output.exists()
    assert "PSM table" in html
    assert "2 rows, 3 columns" in html
    assert "PEPTIDE" in html
    assert "Download CSV" in html
    assert "Columns" in html
    assert "Clear pins" in html
    assert "PEPTIDE, PEP[UNIMOD:21]TIDE" in html
    assert '[&quot;PEPTIDE&quot;' not in html


def test_render_summary_html_writes_simple_summary(tmp_path):
    output = tmp_path / "summary.html"
    summary = {"overview": {"total_psm_count": 10}, "rows": [{"run": "runA", "count": 4}, {"run": "runB", "count": 6}]}
    html = render_summary_html(summary, title="Dashboard summary", output_path=output)
    assert output.exists()
    assert "Dashboard summary" in html
    assert "total_psm_count" in html
    assert "runA" in html

    tool_result = render_summary_html_tool(summary, str(tmp_path / "tool_summary.html"), title="Tool summary")
    assert tool_result["format"] == "html"


def test_render_dashboard_bundle_html_writes_dashboard(tmp_path):
    output = tmp_path / "dashboard.html"
    summary = {
        "dashboard": {
            "overview": {
                "total_psm_count": 10,
                "distinct_peptide_count": 5,
                "distinct_peptidoform_count": 6,
                "distinct_peptidoform_precursor_count": 7,
                "distinct_run_count": 2,
            },
            "ids_per_file": {"counts": {"runA": 4, "runB": 6}},
            "missing_values": {"fields": [{"field": "confidence", "missing_fraction": 0.1}]},
            "histograms": {
                "precursor_mz": {"bins": [{"left": 400.0, "right": 500.0, "count": 3}]},
                "peptide_length": {"bins": [{"left": 8.0, "right": 10.0, "count": 2}]},
                "confidence": {"bins": [{"left": 0.0, "right": 0.5, "count": 2}]},
            },
            "scatter_2d": {
                "confidence_vs_spectral_angle": {
                    "points": [
                        {"x": 0.99, "y": 0.82, "group": "runA", "candidate_id": "c1", "peptidoform": "PEPTIDE", "run_key": "runA", "scan_number": 10},
                        {"x": 0.95, "y": 0.76, "group": "runB", "candidate_id": "c2", "peptidoform": "TESTIDE", "run_key": "runB", "scan_number": 20},
                    ]
                }
            },
        },
        "overlap_peptide": {"pairwise": [{"group_a": "runA", "group_b": "runB", "jaccard": 0.5}]},
        "modifications": {"modified_fraction": 0.2, "distinct_modified_peptidoform_count": 2, "modification_counts": {"UNIMOD:35": 3}, "group_modification_counts": {"runA": {"UNIMOD:35": 1}, "runB": {"UNIMOD:35": 2}}},
        "file_coverage": {
            "files": {
                "runA": {"overview": {"total_psm_count": 4, "distinct_peptide_count": 3, "distinct_peptidoform_count": 3, "modified_psm_count": 1}, "charge_counts": {"2": 2, "3": 2}, "confidence_histogram": {"bins": [{"left": 0.0, "right": 0.5, "count": 1}]}, "rt_histogram": {"bins": [{"left": 10.0, "right": 12.0, "count": 2}]}, "peptide_length_histogram": {"bins": [{"left": 8.0, "right": 10.0, "count": 2}]}},
                "runB": {"overview": {"total_psm_count": 6, "distinct_peptide_count": 4, "distinct_peptidoform_count": 5, "modified_psm_count": 1}, "charge_counts": {"2": 3, "3": 3}, "confidence_histogram": {"bins": [{"left": 0.5, "right": 1.0, "count": 2}]}, "rt_histogram": {"bins": [{"left": 20.0, "right": 22.0, "count": 3}]}, "peptide_length_histogram": {"bins": [{"left": 10.0, "right": 12.0, "count": 3}]}},
            }
        },
    }
    html = render_dashboard_bundle_html(summary, title="Dashboard", output_path=output)
    assert output.exists()
    assert "Peptide overlap (Jaccard)" in html
    assert html.index("Peptide overlap (Jaccard)") < html.index("Per-file coverage")
    assert "IDs per file" not in html
    assert "Modification counts" in html
    assert "Charge distribution by run" in html
    assert "Conf. vs SA" in html
    tool_result = render_dashboard_html_tool(summary, str(tmp_path / "tool_dashboard.html"), title="Tool dashboard")
    assert tool_result["format"] == "html"


def test_build_id_dashboard_tool_enriches_rt_from_feather(tmp_path):
    candidate_path = tmp_path / "preds.csv"
    candidate_path.write_text(
        "experiment_name,scan_number,preds,log_probs,precursor_charge,precursor_mz,predictions_log_probability_beam_1\n"
        "runA,10,PEPTIDE,-0.1,2,500.0,-0.2\n"
        "runB,11,TESTIDE,-0.2,3,650.0,-0.3\n",
        encoding="utf-8",
    )
    spectra_path = tmp_path / "spectra.feather"
    pd.DataFrame(
        [
            {"source_file": "data/mzml/runA.mzML", "scan_number": 10, "rt": 600.0, "filter_string": "f1", "scan_resolution": 60000.0},
            {"source_file": "data/mzml/runB.mzML", "scan_number": 11, "rt": 1200.0, "filter_string": "f2", "scan_resolution": 120000.0},
        ]
    ).to_feather(spectra_path)

    result = build_id_dashboard_tool(
        path=str(candidate_path),
        spectra_path=str(spectra_path),
        output_prefix=str(tmp_path / "dashboard"),
        adapter="instanovo",
        top_n=2,
    )

    assert result["metadata_join"]["matched_rows"] == 2
    assert Path(result["output_files"]["dashboard_html"]).exists()
    assert Path(result["output_files"]["summary_json"]).exists()
    enriched = pd.read_csv(result["output_files"]["enriched_csv"])
    assert enriched["retention_time_minutes"].tolist() == [10.0, 20.0]


def test_build_id_dashboard_fills_missing_instanovo_metadata_from_feather(tmp_path):
    candidate_path = tmp_path / "preds.csv"
    candidate_path.write_text("spectrum_index,preds,log_probs\n10,PEPTIDE,-0.1\n", encoding="utf-8")
    spectra_path = tmp_path / "spectra.feather"
    pd.DataFrame(
        [{"source_file": "data/mzml/runA.mzML", "scan_number": 10, "rt": 600.0, "precursor_mz": 500.2, "precursor_charge": 2}]
    ).to_feather(spectra_path)

    result = build_id_dashboard_tool(
        path=str(candidate_path),
        spectra_path=str(spectra_path),
        output_prefix=str(tmp_path / "dashboard"),
        adapter="instanovo",
        default_experiment_name="runA",
    )

    enriched = pd.read_csv(result["output_files"]["enriched_csv"])
    assert result["metadata_join"]["matched_rows"] == 1
    assert enriched["precursor_mz"].tolist() == [500.2]
    assert enriched["precursor_charge"].tolist() == [2]


def test_build_id_dashboard_tool_works_without_spectra_path(tmp_path):
    candidate_path = tmp_path / "preds.csv"
    candidate_path.write_text(
        "prediction_id,precursor_mz,precursor_charge,prediction,log_probability\n"
        "0,500.0,2,PEPTIDE,-0.1\n"
        "1,650.0,3,TESTIDE,-0.2\n",
        encoding="utf-8",
    )

    result = build_id_dashboard_tool(
        path=str(candidate_path),
        output_prefix=str(tmp_path / "dashboard_id_only"),
        adapter="instanovo",
    )

    assert result["metadata_join"] is None
    dashboard_html = Path(result["output_files"]["dashboard_html"])
    assert dashboard_html.exists()
    html = dashboard_html.read_text(encoding="utf-8")
    assert "<div class='value'>2</div>" in html
    assert "NA" not in html.split("<div class='tile'>")[1]
    summary = json.loads(Path(result["output_files"]["summary_json"]).read_text(encoding="utf-8"))
    assert summary["dashboard"]["overview"]["total_psm_count"] == 2


def test_render_dashboard_html_tool_accepts_flat_summarize_dashboard_output(tmp_path):
    flat_summary = {
        "overview": {
            "total_psm_count": 12,
            "distinct_peptide_count": 7,
            "distinct_peptidoform_count": 8,
            "distinct_peptidoform_precursor_count": 9,
            "distinct_run_count": 0,
        },
        "ids_per_file": {"counts": {}},
        "missing_values": {"fields": []},
        "histograms": {
            "precursor_mz": {"bins": [{"left": 400.0, "right": 500.0, "count": 6}]},
            "confidence": {"bins": [{"left": 0.0, "right": 0.5, "count": 6}]},
        },
        "scatter_2d": {},
    }
    output = tmp_path / "flat_dashboard.html"
    result = render_dashboard_html_tool(flat_summary, str(output), title="Flat dashboard")
    assert result["format"] == "html"
    html = output.read_text(encoding="utf-8")
    assert "<div class='value'>12</div>" in html
    assert "<div class='value'>7</div>" in html
    assert "Plotly.newPlot" in html


def test_build_and_inspect_prediction_batch_tool(tmp_path):
    candidate_path = tmp_path / "preds.csv"
    candidate_path.write_text(
        "experiment_name,scan_number,preds,log_probs,precursor_charge,predictions_log_probability_beam_1\n"
        "runA,1,PEPTIDE,-0.1,2,-0.2\n"
        "runA,2,PEPTIDE,-0.2,2,-0.3\n"
        "runB,3,TESTIDE,-0.3,3,-0.4\n",
        encoding="utf-8",
    )

    result = build_prediction_batch_tool(
        path=str(candidate_path),
        output_prefix=str(tmp_path / "batch"),
        adapter="instanovo",
        predictor_backend="mock",
    )

    assert result["unique_prediction_requests"] == 2
    assert result["predictions_succeeded"] == 2
    assert Path(result["output_files"]["predictions_csv"]).exists()
    assert Path(result["output_files"]["fragments_csv"]).exists()
    assert Path(result["output_files"]["prediction_map_csv"]).exists()
    assert Path(result["output_files"]["log_jsonl"]).exists()
    assert Path(result["output_files"]["status_json"]).exists()

    predictions = pd.read_csv(result["output_files"]["predictions_csv"])
    assert len(predictions) == 2
    assert set(predictions["peptidoform"]) == {"PEPTIDE", "TESTIDE"}

    prediction_map = pd.read_csv(result["output_files"]["prediction_map_csv"])
    assert len(prediction_map) == 3
    assert prediction_map["prediction_id"].nunique() == 2

    inspection = inspect_prediction_batch_tool(result["output_files"]["manifest_json"])
    assert inspection["overview"]["unique_prediction_requests"] == 2
    assert inspection["overview"]["mapped_candidate_rows"] == 3
    assert inspection["fragment_count_summary"]["count"] == 2

    status = inspect_prediction_batch_status_tool(result["output_files"]["manifest_json"])
    assert status["phase"] == "completed"
    assert status["completed"] == 2
    assert status["succeeded"] == 2
    assert status["failed"] == 0


def test_build_prediction_batch_uses_native_ms2pip_batch(monkeypatch, tmp_path):
    candidate_path = tmp_path / "preds.csv"
    candidate_path.write_text(
        "experiment_name,scan_number,preds,log_probs,precursor_charge\n"
        "runA,1,PEPTIDE,-0.1,2\n"
        "runA,2,TESTIDE,-0.2,3\n",
        encoding="utf-8",
    )

    def fake_predict_batch(requests, **kwargs):
        assert requests == [("PEPTIDE", 2), ("TESTIDE", 3)]
        return [mock_predicted_spectrum(peptidoform, charge) for peptidoform, charge in requests]

    monkeypatch.setattr("proteomics_mcp.tools.workflows.predict_fragments_ms2pip_batch", fake_predict_batch)
    result = build_prediction_batch_tool(
        path=str(candidate_path),
        output_prefix=str(tmp_path / "batch"),
        adapter="instanovo",
        predictor_backend="ms2pip",
    )

    assert result["predictions_succeeded"] == 2


def test_score_from_prediction_batch_tool(tmp_path):
    candidate_path = tmp_path / "preds.csv"
    candidate_path.write_text(
        "experiment_name,scan_number,preds,log_probs,precursor_charge,precursor_mz,predictions_log_probability_beam_1\n"
        "runA,10,PEPTIDE,-0.1,2,500.0,-0.2\n",
        encoding="utf-8",
    )
    batch = build_prediction_batch_tool(
        path=str(candidate_path),
        output_prefix=str(tmp_path / "batch"),
        adapter="instanovo",
        predictor_backend="mock",
    )
    spectra_path = tmp_path / "spectra.feather"
    pd.DataFrame(
        [
            {
                "source_file": "data/mzml/runA.mzML",
                "scan_number": 10,
                "spectrum_id": "runA:10",
                "precursor_mz": 500.0,
                "precursor_charge": 2,
                "rt": 600.0,
                "mz_array": [100.0, 200.0, 300.2, 400.2, 500.0, 600.2],
                "intensity_array": [10.0, 30.0, 40.0, 25.0, 12.0, 18.0],
                "metadata": {},
            }
        ]
    ).to_feather(spectra_path)

    scored = score_from_prediction_batch_tool(
        manifest_path=batch["output_files"]["manifest_json"],
        spectra_path=str(spectra_path),
        output_prefix=str(tmp_path / "scored"),
    )
    assert scored["rows_scored"] == 1
    assert Path(scored["output_files"]["scored_csv"]).exists()
    assert Path(scored["output_files"]["status_json"]).exists()

    inspection = inspect_scored_batch_tool(scored["output_files"]["manifest_json"])
    assert inspection["overview"]["rows_scored"] == 1
    assert inspection["spectral_angle_summary"]["count"] == 1

    status = inspect_scored_batch_status_tool(scored["output_files"]["manifest_json"])
    assert status["phase"] == "completed"
    assert status["completed"] == 1
    assert status["succeeded"] == 1
    assert status["failed"] == 0


def test_score_from_prediction_batch_tool_ignores_manifest_top_n_for_candidate_reload(tmp_path):
    candidate_path = tmp_path / "preds.csv"
    rows = ["experiment_name,scan_number,preds,log_probs,precursor_charge,precursor_mz,predictions_log_probability_beam_1"]
    for idx in range(5):
        rows.append(f"runA,{idx+1},LOW{idx},-10.0,2,500.0,-11.0")
    rows.append("runA,10,PEPTIDE,-0.1,2,500.0,-0.2")
    candidate_path.write_text("\n".join(rows) + "\n", encoding="utf-8")

    batch = build_prediction_batch_tool(
        path=str(candidate_path),
        output_prefix=str(tmp_path / "batch"),
        adapter="instanovo",
        predictor_backend="mock",
        top_n=1,
    )
    spectra_path = tmp_path / "spectra.feather"
    pd.DataFrame(
        [
            {
                "source_file": "data/mzml/runA.mzML",
                "scan_number": 10,
                "spectrum_id": "runA:10",
                "precursor_mz": 500.0,
                "precursor_charge": 2,
                "rt": 600.0,
                "mz_array": [100.0, 200.0, 300.2, 400.2, 500.0, 600.2],
                "intensity_array": [10.0, 30.0, 40.0, 25.0, 12.0, 18.0],
                "metadata": {},
            }
        ]
    ).to_feather(spectra_path)

    scored = score_from_prediction_batch_tool(
        manifest_path=batch["output_files"]["manifest_json"],
        spectra_path=str(spectra_path),
        output_prefix=str(tmp_path / "scored"),
    )
    assert scored["rows_scored"] == 1
    assert scored["rows_failed"] == 0


def test_candidate_inspection_and_dashboard_accept_scored_enrichment(tmp_path):
    candidate_path = tmp_path / "preds.csv"
    candidate_path.write_text(
        "experiment_name,scan_number,preds,log_probs,precursor_charge,precursor_mz,predictions_log_probability_beam_1\n"
        "runA,10,PEPTIDE,-0.1,2,500.0,-0.2\n"
        "runB,11,TESTIDE,-0.2,3,650.0,-0.3\n",
        encoding="utf-8",
    )
    scored_path = tmp_path / "scored.csv"
    pd.DataFrame(
        [
                {
                    "candidate_id": "preds:0",
                    "peptidoform": "PEPTIDE",
                    "spectral_angle": 0.83,
                    "pearson_correlation": 0.79,
                    "xcorr": 4.2,
                    "matched_ion_fraction": 0.61,
                    "explained_intensity_fraction": 0.57,
                    "matched_ion_count": 12,
                    "precursor_mass_error_ppm": 3.5,
                    "fragment_mass_error_abs_mean_ppm": 4.2,
                    "predictor_name": "ms2pip",
                    "predictor_model": "HCD",
                    "prediction_id": "pred:1",
                },
                {
                    "candidate_id": "preds:1",
                    "peptidoform": "TESTIDE",
                    "spectral_angle": 0.72,
                    "pearson_correlation": 0.68,
                    "xcorr": 3.8,
                    "matched_ion_fraction": 0.52,
                    "explained_intensity_fraction": 0.48,
                    "matched_ion_count": 9,
                    "precursor_mass_error_ppm": -2.1,
                    "fragment_mass_error_abs_mean_ppm": 5.8,
                    "predictor_name": "ms2pip",
                    "predictor_model": "HCD",
                    "prediction_id": "pred:2",
                },
        ]
    ).to_csv(scored_path, index=False)
    spectra_path = tmp_path / "spectra.feather"
    pd.DataFrame(
        [
            {"source_file": "data/mzml/runA.mzML", "scan_number": 10, "rt": 600.0},
            {"source_file": "data/mzml/runB.mzML", "scan_number": 11, "rt": 1200.0},
        ]
    ).to_feather(spectra_path)

    inspection = inspect_psm_table_tool(path=str(candidate_path), adapter="instanovo", scored_path=str(scored_path))
    assert inspection["scoring_summary"]["scored_psm_count"] == 2
    assert "spectral_angle" in inspection["psm_preview"][0]
    assert "best_spectral_angle" in inspection["peptide_preview"][0]

    dashboard = summarize_dashboard_tool(
        path=str(candidate_path),
        adapter="instanovo",
        scored_path=str(scored_path),
    )
    assert "spectral_angle" in dashboard["histograms"]
    assert "confidence_vs_spectral_angle" in dashboard["grids_2d"]
    assert "confidence_vs_spectral_angle" in dashboard["scatter_2d"]
    assert "spectral_angle_vs_pearson" in dashboard["scatter_2d"]
    assert "spectral_angle_vs_xcorr" in dashboard["scatter_2d"]
    assert "confidence_vs_matched_ion_count" in dashboard["scatter_2d"]
    assert "spectral_angle_vs_precursor_mass_error_ppm" in dashboard["scatter_2d"]
    assert "spectral_angle_vs_fragment_mass_error_abs_mean_ppm" in dashboard["scatter_2d"]

    built = build_id_dashboard_tool(
        path=str(candidate_path),
        spectra_path=str(spectra_path),
        output_prefix=str(tmp_path / "dashboard"),
        adapter="instanovo",
        scored_path=str(scored_path),
        top_n=2,
    )
    summary = json.loads(Path(built["output_files"]["summary_json"]).read_text(encoding="utf-8"))
    assert "spectral_angle" in summary["dashboard"]["histograms"]


def test_distribution_field_catalog_and_numeric_summary_from_source_columns():
    candidates = [
        CandidatePSM(
            candidate_id="c1",
            source_engine="instanovo",
            peptidoform="PEPTIDE",
            experiment_name="runA",
            scan_number=1,
            confidence=0.9,
            precursor_mz=500.0,
            precursor_charge=2,
            source_columns={"p_value": 0.01, "abundance": 1000, "sample_type": "case"},
        ),
        CandidatePSM(
            candidate_id="c2",
            source_engine="instanovo",
            peptidoform="PEPTIDE",
            experiment_name="runB",
            scan_number=2,
            confidence=0.8,
            precursor_mz=700.0,
            precursor_charge=3,
            source_columns={"p_value": 0.02, "abundance": 1500, "sample_type": "control"},
        ),
    ]
    frame = candidate_models_to_frame(candidates)

    catalog = distribution_field_catalog(frame)
    assert "confidence" in catalog["numeric_fields"]
    assert "source_columns.p_value" in catalog["source_numeric_fields"]
    assert "source_columns.abundance" in catalog["source_numeric_fields"]
    assert "source_columns.sample_type" in catalog["source_categorical_fields"]
    assert "source_columns.p_value" in catalog["semantic_hints"]["p_value"]

    summary = summarize_numeric_distribution(frame, field="source_columns.p_value")
    assert summary["count"] == 2
    assert summary["min"] == 0.01
    assert summary["max"] == 0.02


def test_numeric_distribution_grouping_and_binning():
    candidates = [
        CandidatePSM(candidate_id="c1", source_engine="instanovo", peptidoform="PEPTIDE", experiment_name="runA", scan_number=1, confidence=0.9, precursor_mz=500.0, precursor_charge=2),
        CandidatePSM(candidate_id="c2", source_engine="instanovo", peptidoform="PEPTIDE", experiment_name="runA", scan_number=2, confidence=0.8, precursor_mz=550.0, precursor_charge=2),
        CandidatePSM(candidate_id="c3", source_engine="instanovo", peptidoform="TESTIDE", experiment_name="runB", scan_number=3, confidence=0.7, precursor_mz=650.0, precursor_charge=3),
    ]
    frame = candidate_models_to_frame(candidates)

    grouped = summarize_numeric_distribution(frame, field="confidence", group_by="experiment_name")
    assert grouped["groups"]["runA"]["count"] == 2
    assert grouped["groups"]["runB"]["mean"] == 0.7

    binned = bin_numeric_distribution(frame, field="precursor_mz", bin_edges=[500.0, 600.0, 700.0])
    assert binned["bins"][0]["count"] == 2
    assert binned["bins"][1]["count"] == 1

    weighted = summarize_numeric_distribution(frame, field="precursor_mz", weight_field="source_columns.abundance")
    assert weighted["weight_field"] == "source_columns.abundance"
    assert weighted["weighted_count"] == 0


def test_weighted_numeric_distribution_and_categorical_summary():
    candidates = [
        CandidatePSM(candidate_id="c1", source_engine="instanovo", peptidoform="PEPTIDE", experiment_name="runA", scan_number=1, precursor_mz=500.0, precursor_charge=2, source_columns={"abundance": 1.0, "sample_type": "case"}),
        CandidatePSM(candidate_id="c2", source_engine="instanovo", peptidoform="PEPTIDE", experiment_name="runA", scan_number=2, precursor_mz=700.0, precursor_charge=3, source_columns={"abundance": 3.0, "sample_type": "control"}),
        CandidatePSM(candidate_id="c3", source_engine="instanovo", peptidoform="TESTIDE", experiment_name="runB", scan_number=3, precursor_mz=900.0, precursor_charge=3, source_columns={"abundance": 2.0, "sample_type": "case"}),
    ]
    frame = candidate_models_to_frame(candidates)

    weighted = summarize_numeric_distribution(frame, field="precursor_mz", weight_field="source_columns.abundance")
    assert weighted["weighted_mean"] == 733.3333333333334
    assert weighted["total_weight"] == 6.0

    weighted_bins = bin_numeric_distribution(
        frame,
        field="precursor_mz",
        bin_edges=[400.0, 800.0, 1000.0],
        weight_field="source_columns.abundance",
    )
    assert weighted_bins["bins"][0]["count"] == 2
    assert weighted_bins["bins"][0]["weight_sum"] == 4.0
    assert weighted_bins["bins"][1]["weight_sum"] == 2.0

    categorical = summarize_categorical_distribution(frame, field="precursor_charge")
    assert categorical["counts"] == {"3": 2, "2": 1}

    grouped_categorical = summarize_categorical_distribution(frame, field="source_columns.sample_type", group_by="experiment_name")
    assert grouped_categorical["groups"]["runA"]["counts"] == {"case": 1, "control": 1}


def test_distribution_tools_on_candidate_table(tmp_path):
    candidate_path = tmp_path / "preds.csv"
    candidate_path.write_text(
        "experiment_name,scan_number,preds,log_probs,precursor_charge,precursor_mz,predictions_log_probability_beam_1,p_value,abundance\n"
        "runA,1,PEPTIDE,-0.1,2,500.0,-0.2,0.01,1000\n"
        "runB,2,TESTIDE,-0.2,3,650.0,-0.3,0.02,1200\n",
        encoding="utf-8",
    )

    catalog = list_distribution_fields_tool(path=str(candidate_path), adapter="instanovo")
    assert "source_columns.p_value" in catalog["source_numeric_fields"]

    summary = summarize_numeric_field_tool(path=str(candidate_path), adapter="instanovo", field="confidence")
    assert summary["count"] == 2

    grouped = summarize_numeric_field_tool(path=str(candidate_path), adapter="instanovo", field="source_columns.abundance", group_by="experiment_name")
    assert grouped["groups"]["runA"]["mean"] == 1000.0

    weighted = summarize_numeric_field_tool(
        path=str(candidate_path),
        adapter="instanovo",
        field="precursor_mz",
        weight_field="source_columns.abundance",
    )
    assert weighted["weight_field"] == "source_columns.abundance"
    assert weighted["total_weight"] == 2200.0

    binned = bin_numeric_field_tool(path=str(candidate_path), adapter="instanovo", field="precursor_mz", bin_edges=[400.0, 600.0, 800.0])
    assert binned["bins"][0]["count"] == 1
    assert binned["bins"][1]["count"] == 1

    weighted_binned = bin_numeric_field_tool(
        path=str(candidate_path),
        adapter="instanovo",
        field="precursor_mz",
        bin_edges=[400.0, 600.0, 800.0],
        weight_field="source_columns.abundance",
    )
    assert weighted_binned["bins"][0]["weight_sum"] == 1000.0

    categorical = summarize_categorical_field_tool(path=str(candidate_path), adapter="instanovo", field="precursor_charge")
    assert categorical["counts"] == {"2": 1, "3": 1}


def test_missing_cv_2d_and_dashboard_summaries():
    candidates = [
        CandidatePSM(candidate_id="c1", source_engine="instanovo", peptidoform="PEPTIDE", experiment_name="runA", scan_number=1, precursor_mz=500.0, precursor_charge=2, confidence=0.9, source_columns={"retention_time_minutes": 10.0, "abundance": 100.0}),
        CandidatePSM(candidate_id="c2", source_engine="instanovo", peptidoform="PEPTIDE", experiment_name="runA", scan_number=2, precursor_mz=700.0, precursor_charge=3, confidence=0.7, source_columns={"retention_time_minutes": 20.0, "abundance": 120.0}),
        CandidatePSM(candidate_id="c3", source_engine="instanovo", peptidoform="TESTIDE", experiment_name="runB", scan_number=3, precursor_mz=900.0, precursor_charge=3, confidence=0.5, source_columns={"retention_time_minutes": None, "abundance": 80.0}),
    ]
    frame = candidate_models_to_frame(candidates)

    missing = summarize_missing_values(frame, fields=["confidence", "source_columns.retention_time_minutes"])
    assert missing["fields"][0]["missing_count"] == 0
    assert missing["fields"][1]["missing_count"] == 1

    cv = summarize_cv_by_category(frame, value_field="source_columns.abundance", category_field="experiment_name")
    run_a = next(row for row in cv["rows"] if row["category"] == "runA")
    assert run_a["count"] == 2
    assert run_a["cv"] is not None

    grid = bin_2d_numeric_distribution(
        frame,
        x_field="precursor_mz",
        y_field="confidence",
        x_bin_edges=[400.0, 800.0, 1000.0],
        y_bin_edges=[0.0, 0.6, 1.0],
        weight_field="source_columns.abundance",
    )
    assert len(grid["cells"]) == 4
    assert grid["weight_field"] == "source_columns.abundance"

    dashboard = build_dashboard_summary(
        frame,
        rt_field="source_columns.retention_time_minutes",
        quant_field="source_columns.abundance",
        category_field="experiment_name",
    )
    assert dashboard["overview"]["total_psm_count"] == 3
    assert dashboard["ids_per_file"]["counts"] == {"runA": 2, "runB": 1}
    assert "confidence_vs_mz" in dashboard["grids_2d"]
    assert "mz_vs_rt" in dashboard["grids_2d"]


def test_dashboard_tools_on_candidate_table(tmp_path):
    candidate_path = tmp_path / "preds.csv"
    candidate_path.write_text(
        "experiment_name,scan_number,preds,log_probs,precursor_charge,precursor_mz,predictions_log_probability_beam_1,retention_time_minutes,abundance\n"
        "runA,1,PEPTIDE,-0.1,2,500.0,-0.2,10.0,100\n"
        "runA,2,TESTIDE,-0.2,3,700.0,-0.3,20.0,120\n"
        "runB,3,ALTIDE,-0.3,3,900.0,-0.4,,80\n",
        encoding="utf-8",
    )

    missing = summarize_missing_values_tool(path=str(candidate_path), adapter="instanovo", fields=["source_columns.retention_time_minutes"])
    assert missing["fields"][0]["missing_count"] == 1

    cv = summarize_cv_by_category_tool(
        path=str(candidate_path),
        adapter="instanovo",
        value_field="source_columns.abundance",
        category_field="experiment_name",
    )
    assert cv["rows"][0]["category"] == "runA"

    grid = bin_2d_numeric_fields_tool(
        path=str(candidate_path),
        adapter="instanovo",
        x_field="precursor_mz",
        y_field="confidence",
        x_bin_edges=[400.0, 800.0, 1000.0],
        y_bin_edges=[0.0, 0.6, 1.0],
        weight_field="source_columns.abundance",
    )
    assert grid["weight_field"] == "source_columns.abundance"

    dashboard = summarize_dashboard_tool(
        path=str(candidate_path),
        adapter="instanovo",
        rt_field="source_columns.retention_time_minutes",
        quant_field="source_columns.abundance",
    )
    assert dashboard["overview"]["distinct_run_count"] == 2
    assert "retention_time" in dashboard["histograms"]


def test_summarize_and_render_prediction_accuracy(tmp_path):
    candidates = [
        # exact match
        CandidatePSM(candidate_id="c1", source_engine="instanovo", peptidoform="PEPTIDE", precursor_charge=2, confidence=0.95, source_columns={"target": "PEPTIDE", "mass_error": 0.01}),
        # I/L-equivalent only (PEPTLDE vs PEPTIDE)
        CandidatePSM(candidate_id="c2", source_engine="instanovo", peptidoform="PEPTLDE", precursor_charge=2, confidence=0.7, source_columns={"target": "PEPTIDE", "mass_error": -0.02}),
        # mod-only difference: stripped sequence matches
        CandidatePSM(candidate_id="c3", source_engine="instanovo", peptidoform="M[UNIMOD:35]ACK", precursor_charge=3, confidence=0.85, source_columns={"target": "MACK", "mass_error": 0.0}),
        # wrong
        CandidatePSM(candidate_id="c4", source_engine="instanovo", peptidoform="WRONG", precursor_charge=3, confidence=0.2, source_columns={"target": "RIGHT", "mass_error": 1.5}),
        # missing target → skipped
        CandidatePSM(candidate_id="c5", source_engine="instanovo", peptidoform="ANYK", precursor_charge=2, confidence=0.5, source_columns={"mass_error": 0.0}),
    ]
    frame = candidate_models_to_frame(candidates)

    result = summarize_prediction_accuracy(
        frame,
        target_field="source_columns.target",
        mass_error_field="source_columns.mass_error",
    )

    overview = result["overview"]
    assert overview["psm_count"] == 4
    assert overview["skipped_missing_pred_or_target"] == 1
    assert overview["exact_accuracy"] == 0.25
    assert overview["stripped_sequence_accuracy"] == 0.5
    assert overview["il_equivalent_accuracy"] == 0.75

    by_charge = {row["charge"]: row for row in result["by_charge"]}
    assert by_charge[2]["accuracy_il"] == 1.0
    assert by_charge[3]["accuracy_il"] == 0.5

    op_threshold = {row["threshold"]: row for row in result["operating_points"]}
    high = op_threshold[0.9]
    assert high["kept"] == 1 and high["accuracy_il"] == 1.0
    assert op_threshold[0.0]["kept"] == 4 and op_threshold[0.0]["accuracy_il"] == 0.75

    by_mod = {row["label"]: row for row in result["by_target_modification_state"]}
    assert "modified" not in by_mod  # nothing in this fixture has a target with [UNIMOD:..]; targets are plain
    assert by_mod["unmodified"]["count"] == 4

    assert result["confidence_distribution_by_correctness"]["correct_il"]["bins"]
    assert result["mass_error_distribution_by_correctness"]["correct_il"]["bins"]

    output = tmp_path / "accuracy.html"
    html = render_prediction_accuracy_html(result, title="Accuracy", output_path=output)
    assert output.exists()
    assert "<div class='value'>4</div>" in html
    assert "Reliability diagram" in html
    assert "Operating points" in html
    assert "Accuracy by precursor charge" in html

    tool_result = render_prediction_accuracy_html_tool(result, str(tmp_path / "accuracy_tool.html"), title="Tool")
    assert tool_result["format"] == "html"


def test_summarize_prediction_accuracy_tool_loads_from_csv(tmp_path):
    csv_path = tmp_path / "preds.csv"
    csv_path.write_text(
        "experiment_name,scan_number,preds,log_probs,target,precursor_charge\n"
        "run,1,PEPTIDE,-0.05,PEPTIDE,2\n"
        "run,2,PEPTLDE,-0.30,PEPTIDE,2\n"
        "run,3,WRONG,-3.00,RIGHT,3\n",
        encoding="utf-8",
    )
    summary = summarize_prediction_accuracy_tool(
        path=str(csv_path),
        target_field="source_columns.target",
        adapter="instanovo",
    )
    assert summary["overview"]["psm_count"] == 3
    assert summary["overview"]["il_equivalent_accuracy"] == 2 / 3


def test_overlap_modification_quant_and_file_coverage_summaries():
    candidates = [
        CandidatePSM(candidate_id="c1", source_engine="instanovo", peptidoform="PEPTIDE", experiment_name="runA", scan_number=1, precursor_charge=2, confidence=0.9, precursor_mz=500.0, source_columns={"abundance": 100.0, "rt_minutes": 10.0}),
        CandidatePSM(candidate_id="c2", source_engine="instanovo", peptidoform="M[UNIMOD:35]PEPTIDE", experiment_name="runA", scan_number=2, precursor_charge=2, confidence=0.8, precursor_mz=520.0, source_columns={"abundance": 120.0, "rt_minutes": 11.0}),
        CandidatePSM(candidate_id="c3", source_engine="instanovo", peptidoform="PEPTIDE", experiment_name="runB", scan_number=3, precursor_charge=3, confidence=0.7, precursor_mz=700.0, source_columns={"abundance": 110.0, "rt_minutes": 20.0}),
        CandidatePSM(candidate_id="c4", source_engine="instanovo", peptidoform="ALTIDE", experiment_name="runB", scan_number=4, precursor_charge=2, confidence=0.6, precursor_mz=800.0, source_columns={"abundance": 130.0, "rt_minutes": 21.0}),
    ]
    frame = candidate_models_to_frame(candidates)

    overlap = summarize_overlap(frame, entity="peptide", group_field="run_key")
    assert overlap["total_distinct_entities"] == 3
    assert overlap["pairwise"][0]["intersection_count"] == 1

    mods = summarize_modifications(frame, group_field="run_key")
    assert mods["modified_psm_count"] == 1
    assert mods["modification_counts"]["UNIMOD:35"] == 1
    assert mods["residue_counts"]["M"] == 1

    quant = summarize_quant_quality(frame, quant_field="abundance", entity="peptide", group_field="run_key")
    assert quant["overview"]["entity_count"] == 3
    assert quant["overview"]["complete_entity_count"] == 1

    coverage = summarize_file_coverage(frame, group_field="run_key", rt_field="rt_minutes", quant_field="abundance")
    assert coverage["files"]["runA"]["overview"]["total_psm_count"] == 2
    assert "rt_histogram" in coverage["files"]["runA"]
    assert coverage["files"]["runB"]["charge_counts"] == {"2": 1, "3": 1}


def test_new_summary_tools_on_candidate_table(tmp_path):
    candidate_path = tmp_path / "preds.csv"
    candidate_path.write_text(
        "experiment_name,scan_number,preds,log_probs,precursor_charge,precursor_mz,predictions_log_probability_beam_1,abundance,rt_minutes\n"
        "runA,1,PEPTIDE,-0.1,2,500.0,-0.2,100,10.0\n"
        "runA,2,M[UNIMOD:35]PEPTIDE,-0.2,2,-0.3,520.0,120,11.0\n"
        "runB,3,PEPTIDE,-0.3,3,700.0,-0.4,110,20.0\n"
        "runB,4,ALTIDE,-0.4,2,800.0,-0.5,130,21.0\n",
        encoding="utf-8",
    )

    overlap = summarize_overlap_tool(path=str(candidate_path), adapter="instanovo", entity="peptide")
    assert overlap["pairwise"][0]["intersection_count"] == 1

    mods = summarize_modifications_tool(path=str(candidate_path), adapter="instanovo")
    assert mods["modification_counts"]["UNIMOD:35"] == 1

    quant = summarize_quant_quality_tool(path=str(candidate_path), adapter="instanovo", quant_field="abundance", entity="peptide")
    assert quant["overview"]["complete_entity_count"] == 1

    coverage = summarize_file_coverage_tool(path=str(candidate_path), adapter="instanovo", rt_field="rt_minutes", quant_field="abundance")
    assert coverage["files"]["runA"]["overview"]["total_psm_count"] == 2


def test_dataset_summary_with_a_multiplexed_sdrf(tmp_path):
    # A multiplexed run's PSMs belong to all of its channels; the sample metadata lists them all (#81).
    scored_path = tmp_path / "scores.csv"
    _scored_frame().to_csv(scored_path, index=False)
    sdrf_path = tmp_path / "tmt.sdrf.tsv"
    sdrf_path.write_text(
        "assay name\tcomment[data file]\tsource name\tcomment[label]\n"
        "runA_126\trunA.raw\tsampleA1\tTMT126\n"
        "runA_127N\trunA.raw\tsampleA2\tTMT127N\n"
        "runB_126\trunB.raw\tsampleB1\tTMT126\n",
        encoding="utf-8",
    )

    summary = summarize_dataset_tool(scored_path=str(scored_path), output_prefix=str(tmp_path / "dataset"), sdrf_path=str(sdrf_path))
    samples = pd.read_csv(summary["output_files"]["sample_summary_csv"])

    by_sample = samples.set_index("sample_key")
    assert by_sample.loc["runA_126; runA_127N", "label"] == "TMT126; TMT127N"
    assert by_sample.loc["runA_126; runA_127N", "total_psm_count"] == 2
    assert by_sample.loc["runB_126", "label"] == "TMT126"


def _best_row_by_sorting_each_group(group: pd.DataFrame) -> pd.Series:
    """The per-group sort that _best_row used before #68, kept as the reference."""
    sortable = group.copy()
    keys = {
        "__has_confidence": sortable["confidence"].notna().astype(int),
        "__confidence": sortable["confidence"].fillna(float("-inf")),
        "__has_raw_score": sortable["raw_score"].notna().astype(int),
        "__raw_score": sortable["raw_score"].fillna(float("-inf")),
        "__spectral_angle": sortable["spectral_angle"].fillna(float("-inf")),
        "__matched_ion_count": sortable["matched_ion_count"].fillna(float("-inf")),
        "__candidate_id_key": sortable["candidate_id"].fillna("").astype(str),
    }
    for name, values in keys.items():
        sortable[name] = values
    ordered = sortable.sort_values(
        by=[*list(keys)[:6], "__row_order", "__candidate_id_key"],
        ascending=[False, False, False, False, False, False, True, True],
        kind="mergesort",
    )
    return ordered.iloc[0]


def test_best_rows_match_a_per_group_sort_with_ties_and_missing_scores():
    import numpy as np

    from proteomics_mcp.core.dataset import _best_row, _summarize_group, _with_best_rank

    rng = np.random.default_rng(68)
    size = 600
    frame = pd.DataFrame(
        {
            "candidate_id": rng.choice(["a", "b", "c", None], size=size),
            "peptidoform": rng.choice(["PEPTIDE", "PEPM[UNIMOD:35]TIDE", "PEPTIDEK", "PEPTIDEK[UNIMOD:21]"], size=size),
            "confidence": rng.choice([np.nan, 0.5, 0.9], size=size),
            "raw_score": rng.choice([np.nan, 1.0, 2.0], size=size),
            "spectral_angle": rng.choice([np.nan, 0.3, 0.6], size=size),
            "matched_ion_count": rng.choice([np.nan, 4.0, 7.0], size=size),
            "experiment_name": "run",
            "group": rng.integers(0, 40, size=size),
        }
    )
    frame["__row_order"] = rng.permutation(size)  # row order need not follow the index
    frame["stripped_sequence"] = frame["peptidoform"].str.replace(r"\[[^\]]*\]", "", regex=True)
    frame["is_modified"] = frame["peptidoform"].str.contains("[", regex=False)
    frame["run_key"] = "run"
    ranked = _with_best_rank(frame)

    for _, group in ranked.groupby("group"):
        for subset in (group, group[group["is_modified"]], group[~group["is_modified"]]):
            if subset.empty:
                continue
            expected = _best_row_by_sorting_each_group(subset)
            assert _best_row(subset).name == expected.name
            assert _best_row(subset.drop(columns="__best_rank")).name == expected.name
        summary = _summarize_group(group, include_peptide_context=True)
        modified, unmodified = group[group["is_modified"]], group[~group["is_modified"]]
        assert summary["best_peptidoform"] == _best_row_by_sorting_each_group(group)["peptidoform"]
        assert summary["best_modified_peptidoform"] == (_best_row_by_sorting_each_group(modified)["peptidoform"] if not modified.empty else None)
        assert summary["best_unmodified_peptidoform"] == (_best_row_by_sorting_each_group(unmodified)["peptidoform"] if not unmodified.empty else None)


SCRIPT_BREAKOUT = "</script><script>alert(1)</script>"
IMG_ONERROR = '<img src=x onerror="alert(2)">'


def _assert_no_injected_markup(html: str) -> None:
    assert SCRIPT_BREAKOUT not in html
    assert IMG_ONERROR not in html
    assert "<script>alert(1)" not in html
    assert "<img src=x" not in html


def test_html_reports_escape_untrusted_strings(tmp_path):
    candidates = [
        CandidatePSM(candidate_id="c1", source_engine="instanovo", peptidoform="PEPTIDE", experiment_name=SCRIPT_BREAKOUT, scan_number=1, confidence=0.9, precursor_charge=2),
        CandidatePSM(candidate_id="c2", source_engine="instanovo", peptidoform="PEPTIDEK", experiment_name=IMG_ONERROR, scan_number=2, confidence=0.8, precursor_charge=3),
    ]
    inspection = summarize_candidate_inspection_frame(candidate_models_to_frame(candidates))
    for render in (render_psm_table_html, render_peptide_table_html):
        output = tmp_path / f"{render.__name__}.html"
        render(inspection, output)
        html = output.read_text(encoding="utf-8")
        _assert_no_injected_markup(html)
    assert "&lt;img src=x onerror=&quot;alert(2)&quot;&gt;" in (tmp_path / "render_psm_table_html.html").read_text(encoding="utf-8")

    frame = pd.DataFrame([{"run": SCRIPT_BREAKOUT, "note": IMG_ONERROR}])
    html = render_table_browse_html(frame, title=IMG_ONERROR, output_path=tmp_path / "browse.html", csv_filename=SCRIPT_BREAKOUT)
    _assert_no_injected_markup(html)
    assert "<title>&lt;img src=x" in html
    assert "\\u003c/script\\u003e\\u003cscript\\u003ealert(1)" in html  # embedded CSV and file name stay intact as JSON

    summary = {
        IMG_ONERROR: {SCRIPT_BREAKOUT: IMG_ONERROR},
        "rows": [{"run": SCRIPT_BREAKOUT, "count": 1}],
        "values": [IMG_ONERROR],
        "scalar": SCRIPT_BREAKOUT,
    }
    html = render_summary_html(summary, title=SCRIPT_BREAKOUT, output_path=tmp_path / "summary.html")
    _assert_no_injected_markup(html)
    assert "<h2>&lt;img src=x" in html
    assert "<td>&lt;/script&gt;&lt;script&gt;alert(1)&lt;/script&gt;</td>" in html


def test_dashboard_and_accuracy_html_escape_untrusted_strings(tmp_path):
    summary = {
        "dashboard": {
            "overview": {"total_psm_count": 2, "distinct_run_count": 2},
            "ids_per_file": {"counts": {SCRIPT_BREAKOUT: 1, IMG_ONERROR: 1}},
            "scatter_2d": {
                "confidence_vs_spectral_angle": {
                    "points": [
                        {"x": 0.9, "y": 0.8, "group": SCRIPT_BREAKOUT, "candidate_id": "c1", "peptidoform": IMG_ONERROR, "run_key": SCRIPT_BREAKOUT, "scan_number": 1},
                    ]
                }
            },
        },
        "modifications": {"modified_fraction": 0.5, "distinct_modified_peptidoform_count": 1, "modification_counts": {IMG_ONERROR: 1}, "group_modification_counts": {SCRIPT_BREAKOUT: {IMG_ONERROR: 1}}},
    }
    html = render_dashboard_bundle_html(summary, title=SCRIPT_BREAKOUT, output_path=tmp_path / "dashboard.html")
    _assert_no_injected_markup(html)
    assert "<h1>&lt;/script&gt;" in html
    assert "\\u003cimg src=x" in html  # Plotly labels are kept, as escaped JSON

    candidate_path = tmp_path / "preds.csv"
    pd.DataFrame(
        [
            {"experiment_name": SCRIPT_BREAKOUT, "scan_number": 1, "preds": "PEPTIDE", "log_probs": -0.1, "target": "PEPTIDE", "precursor_charge": 2},
            {"experiment_name": IMG_ONERROR, "scan_number": 2, "preds": "PEPTIDEK", "log_probs": -0.2, "target": "PEPTLDEK", "precursor_charge": 3},
        ]
    ).to_csv(candidate_path, index=False)
    result = summarize_prediction_accuracy_tool(path=str(candidate_path), target_field="target", adapter="instanovo")
    html = render_prediction_accuracy_html(result, title=IMG_ONERROR, output_path=tmp_path / "accuracy.html")
    _assert_no_injected_markup(html)
    assert "<title>&lt;img src=x" in html


def test_binning_counts_the_maximum_when_the_range_comes_from_the_data():
    """The last bin is closed when the edges come from the data (numpy.histogram's convention): before, the values at
    the maximum fell outside every half-open bin, e.g. 1 of 100 PSMs in the Figure 4 spectral-angle histogram."""
    frame = pd.DataFrame({"value": [0.1, 0.3, 0.7]})
    binned = bin_numeric_distribution(frame, field="value", bin_count=20)
    assert sum(row["count"] for row in binned["bins"]) == 3
    assert binned["bin_edges"][-1] == 0.7 and binned["bins"][-1]["count"] == 1

    figure4_extremes = pd.DataFrame({"value": [0.2025552246253021, 0.9617314508392228]})
    assert sum(row["count"] for row in bin_numeric_distribution(figure4_extremes, field="value", bin_count=20)["bins"]) == 2

    ties = pd.DataFrame({"value": [0.5, 1.0, 1.0, 1.0]})  # e.g. many PSMs with matched fraction 1.0
    assert bin_numeric_distribution(ties, field="value", bin_count=10)["bins"][-1]["count"] == 3
    assert bin_numeric_distribution(ties, field="value", bin_count=10, value_range=(0.0, 1.0))["bins"][-1]["count"] == 3


@pytest.mark.parametrize("seed", range(20))
def test_binning_counts_every_finite_value(seed):
    rng = np.random.default_rng(seed)
    values = np.concatenate([rng.normal(size=200) * rng.uniform(1e-6, 1e6), [np.nan, np.inf, -np.inf]])
    frame = pd.DataFrame({"x": values, "y": rng.uniform(0.0, 1.0, size=len(values))})
    finite = int(np.isfinite(values).sum())
    binned = bin_numeric_distribution(frame.replace([np.inf, -np.inf], np.nan), field="x", bin_count=int(rng.integers(1, 40)))
    assert sum(row["count"] for row in binned["bins"]) == finite
    cells = bin_2d_numeric_distribution(frame.replace([np.inf, -np.inf], np.nan), x_field="x", y_field="y", x_bin_count=7, y_bin_count=5)["cells"]
    assert sum(cell["count"] for cell in cells) == finite


def test_explicit_bin_edges_stay_half_open():
    frame = pd.DataFrame({"value": [1.0, 2.0, 3.0]})
    binned = bin_numeric_distribution(frame, field="value", bin_edges=[1.0, 2.0, 3.0])
    assert [row["count"] for row in binned["bins"]] == [1, 1]  # 3.0 is outside [2, 3): the caller chose the edges


# A range given at one end only takes the other from the data (#175)


@pytest.mark.parametrize(
    ("value_range", "edges"),
    [
        ((0.0, None), (0.0, 0.7)),
        ((None, 1.0), (0.1, 1.0)),
        ((0.0, 1.0), (0.0, 1.0)),
        (None, (0.1, 0.7)),
    ],
)
def test_a_one_sided_range_takes_the_other_end_from_the_data(value_range, edges):
    frame = pd.DataFrame({"value": [0.1, 0.3, 0.7]})

    binned = bin_numeric_distribution(frame, field="value", bin_count=4, value_range=value_range)

    assert (binned["bin_edges"][0], binned["bin_edges"][-1]) == pytest.approx(edges)
    assert sum(row["count"] for row in binned["bins"]) == 3


@pytest.mark.parametrize(("value_range", "message"), [((0.9, None), "start at 0.9 and end at 0.7"), ((0.5, 0.2), "start at 0.5 and end at 0.2"), ((0.0,), "must be \\[minimum, maximum\\]")])
def test_a_range_that_cannot_be_binned_is_refused(value_range, message):
    frame = pd.DataFrame({"value": [0.1, 0.3, 0.7]})

    with pytest.raises(ValueError, match=message):
        bin_numeric_distribution(frame, field="value", bin_count=4, value_range=value_range)


def test_range_min_alone_reaches_the_bins_through_mcp_and_the_cli(tmp_path, capsys):
    import asyncio

    from proteomics_mcp.cli import main
    from proteomics_mcp.server import mcp

    path = tmp_path / "preds.csv"
    pd.DataFrame({"scan_number": [1, 2, 3], "preds": ["PEPTIDEK"] * 3, "log_probs": [-0.1, -0.2, -0.3],
                  "precursor_mz": [450.0, 500.0, 550.0]}).to_csv(path, index=False)

    async def call():
        from mcp import Client

        async with Client(mcp) as client:
            return await client.call_tool("bin_numeric_field", {"path": str(path), "field": "precursor_mz", "bin_count": 4, "range_min": 400.0})

    result = asyncio.run(call())
    assert not result.is_error, result.content[0].text
    assert json.loads(result.content[0].text)["bin_edges"] == [400.0, 437.5, 475.0, 512.5, 550.0]

    main(["bin-numeric-field", "--path", str(path), "--field", "precursor_mz", "--bin-count", "4", "--range-json", "[null, 600]"])
    assert json.loads(capsys.readouterr().out)["bin_edges"] == [450.0, 487.5, 525.0, 562.5, 600.0]
