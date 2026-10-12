from __future__ import annotations

import json
from pathlib import Path

import pytest

import proteomics_mcp.core.visualization as visualization
from proteomics_mcp.core.visualization import (
    build_spectrum_payload,
    observed_spectrum_metrics,
    render_annotated_observed_spectrum_html,
    render_observed_spectrum_html,
    render_spectrum_comparison_html,
    render_static_spectrum,
)
from proteomics_mcp.core.scoring import score_candidate_spectrum
from proteomics_mcp.models.schemas import CandidatePSM, ObservedSpectrum, PredictedFragment, PredictedSpectrum, SpectrumMetricResult
from proteomics_mcp.predictors.mock import mock_predicted_spectrum
from proteomics_mcp.predictors.ms2pip import predict_fragments_ms2pip, predict_fragments_ms2pip_batch
from proteomics_mcp.tools.workflows import render_spectra_batch_tool


def test_ms2pip_wrapper_requires_dependency_when_absent():
    try:
        import ms2pip  # noqa: F401
    except ImportError:
        try:
            predict_fragments_ms2pip("PEPTIDE", 2)
        except RuntimeError as exc:
            assert "ms2pip" in str(exc).lower()
        else:
            raise AssertionError("Expected RuntimeError when ms2pip is absent")


@pytest.mark.network
def test_ms2pip_auto_uses_charge2_fragments_for_charge3_precursor():
    try:
        import ms2pip  # noqa: F401
    except ImportError:
        return

    predicted = predict_fragments_ms2pip("PEPTIDE", 3, model="HCD")

    assert predicted.model == "HCDch2"
    assert predicted.metadata["max_fragment_charge"] == 2
    assert max(fragment.charge for fragment in predicted.fragments) == 2
    assert {fragment.charge for fragment in predicted.fragments} == {1, 2}


@pytest.mark.network
def test_ms2pip_explicit_max_fragment_charge_filters_charge2_fragments():
    try:
        import ms2pip  # noqa: F401
    except ImportError:
        return

    predicted = predict_fragments_ms2pip("PEPTIDE", 3, model="HCD", max_fragment_charge=1)

    assert predicted.model == "HCD"
    assert predicted.metadata["max_fragment_charge"] == 1
    assert {fragment.charge for fragment in predicted.fragments} == {1}


@pytest.mark.network
def test_ms2pip_batch_wrapper_retains_request_order():
    try:
        import ms2pip  # noqa: F401
    except ImportError:
        return

    predictions = predict_fragments_ms2pip_batch([("PEPTIDE", 2), ("PEPTIDER", 3)], model="HCD", processes=1)

    assert [prediction.peptidoform for prediction in predictions] == ["PEPTIDE", "PEPTIDER"]
    assert [prediction.precursor_charge for prediction in predictions] == [2, 3]
    assert max(fragment.charge for fragment in predictions[1].fragments) == 2


def test_render_html_contains_plotly_download_config(tmp_path):
    observed = ObservedSpectrum(
        mz=[100.0],
        intensity=[10.0],
        rt_seconds=125.0,
        metadata={
            "attributes": [
                {"accession": "MS:1000422", "name": "higher energy beam-type collision-induced dissociation"},
                {"accession": "MS:1000045", "name": "collision energy", "value": 30, "unitName": "electronvolt"},
                {"accession": "MS:1000463", "name": "instrument", "value": "Orbitrap Fusion Lumos"},
                {"accession": "MS:1000484", "name": "orbitrap"},
            ]
        },
    )
    predicted = mock_predicted_spectrum("PEPTIDE", 2)
    payload = build_spectrum_payload(observed, predicted)
    output = tmp_path / "plot.html"

    html = render_spectrum_comparison_html(payload, output)

    assert "Plotly.newPlot" in html
    assert "toImageButtonOptions" in html
    assert "modeBarButtons: [" in html
    assert "['zoom2d', 'pan2d']" in html
    assert "['zoomIn2d', 'zoomOut2d', 'autoScale2d']" in html
    assert "'select2d'" not in html
    assert "'lasso2d'" not in html
    assert "'resetScale2d'" not in html
    assert "document.getElementById('download-svg')" not in html
    assert "document.getElementById('download-pdf')" not in html
    assert "name: 'Download SVG'" not in html
    assert "name: 'Print / Save PDF'" not in html
    assert 'id="save-spectrum-pdf"' not in html
    assert "Save plot as PDF" not in html
    assert "name: 'Print page / Save PDF'" in html
    assert "name: 'Toggle dark mode'" in html
    assert "prefers-color-scheme: dark" in html
    assert "#111827" in html
    assert "#9792E3" in html
    assert "#B8336A" in html
    assert "showgrid: false" in html
    assert "tickformat: '.2f'" in html
    assert "Relative intensity; observed positive, predicted mirrored below zero" in html
    assert "predicted intensity=%{customdata:.3e}" in html
    assert "Predicted unmatched labels" not in html
    assert "hover area" in html
    assert '"marker": {"color": "rgba(0,0,0,0)", "size": 12}' in html
    assert "hoverlabel:" in html
    assert "hoverBg" in html
    assert "rgba(255, 255, 255, 0.96)" in html
    assert output.read_text(encoding="utf-8") == html


def test_fragment_sequence_toggle_uses_compact_modified_sequence(tmp_path):
    observed = ObservedSpectrum(mz=[100.0], intensity=[10.0])
    predicted = PredictedSpectrum(
        peptidoform="PEPM[UNIMOD:35]TIDE",
        precursor_charge=2,
        fragments=[PredictedFragment(mz=100.0, intensity=1.0, annotation="b4+", ion_type="b", ordinal=4, charge=1)],
        predictor_name="mock",
        model="mock",
    )
    candidate = CandidatePSM(candidate_id="c1", source_engine="test", peptidoform=predicted.peptidoform, precursor_charge=2)
    result = score_candidate_spectrum(candidate=candidate, observed=observed, predicted=predicted)
    html = render_spectrum_comparison_html(build_spectrum_payload(observed, predicted, result), tmp_path / "modified.html")

    assert "PEPM(ox)" in html


def test_render_observed_spectrum_contains_observed_only_metrics(tmp_path):
    observed = ObservedSpectrum(
        mz=[100.0, 101.0, 200.0],
        intensity=[10.0, 30.0, 5.0],
        precursor_mz=100.5,
        metadata={"attributes": [{"accession": "MS:1000024", "name": "resolving power", "value": 60000}]},
    )

    metrics = observed_spectrum_metrics(observed.model_dump())
    html = render_observed_spectrum_html(observed.model_dump(), tmp_path / "observed.html")

    assert 'id="save-spectrum-pdf"' not in html
    assert "Save plot as PDF" not in html
    assert "name: 'Print page / Save PDF'" in html

    assert metrics["peak_count"] == 3
    assert metrics["tic"] == 45.0
    assert metrics["base_peak_mz"] == 101.0
    assert "Observed Spectrum" in html
    assert "Prediction Metadata" in html
    assert "Fragmentation method" in html
    assert "Isolation target m/z" in html
    assert "Scan resolution" in html
    assert "60000" in html
    assert "Spectrum USI" in html
    assert "Predicted sequence" not in html
    assert "Predictor" not in html
    assert "Acquisition" not in html
    assert "Base peak m/z" in html
    assert "Precursor region" not in html
    assert "m/z min" not in html
    assert "m/z max" not in html
    assert "Mean intensity" not in html
    assert "Predicted spectrum" not in html
    assert "#068D9D" in html
    assert "observed hover area" in html


def test_render_static_spectrum_writes_mirror_svg(tmp_path):
    observed = ObservedSpectrum(mz=[100.0, 150.0], intensity=[10.0, 20.0], precursor_charge=2)
    predicted = mock_predicted_spectrum("PEPTIDE", 2)
    output = tmp_path / "mirror.svg"

    rendered = render_static_spectrum(build_spectrum_payload(observed, predicted), output, view="mirror", label_top_n=0)

    assert rendered == output
    text = output.read_text(encoding="utf-8")
    assert "<svg" in text
    assert "b-ions observed" in text
    assert "y-ions observed" in text


def test_render_static_spectrum_writes_annotated_svg(monkeypatch, tmp_path):
    def fake_annotate(*_, **kwargs):
        return {
            "annotations": ["b1+", "a1+"],
            "matches": [
                {"observed_index": 0, "annotation": "b1+", "ion_type": "b"},
                {"observed_index": 1, "annotation": "a1+", "ion_type": "a"},
            ],
        }

    monkeypatch.setattr(visualization, "annotate_observed_ions", fake_annotate)
    output = tmp_path / "annotated.svg"
    render_static_spectrum(
        {"observed": ObservedSpectrum(mz=[100.0, 150.0], intensity=[10.0, 20.0], precursor_charge=2).model_dump(), "metrics": {"candidate": {"peptidoform": "PEPTIDE"}}},
        output,
        view="annotated",
        include_fragment_sequences=True,
    )

    text = output.read_text(encoding="utf-8")
    assert "Annotated observed spectrum" in text
    assert "b-ions" in text
    assert "a1+" in text
    assert "b1+ | P" in text
    assert "a1+ | P" not in text


def test_render_annotated_observed_spectrum_supports_configurable_ion_series(monkeypatch, tmp_path):
    def fake_annotate(*_, **kwargs):
        assert kwargs["max_isotope"] == 2
        return {
            "annotations": ["b1+", "y1+", "a1+", ""],
            "annotated_peak_count": 3,
            "ion_type_counts": {"a": 1, "b": 1, "y": 1},
            "ion_types": kwargs["ion_types"],
            "total_theoretical_ions": 8,
        }

    monkeypatch.setattr(visualization, "annotate_observed_ions", fake_annotate)
    observed = ObservedSpectrum(mz=[100.0, 150.0, 200.0, 250.0], intensity=[10.0, 20.0, 30.0, 5.0], precursor_mz=300.0, precursor_charge=2)

    html = render_annotated_observed_spectrum_html(
        {"observed": observed.model_dump(), "metrics": {"candidate": {"peptidoform": "PEPTIDE"}}},
        tmp_path / "annotated.html",
        ion_types="abcxyz",
    )

    assert "Annotated observed spectrum" in html
    assert "Prediction Metadata" in html
    assert "Predicted sequence" in html
    assert "Predictor" not in html
    assert "Fragmentation method" in html
    assert "Isolation target m/z" in html
    assert "Precursor charge" in html
    assert "Spectrum USI" in html
    assert "Scan resolution" not in html
    assert "Spectrum Metadata" not in html
    assert "Observed Annotation" in html
    assert "Peptidoform" not in html
    assert "abcxyz" in html
    assert "Max isotope" not in html
    assert "3 / 8" in html
    assert "b-ions" in html
    assert "y-ions" in html
    assert "a-ions" in html
    assert 'class="ion-ladder-controls"' in html
    assert 'data-ion-type="b"' not in html
    assert 'data-ion-type="y"' not in html
    assert 'data-ion-type="a"' in html
    assert '"visible": false' in html
    assert "b1+" in html
    assert "Predicted spectrum" not in html


def test_render_mirror_spectrum_exposes_fragment_sequence_toggle(tmp_path):
    candidate = CandidatePSM(candidate_id="c1", source_engine="instanovo", peptidoform="PEPTIDE", precursor_charge=2, confidence=0.9)
    observed = ObservedSpectrum(mz=[100.0], intensity=[10.0], precursor_charge=2)
    predicted = PredictedSpectrum(
        peptidoform="PEPTIDE",
        precursor_charge=2,
        predictor_name="mock",
        fragments=[PredictedFragment(mz=100.0, intensity=1.0, annotation="b3+", ion_type="b", ordinal=3)],
    )
    result = SpectrumMetricResult(
        candidate=candidate,
        spectral_angle=0.0,
        cosine_similarity=0.0,
        pearson_correlation=0.0,
        xcorr=0.0,
        log10_xcorr=0.0,
        log2_xcorr=0.0,
        hyperscore=0.0,
        log2_hyperscore=0.0,
        matched_ion_count=1,
        total_predicted_ions=1,
        matched_ion_fraction=1.0,
        explained_intensity_fraction=1.0,
        precursor_signal_fraction=0.0,
        fragment_mass_error_mean_ppm=0.0,
        fragment_mass_error_abs_mean_ppm=0.0,
        fragment_mass_error_mean_da=0.0,
        fragment_mass_error_abs_mean_da=0.0,
        precursor_mass_error_ppm=0.0,
        precursor_mass_error_da=0.0,
        matched_peaks=[{"observed_index": 0, "predicted_index": 0, "annotation": "b3+", "ion_type": "b"}],
    )

    html = render_spectrum_comparison_html(build_spectrum_payload(observed, predicted, result), tmp_path / "mirror.html")

    assert "Show fragment ion sequences" in html
    assert "formatSequenceLabel" in html
    assert "fragmentSequenceText" in html
    assert "PEP" in html


def test_render_annotated_spectrum_exposes_fragment_sequence_toggle(monkeypatch, tmp_path):
    def fake_annotate(*_, **kwargs):
        return {
            "annotations": ["b2+", "y2+"],
            "matches": [
                {"observed_index": 0, "annotation": "b2+", "ion_type": "b"},
                {"observed_index": 1, "annotation": "y2+", "ion_type": "y"},
            ],
            "annotated_peak_count": 2,
            "ion_type_counts": {"b": 1, "y": 1},
            "ion_types": kwargs["ion_types"],
            "total_theoretical_ions": 4,
        }

    monkeypatch.setattr(visualization, "annotate_observed_ions", fake_annotate)
    observed = ObservedSpectrum(mz=[100.0, 150.0], intensity=[10.0, 20.0], precursor_mz=300.0, precursor_charge=2)

    html = render_annotated_observed_spectrum_html(
        {"observed": observed.model_dump(), "metrics": {"candidate": {"peptidoform": "PEPTIDE"}}},
        tmp_path / "annotated_toggle.html",
        ion_types="by",
    )

    assert "Show fragment ion sequences" in html
    assert "formatSequenceLabel" in html
    assert "fragmentSequenceText" in html
    assert "PE" in html
    assert "DE" in html


def test_render_html_contains_metric_tiles_when_metrics_present(tmp_path):
    candidate = CandidatePSM(
        candidate_id="c1",
        source_engine="instanovo",
        peptidoform="PEPTIDE",
        confidence=0.9,
        raw_score=-0.1,
        raw_score_name="log_probability",
        next_beam=-0.5,
        next_beam_margin=0.4,
        token_log_probabilities=[-0.2, -0.4],
        token_log_probability_mean=-0.3,
        token_log_probability_min=-0.4,
        token_confidences=[0.8187307530779818, 0.6703200460356393],
        token_confidence_mean=0.7445253995568105,
        token_confidence_min=0.6703200460356393,
        token_confidence_min_position=2,
        token_confidence_min_token="E",
    )
    observed = ObservedSpectrum(
        mz=[100.0],
        intensity=[10.0],
        rt_seconds=125.0,
        metadata={
            "attributes": [
                {"accession": "MS:1000422", "name": "higher energy beam-type collision-induced dissociation"},
                {"accession": "MS:1000045", "name": "collision energy", "value": 30, "unitName": "electronvolt"},
                {"accession": "MS:1000024", "name": "resolving power", "value": 120000},
                {"accession": "MS:1000463", "name": "instrument", "value": "Orbitrap Fusion Lumos"},
                {"accession": "MS:1000484", "name": "orbitrap"},
            ]
        },
    )
    predicted = PredictedSpectrum(
        peptidoform="PEPTIDE",
        precursor_charge=2,
        predictor_name="mock",
        fragments=[PredictedFragment(mz=100.0, intensity=1.0, annotation="b1+", ion_type="b", ordinal=1)],
    )
    result = SpectrumMetricResult(
        candidate=candidate,
        spectral_angle=0.75,
        cosine_similarity=0.8,
        pearson_correlation=0.7,
        xcorr=10.0,
        log10_xcorr=1.0,
        log2_xcorr=3.321928094887362,
        hyperscore=4.2,
        log2_hyperscore=2.070389327891398,
        matched_ion_count=1,
        total_predicted_ions=2,
        matched_ion_fraction=0.5,
        explained_intensity_fraction=0.6,
        precursor_signal_fraction=0.1,
        fragment_mass_error_mean_ppm=1.0,
        fragment_mass_error_abs_mean_ppm=1.2,
        fragment_mass_error_mean_da=0.001,
        fragment_mass_error_abs_mean_da=0.002,
        precursor_mass_error_ppm=-2.0,
        precursor_mass_error_da=-0.001,
        matched_peaks=[],
        metric_profile="legacy_v1",
        predictor_metadata={},
        annotation_summary={"annotated_peak_count": 6, "ion_types": "abcxyz", "ion_type_counts": {"a": 1, "b": 2, "y": 3}},
    )
    payload = build_spectrum_payload(observed, predicted, result)
    output = tmp_path / "plot.html"

    html = render_spectrum_comparison_html(payload, output)

    assert 'id="metrics-panel"' in html
    assert "Match Quality" in html
    assert "Ion Evidence" in html
    assert "Mass Error" in html
    assert "Candidate" in html
    assert "Prediction Metadata" in html
    assert "Predicted sequence" in html
    assert "Fragmentation method" in html
    assert "Isolation target m/z" in html
    assert "Precursor charge" in html
    assert "Spectrum USI" in html
    assert "metadata-usi" not in html
    assert "word-break: break-word" in html
    assert "Retention time" in html
    assert "2.08 min" in html
    assert "Collision energy" in html
    assert "30 electronvolt" in html
    assert "Scan resolution" in html
    assert "120000" in html
    assert "Instrument" in html
    assert "Orbitrap Fusion Lumos" in html
    assert "Detector" in html
    assert "orbitrap" in html
    assert "Spectral angle" in html
    assert "spectrum similarity" in html
    assert "higher is better" not in html
    assert "rank-supporting" not in html
    assert "carryover window" not in html
    assert "log2 XCorr" in html
    assert "cross-correlation-like" not in html
    assert "log2 Hyperscore" in html
    assert "Matched ions" in html
    assert "1 / 2" in html
    assert "Annotated peaks" not in html
    assert "Observed ion series" in html
    assert "a:1, b:2, y:3" in html
    assert "90.0%" in html
    assert "Next beam" in html
    assert "Next beam margin" in html
    assert "Next score" not in html
    assert "0.400" in html
    assert "67.0%" in html
    assert "60.0%" in html
    assert "instanovo" in html


def test_render_spectra_batch_tool_writes_parallel_category_outputs(tmp_path):
    observed = ObservedSpectrum(mz=[100.0, 150.0], intensity=[10.0, 25.0], precursor_mz=400.0, precursor_charge=2)
    predicted = mock_predicted_spectrum("PEPTIDE", 2)
    payload = build_spectrum_payload(observed, predicted)

    result = render_spectra_batch_tool(
        [
            {"name": "scan_a", "payload": payload},
            {"name": "scan_b", "payload": payload, "annotated_options": {"peptidoform": "PEPTIDE", "ion_types": "by"}},
        ],
        str(tmp_path / "batch"),
        views=["comparison", "observed", "annotated"],
        max_workers=2,
    )

    assert result["items_requested"] == 2
    assert result["items_per_view"] == {"comparison": 2, "observed": 2, "annotated": 2}
    assert result["failures"] == []
    assert Path(result["status_path"]).exists()
    for view in ("comparison", "observed", "annotated"):
        for item in result["results"][view]:
            assert (tmp_path / "batch" / view / f'{item["output_stem"]}.html').exists()


def test_render_spectra_batch_tool_accepts_input_path(tmp_path):
    observed = ObservedSpectrum(mz=[100.0], intensity=[10.0], precursor_mz=400.0, precursor_charge=2)
    predicted = mock_predicted_spectrum("PEPTIDE", 2)
    payload = build_spectrum_payload(observed, predicted)
    input_path = tmp_path / "batch_items.json"
    input_path.write_text(json.dumps({"items": [{"name": "scan_a", "payload": payload}], "views": ["observed"]}), encoding="utf-8")

    result = render_spectra_batch_tool(
        input_path=str(input_path),
        output_dir=str(tmp_path / "batch_out"),
    )

    assert result["input_path"] == str(input_path)
    assert result["views"] == ["observed"]
    assert result["items_per_view"] == {"observed": 1}
    assert (tmp_path / "batch_out" / "observed" / "scan_a.html").exists()


def test_render_spectra_batch_tool_validate_only_reports_schema_failures(tmp_path):
    input_path = tmp_path / "invalid_batch_items.json"
    input_path.write_text(json.dumps({"items": [{"name": "bad_item"}], "views": ["comparison"]}), encoding="utf-8")

    result = render_spectra_batch_tool(
        input_path=str(input_path),
        validate_only=True,
    )

    assert result["format"] == "validation"
    assert result["items_requested"] == 1
    assert result["valid_items"] == 0
    assert result["invalid_items"] == 1
    assert "payload" in result["failures"][0]["error"]


def test_render_html_falls_back_to_precursor_target_metadata(tmp_path):
    candidate = CandidatePSM(
        candidate_id="c1",
        source_engine="instanovo",
        peptidoform="PEPTIDE",
        precursor_charge=2,
        confidence=0.9,
    )
    observed = ObservedSpectrum(
        mz=[100.0],
        intensity=[10.0],
        precursor_charge=2,
        metadata={"attributes": [{"accession": "MS:1000827", "name": "isolation window target m/z", "value": 400.0}]},
    )
    predicted = PredictedSpectrum(
        peptidoform="PEPTIDE",
        precursor_charge=2,
        predictor_name="mock",
        fragments=[PredictedFragment(mz=100.0, intensity=1.0, annotation="b1+", ion_type="b", ordinal=1)],
    )
    result = SpectrumMetricResult(
        candidate=candidate,
        spectral_angle=0.75,
        cosine_similarity=0.8,
        pearson_correlation=0.7,
        xcorr=10.0,
        log10_xcorr=1.0,
        log2_xcorr=3.321928094887362,
        hyperscore=4.2,
        log2_hyperscore=2.070389327891398,
        matched_ion_count=1,
        total_predicted_ions=2,
        matched_ion_fraction=0.5,
        explained_intensity_fraction=0.6,
        precursor_signal_fraction=0.1,
        fragment_mass_error_mean_ppm=1.0,
        fragment_mass_error_abs_mean_ppm=1.2,
        fragment_mass_error_mean_da=0.001,
        fragment_mass_error_abs_mean_da=0.002,
        precursor_mass_error_ppm=None,
        precursor_mass_error_da=None,
        matched_peaks=[],
        metric_profile="legacy_v1",
        predictor_metadata={},
    )
    html = render_spectrum_comparison_html(build_spectrum_payload(observed, predicted, result), tmp_path / "plot.html")

    assert "Isolation target m/z" in html
    assert "400.00000" in html
    assert "Precursor charge" in html
    assert ">2</div>" in html
    assert "Precursor ppm" in html
    assert "Precursor Da" in html
    assert "target minus theoretical" in html
