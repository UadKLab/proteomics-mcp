from __future__ import annotations

import math

from proteomics_mcp.adapters.instanovo import load_instanovo_predictions, normalize_peptidoform_to_instanovo_style


def test_legacy_instanovo_confidence(tmp_path):
    path = tmp_path / "legacy.csv"
    path.write_text(
        "experiment_name,scan_number,preds,log_probs\n"
        "run,42,PEPTIDE,-0.5\n",
        encoding="utf-8",
    )

    [candidate] = load_instanovo_predictions(path)

    assert candidate.source_engine == "instanovo"
    assert candidate.peptidoform == "PEPTIDE"
    assert candidate.scan_number == 42
    assert candidate.raw_score == -0.5
    assert abs(candidate.confidence - math.exp(-0.5)) < 1e-12
    assert candidate.source_columns["schema"] == "legacy_simple"


def test_instanovo_default_experiment_name_for_minimal_export(tmp_path):
    path = tmp_path / "minimal.csv"
    path.write_text("spectrum_index,preds,log_probs\n42,PEPTIDE,-0.5\n", encoding="utf-8")

    [candidate] = load_instanovo_predictions(path, default_experiment_name="runA")

    assert candidate.experiment_name == "runA"


def test_singular_prediction_log_probability_schema(tmp_path):
    path = tmp_path / "singular.csv"
    path.write_text(
        "prediction_id,precursor_mz,precursor_charge,prediction,log_probability\n"
        "0,500.0,2,PEPTIDE,-0.25\n",
        encoding="utf-8",
    )

    [candidate] = load_instanovo_predictions(path)

    assert candidate.source_engine == "instanovo"
    assert candidate.peptidoform == "PEPTIDE"
    assert candidate.precursor_mz == 500.0
    assert candidate.precursor_charge == 2
    assert candidate.raw_score == -0.25
    assert abs(candidate.confidence - math.exp(-0.25)) < 1e-12


def test_instanovo_normalizes_parenthesized_oxidation(tmp_path):
    path = tmp_path / "oxidized.csv"
    path.write_text("spectrum_index,preds,log_probs\n1,PEPM(ox)TIDE,-0.1\n", encoding="utf-8")

    [candidate] = load_instanovo_predictions(path)

    assert candidate.peptidoform == "PEPM[UNIMOD:35]TIDE"


def test_beam_output_token_confidence_summary(tmp_path):
    path = tmp_path / "beam.csv"
    path.write_text(
        "experiment_name,scan_number,spectrum_id,precursor_mz,precursor_charge,prediction_id,predictions,log_probs,token_log_probs,predictions_tokenised,predictions_beam_0,predictions_log_probability_beam_0\n"
        'run,7,run:7,500.1,2,abc,ACD,-0.1,"[-0.1, -2.0, -0.3]","A, C, D",ACD,-0.1\n',
        encoding="utf-8",
    )

    [candidate] = load_instanovo_predictions(path)

    assert candidate.source_columns["schema"] == "instanovo_beam"
    assert candidate.precursor_mz == 500.1
    assert candidate.precursor_charge == 2
    assert math.isclose(candidate.token_log_probability_mean, (-0.1 - 2.0 - 0.3) / 3)
    assert candidate.token_log_probability_min == -2.0
    assert candidate.token_confidences == [math.exp(-0.1), math.exp(-2.0), math.exp(-0.3)]
    assert candidate.token_confidence_min_position == 2
    assert candidate.token_confidence_min_token == "C"
    assert candidate.token_confidence_mean is not None


def test_beam_output_derives_next_beam_from_second_beam_score(tmp_path):
    path = tmp_path / "beam_next.csv"
    path.write_text(
        "experiment_name,scan_number,predictions,log_probs,predictions_beam_0,predictions_log_probability_beam_0,predictions_beam_1,predictions_log_probability_beam_1\n"
        "run,7,ACD,-0.1,ACD,-0.1,ACE,-0.8\n",
        encoding="utf-8",
    )

    [candidate] = load_instanovo_predictions(path)

    assert math.isclose(candidate.next_beam or 0.0, math.exp(-0.8))
    assert math.isclose(candidate.next_beam_margin or 0.0, math.exp(-0.1) - math.exp(-0.8))


def test_documented_combined_output_defaults_to_predictions(tmp_path):
    path = tmp_path / "combined.csv"
    path.write_text(
        "scan_number,precursor_mz,precursor_charge,experiment_name,spectrum_id,predictions,log_probs,instanovo_predictions,instanovo_log_probabilities,diffusion_predictions,diffusion_log_probabilities\n"
        "1,451.2,2,spectra,spectra:1,FINAL,-0.01,BASE,-1.0,REFINED,-0.2\n",
        encoding="utf-8",
    )

    [candidate] = load_instanovo_predictions(path)

    assert candidate.source_columns["schema"] == "documented_combined"
    assert candidate.peptidoform == "FINAL"
    assert candidate.raw_score_name == "log_probs"


def test_empty_predictions_are_skipped(tmp_path):
    path = tmp_path / "empty.csv"
    path.write_text(
        "experiment_name,scan_number,preds,log_probs\n"
        "run,42,,-1.0\n",
        encoding="utf-8",
    )

    assert load_instanovo_predictions(path) == []


def test_legacy_instanovo_accepts_spectrum_index_and_global_index(tmp_path):
    path = tmp_path / "hela_labelled.csv"
    path.write_text(
        "spectrum_index,global_index,targets,preds,log_probs\n"
        "28599,2,ARPFPDGLAEDIDKGEVSAR,ARPFPDGLAEDLDKGEVSAR,-0.0144759882241487\n",
        encoding="utf-8",
    )

    [candidate] = load_instanovo_predictions(path)

    assert candidate.scan_number == 28599
    assert candidate.spectrum_id == "28599"
    assert candidate.candidate_id == "2"
    assert math.isclose(candidate.confidence or 0.0, math.exp(-0.0144759882241487))


def test_legacy_sample_export_maps_zero_based_index_to_mzml_scan_and_run(tmp_path):
    path = tmp_path / "wound.csv"
    path.write_text(
        "spectrum_index,global_index,preds,log_probs,sample\n"
        "57947,2739,PEPTIDE,-0.1,20190905_KMJ_JUT_Louise_HYC_PSM\n",
        encoding="utf-8",
    )

    [candidate] = load_instanovo_predictions(path)

    assert candidate.scan_number == 57948
    assert candidate.spectrum_id == "57948"
    assert candidate.experiment_name == "20190905_KMJ_JUT_Louise_HYC"


def test_legacy_instanovo_accepts_native_id_column(tmp_path):
    path = tmp_path / "hela_raw.csv"
    path.write_text(
        "id,preds,log_probs\n"
        "\"controllerType=0 controllerNumber=1 scan=1235\",PPPEEGPPAHLPRPR,-17.894763946533203\n",
        encoding="utf-8",
    )

    [candidate] = load_instanovo_predictions(path)

    assert candidate.scan_number == 1235
    assert candidate.spectrum_id == "controllerType=0 controllerNumber=1 scan=1235"
    assert candidate.candidate_id == "controllerType=0 controllerNumber=1 scan=1235"


def test_normalize_peptidoform_to_instanovo_style_maps_named_modifications():
    assert normalize_peptidoform_to_instanovo_style("C[Carbamidomethyl]M[Oxidation]") == "C[UNIMOD:4]M[UNIMOD:35]"
    assert normalize_peptidoform_to_instanovo_style("[Acetyl]PEPTIDE") == "[UNIMOD:1]-PEPTIDE"
    assert normalize_peptidoform_to_instanovo_style("M[15.994915]PEPTIDE") == "M[UNIMOD:35]PEPTIDE"


def test_normalize_peptidoform_to_instanovo_style_resolves_unimod_catalog_entries():
    assert normalize_peptidoform_to_instanovo_style("K[Methyl]PEPTIDE") == "K[UNIMOD:34]PEPTIDE"
    assert normalize_peptidoform_to_instanovo_style("[Formyl]PEPTIDE") == "[UNIMOD:122]-PEPTIDE"
