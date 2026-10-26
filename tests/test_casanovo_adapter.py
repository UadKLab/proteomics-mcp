from __future__ import annotations

import math

import pytest

from proteomics_mcp.adapters.casanovo import load_casanovo_predictions
from proteomics_mcp.tools.workflows import load_candidates_tool


def test_load_casanovo_predictions_normalizes_sequence_and_confidence(tmp_path):
    path = tmp_path / "casanovo.csv"
    path.write_text(
        "prediction_id,predictions,targets,log_probs,predictable,precursor_mz,precursor_charge,group\n"
        "0,DM[Oxidation]HGVTSHLSNQELQDLVEFMK,DM[UNIMOD:35]HGVTSHLSNQELQDLVEFMK,-0.24043669456534025,True,858.73639,3,runA\n",
        encoding="utf-8",
    )

    [candidate] = load_casanovo_predictions(path)

    assert candidate.candidate_id == "0"
    assert candidate.source_engine == "casanovo"
    assert candidate.peptidoform == "DM[UNIMOD:35]HGVTSHLSNQELQDLVEFMK"
    assert candidate.experiment_name == "runA"
    assert candidate.precursor_mz == 858.73639
    assert candidate.precursor_charge == 3
    assert math.isclose(candidate.raw_score or 0.0, -0.24043669456534025)
    assert math.isclose(candidate.confidence or 0.0, math.exp(-0.24043669456534025))
    assert candidate.confidence_type == "sequence_probability_from_log_probability"
    assert candidate.source_columns["schema"] == "casanovo_predictions"


def test_load_casanovo_predictions_keeps_missing_scan_number_empty(tmp_path):
    path = tmp_path / "casanovo.csv"
    path.write_text(
        "prediction_id,predictions,log_probs,precursor_mz,precursor_charge,group\n"
        "1,GQLYFC[Carbamidomethyl]ELC[Carbamidomethyl]R,-0.12892728609586385,673.30444,2,runB\n",
        encoding="utf-8",
    )

    [candidate] = load_casanovo_predictions(path)

    assert candidate.scan_number is None
    assert candidate.spectrum_id is None
    assert candidate.peptidoform == "GQLYFC[UNIMOD:4]ELC[UNIMOD:4]R"


def _table(tmp_path, header, *rows):
    path = tmp_path / "casanovo.csv"
    path.write_text("\n".join([header, *rows]) + "\n", encoding="utf-8")
    return path


# Casanovo's own peptide score (#157): a mean of amino-acid scores in 4.x, their product in 5.x, minus 1 up to 5.1 when
# the peptide does not fit the precursor m/z


def test_a_precursor_mismatch_gets_no_confidence(tmp_path):
    path = _table(tmp_path, "predictions,peptide_score", "PEPTIDEK,-0.2")

    [candidate] = load_casanovo_predictions(path, casanovo_version="4.2.1")

    # read as a log-probability, -0.2 became confidence 0.82
    assert candidate.confidence is None and candidate.confidence_type is None
    assert candidate.raw_score == -0.2
    assert candidate.raw_score_name == "Casanovo peptide score (mean amino-acid score)"
    assert "does not fit the precursor m/z" in candidate.warnings[0]


def test_a_v5_peptide_score_is_a_sequence_probability(tmp_path):
    path = _table(tmp_path, "predictions,search_engine_score[1]", "PEPTIDEK,0.73")

    [candidate] = load_casanovo_predictions(path, casanovo_version="v5.2.1")

    assert candidate.confidence == 0.73
    assert candidate.confidence_type == "casanovo_sequence_probability"
    assert candidate.raw_score_name == "Casanovo peptide score (amino-acid score product)"
    assert candidate.source_engine_version == "5.2.1"
    assert candidate.warnings == []


def test_a_v4_mean_score_is_not_a_probability(tmp_path):
    path = _table(tmp_path, "predictions,peptide_score", "PEPTIDEK,0.91")

    [candidate] = load_casanovo_predictions(path, casanovo_version="4.3.0")

    assert candidate.confidence is None and candidate.raw_score == 0.91
    assert candidate.warnings == []


def test_without_a_version_the_peptide_score_gets_no_confidence(tmp_path):
    path = _table(tmp_path, "predictions,peptide_score", "PEPTIDEK,0.91")

    [candidate] = load_casanovo_predictions(path)

    assert candidate.confidence is None and candidate.raw_score == 0.91
    assert "Pass casanovo_version" in candidate.warnings[0]


def test_log_probs_are_read_as_log_probabilities_whatever_their_sign(tmp_path):
    # InstaNovo-style noise just above 0 is a probability of 1, as is 0 itself
    path = _table(tmp_path, "predictions,log_probs", "PEPTIDEK,-23.025850929940457", "SAMPLER,1.9e-08", "PEPTIDER,0.0")

    candidates = load_casanovo_predictions(path)

    assert math.isclose(candidates[0].confidence, 1e-10)
    assert [candidate.confidence for candidate in candidates[1:]] == [1.0, 1.0]
    assert {candidate.confidence_type for candidate in candidates} == {"sequence_probability_from_log_probability"}


def test_a_log_probs_column_with_positive_scores_is_refused(tmp_path):
    path = _table(tmp_path, "predictions,log_probs", "PEPTIDEK,-0.2", "SAMPLER,0.87")

    with pytest.raises(ValueError, match="name the column 'peptide_score'"):
        load_casanovo_predictions(path)


def test_two_score_columns_are_refused(tmp_path):
    path = _table(tmp_path, "predictions,log_probs,peptide_score", "PEPTIDEK,-0.2,0.8")

    with pytest.raises(ValueError, match="several score columns"):
        load_casanovo_predictions(path)


def test_casanovo_version_is_checked(tmp_path):
    path = _table(tmp_path, "predictions,peptide_score", "PEPTIDEK,0.8")

    with pytest.raises(ValueError, match="casanovo_version must look like"):
        load_casanovo_predictions(path, casanovo_version="latest")
    with pytest.raises(ValueError, match="applies only to adapter='casanovo'"):
        load_candidates_tool(str(path), "instanovo", casanovo_version="5.2.1")
    [candidate] = load_candidates_tool(str(path), "casanovo", casanovo_version="5.0.0")
    assert candidate["confidence"] == 0.8
