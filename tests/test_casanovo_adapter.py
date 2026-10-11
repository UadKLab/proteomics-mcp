from __future__ import annotations

import math

from proteomics_mcp.adapters.casanovo import load_casanovo_predictions


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
