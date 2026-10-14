"""Joining a scored table onto candidate rows: candidate_id first, then run + scan + peptidoform (#39)."""

from __future__ import annotations

import pandas as pd

from proteomics_mcp.core.dataset import join_scored_enrichment
from proteomics_mcp.tools.workflows import inspect_psm_table_tool, summarize_numeric_field_tool

FRAGPIPE_COLUMNS = ["Spectrum", "Peptide", "Modified Peptide", "Charge", "Calibrated Observed M/Z", "Hyperscore", "Probability",
                    "Retention", "Protein"]


def _fragpipe_table(path):
    rows = [
        ["runA.00010.00010.2", "PEPTIDE", "", 2, 400.2, 20.0, 0.99, 600.0, "P1"],
        ["runA.00011.00011.2", "PEPTIDEK", "", 2, 464.2, 18.0, 0.98, 610.0, "P1"],
        ["runA.00012.00012.2", "SAMPLER", "", 2, 401.2, 15.0, 0.95, 620.0, "P2"],
    ]
    pd.DataFrame(rows, columns=FRAGPIPE_COLUMNS).to_csv(path, sep="\t", index=False)


def _scored_table(path, candidate_ids):
    """Scored rows as score_candidates writes them: the ids come from the file the candidates were scored from."""
    pd.DataFrame([
        {"candidate_id": candidate_ids[0], "experiment_name": "runA", "scan_number": 10, "peptidoform": "PEPTIDE", "spectral_angle": 0.81},
        {"candidate_id": candidate_ids[1], "experiment_name": "runA", "scan_number": 11, "peptidoform": "PEPTIDEK", "spectral_angle": 0.62},
    ]).to_csv(path, index=False)


def test_scored_table_joins_by_candidate_id(tmp_path):
    candidates = tmp_path / "psm.tsv"
    _fragpipe_table(candidates)
    scored = tmp_path / "scored.csv"
    _scored_table(scored, ["psm:0", "psm:1"])
    result = inspect_psm_table_tool(path=str(candidates), adapter="fragpipe", scored_path=str(scored))
    assert {key: result["scored_join"][key] for key in ("scored_rows", "matched_by_candidate_id", "matched_by_run_scan_peptidoform",
                                                      "unmatched_scored_rows")} == {
        "scored_rows": 2, "matched_by_candidate_id": 2, "matched_by_run_scan_peptidoform": 0, "unmatched_scored_rows": 0}
    assert "warnings" not in result


def test_renamed_candidate_table_joins_by_run_scan_and_peptidoform(tmp_path):
    """The Zenodo case: the candidate table was renamed after scoring, so no candidate_id matches."""
    candidates = tmp_path / "fragpipe_psm.tsv"
    _fragpipe_table(candidates)
    scored = tmp_path / "scored.csv"
    _scored_table(scored, ["psm:0", "psm:1"])
    result = summarize_numeric_field_tool(path=str(candidates), adapter="fragpipe", scored_path=str(scored), field="spectral_angle")
    assert result["scored_join"]["matched_by_candidate_id"] == 0
    assert result["scored_join"]["matched_by_run_scan_peptidoform"] == 2
    assert result["scored_join"]["unmatched_scored_rows"] == 0
    assert "warnings" not in result
    assert (result["count"], result["min"], result["max"]) == (2, 0.62, 0.81)  # the scores reached the candidate rows


def test_unmatched_scored_table_is_reported(tmp_path):
    candidates = tmp_path / "fragpipe_psm.tsv"
    _fragpipe_table(candidates)
    scored = tmp_path / "scored.csv"
    pd.DataFrame([{"candidate_id": "other:0", "experiment_name": "runB", "scan_number": 99, "peptidoform": "OTHERPEP",
                   "spectral_angle": 0.5}]).to_csv(scored, index=False)
    result = inspect_psm_table_tool(path=str(candidates), adapter="fragpipe", scored_path=str(scored))
    assert result["scored_join"]["unmatched_scored_rows"] == 1
    assert len(result["warnings"]) == 1 and "matched none of its 1 rows" in result["warnings"][0]


def test_ambiguous_fallback_keys_are_not_guessed():
    frame = pd.DataFrame([
        {"candidate_id": "a:0", "run_key": "runA", "scan_number": 10, "peptidoform": "PEPTIDE"},
        {"candidate_id": "a:1", "run_key": "runA", "scan_number": 10, "peptidoform": "PEPTIDE"},  # same key twice
        {"candidate_id": "a:2", "run_key": "runA", "scan_number": 11, "peptidoform": "PEPTIDEK"},
    ])
    scored = pd.DataFrame([
        {"candidate_id": "b:0", "run_key": "runA", "scan_number": 10, "peptidoform": "PEPTIDE", "spectral_angle": 0.9},
        {"candidate_id": "b:1", "run_key": "runA", "scan_number": 11, "peptidoform": "PEPTIDEK", "spectral_angle": 0.7},
    ])
    merged, summary = join_scored_enrichment(frame, scored)
    assert summary["matched_by_run_scan_peptidoform"] == 1 and summary["unmatched_scored_rows"] == 1
    angles = merged.set_index("candidate_id")["spectral_angle"]
    assert angles.isna().to_dict() == {"a:0": True, "a:1": True, "a:2": False}
    assert angles["a:2"] == 0.7


def test_candidate_id_match_takes_precedence_over_the_fallback():
    frame = pd.DataFrame([
        {"candidate_id": "a:0", "run_key": "runA", "scan_number": 10, "peptidoform": "PEPTIDE"},
        {"candidate_id": "a:1", "run_key": "runA", "scan_number": 11, "peptidoform": "PEPTIDEK"},
    ])
    scored = pd.DataFrame([
        {"candidate_id": "a:1", "run_key": "runA", "scan_number": 11, "peptidoform": "PEPTIDEK", "spectral_angle": 0.7},
        # same run/scan/peptidoform as a:1, but a:1 is already taken by its own id
        {"candidate_id": "x:9", "run_key": "runA", "scan_number": 11, "peptidoform": "PEPTIDEK", "spectral_angle": 0.1},
    ])
    merged, summary = join_scored_enrichment(frame, scored)
    assert summary["matched_by_candidate_id"] == 1 and summary["matched_by_run_scan_peptidoform"] == 0
    assert merged.set_index("candidate_id").loc["a:1", "spectral_angle"] == 0.7
