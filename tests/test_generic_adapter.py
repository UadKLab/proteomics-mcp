from __future__ import annotations

import pandas as pd
import pytest

from proteomics_mcp.adapters.generic import load_generic_candidates


def test_load_generic_candidates_supports_xlsx(tmp_path):
    path = tmp_path / "generic.xlsx"
    pd.DataFrame(
        [
            {
                "sequence": "PEPTIDE",
                "scan": 42,
                "score": 0.9,
                "charge": 2,
                "precursor_mz": 500.2,
            }
        ]
    ).to_excel(path, index=False)

    [candidate] = load_generic_candidates(
        path,
        peptide_column="sequence",
        scan_column="scan",
        score_column="score",
        score_type="confidence_0_1",
        charge_column="charge",
        precursor_mz_column="precursor_mz",
        source_engine="mapped",
    )

    assert candidate.source_engine == "mapped"
    assert candidate.peptidoform == "PEPTIDE"
    assert candidate.scan_number == 42
    assert candidate.confidence == 0.9


# What each declared score_type makes of the score (#166)


def _generic(tmp_path, score, **mapping):
    path = tmp_path / "generic.csv"
    pd.DataFrame([{"Sequence": "PEPTIDEK", "Score": score, "Next": 0.25}, {"Sequence": " ", "Score": 0.5, "Next": 0.1}]).to_csv(path, index=False)
    return load_generic_candidates(path, peptide_column="Sequence", score_column="Score", **mapping)


@pytest.mark.parametrize(
    ("score_type", "score", "confidence"),
    [
        ("confidence_0_1", 0.8, 0.8),
        ("confidence_0_1", 1.7, 1.0),  # clipped to [0, 1]
        ("posterior_error_probability", 0.05, 0.95),
        ("posterior_error_probability", 1.5, 0.0),
        ("evalue", 1e-5, None),
        ("score_unscaled", 42.0, None),
        ("unknown", 0.8, None),
    ],
)
def test_score_type_sets_the_confidence(tmp_path, score_type, score, confidence):
    [candidate] = _generic(tmp_path, score, score_type=score_type)  # the blank peptide row is skipped

    assert candidate.raw_score == score
    assert candidate.confidence == (pytest.approx(confidence) if confidence is not None else None)
    assert candidate.confidence_type == score_type


def test_an_unsupported_score_type_is_refused(tmp_path):
    with pytest.raises(ValueError, match="Unsupported score_type: probability"):
        _generic(tmp_path, 0.8, score_type="probability")


def test_next_beam_margin_is_derived_from_the_next_score(tmp_path):
    [candidate] = _generic(tmp_path, 0.75, next_score_column="Next")

    assert candidate.next_beam == 0.25
    assert candidate.next_beam_margin == pytest.approx(0.5)
