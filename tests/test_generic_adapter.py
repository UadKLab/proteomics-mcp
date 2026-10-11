from __future__ import annotations

import pandas as pd

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
