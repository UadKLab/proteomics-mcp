"""Modification notation is normalised in the generic adapter and in accuracy targets (#124)."""

from __future__ import annotations

import json

import pandas as pd
import pytest

from proteomics_mcp.core.dataset import candidate_models_to_frame, summarize_prediction_accuracy
from proteomics_mcp.models.schemas import CandidatePSM
from proteomics_mcp.tools.workflows import load_candidates_tool


def _generic(tmp_path, peptides):
    table = tmp_path / "psms.csv"
    pd.DataFrame({"peptide": peptides, "charge": [2] * len(peptides)}).to_csv(table, index=False)
    mapping = tmp_path / "mapping.json"
    mapping.write_text(json.dumps({"peptide_column": "peptide", "charge_column": "charge"}), encoding="utf-8")
    return load_candidates_tool(str(table), adapter="generic", column_mapping_path=str(mapping))


def test_generic_adapter_writes_every_notation_of_a_modification_the_same_way(tmp_path):
    candidates = _generic(tmp_path, ["ALPM[+15.9949]HIR", "ALPM(+15.99)HIR", "ALPM[Oxidation]HIR", "ALPM[UNIMOD:35]HIR",
                                     "PEPS[+79.9663]TIDE", "PEPC[+57.0215]TIDE", "[+42.0106]PEPTIDE"])

    assert [c["peptidoform"] for c in candidates] == ["ALPM[UNIMOD:35]HIR"] * 4 + ["PEPS[UNIMOD:21]TIDE", "PEPC[UNIMOD:4]TIDE", "[UNIMOD:1]-PEPTIDE"]
    assert all(c["warnings"] == [] for c in candidates)
    assert candidates[1]["source_columns"]["peptide"] == "ALPM(+15.99)HIR"  # as written, kept


def test_generic_adapter_warns_about_a_modification_it_cannot_resolve(tmp_path):
    unknown, unmatched_mass = _generic(tmp_path, ["PEPM[Foo]K", "PEPK[+123.4567]"])

    assert unknown["peptidoform"] == "PEPM[Foo]K"
    assert len(unknown["warnings"]) == 1 and "[Foo]" in unknown["warnings"][0]
    # a mass that matches no UNIMOD entry is still a usable mass shift
    assert unmatched_mass["peptidoform"] == "PEPK[+123.4567]" and unmatched_mass["warnings"] == []


def _accuracy(pairs):
    candidates = [
        CandidatePSM(candidate_id=f"c{i}", source_engine="instanovo", peptidoform=prediction, precursor_charge=2, confidence=0.9,
                     source_columns={"target": target})
        for i, (prediction, target) in enumerate(pairs)
    ]
    return summarize_prediction_accuracy(candidate_models_to_frame(candidates), target_field="source_columns.target")


def test_a_target_in_another_notation_counts_as_exactly_correct():
    result = _accuracy([("ALPM[UNIMOD:35]HIR", "ALPM[+15.9949]HIR"), ("GQLYFC[UNIMOD:4]ELC[UNIMOD:4]R", "GQLYFC[Carbamidomethyl]ELC[Carbamidomethyl]R"),
                        ("PEPTIDE", "PEPTIDE")])

    overview = result["overview"]
    assert overview["exact_accuracy"] == 1.0
    assert overview["targets_renotated"] == 2
    assert overview["targets_with_unresolved_modifications"] == 0
    by_mod = {row["label"]: row for row in result["by_target_modification_state"]}
    assert by_mod["modified"]["count"] == 2 and by_mod["modified"]["accuracy_exact"] == 1.0


def test_targets_that_keep_an_unknown_modification_are_counted():
    result = _accuracy([("PEPM[UNIMOD:35]K", "PEPM[Foo]K")])

    assert result["overview"]["exact_accuracy"] == 0.0
    assert result["overview"]["targets_with_unresolved_modifications"] == 1


def test_a_wrong_modification_is_still_wrong():
    assert _accuracy([("ALPM[UNIMOD:35]HIR", "ALPM[+0.9840]HIR")])["overview"]["exact_accuracy"] == 0.0
