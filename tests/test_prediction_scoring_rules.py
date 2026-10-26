"""One set of prediction-side rules for every path to the scorer, and a record of what was scored (#137)."""

from __future__ import annotations

import json
import math
from pathlib import Path

import pandas as pd
import pytest

import proteomics_mcp.tools.workflows as workflows
from proteomics_mcp.core.scoring import prepare_predicted_for_scoring, score_candidate_spectrum
from proteomics_mcp.models.schemas import CandidatePSM, ObservedSpectrum, PredictedFragment, PredictedSpectrum, model_dump
from proteomics_mcp.predictors import koina
from proteomics_mcp.tools.workflows import (
    build_prediction_batch_tool,
    score_candidates_tool,
    score_from_prediction_batch_tool,
    score_spectrum_tool,
)

FIXTURE = Path(__file__).parent / "fixtures" / "koina_prosit2020_FFLVLR_z2.json"
SIMILARITY = ("spectral_angle", "cosine_similarity", "pearson_correlation")


class _Response:
    def __init__(self, payload: dict) -> None:
        self._payload = payload

    def read(self) -> bytes:
        return json.dumps(self._payload).encode("utf-8")

    def __enter__(self) -> "_Response":
        return self

    def __exit__(self, *args) -> None:
        return None


@pytest.fixture
def recorded_koina(monkeypatch, tmp_path):
    """Serve the recorded Koina response for FFLVLR 2+ (Prosit 2020, collision energy 30) instead of the network."""
    recorded = json.loads(FIXTURE.read_text(encoding="utf-8"))
    cache = tmp_path / "koina_cache"
    cache.mkdir()
    monkeypatch.setattr(koina, "_koina_cache_dir", lambda cache_dir: cache)

    def fake_urlopen(req, timeout):
        return _Response(recorded["model_metadata"] if req.data is None else recorded["infer_response"])

    monkeypatch.setattr(koina, "_urlopen", fake_urlopen)
    koina._describe_koina_model_cached.cache_clear()
    return recorded


def _raw_koina_prediction(recorded: dict) -> dict:
    """The recorded response as a caller would pass it in: every fragment charge, intensities as Koina returned them."""
    outputs = {output["name"]: output["data"] for output in recorded["infer_response"]["outputs"]}
    fragments = []
    for mz, intensity, annotation in zip(outputs["mz"], outputs["intensities"], outputs["annotation"]):
        if mz <= 0:
            continue
        ion, rest = annotation[0], annotation[1:]
        ordinal, charge = rest.split("+")
        fragments.append({"mz": mz, "intensity": intensity, "annotation": annotation, "ion_type": ion, "ordinal": int(ordinal), "charge": int(charge)})
    return {"peptidoform": "FFLVLR", "precursor_charge": 2, "predictor_name": "koina", "model": "Prosit_2020_intensity_HCD", "fragments": fragments}


def _fflvlr_observed(recorded: dict) -> dict:
    # Peaks on six of the predicted 1+ ions and on y5++, plus noise; enough for every metric to be defined.
    peaks = {175.119: 20.0, 288.2030: 8.0, 295.1441: 30.0, 387.2714: 25.0, 500.3555: 60.0, 647.4239: 100.0, 324.2156: 3.0, 410.0: 5.0, 512.3: 4.0}
    return {"mz": list(peaks), "intensity": list(peaks.values()), "precursor_mz": 397.7470, "precursor_charge": 2}


CANDIDATE = {"candidate_id": "c1", "source_engine": "test", "peptidoform": "FFLVLR", "precursor_charge": 2}
KOINA = {"predictor_backend": "koina", "predictor_model": "Prosit_2020_intensity_HCD"}


def test_supplied_koina_prediction_scores_like_the_internal_one(recorded_koina):
    observed = _fflvlr_observed(recorded_koina)
    internal = score_spectrum_tool(CANDIDATE, observed, **KOINA)
    supplied = score_spectrum_tool(CANDIDATE, observed, _raw_koina_prediction(recorded_koina), **KOINA)

    for metric in (*SIMILARITY, "matched_ion_count", "total_predicted_ions", "explained_intensity_fraction"):
        assert supplied[metric] == pytest.approx(internal[metric]), metric
    assert internal["total_predicted_ions"] == 7
    assert [ion["annotation"] for ion in supplied["scored_ions"]] == [ion["annotation"] for ion in internal["scored_ions"]]
    above_limit = [ion for ion in supplied["excluded_predicted_ions"] if ion["reason"] == "fragment_charge_above_limit"]
    assert above_limit and all(ion["charge"] == 2 for ion in above_limit)
    assert "y5+2" in [ion["annotation"] for ion in above_limit]


def test_supplied_prediction_can_be_scored_as_given(recorded_koina):
    observed = _fflvlr_observed(recorded_koina)
    result = score_spectrum_tool(CANDIDATE, observed, _raw_koina_prediction(recorded_koina), score_predictions_as_given=True)

    assert result["total_predicted_ions"] == 8
    assert "y5+2" in [ion["annotation"] for ion in result["scored_ions"]]
    assert result["scoring_settings"]["score_predictions_as_given"] is True
    assert result["scoring_settings"]["fragment_charge_rule"] == "none (score_predictions_as_given)"


def test_internal_ms2pip_prediction_scores_the_same_when_passed_back(ms2pip_runtime):
    observed = {"mz": [175.119, 288.203, 401.287, 500.36, 613.44, 726.53], "intensity": [5.0, 20.0, 10.0, 40.0, 30.0, 8.0],
                "precursor_mz": 640.85, "precursor_charge": 3}
    candidate = {**CANDIDATE, "peptidoform": "LGEYGFQNALIVR", "precursor_charge": 3}
    internal = score_spectrum_tool(candidate, observed)
    predicted = workflows._predict_fragment_spectrum(workflows.load_config(), "LGEYGFQNALIVR", 3)
    supplied = score_spectrum_tool(candidate, observed, model_dump(predicted))

    for metric in (*SIMILARITY, "matched_ion_count", "total_predicted_ions", "xcorr", "hyperscore"):
        assert supplied[metric] == internal[metric], metric
    assert supplied["scored_ions"] == internal["scored_ions"]


def _spectrum(fragments: list[PredictedFragment], charge: int = 2) -> PredictedSpectrum:
    return PredictedSpectrum(peptidoform="PEPTIDEK", precursor_charge=charge, predictor_name="test", fragments=fragments)


def test_each_left_out_ion_records_its_reason():
    fragments = [
        PredictedFragment(mz=200.0, intensity=2.0, annotation="y1+", charge=1),
        PredictedFragment(mz=300.0, intensity=4.0, annotation="y2+", charge=1),
        PredictedFragment(mz=150.5, intensity=3.0, annotation="y2++", charge=2),
        PredictedFragment(mz=-1.0, intensity=0.5, annotation="y9+", charge=1),
        PredictedFragment(mz=400.0, intensity=math.nan, annotation="y3+", charge=1),
        PredictedFragment(mz=500.0, intensity=0.0, annotation="y4+", charge=1),
        PredictedFragment(mz=600.0, intensity=0.001, annotation="y5+", charge=1),
    ]

    prepared = prepare_predicted_for_scoring(_spectrum(fragments), precursor_charge=2, min_predicted_intensity=0.001)

    assert prepared.indices == [0, 1]
    assert [fragment.intensity for fragment in prepared.fragments] == [0.5, 1.0]
    assert {item["annotation"]: item["reason"] for item in prepared.exclusions} == {
        "y2++": "fragment_charge_above_limit",
        "y9+": "invalid_mz",
        "y3+": "non_finite_intensity",
        "y4+": "non_positive_intensity",
        "y5+": "below_min_predicted_intensity",
    }
    assert (prepared.max_fragment_charge, prepared.fragment_charge_rule) == (1, "charge <= precursor_charge - 1")


def test_explicit_max_fragment_charge_wins_over_the_default_rule():
    fragments = [PredictedFragment(mz=200.0, intensity=1.0, annotation="y1+", charge=1), PredictedFragment(mz=150.5, intensity=1.0, annotation="y2++", charge=2)]

    prepared = prepare_predicted_for_scoring(_spectrum(fragments), precursor_charge=2, max_fragment_charge=2)

    assert prepared.indices == [0, 1]
    assert prepared.fragment_charge_rule == "explicit max_fragment_charge"


def test_legacy_v1_scores_a_supplied_prediction_as_given():
    fragments = [PredictedFragment(mz=200.0, intensity=1.0, annotation="y1+", charge=1), PredictedFragment(mz=150.5, intensity=1.0, annotation="y2++", charge=2)]
    candidate = CandidatePSM(candidate_id="c1", source_engine="test", peptidoform="PEPTIDEK", precursor_charge=2)
    observed = ObservedSpectrum(mz=[150.5, 200.0], intensity=[10.0, 20.0])

    result = score_candidate_spectrum(candidate=candidate, observed=observed, predicted=_spectrum(fragments), metric_profile="legacy_v1", use_observed_annotation=False)

    assert result.total_predicted_ions == 2
    assert result.excluded_predicted_ions == []


def test_result_records_what_was_scored(recorded_koina):
    result = score_spectrum_tool(CANDIDATE, _fflvlr_observed(recorded_koina), **KOINA)

    assert result["metric_n"] == {"spectral_angle": 7, "cosine_similarity": 7, "pearson_correlation": 7, "matched_ions": 6, "total_predicted_ions": 7}
    settings = result["scoring_settings"]
    assert settings["collision_energy"] == 30.0
    assert settings["model_version"] == "1"
    assert settings["tolerance"] == 0.05 and settings["tolerance_unit"] == "Da"
    assert settings["metric_profile"] == "v2" and settings["min_predicted_intensity"] == 0.001
    unmatched = [ion for ion in result["scored_ions"] if ion["observed_index"] is None]
    assert [ion["observed_intensity"] for ion in unmatched] == [0.0]
    assert sum(1 for ion in result["excluded_predicted_ions"] if ion["reason"] == "below_min_predicted_intensity") > 0


@pytest.fixture
def candidates(tmp_path):
    path = tmp_path / "preds.csv"
    path.write_text("experiment_name,scan_number,precursor_charge,preds,log_probs\nrun,10,2,PEPTIDE,-0.1\nrun,11,2,PEPTIDEK,-0.2\n", encoding="utf-8")
    return path


def _fake_spectra(monkeypatch):
    spectrum = workflows.ObservedSpectrum(mz=[200.0, 300.0, 400.0], intensity=[10.0, 20.0, 5.0], precursor_mz=500.0, precursor_charge=2)
    monkeypatch.setattr(workflows, "load_spectrum_from_table", lambda *_, **__: spectrum)


def test_scored_batch_records_settings_and_the_scored_ions(candidates, tmp_path, monkeypatch):
    _fake_spectra(monkeypatch)
    manifest = build_prediction_batch_tool(path=str(candidates), output_prefix=str(tmp_path / "batch"), predictor_backend="mock")
    scored = score_from_prediction_batch_tool(manifest_path=manifest["output_files"]["manifest_json"], spectra_path=str(tmp_path / "x.feather"),
                                              output_prefix=str(tmp_path / "scored"))

    scoring_manifest = json.loads(Path(scored["output_files"]["manifest_json"]).read_text(encoding="utf-8"))
    assert scoring_manifest["scoring_settings"]["fragment_charge_rule"] == "charge <= precursor_charge - 1"
    rows = pd.read_csv(scored["output_files"]["scored_csv"])
    assert {"tolerance", "min_predicted_intensity", "similarity_n", "excluded_predicted_ion_count"} <= set(rows.columns)
    details = [json.loads(line) for line in Path(scored["output_files"]["scored_ions_jsonl"]).read_text(encoding="utf-8").splitlines()]
    assert [detail["candidate_id"] for detail in details] == rows["candidate_id"].tolist()
    assert all(detail["scored_ions"] for detail in details)


def test_score_candidates_writes_the_scored_ions_next_to_the_table(candidates, tmp_path, monkeypatch):
    _fake_spectra(monkeypatch)
    summary = score_candidates_tool(candidates_path=str(candidates), spectra_path=str(tmp_path / "x.feather"),
                                    output_path=str(tmp_path / "scores.csv"), predictor_backend="mock")

    assert summary["scored_ions_path"] == str(tmp_path / "scores_scored_ions.jsonl")
    assert len(Path(summary["scored_ions_path"]).read_text(encoding="utf-8").splitlines()) == 2
    assert summary["scoring_settings"]["metric_profile"] == "v2"


def test_batch_key_depends_on_the_koina_settings():
    common = {"peptidoform": "PEPTIDEK", "precursor_charge": 2, "predictor_backend": "koina", "model": "Prosit_2020_intensity_HCD", "max_fragment_charge": None}
    at_27 = workflows._prediction_batch_key(**common, koina_settings={"collision_energy": 27.0})
    at_30 = workflows._prediction_batch_key(**common, koina_settings={"collision_energy": 30.0})
    assert at_27 != at_30
    assert workflows._prediction_batch_key(**common) == workflows._prediction_batch_key(**common, koina_settings=None)
