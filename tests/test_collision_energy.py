"""Koina predictions use the spectrum's own normalized collision energy, or one given per call (#128)."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

from proteomics_mcp.core.acquisition import normalized_collision_energy
from proteomics_mcp.predictors import koina
from proteomics_mcp.tools.workflows import build_prediction_batch_tool, predict_fragments_tool, score_spectrum_tool

THERMO_FILTER = "FTMS + c NSI cv=-50.00 d Full ms2 498.7812@hcd28.00 [120.0000-1039.1456]"


@pytest.mark.parametrize(
    ("metadata", "expected"),
    [
        ({"filter_string": THERMO_FILTER}, 28.0),
        ({"filter_string": "ITMS + c NSI d Full ms2 445.12@cid35.00 [110.00-900.00]"}, 35.0),
        ({"filter_string": "FTMS + p NSI d Full ms2 500.0", "collision_energy": 27}, 27.0),  # Thermo, no @ token
        ({"collision_energy": 27, "instrument_model": "Orbitrap Exploris 480"}, 27.0),
        ({"collision_energy": 42.0}, None),  # eV from a timsTOF: not an NCE
        ({}, None),
    ],
)
def test_normalized_collision_energy(metadata, expected):
    assert normalized_collision_energy({"metadata": metadata})[0] == expected


def test_the_pride_response_states_its_nce_in_the_filter_string():
    record = json.loads((Path(__file__).parent / "fixtures" / "proxi_pride_PXD044934_scan39484.json").read_text(encoding="utf-8"))
    assert normalized_collision_energy({"metadata": {"attributes": record["attributes"]}}) == (28.0, "scan filter @hcd28.00")


class _Response:
    def __init__(self, payload):
        self._payload = payload

    def read(self):
        return json.dumps(self._payload).encode("utf-8")

    def __enter__(self):
        return self

    def __exit__(self, *args):
        return None


@pytest.fixture
def koina_requests(monkeypatch, tmp_path):
    """A fake Koina server; the list holds the collision energies of every inference request, one list per request."""
    sent: list[list[float]] = []

    def fake_urlopen(req, timeout):
        if req.data is None:
            return _Response({"name": "Prosit_2020_intensity_HCD", "inputs": [
                {"name": "peptide_sequences", "datatype": "BYTES", "shape": [-1, -1]},
                {"name": "precursor_charges", "datatype": "INT32", "shape": [-1, 1]},
                {"name": "collision_energies", "datatype": "FP32", "shape": [-1, 1]}]})
        inputs = {item["name"]: item["data"] for item in json.loads(req.data.decode("utf-8"))["inputs"]}
        energies = [row[0] for row in inputs["collision_energies"]]
        sent.append(energies)
        n = len(energies)
        return _Response({"model_version": "1", "outputs": [
            {"name": "mz", "shape": [n, 2], "data": [175.119, 288.203] * n},
            {"name": "intensities", "shape": [n, 2], "data": [1.0, 0.5] * n},
            {"name": "annotation", "shape": [n, 2], "data": ["y1+1", "y2+1"] * n}]})

    monkeypatch.setattr(koina, "_urlopen", fake_urlopen)
    monkeypatch.setattr(koina, "_koina_cache_dir", lambda cache_dir: tmp_path)
    koina._describe_koina_model_cached.cache_clear()
    return sent


CANDIDATE = {"candidate_id": "c", "source_engine": "test", "peptidoform": "PEPTIDEK", "precursor_charge": 2}
KOINA = {"predictor_backend": "koina", "predictor_model": "Prosit_2020_intensity_HCD"}


def _observed(**metadata):
    return {"mz": [175.119, 288.203], "intensity": [10.0, 5.0], "precursor_mz": 464.7, "precursor_charge": 2, "metadata": metadata}


def test_score_spectrum_predicts_at_the_spectrum_nce(koina_requests):
    result = score_spectrum_tool(CANDIDATE, _observed(filter_string=THERMO_FILTER), **KOINA)

    assert koina_requests == [[28.0]]
    assert result["predictor_metadata"]["collision_energy"] == 28.0
    assert result["predictor_metadata"]["collision_energy_source"] == "spectrum: scan filter @hcd28.00"


def test_an_explicit_collision_energy_wins(koina_requests):
    result = score_spectrum_tool(CANDIDATE, _observed(filter_string=THERMO_FILTER), collision_energy=33.0, **KOINA)

    assert koina_requests == [[33.0]]
    assert result["predictor_metadata"]["collision_energy_source"] == "collision_energy argument"


def test_an_ev_collision_energy_falls_back_to_the_config(koina_requests):
    result = score_spectrum_tool(CANDIDATE, _observed(collision_energy=42.0), **KOINA)

    assert koina_requests == [[30.0]]
    assert result["predictor_metadata"]["collision_energy_source"].startswith("config koina_collision_energy (collision energy 42")


def test_predict_fragments_takes_a_collision_energy(koina_requests):
    predicted = predict_fragments_tool("PEPTIDEK", 2, "Prosit_2020_intensity_HCD", "koina", collision_energy=25.0)

    assert koina_requests == [[25.0]]
    assert predicted["metadata"]["collision_energy"] == 25.0


def test_prediction_batch_groups_requests_by_each_spectrum_nce(koina_requests, tmp_path):
    spectra = tmp_path / "spectra.feather"
    pd.DataFrame([
        {"source_name": "run", "scan_number": scan, "precursor_mz": 464.7, "precursor_charge": 2, "mz_array": [175.1], "intensity_array": [1.0],
         "filter_string": f"FTMS + c NSI d Full ms2 464.7@hcd{nce:.2f} [100.0-1000.0]"}
        for scan, nce in ((1, 28.0), (2, 32.0), (3, 28.0))
    ]).to_feather(spectra)
    candidates = tmp_path / "preds.csv"
    candidates.write_text("experiment_name,scan_number,precursor_charge,preds,log_probs\nrun,1,2,PEPTIDEK,-0.1\nrun,2,2,PEPTIDEK,-0.1\nrun,3,2,SAMPLER,-0.2\n",
                          encoding="utf-8")

    manifest = build_prediction_batch_tool(path=str(candidates), output_prefix=str(tmp_path / "batch"), spectra_path=str(spectra),
                                           predictor_backend="koina", predictor_model="Prosit_2020_intensity_HCD")

    # PEPTIDEK at 28 and at 32 are two predictions; one Koina request per collision energy
    assert manifest["unique_prediction_requests"] == 3
    assert sorted(sorted(set(request)) for request in koina_requests) == [[28.0], [32.0]]
    predictions = pd.read_csv(manifest["output_files"]["predictions_csv"])
    assert sorted(predictions["collision_energy"]) == [28.0, 28.0, 32.0]
    assert predictions["collision_energy_source"].str.startswith("spectrum: scan filter").all()
