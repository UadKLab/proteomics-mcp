from __future__ import annotations

import io
import json

from proteomics_mcp.predictors import koina
from proteomics_mcp.tools import workflows


class _FakeResponse:
    def __init__(self, payload: dict[str, object]) -> None:
        self._payload = payload

    def read(self) -> bytes:
        return json.dumps(self._payload).encode("utf-8")

    def __enter__(self) -> "_FakeResponse":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        return None


def test_koina_predictor_shapes_triton_request_and_normalizes_response(monkeypatch, tmp_path):
    captured: dict[str, object] = {}

    def fake_urlopen(req, timeout):
        captured["url"] = req.full_url
        captured["timeout"] = timeout
        if req.data is None:
            return _FakeResponse(
                {
                    "name": "Prosit_2020_intensity_HCD",
                    "inputs": [
                        {"name": "peptide_sequences", "datatype": "BYTES", "shape": [-1, -1]},
                        {"name": "precursor_charges", "datatype": "INT32", "shape": [-1, 1]},
                        {"name": "collision_energies", "datatype": "FP32", "shape": [-1, 1]},
                    ],
                }
            )
        captured["payload"] = json.loads(req.data.decode("utf-8"))
        return _FakeResponse(
            {
                "outputs": [
                    {"name": "mz", "shape": [1, 2], "data": [100.0, 200.0]},
                    {"name": "intensities", "shape": [1, 2], "data": [0.2, 1.0]},
                    {"name": "annotation", "shape": [1, 2], "data": ["b1+", "y2+2"]},
                ]
            }
        )

    monkeypatch.setattr(koina, "_urlopen", fake_urlopen)
    koina._describe_koina_model_cached.cache_clear()

    predicted = koina.predict_fragments_koina(
        "PEPM[UNIMOD:35]TIDE",
        3,
        model="Prosit_2020_intensity_HCD",
        collision_energy=27.0,
        fragmentation_type="HCD",
        instrument_type="QE",
        cache_dir=str(tmp_path),
    )

    assert captured["url"] == "https://koina.wilhelmlab.org/v2/models/Prosit_2020_intensity_HCD/infer"
    assert captured["timeout"] == 30.0
    payload = captured["payload"]
    assert isinstance(payload, dict)
    assert payload["inputs"] == [
        {"name": "peptide_sequences", "shape": [1, 1], "datatype": "BYTES", "data": [["PEPM[UNIMOD:35]TIDE"]]},
        {"name": "precursor_charges", "shape": [1, 1], "datatype": "INT32", "data": [[3]]},
        {"name": "collision_energies", "shape": [1, 1], "datatype": "FP32", "data": [[27.0]]},
    ]
    assert predicted.predictor_name == "koina"
    assert predicted.model == "Prosit_2020_intensity_HCD"
    assert [fragment.annotation for fragment in predicted.fragments] == ["b1+", "y2+2"]
    assert [fragment.charge for fragment in predicted.fragments] == [1, 2]
    assert max(fragment.intensity for fragment in predicted.fragments) == 1.0


def test_koina_predictor_uses_cache(monkeypatch, tmp_path):
    calls = {"count": 0}

    def fake_urlopen(req, timeout):
        calls["count"] += 1
        if req.data is None:
            return _FakeResponse(
                {
                    "name": "Prosit_2020_intensity_HCD",
                    "inputs": [
                        {"name": "peptide_sequences", "datatype": "BYTES", "shape": [-1, -1]},
                        {"name": "precursor_charges", "datatype": "INT32", "shape": [-1, 1]},
                        {"name": "collision_energies", "datatype": "FP32", "shape": [-1, 1]},
                    ],
                }
            )
        return _FakeResponse(
            {
                "outputs": [
                    {"name": "mz", "shape": [1, 1], "data": [150.0]},
                    {"name": "intensities", "shape": [1, 1], "data": [1.0]},
                    {"name": "annotation", "shape": [1, 1], "data": ["b1+"]},
                ]
            }
        )

    monkeypatch.setattr(koina, "_urlopen", fake_urlopen)
    koina._describe_koina_model_cached.cache_clear()

    first = koina.predict_fragments_koina("PEPTIDE", 2, cache_dir=str(tmp_path))
    second = koina.predict_fragments_koina("PEPTIDE", 2, cache_dir=str(tmp_path))

    assert calls["count"] == 2
    assert first.fragments[0].mz == second.fragments[0].mz
    assert second.metadata["cache_hit"] is True


def test_koina_prediction_rejects_missing_fragment_outputs(monkeypatch, tmp_path):
    def fake_urlopen(req, timeout):
        if req.data is None:
            return _FakeResponse(
                {
                    "inputs": [
                        {"name": "peptide_sequences", "datatype": "BYTES", "shape": [-1, -1]},
                        {"name": "precursor_charges", "datatype": "INT32", "shape": [-1, 1]},
                    ]
                }
            )
        return _FakeResponse({"outputs": [{"name": "mz", "shape": [1, 1], "data": [150.0]}]})

    monkeypatch.setattr(koina, "_urlopen", fake_urlopen)
    koina._describe_koina_model_cached.cache_clear()

    try:
        koina.predict_fragments_koina("PEPTIDE", 2, cache_dir=str(tmp_path))
    except ValueError as exc:
        assert "m/z and intensity" in str(exc)
    else:  # pragma: no cover - defensive assertion for a required response contract
        raise AssertionError("Expected an invalid Koina response to be rejected.")


def test_predict_fragments_tool_preserves_model_from_config(tmp_path):
    config_path = tmp_path / "predictor.toml"
    config_path.write_text('predictor_backend = "mock"\npredictor_model = "configured-mock"\n', encoding="utf-8")

    result = workflows.predict_fragments_tool("PEPTIDE", 2, config_path=str(config_path))

    assert result["predictor_name"] == "mock"


def test_describe_koina_model_returns_live_schema_shape(monkeypatch):
    def fake_urlopen(req, timeout):
        assert req.data is None
        return _FakeResponse(
            {
                "name": "Prosit_2020_intensity_HCD",
                "inputs": [{"name": "peptide_sequences", "datatype": "BYTES", "shape": [-1, -1]}],
                "outputs": [{"name": "intensities", "datatype": "FP32", "shape": [-1, 174]}],
            }
        )

    monkeypatch.setattr(koina, "_urlopen", fake_urlopen)
    koina._describe_koina_model_cached.cache_clear()
    result = workflows.describe_predictor_model_tool("Prosit_2020_intensity_HCD")

    assert result["backend"] == "koina"
    assert result["supported_by_proteomics_mcp"] is True
    assert result["metadata"]["outputs"][0]["name"] == "intensities"


def test_koina_batch_prediction_integrates_with_prediction_batch_workflow(monkeypatch, tmp_path):
    candidate_path = tmp_path / "preds.csv"
    candidate_path.write_text(
        "experiment_name,scan_number,preds,log_probs,precursor_charge,predictions_log_probability_beam_1\n"
        "runA,1,PEPTIDE,-0.1,2,-0.2\n"
        "runB,2,TESTIDE,-0.2,3,-0.3\n",
        encoding="utf-8",
    )

    def fake_predict_batch(items, **kwargs):
        assert kwargs["model"] == "Prosit_2020_intensity_HCD"
        return [
            koina._predicted_spectrum_from_rows(
                peptidoform=peptidoform,
                precursor_charge=charge,
                model=kwargs["model"],
                max_fragment_charge=kwargs["max_fragment_charge"],
                mz_row=[100.0, 200.0],
                intensity_row=[0.1, 1.0],
                annotation_row=["b1+", "y2+"],
                ion_type_row=None,
                ordinal_row=None,
                charge_row=None,
                collision_energy=kwargs["collision_energy"],
                fragmentation_type=kwargs["fragmentation_type"],
                instrument_type=kwargs["instrument_type"],
                server_url=kwargs["server_url"],
                cache_hit=False,
            )
            for peptidoform, charge in items
        ]

    monkeypatch.setattr(workflows, "predict_fragments_koina_batch", fake_predict_batch)

    result = workflows.build_prediction_batch_tool(
        path=str(candidate_path),
        output_prefix=str(tmp_path / "batch"),
        adapter="instanovo",
        predictor_backend="koina",
    )

    assert result["unique_prediction_requests"] == 2
    predictions = workflows.pd.read_csv(result["output_files"]["predictions_csv"])
    assert set(predictions["predictor_name"]) == {"koina"}
    assert set(predictions["resolved_model"]) == {"Prosit_2020_intensity_HCD"}


def test_koina_batch_failure_is_written_as_row_failures(monkeypatch, tmp_path):
    candidate_path = tmp_path / "preds.csv"
    candidate_path.write_text(
        "experiment_name,scan_number,preds,log_probs,precursor_charge\nrunA,1,PEPTIDE,-0.1,2\n",
        encoding="utf-8",
    )

    def fail_predict_batch(*args, **kwargs):
        raise TimeoutError("Koina timed out")

    monkeypatch.setattr(workflows, "predict_fragments_koina_batch", fail_predict_batch)
    result = workflows.build_prediction_batch_tool(
        path=str(candidate_path),
        output_prefix=str(tmp_path / "failed_batch"),
        adapter="instanovo",
        predictor_backend="koina",
    )

    assert result["predictions_succeeded"] == 0
    assert result["predictions_failed"] == 1
    failures = workflows.pd.read_csv(result["output_files"]["failures_csv"])
    assert "Koina prediction failed" in failures.loc[0, "error"]


def test_koina_prediction_records_the_model_version_that_answered(monkeypatch, tmp_path):
    """predictor_version is "<model>@<Triton model_version>" from the inference response (#22)."""

    def fake_urlopen(req, timeout):
        if req.data is None:
            return _FakeResponse({"name": "Prosit_2020_intensity_HCD", "inputs": [
                {"name": "peptide_sequences", "datatype": "BYTES", "shape": [-1, -1]},
                {"name": "precursor_charges", "datatype": "INT32", "shape": [-1, 1]},
                {"name": "collision_energies", "datatype": "FP32", "shape": [-1, 1]}]})
        return _FakeResponse({
            "model_name": "Prosit_2020_intensity_HCD",
            "model_version": "1",
            "outputs": [
                {"name": "mz", "shape": [1, 2], "data": [100.0, 200.0]},
                {"name": "intensities", "shape": [1, 2], "data": [0.2, 1.0]},
                {"name": "annotation", "shape": [1, 2], "data": ["b1+1", "y2+1"]},
            ],
        })

    monkeypatch.setattr(koina, "_urlopen", fake_urlopen)
    koina._describe_koina_model_cached.cache_clear()

    predicted = koina.predict_fragments_koina("PEPTIDE", 2, model="Prosit_2020_intensity_HCD", cache_dir=str(tmp_path))
    cached = koina.predict_fragments_koina("PEPTIDE", 2, model="Prosit_2020_intensity_HCD", cache_dir=str(tmp_path))

    assert predicted.predictor_version == "Prosit_2020_intensity_HCD@1"
    assert cached.metadata["cache_hit"] is True and cached.predictor_version == "Prosit_2020_intensity_HCD@1"


def test_koina_cache_entries_without_a_version_are_predicted_again(monkeypatch, tmp_path):
    """Predictions cached by releases that did not record the model version are refreshed."""
    calls = {"infer": 0}

    def fake_urlopen(req, timeout):
        if req.data is None:
            return _FakeResponse({"name": "Prosit_2020_intensity_HCD", "inputs": [
                {"name": "peptide_sequences", "datatype": "BYTES", "shape": [-1, -1]},
                {"name": "precursor_charges", "datatype": "INT32", "shape": [-1, 1]},
                {"name": "collision_energies", "datatype": "FP32", "shape": [-1, 1]}]})
        calls["infer"] += 1
        return _FakeResponse({"model_version": "1", "outputs": [
            {"name": "mz", "shape": [1, 1], "data": [100.0]},
            {"name": "intensities", "shape": [1, 1], "data": [1.0]},
            {"name": "annotation", "shape": [1, 1], "data": ["b1+1"]}]})

    monkeypatch.setattr(koina, "_urlopen", fake_urlopen)
    koina._describe_koina_model_cached.cache_clear()
    first = koina.predict_fragments_koina("PEPTIDE", 2, model="Prosit_2020_intensity_HCD", cache_dir=str(tmp_path))
    [entry] = list(tmp_path.glob("*.json"))
    old = json.loads(entry.read_text(encoding="utf-8"))
    old["predictor_version"] = None  # as written by 0.1.9: no version, and no model_version in the metadata
    del old["metadata"]["model_version"]
    entry.write_text(json.dumps(old), encoding="utf-8")

    refreshed = koina.predict_fragments_koina("PEPTIDE", 2, model="Prosit_2020_intensity_HCD", cache_dir=str(tmp_path))

    assert first.predictor_version == refreshed.predictor_version == "Prosit_2020_intensity_HCD@1"
    assert calls["infer"] == 2 and refreshed.metadata["cache_hit"] is False


def _koina_returning(mz: list[float], intensities: list[float], annotations: list[str]):
    def fake_urlopen(req, timeout):
        if req.data is None:
            return _FakeResponse(
                {
                    "name": "Prosit_2020_intensity_HCD",
                    "inputs": [
                        {"name": "peptide_sequences", "datatype": "BYTES", "shape": [-1, -1]},
                        {"name": "precursor_charges", "datatype": "INT32", "shape": [-1, 1]},
                        {"name": "collision_energies", "datatype": "FP32", "shape": [-1, 1]},
                    ],
                }
            )
        n = len(mz)
        return _FakeResponse(
            {
                "model_version": "1",
                "outputs": [
                    {"name": "mz", "shape": [1, n], "data": mz},
                    {"name": "intensities", "shape": [1, n], "data": intensities},
                    {"name": "annotation", "shape": [1, n], "data": annotations},
                ],
            }
        )

    return fake_urlopen


def test_koina_drops_fragments_without_predicted_intensity(monkeypatch, tmp_path):
    # Prosit predicts an intensity for every annotated ion; 0 or less means the ion is not expected (#83).
    monkeypatch.setattr(
        koina,
        "_urlopen",
        _koina_returning(
            [100.0, 200.0, 300.0, 400.0, -1.0, 500.0],
            [0.5, 0.0, 1.0, -0.01, -1.0, float("nan")],
            ["b1+", "b2+", "y1+", "y2+", "y3+", "y4+"],
        ),
    )
    koina._describe_koina_model_cached.cache_clear()

    predicted = koina.predict_fragments_koina("PEPTIDE", 2, model="Prosit_2020_intensity_HCD", cache_dir=str(tmp_path))

    assert [fragment.annotation for fragment in predicted.fragments] == ["b1+", "y1+"]
    assert all(fragment.intensity > 0 for fragment in predicted.fragments)
    assert predicted.metadata["dropped_zero_intensity_fragments"] == 3


def test_koina_cache_entries_with_zero_intensity_fragments_are_filtered(monkeypatch, tmp_path):
    def no_network(req, timeout):
        raise AssertionError("cache hit expected")

    key = koina._cache_key(
        peptidoform="PEPTIDE",
        precursor_charge=2,
        model="Prosit_2020_intensity_HCD",
        collision_energy=30.0,
        fragmentation_type="HCD",
        instrument_type="QE",
        max_fragment_charge=None,
    )
    # An entry as written by <= 0.1.10, which kept zero-intensity fragments.
    koina._cache_path(str(tmp_path), key).write_text(
        json.dumps(
            {
                "peptidoform": "PEPTIDE",
                "precursor_charge": 2,
                "predictor_name": "koina",
                "model": "Prosit_2020_intensity_HCD",
                "metadata": {"model_version": "1"},
                "fragments": [
                    {"mz": 100.0, "intensity": 1.0, "annotation": "b1+", "ion_type": "b", "ordinal": 1, "charge": 1},
                    {"mz": 200.0, "intensity": 0.0, "annotation": "b2+", "ion_type": "b", "ordinal": 2, "charge": 1},
                ],
            }
        ),
        encoding="utf-8",
    )
    monkeypatch.setattr(koina, "_urlopen", no_network)
    monkeypatch.setattr(koina, "_validate_prosit_constraints", lambda *args: None)

    predicted = koina.predict_fragments_koina("PEPTIDE", 2, model="Prosit_2020_intensity_HCD", cache_dir=str(tmp_path))

    assert predicted.metadata["cache_hit"] is True
    assert [fragment.annotation for fragment in predicted.fragments] == ["b1+"]
    assert predicted.metadata["dropped_zero_intensity_fragments"] == 1
