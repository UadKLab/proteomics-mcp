"""Koina request building and response parsing, on hand-written Triton payloads (#166). No network."""

from __future__ import annotations

import json

import pytest

from proteomics_mcp.predictors import koina

SPECTRUM_SETTINGS = {
    "collision_energy": 30.0,
    "fragmentation_type": "HCD",
    "instrument_type": "QE",
    "server_url": "https://koina.example",
    "cache_hit": False,
}


@pytest.mark.parametrize(
    ("data", "shape"),
    [([[1, 2, 3], [4, 5, 6]], [2, 3]), ([1, 2], [2]), ([], [0]), ("x", [1])],
)
def test_infer_triton_shape(data, shape):
    assert koina._infer_triton_shape(data) == shape


@pytest.mark.parametrize(
    ("values", "batch_size", "rows"),
    [
        (5, 2, [[5], [5]]),  # a scalar is repeated for every item
        ([1, 2, 3], 1, [[1, 2, 3]]),  # one item: the flat list is its row
        ([[1, 2], [3, 4]], 2, [[1, 2], [3, 4]]),
        ([1, 2], 2, [[1], [2]]),  # one value per item
        ([[1, 2], 3, [4]], 2, [[1, 2], [3]]),  # mixed: keep lists, wrap scalars, stop at the batch size
    ],
)
def test_ensure_row_lists(values, batch_size, rows):
    assert koina._ensure_row_lists(values, batch_size) == rows


def test_decode_outputs_reshapes_flat_tensors_and_keeps_the_rest_as_sent():
    payload = {
        "outputs": [
            {"name": "mz", "shape": [2, 2], "data": [1.0, 2.0, 3.0, 4.0]},
            {"name": "intensities", "shape": [2, 3], "data": [1.0, 2.0]},  # shape does not match: kept as sent
            {"name": "", "shape": [1], "data": [0]},  # unnamed: dropped
            "not an output",
        ]
    }

    assert koina._decode_outputs(payload) == {"mz": [[1.0, 2.0], [3.0, 4.0]], "intensities": [1.0, 2.0]}


def test_reshape_flat_values_of_a_scalar_shape():
    assert koina._reshape_flat_values([7], []) == 7
    assert koina._reshape_flat_values([], []) is None


def test_annotations_are_rebuilt_from_ion_types_ordinals_and_charges():
    predicted = koina._predicted_spectrum_from_rows(
        peptidoform="PEPTIDEK",
        precursor_charge=3,
        model="Prosit_2020_intensity_HCD",
        max_fragment_charge=None,
        mz_row=[175.119, -1.0, 300.5, 450.2],
        intensity_row=[0.5, 1.0, 1.0, 2.0],
        annotation_row=None,
        ion_type_row=["y", "y", "b", "y"],
        ordinal_row=[1, 2, 3, 4],
        charge_row=[1, 1, 2, 3],
        **SPECTRUM_SETTINGS,
    )

    # the negative m/z is Prosit's "not available" mark; 3+ is above the precursor charge - 1 limit
    assert [fragment.annotation for fragment in predicted.fragments] == ["y1+", "b3++"]
    assert [fragment.intensity for fragment in predicted.fragments] == [0.5, 1.0]  # divided by the strongest kept
    assert predicted.metadata["available_fragment_charges"] == [1, 2]
    assert predicted.predictor_version is None  # no model version reported


def test_missing_charges_default_to_one():
    predicted = koina._predicted_spectrum_from_rows(
        peptidoform="PEPTIDEK", precursor_charge=2, model="m", max_fragment_charge=1, mz_row=[175.1], intensity_row=[0.0],
        annotation_row=None, ion_type_row=["y"], ordinal_row=[1], charge_row=None, **SPECTRUM_SETTINGS,
    )

    assert predicted.fragments[0].annotation == "y1+"
    assert predicted.fragments[0].intensity == 0.0  # an all-zero spectrum stays zero
    assert predicted.metadata["fragment_charge_rule"] == "explicit max_fragment_charge"


def test_a_response_without_annotations_or_their_parts_is_refused():
    with pytest.raises(ValueError, match="did not include fragment annotations"):
        koina._predicted_spectrum_from_rows(
            peptidoform="PEPTIDEK", precursor_charge=2, model="m", max_fragment_charge=None, mz_row=[175.1], intensity_row=[1.0],
            annotation_row=None, ion_type_row=["y"], ordinal_row=None, charge_row=None, **SPECTRUM_SETTINGS,
        )


def test_arrays_of_different_lengths_are_refused():
    with pytest.raises(ValueError, match="equal lengths"):
        koina._predicted_spectrum_from_rows(
            peptidoform="PEPTIDEK", precursor_charge=2, model="m", max_fragment_charge=None, mz_row=[175.1, 200.0], intensity_row=[1.0],
            annotation_row=["y1+1"], ion_type_row=None, ordinal_row=None, charge_row=None, **SPECTRUM_SETTINGS,
        )


def test_request_payload_follows_the_advertised_inputs():
    payload = koina._koina_request_payload(
        [("PEPTIDEK", 2), ("SAMPLER", 3)],
        model_inputs=[
            {"name": "peptide_sequences", "datatype": "BYTES", "shape": [-1, 1]},
            {"name": "precursor_charges", "datatype": "INT32", "shape": [-1]},
        ],
        collision_energy=28,
        fragmentation_type="HCD",
        instrument_type="QE",
    )

    assert payload["inputs"][0] == {"name": "peptide_sequences", "shape": [2, 1], "datatype": "BYTES", "data": [["PEPTIDEK"], ["SAMPLER"]]}
    assert payload["inputs"][1] == {"name": "precursor_charges", "shape": [2], "datatype": "INT32", "data": [2, 3]}


@pytest.mark.parametrize(
    ("model_inputs", "message"),
    [
        ([{"name": "retention_times", "shape": [-1, 1]}], "unsupported input 'retention_times'"),
        ([], "did not advertise any inputs"),
        ([{"name": "peptide_sequences", "shape": [-1, 30]}], "Only scalar per-spectrum Koina inputs"),
    ],
)
def test_request_payload_refuses_inputs_it_cannot_fill(model_inputs, message):
    with pytest.raises(ValueError, match=message):
        koina._koina_request_payload([("PEPTIDEK", 2)], model_inputs=model_inputs, collision_energy=30, fragmentation_type="HCD", instrument_type="QE")


@pytest.mark.parametrize(
    ("peptidoform", "charge", "message"),
    [("A" * 31, 2, "length <= 30 residues; got 31"), ("PEPTIDEK", 7, "precursor charge <= 6; got 7")],
)
def test_prosit_limits_are_checked_before_sending(peptidoform, charge, message):
    with pytest.raises(ValueError, match=message):
        koina._validate_prosit_constraints(peptidoform, charge, "Prosit_2020_intensity_HCD")
    koina._validate_prosit_constraints(peptidoform, charge, "ms2pip_HCD2021")  # other models are not limited here


def test_the_default_cache_is_in_the_temporary_directory(monkeypatch, tmp_path):
    monkeypatch.setattr(koina.tempfile, "gettempdir", lambda: str(tmp_path))

    assert koina._koina_cache_dir(None) == tmp_path / "proteomics_mcp_koina_cache"
    assert (tmp_path / "proteomics_mcp_koina_cache").is_dir()


class _Response:
    def __init__(self, payload) -> None:
        self._payload = payload

    def read(self) -> bytes:
        return json.dumps(self._payload).encode("utf-8")

    def __enter__(self):
        return self

    def __exit__(self, *_) -> None:
        return None


def _serve(monkeypatch, metadata, outputs):
    def fake_urlopen(req, timeout):
        return _Response(metadata if req.data is None else {"model_version": "3", "outputs": outputs})

    monkeypatch.setattr(koina, "_urlopen", fake_urlopen)
    koina._describe_koina_model_cached.cache_clear()


INPUTS = {"inputs": [{"name": "peptide_sequences", "datatype": "BYTES", "shape": [-1, 1]}, {"name": "precursor_charges", "datatype": "INT32", "shape": [-1, 1]}]}


def test_a_batch_answered_with_ion_types_ordinals_and_charges(monkeypatch, tmp_path):
    _serve(monkeypatch, INPUTS, [
        {"name": "mz", "shape": [2, 2], "data": [175.119, 300.2, 147.1, 250.3]},
        {"name": "intensities", "shape": [2, 2], "data": [1.0, 0.5, 0.25, 1.0]},
        {"name": "ion_types", "shape": [2, 2], "data": ["y", "b", "y", "b"]},
        {"name": "fragment_numbers", "shape": [2, 2], "data": [1, 3, 1, 2]},
        {"name": "charges", "shape": [2, 2], "data": [1, 1, 1, 2]},
    ])

    first, second = koina.predict_fragments_koina_batch([("PEPTIDEK", 2), ("SAMPLER", 3)], model="Prosit_2020_intensity_HCD", cache_dir=str(tmp_path))

    assert [fragment.annotation for fragment in first.fragments] == ["y1+", "b3+"]
    assert [fragment.annotation for fragment in second.fragments] == ["y1+", "b2++"]
    assert first.predictor_version == "Prosit_2020_intensity_HCD@3"


def test_a_batch_answer_for_fewer_items_is_refused(monkeypatch, tmp_path):
    _serve(monkeypatch, INPUTS, [
        {"name": "mz", "shape": [1, 2], "data": [175.119, 300.2]},
        {"name": "intensities", "shape": [1, 2], "data": [1.0, 0.5]},
        {"name": "annotation", "shape": [1, 2], "data": ["y1+1", "b3+1"]},
    ])

    with pytest.raises(ValueError, match="batch size does not match"):
        koina.predict_fragments_koina_batch([("PEPTIDEK", 2), ("SAMPLER", 3)], model="Prosit_2020_intensity_HCD", cache_dir=str(tmp_path))


def test_model_metadata_that_is_not_an_object_is_refused(monkeypatch):
    _serve(monkeypatch, ["not", "an", "object"], [])

    with pytest.raises(ValueError, match="invalid metadata for model 'Prosit_2020_intensity_HCD'"):
        koina.describe_koina_model(model="Prosit_2020_intensity_HCD")
