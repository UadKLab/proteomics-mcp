"""Every InstaNovo output format of Zenodo 10.5281/zenodo.20756892, read by the InstaNovo adapter."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from proteomics_mcp.adapters.instanovo import _split_tokens, load_instanovo_predictions
from proteomics_mcp.core.spectrum_loader import load_spectrum_from_table

FIXTURES = Path(__file__).parent / "fixtures" / "instanovo_formats"
ALL_FILES = sorted(FIXTURES.glob("*.csv"))
MGF_FILES = [path for path in ALL_FILES if path.stem.endswith("_mgf")]


def _raw(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def _by_index(path: Path) -> dict[int, object]:
    return {candidate.spectrum_index: candidate for candidate in load_instanovo_predictions(path)}


def test_fixtures_cover_every_format():
    assert len(ALL_FILES) == 17


@pytest.mark.parametrize("path", ALL_FILES, ids=[path.stem for path in ALL_FILES])
def test_every_format_loads(path):
    candidates = load_instanovo_predictions(path)

    assert candidates
    assert all(0.0 <= candidate.confidence <= 1.0 for candidate in candidates if candidate.confidence is not None)
    for candidate in candidates:
        if candidate.token_log_probabilities is not None and candidate.token_confidence_min_position is not None:
            assert len(candidate.token_log_probabilities) == len(_split_tokens(candidate.peptidoform, candidate.source_columns.get(
                f"{candidate.source_columns['prediction_column']}_tokenised") or None))


@pytest.mark.parametrize("path", MGF_FILES, ids=[path.stem for path in MGF_FILES])
def test_mgf_scan_numbers_are_positions(path):
    # InstaNovo writes the 0-based position of the spectrum in the MGF file as scan_number.
    raw = _raw(path)
    candidates = load_instanovo_predictions(path)

    assert all(candidate.scan_number is None for candidate in candidates)
    assert {candidate.spectrum_index for candidate in candidates} <= {int(value) for value in raw["scan_number"]}


def test_mzml_scan_numbers_are_kept():
    candidates = load_instanovo_predictions(FIXTURES / "v1.1.3_transformer_greedy_mzml.csv")

    assert [candidate.scan_number for candidate in candidates] == [1501, 1727, 1738]
    assert all(candidate.spectrum_index is None for candidate in candidates)


@pytest.mark.parametrize("decoding", ["greedy", "beam5"])
def test_v110_transformer_columns_are_read(decoding):
    # 1.1.0 writes predictions with log_probabilities and token_log_probabilities, which was not recognised.
    candidates = load_instanovo_predictions(FIXTURES / f"v1.1.0_transformer_{decoding}_mgf.csv")

    assert candidates
    assert {candidate.raw_score_name for candidate in candidates} == {"log_probabilities"}
    assert all(candidate.token_log_probabilities for candidate in candidates)


def test_refined_instanovoplus_uses_final_prediction():
    path = FIXTURES / "v1.1.3_instanovoplus_refined_mgf.csv"
    raw = _raw(path).set_index("scan_number")
    candidates = load_instanovo_predictions(path)

    differing = [c for c in candidates if raw.loc[str(c.spectrum_index), "final_prediction"] != raw.loc[str(c.spectrum_index), "diffusion_predictions"]]
    assert differing
    for candidate in candidates:
        row = raw.loc[str(candidate.spectrum_index)]
        expected = "final_prediction" if row["final_prediction"] else "transformer_predictions"
        assert candidate.source_columns["prediction_column"] == expected


def test_empty_refinement_falls_back_to_the_transformer_prediction():
    path = FIXTURES / "v1.1.0_instanovoplus_alpha_refined_mgf.csv"
    raw = _raw(path).set_index("scan_number")
    candidates = _by_index(path)

    empty = [int(index) for index, row in raw.iterrows() if not row["diffusion_predictions"]]
    assert empty
    for index in empty:
        assert candidates[index].source_columns["prediction_column"] == "transformer_predictions"
        assert candidates[index].raw_score_name == "transformer_log_probabilities"


@pytest.mark.parametrize("path", [p for p in ALL_FILES if "instanovoplus" in p.stem], ids=lambda p: p.stem)
def test_log_probability_just_above_zero_is_certainty(path):
    candidates = load_instanovo_predictions(path)

    positive = [c for c in candidates if c.raw_score is not None and c.raw_score > 0]
    for candidate in positive:
        assert candidate.confidence == pytest.approx(1.0)


def test_empty_beam_prediction_is_skipped():
    path = FIXTURES / "v1.0.0_transformer_beam5_mgf.csv"
    raw = _raw(path)

    candidates = load_instanovo_predictions(path)

    assert raw.loc[0, "preds"] == "" and raw.loc[0, "log_probs"] == "-inf"
    assert len(candidates) == int((raw["preds"] != "").sum())


@pytest.mark.parametrize("version", ["v1.0.0", "v1.1.0", "v1.1.3"])
def test_beam_token_log_probabilities_follow_the_sequence(version):
    # Legacy beam files store one extra trailing value, and 1.1.x the rest in reverse order. Compare with greedy
    # decoding of the same spectrum to the same peptide.
    greedy = _by_index(FIXTURES / f"{version}_transformer_greedy_mgf.csv")
    beam = _by_index(FIXTURES / f"{version}_transformer_beam5_mgf.csv")

    same = [index for index in greedy.keys() & beam.keys() if greedy[index].peptidoform == beam[index].peptidoform]
    assert same
    for index in same:
        assert np.allclose(beam[index].token_log_probabilities, greedy[index].token_log_probabilities, atol=1e-2)
        assert beam[index].token_confidence_min_position == greedy[index].token_confidence_min_position


@pytest.mark.parametrize("path", [p for p in ALL_FILES if "beam5" in p.stem], ids=lambda p: p.stem)
def test_beam_files_report_the_next_beam(path):
    candidates = load_instanovo_predictions(path)

    assert any(candidate.next_beam is not None for candidate in candidates)


def test_spectrum_index_looks_up_by_position_within_the_run(tmp_path):
    table = tmp_path / "spectra.feather"
    pd.DataFrame(
        {
            "source_file": ["runA.mgf", "runA.mgf", "runB.mgf"],
            "spectrum_id": ["runA.1501.1501.2", "runA.1727.1727.4", "runB.10.10.2"],
            "scan_number": [1501, 1727, 10],
            "precursor_mz": [500.0, 600.0, 700.0],
            "precursor_charge": [2, 4, 2],
            "mz_array": [[1.0], [2.0], [3.0]],
            "intensity_array": [[1.0], [1.0], [1.0]],
        }
    ).to_feather(table)

    assert load_spectrum_from_table(table, spectrum_index=1, experiment_name="runA").precursor_mz == 600.0
    assert load_spectrum_from_table(table, spectrum_index=0, experiment_name="runB").precursor_mz == 700.0
    with pytest.raises(LookupError, match="run name"):
        load_spectrum_from_table(table, spectrum_index=0, experiment_name="typo")
    with pytest.raises(LookupError, match="outside"):
        load_spectrum_from_table(table, spectrum_index=2, experiment_name="runA")


def test_spectrum_index_matches_mgf_order_after_conversion(tmp_path):
    # spectrum_index relies on mgf_to_feather keeping spectra in file order. The scans are deliberately not sorted,
    # so a table ordered by scan number would return the wrong spectra.
    from proteomics_mcp.core.spectrum_loader import mgf_to_feather

    run = "SF_test_run"
    spectra = [(1727, 450.25, 4), (1501, 893.91, 2), (1738, 610.80, 3)]
    mgf = tmp_path / f"{run}.mgf"
    mgf.write_text(
        "".join(
            f"BEGIN IONS\nTITLE={run}.{scan}.{scan}.{charge}\nPEPMASS={mz}\nCHARGE={charge}+\nSCANS={scan}\n100.0 10.0\n{200.0 + i} 20.0\nEND IONS\n"
            for i, (scan, mz, charge) in enumerate(spectra)
        ),
        encoding="utf-8",
    )
    table = tmp_path / "spectra.feather"
    mgf_to_feather(mgf, table)

    for position, (scan, mz, charge) in enumerate(spectra):
        observed = load_spectrum_from_table(table, spectrum_index=position, experiment_name=run)
        assert (observed.precursor_mz, observed.precursor_charge, observed.scan_number) == (mz, charge, scan)

    # The same through the InstaNovo adapter: MGF-based output writes positions 0, 1, 2 as scan_number.
    predictions = tmp_path / "predictions.csv"
    pd.DataFrame(
        {
            "scan_number": [0, 1, 2],
            "precursor_mz": [mz for _, mz, _ in spectra],
            "precursor_charge": [charge for _, _, charge in spectra],
            "experiment_name": [run] * 3,
            "spectrum_id": [f"{run}:{position}" for position in range(3)],
            "preds": ["PEPTIDEK", "LESLIEK", "AAAAK"],
            "log_probs": [-0.1, -0.2, -0.3],
        }
    ).to_csv(predictions, index=False)
    for candidate in load_instanovo_predictions(predictions):
        observed = load_spectrum_from_table(
            table, spectrum_index=candidate.spectrum_index, scan_number=candidate.scan_number, experiment_name=candidate.experiment_name
        )
        assert observed.precursor_mz == candidate.precursor_mz
