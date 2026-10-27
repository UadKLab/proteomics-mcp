"""Spectrum loader paths without real instrument files: value conversion, table formats, and the Bruker .d conversion
against a stand-in AlphaTims reader (#166)."""

from __future__ import annotations

import sys

import numpy as np
import pandas as pd
import pytest

from proteomics_mcp.core import spectrum_loader
from proteomics_mcp.core.spectrum_loader import bruker_d_to_feather, load_spectrum_table, to_float_list

EMBEDDED_MGF = (
    "BEGIN IONS\nTITLE=Cmpd 2, +PASEF(922.0056)\nRTINSECONDS=2.74687\nPEPMASS=922.00560 8713\nCHARGE=1+\n"
    "921.91941 71\n922.00624 7624\nEND IONS\n"
)


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (None, []),
        ([1, "2.5"], [1.0, 2.5]),
        (np.array([1.0, 2.0]), [1.0, 2.0]),
        ("[100.5, 200]", [100.5, 200.0]),
        ("100.5 200 3e2", [100.5, 200.0, 300.0]),  # not JSON: every number in the text
        ("  ", []),
        ((1, 2), [1.0, 2.0]),
    ],
)
def test_to_float_list(value, expected):
    assert to_float_list(value) == expected


class _Undecodable(bytes):
    def decode(self, *args, **kwargs):
        raise UnicodeDecodeError("utf-8", b"", 0, 1, "bad")


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (b"HCD", "HCD"),
        (np.float32(1.5), 1.5),
        (np.int64(3), 3.0),
        ([1, 2], "[1, 2]"),
        (None, None),
        (True, True),
    ],
)
def test_normalize_meta_value(value, expected):
    assert spectrum_loader._normalize_meta_value(value) == expected


def test_undecodable_bytes_become_their_text():
    assert spectrum_loader._normalize_meta_value(_Undecodable(b"\xff")) == str(_Undecodable(b"\xff"))


def test_a_spectrum_table_must_be_feather_or_arrow(tmp_path):
    path = tmp_path / "spectra.csv"
    path.write_text("a\n1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Unsupported spectrum table format: .csv"):
        load_spectrum_table(path)


class _TimsTOF:
    """The AlphaTims tables bruker_d_to_feather reads; each test removes or changes what it needs."""

    frames = pd.DataFrame({"Id": [1, 2], "Time": [10.0, 10.1]})
    precursors = pd.DataFrame({"Id": [0, 1, 2, 3], "MonoisotopicMz": [0.0, 500.0, 600.0, 700.0], "Charge": [0, 2, 3, 2],
                               "Parent": [1, 1, 1, 1], "OneOverK0": [np.nan, 0.9, 1.1, 1.2]})
    fragment_frames = pd.DataFrame({"Frame": [2, 2, 2], "Precursor": [1, 1, 2], "IsolationMz": [500.0, 500.1, 600.0],
                                    "CollisionEnergy": [30.0, 31.0, 35.0]})
    mz_values = np.array([100.0, 200.0, 300.0])

    def __init__(self, path):
        self.path = path

    def index_precursors(self):
        # precursor 1: two peaks at m/z 100 (summed) and one at 200; precursor 2: one peak; precursor 3: none
        return np.array([0, 0, 3, 4, 4]), np.array([0, 0, 1, 2]), np.array([5.0, 5.0, 20.0, 30.0])


@pytest.fixture
def d_path(tmp_path):
    path = tmp_path / "run.d"
    path.mkdir()
    return path


@pytest.fixture
def alphatims(monkeypatch):
    alphatims_bruker = pytest.importorskip("alphatims.bruker")

    def use(reader):
        monkeypatch.setattr(alphatims_bruker, "TimsTOF", reader)

    return use


def test_bruker_precursors_become_ms2_spectra(alphatims, d_path, tmp_path):
    alphatims(_TimsTOF)
    output = tmp_path / "run.feather"

    result = bruker_d_to_feather(d_path, output)
    table = pd.read_feather(output)

    # precursor 0 is not a precursor, and precursor 3 has no peaks
    assert result == {"output_path": str(output), "spectra_written": 2, "source_format": "bruker_d"}
    first = table.iloc[0]
    assert first["spectrum_id"] == "run:precursor:1"
    assert list(first["mz_array"]) == [100.0, 200.0] and list(first["intensity_array"]) == [10.0, 20.0]
    assert first["mobility"] == 0.9  # the precursor's own 1/K0
    assert first["rt"] == 10.0  # its parent frame's time
    assert first["isolation_target_mz"] == 500.0 and first["collision_energy"] == 30.0  # the first fragment frame
    assert table["precursor_charge"].tolist() == [2, 3]


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"frames": None}, "did not expose a frames table"),
        ({"precursors": pd.DataFrame()}, "no precursors table"),
        ({"fragment_frames": pd.DataFrame()}, "no fragment_frames table"),
    ],
)
def test_bruker_tables_alphatims_did_not_give_are_reported(alphatims, d_path, tmp_path, change, message):
    alphatims(type("Reader", (_TimsTOF,), change))

    with pytest.raises(RuntimeError, match=message):
        bruker_d_to_feather(d_path, tmp_path / "run.feather")


def test_without_fragment_frames_the_embedded_mgf_is_used(alphatims, d_path, tmp_path):
    alphatims(type("Reader", (_TimsTOF,), {"fragment_frames": pd.DataFrame()}))
    (d_path / "ProteinAnalysisResults.mgf").write_text(EMBEDDED_MGF, encoding="utf-8")
    output = tmp_path / "run.feather"

    result = bruker_d_to_feather(d_path, output)

    assert result["source_format"] == "bruker_d_embedded_mgf"
    table = pd.read_feather(output)
    assert table["source_name"].tolist() == ["run"]
    assert table["metadata"].iloc[0]["bruker_d_fallback"] is True


def test_an_unreadable_d_folder_without_an_embedded_mgf_raises_alphatims_error(alphatims, d_path, tmp_path):
    def unreadable(path):
        raise OSError("not a timsTOF run")

    alphatims(unreadable)

    with pytest.raises(OSError, match="not a timsTOF run"):
        bruker_d_to_feather(d_path, tmp_path / "run.feather")


def test_without_alphatims_the_embedded_mgf_is_used(monkeypatch, d_path, tmp_path):
    monkeypatch.setitem(sys.modules, "alphatims.bruker", None)
    output = tmp_path / "run.feather"

    with pytest.raises(RuntimeError, match="requires AlphaTims"):
        bruker_d_to_feather(d_path, output)

    (d_path / "ProteinAnalysisResults.mgf").write_text(EMBEDDED_MGF, encoding="utf-8")
    result = bruker_d_to_feather(d_path, output)
    assert result["fallback_used"] == "ProteinAnalysisResults.mgf"
    assert pd.read_feather(output)["metadata"].iloc[0]["embedded_mgf_path"] == str(d_path / "ProteinAnalysisResults.mgf")


def test_the_d_path_must_be_an_existing_folder(tmp_path):
    with pytest.raises(FileNotFoundError):
        bruker_d_to_feather(tmp_path / "missing.d", tmp_path / "out.feather")
    (tmp_path / "file.d").write_text("", encoding="utf-8")
    with pytest.raises(ValueError, match="must be a directory"):
        bruker_d_to_feather(tmp_path / "file.d", tmp_path / "out.feather")


class _Sliceable:
    """AlphaTims' raw[:, :, precursor, :, "df"] access, used when no precursor index is given."""

    def __init__(self, result):
        self.result = result

    def __getitem__(self, key):
        if isinstance(self.result, Exception):
            raise self.result
        return self.result


def test_peaks_by_slicing_are_summed_per_mz():
    raw = _Sliceable(pd.DataFrame({"mz_values": [200.0, 100.0, 200.0], "intensity_values": [1.0, 2.0, 3.0]}))

    assert spectrum_loader._extract_alphatims_precursor_peaks(raw, 1) == ([100.0, 200.0], [2.0, 4.0])


@pytest.mark.parametrize("result", [pd.DataFrame(), pd.DataFrame({"other": [1.0]}), None])
def test_slicing_without_peaks_gives_none(result):
    assert spectrum_loader._extract_alphatims_precursor_peaks(_Sliceable(result), 1) == ([], [])


def test_slicing_that_fails_is_reported():
    with pytest.raises(RuntimeError, match="AlphaTims precursor extraction failed"):
        spectrum_loader._extract_alphatims_precursor_peaks(_Sliceable(KeyError("precursor")), 1)


def test_inverse_mobility_needs_a_calibration_and_a_scan_position():
    assert spectrum_loader._one_over_k0(type("Raw", (), {"mobility_values": np.array([])})(), 10.0) is None
    assert spectrum_loader._one_over_k0(type("Raw", (), {})(), 10.0) is None
    assert spectrum_loader._one_over_k0(type("Raw", (), {"mobility_values": np.array([1.6, 0.6])})(), None) is None
