"""mzml_to_feather keeps the isolation window and ion injection time (#77).

The fixture is a minimal mzML with two MS2 spectra: one whose isolation window target equals the selected-ion m/z,
and one where they differ (a 13C pick), as ThermoRawFileParser writes them. Expected values are read from the mzML's
cvParams with ElementTree, independently of pyOpenMS.
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("pyopenms")

from proteomics_mcp.core.spectrum_loader import mzml_to_feather

FIXTURE = Path(__file__).parent / "fixtures" / "mzml_isolation_injection.mzML"
NS = {"mz": "http://psi.hupo.org/ms/mzml"}


def _cv_params(element) -> dict[str, float]:
    return {param.get("accession"): float(param.get("value")) for param in element.iter("{http://psi.hupo.org/ms/mzml}cvParam")
            if param.get("value") not in (None, "")}


def _reference_from_xml() -> dict[int, dict[str, float]]:
    reference = {}
    for spectrum in ET.parse(FIXTURE).getroot().iterfind(".//mz:spectrum", NS):
        scan = int(spectrum.get("id").rsplit("scan=", 1)[1])
        scan_params = _cv_params(spectrum.find(".//mz:scan", NS))
        window = _cv_params(spectrum.find(".//mz:isolationWindow", NS))
        selected = _cv_params(spectrum.find(".//mz:selectedIon", NS))
        reference[scan] = {
            "precursor_mz": selected["MS:1000744"],
            "isolation_target_mz": window["MS:1000827"],
            "isolation_lower_offset": window["MS:1000828"],
            "isolation_upper_offset": window["MS:1000829"],
            "ion_injection_time": scan_params["MS:1000927"],
            "scan_resolution": scan_params["MS:1000800"],
        }
    return reference


@pytest.fixture(scope="module")
def converted(tmp_path_factory) -> pd.DataFrame:
    output = tmp_path_factory.mktemp("mzml") / "spectra.feather"
    mzml_to_feather(FIXTURE, output)
    return pd.read_feather(output).set_index("scan_number")


def test_isolation_window_and_injection_time_match_the_mzml(converted):
    reference = _reference_from_xml()
    assert set(converted.index) == set(reference) == {101, 102}
    for scan, expected in reference.items():
        row = converted.loc[scan]
        for column, value in expected.items():
            assert float(row[column]) == pytest.approx(value), (scan, column)


def test_isolation_target_differs_from_the_precursor_when_the_instrument_picked_another_peak(converted):
    """Scan 101: target = selected ion. Scan 102: the selected ion is the 13C peak, 0.5 m/z (1/z) above the target."""
    assert converted.loc[101, "isolation_target_mz"] == converted.loc[101, "precursor_mz"]
    assert converted.loc[102, "precursor_mz"] - converted.loc[102, "isolation_target_mz"] == pytest.approx(0.5015, abs=1e-4)


def test_scan_level_metadata_is_kept_in_the_metadata_column(converted):
    metadata = converted.loc[101, "metadata"]
    assert metadata["acquisition_meta"]["MS:1000927"] == pytest.approx(22.5)
