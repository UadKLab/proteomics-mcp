from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

import pytest

from proteomics_mcp.core import spectrum_loader
from proteomics_mcp.core.spectrum_loader import load_spectrum_from_table, retrieve_spectrum_from_usi


def test_load_spectrum_from_table_prefers_experiment_name_with_duplicate_scans(tmp_path):
    table = pd.DataFrame(
        [
            {
                "source_file": "/tmp/run_a.mzML",
                "spectrum_id": "controllerType=0 controllerNumber=1 scan=100",
                "scan_number": 100,
                "precursor_mz": 500.0,
                "precursor_charge": 2,
                "mz_array": [100.0, 200.0],
                "intensity_array": [10.0, 20.0],
            },
            {
                "source_file": "/tmp/run_b.mzML",
                "spectrum_id": "controllerType=0 controllerNumber=1 scan=100",
                "scan_number": 100,
                "precursor_mz": 700.0,
                "precursor_charge": 3,
                "mz_array": [300.0, 400.0],
                "intensity_array": [30.0, 40.0],
            },
        ]
    )
    path = tmp_path / "spectra.feather"
    table.to_feather(path)

    observed = load_spectrum_from_table(path, scan_number=100, experiment_name="run_b")

    assert observed.precursor_mz == 700.0
    assert observed.precursor_charge == 3
    assert observed.mz == [300.0, 400.0]


def test_mgf_to_feather_converts_basic_mgf(tmp_path):
    mgf_path = tmp_path / "demo.mgf"
    mgf_path.write_text(
        "BEGIN IONS\n"
        "TITLE=demo.42.42.2.dta\n"
        "CHARGE=2+\n"
        "PEPMASS=500.2 1234\n"
        "RTINSECONDS=12.5\n"
        "100.0 10.0\n"
        "200.0 20.0\n"
        "END IONS\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "demo.feather"

    result = spectrum_loader.mgf_to_feather(mgf_path, output_path)
    observed = load_spectrum_from_table(output_path, scan_number=42)

    assert result["source_format"] == "mgf"
    assert observed.scan_number == 42
    assert observed.precursor_mz == 500.2
    assert observed.precursor_charge == 2
    assert observed.mz == [100.0, 200.0]
    assert observed.intensity == [10.0, 20.0]


def test_bruker_d_to_feather_falls_back_to_embedded_mgf(tmp_path):
    d_path = tmp_path / "sample.d"
    d_path.mkdir()
    mgf_path = d_path / "ProteinAnalysisResults.mgf"
    mgf_path.write_text(
        "BEGIN IONS\n"
        "TITLE=Cmpd 2, +PASEF(922.0056), 41.4eV, 0.0min, 1/K0=1.1894, #21-25\n"
        "RTINSECONDS=2.74687\n"
        "RAWSCANS=19,21-23,25\n"
        "PEPMASS=922.00560 8713\n"
        "CHARGE=1+\n"
        "ION_MOBILITY=922.00560 8713 1.18943\n"
        "921.91941 71\n"
        "922.00624 7624\n"
        "END IONS\n",
        encoding="utf-8",
    )
    output_path = tmp_path / "bruker.feather"

    result = spectrum_loader.bruker_d_to_feather(d_path, output_path)
    observed = load_spectrum_from_table(output_path, spectrum_id="Cmpd 2, +PASEF(922.0056), 41.4eV, 0.0min, 1/K0=1.1894, #21-25")

    assert result["source_format"] == "bruker_d_embedded_mgf"
    assert result["fallback_used"] == "ProteinAnalysisResults.mgf"
    assert observed.precursor_mz == pytest.approx(922.00560)
    assert observed.precursor_charge == 1
    assert observed.metadata["bruker_d_fallback"] is True


def _recorded_proxi_response() -> dict:
    """The PRIDE PROXI response for the HeLa TLFPLFMK spectrum (scan 39484), recorded on 2026-10-14."""
    record = json.loads((Path(__file__).parent / "fixtures" / "proxi_pride_PXD044934_scan39484.json").read_text(encoding="utf-8"))
    return {"usi": record["usi"], "status": "READABLE", "attributes": record["attributes"],
            "m/z array": record["mzs"], "intensity array": record["intensities"]}


def test_usi_spectrum_carries_the_precursor_from_proxi_attributes(monkeypatch):
    """PRIDE sends the precursor as MS:1000744 selected ion m/z, not as "precursor m/z" (#41)."""
    from pyteomics import usi as pyteomics_usi

    response = _recorded_proxi_response()
    monkeypatch.setattr(pyteomics_usi, "proxi", lambda *_, **__: response)

    spectrum = retrieve_spectrum_from_usi(response["usi"])

    assert spectrum.precursor_mz == 498.781219482422
    assert spectrum.precursor_charge == 2
    assert spectrum.scan_number == 39484
    assert len(spectrum.mz) == 77


def test_usi_precursor_falls_back_to_the_isolation_target(monkeypatch):
    from pyteomics import usi as pyteomics_usi

    response = _recorded_proxi_response()
    response["attributes"] = [a for a in response["attributes"] if a.get("accession") != "MS:1000744"]
    response["attributes"].append({"accession": "MS:1000827", "name": "isolation window target m/z", "value": 498.5})
    response["attributes"] = [a for a in response["attributes"] if a.get("value") != 498.781219482422]
    monkeypatch.setattr(pyteomics_usi, "proxi", lambda *_, **__: response)

    assert retrieve_spectrum_from_usi(response["usi"]).precursor_mz == 498.5


def test_usi_spectrum_scores_with_precursor_signal_fraction(monkeypatch):
    """With the precursor m/z, USI-based scores include precursor_signal_fraction again."""
    from pyteomics import usi as pyteomics_usi

    from proteomics_mcp.tools.workflows import score_spectrum_tool

    payload_path = Path(__file__).resolve().parents[1] / "examples" / "22_hela_instanovo_top5" / "rank_4_scan_39484_TLFPLFMK_payload.json"
    if not payload_path.is_file():
        pytest.skip("examples/ is only available in a source checkout")
    payload = json.loads(payload_path.read_text(encoding="utf-8"))
    response = _recorded_proxi_response()
    monkeypatch.setattr(pyteomics_usi, "proxi", lambda *_, **__: response)
    observed = retrieve_spectrum_from_usi(response["usi"]).model_dump()

    result = score_spectrum_tool(payload["metrics"]["candidate"], observed, payload["predicted"])

    assert result["precursor_signal_fraction"] == pytest.approx(payload["metrics"]["precursor_signal_fraction"], abs=1e-3)
