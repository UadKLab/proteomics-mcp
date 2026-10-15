from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from threading import Event, Thread
from time import monotonic, sleep

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


def _start_local_proxi_server(payload: dict, response_gate: Event | None = None):
    request_started = Event()
    body = json.dumps(payload).encode("utf-8")

    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            request_started.set()
            if response_gate is not None:
                response_gate.wait(timeout=5)
            self.send_response(200)
            self.send_header("Content-Type", "application/json")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            try:
                self.wfile.write(body)
            except OSError:
                pass

        def log_message(self, *args):
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    server.daemon_threads = True
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    return server, thread, request_started


def _route_proxi_backends_to_local_servers(monkeypatch, servers):
    from pyteomics import usi as pyteomics_usi

    original_backends = pyteomics_usi._proxies

    def make_local_backend(backend_type, server):
        url_template = f"http://127.0.0.1:{server.server_port}/?usi={{usi}}"

        class LocalBackend(backend_type):
            def __init__(self, **kwargs):
                super().__init__(**kwargs)
                self.url_template = url_template

        return LocalBackend

    monkeypatch.setattr(
        pyteomics_usi,
        "_proxies",
        {
            name: make_local_backend(original_backends[name], server)
            for name, server in servers.items()
        },
    )


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


def test_usi_spectrum_parses_a_delayed_proxi_response_from_the_selected_backend(monkeypatch):
    response_gate = Event()
    payload = {
        "status": "READABLE",
        "attributes": [
            {"accession": "MS:1000744", "name": "selected ion m/z", "value": "500.25"},
            {"name": "charge state", "value": "2"},
        ],
        "mzs": [100.0, 200.0],
        "intensities": [10.0, 20.0],
    }
    server, thread, request_started = _start_local_proxi_server(payload, response_gate)
    _route_proxi_backends_to_local_servers(monkeypatch, {"pride": server})

    def release_delayed_response() -> None:
        if request_started.wait(timeout=1):
            sleep(0.1)
            response_gate.set()

    response_thread = Thread(target=release_delayed_response, daemon=True)
    response_thread.start()

    try:
        from proteomics_mcp.tools.workflows import retrieve_usi_spectrum_tool

        result = retrieve_usi_spectrum_tool("mzspec:PRIDE:sample:scan:42", timeout_seconds=1)
    finally:
        response_gate.set()
        server.shutdown()
        server.server_close()
        thread.join(timeout=1)
        response_thread.join(timeout=1)

    assert result["mz"] == [100.0, 200.0]
    assert result["intensity"] == [10.0, 20.0]
    assert result["precursor_mz"] == 500.25
    assert result["precursor_charge"] == 2
    assert result["scan_number"] == 42


def test_usi_spectrum_single_backend_request_honors_timeout(monkeypatch):
    response_gate = Event()
    server, thread, request_started = _start_local_proxi_server(
        {"attributes": [], "mzs": [100.0], "intensities": [10.0]}, response_gate
    )
    _route_proxi_backends_to_local_servers(monkeypatch, {"pride": server})

    try:
        from proteomics_mcp.tools.workflows import retrieve_usi_spectrum_tool

        with ThreadPoolExecutor(max_workers=1) as executor:
            future = executor.submit(
                retrieve_usi_spectrum_tool,
                "mzspec:PRIDE:sample:scan:42",
                "pride",
                0.05,
            )
            try:
                assert request_started.wait(timeout=1)
                started_at = monotonic()
                with pytest.raises(TimeoutError):
                    future.result(timeout=0.5)
                assert future.done()
                assert monotonic() - started_at < 0.5
            finally:
                response_gate.set()
    finally:
        response_gate.set()
        server.shutdown()
        server.server_close()
        thread.join(timeout=1)


def test_usi_aggregator_honors_timeout_for_its_backend_requests(monkeypatch):
    response_gate = Event()
    server, thread, request_started = _start_local_proxi_server(
        {"attributes": [], "mzs": [100.0], "intensities": [10.0]}, response_gate
    )
    _route_proxi_backends_to_local_servers(monkeypatch, {"pride": server})

    try:
        from proteomics_mcp.tools.workflows import retrieve_usi_spectrum_tool

        with ThreadPoolExecutor(max_workers=1) as executor:
            future = executor.submit(
                retrieve_usi_spectrum_tool,
                "mzspec:PRIDE:sample:scan:42",
                "aggregator",
                0.05,
            )
            try:
                assert request_started.wait(timeout=1)
                started_at = monotonic()
                with pytest.raises(ValueError, match="No valid responses found"):
                    future.result(timeout=0.5)
                assert future.done()
                assert monotonic() - started_at < 0.5
            finally:
                response_gate.set()
    finally:
        response_gate.set()
        server.shutdown()
        server.server_close()
        thread.join(timeout=1)


def test_usi_aggregator_uses_other_backends_when_one_times_out(monkeypatch):
    response_gate = Event()
    slow_server, slow_thread, slow_request_started = _start_local_proxi_server(
        {"attributes": [], "mzs": [100.0], "intensities": [10.0]}, response_gate
    )
    fast_server, fast_thread, fast_request_started = _start_local_proxi_server(
        {"attributes": [], "mzs": [200.0], "intensities": [20.0]}
    )
    _route_proxi_backends_to_local_servers(
        monkeypatch, {"pride": slow_server, "massive": fast_server}
    )

    try:
        from proteomics_mcp.tools.workflows import retrieve_usi_spectrum_tool

        with ThreadPoolExecutor(max_workers=1) as executor:
            future = executor.submit(
                retrieve_usi_spectrum_tool,
                "mzspec:PRIDE:sample:scan:42",
                "aggregator",
                0.05,
            )
            try:
                assert slow_request_started.wait(timeout=1)
                assert fast_request_started.wait(timeout=1)
                response = future.result(timeout=0.5)
            finally:
                response_gate.set()

    finally:
        response_gate.set()
        for server, thread in ((slow_server, slow_thread), (fast_server, fast_thread)):
            server.shutdown()
            server.server_close()
            thread.join(timeout=1)

    assert response["mz"] == [200.0]
    assert response["intensity"] == [20.0]


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
