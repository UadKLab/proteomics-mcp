"""Batch 2 of #30: live prediction and the batch prediction -> scoring workflow on the published package.

- Live MS2PIP re-scoring of all 25 reference payloads (every spectrum shown in Figures 2 and 3 and
  Supplementary Figures 1, 2 and 5), compared with the stored metrics.
- The Figure 2 legend values, at the precision printed in the manuscript.
- Koina: `describe_predictor_model` and the §2.12 live check (PEPTIDE 2+, Prosit_2020_intensity_HCD).
- Batch workflow (§2.8, Supplementary Note 1 example 3): build_prediction_batch -> inspect -> score_from_prediction_batch,
  reuse of the prediction artifacts, isolation of a malformed row, and what an interrupted batch leaves behind.
"""

from __future__ import annotations

import csv
import json
import os
import signal
import subprocess
import sys
import time
from importlib.metadata import version
from pathlib import Path

import pytest
from packaging.version import Version

pytestmark = pytest.mark.e2e

INSTALLED = Version(version("proteomics-mcp"))

REPO_EXAMPLES = Path(__file__).resolve().parents[2] / "examples"
ALL_PAYLOADS = sorted(p.relative_to(REPO_EXAMPLES).as_posix() for p in REPO_EXAMPLES.glob("*_top5/rank_*_payload.json"))
FIGURE2 = "22_hela_instanovo_top5/rank_2_scan_11648_EFHLNESGDPSSK_payload.json"
HELA_RUN = "20230606_MK_E1200_MVL_FAIMS_2CV_15cm_70min_0000_DDA_200ng_QC1_Hela_QC"

# MS2PIP 4.2.0 vs the 4.1.2 used for the stored payloads: identical for 2+, <= 1e-7 for 3+ (HCDch2), see #23.
LIVE_ATOL = 1e-6
EXACT = ("matched_ion_count", "total_predicted_ions")
CLOSE = ("spectral_angle", "cosine_similarity", "pearson_correlation", "explained_intensity_fraction", "matched_ion_fraction")
RELATIVE = ("log2_xcorr", "log2_hyperscore")

# Figure 2 legend, at the printed precision.
FIGURE2_LEGEND = {
    "spectral_angle": (0.833, 3), "cosine_similarity": (0.966, 3), "pearson_correlation": (0.962, 3),
    "log2_xcorr": (18.765, 3), "log2_hyperscore": (5.983, 3),
    "matched_ion_count": (22, 0), "total_predicted_ions": (24, 0),
    "matched_ion_fraction": (91.7, 1, 100), "explained_intensity_fraction": (48.0, 1, 100),
    "precursor_signal_fraction": (10.3, 1, 100), "fragment_mass_error_abs_mean_ppm": (15.822, 3),
    "fragment_mass_error_abs_mean_da": (0.004, 3), "precursor_mass_error_ppm": (-3.055, 3),
}
# Fragment mass errors are measured against the predicted fragment m/z. MS2PIP 4.2.0 moved those by ~4.6e-5 Da (b1 of E:
# 130.04982 -> 130.04987, the theoretical value), so live values match the legend only with the 4.1.2 of the figures (#23).
PREDICTED_MZ_DEPENDENT = ("fragment_mass_error_abs_mean_ppm", "fragment_mass_error_abs_mean_da")
# Since #73 the theoretical precursor mass comes from exact residue masses instead of pyOpenMS, which moves the error by
# about +0.006 ppm: -3.0555 -> -3.0495, printed -3.050 instead of -3.055. Packages with #73 report precursor_isotope_error.
PRECURSOR_MASS_DEPENDENT = ("precursor_mass_error_ppm",)
PRECURSOR_MASS_ATOL_PPM = 0.01

LIBOMP_HINT = "brew install libomp"


def ms2pip_available() -> bool:
    """MS2PIP needs a loadable xgboost (on macOS: Homebrew libomp)."""
    try:
        import xgboost  # noqa: F401
    except Exception:
        return False
    return True


def _load(name: str) -> dict:
    path = REPO_EXAMPLES / name
    if not path.is_file():
        pytest.skip("examples/ (reference payloads) is only available in a repository checkout")
    return json.loads(path.read_text())


def _settings(payload: dict) -> dict:
    s = payload["metrics"]["annotation_settings"]
    return {"tolerance": s["tolerance"], "tolerance_unit": s["tolerance_unit"], "annotation_ion_types": s["annotation_ion_types"],
            "annotation_max_isotope": s["annotation_max_isotope"], "annotation_max_ion_charge": s["annotation_max_ion_charge"],
            "annotation_neutral_losses": s["annotation_neutral_losses"]}


def _deltas(result: dict, stored: dict) -> dict[str, float]:
    deltas = {}
    for m in EXACT + CLOSE:
        if stored.get(m) is not None:
            deltas[m] = abs(result[m] - stored[m])
    for m in RELATIVE:
        if stored.get(m):
            deltas[m] = abs(result[m] - stored[m]) / abs(stored[m])
    return deltas


def _printed(value: float, spec: tuple) -> float:
    printed, digits, *scale = spec
    return round(value * (scale[0] if scale else 1), digits)


def _require_ms2pip(tool, report, item: str) -> None:
    """Without a loadable MS2PIP, the tool must fail with the actionable libomp error (#25), then skip."""
    if ms2pip_available():
        return
    is_error, message = tool("predict_fragments", {"peptidoform": "PEPTIDE", "charge": 2})
    report("prediction", f"{item}: MS2PIP unavailable", "info", {"is_error": is_error, "message": str(message)[:300]})
    if sys.platform == "darwin":
        assert is_error and LIBOMP_HINT in str(message), message
    pytest.skip("MS2PIP cannot load xgboost here (macOS without Homebrew libomp, see #17)")


# --- live MS2PIP on the reference payloads -------------------------------------------------------------------------


@pytest.mark.parametrize("name", ALL_PAYLOADS, ids=[n.replace("/", ":")[:48] for n in ALL_PAYLOADS])
def test_live_ms2pip_reproduces_stored_metrics(tool, report, name):
    payload = _load(name)
    _require_ms2pip(tool, report, "live re-scoring")
    is_error, result = tool("score_spectrum", {"candidate": payload["metrics"]["candidate"], "observed": payload["observed"],
                                               "predictor_backend": "ms2pip", "predictor_model": "HCD", **_settings(payload)})
    assert not is_error, result
    stored = payload["metrics"]
    deltas = _deltas(result, stored)
    printed_equal = all(round(result[m], 3) == round(stored[m], 3) for m in CLOSE if stored.get(m) is not None)
    if INSTALLED >= Version("0.1.10"):  # #22: the MS2PIP release is recorded with the score
        assert (result.get("predictor_metadata") or {}).get("predictor_version") == version("ms2pip")
    report("prediction", f"live MS2PIP re-score {name}", "pass" if max(deltas.values()) <= LIVE_ATOL else "fail",
           {"resolved_model": (result.get("predictor_metadata") or {}).get("model"), "max_delta": deltas,
            "equal_at_3_decimals": printed_equal})
    for m in EXACT:
        assert deltas.get(m, 0) == 0, (m, result[m], stored[m])
    assert max(deltas.values()) <= LIVE_ATOL, deltas
    assert printed_equal


@pytest.mark.parametrize("source", ["stored", "live"])
def test_figure2_legend_values(tool, report, source):
    """Every number in the Figure 2 legend, from the stored payload and from a live MS2PIP prediction."""
    payload = _load(FIGURE2)
    if source == "stored":
        metrics = payload["metrics"]
    else:
        _require_ms2pip(tool, report, "Figure 2 live")
        is_error, metrics = tool("score_spectrum", {"candidate": payload["metrics"]["candidate"], "observed": payload["observed"],
                                                    "predictor_backend": "ms2pip", "predictor_model": "HCD", **_settings(payload)})
        assert not is_error, metrics
    mismatches = {m: (_printed(metrics[m], spec), spec[0]) for m, spec in FIGURE2_LEGEND.items()
                  if _printed(metrics[m], spec) != spec[0]}
    version_dependent = {}
    if source == "live" and Version(version("ms2pip")) != Version("4.1.2"):
        version_dependent = {m: mismatches.pop(m) for m in PREDICTED_MZ_DEPENDENT if m in mismatches}
    if source == "live" and "precursor_isotope_error" in metrics:
        for m in PRECURSOR_MASS_DEPENDENT:
            if m in mismatches:
                assert abs(metrics[m] - FIGURE2_LEGEND[m][0]) <= PRECURSOR_MASS_ATOL_PPM, (m, metrics[m])
                version_dependent[m] = mismatches.pop(m)
    candidate = payload["metrics"]["candidate"]
    if round(candidate["confidence"] * 100, 1) != 100.0:
        mismatches["confidence"] = (candidate["confidence"], 100.0)
    if f"{candidate['raw_score']:.2e}" != "-2.80e-04":
        mismatches["raw_score"] = (candidate["raw_score"], -2.80e-4)
    report("figures", f"Figure 2 legend values ({source})", "fail" if mismatches else "differs (#23, #73)" if version_dependent else "pass",
           {"mismatches": mismatches, "ms2pip": version("ms2pip"), "version_dependent": version_dependent})
    assert not mismatches, mismatches


# --- Koina ---------------------------------------------------------------------------------------------------------


def test_koina_describe_and_predict(tool, cli, report):
    """§2.4 describe_predictor_model, and the §2.12 live check: PEPTIDE 2+ gives b/y charge-1 ions only."""
    is_error, described = tool("describe_predictor_model", {"model": "Prosit_2020_intensity_HCD"})
    assert not is_error, described
    inputs = {i["name"] for i in described["metadata"]["inputs"]}
    is_error, predicted = tool("predict_fragments", {"peptidoform": "PEPTIDE", "charge": 2, "backend": "koina",
                                                     "model": "Prosit_2020_intensity_HCD"})
    assert not is_error, predicted
    kinds = sorted({(f["ion_type"], f["charge"]) for f in predicted["fragments"]})
    completed = cli("predict-fragments", "--peptidoform", "PEPTIDE", "--charge", "2", "--backend", "koina",
                    "--model", "Prosit_2020_intensity_HCD")
    cli_fragments = json.loads(completed.stdout)["fragments"] if completed.returncode == 0 else None
    report("prediction", "Koina describe + predict (PEPTIDE 2+)", "pass", {
        "inputs": sorted(inputs), "supported": described.get("supported_by_proteomics_mcp"),
        "fragments": len(predicted["fragments"]), "ion_kinds": kinds, "predictor_version": predicted.get("predictor_version"),
        "cli_returncode": completed.returncode})
    assert described.get("supported_by_proteomics_mcp") is True
    assert {"peptide_sequences", "precursor_charges"} <= inputs
    assert kinds == [("b", 1), ("y", 1)]
    assert all(f["mz"] > 0 for f in predicted["fragments"])
    assert completed.returncode == 0, completed.stderr[-2000:]
    assert [(f["annotation"], f["mz"]) for f in cli_fragments] == [(f["annotation"], f["mz"]) for f in predicted["fragments"]]
    if INSTALLED >= Version("0.1.10"):  # #22
        assert predicted["predictor_version"] == "Prosit_2020_intensity_HCD@1", predicted["predictor_version"]


# --- batch prediction -> scoring ------------------------------------------------------------------------------------


def _write_mgf(path: Path, payloads: list[dict]) -> None:
    """The observed spectra of the HeLa reference payloads, as an MGF named after their run."""
    lines = []
    for payload in payloads:
        o = payload["observed"]
        lines += ["BEGIN IONS", f"TITLE={o['spectrum_id']}", f"SCANS={o['scan_number']}", f"PEPMASS={o['precursor_mz']!r}",
                  f"CHARGE={o['precursor_charge']}+"]
        lines += [f"{mz!r} {intensity!r}" for mz, intensity in zip(o["mz"], o["intensity"])]
        lines.append("END IONS")
    path.write_text("\n".join(lines) + "\n")


def _write_candidates(path: Path, payloads: list[dict]) -> None:
    """InstaNovo-style candidate table for the same PSMs, plus one row with an unknown modification."""
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["experiment_name", "scan_number", "precursor_charge", "predictions", "log_probs"])
        for payload in payloads:
            c = payload["metrics"]["candidate"]
            writer.writerow([HELA_RUN, c["scan_number"], c["precursor_charge"], c["peptidoform"], c["raw_score"]])
        # malformed: the peptidoform cannot be predicted; it must fail on its own without stopping the batch
        writer.writerow([HELA_RUN, 99999999, 2, "PEPT[UNIMOD:999999]IDEK", -0.5])


@pytest.fixture(scope="module")
def hela_batch_inputs(tmp_path_factory, tool):
    payloads = [_load(p) for p in ALL_PAYLOADS if p.startswith("22_hela_instanovo_top5/")]
    root = tmp_path_factory.mktemp("batch2")
    mgf = root / f"{HELA_RUN}.mgf"
    _write_mgf(mgf, payloads)
    candidates = root / "hela_candidates.csv"
    _write_candidates(candidates, payloads)
    is_error, converted = tool("mgf_to_feather", {"mgf_path": str(mgf), "output_path": str(root / "hela_spectra.feather")})
    assert not is_error, converted
    return {"root": root, "payloads": payloads, "candidates": candidates, "spectra": root / "hela_spectra.feather"}


def _backend() -> dict:
    # MS2PIP where it loads (the figures' predictor); Koina otherwise, so the workflow is still exercised.
    if ms2pip_available():
        return {"predictor_backend": "ms2pip", "predictor_model": "HCD"}
    return {"predictor_backend": "koina", "predictor_model": "Prosit_2020_intensity_HCD"}


def _read_csv(path: str) -> list[dict]:
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def test_prediction_and_scoring_batch(tool, report, hela_batch_inputs):
    inputs = hela_batch_inputs
    backend = _backend()
    prefix = inputs["root"] / "run1" / "hela"
    is_error, built = tool("build_prediction_batch", {"path": str(inputs["candidates"]), "output_prefix": str(prefix),
                                                      "spectra_path": str(inputs["spectra"]), "adapter": "instanovo", **backend})
    assert not is_error, built
    manifest = built["output_files"]["manifest_json"]

    is_error, inspected = tool("inspect_prediction_batch", {"manifest_path": manifest})
    assert not is_error, inspected
    is_error, status = tool("inspect_prediction_batch_status", {"manifest_path": manifest})
    assert not is_error, status
    failures = _read_csv(json.loads(Path(manifest).read_text())["output_files"]["failures_csv"])

    settings = _settings(inputs["payloads"][0])
    score_args = {"manifest_path": manifest, "spectra_path": str(inputs["spectra"]), **{k: v for k, v in settings.items()}}
    is_error, scored = tool("score_from_prediction_batch", {**score_args, "output_prefix": str(inputs["root"] / "run1" / "scored")})
    assert not is_error, scored
    scored_manifest = scored["output_files"]["manifest_json"]
    is_error, scored_status = tool("inspect_scored_batch_status", {"manifest_path": scored_manifest})
    assert not is_error, scored_status
    is_error, scored_inspected = tool("inspect_scored_batch", {"manifest_path": scored_manifest})
    assert not is_error, scored_inspected
    scored_rows = _read_csv(json.loads(Path(scored_manifest).read_text())["output_files"]["scored_csv"])

    # Reuse: scoring again from the same prediction artifacts must give the same scores, without predicting again.
    is_error, rescored = tool("score_from_prediction_batch", {**score_args, "output_prefix": str(inputs["root"] / "run2" / "scored")})
    assert not is_error, rescored
    rescored_manifest = rescored["output_files"]["manifest_json"]
    rescored_rows = _read_csv(json.loads(Path(rescored_manifest).read_text())["output_files"]["scored_csv"])
    key = lambda rows: {(r["candidate_id"], r.get("spectral_angle")) for r in rows}  # noqa: E731

    # With MS2PIP, the batch scores must match the stored per-PSM metrics.
    deltas = {}
    if backend["predictor_backend"] == "ms2pip":
        by_scan = {int(float(r["scan_number"])): r for r in scored_rows if r.get("spectral_angle")}
        for payload in inputs["payloads"]:
            row = by_scan.get(payload["observed"]["scan_number"])
            assert row is not None, f"scan {payload['observed']['scan_number']} not scored"
            deltas[payload["observed"]["scan_number"]] = max(abs(float(row[m]) - payload["metrics"][m]) for m in
                                                             ("spectral_angle", "cosine_similarity", "pearson_correlation"))

    report("batch", "build_prediction_batch -> score_from_prediction_batch (HeLa reference PSMs)", "pass", {
        "backend": backend, "prediction_status": {k: status.get(k) for k in ("phase", "completed", "total", "succeeded", "failed")},
        "prediction_failures": [(f["peptidoform"], f["error_type"]) for f in failures],
        "scored_status": {k: scored_status.get(k) for k in ("phase", "completed", "total", "succeeded", "failed")},
        "scored_rows": len(scored_rows), "rescored_identical": key(scored_rows) == key(rescored_rows),
        "max_delta_vs_payload_by_scan": deltas})
    assert status["phase"] == "completed" and status["succeeded"] == len(inputs["payloads"]), status
    assert [f["peptidoform"] for f in failures] == ["PEPT[UNIMOD:999999]IDEK"], failures
    assert key(scored_rows) == key(rescored_rows)
    assert all(d <= LIVE_ATOL for d in deltas.values()), deltas


@pytest.fixture(scope="module")
def interrupted_batch(hela_batch_inputs):
    """Start build-prediction-batch in a subprocess and stop it as soon as its status file appears."""
    inputs = hela_batch_inputs
    backend = _backend()
    prefix = inputs["root"] / "interrupted" / "hela"
    prefix.parent.mkdir(parents=True, exist_ok=True)
    status_json = prefix.with_name(prefix.name + "_prediction_status.json")
    manifest_json = prefix.with_name(prefix.name + "_prediction_manifest.json")
    bindir = Path(sys.executable).parent
    script = next(bindir / n for n in ("proteomics-mcp", "proteomics-mcp.exe") if (bindir / n).is_file())
    process = subprocess.Popen([str(script), "build-prediction-batch", "--path", str(inputs["candidates"]),
                                "--output-prefix", str(prefix), "--spectra-path", str(inputs["spectra"]),
                                "--predictor-backend", backend["predictor_backend"], "--predictor-model", backend["predictor_model"]],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    deadline = time.time() + 120
    while not status_json.exists() and process.poll() is None and time.time() < deadline:
        time.sleep(0.01)
    if process.poll() is None:
        if os.name == "nt":
            process.kill()
        else:
            process.send_signal(signal.SIGINT)
    process.wait(timeout=60)
    if manifest_json.exists() and json.loads(manifest_json.read_text()).get("phase", "completed") == "completed":
        pytest.skip("the batch finished before it could be interrupted")
    return {"dir": prefix.parent, "status_json": status_json, "manifest_json": manifest_json, "backend": backend}


def test_interrupted_batch_keeps_only_its_status(report, interrupted_batch):
    """No resume is claimed (§2.8 says reusable artifacts); record what an interrupted batch leaves behind."""
    left = sorted(p.name for p in interrupted_batch["dir"].iterdir())
    status = json.loads(interrupted_batch["status_json"].read_text()) if interrupted_batch["status_json"].exists() else None
    report("batch", "interrupted build_prediction_batch", "info", {
        "backend": interrupted_batch["backend"], "files_left": left,
        "status": {k: (status or {}).get(k) for k in ("phase", "completed", "total")}})
    assert status is not None and status["phase"] == "predicting", status
    # results are written at the end; only the status (and, from 0.1.10, the running manifest) survive
    expected = [interrupted_batch["status_json"].name]
    if INSTALLED >= Version("0.1.10"):
        expected = sorted(expected + [interrupted_batch["manifest_json"].name])
    assert left == expected, left


@pytest.mark.xfail(INSTALLED < Version("0.1.10"), strict=True,
                   reason="#46 (fixed in 0.1.10): the status tool reads the manifest, which is written only when the batch finishes")
def test_interrupted_batch_status_can_be_inspected(tool, interrupted_batch):
    is_error, result = tool("inspect_prediction_batch_status", {"manifest_path": str(interrupted_batch["manifest_json"])})
    assert not is_error, result
    assert result["phase"] != "completed" and result["seconds_since_update"] >= 0, result
