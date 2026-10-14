"""Batch 1, item 4 of #30: PRIDE USI spectra vs the stored reference payloads.

Step a compares observed peaks (data source); step b compares metrics, and only where
the peaks match. Peak differences are data-source differences, not bugs.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import pytest

pytestmark = [pytest.mark.e2e, pytest.mark.network]

# examples/ folder -> dataset key in the Zenodo raw_file_inventory.tsv
EXAMPLE_DATASETS = {
    "17_maxquant_top5": "tum",
    "18_fragpipe_top5": "arabidopsis",
    "19_proteome_discoverer_top5": "sbrodae",
    "22_hela_instanovo_top5": "hela",
    "23_sbrodae_casanovo_top5": "sbrodae",
}
MZ_TOLERANCE = 1e-4  # Da
INTENSITY_RTOL = 1e-4  # relative
METRIC_ATOL = 1e-6  # absolute; MS2PIP 4.2.0 vs 4.1.2 differs by < 1e-7 for 3+ precursors (#23)
METRICS = ("spectral_angle", "cosine_similarity", "pearson_correlation", "explained_intensity_fraction")
REPO_EXAMPLES = Path(__file__).resolve().parents[2] / "examples"


def _payload_cases() -> list[tuple[str, Path]]:
    cases = []
    for folder in sorted(EXAMPLE_DATASETS):
        cases += [(folder, path) for path in sorted((REPO_EXAMPLES / folder).glob("*_payload.json"))]
    return cases


CASES = _payload_cases()


def _usi(zenodo: Path, folder: str, payload: dict) -> str:
    inventory = {row["dataset"]: row for row in csv.DictReader((zenodo / "raw_file_inventory.tsv").open(), delimiter="\t")}
    row = inventory[EXAMPLE_DATASETS[folder]]
    stem = row["pride_raw_filename"].removesuffix(".raw")
    return f"mzspec:{row['pride_accession']}:{stem}:scan:{payload['observed']['scan_number']}"


def compare_peaks(reference: dict, retrieved: dict) -> dict:
    """Match peaks by m/z within MZ_TOLERANCE; report counts and intensity agreement."""
    ref = sorted(zip(reference["mz"], reference["intensity"]))
    got = sorted(zip(retrieved["mz"], retrieved["intensity"]))
    matched, unmatched_ref, i = [], [], 0
    used = set()
    for mz, intensity in ref:
        while i < len(got) and got[i][0] < mz - MZ_TOLERANCE:
            i += 1
        j = i
        hit = None
        while j < len(got) and got[j][0] <= mz + MZ_TOLERANCE:
            if j not in used:
                hit = j
                break
            j += 1
        if hit is None:
            unmatched_ref.append(mz)
        else:
            used.add(hit)
            matched.append((intensity, got[hit][1]))
    intensity_ok = all(math.isclose(a, b, rel_tol=INTENSITY_RTOL) for a, b in matched)
    return {
        "reference_peaks": len(ref), "retrieved_peaks": len(got), "matched": len(matched),
        "missing_in_usi": len(unmatched_ref), "extra_in_usi": len(got) - len(used),
        "intensities_within_rtol": intensity_ok,
        "identical": len(matched) == len(ref) == len(got) and intensity_ok,
        "precursor_mz": [reference.get("precursor_mz"), retrieved.get("precursor_mz")],
        "precursor_charge": [reference.get("precursor_charge"), retrieved.get("precursor_charge")],
    }


def _score(tool, payload: dict, observed: dict, predicted: dict | None) -> dict:
    settings = payload["metrics"]["annotation_settings"]
    arguments = {
        "candidate": payload["metrics"]["candidate"], "observed": observed,
        "tolerance": settings["tolerance"], "tolerance_unit": settings["tolerance_unit"],
        "annotation_ion_types": settings["annotation_ion_types"], "annotation_max_isotope": settings["annotation_max_isotope"],
        "annotation_max_ion_charge": settings["annotation_max_ion_charge"],
        "annotation_neutral_losses": settings["annotation_neutral_losses"],
    }
    if predicted is not None:
        arguments["predicted"] = predicted
    else:
        arguments.update({"predictor_backend": "ms2pip", "ms2pip_model": payload["predicted"]["metadata"]["requested_model"]})
    is_error, result = tool("score_spectrum", arguments)
    assert not is_error, result
    return result


@pytest.fixture(scope="module")
def retrieved(zenodo, tool):
    """Retrieve each USI once per module."""
    cache = {}
    for folder, path in CASES:
        payload = json.loads(path.read_text())
        usi = _usi(zenodo, folder, payload)
        cache[path.name] = (usi, tool("retrieve_usi_spectrum", {"usi": usi, "timeout_seconds": 60}))
    return cache


@pytest.mark.parametrize(("folder", "path"), CASES, ids=[p.name.removesuffix("_payload.json") for _, p in CASES])
def test_usi_peaks_vs_payload(retrieved, report, folder, path):
    """Step a: observed peaks from PRIDE vs the payload (data-source comparison)."""
    payload = json.loads(path.read_text())
    usi, (is_error, spectrum) = retrieved[path.name]
    if is_error:
        report("usi", f"step a {path.name}", "fail", {"usi": usi, "error": spectrum})
        pytest.fail(f"retrieve_usi_spectrum failed for {usi}: {spectrum}")
    comparison = compare_peaks(payload["observed"], spectrum)
    report("usi", f"step a {path.name}", "identical" if comparison["identical"] else "differs", {"usi": usi, **comparison})


@pytest.mark.xfail(strict=True, reason="#41: the PRIDE precursor m/z (MS:1000744 / MS:1000827) is not parsed")
def test_usi_spectra_carry_precursor_mz(retrieved, report):
    missing = [usi for usi, (is_error, spectrum) in retrieved.values() if not is_error and spectrum.get("precursor_mz") is None]
    report("usi", "precursor_mz present in USI spectra", "fail (#41)" if missing else "pass", {"missing": len(missing), "of": len(retrieved)})
    assert not missing, f"{len(missing)} of {len(retrieved)} USI spectra have no precursor_mz"


@pytest.mark.parametrize(("folder", "path"), CASES, ids=[p.name.removesuffix("_payload.json") for _, p in CASES])
def test_usi_metrics_with_stored_prediction(retrieved, tool, report, folder, path):
    """Step b (scoring engine): only where peaks match; predictions taken from the payload."""
    payload = json.loads(path.read_text())
    usi, (is_error, spectrum) = retrieved[path.name]
    if is_error or not compare_peaks(payload["observed"], spectrum)["identical"]:
        pytest.skip("peaks differ from the reference payload (data source); metrics not compared")
    result = _score(tool, payload, spectrum, payload["predicted"])
    deltas = {name: abs(result[name] - payload["metrics"][name]) for name in METRICS}
    status = "pass" if max(deltas.values()) <= METRIC_ATOL else "fail"
    report("usi", f"step b stored prediction {path.name}", status, {"usi": usi, "max_abs_delta": deltas})
    assert status == "pass", deltas


@pytest.mark.usefixtures("ms2pip_runtime")
@pytest.mark.parametrize(("folder", "path"), CASES, ids=[p.name.removesuffix("_payload.json") for _, p in CASES])
def test_usi_metrics_with_live_ms2pip(retrieved, tool, report, folder, path):
    """Step b (live MS2PIP 4.2.0 in 0.1.9): only where peaks match."""
    payload = json.loads(path.read_text())
    usi, (is_error, spectrum) = retrieved[path.name]
    if is_error or not compare_peaks(payload["observed"], spectrum)["identical"]:
        pytest.skip("peaks differ from the reference payload (data source); metrics not compared")
    result = _score(tool, payload, spectrum, None)
    deltas = {name: abs(result[name] - payload["metrics"][name]) for name in METRICS}
    status = "pass" if max(deltas.values()) <= METRIC_ATOL else "fail"
    report("usi", f"step b live MS2PIP {path.name}", status, {"usi": usi, "charge": payload["predicted"]["precursor_charge"],
                                                               "max_abs_delta": deltas})
    assert status == "pass", deltas
