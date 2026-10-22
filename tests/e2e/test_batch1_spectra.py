"""Batch 1, item 4 of #30 (offline part): scoring, annotation and rendering of the reference payloads."""

from __future__ import annotations

import json
import sys
from importlib.metadata import version
from pathlib import Path

import pytest
from packaging.version import Version

pytestmark = pytest.mark.e2e

INSTALLED = Version(version("proteomics-mcp"))

REPO_EXAMPLES = Path(__file__).resolve().parents[2] / "examples"
# One PSM per dataset/engine; rank 2 of the HeLa set is the Figure 2 PSM (EFHLNESGDPSSK, scan 11648).
PAYLOADS = [
    "17_maxquant_top5/rank_1_scan_6637_AQNSCPICQK_payload.json",
    "18_fragpipe_top5/rank_3_scan_196914_ISM_UNIMOD_35_DEQDLLNGAGETYEVALTEPGTYSFYC_UNIMOD_4_APHQGAGM_UNIMOD_35_VGK_payload.json",
    "19_proteome_discoverer_top5/rank_1_scan_24381_DM_UNIMOD_35_HGVTSHLSNQELQDLVEFMK_payload.json",
    "22_hela_instanovo_top5/rank_2_scan_11648_EFHLNESGDPSSK_payload.json",
    "23_sbrodae_casanovo_top5/rank_2_scan_15106_EYLM_UNIMOD_35_FR_payload.json",
]
RESCORE_ATOL = 1e-9  # same inputs; allows cross-platform floating-point rounding (observed <= 4e-15)
METRICS = ("spectral_angle", "cosine_similarity", "pearson_correlation", "explained_intensity_fraction",
           "matched_ion_count", "log2_xcorr", "log2_hyperscore")
IDS = [name.split("/")[0] for name in PAYLOADS]


def _metric(result: dict, name: str):
    """A result metric under the name the 0.1.x payloads use: 0.2.0 renamed the legacy_v1 XCorr and hyperscore
    variants to *_like (#51), so log2_xcorr is log2_xcorr_like in a 0.2.0 legacy_v1 result."""
    return result[name] if name in result else result[f"{name}_like"]


def _load(name: str) -> dict:
    path = REPO_EXAMPLES / name
    if not path.is_file():
        pytest.skip("examples/ (reference payloads) is only available in a repository checkout")
    return json.loads(path.read_text())


def _settings(payload: dict) -> dict:
    settings = payload["metrics"]["annotation_settings"]
    return {"tolerance": settings["tolerance"], "tolerance_unit": settings["tolerance_unit"],
            "annotation_ion_types": settings["annotation_ion_types"], "annotation_max_isotope": settings["annotation_max_isotope"],
            "annotation_max_ion_charge": settings["annotation_max_ion_charge"],
            "annotation_neutral_losses": settings["annotation_neutral_losses"]}


@pytest.mark.parametrize("name", PAYLOADS, ids=IDS)
def test_rescoring_reproduces_stored_metrics(tool, cli, report, tmp_path, legacy_metrics, name):
    """Scoring engine only: the stored observed spectrum and stored prediction must reproduce the stored metrics."""
    payload = _load(name)
    is_error, result = tool("score_spectrum", {"candidate": payload["metrics"]["candidate"], "observed": payload["observed"],
                                               "predicted": payload["predicted"], "config_path": str(legacy_metrics),
                                               **_settings(payload)})
    assert not is_error, result
    deltas = {m: abs(_metric(result, m) - payload["metrics"][m]) for m in METRICS if payload["metrics"].get(m) is not None}

    # CLI parity for the same inputs
    files = {}
    for key, data in (("candidate", payload["metrics"]["candidate"]), ("spectrum", payload["observed"]), ("predicted", payload["predicted"])):
        files[key] = tmp_path / f"{key}.json"
        files[key].write_text(json.dumps(data))
    s = _settings(payload)
    args = ["score-spectrum", "--candidate-json", str(files["candidate"]), "--spectrum-json", str(files["spectrum"]),
            "--predicted-json", str(files["predicted"]), "--tolerance", str(s["tolerance"]), "--tolerance-unit", s["tolerance_unit"],
            "--annotation-ion-types", s["annotation_ion_types"], "--annotation-max-isotope", str(s["annotation_max_isotope"]),
            "--config", str(legacy_metrics)]
    if s["annotation_neutral_losses"]:
        args.append("--annotation-neutral-losses")
    completed = cli(*args)
    cli_ok = completed.returncode == 0
    cli_deltas = {}
    if cli_ok:
        cli_result = json.loads(completed.stdout)
        cli_deltas = {m: abs(_metric(cli_result, m) - _metric(result, m)) for m in deltas}
    report("spectra", f"score_spectrum (stored prediction) {name.split('/')[0]}", "pass" if max(deltas.values()) <= RESCORE_ATOL else "fail",
           {"max_abs_delta_vs_payload": deltas, "cli_returncode": completed.returncode, "cli_vs_mcp": cli_deltas,
            "cli_stderr": completed.stderr[-300:] if not cli_ok else ""})
    assert max(deltas.values()) <= RESCORE_ATOL, deltas
    assert cli_ok, completed.stderr[-2000:]
    assert max(cli_deltas.values(), default=0) <= RESCORE_ATOL, cli_deltas


@pytest.mark.parametrize("name", PAYLOADS, ids=IDS)
def test_annotation_and_score_summary(tool, report, name):
    payload = _load(name)
    settings = payload["metrics"]["annotation_settings"]
    is_error, annotation = tool("annotate_observed_spectrum", {
        "observed": payload["observed"], "peptidoform": payload["predicted"]["peptidoform"],
        "ion_types": settings["annotation_ion_types"], "tolerance": settings["tolerance"],
        "tolerance_unit": settings["tolerance_unit"], "max_isotope": settings["annotation_max_isotope"],
        "neutral_losses": settings["annotation_neutral_losses"]})
    assert not is_error, annotation
    stored = payload["metrics"]["annotation_summary"]
    same = (annotation["annotated_peak_count"] == stored["annotated_peak_count"]
            and abs(annotation["annotated_intensity_fraction"] - stored["annotated_intensity_fraction"]) <= 1e-12)
    is_error, summary = tool("summarize_score_metrics", {"score_result": payload["metrics"]})
    assert not is_error, summary
    report("spectra", f"annotate_observed_spectrum / summarize_score_metrics {name.split('/')[0]}", "pass" if same else "fail",
           {"annotated_peak_count": [annotation["annotated_peak_count"], stored["annotated_peak_count"]],
            "annotated_intensity_fraction": [annotation["annotated_intensity_fraction"], stored["annotated_intensity_fraction"]]})
    assert same


@pytest.mark.parametrize("name", PAYLOADS, ids=IDS)
def test_renderers(tool, cli, report, tmp_path, name):
    payload = _load(name)
    outputs = {}
    calls = {
        "render_spectrum_comparison": {"payload": payload},
        "render_observed_spectrum": {"payload_or_observed": payload},
        "render_annotated_observed_spectrum": {"payload_or_observed": payload, "peptidoform": payload["predicted"]["peptidoform"]},
    }
    for tool_name, arguments in calls.items():
        target = tmp_path / f"{tool_name}.html"
        is_error, result = tool(tool_name, {**arguments, "output_path": str(target)})
        assert not is_error, f"{tool_name}: {result}"
        assert target.is_file() and target.stat().st_size > 1000
        outputs[tool_name] = target.stat().st_size

    static = {}
    payload_file = tmp_path / "payload.json"
    payload_file.write_text(json.dumps(payload))
    cli_target = tmp_path / "static_cli.svg"
    completed = cli("render-static-spectrum", "--input-json", str(payload_file), "--output-path", str(cli_target), "--view", "annotated",
                    "--peptidoform", payload["predicted"]["peptidoform"])
    static["cli:annotated.svg"] = cli_target.stat().st_size if completed.returncode == 0 else "error: " + completed.stderr[-200:]
    report("spectra", f"renderers {name.split('/')[0]}", "pass", {"html_bytes": outputs, "static": static})
    assert completed.returncode == 0, completed.stderr[-2000:]


def _tk_unusable() -> bool:
    """Whether pyplot's Tk backend cannot create a figure here (broken Tcl/Tk on some Windows runners)."""
    import subprocess

    probe = "import matplotlib; matplotlib.use('TkAgg'); import matplotlib.pyplot as plt; plt.figure()"
    return subprocess.run([sys.executable, "-c", probe], capture_output=True).returncode != 0


@pytest.mark.skipif(sys.platform == "win32" and INSTALLED < Version("0.1.10") and _tk_unusable(),
                    reason="#35 (fixed in 0.1.10): with a broken Tcl/Tk, pyplot's TkAgg aborts the whole interpreter")
@pytest.mark.xfail(sys.platform == "darwin" and INSTALLED < Version("0.1.10"), strict=True,
                   reason="#35 (fixed in 0.1.10): pyplot selects the macosx backend, which fails in the MCP worker thread")
@pytest.mark.parametrize("name", PAYLOADS, ids=IDS)
def test_static_spectrum_via_mcp(tool, report, tmp_path, name):
    payload = _load(name)
    static = {}
    for view in ("mirror", "annotated", "plain"):
        for suffix in (".svg", ".png", ".pdf"):
            target = tmp_path / f"static_{view}{suffix}"
            is_error, result = tool("render_static_spectrum", {"payload_or_observed": payload, "output_path": str(target), "view": view,
                                                               "peptidoform": payload["predicted"]["peptidoform"]})
            static[f"{view}{suffix}"] = "error: " + str(result)[:120] if is_error else target.stat().st_size
    failed = {key: value for key, value in static.items() if isinstance(value, str)}
    report("spectra", f"render_static_spectrum via MCP {name.split('/')[0]}", "fail (#35)" if failed else "pass", {"outputs": static})
    assert not failed, failed


@pytest.mark.xfail(INSTALLED < Version("0.1.10"), strict=True,
                   reason="#45 (fixed in 0.1.10): render_spectra_batch reads only plain JSON")
def test_render_spectra_batch_reads_gzip_json_lines(tool, report, tmp_path):
    """§2.8: batch rendering from a file-backed JSON Lines manifest, gzip-compressed."""
    import gzip

    spec = tmp_path / "items.jsonl.gz"
    with gzip.open(spec, "wt", encoding="utf-8") as handle:
        for name in PAYLOADS:
            handle.write(json.dumps({"name": name.split("/")[0], "payload": _load(name)}) + "\n")
    is_error, result = tool("render_spectra_batch", {"input_path": str(spec), "validate_only": True})
    report("spectra", "render_spectra_batch from .jsonl.gz", "pass" if not is_error else "fail (#45)", {"result": str(result)[:300]})
    assert not is_error, result
    assert (result["items_requested"], result["valid_items"]) == (len(PAYLOADS), len(PAYLOADS))


def test_default_profile_scores_a_single_matched_peak_below_one(tool, report):
    """The default metric profile end to end (#72): keep only the observed peak matched to the strongest predicted ion.

    Under v2 the spectral angle is 1 - 2*arccos(p_max / ||p||)/pi; under legacy_v1 it would be 1.0.
    """
    import math

    payload = _load(next(name for name in PAYLOADS if "scan_11648" in name))  # the Figure 2 PSM
    stored = payload["metrics"]
    fragments = payload["predicted"]["fragments"]
    strongest = max(range(len(fragments)), key=lambda index: fragments[index]["intensity"])
    match = next((m for m in stored["matched_peaks"] if m["predicted_index"] == strongest), None)
    assert match is not None, "the strongest predicted ion of the reference payload is matched"
    observed = {**payload["observed"], "mz": [match["observed_mz"]], "intensity": [match["observed_intensity"]]}

    is_error, result = tool("score_spectrum", {"candidate": stored["candidate"], "observed": observed, "predicted": payload["predicted"],
                                               **_settings(payload)})
    assert not is_error, result
    if result.get("metric_profile") != "v2":
        pytest.skip(f"installed package defaults to metric profile {result.get('metric_profile')!r}")
    norm = math.sqrt(sum(f["intensity"] ** 2 for f in fragments))
    expected = 1 - 2 * math.acos(fragments[strongest]["intensity"] / norm) / math.pi
    report("spectra", "default profile: one matched peak (#72)", "pass" if abs(result["spectral_angle"] - expected) < 1e-9 else "fail",
           {"spectral_angle": result["spectral_angle"], "expected": expected, "matched_ion_count": result["matched_ion_count"]})
    assert result["matched_ion_count"] == 1
    assert result["spectral_angle"] == pytest.approx(expected, abs=1e-9)
    assert result["spectral_angle"] < 0.9
