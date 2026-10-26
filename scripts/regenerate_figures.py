"""Regenerate the manuscript figures that change with 0.2.0 (Figures 2, 3, 4 and Supplementary Figures 2 and 5).

Inputs are the manuscript's Zenodo supplement (record 23295902, md5-checked) and, for Figure 4, PRIDE spectra by USI.
Every PSM is resolved from the supplement's identification tables with the installed package's adapters; observed
spectra come from the supplement's reference payloads, so the figures show the same spectra as the published ones.
Scoring uses the package defaults (metric profile v2) with the published tolerance and annotation settings, and the
panels use the published render settings. Nothing is written inside the repository.

Run from a checkout (the release candidate or the tag):

    uv run python scripts/regenerate_figures.py --out ~/pmcp/figures-v0.2.0 --cache-dir ~/pmcp/cache

The figures use MS2PIP 4.1.2, as the published ones did (#23). The package itself does not pin MS2PIP, so the script
runs itself again in a uv environment with this checkout, every other dependency at its uv.lock version, MS2PIP 4.1.2
and the figure-assembly tools (FIGURE_ENVIRONMENT), unless it already runs in one.

How the figures are made: every number, score, annotation and plot comes from the package, through the functions behind
its MCP tools (load_candidates, score_spectrum, predict_fragments, the prediction batch and its scoring, build_id_dashboard,
render_static_spectrum, render_spectrum_comparison, the observed-spectrum views, render_psm_table, render_peptide_table,
inspect_peptide_table). The script computes no score and draws no data. It only assembles: it places the package's panels
on the published page widths and adds the panel letters (in place of the manual Illustrator step), zooms and crops the
package's interactive plots (Supplementary Figure 2B/2C) and captures them, and the Figure 4 dashboard plots, as
screenshots of at least MIN_RASTER_DPI in the figure, and writes the PRIDE spectra of Figure 4 to an MGF at full precision.

HTML panels (Figure 2B, Supplementary Figure 2, Figure 4) are captured with Playwright in the local Google Chrome
(`--browser-channel chrome`) or Playwright's own Chromium (`--browser-channel chromium`). Exit codes: 0 done,
3 a verification against the expected 0.2.0 values failed (nothing is assembled), 4 an input checksum differs.
"""

from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.request
import warnings
import zipfile
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parents[1]
ZENODO_RECORD = "23295902"
ZENODO_URL = f"https://zenodo.org/records/{ZENODO_RECORD}/files/supplement.zip?download=1"
ZENODO_MD5 = "1445e7af348b12c703f477d50fb01130"
EXAMPLES = "figure_assets/examples"
PUBLISHED_PDFS = "figure_assets/pdfs"

# --- the published figures: PSMs, settings and layout ---------------------------------------------------------------

# One PSM per spectrum panel: the identification table it is resolved from, and the reference payload whose observed
# spectrum the published panel shows. Casanovo reports no scan; its published scan link came from the matched
# Proteome Discoverer export, so the scan is taken from the payload and the candidate from the adapter.
PSMS: dict[str, dict[str, Any]] = {
    "fig2": {"adapter": "instanovo", "table": "processed/hela/instanovo_labelled_kpreds.csv", "scan": 11648,
             "payload": "22_hela_instanovo_top5/rank_2_scan_11648_EFHLNESGDPSSK_payload.json"},
    "fig3a": {"adapter": "proteome_discoverer", "table": "processed/sbrodae/proteome_discoverer_PSMs.xlsx", "scan": 24381,
              "payload": "19_proteome_discoverer_top5/rank_1_scan_24381_DM_UNIMOD_35_HGVTSHLSNQELQDLVEFMK_payload.json"},
    "fig3b": {"adapter": "maxquant", "table": "processed/proteome_tools/maxquant_msms.txt", "scan": 6637,
              "payload": "17_maxquant_top5/rank_1_scan_6637_AQNSCPICQK_payload.json"},
    "s5a": {"adapter": "maxquant", "table": "processed/proteome_tools/maxquant_msms.txt", "scan": 38257, "title": "MaxQuant",
            "payload": "17_maxquant_top5/rank_2_scan_38257_HTYDDEFFTFK_payload.json"},
    "s5b": {"adapter": "fragpipe", "table": "processed/arabidopsis/fragpipe_psm.tsv", "scan": 183225, "title": "FragPipe",
            "payload": "18_fragpipe_top5/rank_1_scan_183225_FIHVSTDEVYGETDEDALVGNHEASQLLPTNPYSATK_payload.json"},
    "s5c": {"adapter": "proteome_discoverer", "table": "processed/sbrodae/proteome_discoverer_PSMs.xlsx", "scan": 17557,
            "title": "Proteome Discoverer", "payload": "19_proteome_discoverer_top5/rank_3_scan_17557_YVSRPGIISEDDDTEILGSDR_payload.json"},
    "s5d": {"adapter": "casanovo", "table": "processed/sbrodae/casanovo_predictions.csv", "candidate_id": "8352", "title": "Casanovo",
            "payload": "23_sbrodae_casanovo_top5/rank_3_scan_10031_FLYNTR_payload.json"},
}
PREDICTOR = {"predictor_backend": "ms2pip", "predictor_model": "HCD"}
FIGURE_MS2PIP = "4.1.2"  # the release behind the published figures (#23); pinned here only, not in the package
FIGURE_ENVIRONMENT = [f"ms2pip=={FIGURE_MS2PIP}", "playwright==1.55.0", "pypdf==6.21.1", "pypdfium2==5.14.0"]
MIRROR_LABELS = 14  # the published mirror plots label the 14 most intense matched peaks
ANNOTATED_LABELS = 18  # the published Figure 3 panels label the 18 most intense annotated peaks
FIG3A_ION_TYPES = "a,b,c,x,y,z,z-dot,immonium,precursor"  # the ion families in the published Figure 3A legend
FIG3B_ION_TYPES = "all"
METRIC_STRIP_VIEWPORT = 1400  # five metric tiles per row, as in Figure 2B
# Panels that are screenshots of the package's interactive plots; the rest are vector PDFs from the package's renderers.
RASTER_PANELS = {"Suppl_figure2B", "Suppl_figure2C", *(f"Figure4{letter}" for letter in "ABCDE")}
MIN_RASTER_DPI = 300
METADATA_VIEWPORT = 960  # four metadata columns, as in Supplementary Figure 2A
S2B_MZ_RANGE = (195.0, 1085.0)
S2C_MZ_RANGE = (413.6, 420.4)
S2C_HOVER_MZ = 418.228  # y4+ | PSSK
# Page widths of the published PDFs (points); heights follow from the panels.
PAGE_WIDTH = {"Figure2": 527.33, "Figure3": 527.33, "Figure4": 543.81, "Suppl_figure2": 611.47, "Suppl_figure5": 1085.82}
FIGURE4_PLOTS = {"A": "scatter_spectral_angle_vs_fragment_mass_error", "B": "hist_rt", "C": "hist_spectral_angle",
                 "D": "hist_matched_fraction", "E": "hist_explained_fraction"}
FIGURE4_SUBSET_ROWS = 100  # the Figure 4 subset (#42, #47): the first 100 PSMs in the deposited FragPipe PSM table, retained in source-file order
FIGURE4_RUN = "20180905_QX4_JoMu_SA_EASY12-7_uPAC_500ng_ArabidopsisiRT-F5"
FIGURE4_PRIDE = "PXD014877"
FIGURE4_RT_TOLERANCE_MINUTES = 1e-5  # FragPipe prints Retention in seconds to 4 decimals
# sha256 of the Figure 4 MGF (the 100 PRIDE spectra with their retention times, full precision, 2026-10-26): PRIDE serving other
# spectra, or another --figure4-mgf, fails the build's checks
FIGURE4_MGF_SHA256 = "e1d2ee5b8632f36b4ae8291196735b4ef7939c045ff6ffba57d2c3ee8ceaf059"
# Where the Figure 4 build lands in the Zenodo deposit (deposit/ in the output folder; scripts/build_zenodo_package.py)
DEPOSIT_ARABIDOPSIS = "processed/arabidopsis"
DEPOSIT_INTERACTIVE = "figure_assets/interactive"
FIGURE4_DEPOSIT = {
    "spectra_mgf": f"{DEPOSIT_ARABIDOPSIS}/figure4_usi_spectra.mgf",  # new: the PRIDE spectra, so Figure 4 rebuilds from Zenodo alone
    "scored_psms_full": f"{DEPOSIT_ARABIDOPSIS}/figure4_scored_psms_full.csv",  # new: the batch scoring output, every column
    "scored_ions": f"{DEPOSIT_ARABIDOPSIS}/figure4_scored_ions.jsonl",  # new: the ions behind each score (#137)
    "prediction_manifest": "regeneration/figure4_prediction_manifest.json",  # new: build_prediction_batch's manifest
    "scoring_manifest": "regeneration/figure4_scoring_manifest.json",  # new: score_from_prediction_batch's manifest
    "scored_psms": f"{DEPOSIT_ARABIDOPSIS}/supplementary_scored_psms.csv",
    "scored_peptides": f"{DEPOSIT_ARABIDOPSIS}/supplementary_scored_peptides.csv",
    "plot_data": f"{DEPOSIT_ARABIDOPSIS}/figure4_dashboard_plot_data.json",
    "dashboard_html": f"{DEPOSIT_INTERACTIVE}/Figure_3D_arabidopsis_fragpipe_score_aware_dashboard.html",
    "psm_table_html": f"{DEPOSIT_INTERACTIVE}/Supplementary_Figure_3_fragpipe_scored_psm_table.html",
    "peptide_table_html": f"{DEPOSIT_INTERACTIVE}/Supplementary_Figure_4_fragpipe_scored_peptide_table.html",
}

STRIP_METRICS = [  # Figure 2B: metric, label, printed format
    ("spectral_angle", "Spectral angle", "{:.3f}"), ("cosine_similarity", "Cosine", "{:.3f}"),
    ("pearson_correlation", "Pearson", "{:.3f}"), ("xcorr", "XCorr (Comet)", "{:.3f}"),
    ("hyperscore", "Hyperscore (X! Tandem)", "{:.1f}"), ("matched_ion_count", "Matched ions (observed)", "{:d}"),
    ("total_predicted_ions", "Matched ions (predicted)", "{:d}"), ("matched_ion_fraction", "Matched fraction", "pct"),
    ("explained_intensity_fraction", "Explained intensity", "pct"), ("precursor_signal_fraction", "Precursor signal", "pct"),
    ("fragment_mass_error_abs_mean_ppm", "Fragment abs ppm", "{:.3f}"), ("fragment_mass_error_abs_mean_da", "Fragment abs Da", "{:.3f}"),
    ("precursor_mass_error_ppm", "Precursor ppm", "{:.3f}"), ("precursor_mass_error_da", "Precursor Da", "{:.3f}"),
]
# The published strip showed the 0.1.x variants in these places (#51): CHANGES.md compares the new scores with them.
PUBLISHED_AS = {"xcorr": ("log2_xcorr", "log2 XCorr, 0.1.x variant"), "hyperscore": ("log2_hyperscore", "log2 Hyperscore, 0.1.x variant")}

# What 0.2.0 must give (issue #147). A mismatch stops the run before anything is assembled.
EXPECTED = {
    # xcorr and hyperscore: the values Comet 2026.01.1 (comet -p defaults) and X! Tandem 2015.12.15.2 (default input,
    # 0.4 Da) report for this PSM, run from the BioContainers images pinned by digest in #154 (#51)
    "fig2": {"spectral_angle": "0.834", "pearson_correlation": "0.960", "xcorr": "3.336", "hyperscore": "59.5", "matched_ion_count": "18",
             "total_predicted_ions": "18", "matched_ion_fraction": "100.0%", "explained_intensity_fraction": "47.1%",
             "precursor_mass_error_ppm": "-3.050", "ion_counts": {"x": 8, "precursor": 5},
             "fragment_mass_error_abs_mean_ppm": "2.796",  # MS2PIP 4.1.2, the figures' release
             "fragment_mass_error_abs_mean_ppm_by_ms2pip": {"4.2.0": "2.867"},  # the --compare-ms2pip side-by-side
             "mirror_labels_removed": ["b11+", "b1+", "b2+", "b8+"]},
    "fig3a": {"labelled_peaks": 237},
    "fig3b": {"peptidoform": "AQNSC[UNIMOD:4]PIC[UNIMOD:4]QK", "labelled_peaks": 110, "precursor_abs_ppm_max": 0.05,
              "legend_without": ["a+1", "d-alt"], "legend_with": ["precursor"]},
    "mirror_labels_removed": {"s5a": [], "s5b": ["y22++", "y23++", "y24++"], "s5c": ["b3+", "b6+", "y14+", "y15+", "y15++", "y17++"], "s5d": []},
}
# sha256 of the Figure 3A peak labels (ion types of its legend, 0.05 Da, max isotope 2, neutral losses) from 0.1.10,
# the version of the published figure: 3A must come out with exactly these labels.
FIG3A_LABELS_SHA256 = "a02da17e233777f6ada441a933c58a533aa7bbe1f192ed36d2916f67c623c44c"

CAUSES = {  # metric -> issue/PR, for CHANGES.md
    "spectral_angle": "#72 (PR #99), #74 (PR #103), #83 (PR #106): metric profile v2",
    "cosine_similarity": "#72, #74, #83: metric profile v2",
    "pearson_correlation": "#72, #74, #83: metric profile v2",
    "xcorr": "#51 (PR #152): Comet's XCorr replaces the 0.1.x variant (log2 of a binned dot product)",
    "hyperscore": "#51 (PR #152): X! Tandem's hyperscore replaces the 0.1.x variant (log2 of ln n! plus summed intensities)",
    "matched_ion_count": "#83 (PR #106): predicted ions below 0.1% of the strongest are not scored; #74 (PR #103)",
    "total_predicted_ions": "#83 (PR #106): predicted ions below 0.1% of the strongest are not scored",
    "matched_ion_fraction": "#83 (PR #106), #74 (PR #103)",
    "explained_intensity_fraction": "#74 (PR #103): one peak explains one ion; #83 (PR #106)",
    "precursor_signal_fraction": "-",
    "fragment_mass_error_abs_mean_ppm": "#83 (PR #106): faint predicted ions no longer matched to noise peaks up to 0.05 Da away; #74",
    "fragment_mass_error_abs_mean_da": "#83 (PR #106), #74 (PR #103)",
    "precursor_mass_error_ppm": "#73 (PR #100): exact residue masses, every modification",
    "precursor_mass_error_da": "#73 (PR #100)",
}


# --- small helpers -----------------------------------------------------------------------------------------------------


def log(message: str) -> None:
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {message}", flush=True)


def file_digest(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def printed(value: Any, spec: str) -> str:
    if value is None:
        return "NA"
    if spec == "pct":
        return f"{100 * value:.1f}%"
    if spec == "{:d}":
        return str(int(value))
    return spec.format(value)


def installed_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not installed"


def package_versions() -> dict[str, str]:
    names = ("proteomics-mcp", "ms2pip", "matplotlib", "seaborn", "numpy", "pandas", "pyopenms", "playwright", "pypdf", "pypdfium2")
    return {name: installed_version(name) for name in names}


def user_command() -> str:
    """The command as run from the checkout, without the figure environment's re-run (--in-figure-environment)."""
    script = Path(sys.argv[0]).resolve()
    shown = str(script.relative_to(REPO_ROOT)) if REPO_ROOT in script.parents else str(script)
    return shlex.join(["uv", "run", "python", shown, *(a for a in sys.argv[1:] if a != "--in-figure-environment")])


def git_state() -> dict[str, Any]:
    def run(*args: str) -> str:
        try:
            return subprocess.run(["git", *args], cwd=REPO_ROOT, capture_output=True, text=True, check=True).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            return ""

    return {"sha": run("rev-parse", "HEAD"), "describe": run("describe", "--tags", "--always", "--dirty"),
            "branch": run("rev-parse", "--abbrev-ref", "HEAD"), "dirty": bool(run("status", "--porcelain", "--untracked-files=no"))}


# --- inputs ----------------------------------------------------------------------------------------------------------


def zenodo_supplement(zip_path: Path | None, cache_dir: Path) -> tuple[Path, Path]:
    """The supplement archive (downloaded once into the cache unless given), md5-checked, and its extracted folder."""
    archive = zip_path or cache_dir / f"zenodo-{ZENODO_RECORD}-supplement.zip"
    if not archive.is_file():
        log(f"downloading the Zenodo supplement (record {ZENODO_RECORD})")
        partial = archive.with_suffix(".part")
        with urllib.request.urlopen(ZENODO_URL, timeout=300) as response, partial.open("wb") as handle:
            shutil.copyfileobj(response, handle)
        partial.replace(archive)
    checksum = file_digest(archive, "md5")
    if checksum != ZENODO_MD5:
        log(f"Zenodo supplement md5 {checksum} != {ZENODO_MD5}: not the deposited supplement")
        sys.exit(4)
    target = cache_dir / f"zenodo-{ZENODO_RECORD}"
    if not (target / "supplement" / "raw_file_inventory.tsv").is_file():
        with zipfile.ZipFile(archive) as bundle:
            members = [m for m in bundle.namelist() if not m.startswith("__MACOSX/") and not m.endswith(".DS_Store")]
            bundle.extractall(target, members=members)
    return archive, target / "supplement"


def resolve_candidate(supplement: Path, spec: dict[str, Any], cache: dict[str, list[dict[str, Any]]]) -> dict[str, Any]:
    """The candidate as the installed package's adapter reads it from the published identification table."""
    from proteomics_mcp.tools.workflows import load_candidates_tool

    table = str(supplement / spec["table"])
    if table not in cache:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cache[table] = load_candidates_tool(table, adapter=spec["adapter"])
    rows = cache[table]
    if "candidate_id" in spec:
        matches = [row for row in rows if str(row.get("candidate_id")) == spec["candidate_id"]]
    else:
        matches = [row for row in rows if row.get("scan_number") == spec["scan"]]
    peptidoforms = {row["peptidoform"] for row in matches}
    if len(matches) != 1 and len(peptidoforms) != 1:
        raise SystemExit(f"{spec['table']}: expected one PSM for {spec.get('scan') or spec.get('candidate_id')}, found {len(matches)}")
    return dict(matches[0])


def score_psm(candidate: dict[str, Any], observed: dict[str, Any], settings: dict[str, Any]) -> dict[str, Any]:
    """Score one PSM with the package defaults; returns the render payload (observed, predicted, metrics)."""
    from proteomics_mcp.tools.workflows import predict_fragments_tool, score_spectrum_tool

    charge = candidate.get("precursor_charge") or observed.get("precursor_charge") or 2
    # The package predicts inside score_spectrum, which applies the fragment-charge rule (#137 concerns only
    # predictions passed in from outside). The returned spectrum is the same prediction, kept for the mirror panel.
    metrics = score_spectrum_tool(candidate, observed, **PREDICTOR, **settings)
    predicted = predict_fragments_tool(candidate["peptidoform"], int(charge), model=PREDICTOR["predictor_model"], backend=PREDICTOR["predictor_backend"])
    supplied = score_spectrum_tool(candidate, observed, predicted=predicted, **settings)
    for key in ("spectral_angle", "pearson_correlation", "matched_ion_count", "total_predicted_ions", "explained_intensity_fraction"):
        if supplied.get(key) != metrics.get(key):
            raise SystemExit(f"{candidate['peptidoform']}: {key} differs between predicted and supplied fragments (#137)")
    if metrics.get("metric_profile") != "v2":
        raise SystemExit(f"metric profile is {metrics.get('metric_profile')}, expected the 0.2.0 default v2")
    return {"observed": observed, "predicted": predicted, "metrics": metrics}


def annotation(payload: dict[str, Any], ion_types: str) -> dict[str, Any]:
    from proteomics_mcp.core.annotation import annotate_observed_ions
    from proteomics_mcp.models import ObservedSpectrum

    s = payload["metrics"]["annotation_settings"]
    return annotate_observed_ions(ObservedSpectrum(**payload["observed"]), peptidoform=payload["metrics"]["candidate"]["peptidoform"],
                                  tolerance=s["tolerance"], tolerance_unit=s["tolerance_unit"], ion_types=ion_types, max_isotope=2,
                                  neutral_losses=True)


def mirror_labels(payload: dict[str, Any]) -> set[tuple[int, str]]:
    return {(m["observed_index"], m["annotation"]) for m in payload["metrics"]["matched_peaks"]}


# --- verification ------------------------------------------------------------------------------------------------------


def verify(results: dict[str, dict[str, Any]], ms2pip_version: str, fig2_by_version: dict[str, dict[str, Any]]) -> list[str]:
    problems: list[str] = []
    fig2 = results["fig2"]["new"]["metrics"]
    specs = {key: spec for key, _, spec in STRIP_METRICS}
    for key in ("spectral_angle", "pearson_correlation", "xcorr", "hyperscore", "matched_ion_count", "total_predicted_ions",
                "matched_ion_fraction", "explained_intensity_fraction", "precursor_mass_error_ppm"):
        got = printed(fig2.get(key), specs[key])
        if got != EXPECTED["fig2"][key]:
            problems.append(f"Figure 2 {key}: {got}, expected {EXPECTED['fig2'][key]}")
    counts = (fig2.get("annotation_summary") or {}).get("ion_type_counts") or {}
    for ion, count in EXPECTED["fig2"]["ion_counts"].items():
        if counts.get(ion) != count:
            problems.append(f"Figure 2 ion tile {ion}: {counts.get(ion)}, expected {count}")
    if ms2pip_version != FIGURE_MS2PIP:
        problems.append(f"MS2PIP {ms2pip_version}, expected {FIGURE_MS2PIP} for the figures (#23)")
    got = printed(fig2.get("fragment_mass_error_abs_mean_ppm"), "{:.3f}")
    if got != EXPECTED["fig2"]["fragment_mass_error_abs_mean_ppm"]:
        problems.append(f"Figure 2 fragment abs ppm: {got}, expected {EXPECTED['fig2']['fragment_mass_error_abs_mean_ppm']}")
    for ms2pip, metrics in fig2_by_version.items():
        want = EXPECTED["fig2"]["fragment_mass_error_abs_mean_ppm_by_ms2pip"].get(ms2pip)
        got = printed(metrics.get("fragment_mass_error_abs_mean_ppm"), "{:.3f}")
        if want and got != want:
            problems.append(f"Figure 2 fragment abs ppm with MS2PIP {ms2pip}: {got}, expected {want}")
    for key in ("fig2", "s5a", "s5b", "s5c", "s5d"):
        old, new = mirror_labels(results[key]["published"]), mirror_labels(results[key]["new"])
        removed, added = sorted(label for _, label in old - new), sorted(label for _, label in new - old)
        want = EXPECTED["fig2"]["mirror_labels_removed"] if key == "fig2" else EXPECTED["mirror_labels_removed"][key]
        if removed != sorted(want) or added:
            problems.append(f"{key} mirror labels: removed {removed}, added {added}; expected removed {sorted(want)}, none added")
    fig3a = results["fig3a"]
    if fig3a["new"]["metrics"]["candidate"]["peptidoform"] != fig3a["published"]["metrics"]["candidate"]["peptidoform"]:
        problems.append("Figure 3A peptidoform differs from the published payload")
    labels = annotation(fig3a["new"], FIG3A_ION_TYPES)
    if labels["annotated_peak_count"] != EXPECTED["fig3a"]["labelled_peaks"]:
        problems.append(f"Figure 3A labelled peaks: {labels['annotated_peak_count']}, expected {EXPECTED['fig3a']['labelled_peaks']}")
    if hashlib.sha256(json.dumps(labels["annotations"]).encode()).hexdigest() != FIG3A_LABELS_SHA256:
        problems.append("Figure 3A peak labels differ from the published (0.1.10) labels")
    fig3b = results["fig3b"]["new"]
    expected_3b = EXPECTED["fig3b"]
    if fig3b["metrics"]["candidate"]["peptidoform"] != expected_3b["peptidoform"]:
        problems.append(f"Figure 3B peptidoform {fig3b['metrics']['candidate']['peptidoform']}, expected {expected_3b['peptidoform']}")
    labels_3b = annotation(fig3b, FIG3B_ION_TYPES)
    if labels_3b["annotated_peak_count"] != expected_3b["labelled_peaks"]:
        problems.append(f"Figure 3B labelled peaks: {labels_3b['annotated_peak_count']}, expected {expected_3b['labelled_peaks']}")
    precursor = fig3b["metrics"].get("precursor_mass_error_ppm")
    if precursor is None or abs(precursor) > expected_3b["precursor_abs_ppm_max"]:
        problems.append(f"Figure 3B precursor error {precursor} ppm, expected ~0")
    present = set(labels_3b["ion_type_counts"])
    for ion in expected_3b["legend_without"]:
        if ion in present:
            problems.append(f"Figure 3B legend still has {ion}")
    for ion in expected_3b["legend_with"]:
        if ion not in present:
            problems.append(f"Figure 3B legend lacks {ion}")
    return problems


# --- pinned environments -------------------------------------------------------------------------------------------------


def pinned_environment(packages: list[str], python: str, workdir: Path, name: str) -> list[str]:
    """`uv run` with this checkout, every other dependency at its uv.lock version, and `packages` (pins that replace or
    add to the locked ones); returns the command prefix up to and including `python`."""
    exported = subprocess.run(["uv", "export", "--frozen", "--no-hashes", "--no-emit-project", "--no-dev", "--no-header",
                               "--python", python], cwd=REPO_ROOT, capture_output=True, text=True, check=True).stdout
    replaced = tuple(f"{package.split('==')[0]}==" for package in packages)
    lines = [line for line in exported.splitlines() if line.strip() and not line.lstrip().startswith(("#", *replaced))]
    workdir.mkdir(parents=True, exist_ok=True)
    requirements = workdir / f"requirements_{name}.txt"
    requirements.write_text("\n".join([*lines, *packages, ""]))
    return ["uv", "run", "--no-project", "--python", python, "--with-requirements", str(requirements), "--with", str(REPO_ROOT),
            "--reinstall-package", "proteomics-mcp",  # build this checkout, never a cached build of another commit
            "python"]


def rerun_in_figure_environment(args: argparse.Namespace, argv: list[str]) -> int:
    """Run this script again in the figure environment (MS2PIP FIGURE_MS2PIP and the assembly tools)."""
    out = args.out.expanduser().resolve()
    command = [*pinned_environment(FIGURE_ENVIRONMENT, args.python, out / "work", "figures"), str(Path(__file__).resolve()),
               *argv, "--in-figure-environment"]
    log(f"MS2PIP {installed_version('ms2pip')} installed; running in the figure environment (MS2PIP {FIGURE_MS2PIP}): {shlex.join(command)}")
    return subprocess.run(command).returncode


# --- MS2PIP comparison for Figure 2 -------------------------------------------------------------------------------------


def figure2_metrics_with(ms2pip: str, args: argparse.Namespace, workdir: Path) -> dict[str, Any]:
    """Figure 2 metrics in a separate environment with another MS2PIP release: the package from this checkout and every
    other dependency at its uv.lock version, so MS2PIP is the only difference."""
    out = workdir / f"figure2_ms2pip_{ms2pip}.json"
    command = [*pinned_environment([f"ms2pip=={ms2pip}"], args.python, workdir, f"ms2pip_{ms2pip}"), str(Path(__file__).resolve()),
               "--figure2-metrics-only", str(out), "--cache-dir", str(args.cache_dir)]
    if args.zenodo_zip:
        command += ["--zenodo-zip", str(args.zenodo_zip)]
    log(f"Figure 2 with MS2PIP {ms2pip}: {shlex.join(command)}")
    completed = subprocess.run(command, capture_output=True, text=True)
    if completed.returncode != 0:
        raise SystemExit(f"Figure 2 with MS2PIP {ms2pip} failed:\n{completed.stderr[-3000:]}")
    result = json.loads(out.read_text())
    result["command"] = shlex.join(command)
    return result


# --- HTML panels (Playwright) ---------------------------------------------------------------------------------------


class Browser:
    def __init__(self, channel: str):
        from playwright.sync_api import sync_playwright

        self._playwright = sync_playwright().start()
        self.browser = self._playwright.chromium.launch(**({"channel": channel} if channel != "chromium" else {}))

    def page(self, html: Path, width: int, height: int = 900):
        page = self.browser.new_page(viewport={"width": width, "height": height}, device_scale_factor=3)
        page.goto(html.resolve().as_uri())
        page.wait_for_function("window.Plotly !== undefined || !document.querySelector('.js-plotly-plot, #plot')", timeout=60000)
        page.wait_for_timeout(1500)
        return page

    def close(self) -> None:
        self.browser.close()
        self._playwright.stop()


def isolate_and_capture(page, selector_js: str, stem: Path) -> None:
    """Capture the elements returned by `selector_js` as PNG (screenshot) and vector PDF (the elements alone)."""
    page.evaluate(f"""() => {{
        const nodes = ({selector_js})();
        const width = Math.max(...nodes.map(n => n.getBoundingClientRect().width));  // keep the page's column layout
        const box = document.createElement('div'); box.id = 'capture-box';
        box.style.cssText = `display:inline-block;padding:8px;background:#fff;width:${{width}}px`;
        nodes.forEach(n => {{ const c = n.cloneNode(true); c.style.width = `${{width}}px`; c.style.boxSizing = 'border-box'; box.appendChild(c); }});
        document.body.innerHTML = ''; document.body.appendChild(box);
        document.body.style.margin = '0';
    }}""")
    element = page.locator("#capture-box")
    element.screenshot(path=str(stem.with_suffix(".png")))
    box = element.bounding_box()
    page.pdf(path=str(stem.with_suffix(".pdf")), width=f"{box['width'] + 1}px", height=f"{box['height'] + 1}px", print_background=True,
             margin={"top": "0", "bottom": "0", "left": "0", "right": "0"}, page_ranges="1")


def capture_plot(page, stem: Path, relayout: dict[str, Any] | None = None, size: tuple[int, int] | None = None,
                 hover_mz: float | None = None) -> None:
    """Screenshot the Plotly plot (#plot), optionally zoomed, resized and with a peak hover card."""
    if size:  # the div keeps the plot's size, so the screenshot holds the plot and nothing below it
        page.evaluate(f"() => {{ const gd = document.getElementById('plot'); gd.style.width = '{size[0]}px'; gd.style.height = '{size[1]}px';"
                      f" return Plotly.relayout(gd, {{width: {size[0]}, height: {size[1]}}}); }}")
    if relayout:
        page.evaluate(f"() => Plotly.relayout('plot', {json.dumps(relayout)})")
    page.evaluate("() => document.querySelectorAll('.plot-controls').forEach(e => { e.style.display = 'none'; })")
    page.wait_for_timeout(800)
    if hover_mz is not None:
        found = page.evaluate(f"""() => {{
            const gd = document.getElementById('plot');
            for (let c = 0; c < gd.data.length; c++) {{
                const t = gd.data[c];
                if (!String(t.name || '').endsWith('hover area')) continue;
                for (let i = 0; i < (t.x || []).length; i++) {{
                    if (t.x[i] !== null && Math.abs(t.x[i] - {hover_mz}) < 0.002 && t.y[i] > 0) {{
                        Plotly.Fx.hover(gd, [{{curveNumber: c, pointNumber: i}}]); return t.name;
                    }}
                }}
            }}
            return null;
        }}""")
        if not found:
            raise SystemExit(f"no observed peak at m/z {hover_mz} to hover")
        page.wait_for_timeout(500)
    page.locator("#plot").screenshot(path=str(stem.with_suffix(".png")))
    png_to_pdf(stem.with_suffix(".png"), stem.with_suffix(".pdf"), dpi=288)


def png_to_pdf(png: Path, pdf: Path, dpi: int) -> None:
    from PIL import Image

    Image.open(png).convert("RGB").save(pdf, "PDF", resolution=dpi)


# --- figure assembly ---------------------------------------------------------------------------------------------------


def compose(rows: list[list[tuple[str, Path, float]]], width: float, output: Path, gap: float = 10.0, gutter: float = 16.0) -> dict[str, float]:
    """Lay out panel PDFs in rows ((letter, pdf, fraction of the width)), letter at each panel's top left, as the
    published figures are (letters in a left gutter of each panel); the vector content of each panel is kept.
    Writes output (.pdf) and a PNG next to it; returns each panel's placed width in points, by panel file stem."""
    from matplotlib.figure import Figure
    from pypdf import PageObject, PdfReader, PdfWriter, Transformation

    placed: list[tuple[str, Any, float, float, float]] = []  # letter, page, x, y_top, scale
    widths: dict[str, float] = {}
    y_top = 0.0
    for row in rows:
        x = 0.0
        row_height = 0.0
        for letter, pdf, fraction in row:
            page = PdfReader(str(pdf)).pages[0]
            panel_width = width * fraction - gutter - (gap if len(row) > 1 else 0)
            scale = panel_width / float(page.mediabox.width)
            placed.append((letter, page, x + (gap / 2 if len(row) > 1 and x > 0 else 0), y_top, scale))
            widths[Path(pdf).stem] = panel_width
            row_height = max(row_height, float(page.mediabox.height) * scale)
            x += width * fraction
        y_top += row_height + gap
    height = y_top - gap
    letters = Figure(figsize=(width / 72, height / 72))
    for letter, _, x, top, _ in placed:
        letters.text(x / width, 1 - top / height, letter, fontsize=16, va="top", ha="left", family="DejaVu Sans")
    buffer = io.BytesIO()
    letters.savefig(buffer, format="pdf", transparent=True)
    overlay = PdfReader(buffer).pages[0]
    page = PageObject.create_blank_page(width=width, height=height)
    for letter, panel, x, top, scale in placed:
        y = height - top - float(panel.mediabox.height) * scale
        page.merge_transformed_page(panel, Transformation().translate(-float(panel.mediabox.left), -float(panel.mediabox.bottom)).scale(scale).translate(x + gutter, y))
    page.merge_page(overlay)
    writer = PdfWriter()
    writer.add_page(page)
    with output.open("wb") as handle:
        writer.write(handle)
    pdf_to_png(output, output.with_suffix(".png"))
    return widths


def raster_dpi(panels: Path, widths: dict[str, float]) -> dict[str, float]:
    """Effective resolution of each screenshot panel in its figure: PNG pixels per inch of placed width."""
    from PIL import Image

    return {stem: round(Image.open(panels / f"{stem}.png").width / (widths[stem] / 72), 1) for stem in RASTER_PANELS if stem in widths}


def pdf_to_png(pdf: Path, png: Path, longest: int = 2400) -> None:
    import pypdfium2

    document = pypdfium2.PdfDocument(str(pdf))
    page = document[0]
    width, height = page.get_size()
    page.render(scale=longest / max(width, height)).to_pil().save(png)


def side_by_side(published_pdf: Path, new_png: Path, output: Path, label_new: str) -> None:
    """Published figure (left) next to the regenerated one (right)."""
    from matplotlib import font_manager
    from PIL import Image, ImageDraw, ImageFont

    with tempfile.TemporaryDirectory() as tmp:
        published_png = Path(tmp) / "published.png"
        pdf_to_png(published_pdf, published_png)
        left, right = Image.open(published_png).convert("RGB"), Image.open(new_png).convert("RGB")
    height = max(left.height, right.height)
    left = left.resize((round(left.width * height / left.height), height)) if left.height != height else left
    right = right.resize((round(right.width * height / right.height), height)) if right.height != height else right
    header, margin = 90, 40
    canvas = Image.new("RGB", (left.width + right.width + 3 * margin, height + header + margin), "white")
    canvas.paste(left, (margin, header))
    canvas.paste(right, (2 * margin + left.width, header))
    draw = ImageDraw.Draw(canvas)
    font = ImageFont.truetype(font_manager.findfont("DejaVu Sans"), 44)
    draw.text((margin, 20), "Published (Zenodo 23295902, proteomics-mcp 0.1.10 era)", fill="black", font=font)
    draw.text((2 * margin + left.width, 20), label_new, fill="black", font=font)
    draw.line([(1.5 * margin + left.width, header), (1.5 * margin + left.width, height + header)], fill="#999999", width=3)
    canvas.save(output)


# --- Figure 4 (PRIDE USI) ----------------------------------------------------------------------------------------------


def fetch_usi_spectra(scans: list[tuple[int, int]], cache_dir: Path, timeout: float, retries: int) -> tuple[dict[int, dict[str, Any]], list[int]]:
    """PRIDE spectra for (scan, charge), cached per scan; returns spectra and the scans PRIDE did not serve."""
    from proteomics_mcp.tools.workflows import retrieve_usi_spectrum_tool

    cache = cache_dir / "usi" / FIGURE4_PRIDE
    cache.mkdir(parents=True, exist_ok=True)
    spectra, missing = {}, []
    for number, (scan, _) in enumerate(scans, 1):
        path = cache / f"{FIGURE4_RUN}.{scan}.json"
        if path.is_file():
            cached = json.loads(path.read_text())
            if cached.get("rt_seconds") is not None:  # cached before #135 (PR #156): fetched again for its retention time
                spectra[scan] = cached
                continue
        usi = f"mzspec:{FIGURE4_PRIDE}:{FIGURE4_RUN}:scan:{scan}"
        for attempt in range(1, retries + 1):
            try:
                spectrum = retrieve_usi_spectrum_tool(usi, "pride", timeout_seconds=timeout)
                if spectrum.get("mz"):
                    path.write_text(json.dumps(spectrum))
                    spectra[scan] = spectrum
                    log(f"USI {number}/{len(scans)} scan {scan}: {len(spectrum['mz'])} peaks")
                    break
            except Exception as exc:  # noqa: BLE001 - PRIDE answers 500 or times out intermittently; retried
                log(f"USI {number}/{len(scans)} scan {scan}, attempt {attempt}/{retries}: {type(exc).__name__}: {str(exc)[:120]}")
                time.sleep(min(60, 5 * attempt))
        else:
            missing.append(scan)
    return spectra, missing


def write_figure4_mgf(scans: list[tuple[int, int]], spectra: dict[int, dict[str, Any]], mgf: Path) -> list[int]:
    """The PRIDE spectra as one MGF, with their retention times (#135, PR #156); returns the scans without one.

    Values are written at full precision (repr round-trips every float), so the MGF holds exactly what PRIDE served;
    the package has no MGF writer."""
    without_rt = []
    with mgf.open("w", newline="\n") as handle:
        for scan, charge in scans:
            spectrum = spectra[scan]
            handle.write(f"BEGIN IONS\nTITLE={FIGURE4_RUN}.{scan}.{scan}.{charge}\nSCANS={scan}\n")
            if spectrum.get("rt_seconds") is None:
                without_rt.append(scan)
            else:
                handle.write(f"RTINSECONDS={float(spectrum['rt_seconds'])!r}\n")
            handle.write(f"PEPMASS={float(spectrum['precursor_mz'])!r}\nCHARGE={spectrum.get('precursor_charge') or charge}+\n")
            handle.writelines(f"{float(mz)!r} {float(intensity)!r}\n" for mz, intensity in zip(spectrum["mz"], spectrum["intensity"]))
            handle.write("END IONS\n")
    return without_rt


def batch_vs_direct_scores(candidates: Path, scored_csv: str, feather: Path) -> dict[str, Any]:
    """#137: the batch scores the predictions of its prediction step. Re-score every PSM with score_spectrum, where the
    package predicts, against the same spectrum table, and report rows whose metrics differ."""
    import pandas as pd

    from proteomics_mcp.core.spectrum_loader import load_spectrum_from_table
    from proteomics_mcp.tools.workflows import load_candidates_tool, score_spectrum_tool

    scored = pd.read_csv(scored_csv).set_index("candidate_id")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        rows = {row["candidate_id"]: row for row in load_candidates_tool(str(candidates), adapter="fragpipe")}
    keys = ("spectral_angle", "cosine_similarity", "pearson_correlation", "xcorr", "hyperscore", "matched_ion_count", "total_predicted_ions",
            "explained_intensity_fraction")
    differing, compared, largest = [], 0, 0.0
    for candidate_id, row in rows.items():
        if candidate_id not in scored.index:
            continue
        observed = load_spectrum_from_table(feather, scan_number=row["scan_number"]).model_dump()
        direct = score_spectrum_tool(row, observed, **PREDICTOR, tolerance=0.05, tolerance_unit="Da")
        compared += 1
        batch = scored.loc[candidate_id]
        pairs = {k: (float(batch[k]), float(direct[k])) for k in keys if pd.notna(batch.get(k)) and direct.get(k) is not None}
        largest = max([largest, *(abs(a - b) for a, b in pairs.values())])
        deltas = {k: v for k, v in pairs.items() if abs(v[0] - v[1]) > 1e-6}  # the batch CSV round trip moves floats by ~1e-9
        if deltas:
            differing.append({"candidate_id": candidate_id, "deltas": deltas})
    return {"compared": compared, "differing_above_1e-6": len(differing), "largest_abs_difference": largest, "examples": differing[:5]}


def deposited_column_order(new: Any, deposited: Path) -> tuple[Any, dict[str, list[str]]]:
    """The regenerated table with the deposited table's column order (new columns last), and the columns that differ."""
    import pandas as pd

    old_columns = list(pd.read_csv(deposited, nrows=0).columns)
    kept = [c for c in old_columns if c in new.columns]
    added = [c for c in new.columns if c not in old_columns]
    return new[kept + added], {"removed": [c for c in old_columns if c not in new.columns], "added": added}


def figure4(supplement: Path, args: argparse.Namespace, work: Path, panels: Path, browser: Browser | None, deposit: Path) -> dict[str, Any]:
    """Figure 4 from the first 100 PSMs in the deposited FragPipe PSM table, retained in source-file order (#42, #47),
    and their PRIDE spectra by USI, with the tables and interactive files of the deposit. The spectra go through one
    MGF, which the deposit carries, so `--figure4-mgf` rebuilds the figure from the Zenodo files alone."""
    import pandas as pd

    from proteomics_mcp.tools.workflows import (
        build_id_dashboard_tool,
        build_prediction_batch_tool,
        inspect_peptide_table_tool,
        mgf_to_feather_tool,
        render_peptide_table_tool,
        render_psm_table_tool,
        score_from_prediction_batch_tool,
    )

    status: dict[str, Any] = {"status": "PROVISIONAL", "reasons": [], "notes": [
        "PSMs: the first 100 PSMs in the deposited FragPipe PSM table, retained in source-file order (fragpipe_psm.tsv, "
        "the subset of Supplementary Figures 3 and 4), as decided for 0.2.0 (#42, #47); the published Figure 4 showed "
        "~960 PSMs of an unrecorded selection",
        "spectra: PRIDE's USI rendering of the scans (PXD014877), not the local mzML conversion; PRIDE serves fewer peaks than "
        "the local files for 23 of 25 reference PSMs (#31). The MGF of these spectra is in the deposit",
        "histograms: the published dashboard left out the PSMs at each histogram's maximum (1-2 per panel); fixed in #169, and "
        "checked here: every histogram panel counts all PSMs",
    ]}
    source = supplement / "processed/arabidopsis/fragpipe_psm.tsv"
    table = pd.read_csv(source, sep="\t", low_memory=False)
    subset = table.head(FIGURE4_SUBSET_ROWS)
    deposited = pd.read_csv(supplement / "processed/arabidopsis/supplementary_scored_psms.csv", usecols=["candidate_id", "scan_number"])
    scans = [(int(s.split(".")[-2]), int(c)) for s, c in zip(subset["Spectrum"], subset["Charge"])]
    if [s for s, _ in scans] != deposited["scan_number"].astype(int).tolist():
        raise SystemExit("the first 100 FragPipe PSMs are not the deposited Figure 4 subset")
    work4 = work / "figure4"
    work4.mkdir(parents=True, exist_ok=True)
    mgf = work4 / f"{FIGURE4_RUN}.mgf"
    if args.figure4_mgf:
        shutil.copyfile(args.figure4_mgf.expanduser(), mgf)
        status["spectra"] = {"source": f"--figure4-mgf {args.figure4_mgf}"}
    else:
        spectra, missing = fetch_usi_spectra(scans, args.cache_dir, args.usi_timeout, args.usi_retries)
        status["spectra"] = {"source": f"PRIDE USI ({FIGURE4_PRIDE})", "requested": len(scans), "retrieved": len(spectra), "missing_scans": missing,
                             "cache_sha256": {str(s): file_digest(args.cache_dir / "usi" / FIGURE4_PRIDE / f"{FIGURE4_RUN}.{s}.json") for s in sorted(spectra)}}
        if missing:
            status["status"] = "SKIPPED"
            status["reasons"].append(f"PRIDE did not serve {len(missing)} of {len(scans)} spectra by USI")
            return status
        without_rt = write_figure4_mgf(scans, spectra, mgf)
        if without_rt:
            status["reasons"].append(f"{len(without_rt)} PRIDE spectra came without a retention time (#135): {without_rt[:10]}")
    status["spectra"]["mgf_sha256"] = file_digest(mgf)
    # The package's own MGF reader builds the spectrum table; every score and retention time below comes from it.
    feather = work4 / f"{FIGURE4_RUN}.feather"
    mgf_to_feather_tool(str(mgf), str(feather))
    candidates = work4 / "psm.tsv"  # adapters derive candidate_id from the file name; the deposited subset uses "psm:<row>"
    subset.to_csv(candidates, sep="\t", index=False)
    batch = build_prediction_batch_tool(path=str(candidates), output_prefix=str(work4 / "batch"), spectra_path=str(feather), adapter="fragpipe", **PREDICTOR)
    scored = score_from_prediction_batch_tool(manifest_path=batch["output_files"]["manifest_json"],
                                              spectra_path=str(feather), output_prefix=str(work4 / "scored"), tolerance=0.05, tolerance_unit="Da")
    scored_csv = scored["output_files"]["scored_csv"]
    dashboard = build_id_dashboard_tool(path=str(candidates), adapter="fragpipe", scored_path=str(scored_csv), spectra_path=str(feather),
                                        output_prefix=str(work4 / "figure4"))
    enriched = pd.read_csv(dashboard["output_files"]["enriched_csv"], low_memory=False)
    status["issue_137"] = batch_vs_direct_scores(candidates, scored_csv, feather)
    rt = enriched.set_index("scan_number")["retention_time_minutes"] if "retention_time_minutes" in enriched else pd.Series(dtype=float)
    fragpipe_rt = pd.Series(subset["Retention"].to_numpy() / 60, index=[s for s, _ in scans])
    common = rt.dropna().index.intersection(fragpipe_rt.index)
    status["retention_time_check"] = {"compared": len(common),
                                      "largest_abs_difference_minutes": float((rt.loc[common] - fragpipe_rt.loc[common]).abs().max()) if len(common) else None}
    join = dashboard.get("metadata_join") or {}
    status["batch"] = {"prediction": {k: v for k, v in batch.items() if isinstance(v, (int, float, str))},
                       "scoring": {k: v for k, v in scored.items() if isinstance(v, (int, float, str))},
                       "dashboard_rows": len(enriched), "rows_with_spectral_angle": int(enriched["spectral_angle"].notna().sum()),
                       "rows_with_rt": int(enriched.get("retention_time_minutes", pd.Series(dtype=float)).notna().sum()),
                       "metadata_join": join}
    # The deposit: the spectra, the scored PSM and peptide tables, the dashboard and the two interactive tables.
    targets = {key: deposit / path for key, path in FIGURE4_DEPOSIT.items()}
    for path in targets.values():
        path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(mgf, targets["spectra_mgf"])
    psms, status["scored_psms_columns"] = deposited_column_order(enriched, supplement / FIGURE4_DEPOSIT["scored_psms"])
    psms.to_csv(targets["scored_psms"], index=False)
    shutil.copyfile(scored_csv, targets["scored_psms_full"])
    shutil.copyfile(scored["output_files"]["scored_ions_jsonl"], targets["scored_ions"])
    shutil.copyfile(batch["output_files"]["manifest_json"], targets["prediction_manifest"])
    shutil.copyfile(scored["output_files"]["manifest_json"], targets["scoring_manifest"])
    # the inspect_peptide_table tool: every peptide of the scored subset, as the package tabulates them
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        inspected = inspect_peptide_table_tool(path=str(candidates), adapter="fragpipe", scored_path=str(scored_csv), limit=len(subset))
    peptides, status["scored_peptides_columns"] = deposited_column_order(pd.DataFrame(inspected["peptide_preview"]),
                                                                         supplement / FIGURE4_DEPOSIT["scored_peptides"])
    status["scored_peptides_rows"] = len(peptides)
    peptides.to_csv(targets["scored_peptides"], index=False)
    shutil.copyfile(dashboard["output_files"]["dashboard_html"], targets["dashboard_html"])
    sdrf = supplement / DEPOSIT_ARABIDOPSIS / "sdrf.tsv"
    table_options = {"path": str(candidates), "adapter": "fragpipe", "scored_path": str(scored_csv), "sdrf_path": str(sdrf) if sdrf.is_file() else None}
    render_psm_table_tool(output_path=str(targets["psm_table_html"]), **table_options)
    render_peptide_table_tool(output_path=str(targets["peptide_table_html"]), **table_options)

    checks = {
        "the spectra are the recorded ones (FIGURE4_MGF_SHA256)": status["spectra"]["mgf_sha256"] == FIGURE4_MGF_SHA256,
        "100 PSMs scored": status["batch"]["rows_with_spectral_angle"] == FIGURE4_SUBSET_ROWS == len(enriched),
        "every PSM joined to its spectrum (#146, PR #158)": join.get("matched_rows") == FIGURE4_SUBSET_ROWS and not join.get("unmatched_rows"),
        "every PSM has the PRIDE retention time (#135, PR #156)": status["batch"]["rows_with_rt"] == FIGURE4_SUBSET_ROWS,
        "retention times equal FragPipe's": status["retention_time_check"]["compared"] == FIGURE4_SUBSET_ROWS
        and status["retention_time_check"]["largest_abs_difference_minutes"] <= FIGURE4_RT_TOLERANCE_MINUTES,
        "batch scores equal score_spectrum (#137, PR #150)": status["issue_137"]["differing_above_1e-6"] == 0
        and status["issue_137"]["compared"] == FIGURE4_SUBSET_ROWS,
    }
    status["checks"] = checks
    status["reasons"] += [f"check failed: {name}" for name, passed in checks.items() if not passed]
    if not status["reasons"]:
        status["status"] = "REGENERATED"
    if browser:
        page = browser.page(Path(dashboard["output_files"]["dashboard_html"]), 1400)
        ids = page.evaluate("() => [...document.querySelectorAll('.js-plotly-plot')].map(e => e.id)")
        plot_data = page.evaluate("() => Object.fromEntries([...document.querySelectorAll('.js-plotly-plot')]"
                                  ".map(gd => [gd.id, JSON.parse(Plotly.Plots.graphJson(gd, false, 'keepdata'))]))")
        targets["plot_data"].write_text(json.dumps(plot_data))
        status["plot_data_ids"] = sorted(plot_data)
        # #169: each histogram panel counts every PSM, including those at its maximum
        counted = {letter: sum(sum(v for v in (trace.get("y") or []) if isinstance(v, (int, float))) for trace in plot_data[plot_id]["data"])
                   for letter, plot in FIGURE4_PLOTS.items() if plot.startswith("hist_")
                   for plot_id in [next(i for i in plot_data if i == plot or i.endswith(plot))]}
        status["histogram_counts"] = counted
        checks["every histogram panel counts all 100 PSMs (#169)"] = all(count == FIGURE4_SUBSET_ROWS for count in counted.values())
        if not checks["every histogram panel counts all 100 PSMs (#169)"]:
            status["status"] = "PROVISIONAL"
            status["reasons"].append(f"check failed: every histogram panel counts all 100 PSMs (#169): {counted}")
        for letter, plot in FIGURE4_PLOTS.items():
            match = next((i for i in ids if i == plot or i.endswith(plot) or plot in i), None)
            if not match:
                raise SystemExit(f"Figure 4 panel {letter}: plot {plot} not in the dashboard ({ids})")
            stem = panels / f"Figure4{letter}"
            page.locator(f"#{match}").screenshot(path=str(stem.with_suffix(".png")))
            png_to_pdf(stem.with_suffix(".png"), stem.with_suffix(".pdf"), dpi=288)
        page.close()
    status["files"] = {k: str(v) for k, v in dashboard["output_files"].items()}
    status["scored_csv"] = str(scored_csv)
    status["deposit"] = {key: str(path.relative_to(deposit)) for key, path in targets.items() if path.is_file()}
    return status


# --- the deposit's reference examples ----------------------------------------------------------------------------------

# Each deposited example folder: the identification table and adapter its five PSMs come from (as in PSMS).
EXAMPLE_SETS = {
    "17_maxquant_top5": {"adapter": "maxquant", "table": "processed/proteome_tools/maxquant_msms.txt"},
    "18_fragpipe_top5": {"adapter": "fragpipe", "table": "processed/arabidopsis/fragpipe_psm.tsv"},
    "19_proteome_discoverer_top5": {"adapter": "proteome_discoverer", "table": "processed/sbrodae/proteome_discoverer_PSMs.xlsx"},
    "22_hela_instanovo_top5": {"adapter": "instanovo", "table": "processed/hela/instanovo_labelled_kpreds.csv"},
    "23_sbrodae_casanovo_top5": {"adapter": "casanovo", "table": "processed/sbrodae/casanovo_predictions.csv", "by": "candidate_id"},
}
EXAMPLE_HTML = {"annotated": {"ion_types": "all", "max_isotope": 2, "neutral_losses": False}}  # as in examples/README.md


def regenerate_examples(supplement: Path, deposit: Path, candidate_cache: dict[str, list[dict[str, Any]]]) -> dict[str, Any]:
    """Re-score the deposit's reference examples: each folder's top5_summary.csv (all five PSMs), and the payloads and
    HTML views the deposit carries. The deposit holds 7 of the 25 payloads; the observed spectra of the other 18 come
    from the repository's examples/, whose payloads equal the deposited ones where both exist."""
    import pandas as pd

    from proteomics_mcp.core.visualization import (
        render_annotated_observed_spectrum_html,
        render_observed_spectrum_html,
        render_spectrum_comparison_html,
    )

    report: dict[str, Any] = {}
    for folder, spec in EXAMPLE_SETS.items():
        published = supplement / EXAMPLES / folder
        target = deposit / EXAMPLES / folder
        target.mkdir(parents=True, exist_ok=True)
        summary = pd.read_csv(published / "top5_summary.csv")
        rows, written = [], []
        for row in summary.itertuples():
            name = f"rank_{row.rank}_scan_{row.scan_number}_"
            deposited = sorted(published.glob(f"{name}*_payload.json"))
            source = deposited[0] if deposited else next((REPO_ROOT / "examples" / folder).glob(f"{name}*_payload.json"), None)
            if source is None:
                raise SystemExit(f"{folder} rank {row.rank}: no payload with the observed spectrum (deposit or examples/)")
            payload = json.loads(source.read_text())
            lookup = {"candidate_id": str(row.candidate_id)} if spec.get("by") == "candidate_id" else {"scan": int(row.scan_number)}
            candidate = resolve_candidate(supplement, {**spec, **lookup}, candidate_cache)
            if candidate.get("scan_number") is None:  # Casanovo: the scan link of the published payload
                candidate["scan_number"] = payload["metrics"]["candidate"]["scan_number"]
                candidate["spectrum_id"] = payload["metrics"]["candidate"]["spectrum_id"]
            s = payload["metrics"]["annotation_settings"]
            settings = {"tolerance": s["tolerance"], "tolerance_unit": s["tolerance_unit"], "annotation_ion_types": s["annotation_ion_types"],
                        "annotation_max_isotope": s["annotation_max_isotope"], "annotation_max_ion_charge": s["annotation_max_ion_charge"],
                        "annotation_neutral_losses": s["annotation_neutral_losses"]}
            new = score_psm(candidate, payload["observed"], settings)
            metrics, annotation_summary = new["metrics"], new["metrics"].get("annotation_summary") or {}
            rows.append({"rank": row.rank, "candidate_id": candidate.get("candidate_id"), "scan_number": candidate.get("scan_number"),
                         "spectrum_id": candidate.get("spectrum_id"), "experiment_name": candidate.get("experiment_name"),
                         "peptidoform": candidate.get("peptidoform"), "confidence": candidate.get("confidence"), "raw_score": candidate.get("raw_score"),
                         **{k: metrics.get(k) for k in ("spectral_angle", "cosine_similarity", "pearson_correlation", "xcorr", "hyperscore",
                                                        "matched_ion_count", "matched_ion_fraction")},
                         "annotated_peak_count": annotation_summary.get("annotated_peak_count"),
                         "annotated_intensity_fraction": annotation_summary.get("annotated_intensity_fraction")})
            if deposited:  # the deposit keeps its file names, so the stem (with the 0.1.10 peptide text) stays
                stem = deposited[0].name.removesuffix("_payload.json")
                (target / f"{stem}_payload.json").write_text(json.dumps(new, default=str))
                render_spectrum_comparison_html(new, target / f"{stem}_mirror.html")
                render_annotated_observed_spectrum_html(new, target / f"{stem}_annotated.html", tolerance=s["tolerance"],
                                                        tolerance_unit=s["tolerance_unit"], **EXAMPLE_HTML["annotated"])
                render_observed_spectrum_html(new, target / f"{stem}_plain.html")
                written += [f"{stem}_{kind}" for kind in ("payload.json", "mirror.html", "annotated.html", "plain.html")]
        table, columns = deposited_column_order(pd.DataFrame(rows), published / "top5_summary.csv")
        table.to_csv(target / "top5_summary.csv", index=False)
        report[folder] = {"psms": len(rows), "payload_files": written, "top5_summary_columns": columns,
                          "peptidoform_changes": {str(r["rank"]): (o, r["peptidoform"]) for r, o in zip(rows, summary["peptidoform"]) if o != r["peptidoform"]}}
    return report


# --- CHANGES.md ----------------------------------------------------------------------------------------------------------


def changes_markdown(results: dict[str, dict[str, Any]], fig2_by_version: dict[str, dict[str, Any]], ms2pip_version: str,
                     fig4: dict[str, Any], manifest: dict[str, Any]) -> str:
    specs = {key: (label, spec) for key, label, spec in STRIP_METRICS}
    old2, new2 = results["fig2"]["published"]["metrics"], results["fig2"]["new"]["metrics"]
    lines = [f"# Figure changes for proteomics-mcp 0.2.0 (release-candidate run, {manifest['git']['describe']})", "",
             f"Generated by `scripts/regenerate_figures.py` on {manifest['started']}; MS2PIP {ms2pip_version}, metric profile v2, "
             "0.05 Da. Published values: the reference payloads of Zenodo record 23295902 (proteomics-mcp 0.1.10 era, MS2PIP 4.1.2).", "",
             "## Figure 2 → image changes (panel B numbers; panel A recolouring)", "",
             "| Value | Published | New | " + " | ".join(f"MS2PIP {v}" for v in sorted(fig2_by_version)) + " | Cause |",
             "|---|---|---|" + "---|" * len(fig2_by_version) + "---|"]
    for key, (label, spec) in specs.items():
        published_key, published_label = PUBLISHED_AS.get(key, (key, None))
        old = printed(old2.get(published_key), "{:.3f}" if published_label else spec) + (f" ({published_label})" if published_label else "")
        new = printed(new2.get(key), spec)
        others = [printed(fig2_by_version[v]["metrics"].get(key), spec) for v in sorted(fig2_by_version)]
        changed = old != new or any(o != new for o in others)
        lines.append(f"| {label} | {old} | {'**' + new + '**' if old != new else new} | " + " | ".join(others) + f" | {CAUSES[key] if changed else '-'} |")
    old_counts = (old2.get("annotation_summary") or {}).get("ion_type_counts") or {}
    new_counts = (new2.get("annotation_summary") or {}).get("ion_type_counts") or {}
    diff = {k: (old_counts.get(k), new_counts.get(k)) for k in sorted(set(old_counts) | set(new_counts)) if old_counts.get(k) != new_counts.get(k)}
    lines += ["", f"Observed ion series tile: " + ", ".join(f"{k}:{o} → {n}" for k, (o, n) in diff.items()) + " (#79, PR #104: w-ion masses).", ""]
    removed = sorted(label for _, label in mirror_labels(results["fig2"]["published"]) - mirror_labels(results["fig2"]["new"]))
    lines += [f"Panel A: matched peaks {', '.join(removed)} are no longer matched (grey instead of blue); the 14 labelled peaks are "
              "unchanged. Their predicted ions fall below the 0.1% cut-off (#83, PR #106). Unchanged: confidence 100.0%, raw score "
              "−2.80e-04, next beam NA.", ""]
    lines += ["## Figure 3", "", "**3A → unchanged.** Same peptidoform (`DM[UNIMOD:35]HGVTSHLSNQELQDLVEFMK`), same 237 peak labels as 0.1.10. "
              "Rendering note: the sequence stubs read `DM(ox)HGVTSH…`; the published panel shows `DMHGVTSH…` (pre-release renderer; "
              "every release since 0.1.8 shows the modification). Keep the published 3A or accept the stub notation.", ""]
    b_old, b_new = results["fig3b"]["published"], results["fig3b"]["new"]
    ann_old = b_old["metrics"]["annotation_summary"]  # the published (0.1.10) labels, stored with the panel's settings
    ann_new = annotation(b_new, FIG3B_ION_TYPES)
    changed_labels = sum(1 for a, b in zip(ann_old["annotations"], ann_new["annotations"]) if a != b)
    lines += [f"**3B → image changes (published image wrong).** Peptidoform `{b_old['metrics']['candidate']['peptidoform']}` → "
              f"`{b_new['metrics']['candidate']['peptidoform']}`: the MaxQuant parameters fix carbamidomethyl on C (#101, PR #112). "
              f"Labelled peaks {ann_old['annotated_peak_count']} → {ann_new['annotated_peak_count']}; {changed_labels} peak labels change. "
              f"Legend: {', '.join(sorted(set(ann_old['ion_type_counts']) - set(ann_new['ion_type_counts'])))} removed, "
              f"{', '.join(sorted(set(ann_new['ion_type_counts']) - set(ann_old['ion_type_counts']))) or 'nothing'} added. Precursor error "
              f"{printed(b_old['metrics'].get('precursor_mass_error_ppm'), '{:.1f}')} ppm → {printed(b_new['metrics'].get('precursor_mass_error_ppm'), '{:.3f}')} ppm "
              f"(not printed in the figure). w-ion labels also follow #79 (PR #104).", ""]
    lines += ["## Supplementary Figure 2 → image changes (panel B only)", "",
              f"Panel B: {', '.join(removed)} are unmatched under v2, so the `b2+ | EF` and `b8+ | EFHLNESG` labels disappear and those "
              "peaks turn grey (#83, PR #106). Panels A and C are unchanged; no caption values change.", "",
              "## Supplementary Figure 5 → image changes (panels B and C)", ""]
    for key, letter in (("s5a", "A"), ("s5b", "B"), ("s5c", "C"), ("s5d", "D")):
        removed5 = sorted(label for _, label in mirror_labels(results[key]["published"]) - mirror_labels(results[key]["new"]))
        lines.append(f"- {letter} ({PSMS[key]['title']}): " + (f"no longer matched: {', '.join(removed5)} (#74, #83)" if removed5 else "unchanged"))
    lines += ["", "No values are printed in Supplementary Figure 5.", "", f"## Figure 4 → {fig4.get('status')}", ""]
    lines += [f"- {note}" for note in fig4.get("notes", [])]
    lines += [f"- **{reason}**" for reason in fig4.get("reasons", [])]
    lines += [f"- check: {name}: {'passed' if passed else 'FAILED'}" for name, passed in (fig4.get("checks") or {}).items()]
    if fig4.get("batch"):
        lines += ["", f"Built from {fig4['spectra']['source']} (MGF sha256 `{fig4['spectra']['mgf_sha256']}`); "
                  f"{fig4['batch']['rows_with_spectral_angle']} of {fig4['batch']['dashboard_rows']} dashboard rows scored; "
                  f"{fig4['batch']['rows_with_rt']} with retention time."]
        issue_137 = fig4.get("issue_137") or {}
        lines += ["", f"#137 check: {issue_137.get('compared')} PSMs re-scored with `score_spectrum` (the package predicts); "
                  f"{issue_137.get('differing_above_1e-6')} differ by more than 1e-6 from the batch scores, which score the batch's stored "
                  f"predictions (largest difference {issue_137.get('largest_abs_difference', 0):.1e})."]
        rt = fig4.get("retention_time_check") or {}
        join = fig4["batch"]["metadata_join"] or {}
        lines += ["", f"#146 check: the metadata join matched {join.get('matched_rows')} rows, {join.get('unmatched_rows')} unmatched; "
                  f"dashboard retention times vs FragPipe `Retention`: largest difference "
                  f"{rt.get('largest_abs_difference_minutes') or float('nan'):.6f} min over {rt.get('compared')} PSMs."]
        for key in ("scored_psms_columns", "scored_peptides_columns"):
            columns = fig4.get(key) or {}
            lines += ["", f"`{FIGURE4_DEPOSIT[key.removesuffix('_columns')].split('/')[-1]}` columns: removed {columns.get('removed') or 'none'}, "
                      f"added {columns.get('added') or 'none'}."]
    lines += ["", "## Checks", "",
              f"- Verification against the expected 0.2.0 values: {'passed' if manifest['verification']['passed'] else 'FAILED: ' + '; '.join(manifest['verification']['problems'])}.",
              "- #137 (supplied predictions skip the fragment-charge rule): every spectrum panel is scored by `score_spectrum` with the "
              "package predicting; the same PSMs scored with that prediction passed in give identical metrics (checked for all 7).",
              "- Figure 3A: its 237 peak labels equal those of 0.1.10, the version of the published figure (checked in this run "
              "against `FIG3A_LABELS_SHA256`)."]
    lines += ["", "## Manuscript text and captions", "",
              f"- **§2.5 Matching and metrics:** add that each observed peak explains at most one predicted ion (closest first), and "
              "that predicted ions below 0.1% of the most intense predicted ion are not scored (metric profile v2, the default); "
              "spectral angle, cosine and Pearson compare every scored predicted ion, unmatched ones with zero observed intensity.",
              f"- **§3.2:** spectral angle 0.833 → {printed(new2['spectral_angle'], '{:.3f}')}, cosine {printed(new2['cosine_similarity'], '{:.3f}')}, "
              f"Pearson 0.962 → {printed(new2['pearson_correlation'], '{:.3f}')}, \"22 of 24 predicted b/y ions\" → "
              f"\"{new2['matched_ion_count']} of {new2['total_predicted_ions']}\", 91.7% → {printed(new2['matched_ion_fraction'], 'pct')} matched fraction, "
              f"48.0% → {printed(new2['explained_intensity_fraction'], 'pct')} explained intensity.",
              f"- **Figure 2 legend:** every changed value in the table above (spectral angle, Pearson, matched ions, matched fraction, "
              "explained intensity, fragment ppm and Da, precursor ppm). The strip now shows Comet's XCorr "
              f"({printed(new2.get('xcorr'), '{:.3f}')}) and X! Tandem's hyperscore ({printed(new2.get('hyperscore'), '{:.1f}')}) instead of "
              f"log2 XCorr and log2 hyperscore (#51); name them so, and state MS2PIP {FIGURE_MS2PIP} (#23).",
              "- **§3.3 and Figure 3 legend:** write the peptide as `AQNSC[UNIMOD:4]PIC[UNIMOD:4]QK` (carbamidomethylated cysteines, "
              f"from the MaxQuant fixed modifications); ion vocabulary in panel B: {', '.join(sorted(ann_new['ion_type_counts']))} "
              "(no a+1 or d-alt; precursor ions now labelled).",
              "- **Data availability:** cite proteomics-mcp 0.2.0 and the new Zenodo version.", ""]
    return "\n".join(lines)


# --- main --------------------------------------------------------------------------------------------------------------


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", type=Path, help="output folder (outside the repository)")
    parser.add_argument("--cache-dir", type=Path, required=True, help="downloads: the Zenodo supplement and PRIDE USI spectra")
    parser.add_argument("--zenodo-zip", type=Path, help="a local copy of the Zenodo supplement.zip (md5-checked); downloaded otherwise")
    parser.add_argument("--compare-ms2pip", default="4.2.0", help="other MS2PIP releases for the Figure 2 side-by-side ('' to skip)")
    parser.add_argument("--python", default="3.12", help="Python for the figure and MS2PIP comparison environments (4.1.2 has no 3.14 wheels)")
    parser.add_argument("--browser-channel", default="chrome", help="Playwright browser: chrome (local Google Chrome) or chromium")
    parser.add_argument("--skip-figure4", action="store_true", help="do not fetch PRIDE spectra or build Figure 4")
    parser.add_argument("--figure4-mgf", type=Path, help="build Figure 4 from this MGF (the deposit's figure4_usi_spectra.mgf) instead of PRIDE")
    parser.add_argument("--usi-timeout", type=float, default=120.0, help="seconds per PRIDE USI request")
    parser.add_argument("--usi-retries", type=int, default=4)
    parser.add_argument("--no-verify", action="store_true", help="assemble even if the expected 0.2.0 values differ")
    parser.add_argument("--figure2-metrics-only", type=Path, help=argparse.SUPPRESS)
    parser.add_argument("--in-figure-environment", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    if not args.out and not args.figure2_metrics_only:
        parser.error("--out is required")
    return args


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    args = parse_args(argv)
    if not args.figure2_metrics_only and not args.in_figure_environment:
        return rerun_in_figure_environment(args, argv)
    args.cache_dir = args.cache_dir.expanduser().resolve()
    args.cache_dir.mkdir(parents=True, exist_ok=True)
    archive, supplement = zenodo_supplement(args.zenodo_zip and args.zenodo_zip.expanduser().resolve(), args.cache_dir)
    examples = supplement / EXAMPLES
    candidate_cache: dict[str, list[dict[str, Any]]] = {}

    def published_payload(key: str) -> dict[str, Any]:
        return json.loads((examples / PSMS[key]["payload"]).read_text())

    def new_payload(key: str) -> dict[str, Any]:
        published = published_payload(key)
        candidate = resolve_candidate(supplement, PSMS[key], candidate_cache)
        candidate.setdefault("scan_number", None)
        if candidate["scan_number"] is None:  # Casanovo: the published scan link (from the Proteome Discoverer export)
            candidate["scan_number"] = published["metrics"]["candidate"]["scan_number"]
            candidate["spectrum_id"] = published["metrics"]["candidate"]["spectrum_id"]
        s = published["metrics"]["annotation_settings"]
        settings = {"tolerance": s["tolerance"], "tolerance_unit": s["tolerance_unit"], "annotation_ion_types": s["annotation_ion_types"],
                    "annotation_max_isotope": s["annotation_max_isotope"], "annotation_max_ion_charge": s["annotation_max_ion_charge"],
                    "annotation_neutral_losses": s["annotation_neutral_losses"]}
        return score_psm(candidate, published["observed"], settings)

    if args.figure2_metrics_only:
        payload = new_payload("fig2")
        args.figure2_metrics_only.write_text(json.dumps({"metrics": payload["metrics"], "versions": package_versions(),
                                                         "python": platform.python_version()}, default=str))
        return 0

    out = args.out.expanduser().resolve()
    if REPO_ROOT in out.parents or out == REPO_ROOT:
        raise SystemExit("--out must be outside the repository")
    panels, comparison, work, deposit = out / "panels", out / "comparison", out / "work", out / "deposit"
    for folder in (out, panels, comparison, work, deposit):
        folder.mkdir(parents=True, exist_ok=True)
    started = datetime.now(timezone.utc).isoformat(timespec="seconds")
    ms2pip_version = version("ms2pip")
    log(f"proteomics-mcp {version('proteomics-mcp')} ({git_state()['describe']}), MS2PIP {ms2pip_version}, Python {platform.python_version()}")

    results: dict[str, dict[str, Any]] = {}
    for key in PSMS:
        log(f"scoring {key}")
        results[key] = {"published": published_payload(key), "new": new_payload(key)}
        (work / f"{key}_payload.json").write_text(json.dumps(results[key]["new"], default=str))

    fig2_by_version: dict[str, dict[str, Any]] = {}
    for release in [v.strip() for v in args.compare_ms2pip.split(",") if v.strip()]:
        fig2_by_version[release] = ({"metrics": results["fig2"]["new"]["metrics"], "versions": package_versions(), "python": platform.python_version(),
                                     "command": "this run"} if release == ms2pip_version else figure2_metrics_with(release, args, work))
    problems = verify(results, ms2pip_version, {v: r["metrics"] for v, r in fig2_by_version.items()})
    (work / "verification.json").write_text(json.dumps({"problems": problems}, indent=1))
    if problems and not args.no_verify:
        for problem in problems:
            log(f"MISMATCH {problem}")
        log("verification failed: nothing assembled (see work/verification.json)")
        return 3
    log("verification passed" if not problems else f"verification: {len(problems)} mismatches (--no-verify)")
    log("re-scoring the deposit's reference examples")
    examples_report = regenerate_examples(supplement, deposit, candidate_cache)

    from proteomics_mcp.core.visualization import render_spectrum_comparison_html, render_static_spectrum

    def mirror(key: str, stem: str, title: str | None = None) -> Path:
        for suffix in (".pdf", ".png"):
            render_static_spectrum(results[key]["new"], panels / f"{stem}{suffix}", view="mirror", label_top_n=MIRROR_LABELS, title=title)
        return panels / f"{stem}.pdf"

    def annotated(key: str, stem: str, ion_types: str, title: str) -> Path:
        for suffix in (".pdf", ".png"):
            render_static_spectrum(results[key]["new"], panels / f"{stem}{suffix}", view="annotated", ion_types=ion_types, max_isotope=2,
                                   neutral_losses=True, include_fragment_sequences=True, label_top_n=ANNOTATED_LABELS, title=title)
        return panels / f"{stem}.pdf"

    log("rendering spectrum panels")
    fig2a = mirror("fig2", "Figure2A")
    fig3a = annotated("fig3a", "Figure3A", FIG3A_ION_TYPES, "Modified peptidoform annotation")
    fig3b = annotated("fig3b", "Figure3B", FIG3B_ION_TYPES, "Observed spectrum with extended ion series")
    s5 = {letter: mirror(key, f"Suppl_figure5{letter}", PSMS[key]["title"]) for key, letter in (("s5a", "A"), ("s5b", "B"), ("s5c", "C"), ("s5d", "D"))}
    html = work / "fig2_mirror.html"
    render_spectrum_comparison_html(results["fig2"]["new"], html)

    browser = Browser(args.browser_channel)
    try:
        log("capturing HTML panels")
        page = browser.page(html, METRIC_STRIP_VIEWPORT)
        isolate_and_capture(page, "() => [...document.querySelectorAll('section.metric-group')]", panels / "Figure2B")
        page.close()
        page = browser.page(html, METADATA_VIEWPORT)
        isolate_and_capture(page, "() => [[...document.querySelectorAll('h2')].find(h => h.textContent.trim() === 'Prediction Metadata').parentElement]",
                            panels / "Suppl_figure2A")
        page.close()
        page = browser.page(html, 1950)
        page.locator("input[type=checkbox]").first.check()
        capture_plot(page, panels / "Suppl_figure2B", relayout={"xaxis.range": list(S2B_MZ_RANGE), "showlegend": False}, size=(1900, 655))
        page.close()
        page = browser.page(html, METRIC_STRIP_VIEWPORT)
        page.locator("input[type=checkbox]").first.check()
        # Supplementary Figure 2C is a zoom without title, y axis or legend, with the hover card of y4+ | PSSK
        capture_plot(page, panels / "Suppl_figure2C", relayout={"xaxis.range": list(S2C_MZ_RANGE), "yaxis.autorange": True, "yaxis.showticklabels": False, "yaxis.title.text": "",
                                                                "title.text": "", "showlegend": False, "margin": {"t": 30, "l": 30, "r": 30, "b": 70}},
                     size=(964, 518), hover_mz=S2C_HOVER_MZ)
        page.close()
        fig4 = {"status": "SKIPPED", "reasons": ["--skip-figure4"]} if args.skip_figure4 else figure4(supplement, args, work, panels, browser, deposit)
    finally:
        browser.close()

    log("assembling figures")
    widths: dict[str, float] = {}
    widths |= compose([[("A", fig2a, 1.0)], [("B", panels / "Figure2B.pdf", 1.0)]], PAGE_WIDTH["Figure2"], out / "Figure2.pdf")
    compose([[("A", fig3a, 1.0)], [("B", fig3b, 1.0)]], PAGE_WIDTH["Figure3"], out / "Figure3.pdf")
    widths |= compose([[("A", panels / "Suppl_figure2A.pdf", 1.0)], [("B", panels / "Suppl_figure2B.pdf", 1.0)], [("C", panels / "Suppl_figure2C.pdf", 0.66)]],
            PAGE_WIDTH["Suppl_figure2"], out / "Suppl_figure2.pdf")
    widths |= compose([[("A", s5["A"], 0.5), ("C", s5["C"], 0.5)], [("B", s5["B"], 0.5), ("D", s5["D"], 0.5)]], PAGE_WIDTH["Suppl_figure5"], out / "Suppl_figure5.pdf")
    names = ["Figure2", "Figure3", "Suppl_figure2", "Suppl_figure5"]
    if fig4["status"] != "SKIPPED":
        p = {letter: panels / f"Figure4{letter}.pdf" for letter in FIGURE4_PLOTS}
        widths |= compose([[("A", p["A"], 0.5), ("B", p["B"], 0.5)], [("C", p["C"], 1.0)], [("D", p["D"], 0.5), ("E", p["E"], 0.5)]], PAGE_WIDTH["Figure4"],
                out / "Figure4.pdf")
        names.append("Figure4")
    dpi = raster_dpi(panels, widths)
    low = {stem: value for stem, value in dpi.items() if value < MIN_RASTER_DPI}
    if low:
        problems.append(f"screenshot panels below {MIN_RASTER_DPI} dpi: {low}")
        log(f"MISMATCH screenshot panels below {MIN_RASTER_DPI} dpi: {low}")
    (deposit / PUBLISHED_PDFS).mkdir(parents=True, exist_ok=True)
    for name in names:
        shutil.copyfile(out / f"{name}.pdf", deposit / PUBLISHED_PDFS / f"{name}.pdf")
    git = git_state()
    label = f"Regenerated ({git['describe']}, MS2PIP {ms2pip_version}, v2)" + (" PROVISIONAL" if fig4["status"] == "PROVISIONAL" else "")
    for name in names:
        side_by_side(supplement / PUBLISHED_PDFS / f"{name}.pdf", out / f"{name}.png", comparison / f"{name}_published_vs_new.png",
                     label if name == "Figure4" else label.replace(" PROVISIONAL", ""))

    inputs = {"zenodo_supplement": {"path": str(archive), "md5": file_digest(archive, "md5"), "sha256": file_digest(archive), "record": ZENODO_RECORD}}
    for key, spec in PSMS.items():
        for rel in (spec["table"], f"{EXAMPLES}/{spec['payload']}"):
            inputs[rel] = {"sha256": file_digest(supplement / rel)}
    manifest = {
        "started": started, "finished": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "command": user_command(), "interpreter": sys.executable, "cwd": os.getcwd(),
        "package": {"proteomics-mcp": version("proteomics-mcp"), "source": str(Path(__import__("proteomics_mcp").__file__).parent)},
        "git": git, "ms2pip": ms2pip_version, "python": platform.python_version(), "platform": platform.platform(), "versions": package_versions(),
        "figure_environment": {"pinned": FIGURE_ENVIRONMENT, "requirements": "work/requirements_figures.txt",
                               "requirements_sha256": file_digest(work / "requirements_figures.txt") if (work / "requirements_figures.txt").is_file() else None},
        "settings": {"metric_profile": "v2 (package default)", "tolerance": "0.05 Da", "predictor": PREDICTOR, "mirror_label_top_n": MIRROR_LABELS,
                     "annotated_label_top_n": ANNOTATED_LABELS, "figure3a_ion_types": FIG3A_ION_TYPES, "figure3b_ion_types": FIG3B_ION_TYPES},
        "psms": {key: {"table": spec["table"], "adapter": spec["adapter"], "peptidoform": results[key]["new"]["metrics"]["candidate"]["peptidoform"],
                       "scan_number": results[key]["new"]["metrics"]["candidate"]["scan_number"],
                       "run": results[key]["new"]["metrics"]["candidate"].get("experiment_name"),
                       "observed_spectrum": f"{EXAMPLES}/{spec['payload']}",
                       "warnings": results[key]["new"]["metrics"]["candidate"].get("warnings")} for key, spec in PSMS.items()},
        "inputs": inputs, "figure2_ms2pip": {v: {"metrics": {k: r["metrics"].get(k) for k, _, _ in STRIP_METRICS}, "versions": r.get("versions"),
                                                "python": r.get("python"), "command": r.get("command")} for v, r in fig2_by_version.items()},
        "verification": {"problems": problems, "passed": not problems, "figure4": fig4.get("status")}, "figure4": fig4,
        "raster_panels_dpi": dpi,
        "examples": examples_report,
    }
    # The environment of this run, for the deposit: uv.lock of this checkout, the figure environment's requirements
    # (uv.lock with MS2PIP and the assembly tools pinned) and every installed distribution.
    regeneration = deposit / "regeneration"
    regeneration.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(REPO_ROOT / "uv.lock", regeneration / "uv.lock")
    if (work / "requirements_figures.txt").is_file():
        shutil.copyfile(work / "requirements_figures.txt", regeneration / "requirements_figures.txt")
    from importlib.metadata import distributions

    installed = sorted({f"{d.metadata['Name']}=={d.version}" for d in distributions()}, key=str.lower)
    (regeneration / "environment.txt").write_text(f"# Python {platform.python_version()}, {platform.platform()}\n" + "\n".join(installed) + "\n")
    manifest["outputs"] = {str(p.relative_to(out)): file_digest(p) for p in sorted(out.rglob("*")) if p.is_file() and "work" not in p.relative_to(out).parts}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=1, default=str))
    (out / "CHANGES.md").write_text(changes_markdown(results, fig2_by_version, ms2pip_version, fig4, manifest))
    log(f"done: {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
