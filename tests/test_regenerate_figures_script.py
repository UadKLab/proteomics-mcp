"""scripts/regenerate_figures.py: the helpers and the Figure 3A label reference, without network or browser."""

import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

from proteomics_mcp.core.annotation import annotate_observed_ions
from proteomics_mcp.models import ObservedSpectrum

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "regenerate_figures.py"
FIG3A_PAYLOAD = ROOT / "examples" / "19_proteome_discoverer_top5" / "rank_1_scan_24381_DM_UNIMOD_35_HGVTSHLSNQELQDLVEFMK_payload.json"

pytestmark = pytest.mark.skipif(not SCRIPT.is_file(), reason="scripts/ is only available in a repository checkout")


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("regenerate_figures", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_printed_values_use_the_figure_precision(script):
    assert script.printed(0.83406, "{:.3f}") == "0.834"
    assert script.printed(1.0, "pct") == "100.0%"
    assert script.printed(18.0, "{:d}") == "18"
    assert script.printed(None, "pct") == "NA"


def test_out_is_required_unless_only_figure2_metrics(script):
    with pytest.raises(SystemExit):
        script.parse_args(["--cache-dir", "cache"])
    args = script.parse_args(["--cache-dir", "cache", "--out", "figures"])
    assert args.compare_ms2pip == "4.2.0" and not args.skip_figure4 and not args.in_figure_environment


def test_figures_pin_ms2pip_412_in_their_own_environment_only(script):
    """#23: the figures use MS2PIP 4.1.2; the pin lives in the figure environment, not in the package dependencies."""
    assert script.FIGURE_MS2PIP == "4.1.2" and "ms2pip==4.1.2" in script.FIGURE_ENVIRONMENT
    assert script.EXPECTED["fig2"]["fragment_mass_error_abs_mean_ppm"] == "2.796"
    pyproject = (ROOT / "pyproject.toml").read_text()
    assert "ms2pip==" not in pyproject


def test_expected_values_cover_every_supplementary_figure5_panel(script):
    assert set(script.EXPECTED["mirror_labels_removed"]) == {key for key in script.PSMS if key.startswith("s5")}


@pytest.mark.skipif(not FIG3A_PAYLOAD.is_file(), reason="examples/ is only available in a repository checkout")
def test_figure3a_labels_match_the_published_figure(script):
    """Figure 3A must come out identical: the package's labels for it equal those of 0.1.10 (the published figure)."""
    payload = json.loads(FIG3A_PAYLOAD.read_text())
    labels = annotate_observed_ions(ObservedSpectrum(**payload["observed"]), peptidoform=payload["metrics"]["candidate"]["peptidoform"],
                                    tolerance=0.05, ion_types=script.FIG3A_ION_TYPES, max_isotope=2, neutral_losses=True)
    assert labels["annotated_peak_count"] == script.EXPECTED["fig3a"]["labelled_peaks"]
    assert hashlib.sha256(json.dumps(labels["annotations"]).encode()).hexdigest() == script.FIG3A_LABELS_SHA256


def test_figure2_strip_shows_the_engine_scores(script):
    """#51: the strip shows Comet XCorr and the X! Tandem hyperscore, compared in CHANGES.md with the 0.1.x variants."""
    keys = [key for key, _, _ in script.STRIP_METRICS]
    assert "xcorr" in keys and "hyperscore" in keys and not any(key.startswith("log2_") for key in keys)
    assert set(script.PUBLISHED_AS) == {"xcorr", "hyperscore"} and all(key in script.CAUSES for key in keys)
    assert (script.EXPECTED["fig2"]["xcorr"], script.EXPECTED["fig2"]["hyperscore"]) == ("3.336", "59.5")


def test_figure4_mgf_carries_the_retention_time(script, tmp_path):
    """#135: the Figure 4 MGF takes each PRIDE spectrum's rt_seconds at full precision; spectra without one are reported."""
    spectra = {7: {"mz": [100.0, 200.5], "intensity": [10.0, 20.0], "precursor_mz": 500.25, "precursor_charge": 2, "rt_seconds": 699.59159},
               8: {"mz": [150.0], "intensity": [5.0], "precursor_mz": 600.0, "precursor_charge": None, "rt_seconds": None}}
    mgf = tmp_path / "f4.mgf"
    assert script.write_figure4_mgf([(7, 2), (8, 3)], spectra, mgf) == [8]
    text = mgf.read_text()
    assert "SCANS=7\nRTINSECONDS=699.59159\n" in text and "CHARGE=3+" in text and text.count("BEGIN IONS") == 2
    assert "\n200.5 20.0\n" in text  # full precision: what PRIDE served, not rounded


def test_regenerated_tables_keep_the_deposited_column_order(script, tmp_path):
    import pandas as pd

    deposited = tmp_path / "old.csv"
    pd.DataFrame(columns=["candidate_id", "log2_xcorr", "spectral_angle"]).to_csv(deposited, index=False)
    new = pd.DataFrame({"spectral_angle": [0.5], "xcorr": [1.2], "candidate_id": ["psm:0"]})
    ordered, columns = script.deposited_column_order(new, deposited)
    assert list(ordered.columns) == ["candidate_id", "spectral_angle", "xcorr"]
    assert columns == {"removed": ["log2_xcorr"], "added": ["xcorr"]}


def test_example_sets_cover_the_figure_payloads(script):
    """Every figure PSM's payload folder is re-scored for the deposit, from the same table and adapter."""
    for spec in script.PSMS.values():
        folder = spec["payload"].split("/")[0]
        assert script.EXAMPLE_SETS[folder]["table"] == spec["table"] and script.EXAMPLE_SETS[folder]["adapter"] == spec["adapter"]


def test_raster_panels_are_the_screenshots(script):
    """Only the screenshots of interactive plots are raster panels, and they must reach 300 dpi in the figure."""
    assert script.MIN_RASTER_DPI == 300
    assert script.RASTER_PANELS == {"Suppl_figure2B", "Suppl_figure2C", "Figure4A", "Figure4B", "Figure4C", "Figure4D", "Figure4E",
                                    "Suppl_figure3", "Suppl_figure4"}


def test_supplementary_tables_fill_the_published_page_above_300_dpi(script):
    """Supplementary Figures 3 and 4: the viewport has the published page's aspect, and at device scale 3 it is far above 300 dpi."""
    width, height = script.TABLE_VIEWPORT
    for name in script.SUPPLEMENTARY_TABLES:
        page_width = script.PAGE_WIDTH[name]
        assert abs(width / height - page_width / 361.41) < 0.01
        assert width * 3 / (page_width / 72) >= script.MIN_RASTER_DPI
    assert {key for key, _ in script.SUPPLEMENTARY_TABLES.values()} <= set(script.FIGURE4_DEPOSIT)


def test_table_header_reads_the_row_and_column_counts(script, tmp_path):
    html = tmp_path / "table.html"
    html.write_text("<h1>Peptide table</h1><p class='meta'>68 rows, 61 columns</p>")
    assert script.table_header(html) == "68 rows, 61 columns"
    html.write_text("<h1>empty</h1>")
    with pytest.raises(SystemExit):
        script.table_header(html)


def test_source_of_a_tag_must_equal_the_installed_source(script, monkeypatch, capsys):
    monkeypatch.setattr(script, "tag_source_digest", lambda tag: "a" * 64)
    assert script.check_source_of("v0.2.0", "a" * 64) == {"tag": "v0.2.0", "sha256": "a" * 64, "equal": True}
    with pytest.raises(SystemExit) as stopped:
        script.check_source_of("v0.2.0", "b" * 64)
    assert stopped.value.code == 5
    assert "not that of v0.2.0" in capsys.readouterr().err


def test_manifest_records_the_command_as_typed(script, monkeypatch):
    monkeypatch.setattr(script.sys, "argv", [str(script.REPO_ROOT / "scripts" / "regenerate_figures.py"), "--out", "/tmp/f", "--in-figure-environment"])
    assert script.user_command() == "uv run python scripts/regenerate_figures.py --out /tmp/f"
    assert "\\" not in script.user_command()  # forward slashes on Windows too
    monkeypatch.setattr(script.sys, "argv", [script.sys.argv[0], "--out", r"C:\pmcp\figures"])
    assert script.user_command().endswith("--out C:/pmcp/figures")


def test_a_stale_installed_package_stops_the_run_before_scoring(script, tmp_path, capsys):
    """The figures were once made with uv's cached build of an earlier commit: the installed code must be the checkout's."""
    checkout = tmp_path / "src" / "proteomics_mcp"
    installed = tmp_path / "site-packages" / "proteomics_mcp"
    for package in (checkout, installed):
        (package / "core").mkdir(parents=True)
        (package / "__init__.py").write_text("")
        (package / "core" / "scoring.py").write_text("PRECURSOR_ISOTOPE_ERRORS = (0, 1, 2, 3)\n")
    (installed / "__pycache__").mkdir()
    (installed / "__pycache__" / "x.cpython-312.pyc").write_bytes(b"ignored")
    assert script.check_installed_package(installed, checkout) == script.source_digest(checkout)

    (installed / "core" / "scoring.py").write_text("PRECURSOR_ISOTOPE_ERRORS = (0, 1, 2)\n")
    with pytest.raises(SystemExit) as stopped:
        script.check_installed_package(installed, checkout)
    assert stopped.value.code == 5
    assert "not this checkout's code" in capsys.readouterr().err

    (installed / "core" / "scoring.py").write_text("PRECURSOR_ISOTOPE_ERRORS = (0, 1, 2, 3)\n")
    (installed / "core" / "extra.py").write_text("")
    with pytest.raises(SystemExit):
        script.check_installed_package(installed, checkout)


def test_the_installed_package_is_checked_against_the_checkout(script):
    """In this test environment the package is the checkout's, so the guard passes."""
    import proteomics_mcp

    installed = Path(proteomics_mcp.__file__).parent
    if installed != ROOT / "src" / "proteomics_mcp" and script.source_digest(installed) != script.source_digest(ROOT / "src" / "proteomics_mcp"):
        pytest.skip("the installed package is not this checkout (e.g. a wheel from another commit)")
    assert script.check_installed_package() == script.source_digest(ROOT / "src" / "proteomics_mcp")


def test_a_missing_zenodo_zip_is_reported_not_downloaded(script, tmp_path, monkeypatch):
    monkeypatch.setattr(script.urllib.request, "urlopen", lambda *a, **k: pytest.fail("downloaded"))
    with pytest.raises(SystemExit, match="no such file"):
        script.zenodo_supplement(tmp_path / "supplment.zip", tmp_path / "cache")


def test_an_interrupted_extraction_is_redone(script, tmp_path, monkeypatch):
    import zipfile

    archive = tmp_path / "supplement.zip"
    with zipfile.ZipFile(archive, "w") as bundle:
        bundle.writestr("supplement/raw_file_inventory.tsv", "run\n")
        bundle.writestr("supplement/figure_assets/examples/payload.json", "{}")
    monkeypatch.setattr(script, "ZENODO_MD5", hashlib.md5(archive.read_bytes()).hexdigest())
    cache = tmp_path / "cache"
    # what an interrupted extractall leaves: the inventory, not the later members, and no marker
    partial = cache / f"zenodo-{script.ZENODO_RECORD}" / "supplement"
    partial.mkdir(parents=True)
    (partial / "raw_file_inventory.tsv").write_text("run\n")

    _, supplement = script.zenodo_supplement(archive, cache)

    assert (supplement / "figure_assets" / "examples" / "payload.json").is_file()


def test_out_must_be_new_or_empty(script, tmp_path, monkeypatch):
    monkeypatch.setattr(script, "rerun_in_figure_environment", lambda *a: pytest.fail("ran"))
    (tmp_path / "deposit").mkdir()
    (tmp_path / "deposit" / "Figure4.pdf").write_text("an earlier run")
    with pytest.raises(SystemExit, match="not empty"):
        script.main(["--out", str(tmp_path), "--cache-dir", str(tmp_path / "cache")])


def test_a_missing_dashboard_plot_is_named(script):
    assert script.dashboard_plot(["plot_hist_rt_minutes"], "B", "hist_rt") == "plot_hist_rt_minutes"
    with pytest.raises(SystemExit, match="plot hist_mz not in the dashboard"):
        script.dashboard_plot(["plot_hist_rt"], "C", "hist_mz")


def test_a_scan_number_in_two_runs_is_refused(script):
    table = "processed/table.tsv"
    rows = [{"scan_number": 7, "peptidoform": "PEPTIDE", "experiment_name": run} for run in ("run_a", "run_b")]
    cache = {str(Path("supplement") / table): rows}
    with pytest.raises(SystemExit, match="runs"):
        script.resolve_candidate(Path("supplement"), {"adapter": "generic", "table": table, "scan": 7}, cache)
    assert script.resolve_candidate(Path("supplement"), {"adapter": "generic", "table": table, "scan": 7}, {str(Path("supplement") / table): rows[:1]})["experiment_name"] == "run_a"


def test_deposited_examples_are_rescored_and_rendered(script, tmp_path, monkeypatch):
    """The deposit's own payloads get new payload, mirror, annotated and plain files, the annotated view with the payload's
    published tolerance (a NameError here once stopped the run after verification)."""
    folder = "22_hela_instanovo_top5"
    published = tmp_path / "supplement" / script.EXAMPLES / folder
    published.mkdir(parents=True)
    source = ROOT / "examples" / folder
    (published / "top5_summary.csv").write_bytes((source / "top5_summary.csv").read_bytes())
    deposited = next(source.glob("rank_2_scan_11648_*_payload.json"))
    (published / deposited.name).write_bytes(deposited.read_bytes())
    payload = json.loads(deposited.read_text())

    monkeypatch.setattr(script, "EXAMPLE_SETS", {folder: {"adapter": "instanovo", "table": "unused.csv"}})
    monkeypatch.setattr(script, "resolve_candidate", lambda supplement, spec, cache: dict(payload["metrics"]["candidate"]))
    monkeypatch.setattr(script, "score_psm", lambda candidate, observed, settings: {**payload, "settings_used": settings})
    tolerances = []
    import proteomics_mcp.core.visualization as visualization

    render = visualization.render_annotated_observed_spectrum_html
    monkeypatch.setattr(visualization, "render_annotated_observed_spectrum_html",
                        lambda *args, **kwargs: tolerances.append((kwargs["tolerance"], kwargs["tolerance_unit"])) or render(*args, **kwargs))

    report = script.regenerate_examples(tmp_path / "supplement", tmp_path / "deposit", {})

    stem = deposited.name.removesuffix("_payload.json")
    assert report[folder]["payload_files"] == [f"{stem}_{kind}" for kind in ("payload.json", "mirror.html", "annotated.html", "plain.html")]
    assert all((tmp_path / "deposit" / script.EXAMPLES / folder / name).is_file() for name in report[folder]["payload_files"])
    settings = payload["metrics"]["annotation_settings"]
    assert tolerances == [(settings["tolerance"], settings["tolerance_unit"])]
    assert report[folder]["psms"] == 5
