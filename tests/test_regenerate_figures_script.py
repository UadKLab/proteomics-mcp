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
