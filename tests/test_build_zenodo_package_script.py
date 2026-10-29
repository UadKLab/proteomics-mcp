"""scripts/build_zenodo_package.py on a small fake deposit, without network."""

import hashlib
import importlib.util
import json
import subprocess
import zipfile
from pathlib import Path, PureWindowsPath

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "build_zenodo_package.py"

pytestmark = pytest.mark.skipif(not SCRIPT.is_file(), reason="scripts/ is only available in a repository checkout")


@pytest.fixture(scope="module")
def script():
    spec = importlib.util.spec_from_file_location("build_zenodo_package", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SUMMARY = {"comet_high_res": {"engine": "Comet 2026.01.1", "spectra": "high_res", "psms": 3, "decoy_psms": 1, "exact": 3,
                              "within_tolerance": 3, "tolerance": 0.001, "max_difference": 0.0, "not_scored": 0}}


def _fake_validate(returncode: int = 0):
    """validate.py's outputs without running it (the real run is tests/test_engine_reference_scores.py's business)."""

    def run(out: Path) -> subprocess.CompletedProcess:
        if returncode == 0:
            (out / "comparison.tsv").write_text("run\ttitle\n")
            (out / "summary.json").write_text(json.dumps(SUMMARY))
            (out / "agreement.png").write_bytes(b"png")
            (out / "agreement.pdf").write_bytes(b"pdf")
            (out / "high_res.mgf").write_text("BEGIN IONS\nEND IONS\n")
        return subprocess.CompletedProcess(["validate.py"], returncode, "", "FAIL comet_high_res" if returncode else "")

    return run


@pytest.fixture(autouse=True)
def _no_validation_run(script, monkeypatch):
    monkeypatch.setattr(script, "run_validate", _fake_validate())


def _deposit(tmp_path: Path) -> Path:
    archive = tmp_path / "supplement.zip"
    with zipfile.ZipFile(archive, "w") as bundle:
        bundle.writestr("supplement/raw_file_inventory.tsv", "run\n")
        bundle.writestr("supplement/figure_assets/pdfs/Figure2.pdf", "old figure 2")
        bundle.writestr("supplement/figure_assets/pdfs/Figure3.pdf", "figure 3")
        bundle.writestr("supplement/processed/arabidopsis/fragpipe_psm.tsv", "Spectrum\tPeptide\n")
        bundle.writestr("__MACOSX/supplement/._Figure2.pdf", "resource fork")
    return archive


def _figures(tmp_path: Path, status: str = "REGENERATED") -> Path:
    figures = tmp_path / "figures"
    staged = figures / "deposit"
    (staged / "figure_assets/pdfs").mkdir(parents=True)
    (staged / "processed/arabidopsis").mkdir(parents=True)
    (staged / "figure_assets/pdfs/Figure2.pdf").write_text("new figure 2")
    (staged / "figure_assets/pdfs/Figure3.pdf").write_text("figure 3")
    (staged / "processed/arabidopsis/figure4_usi_spectra.mgf").write_text("BEGIN IONS\nEND IONS\n")
    manifest = {"package": {"proteomics-mcp": "0.1.10"}, "git": {"describe": "v0.1.10-9-gabc", "sha": "abc", "dirty": False}, "ms2pip": "4.1.2",
                "python": "3.12.13", "command": f"uv run python scripts/regenerate_figures.py --out {Path.home()}/figures",
                "verification": {"problems": [], "passed": True}, "figure4": {"status": status, "checks": {"100 PSMs scored": True}},
                "examples": {}}
    (figures / "manifest.json").write_text(json.dumps(manifest))
    (figures / "CHANGES.md").write_text("# Figure changes\n")
    return figures


def test_package_keeps_the_deposit_and_records_every_change(script, tmp_path, monkeypatch):
    archive = _deposit(tmp_path)
    monkeypatch.setattr(script, "ZENODO_MD5", hashlib.md5(archive.read_bytes()).hexdigest())
    out = tmp_path / "package"
    assert script.main(["--figures", str(_figures(tmp_path)), "--zenodo-zip", str(archive), "--out", str(out)]) == 0

    supplement = out / "supplement"
    assert (supplement / "figure_assets/pdfs/Figure2.pdf").read_text() == "new figure 2"
    assert (supplement / "processed/arabidopsis/fragpipe_psm.tsv").read_text() == "Spectrum\tPeptide\n"  # input kept
    assert not (out / "__MACOSX").exists()
    changes = (supplement / "CHANGES.md").read_text()
    assert "`figure_assets/pdfs/Figure2.pdf`" in changes.split("## Regenerated (content changed)")[1].split("##")[0]
    assert "`processed/arabidopsis/figure4_usi_spectra.mgf`" in changes.split("## Added")[1].split("##")[0]
    assert "`figure_assets/pdfs/Figure3.pdf`" in changes.split("## Regenerated, identical")[1].split("##")[0]
    assert "[#51](https://github.com/UadKLab/proteomics-mcp/issues/51)" in changes
    listed = {line.split("  ", 1)[1]: line.split("  ", 1)[0] for line in (supplement / "MANIFEST.sha256").read_text().splitlines()}
    files = {p.relative_to(supplement).as_posix() for p in supplement.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256"}
    assert set(listed) == files
    assert listed["CHANGES.md"] == hashlib.sha256((supplement / "CHANGES.md").read_bytes()).hexdigest()


def test_package_refuses_a_provisional_run(script, tmp_path, monkeypatch):
    archive = _deposit(tmp_path)
    monkeypatch.setattr(script, "ZENODO_MD5", hashlib.md5(archive.read_bytes()).hexdigest())
    figures = _figures(tmp_path, status="PROVISIONAL")
    assert script.main(["--figures", str(figures), "--zenodo-zip", str(archive), "--out", str(tmp_path / "package")]) == 2


def test_provenance_records_versions_conversion_and_no_home_paths(script, tmp_path, monkeypatch):
    archive = _deposit(tmp_path)
    monkeypatch.setattr(script, "ZENODO_MD5", hashlib.md5(archive.read_bytes()).hexdigest())
    out = tmp_path / "package"
    assert script.main(["--figures", str(_figures(tmp_path)), "--zenodo-zip", str(archive), "--out", str(out)]) == 0

    supplement = out / "supplement"
    provenance = (supplement / "PROVENANCE.md").read_text()
    assert "**0.2.0 rc**" in provenance and "MS2PIP:** 4.1.2" in provenance and "3.12.13" in provenance
    assert "ThermoRawFileParser 2.0.0-dev" in provenance and "e00b388" in provenance and "PMID 31755270" in provenance
    assert "0.1.8 / 0.1.10" in provenance and "--out ~/figures" in provenance
    for path in [supplement / "PROVENANCE.md", *(supplement / "regeneration").rglob("*")]:
        assert str(Path.home()) not in path.read_text(), path


def test_recorded_paths_use_forward_slashes_and_every_change_has_a_reason(script, tmp_path, monkeypatch):
    """Paths are written with / on every platform (backslashes on Windows broke the reason lookup), and every file
    in CHANGES.md says why it changed."""
    archive = _deposit(tmp_path)
    monkeypatch.setattr(script, "ZENODO_MD5", hashlib.md5(archive.read_bytes()).hexdigest())
    out = tmp_path / "package"
    assert script.main(["--figures", str(_figures(tmp_path)), "--zenodo-zip", str(archive), "--out", str(out)]) == 0

    supplement = out / "supplement"
    changes, manifest = (supplement / "CHANGES.md").read_text(), (supplement / "MANIFEST.sha256").read_text()
    commands = (supplement / "PROVENANCE.md").read_text().split("```bash")[1].split("```")[0]
    for name, text in (("CHANGES.md", changes), ("MANIFEST.sha256", manifest), ("PROVENANCE.md commands", commands)):
        assert "\\" not in text, name
    rows = [line.split("|")[1:-1] for line in changes.splitlines() if line.startswith("| `")]
    assert rows and all(len(row) == 2 and row[1].strip() for row in rows), [row for row in rows if not row[-1].strip()]


def test_home_directory_is_hidden_as_written_on_windows(script, monkeypatch):
    monkeypatch.setattr(script.Path, "home", lambda: PureWindowsPath(r"C:\Users\me"))
    text = json.dumps({"cwd": r"C:\Users\me\pmcp", "posix": "C:/Users/me/pmcp"}) + "\n" + r"C:\Users\me\pmcp"
    assert "Users" not in script.sanitized(text)


@pytest.mark.parametrize(("platform", "noted"), [("macOS-15.6-arm64-arm-64bit", True), ("Linux-6.8.0-x86_64-with-glibc2.39", False)])
def test_provenance_notes_the_macos_openmp_runtime_and_the_packaging_commit(script, platform, noted):
    figures = {"package": {"proteomics-mcp": "0.2.0"}, "git": {"describe": "v0.2.0", "sha": "abc", "dirty": False}, "ms2pip": "4.1.2",
               "python": "3.12.13", "platform": platform, "command": "uv run python scripts/regenerate_figures.py"}
    provenance = script.provenance_markdown(figures, "0.2.0", "uv run python scripts/build_zenodo_package.py",
                                            {"describe": "v0.2.0-3-gdef", "sha": "def", "dirty": False})
    assert ("bundled with pyOpenMS" in provenance and "OpenMP Error #15" in provenance) is noted
    assert "**Packaging:** `scripts/build_zenodo_package.py` at `v0.2.0-3-gdef` (commit `def`)." in provenance
    assert "proteomics-mcp **0.2.0**, built from the tag `v0.2.0`" in provenance


def test_package_adds_the_engine_score_validation(script, tmp_path, monkeypatch):
    """#154: the validation's comparison, summary, plots, versions.json and the commands to rerun it go into validation/."""
    archive = _deposit(tmp_path)
    monkeypatch.setattr(script, "ZENODO_MD5", hashlib.md5(archive.read_bytes()).hexdigest())
    out = tmp_path / "package"
    assert script.main(["--figures", str(_figures(tmp_path)), "--zenodo-zip", str(archive), "--out", str(out)]) == 0

    validation = out / "supplement" / "validation"
    assert sorted(p.name for p in validation.iterdir()) == ["README.md", "agreement.pdf", "agreement.png", "comparison.tsv", "summary.json",
                                                            "versions.json"]
    versions = json.loads((validation / "versions.json").read_text())
    assert versions == json.loads((script.VALIDATION / "golden" / "versions.json").read_text())
    assert "@sha256:" in versions["comet_image"] and "comet_high_res.params" in versions["inputs_sha256"]
    readme = (validation / "README.md").read_text()
    assert "3 of 3 PSMs are equal" in readme and "git switch --detach v0.2.0" in readme
    assert f"docker pull {versions['comet_image']}" in readme and "validate.py --run-engines --out out" in readme
    changes = (out / "supplement" / "CHANGES.md").read_text()
    added = changes.split("## Added")[1].split("\n## ")[0]
    assert "`validation/README.md`" in added and "`validation/versions.json`" in added and "#154" in added
    assert "Engine-score validation outputs" not in changes  # no longer pending
    assert "**Engine-score validation:** `validation/`" in (out / "supplement" / "PROVENANCE.md").read_text()
    listed = {line.split("  ", 1)[1] for line in (out / "supplement" / "MANIFEST.sha256").read_text().splitlines()}
    assert "validation/comparison.tsv" in listed


def test_package_stops_when_the_validation_fails(script, tmp_path, monkeypatch):
    archive = _deposit(tmp_path)
    monkeypatch.setattr(script, "ZENODO_MD5", hashlib.md5(archive.read_bytes()).hexdigest())
    monkeypatch.setattr(script, "run_validate", _fake_validate(returncode=1))
    out = tmp_path / "package"
    assert script.main(["--figures", str(_figures(tmp_path)), "--zenodo-zip", str(archive), "--out", str(out)]) == 5
    assert not (out / "supplement" / "MANIFEST.sha256").exists()


def _tables_run(tmp_path: Path, figures: Path, *, source: str = "s" * 64, passed: bool = True) -> Path:
    """A regenerate_figures.py --supplementary-tables run whose tables are the figures run's."""
    run = tmp_path / "tables"
    (run / "deposit/figure_assets/pdfs").mkdir(parents=True)
    tables = {}
    for name, n, html in (("Suppl_figure3", 3, "psm"), ("Suppl_figure4", 4, "peptide")):
        (run / f"deposit/figure_assets/pdfs/{name}.pdf").write_text(f"new {name}")
        rel = f"figure_assets/interactive/Supplementary_Figure_{n}_fragpipe_scored_{html}_table.html"
        (figures / "deposit" / rel).parent.mkdir(parents=True, exist_ok=True)
        (figures / "deposit" / rel).write_text(f"table {n}")
        tables[name] = {"html": rel, "html_sha256": hashlib.sha256(f"table {n}".encode()).hexdigest(), "header": f"{n} rows, {n} columns"}
    manifest = {"mode": "supplementary-tables", "package": {"source_sha256": source}, "git": {"describe": "v0.2.0-9-gfff", "sha": "fff"},
                "ms2pip": "4.1.2", "source_of": {"tag": "v0.2.0", "sha256": source, "equal": True}, "tables": tables,
                "verification": {"problems": [], "passed": passed}}
    (run / "manifest.json").write_text(json.dumps(manifest))
    return run


def test_package_takes_supplementary_figures_3_and_4_from_a_tables_run(script, tmp_path, monkeypatch):
    archive = _deposit(tmp_path)
    with zipfile.ZipFile(archive, "a") as bundle:
        bundle.writestr("supplement/figure_assets/pdfs/Suppl_figure3.pdf", "old picture 3")
        bundle.writestr("supplement/figure_assets/pdfs/Suppl_figure4.pdf", "old picture 4")
    monkeypatch.setattr(script, "ZENODO_MD5", hashlib.md5(archive.read_bytes()).hexdigest())
    figures = _figures(tmp_path)
    manifest = json.loads((figures / "manifest.json").read_text())
    manifest["package"]["source_sha256"] = "s" * 64
    (figures / "manifest.json").write_text(json.dumps(manifest))
    out = tmp_path / "package"
    assert script.main(["--figures", str(figures), "--zenodo-zip", str(archive), "--out", str(out),
                        "--supplementary-tables", str(_tables_run(tmp_path, figures))]) == 0

    supplement = out / "supplement"
    assert (supplement / "figure_assets/pdfs/Suppl_figure4.pdf").read_text() == "new Suppl_figure4"
    changes = (supplement / "CHANGES.md").read_text()
    regenerated = changes.split("## Regenerated (content changed)")[1].split("\n## ")[0]
    assert "`figure_assets/pdfs/Suppl_figure3.pdf` | screenshot of the 0.2.0 interactive PSM table" in regenerated
    assert "`figure_assets/pdfs/Suppl_figure4.pdf`" in regenerated and "62 columns" in regenerated
    assert "Suppl_figure4.pdf" not in changes.split("## Not regenerated")[1].split("\n## ")[0]
    assert (supplement / "regeneration/supplementary_tables_manifest.json").is_file()
    provenance = (supplement / "PROVENANCE.md").read_text()
    assert "**Supplementary Figures 3 and 4:**" in provenance and "(3 rows, 3 columns; 4 rows, 4 columns)" in provenance
    assert "at `v0.2.0-9-gfff` (commit `fff`)" in provenance and f"identical to `v0.2.0` (sha256 `{'s' * 64}`)" in provenance


@pytest.mark.parametrize("problem", ["source", "failed", "table"])
def test_package_refuses_a_tables_run_that_does_not_fit(script, tmp_path, monkeypatch, problem):
    archive = _deposit(tmp_path)
    monkeypatch.setattr(script, "ZENODO_MD5", hashlib.md5(archive.read_bytes()).hexdigest())
    figures = _figures(tmp_path)
    manifest = json.loads((figures / "manifest.json").read_text())
    manifest["package"]["source_sha256"] = "s" * 64
    (figures / "manifest.json").write_text(json.dumps(manifest))
    run = _tables_run(tmp_path, figures, source="t" * 64 if problem == "source" else "s" * 64, passed=problem != "failed")
    if problem == "table":
        (figures / "deposit/figure_assets/interactive/Supplementary_Figure_4_fragpipe_scored_peptide_table.html").write_text("other table")
    assert script.main(["--figures", str(figures), "--zenodo-zip", str(archive), "--out", str(tmp_path / "package"),
                        "--supplementary-tables", str(run)]) == 6
