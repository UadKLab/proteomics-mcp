"""scripts/build_zenodo_package.py on a small fake deposit, without network."""

import hashlib
import importlib.util
import json
import zipfile
from pathlib import Path

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
    files = {str(p.relative_to(supplement)) for p in supplement.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256"}
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
