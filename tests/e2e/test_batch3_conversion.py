"""Batch 3 of #30 (conversions): mzML, MGF and Bruker .d to Feather, and the InstaNovo export, on public files.

Inputs (downloaded once into PROTEOMICS_MCP_E2E_CACHE and sha256-checked):
- mzML: OpenMS example BSA1.mzML (tag Release3.1.0), a real BSA DDA run with 1,120 MS2 spectra
- MGF:  PRIDE PXD000001, PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz (6,103 spectra)
- .d:   PRIDE PXD079916, 20260311_SL3_DDA_WT_i01_BA5_1_2553.d.zip (timsTOF HT ddaPASEF, 522 MB);
        only with PROTEOMICS_MCP_E2E_BRUKER=1
"""

from __future__ import annotations

import gzip
import os
import re
import shutil
import zipfile
from pathlib import Path

import pandas as pd
import pytest


pytestmark = pytest.mark.e2e

BSA_MZML = ("https://raw.githubusercontent.com/OpenMS/OpenMS/Release3.1.0/share/OpenMS/examples/BSA/BSA1.mzML",
            "openms-Release3.1.0-BSA1.mzML", "dc9ed61d595328d4ef2f1de47d21f41b83e2eae7c9145e1d9b88e910c8cec2f7")
PXD000001_MGF = ("https://ftp.pride.ebi.ac.uk/pride/data/archive/2012/03/PXD000001/generated/PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz",
                 "PXD000001-PRIDE_Exp_Complete_Ac_22134.pride.mgf.gz",
                 "9b35ea49b37c2a8c96d3a00433779675bd16935fcfccdd13b04cc16a8bce3fea")
PXD079916_D = ("https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/09/PXD079916/20260311_SL3_DDA_WT_i01_BA5_1_2553.d.zip",
               "PXD079916-20260311_SL3_DDA_WT_i01_BA5_1_2553.d.zip",
               "9c357b83fc434ba57d2e062a851f2b60a6a861710d853ce518f29d09ca551913")


def _spectra_summary(path: Path) -> dict:
    frame = pd.read_feather(path)
    peaks = frame["mz_array"].map(len)
    sorted_mz = frame["mz_array"].map(lambda mz: all(a <= b for a, b in zip(mz, mz[1:])))
    return {"spectra": len(frame), "with_precursor_mz": int(frame["precursor_mz"].notna().sum()),
            "with_charge": int(frame["precursor_charge"].notna().sum()), "empty": int((peaks == 0).sum()),
            "median_peaks": float(peaks.median()), "mz_sorted": bool(sorted_mz.all()),
            "same_lengths": bool((peaks == frame["intensity_array"].map(len)).all()), "columns": sorted(frame.columns)}


def _convert(tool, cli, tool_name: str, command: str, arguments: dict, cli_args: list[str], tmp_path: Path) -> tuple[dict, Path, Path]:
    mcp_out, cli_out = tmp_path / "mcp.feather", tmp_path / "cli.feather"
    is_error, result = tool(tool_name, {**arguments, "output_path": str(mcp_out)})
    assert not is_error, result
    completed = cli(command, *cli_args, "--output-path", str(cli_out), timeout=3600)
    assert completed.returncode == 0, completed.stderr[-2000:]
    return result, mcp_out, cli_out


def _same_spectra(a: Path, b: Path) -> bool:
    left, right = pd.read_feather(a), pd.read_feather(b)
    columns = ["scan_number", "precursor_mz", "precursor_charge"]
    if len(left) != len(right) or not left[columns].equals(right[columns]):
        return False
    return all(list(x) == list(y) for x, y in zip(left["mz_array"], right["mz_array"])) and \
        all(list(x) == list(y) for x, y in zip(left["intensity_array"], right["intensity_array"]))


@pytest.fixture(scope="module")
def bsa_mzml(fetch) -> Path:
    return fetch(*BSA_MZML)


@pytest.fixture(scope="module")
def pxd000001_mgf(tmp_path_factory, fetch) -> Path:
    archive = fetch(*PXD000001_MGF)
    target = tmp_path_factory.mktemp("mgf") / "PRIDE_Exp_Complete_Ac_22134.mgf"
    with gzip.open(archive, "rb") as source, target.open("wb") as handle:
        shutil.copyfileobj(source, handle)
    return target


def test_mzml_to_feather(tool, cli, report, tmp_path, bsa_mzml):
    expected = len(re.findall(r'name="ms level" value="2"', bsa_mzml.read_text(encoding="utf-8")))
    result, mcp_out, cli_out = _convert(tool, cli, "mzml_to_feather", "mzml-to-feather", {"mzml_path": str(bsa_mzml)},
                                        ["--mzml-path", str(bsa_mzml)], tmp_path)
    summary = _spectra_summary(mcp_out)
    same = _same_spectra(mcp_out, cli_out)
    report("conversion", "mzml_to_feather (OpenMS BSA1.mzML)", "pass" if summary["spectra"] == expected and same else "fail",
           {"expected_ms2": expected, **summary, "cli_equals_mcp": same})
    assert summary["spectra"] == expected
    assert summary["with_precursor_mz"] == expected and summary["same_lengths"] and summary["mz_sorted"]
    assert same


def test_mgf_to_feather(tool, cli, report, tmp_path, pxd000001_mgf):
    text = pxd000001_mgf.read_text(encoding="utf-8")
    expected = text.count("BEGIN IONS")
    first_peaks = re.search(r"BEGIN IONS\n(?:[A-Z_]+=.*\n)*([\d.eE+-]+)\s+([\d.eE+-]+)", text)
    result, mcp_out, cli_out = _convert(tool, cli, "mgf_to_feather", "mgf-to-feather", {"mgf_path": str(pxd000001_mgf)},
                                        ["--mgf-path", str(pxd000001_mgf)], tmp_path)
    summary = _spectra_summary(mcp_out)
    frame = pd.read_feather(mcp_out)
    first = (float(frame["mz_array"].iloc[0][0]), float(frame["intensity_array"].iloc[0][0]))
    same = _same_spectra(mcp_out, cli_out)
    report("conversion", "mgf_to_feather (PXD000001)", "pass" if summary["spectra"] == expected and same else "fail",
           {"expected_spectra": expected, **summary, "first_peak": first, "cli_equals_mcp": same})
    assert summary["spectra"] == expected
    assert first == (float(first_peaks.group(1)), float(first_peaks.group(2)))
    assert summary["same_lengths"] and same


@pytest.mark.parametrize("output_format", ["ipc", "parquet"])
def test_export_to_instanovo(tool, cli, report, tmp_path, bsa_mzml, output_format):
    spectra = tmp_path / "bsa.feather"
    is_error, converted = tool("mzml_to_feather", {"mzml_path": str(bsa_mzml), "output_path": str(spectra)})
    assert not is_error, converted
    source = pd.read_feather(spectra)
    out = tmp_path / f"bsa_instanovo.{output_format}"
    is_error, result = tool("export_to_instanovo", {"input_path": str(spectra), "output_path": str(out), "output_format": output_format})
    assert not is_error, result
    exported = pd.read_parquet(out) if output_format == "parquet" else pd.read_feather(out)
    expected = int(source[["precursor_mz", "precursor_charge"]].notna().all(axis=1).sum())
    required = {"precursor_mz", "precursor_charge", "mz_array", "intensity_array"}
    cli_out = tmp_path / f"cli.{output_format}"
    completed = cli("export-to-instanovo", "--input-path", str(spectra), "--output-path", str(cli_out), "--format", output_format)
    report("conversion", f"export_to_instanovo ({output_format})", "pass" if len(exported) == expected else "fail",
           {"rows": len(exported), "source_spectra": len(source), "with_precursor_and_charge": expected,
            "columns": sorted(exported.columns), "cli_returncode": completed.returncode})
    assert required <= set(exported.columns)
    assert len(exported) == expected
    assert completed.returncode == 0, completed.stderr[-2000:]


def test_bruker_d_embedded_mgf_fallback(tool, report, tmp_path, pxd000001_mgf):
    """A .d directory without a readable analysis.tdf falls back to its embedded ProteinAnalysisResults.mgf (§2.3)."""
    d_dir = tmp_path / "no_tdf.d"
    d_dir.mkdir()
    shutil.copyfile(pxd000001_mgf, d_dir / "ProteinAnalysisResults.mgf")
    out = tmp_path / "fallback.feather"
    is_error, result = tool("bruker_d_to_feather", {"d_path": str(d_dir), "output_path": str(out)})
    assert not is_error, result
    summary = _spectra_summary(out)
    expected = pxd000001_mgf.read_text(encoding="utf-8").count("BEGIN IONS")
    report("conversion", "bruker_d_to_feather embedded-MGF fallback", "pass" if summary["spectra"] == expected else "fail",
           {"source_format": result.get("source_format"), "fallback_used": result.get("fallback_used"), **summary})
    assert result.get("source_format") == "bruker_d_embedded_mgf" and result.get("fallback_used") == "ProteinAnalysisResults.mgf"
    assert summary["spectra"] == expected


@pytest.mark.skipif(os.environ.get("PROTEOMICS_MCP_E2E_BRUKER") != "1", reason="set PROTEOMICS_MCP_E2E_BRUKER=1 for the 522 MB timsTOF run")
def test_bruker_d_alphatims(tool, report, tmp_path_factory, fetch):
    """The direct AlphaTims path on a real ddaPASEF run (PXD079916)."""
    archive = fetch(*PXD079916_D)
    root = tmp_path_factory.mktemp("bruker")
    with zipfile.ZipFile(archive) as bundle:
        bundle.extractall(root)
    d_dir = next(root.glob("*.d"))
    out = root / "timstof.feather"
    is_error, result = tool("bruker_d_to_feather", {"d_path": str(d_dir), "output_path": str(out)})
    report("conversion", "bruker_d_to_feather AlphaTims (PXD079916)", "fail" if is_error else "info",
           {"is_error": is_error, "result": str(result)[:500] if is_error else {k: v for k, v in result.items() if k != "columns"}})
    assert not is_error, result
    summary = _spectra_summary(out)
    report("conversion", "bruker_d_to_feather AlphaTims spectra", "pass" if summary["spectra"] > 0 else "fail", summary)
    assert result.get("source_format") != "bruker_d_embedded_mgf"
    assert summary["spectra"] > 1000 and summary["empty"] == 0
    assert summary["with_precursor_mz"] == summary["spectra"] and summary["same_lengths"]
    if "scan_index" in summary["columns"]:  # packages with #95: precursor Id as scan number, 1/K0 as mobility
        table = pd.read_feather(out, columns=["scan_number", "precursor_id", "mobility", "spectrum_id", "precursor_charge"])
        uncharged = table[table["precursor_charge"].isna()]
        checks = {
            "scan_number_unique": bool(table["scan_number"].is_unique),
            "scan_number_is_precursor_id": bool((table["scan_number"] == table["precursor_id"]).all()),
            "mobility_range": [float(table["mobility"].min()), float(table["mobility"].max())],
            "uncharged": int(len(uncharged)),
        }
        if len(uncharged):
            from proteomics_mcp.core.spectrum_loader import load_spectrum_from_table

            checks["uncharged_loads"] = load_spectrum_from_table(out, spectrum_id=uncharged["spectrum_id"].iloc[0]).precursor_charge is None
        report("conversion", "bruker_d_to_feather timsTOF identifiers (#95)",
               "pass" if checks["scan_number_unique"] and 0.5 <= checks["mobility_range"][0] <= checks["mobility_range"][1] <= 2.0 else "fail", checks)
        assert checks["scan_number_unique"] and checks["scan_number_is_precursor_id"]
        assert 0.5 <= checks["mobility_range"][0] <= checks["mobility_range"][1] <= 2.0
        assert checks.get("uncharged_loads", True)
