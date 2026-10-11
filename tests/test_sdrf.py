from __future__ import annotations

import io
import json
from pathlib import Path
from urllib.error import HTTPError

from proteomics_mcp.core.sdrf import load_sdrf_summary, retrieve_pride_sdrf


def test_load_sdrf_summary_extracts_file_mapping(tmp_path: Path):
    path = tmp_path / "example.sdrf.tsv"
    path.write_text(
        "\t".join(
            [
                "source name",
                "assay name",
                "comment[data file]",
                "characteristics[organism]",
                "comment[label]",
                "factor value[condition]",
            ]
        )
        + "\n"
        + "\t".join(["sample-1", "assay-1", "run1.raw", "Homo sapiens", "TMT126", "treated"])
        + "\n",
        encoding="utf-8",
    )

    summary = load_sdrf_summary(path)

    assert summary["row_count"] == 1
    assert summary["data_files"] == ["run1.raw"]
    assert summary["labels"] == ["TMT126"]
    assert summary["file_metadata"]["run1"]["organism"] == "Homo sapiens"
    assert summary["file_metadata"]["run1.raw"]["factor_values"]["factor value[condition]"] == "treated"


def test_retrieve_pride_sdrf_uses_project_files_listing(monkeypatch):
    project_listing = {
        "_embedded": {
            "files": [
                {
                    "fileName": "PXD000001.sdrf.tsv",
                    "publicFileLocations": ["https://example.org/PXD000001.sdrf.tsv"],
                }
            ]
        }
    }
    sdrf_text = (
        "source name\tassay name\tcomment[data file]\tcharacteristics[organism]\n"
        "sample-1\tassay-1\trun1.raw\tHomo sapiens\n"
    )

    def fake_urlopen(url: str, timeout: float = 20.0):
        if "projects/PXD000001/files" in url:
            return io.BytesIO(json.dumps(project_listing).encode("utf-8"))
        if "PXD000001.sdrf.tsv" in url:
            return io.BytesIO(sdrf_text.encode("utf-8"))
        raise AssertionError(url)

    monkeypatch.setattr("proteomics_mcp.core.sdrf.urlopen", fake_urlopen)

    summary = retrieve_pride_sdrf("PXD000001")

    assert summary["project_accession"] == "PXD000001"
    assert summary["file_name"] == "PXD000001.sdrf.tsv"
    assert summary["data_files"] == ["run1.raw"]
    assert summary["file_metadata"]["run1.raw"]["organism"] == "Homo sapiens"


def test_retrieve_pride_sdrf_reports_pride_http_failure(monkeypatch):
    def fake_urlopen(url: str, timeout: float = 20.0):
        raise HTTPError(url, 404, "Not Found", hdrs=None, fp=io.BytesIO(b'{"message":"unknown project"}'))

    monkeypatch.setattr("proteomics_mcp.core.sdrf.urlopen", fake_urlopen)

    try:
        retrieve_pride_sdrf("pxd999999")
    except LookupError as exc:
        assert "HTTP 404" in str(exc)
        assert "PXD999999" in str(exc)
    else:
        raise AssertionError("Expected a diagnostic LookupError")
