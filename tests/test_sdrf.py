from __future__ import annotations

import io
import json
from pathlib import Path
from urllib.error import HTTPError
from urllib.parse import parse_qs, urlparse

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


def _pride_record(name: str) -> dict:
    # as the PRIDE v2 API lists files: an Aspera location first, then FTP
    path = f"pride/data/archive/2020/06/PXD000002/{name}"
    return {
        "fileName": name,
        "publicFileLocations": [
            {"name": "Aspera Protocol", "value": f"prd_ascp@fasp.ebi.ac.uk:{path}"},
            {"name": "FTP Protocol", "value": f"ftp://ftp.pride.ebi.ac.uk/{path}"},
        ],
    }


def test_retrieve_pride_sdrf_pages_through_the_file_listing(monkeypatch):
    """PRIDE returns at most 100 files per page; the SDRFs of PXD014877/PXD010595 are beyond the first page (#40)."""
    raw_files = [_pride_record(f"run{index:03d}.raw") for index in range(250)]
    listing = raw_files[:120] + [_pride_record("sdrf.tsv")] + raw_files[120:] + [_pride_record("PXD000002_community_annotated.sdrf.tsv")]
    sdrf_text = "source name\tcomment[data file]\nsample-1\trun001.raw\n"
    requested_pages = []

    def fake_urlopen(url: str, timeout: float = 20.0):
        if "projects/PXD000002/files" in url:
            query = parse_qs(urlparse(url).query)
            size, page = int(query["pageSize"][0]), int(query["page"][0])
            requested_pages.append((size, page))
            return io.BytesIO(json.dumps(listing[page * size:(page + 1) * size]).encode("utf-8"))
        if url == "https://ftp.pride.ebi.ac.uk/pride/data/archive/2020/06/PXD000002/PXD000002_community_annotated.sdrf.tsv":
            return io.BytesIO(sdrf_text.encode("utf-8"))
        raise AssertionError(url)

    monkeypatch.setattr("proteomics_mcp.core.sdrf.urlopen", fake_urlopen)

    summary = retrieve_pride_sdrf("PXD000002")

    assert requested_pages == [(100, 0), (100, 1), (100, 2)]
    assert summary["file_name"] == "PXD000002_community_annotated.sdrf.tsv"
    assert summary["sdrf_files_found"] == ["PXD000002_community_annotated.sdrf.tsv", "sdrf.tsv"]
    assert summary["download_url"].startswith("https://ftp.pride.ebi.ac.uk/")
    assert summary["data_files"] == ["run001.raw"]


def test_retrieve_pride_sdrf_reports_the_number_of_files_searched(monkeypatch):
    listing = [_pride_record(f"run{index:03d}.raw") for index in range(150)]

    def fake_urlopen(url: str, timeout: float = 20.0):
        query = parse_qs(urlparse(url).query)
        size, page = int(query["pageSize"][0]), int(query["page"][0])
        return io.BytesIO(json.dumps(listing[page * size:(page + 1) * size]).encode("utf-8"))

    monkeypatch.setattr("proteomics_mcp.core.sdrf.urlopen", fake_urlopen)
    try:
        retrieve_pride_sdrf("PXD000003")
    except LookupError as exc:
        assert "among the 150 files" in str(exc)
    else:
        raise AssertionError("Expected a LookupError")


def _write_sdrf(path: Path, header: list[str], rows: list[list[str]]) -> Path:
    path.write_text("\n".join("\t".join(row) for row in [header, *rows]) + "\n", encoding="utf-8")
    return path


def test_multiplexed_sdrf_keeps_every_channel(tmp_path: Path):
    # One row per TMT channel of the same raw file; the last row used to replace the others (#81).
    path = _write_sdrf(
        tmp_path / "tmt.sdrf.tsv",
        ["source name", "assay name", "comment[data file]", "comment[label]", "characteristics[organism]",
         "characteristics[disease]", "factor value[disease]"],
        [
            ["s1", "run1_126", "run1.raw", "TMT126", "Homo sapiens", "normal", "normal"],
            ["s2", "run1_127N", "run1.raw", "TMT127N", "Homo sapiens", "cancer", "cancer"],
            ["s3", "run2_126", "run2.raw", "TMT126", "Homo sapiens", "normal", "normal"],
        ],
    )

    run1 = load_sdrf_summary(path)["file_metadata"]["run1.raw"]

    assert run1["channel_count"] == 2
    assert [channel["label"] for channel in run1["channels"]] == ["TMT126", "TMT127N"]
    assert run1["label"] == "TMT126; TMT127N"
    assert run1["source_name"] == "s1; s2"
    assert run1["disease"] == "normal; cancer"
    assert run1["factor_values"] == {"factor value[disease]": "normal; cancer"}
    # shared values stay single
    assert run1["organism"] == "Homo sapiens"
    assert run1["data_file"] == "run1.raw"
    assert load_sdrf_summary(path)["file_metadata"]["run2"]["label"] == "TMT126"


def test_label_free_sdrf_is_unchanged(tmp_path: Path):
    path = _write_sdrf(
        tmp_path / "lfq.sdrf.tsv",
        ["source name", "assay name", "comment[data file]", "comment[label]", "factor value[condition]"],
        [["s1", "a1", "run1.raw", "label free sample", "treated"], ["s2", "a2", "run2.raw", "label free sample", "control"]],
    )

    summary = load_sdrf_summary(path)

    assert summary["file_metadata"]["run1"]["channel_count"] == 1
    assert summary["file_metadata"]["run1"]["factor_values"] == {"factor value[condition]": "treated"}
    assert summary["file_metadata"]["run2.raw"]["source_name"] == "s2"
    assert summary["data_files"] == ["run1.raw", "run2.raw"]
