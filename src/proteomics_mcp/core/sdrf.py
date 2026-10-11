from __future__ import annotations

import io
import json
import gzip
from pathlib import Path
from typing import Any
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode
from urllib.request import urlopen

import pandas as pd


PRIDE_API_BASE = "https://www.ebi.ac.uk/pride/ws/archive/v2"


def _read_json_url(url: str, timeout_seconds: float) -> Any:
    with urlopen(url, timeout=timeout_seconds) as handle:  # noqa: S310
        return json.loads(handle.read().decode("utf-8"))


def _read_text_url(url: str, timeout_seconds: float) -> str:
    with urlopen(url, timeout=timeout_seconds) as handle:  # noqa: S310
        content = handle.read()
    if content.startswith(b"\x1f\x8b"):
        content = gzip.decompress(content)
    return content.decode("utf-8")


def load_sdrf_table(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str).fillna("")


def _normalized_column_name(column: str) -> str:
    return str(column).strip().lower()


def _matching_columns(frame: pd.DataFrame, prefix: str) -> list[str]:
    return [col for col in frame.columns if _normalized_column_name(col).startswith(prefix)]


def _first_present(row: dict[str, Any], columns: list[str]) -> str | None:
    for column in columns:
        value = str(row.get(column, "")).strip()
        if value:
            return value
    return None


def _file_stem(path_value: str | None) -> str | None:
    if not path_value:
        return None
    return Path(str(path_value)).stem


def summarize_sdrf(frame: pd.DataFrame, *, source: str | None = None) -> dict[str, Any]:
    data_file_columns = _matching_columns(frame, "comment[data file]") + _matching_columns(frame, "comment[data_file]")
    fraction_file_columns = _matching_columns(frame, "comment[fraction file]")
    label_columns = _matching_columns(frame, "comment[label]")
    instrument_columns = _matching_columns(frame, "comment[instrument]")
    organism_columns = _matching_columns(frame, "characteristics[organism]")
    organism_part_columns = _matching_columns(frame, "characteristics[organism part]")
    disease_columns = _matching_columns(frame, "characteristics[disease]")
    factor_columns = _matching_columns(frame, "factor value[")
    assay_columns = _matching_columns(frame, "assay name")
    source_name_columns = _matching_columns(frame, "source name")
    comment_columns = [col for col in frame.columns if _normalized_column_name(col).startswith("comment[")]
    characteristics_columns = [col for col in frame.columns if _normalized_column_name(col).startswith("characteristics[")]

    file_metadata: dict[str, dict[str, Any]] = {}
    assays: list[str] = []
    for row in frame.to_dict(orient="records"):
        assay_name = _first_present(row, assay_columns)
        if assay_name:
            assays.append(assay_name)
        data_file = _first_present(row, data_file_columns)
        if not data_file:
            continue
        entry = {
            "data_file": data_file,
            "fraction_file": _first_present(row, fraction_file_columns),
            "assay_name": assay_name,
            "source_name": _first_present(row, source_name_columns),
            "label": _first_present(row, label_columns),
            "instrument": _first_present(row, instrument_columns),
            "organism": _first_present(row, organism_columns),
            "organism_part": _first_present(row, organism_part_columns),
            "disease": _first_present(row, disease_columns),
            "factor_values": {col: row.get(col) for col in factor_columns if str(row.get(col, "")).strip()},
            "characteristics": {col: row.get(col) for col in characteristics_columns if str(row.get(col, "")).strip()},
            "comments": {col: row.get(col) for col in comment_columns if str(row.get(col, "")).strip()},
        }
        for key in {data_file, Path(data_file).name, _file_stem(data_file)}:
            if key:
                file_metadata[str(key)] = entry

    unique_values = lambda cols: sorted({str(value).strip() for col in cols for value in frame[col].tolist() if str(value).strip()})
    return {
        "source": source,
        "row_count": int(len(frame)),
        "column_count": int(len(frame.columns)),
        "columns": list(map(str, frame.columns)),
        "data_file_columns": data_file_columns,
        "fraction_file_columns": fraction_file_columns,
        "factor_columns": factor_columns,
        "comment_columns": comment_columns,
        "characteristics_columns": characteristics_columns,
        "assays": sorted(set(assays)),
        "labels": unique_values(label_columns),
        "instruments": unique_values(instrument_columns),
        "organisms": unique_values(organism_columns),
        "organism_parts": unique_values(organism_part_columns),
        "diseases": unique_values(disease_columns),
        "data_files": sorted({entry["data_file"] for entry in file_metadata.values() if entry.get("data_file")}),
        "file_metadata": file_metadata,
        "rows_preview": frame.head(10).to_dict(orient="records"),
    }


def load_sdrf_summary(path: str | Path) -> dict[str, Any]:
    frame = load_sdrf_table(path)
    return summarize_sdrf(frame, source=str(path))


def _extract_pride_page_content(payload: Any) -> list[dict[str, Any]]:
    if isinstance(payload, dict):
        embedded = payload.get("_embedded")
        if isinstance(embedded, dict):
            for key in ("files", "fileList", "projectFiles"):
                values = embedded.get(key)
                if isinstance(values, list):
                    return [value for value in values if isinstance(value, dict)]
        for key in ("list", "content", "files"):
            values = payload.get(key)
            if isinstance(values, list):
                return [value for value in values if isinstance(value, dict)]
    if isinstance(payload, list):
        return [value for value in payload if isinstance(value, dict)]
    return []


def _find_download_url(file_record: dict[str, Any]) -> str | None:
    for key in ("downloadLink", "download_url", "publicFileLocation", "ftpLink", "fileDownloadUrl"):
        value = file_record.get(key)
        if value:
            return str(value)
    locations = file_record.get("publicFileLocations") or file_record.get("publicFileLocation")
    if isinstance(locations, list):
        for value in locations:
            if value:
                return str(value)
    if isinstance(locations, str) and locations:
        return locations
    links = file_record.get("_links")
    if isinstance(links, dict):
        for key in ("download", "self"):
            entry = links.get(key)
            if isinstance(entry, dict) and entry.get("href"):
                return str(entry["href"])
    return None


def _is_sdrf_record(file_record: dict[str, Any]) -> bool:
    name = str(file_record.get("fileName", "")).strip().lower()
    if "sdrf" not in name:
        return False
    return name.endswith((".sdrf", ".sdrf.tsv", ".sdrf.txt", ".tsv", ".txt", ".sdrf.gz", ".tsv.gz", ".txt.gz"))


def _read_pride_file_listing(url: str, timeout_seconds: float, project_accession: str) -> Any:
    try:
        return _read_json_url(url, timeout_seconds)
    except HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace").strip()
        message = f"PRIDE returned HTTP {exc.code} while listing files for {project_accession} at {url}."
        if detail:
            message += f" Response: {detail[:500]}"
        raise LookupError(message) from exc
    except URLError as exc:
        raise ConnectionError(f"Could not reach PRIDE while listing files for {project_accession} at {url}: {exc.reason}") from exc


def retrieve_pride_sdrf(project_accession: str, *, timeout_seconds: float = 20.0) -> dict[str, Any]:
    project_accession = str(project_accession).strip().upper()
    if not project_accession:
        raise ValueError("project_accession must be a non-empty PRIDE accession such as PXD012345.")
    params = urlencode({"pageSize": 1000, "page": 0})
    files_url = f"{PRIDE_API_BASE}/projects/{project_accession}/files?{params}"
    payload = _read_pride_file_listing(files_url, timeout_seconds, project_accession)
    files = _extract_pride_page_content(payload)
    sdrf_record = next((record for record in files if _is_sdrf_record(record)), None)
    if sdrf_record is None:
        names = [str(record.get("fileName", "")) for record in files if record.get("fileName")]
        preview = ", ".join(names[:10]) or "no files returned"
        raise LookupError(f"No SDRF file found in PRIDE project {project_accession}. Files returned: {preview}.")
    download_url = _find_download_url(sdrf_record)
    if not download_url:
        raise LookupError(f"PRIDE project {project_accession} has an SDRF record but no downloadable URL.")
    text = _read_text_url(download_url, timeout_seconds)
    frame = pd.read_csv(io.StringIO(text), sep="\t", dtype=str).fillna("")
    summary = summarize_sdrf(frame, source=download_url)
    summary["project_accession"] = project_accession
    summary["file_name"] = sdrf_record.get("fileName")
    summary["download_url"] = download_url
    return summary
