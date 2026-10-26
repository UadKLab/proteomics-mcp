from __future__ import annotations

import bz2
import csv
import gzip
import io
import lzma
import zipfile
from pathlib import Path
from typing import IO

import pandas as pd

_COMPRESSION_SUFFIXES = {".gz", ".bz2", ".xz", ".zip"}
_EXCEL_SUFFIXES = {".xlsx", ".xlsm", ".xls"}
_TAB_SUFFIXES = {".tsv", ".tab", ".txt"}
_DELIMITERS = (",", "\t", ";")
_DELIMITER_NAMES = {",": "comma", "\t": "tab", ";": "semicolon"}


def _table_suffix(path: Path) -> str:
    """The format suffix, looking through compression: "psm.tsv.gz" -> ".tsv"."""
    suffixes = [suffix.lower() for suffix in path.suffixes]
    while suffixes and suffixes[-1] in _COMPRESSION_SUFFIXES:
        suffixes.pop()
    return suffixes[-1] if suffixes else ""


def _open_text(path: Path) -> IO[str]:
    suffix = path.suffix.lower()
    if suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    if suffix == ".bz2":
        return bz2.open(path, "rt", encoding="utf-8", errors="replace")
    if suffix == ".xz":
        return lzma.open(path, "rt", encoding="utf-8", errors="replace")
    if suffix == ".zip":
        archive = zipfile.ZipFile(path)
        return io.TextIOWrapper(archive.open(archive.namelist()[0]), encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def _header_delimiter(path: Path) -> str | None:
    """The delimiter that splits the header line into the most fields, or None when none splits it."""
    with _open_text(path) as handle:
        header = handle.readline()
    counts = {delimiter: len(next(csv.reader([header], delimiter=delimiter), [])) for delimiter in _DELIMITERS}
    delimiter, fields = max(counts.items(), key=lambda item: item[1])
    return delimiter if fields > 1 else None


def read_candidate_table(path: str | Path, *, limit: int | None = None, sheet_name: str | int = 0) -> pd.DataFrame:
    """Read an identification table, choosing the format from the extension (#140).

    .csv is comma-separated; .tsv, .tab and .txt are tab-separated; Excel files are read as Excel; compressed files
    (.gz, .bz2, .xz, .zip) by the suffix before the compression one. Without a known extension the delimiter is taken
    from the header line. A file whose content does not match its extension (a TSV named .csv) is read again with the
    delimiter its header line uses, so that it loads instead of giving one column or a parser error.
    """
    table_path = Path(path)
    suffix = _table_suffix(table_path)
    if suffix in _EXCEL_SUFFIXES:
        return pd.read_excel(table_path, nrows=limit, sheet_name=sheet_name)
    if suffix == ".csv":
        sep = ","
    elif suffix in _TAB_SUFFIXES:
        sep = "\t"
    else:
        sep = _header_delimiter(table_path) or ","

    try:
        frame = pd.read_csv(table_path, sep=sep, nrows=limit)
    except pd.errors.ParserError as exc:
        detected = _header_delimiter(table_path)
        if detected is None or detected == sep:
            raise ValueError(f"Could not read {table_path.name} as a {_DELIMITER_NAMES[sep]}-separated table: {exc}") from exc
        return pd.read_csv(table_path, sep=detected, nrows=limit)

    if frame.shape[1] == 1:
        detected = _header_delimiter(table_path)
        if detected is not None and detected != sep:
            frame = pd.read_csv(table_path, sep=detected, nrows=limit)
    return frame


def describe_columns(frame: pd.DataFrame, limit: int = 12) -> str:
    """The first columns of a table, for error messages."""
    columns = [str(column) for column in frame.columns[:limit]]
    more = f", … ({len(frame.columns)} in total)" if len(frame.columns) > limit else ""
    return ", ".join(columns) + more
