"""One rule for naming runs and for matching a PSM to its spectrum's run (#142, #126).

The metadata join (build_prediction_batch, build_id_dashboard), the spectrum lookup (score_candidates,
score_from_prediction_batch) and the converters all use these functions, so a PSM finds its spectrum the same way
everywhere.
"""

from __future__ import annotations

import math
from pathlib import PurePosixPath
from typing import Any

import pandas as pd

COMPRESSION_SUFFIXES = (".gz", ".bz2", ".xz", ".zip")
# Only these are removed from a file name to get the run name: a name such as QE-HF-01.5000 keeps its dot
SPECTRUM_FILE_SUFFIXES = (".mzml", ".mzxml", ".mgf", ".ms2", ".raw", ".d", ".wiff2", ".wiff", ".feather", ".arrow", ".ipc", ".parquet")
# The spectrum table's run column: the first of these that the table has and fills
SPECTRA_RUN_COLUMNS = ("experiment_name", "source_name", "source_file")


def run_name(value: Any) -> str | None:
    """The run name for a run or file name: the base name without compression and spectrum-file suffixes.

    "/data/A.mzML.gz", "A.mzML", "A.raw", "A.d" and "A" all give "A" (#126).
    """
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    text = str(value).strip().rstrip("/\\")
    if not text or text.lower() == "nan":
        return None
    name = PurePosixPath(text.replace("\\", "/")).name
    for suffixes in (COMPRESSION_SUFFIXES, SPECTRUM_FILE_SUFFIXES):
        lowered = name.lower()
        suffix = next((s for s in suffixes if lowered.endswith(s) and len(name) > len(s)), None)
        if suffix:
            name = name[: -len(suffix)]
    return name or None


def run_from_spectrum_id(spectrum_id: Any) -> tuple[str | None, int | None]:
    """The run and scan a spectrum identifier states: "<run>:<scan>", or a USI (mzspec:<dataset>:<run>:scan:<scan>)."""
    if spectrum_id is None or (isinstance(spectrum_id, float) and math.isnan(spectrum_id)):
        return None, None
    text = str(spectrum_id).strip()
    if text.startswith("mzspec:"):
        parts = text.split(":")
        if len(parts) >= 5 and parts[3] == "scan":
            return run_name(parts[2]), _int_or_none(parts[4])
        return (run_name(parts[2]) if len(parts) > 2 else None), None
    if ":" not in text or "=" in text:
        # native ids (controllerType=0 ... scan=12) carry no run
        return None, None
    run, rest = text.split(":", 1)
    return run_name(run), _int_or_none(rest)


def candidate_run(experiment_name: Any = None, spectrum_id: Any = None, run_key: Any = None) -> str | None:
    """A PSM's run: experiment_name when set, else run_key, else the run part of its spectrum id."""
    return run_name(experiment_name) or run_name(run_key) or run_from_spectrum_id(spectrum_id)[0]


def spectra_run_column(columns: Any, frame: pd.DataFrame | None = None) -> str | None:
    """The spectrum table's run column (SPECTRA_RUN_COLUMNS order); with a frame, the first that has any value."""
    for column in SPECTRA_RUN_COLUMNS:
        if column in columns and (frame is None or frame[column].notna().any()):
            return column
    return None


def spectra_run_names(frame: pd.DataFrame) -> tuple[pd.Series | None, str | None]:
    """The run name of every spectrum in a spectrum table, and the column it came from."""
    column = spectra_run_column(frame.columns, frame)
    if column is None:
        return None, None
    return frame[column].map(run_name), column


def _int_or_none(value: Any) -> int | None:
    try:
        return int(str(value).strip())
    except (TypeError, ValueError):
        return None
