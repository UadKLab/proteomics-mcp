from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd

from proteomics_mcp.core.spectrum_loader import load_spectrum_table


def export_table_to_instanovo(
    input_path: str | Path,
    output_path: str | Path,
    *,
    output_format: str = "ipc",
    sequence_column: str | None = None,
) -> dict[str, Any]:
    frame = load_spectrum_table(input_path)
    export = _build_instanovo_frame(frame, sequence_column=sequence_column)
    output_path = Path(output_path)
    output_format = str(output_format).strip().lower()
    if output_format == "ipc":
        export.to_feather(output_path)
    elif output_format == "parquet":
        export.to_parquet(output_path, index=False)
    else:
        raise ValueError("output_format must be one of: ipc, parquet")
    return {
        "input_path": str(input_path),
        "output_path": str(output_path),
        "output_format": output_format,
        "spectra_written": int(len(export)),
        "columns": export.columns.tolist(),
        "sequence_included": "sequence" in export.columns,
    }


def _build_instanovo_frame(frame: pd.DataFrame, *, sequence_column: str | None = None) -> pd.DataFrame:
    required = ["precursor_mz", "precursor_charge", "mz_array", "intensity_array"]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"Input spectrum table is missing required columns for InstaNovo export: {missing}")

    export = pd.DataFrame(index=frame.index)
    if "experiment_name" in frame.columns:
        export["experiment_name"] = frame["experiment_name"].astype(str)
    elif "source_name" in frame.columns:
        export["experiment_name"] = frame["source_name"].astype(str)
    elif "source_file" in frame.columns:
        export["experiment_name"] = frame["source_file"].map(lambda value: Path(str(value)).stem if value is not None else "spectra")
    else:
        export["experiment_name"] = "spectra"

    if "scan_number" in frame.columns:
        export["scan_number"] = pd.to_numeric(frame["scan_number"], errors="coerce").astype("Int64")
    else:
        export["scan_number"] = pd.Series(range(len(frame)), index=frame.index, dtype="Int64")

    if "spectrum_id" in frame.columns:
        export["spectrum_id"] = frame["spectrum_id"].astype(str)
    else:
        export["spectrum_id"] = [
            f"{experiment}:{int(scan) if pd.notna(scan) else idx}"
            for idx, (experiment, scan) in enumerate(zip(export["experiment_name"], export["scan_number"]))
        ]

    if "rt" in frame.columns:
        export["rt"] = pd.to_numeric(frame["rt"], errors="coerce")

    export["precursor_mz"] = pd.to_numeric(frame["precursor_mz"], errors="coerce")
    export["precursor_charge"] = pd.to_numeric(frame["precursor_charge"], errors="coerce").astype("Int64")
    export["mz_array"] = frame["mz_array"]
    export["intensity_array"] = frame["intensity_array"]

    if sequence_column:
        if sequence_column not in frame.columns:
            raise ValueError(f"sequence_column not found in input table: {sequence_column}")
        export["sequence"] = frame[sequence_column].astype(str)

    export = export.dropna(subset=["precursor_mz", "precursor_charge"])
    export["precursor_charge"] = export["precursor_charge"].astype(int)
    export["mz_array"] = export["mz_array"].map(_coerce_float_list)
    export["intensity_array"] = export["intensity_array"].map(_coerce_float_list)
    if "sequence" in export.columns:
        export["sequence"] = export["sequence"].map(lambda value: value.strip()).replace("", pd.NA)
    return export.reset_index(drop=True)


def _coerce_float_list(value: Any) -> list[float]:
    if value is None:
        return []
    if isinstance(value, list):
        return [float(item) for item in value]
    if hasattr(value, "tolist"):
        return [float(item) for item in value.tolist()]
    return [float(item) for item in value]
