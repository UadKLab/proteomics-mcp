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
    drop_unlabelled: bool = False,
) -> dict[str, Any]:
    """Write a spectrum table in InstaNovo's input layout.

    With sequence_column, spectra without a sequence get a null sequence; InstaNovo refuses such a table as annotated
    data. drop_unlabelled leaves those spectra out, for a labelled training or evaluation set.
    """
    if drop_unlabelled and not sequence_column:
        raise ValueError("drop_unlabelled needs sequence_column: without sequences, every spectrum is unlabelled.")
    frame = load_spectrum_table(input_path)
    export = _build_instanovo_frame(frame, sequence_column=sequence_column)
    dropped = len(frame) - len(export)
    dropped_unlabelled = 0
    if drop_unlabelled:
        labelled = export["sequence"].notna()
        dropped_unlabelled = int((~labelled).sum())
        export = export[labelled].reset_index(drop=True)
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
        # spectra without precursor m/z or charge, which InstaNovo needs
        "spectra_dropped_without_precursor": int(dropped),
        "spectra_dropped_unlabelled": dropped_unlabelled,
        "columns": export.columns.tolist(),
        "sequence_included": "sequence" in export.columns,
        # written as null, not as a sequence: unlabelled spectra in a labelled export
        "sequences_missing": int(export["sequence"].isna().sum()) if "sequence" in export.columns else None,
    }


def _build_instanovo_frame(frame: pd.DataFrame, *, sequence_column: str | None = None) -> pd.DataFrame:
    required = ["precursor_mz", "precursor_charge", "mz_array", "intensity_array"]
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"Input spectrum table is missing required columns for InstaNovo export: {missing}")

    export = pd.DataFrame(index=frame.index)
    # Run name per row: the first of these that is set (#80: str() turned missing values into "None"/"nan").
    experiment = pd.Series([None] * len(frame), index=frame.index, dtype=object)
    for column, transform in (("experiment_name", None), ("source_name", None), ("source_file", lambda value: Path(value).stem)):
        if column in frame.columns:
            values = _text(frame[column])
            if transform is not None:
                values = values.map(lambda value: transform(value) if value is not None else None)
            experiment = experiment.where(experiment.notna(), values)
    export["experiment_name"] = experiment.where(experiment.notna(), "spectra")

    if "scan_number" in frame.columns:
        export["scan_number"] = pd.to_numeric(frame["scan_number"], errors="coerce").astype("Int64")
    else:
        export["scan_number"] = pd.Series(range(len(frame)), index=frame.index, dtype="Int64")

    generated_ids = pd.Series(
        [
            f"{experiment}:{int(scan) if pd.notna(scan) else idx}"
            for idx, (experiment, scan) in enumerate(zip(export["experiment_name"], export["scan_number"]))
        ],
        index=frame.index,
    )
    if "spectrum_id" in frame.columns:
        spectrum_ids = _text(frame["spectrum_id"])
        export["spectrum_id"] = spectrum_ids.where(spectrum_ids.notna(), generated_ids)
    else:
        export["spectrum_id"] = generated_ids

    if "rt" in frame.columns:
        export["rt"] = pd.to_numeric(frame["rt"], errors="coerce")

    export["precursor_mz"] = pd.to_numeric(frame["precursor_mz"], errors="coerce")
    export["precursor_charge"] = pd.to_numeric(frame["precursor_charge"], errors="coerce").astype("Int64")
    export["mz_array"] = frame["mz_array"]
    export["intensity_array"] = frame["intensity_array"]

    if sequence_column:
        if sequence_column not in frame.columns:
            raise ValueError(f"sequence_column not found in input table: {sequence_column}")
        # Missing sequences stay null; as text they became the peptide labels "None" / "nan" (#80).
        export["sequence"] = _text(frame[sequence_column])

    export = export.dropna(subset=["precursor_mz", "precursor_charge"])
    export["precursor_charge"] = export["precursor_charge"].astype(int)
    export["mz_array"] = export["mz_array"].map(_coerce_float_list)
    export["intensity_array"] = export["intensity_array"].map(_coerce_float_list)
    return export.reset_index(drop=True)


def _text(values: pd.Series) -> pd.Series:
    """Values as stripped text, with missing (None, NaN, pd.NA) and empty values as None."""

    def convert(value: Any) -> str | None:
        if value is None or (not isinstance(value, (list, tuple)) and pd.isna(value)):
            return None
        text = str(value).strip()
        return text or None

    return values.map(convert).astype(object)


def _coerce_float_list(value: Any) -> list[float]:
    if value is None:
        return []
    if isinstance(value, list):
        return [float(item) for item in value]
    if hasattr(value, "tolist"):
        return [float(item) for item in value.tolist()]
    return [float(item) for item in value]
