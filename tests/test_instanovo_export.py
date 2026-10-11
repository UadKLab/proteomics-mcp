from __future__ import annotations

import pandas as pd

from proteomics_mcp.core.instanovo_export import export_table_to_instanovo


def test_export_table_to_instanovo_ipc_uses_required_columns(tmp_path):
    input_path = tmp_path / "spectra.feather"
    pd.DataFrame(
        [
            {
                "source_name": "runA",
                "scan_number": 42,
                "spectrum_id": "runA:42",
                "precursor_mz": 500.2,
                "precursor_charge": 2,
                "mz_array": [100.0, 200.0],
                "intensity_array": [10.0, 20.0],
                "extra_field": "ignored",
            }
        ]
    ).to_feather(input_path)
    output_path = tmp_path / "instanovo.ipc"

    result = export_table_to_instanovo(input_path, output_path, output_format="ipc")
    exported = pd.read_feather(output_path)

    assert result["output_format"] == "ipc"
    assert result["sequence_included"] is False
    assert exported.columns.tolist() == [
        "experiment_name",
        "scan_number",
        "spectrum_id",
        "precursor_mz",
        "precursor_charge",
        "mz_array",
        "intensity_array",
    ]
    assert exported.loc[0, "experiment_name"] == "runA"
    assert exported.loc[0, "precursor_charge"] == 2


def test_export_table_to_instanovo_parquet_can_include_sequence(tmp_path):
    input_path = tmp_path / "spectra.feather"
    pd.DataFrame(
        [
            {
                "source_file": "/tmp/demo_run.mzML",
                "precursor_mz": 650.3,
                "precursor_charge": 3,
                "mz_array": [150.0, 250.0],
                "intensity_array": [30.0, 40.0],
                "target": "PEPTIDE",
            }
        ]
    ).to_feather(input_path)
    output_path = tmp_path / "instanovo.parquet"

    result = export_table_to_instanovo(input_path, output_path, output_format="parquet", sequence_column="target")
    exported = pd.read_parquet(output_path)

    assert result["output_format"] == "parquet"
    assert result["sequence_included"] is True
    assert exported.loc[0, "experiment_name"] == "demo_run"
    assert exported.loc[0, "sequence"] == "PEPTIDE"
    assert exported.loc[0, "precursor_charge"] == 3
