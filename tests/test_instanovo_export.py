from __future__ import annotations

import pandas as pd
import pytest

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


def _spectra(path, rows):
    defaults = {"precursor_mz": 500.0, "precursor_charge": 2, "mz_array": [100.0], "intensity_array": [1.0]}
    pd.DataFrame([{**defaults, **row} for row in rows]).to_feather(path)
    return path


def _read(path, output_format):
    return pd.read_feather(path) if output_format == "ipc" else pd.read_parquet(path)


@pytest.mark.parametrize("output_format", ["ipc", "parquet"])
def test_missing_sequences_stay_null(tmp_path, output_format):
    # As text, missing sequences became the peptide labels "None" and "nan" (#80).
    input_path = _spectra(
        tmp_path / "spectra.feather",
        [{"source_name": "runA", "scan_number": i, "target": target} for i, target in enumerate(["PEPTIDE", None, float("nan"), "  ", " LESLIEK "])],
    )
    output_path = tmp_path / f"instanovo.{output_format}"

    result = export_table_to_instanovo(input_path, output_path, output_format=output_format, sequence_column="target")
    exported = _read(output_path, output_format)

    assert exported["sequence"].tolist()[0] == "PEPTIDE"
    assert exported["sequence"].tolist()[4] == "LESLIEK"
    assert exported["sequence"].iloc[1:4].isna().all()
    assert not exported["sequence"].isin(["None", "nan", "<NA>", ""]).any()
    assert result["sequences_missing"] == 3


def test_missing_run_names_and_spectrum_ids_fall_back_per_row(tmp_path):
    input_path = _spectra(
        tmp_path / "spectra.feather",
        [
            {"experiment_name": "runA", "source_file": "/data/fileA.mzML", "scan_number": 1, "spectrum_id": "runA:1"},
            {"experiment_name": None, "source_file": "/data/fileB.mzML", "scan_number": 2, "spectrum_id": None},
        ],
    )
    output_path = tmp_path / "instanovo.ipc"

    export_table_to_instanovo(input_path, output_path)
    exported = pd.read_feather(output_path)

    assert exported["experiment_name"].tolist() == ["runA", "fileB"]
    assert exported["spectrum_id"].tolist() == ["runA:1", "fileB:2"]


def test_spectra_without_precursor_are_counted(tmp_path):
    input_path = _spectra(
        tmp_path / "spectra.feather",
        [{"source_name": "runA", "scan_number": 1}, {"source_name": "runA", "scan_number": 2, "precursor_charge": None}],
    )

    result = export_table_to_instanovo(input_path, tmp_path / "instanovo.ipc")

    assert result["spectra_written"] == 1
    assert result["spectra_dropped_without_precursor"] == 1
    assert result["sequences_missing"] is None


def test_drop_unlabelled_writes_only_spectra_with_a_sequence(tmp_path):
    input_path = _spectra(
        tmp_path / "spectra.feather",
        [{"source_name": "runA", "scan_number": i, "target": target} for i, target in enumerate(["PEPTIDE", None, " ", "LESLIEK"])],
    )
    output_path = tmp_path / "instanovo.parquet"

    result = export_table_to_instanovo(input_path, output_path, output_format="parquet", sequence_column="target", drop_unlabelled=True)
    exported = pd.read_parquet(output_path)

    assert exported["sequence"].tolist() == ["PEPTIDE", "LESLIEK"]
    assert exported["scan_number"].tolist() == [0, 3]
    assert result["spectra_written"] == 2
    assert result["spectra_dropped_unlabelled"] == 2
    assert result["sequences_missing"] == 0


def test_drop_unlabelled_needs_a_sequence_column(tmp_path):
    input_path = _spectra(tmp_path / "spectra.feather", [{"source_name": "runA", "scan_number": 1}])

    with pytest.raises(ValueError, match="sequence_column"):
        export_table_to_instanovo(input_path, tmp_path / "instanovo.ipc", drop_unlabelled=True)


def test_drop_unlabelled_is_exposed_over_mcp_and_cli(tmp_path, monkeypatch, capsys):
    import asyncio
    import json
    import sys

    from proteomics_mcp import cli
    from proteomics_mcp.server import mcp

    tools = asyncio.run(mcp.list_tools())
    schema = next(tool for tool in getattr(tools, "tools", tools) if tool.name == "export_to_instanovo")
    properties = (getattr(schema, "input_schema", None) or schema.inputSchema)["properties"]
    assert properties["drop_unlabelled"]["default"] is False

    input_path = _spectra(tmp_path / "spectra.feather", [{"source_name": "runA", "scan_number": i, "target": t} for i, t in enumerate(["PEPTIDE", None])])
    monkeypatch.setattr(
        sys,
        "argv",
        ["proteomics-mcp", "export-to-instanovo", "--input-path", str(input_path), "--output-path", str(tmp_path / "out.ipc"),
         "--sequence-column", "target", "--drop-unlabelled"],
    )
    cli.main()
    assert json.loads(capsys.readouterr().out)["spectra_dropped_unlabelled"] == 1
