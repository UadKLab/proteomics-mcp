from __future__ import annotations

from pathlib import Path

import pandas as pd


def read_candidate_table(path: str | Path, *, limit: int | None = None, sheet_name: str | int = 0) -> pd.DataFrame:
    table_path = Path(path)
    suffix = table_path.suffix.lower()
    if suffix in {".xlsx", ".xlsm", ".xls"}:
        return pd.read_excel(table_path, nrows=limit, sheet_name=sheet_name)
    if suffix == ".tsv":
        return pd.read_csv(table_path, sep="\t", nrows=limit)
    return pd.read_csv(table_path, nrows=limit)
