from __future__ import annotations

import ast
import html as html_lib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from proteomics_mcp.adapters.table_io import read_candidate_table
from proteomics_mcp.core.annotation import parse_peptidoform
from proteomics_mcp.models.schemas import CandidatePSM, model_dump

try:  # pragma: no cover - exercised through report rendering integration tests
    import plotly.graph_objects as go
except ImportError:  # pragma: no cover
    go = None  # type: ignore


NUMERIC_COLUMNS = [
    "scan_number",
    "rank",
    "precursor_mz",
    "precursor_charge",
    "confidence",
    "raw_score",
    "next_beam",
    "next_beam_margin",
    "token_log_probability_mean",
    "token_log_probability_min",
    "token_confidence_mean",
    "token_confidence_min",
    "token_confidence_min_position",
    "spectral_angle",
    "cosine_similarity",
    "pearson_correlation",
    "xcorr",
    "log10_xcorr",
    "log2_xcorr",
    "hyperscore",
    "log2_hyperscore",
    "matched_ion_count",
    "total_predicted_ions",
    "matched_ion_fraction",
    "explained_intensity_fraction",
    "precursor_signal_fraction",
    "fragment_mass_error_mean_ppm",
    "fragment_mass_error_abs_mean_ppm",
    "fragment_mass_error_mean_da",
    "fragment_mass_error_abs_mean_da",
    "precursor_mass_error_ppm",
    "precursor_mass_error_da",
    "retention_time_seconds",
    "retention_time_minutes",
    "collision_energy",
    "faims_compensation_voltage",
    "injection_time",
    "scan_resolution",
    "annotated_peak_count",
    "annotated_intensity_fraction",
    "peptide_length",
]

METRIC_COLUMNS = [
    "spectral_angle",
    "cosine_similarity",
    "pearson_correlation",
    "log2_xcorr",
    "log2_hyperscore",
    "matched_ion_count",
    "matched_ion_fraction",
    "explained_intensity_fraction",
    "annotated_peak_count",
    "annotated_intensity_fraction",
]

SEMANTIC_FIELD_PATTERNS: dict[str, tuple[str, ...]] = {
    "confidence": ("confidence", "probability"),
    "p_value": ("pvalue", "p-value", "p_value"),
    "q_value": ("qvalue", "q-value", "q_value"),
    "pep": ("pep", "posterior error probability", "posterior_error_probability"),
    "retention_time": ("retention time", "retention_time", "rt"),
    "precursor_mz": ("precursor_mz", "observed m/z", "observed mz", "m/z [da]", "mz"),
    "quant": ("abundance", "area", "intensity", "quantity", "quant"),
}


@dataclass(slots=True)
class DatasetAggregation:
    frame: pd.DataFrame
    dataset_summary: dict[str, Any]
    peptide_summary: pd.DataFrame
    peptidoform_summary: pd.DataFrame
    run_summary: pd.DataFrame
    sample_summary: pd.DataFrame | None
    charge_summary: pd.DataFrame
    sdrf_summary: dict[str, Any] | None = None


@dataclass(slots=True)
class CandidateInspection:
    frame: pd.DataFrame
    peptide_summary: pd.DataFrame
    dataset_summary: dict[str, Any]
    run_summary: pd.DataFrame
    sample_summary: pd.DataFrame | None
    sdrf_summary: dict[str, Any] | None = None


ID_PSM_COLUMNS = [
    "candidate_id",
    "source_engine",
    "source_engine_version",
    "experiment_name",
    "run_key",
    "scan_number",
    "spectrum_id",
    "usi",
    "peptidoform",
    "stripped_sequence",
    "is_modified",
    "precursor_mz",
    "precursor_charge",
    "rank",
    "raw_score",
    "raw_score_name",
    "next_beam",
    "next_beam_margin",
    "confidence",
    "confidence_type",
    "token_log_probability_mean",
    "token_log_probability_min",
    "token_confidence_mean",
    "token_confidence_min",
    "token_confidence_min_position",
    "token_confidence_min_token",
    "prediction_id",
    "predictor_name",
    "predictor_model",
    "spectral_angle",
    "cosine_similarity",
    "pearson_correlation",
    "log2_xcorr",
    "log2_hyperscore",
    "matched_ion_count",
    "matched_ion_fraction",
    "explained_intensity_fraction",
    "precursor_signal_fraction",
    "annotated_peak_count",
    "annotated_intensity_fraction",
    "retention_time_minutes",
    "fragmentation_method",
    "instrument_model",
    "collision_energy",
    "scan_resolution",
]

ID_PEPTIDE_COLUMNS = [
    "stripped_sequence",
    "total_psm_count",
    "distinct_peptidoform_count",
    "peptidoforms",
    "modified_psm_count",
    "unmodified_psm_count",
    "distinct_modified_peptidoform_count",
    "best_peptidoform",
    "best_modified_peptidoform",
    "best_unmodified_peptidoform",
    "distinct_run_count",
    "per_file_psm_counts",
    "distinct_charge_count",
    "per_charge_counts",
    "best_candidate_id",
    "best_scan_number",
    "best_spectrum_id",
    "best_experiment_name",
    "best_confidence",
    "best_raw_score",
    "best_next_beam",
    "best_next_beam_margin",
    "best_token_log_probability_mean",
    "best_token_log_probability_min",
    "best_token_confidence_mean",
    "best_token_confidence_min",
    "best_token_confidence_min_position",
    "best_token_confidence_min_token",
    "best_spectral_angle",
    "best_cosine_similarity",
    "best_pearson_correlation",
    "best_log2_xcorr",
    "best_log2_hyperscore",
    "best_matched_ion_count",
    "best_matched_ion_fraction",
    "best_explained_intensity_fraction",
    "best_annotated_peak_count",
    "best_annotated_intensity_fraction",
    "spectral_angle_mean",
    "spectral_angle_max",
    "cosine_similarity_mean",
    "cosine_similarity_max",
    "pearson_correlation_mean",
    "pearson_correlation_max",
    "log2_xcorr_mean",
    "log2_xcorr_max",
    "log2_hyperscore_mean",
    "log2_hyperscore_max",
    "matched_ion_count_mean",
    "matched_ion_count_max",
    "matched_ion_fraction_mean",
    "matched_ion_fraction_max",
    "explained_intensity_fraction_mean",
    "explained_intensity_fraction_max",
    "annotated_peak_count_mean",
    "annotated_peak_count_max",
    "annotated_intensity_fraction_mean",
    "annotated_intensity_fraction_max",
    "fragmentation_methods",
    "instrument_counts",
    "detector_type_counts",
    "collision_energies",
    "observed_ion_series_counts",
    "source_engine_counts",
]

OPTIONAL_PEPTIDE_COLUMNS = {
    "best_spectral_angle",
    "best_cosine_similarity",
    "best_pearson_correlation",
    "best_log2_xcorr",
    "best_log2_hyperscore",
    "best_matched_ion_count",
    "best_matched_ion_fraction",
    "best_explained_intensity_fraction",
    "best_annotated_peak_count",
    "best_annotated_intensity_fraction",
    "spectral_angle_mean",
    "spectral_angle_max",
    "cosine_similarity_mean",
    "cosine_similarity_max",
    "pearson_correlation_mean",
    "pearson_correlation_max",
    "log2_xcorr_mean",
    "log2_xcorr_max",
    "log2_hyperscore_mean",
    "log2_hyperscore_max",
    "matched_ion_count_mean",
    "matched_ion_count_max",
    "matched_ion_fraction_mean",
    "matched_ion_fraction_max",
    "explained_intensity_fraction_mean",
    "explained_intensity_fraction_max",
    "annotated_peak_count_mean",
    "annotated_peak_count_max",
    "annotated_intensity_fraction_mean",
    "annotated_intensity_fraction_max",
    "fragmentation_methods",
    "instrument_counts",
    "detector_type_counts",
    "collision_energies",
    "observed_ion_series_counts",
}

SCORING_ENRICHMENT_COLUMNS = [
    "prediction_id",
    "predictor_name",
    "predictor_model",
    "spectral_angle",
    "cosine_similarity",
    "pearson_correlation",
    "xcorr",
    "log10_xcorr",
    "log2_xcorr",
    "hyperscore",
    "log2_hyperscore",
    "matched_ion_count",
    "total_predicted_ions",
    "matched_ion_fraction",
    "explained_intensity_fraction",
    "precursor_signal_fraction",
    "fragment_mass_error_mean_ppm",
    "fragment_mass_error_abs_mean_ppm",
    "fragment_mass_error_mean_da",
    "fragment_mass_error_abs_mean_da",
    "precursor_mass_error_ppm",
    "precursor_mass_error_da",
    "fragmentation_method",
    "isolation_target_mz",
    "retention_time_seconds",
    "retention_time_minutes",
    "collision_energy",
    "faims_compensation_voltage",
    "injection_time",
    "scan_header",
    "scan_resolution",
    "instrument_model",
    "detector_type",
    "annotated_peak_count",
    "ion_type_counts",
    "annotated_intensity_fraction",
    "metric_profile",
]


def _json_safe(value: Any) -> Any:
    if hasattr(value, "item"):
        try:
            return _json_safe(value.item())
        except Exception:
            pass
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    if isinstance(value, dict):
        return {str(k): _json_safe(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_json_safe(v) for v in value]
    return value


def _parse_mapping(value: Any) -> dict[str, Any]:
    if value is None:
        return {}
    if isinstance(value, dict):
        return {str(k): v for k, v in value.items()}
    if isinstance(value, float) and math.isnan(value):
        return {}
    text = str(value).strip()
    if not text:
        return {}
    for loader in (json.loads, ast.literal_eval):
        try:
            parsed = loader(text)
        except Exception:
            continue
        if isinstance(parsed, dict):
            return {str(k): v for k, v in parsed.items()}
    return {}


def _serialize_cell(value: Any) -> Any:
    if isinstance(value, (dict, list)):
        return json.dumps(_json_safe(value), sort_keys=True)
    return value


def _serialize_browse_table_cell(column: str, value: Any) -> Any:
    """Format table-display values without changing the underlying summaries."""
    if column == "peptidoforms" and isinstance(value, str):
        try:
            parsed = json.loads(value)
        except json.JSONDecodeError:
            parsed = value
        if isinstance(parsed, list):
            value = parsed
    if column == "peptidoforms" and isinstance(value, (list, tuple, set)):
        return ", ".join(str(item) for item in value if item is not None)
    return _serialize_cell(value)


def _browse_table_html(frame: pd.DataFrame) -> str:
    columns = [str(column) for column in frame.columns]
    header_cells = "".join(
        f'<th scope="col" data-column-index="{index}">{html_lib.escape(column)}</th>'
        for index, column in enumerate(columns)
    )
    rows: list[str] = []
    for values in frame.itertuples(index=False, name=None):
        cells = []
        for index, value in enumerate(values):
            if value is None or (isinstance(value, float) and math.isnan(value)):
                display = ""
            else:
                display = str(value)
            cells.append(f'<td data-column-index="{index}">{html_lib.escape(display)}</td>')
        rows.append(f"<tr>{''.join(cells)}</tr>")
    return f'<table class="browse-table"><thead><tr>{header_cells}</tr></thead><tbody>{"".join(rows)}</tbody></table>'


def _source_column_value(row: pd.Series, key: str) -> Any:
    mapping = row.get("source_columns")
    if isinstance(mapping, dict):
        return mapping.get(key)
    return None


def _source_column_keys(frame: pd.DataFrame) -> list[str]:
    keys: set[str] = set()
    if "source_columns" not in frame.columns:
        return []
    for value in frame["source_columns"].tolist():
        if isinstance(value, dict):
            keys.update(str(k) for k in value.keys())
    return sorted(keys)


def _resolve_field_name(frame: pd.DataFrame, field: str) -> str:
    if field in frame.columns:
        return field
    if field.startswith("source_columns."):
        return field
    source_keys = _source_column_keys(frame)
    if field in source_keys:
        return f"source_columns.{field}"
    if field in SEMANTIC_FIELD_PATTERNS:
        catalog = distribution_field_catalog(frame)
        matches = catalog["semantic_hints"].get(field) or []
        if matches:
            return str(matches[0])
    return field


def _extract_field_series(frame: pd.DataFrame, field: str) -> pd.Series:
    field = _resolve_field_name(frame, field)
    if field in frame.columns:
        return frame[field]
    source_key = field[len("source_columns.") :] if field.startswith("source_columns.") else field
    if "source_columns" in frame.columns:
        return frame.apply(lambda row: _source_column_value(row, source_key), axis=1)
    raise KeyError(f"Unknown field '{field}'.")


def _numeric_series(frame: pd.DataFrame, field: str) -> pd.Series:
    return pd.to_numeric(_extract_field_series(frame, field), errors="coerce")


def _string_series(frame: pd.DataFrame, field: str) -> pd.Series:
    series = _extract_field_series(frame, field)
    return series.map(_coerce_string)


def _weighted_quantile(values: pd.Series, weights: pd.Series, quantile: float) -> float | None:
    if values.empty:
        return None
    ordered = pd.DataFrame({"value": values, "weight": weights}).sort_values("value", kind="mergesort")
    cumulative = ordered["weight"].cumsum()
    total = float(ordered["weight"].sum())
    if total <= 0.0:
        return None
    threshold = quantile * total
    row = ordered[cumulative >= threshold].head(1)
    if row.empty:
        return float(ordered["value"].iloc[-1])
    return float(row["value"].iloc[0])


def distribution_field_catalog(frame: pd.DataFrame) -> dict[str, Any]:
    numeric_fields: list[str] = []
    categorical_fields: list[str] = []

    for column in frame.columns:
        if column.startswith("__") or column == "source_columns":
            continue
        numeric = pd.to_numeric(frame[column], errors="coerce")
        if numeric.notna().any():
            numeric_fields.append(column)
        elif frame[column].map(lambda value: _coerce_string(value) is not None).any():
            categorical_fields.append(column)

    source_numeric_fields: list[str] = []
    source_categorical_fields: list[str] = []
    for key in _source_column_keys(frame):
        series = frame.apply(lambda row: _source_column_value(row, key), axis=1)
        numeric = pd.to_numeric(series, errors="coerce")
        if numeric.notna().any():
            source_numeric_fields.append(f"source_columns.{key}")
        elif series.map(lambda value: _coerce_string(value) is not None).any():
            source_categorical_fields.append(f"source_columns.{key}")

    semantic_hints: dict[str, list[str]] = {}
    all_fields = numeric_fields + source_numeric_fields + categorical_fields + source_categorical_fields
    lowered = {field: field.lower() for field in all_fields}
    for semantic, patterns in SEMANTIC_FIELD_PATTERNS.items():
        matches = [field for field, lower in lowered.items() if any(pattern in lower for pattern in patterns)]
        if matches:
            semantic_hints[semantic] = sorted(matches)

    return {
        "numeric_fields": sorted(numeric_fields),
        "source_numeric_fields": sorted(source_numeric_fields),
        "categorical_fields": sorted(categorical_fields),
        "source_categorical_fields": sorted(source_categorical_fields),
        "semantic_hints": semantic_hints,
        "semantic_defaults": {semantic: matches[0] for semantic, matches in semantic_hints.items() if matches},
    }


def summarize_numeric_distribution(
    frame: pd.DataFrame,
    *,
    field: str,
    group_by: str | None = None,
    percentiles: list[float] | None = None,
    weight_field: str | None = None,
) -> dict[str, Any]:
    percentile_values = percentiles or [0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99]

    def summarize(series: pd.Series) -> dict[str, Any]:
        numeric = pd.to_numeric(series, errors="coerce")
        valid = numeric.dropna()
        summary = {
            "field": field,
            "count": int(len(valid)),
            "missing_count": int(len(numeric) - len(valid)),
            "distinct_value_count": int(valid.nunique()),
            "min": float(valid.min()) if not valid.empty else None,
            "max": float(valid.max()) if not valid.empty else None,
            "mean": float(valid.mean()) if not valid.empty else None,
            "median": float(valid.median()) if not valid.empty else None,
            "std": float(valid.std(ddof=1)) if len(valid) > 1 else None,
            "percentiles": {str(p): (float(valid.quantile(p)) if not valid.empty else None) for p in percentile_values},
        }
        if weight_field is not None:
            weights = _numeric_series(frame.loc[series.index], weight_field)
            valid_mask = numeric.notna() & weights.notna()
            weighted_values = numeric[valid_mask]
            weighted_weights = weights[valid_mask]
            positive_mask = weighted_weights > 0
            weighted_values = weighted_values[positive_mask]
            weighted_weights = weighted_weights[positive_mask]
            total_weight = float(weighted_weights.sum()) if not weighted_weights.empty else None
            weighted_mean = None
            weighted_percentiles: dict[str, float | None] = {}
            if total_weight is not None and total_weight > 0.0:
                weighted_mean = float((weighted_values * weighted_weights).sum() / total_weight)
                weighted_percentiles = {str(p): _weighted_quantile(weighted_values, weighted_weights, p) for p in percentile_values}
            summary["weight_field"] = weight_field
            summary["weighted_count"] = int(len(weighted_values))
            summary["total_weight"] = total_weight
            summary["weighted_mean"] = weighted_mean
            summary["weighted_percentiles"] = weighted_percentiles
        return _json_safe(summary)

    if group_by is None:
        return summarize(_numeric_series(frame, field))

    group_labels = _string_series(frame, group_by)
    numeric = _numeric_series(frame, field)
    grouped: dict[str, Any] = {}
    for label in sorted({value for value in group_labels.tolist() if value is not None}):
        mask = group_labels == label
        grouped[label] = summarize(numeric[mask])
    result = {"field": field, "group_by": group_by, "groups": grouped}
    if weight_field is not None:
        result["weight_field"] = weight_field
    return result


def bin_numeric_distribution(
    frame: pd.DataFrame,
    *,
    field: str,
    bin_edges: list[float] | None = None,
    bin_count: int | None = None,
    value_range: tuple[float, float] | None = None,
    group_by: str | None = None,
    weight_field: str | None = None,
) -> dict[str, Any]:
    numeric = _numeric_series(frame, field)
    valid = numeric.dropna()
    if valid.empty:
        result = {"field": field, "bin_edges": [], "bins": []}
        return result if group_by is None else {"field": field, "group_by": group_by, "bin_edges": [], "groups": {}}

    if bin_edges is None:
        bins = int(bin_count or 10)
        lower = float(valid.min()) if value_range is None else float(value_range[0])
        upper = float(valid.max()) if value_range is None else float(value_range[1])
        if math.isclose(lower, upper):
            upper = lower + 1.0
        edges = [lower + ((upper - lower) * i / bins) for i in range(bins + 1)]
    else:
        edges = [float(value) for value in bin_edges]
        if len(edges) < 2:
            raise ValueError("bin_edges must contain at least two numeric edges.")

    def summarize_bins(series: pd.Series) -> list[dict[str, Any]]:
        valid_series = pd.to_numeric(series, errors="coerce").dropna()
        if valid_series.empty:
            return [
                {"left": float(edges[i]), "right": float(edges[i + 1]), "count": 0, "weight_sum": None if weight_field is None else 0.0}
                for i in range(len(edges) - 1)
            ]
        categories = pd.cut(valid_series, bins=edges, include_lowest=True, right=False)
        counts = categories.value_counts(sort=False)
        weight_sums: dict[Any, float] = {}
        if weight_field is not None:
            weights = _numeric_series(frame.loc[series.index], weight_field)
            bin_frame = pd.DataFrame({"bin": categories, "weight": weights.loc[valid_series.index]}).dropna(subset=["weight"])
            if not bin_frame.empty:
                grouped_weights = bin_frame.groupby("bin", observed=False)["weight"].sum()
                weight_sums = {interval: float(value) for interval, value in grouped_weights.items()}
        bins_out: list[dict[str, Any]] = []
        for interval, count in counts.items():
            row = {"left": float(interval.left), "right": float(interval.right), "count": int(count)}
            if weight_field is not None:
                row["weight_sum"] = weight_sums.get(interval, 0.0)
            bins_out.append(row)
        return bins_out

    if group_by is None:
        result = {"field": field, "bin_edges": edges, "bins": summarize_bins(numeric)}
        if weight_field is not None:
            result["weight_field"] = weight_field
        return result

    group_labels = _string_series(frame, group_by)
    grouped: dict[str, Any] = {}
    for label in sorted({value for value in group_labels.tolist() if value is not None}):
        mask = group_labels == label
        grouped[label] = summarize_bins(numeric[mask])
    result = {"field": field, "group_by": group_by, "bin_edges": edges, "groups": grouped}
    if weight_field is not None:
        result["weight_field"] = weight_field
    return result


def summarize_categorical_distribution(
    frame: pd.DataFrame,
    *,
    field: str,
    group_by: str | None = None,
    top_n: int | None = None,
) -> dict[str, Any]:
    def summarize(series: pd.Series) -> dict[str, Any]:
        values = series.map(_coerce_string).dropna()
        counts = values.value_counts(sort=True)
        if top_n is not None:
            counts = counts.head(top_n)
        total = int(values.shape[0])
        return {
            "field": field,
            "count": total,
            "missing_count": int(series.shape[0] - total),
            "distinct_value_count": int(values.nunique()),
            "counts": {str(index): int(value) for index, value in counts.items()},
            "fractions": {str(index): (float(value) / total if total else None) for index, value in counts.items()},
        }

    series = _extract_field_series(frame, field)
    if group_by is None:
        return _json_safe(summarize(series))

    group_labels = _string_series(frame, group_by)
    grouped: dict[str, Any] = {}
    for label in sorted({value for value in group_labels.tolist() if value is not None}):
        mask = group_labels == label
        grouped[label] = summarize(series[mask])
    return {"field": field, "group_by": group_by, "groups": _json_safe(grouped)}


def summarize_missing_values(
    frame: pd.DataFrame,
    *,
    fields: list[str] | None = None,
) -> dict[str, Any]:
    selected_fields = fields or sorted(distribution_field_catalog(frame)["numeric_fields"] + distribution_field_catalog(frame)["source_numeric_fields"])
    rows: list[dict[str, Any]] = []
    total_rows = int(len(frame))
    for field in selected_fields:
        series = _extract_field_series(frame, field)
        numeric = pd.to_numeric(series, errors="coerce")
        missing_count = int(numeric.isna().sum()) if numeric.notna().any() else int(series.map(lambda value: _coerce_string(value) is None).sum())
        rows.append(
            {
                "field": field,
                "total_count": total_rows,
                "missing_count": missing_count,
                "missing_fraction": (missing_count / total_rows) if total_rows else None,
            }
        )
    return {"total_rows": total_rows, "fields": rows}


def summarize_cv_by_category(
    frame: pd.DataFrame,
    *,
    value_field: str,
    category_field: str,
    top_n: int | None = None,
) -> dict[str, Any]:
    values = _numeric_series(frame, value_field)
    categories = _string_series(frame, category_field)
    rows: list[dict[str, Any]] = []
    for category in sorted({value for value in categories.tolist() if value is not None}):
        mask = categories == category
        series = values[mask].dropna()
        mean = float(series.mean()) if not series.empty else None
        std = float(series.std(ddof=1)) if len(series) > 1 else None
        cv = (std / mean) if mean not in (None, 0.0) and std is not None else None
        rows.append(
            {
                "category": category,
                "count": int(len(series)),
                "missing_count": int(mask.sum() - len(series)),
                "mean": mean,
                "std": std,
                "cv": cv,
            }
        )
    rows.sort(key=lambda row: (-row["count"], row["category"]))
    if top_n is not None:
        rows = rows[:top_n]
    return {"value_field": value_field, "category_field": category_field, "rows": _json_safe(rows)}


def bin_2d_numeric_distribution(
    frame: pd.DataFrame,
    *,
    x_field: str,
    y_field: str,
    x_bin_edges: list[float] | None = None,
    y_bin_edges: list[float] | None = None,
    x_bin_count: int | None = None,
    y_bin_count: int | None = None,
    weight_field: str | None = None,
) -> dict[str, Any]:
    x_series = _numeric_series(frame, x_field)
    y_series = _numeric_series(frame, y_field)
    valid_mask = x_series.notna() & y_series.notna()
    x_valid = x_series[valid_mask]
    y_valid = y_series[valid_mask]
    if x_valid.empty or y_valid.empty:
        return {"x_field": x_field, "y_field": y_field, "x_bin_edges": [], "y_bin_edges": [], "cells": []}

    def resolve_edges(series: pd.Series, edges: list[float] | None, count: int | None) -> list[float]:
        if edges is not None:
            out = [float(value) for value in edges]
            if len(out) < 2:
                raise ValueError("2D bin edges must contain at least two numeric edges.")
            return out
        bins = int(count or 10)
        lower = float(series.min())
        upper = float(series.max())
        if math.isclose(lower, upper):
            upper = lower + 1.0
        return [lower + ((upper - lower) * i / bins) for i in range(bins + 1)]

    x_edges = resolve_edges(x_valid, x_bin_edges, x_bin_count)
    y_edges = resolve_edges(y_valid, y_bin_edges, y_bin_count)

    x_bins = pd.cut(x_valid, bins=x_edges, include_lowest=True, right=False)
    y_bins = pd.cut(y_valid, bins=y_edges, include_lowest=True, right=False)
    bin_frame = pd.DataFrame({"x_bin": x_bins, "y_bin": y_bins})
    if weight_field is not None:
        weights = _numeric_series(frame.loc[valid_mask], weight_field).fillna(0.0)
        bin_frame["weight"] = weights.values

    cell_counts = bin_frame.groupby(["x_bin", "y_bin"], observed=False).size()
    cell_weights = bin_frame.groupby(["x_bin", "y_bin"], observed=False)["weight"].sum() if weight_field is not None else None

    cells: list[dict[str, Any]] = []
    for (x_bin, y_bin), count in cell_counts.items():
        row = {
            "x_left": float(x_bin.left),
            "x_right": float(x_bin.right),
            "y_left": float(y_bin.left),
            "y_right": float(y_bin.right),
            "count": int(count),
        }
        if cell_weights is not None:
            row["weight_sum"] = float(cell_weights.loc[(x_bin, y_bin)])
        cells.append(row)

    result = {"x_field": x_field, "y_field": y_field, "x_bin_edges": x_edges, "y_bin_edges": y_edges, "cells": cells}
    if weight_field is not None:
        result["weight_field"] = weight_field
    return result


def scatter_2d_distribution(
    frame: pd.DataFrame,
    *,
    x_field: str,
    y_field: str,
    color_field: str | None = "run_key",
    hover_fields: list[str] | None = None,
    max_points: int = 5000,
) -> dict[str, Any]:
    x_series = _numeric_series(frame, x_field)
    y_series = _numeric_series(frame, y_field)
    valid_mask = x_series.notna() & y_series.notna()
    if not valid_mask.any():
        return {"x_field": x_field, "y_field": y_field, "point_count": 0, "points": []}

    plot_frame = pd.DataFrame(
        {
            "x": x_series[valid_mask],
            "y": y_series[valid_mask],
        },
        index=frame.index[valid_mask],
    )
    if color_field:
        plot_frame["group"] = _string_series(frame.loc[plot_frame.index], color_field)
    for field in hover_fields or []:
        plot_frame[field] = _extract_field_series(frame.loc[plot_frame.index], field)

    if len(plot_frame) > max_points:
        step = max(1, math.ceil(len(plot_frame) / max_points))
        plot_frame = plot_frame.iloc[::step].head(max_points).copy()

    points: list[dict[str, Any]] = []
    for _, row in plot_frame.iterrows():
        point: dict[str, Any] = {
            "x": float(row["x"]),
            "y": float(row["y"]),
        }
        if "group" in plot_frame.columns:
            point["group"] = _coerce_string(row.get("group"))
        for field in hover_fields or []:
            point[field] = _json_safe(row.get(field))
        points.append(point)

    result = {
        "x_field": x_field,
        "y_field": y_field,
        "point_count": int(len(points)),
        "points": points,
    }
    if color_field:
        result["color_field"] = color_field
    return _json_safe(result)


def build_dashboard_summary(
    frame: pd.DataFrame,
    *,
    rt_field: str | None = None,
    mz_field: str = "precursor_mz",
    confidence_field: str | None = "confidence",
    quant_field: str | None = None,
    category_field: str = "experiment_name",
) -> dict[str, Any]:
    catalog = distribution_field_catalog(frame)
    rt_field_resolved = rt_field
    if rt_field_resolved is None:
        candidates = catalog["semantic_hints"].get("retention_time", [])
        rt_field_resolved = candidates[0] if candidates else None
    confidence_field_resolved = confidence_field
    if confidence_field_resolved is None:
        candidates = catalog["semantic_hints"].get("confidence", [])
        confidence_field_resolved = candidates[0] if candidates else None

    summary: dict[str, Any] = {
        "overview": {
            "total_psm_count": int(len(frame)),
            "distinct_peptide_count": int(frame["stripped_sequence"].nunique()) if "stripped_sequence" in frame.columns else None,
            "distinct_peptidoform_count": int(frame["peptidoform"].nunique()) if "peptidoform" in frame.columns else None,
            **_distinct_precursor_counts(frame),
            "distinct_run_count": int(frame["run_key"].dropna().nunique()) if "run_key" in frame.columns else None,
        },
        "ids_per_file": summarize_categorical_distribution(frame, field=category_field),
        "missing_values": summarize_missing_values(frame, fields=[field for field in [mz_field, rt_field_resolved, confidence_field_resolved, quant_field] if field]),
        "histograms": {},
        "cv_summaries": {},
        "grids_2d": {},
        "scatter_2d": {},
    }

    summary["histograms"]["precursor_mz"] = bin_numeric_distribution(frame, field=mz_field, bin_count=20)
    peptide_lengths = pd.to_numeric(frame.get("peptide_length", pd.Series(dtype=float)), errors="coerce").dropna()
    if not peptide_lengths.empty:
        min_len = int(peptide_lengths.min())
        max_len = int(peptide_lengths.max())
        summary["histograms"]["peptide_length"] = bin_numeric_distribution(
            frame,
            field="peptide_length",
            bin_edges=[float(value) for value in range(min_len, max_len + 2)],
        )
    if rt_field_resolved:
        summary["histograms"]["retention_time"] = bin_numeric_distribution(frame, field=rt_field_resolved, bin_count=20)
    if confidence_field_resolved:
        summary["histograms"]["confidence"] = bin_numeric_distribution(frame, field=confidence_field_resolved, bin_count=20)
    if "spectral_angle" in frame.columns and pd.to_numeric(frame["spectral_angle"], errors="coerce").notna().any():
        summary["histograms"]["spectral_angle"] = bin_numeric_distribution(frame, field="spectral_angle", bin_count=20)
    if "matched_ion_fraction" in frame.columns and pd.to_numeric(frame["matched_ion_fraction"], errors="coerce").notna().any():
        summary["histograms"]["matched_ion_fraction"] = bin_numeric_distribution(frame, field="matched_ion_fraction", bin_count=20)
    if "explained_intensity_fraction" in frame.columns and pd.to_numeric(frame["explained_intensity_fraction"], errors="coerce").notna().any():
        summary["histograms"]["explained_intensity_fraction"] = bin_numeric_distribution(frame, field="explained_intensity_fraction", bin_count=20)
    if "log2_xcorr" in frame.columns and pd.to_numeric(frame["log2_xcorr"], errors="coerce").notna().any():
        summary["histograms"]["log2_xcorr"] = bin_numeric_distribution(frame, field="log2_xcorr", bin_count=20)
    summary["charge_distribution"] = summarize_categorical_distribution(frame, field="precursor_charge")
    if quant_field:
        summary["histograms"]["quant"] = bin_numeric_distribution(frame, field=quant_field, bin_count=20)
        summary["cv_summaries"]["quant_by_category"] = summarize_cv_by_category(frame, value_field=quant_field, category_field=category_field)

    if confidence_field_resolved:
        summary["grids_2d"]["confidence_vs_mz"] = bin_2d_numeric_distribution(frame, x_field=mz_field, y_field=confidence_field_resolved, x_bin_count=20, y_bin_count=20, weight_field=quant_field)
        if rt_field_resolved:
            summary["grids_2d"]["confidence_vs_rt"] = bin_2d_numeric_distribution(frame, x_field=rt_field_resolved, y_field=confidence_field_resolved, x_bin_count=20, y_bin_count=20, weight_field=quant_field)
        if "spectral_angle" in frame.columns and pd.to_numeric(frame["spectral_angle"], errors="coerce").notna().any():
            summary["grids_2d"]["confidence_vs_spectral_angle"] = bin_2d_numeric_distribution(
                frame,
                x_field=confidence_field_resolved,
                y_field="spectral_angle",
                x_bin_count=20,
                y_bin_count=20,
                weight_field=quant_field,
            )
    if rt_field_resolved:
        summary["grids_2d"]["mz_vs_rt"] = bin_2d_numeric_distribution(frame, x_field=mz_field, y_field=rt_field_resolved, x_bin_count=20, y_bin_count=20, weight_field=quant_field)

    scatter_hover_fields = ["candidate_id", "peptidoform", "run_key", "scan_number"]
    if confidence_field_resolved and "spectral_angle" in frame.columns and pd.to_numeric(frame["spectral_angle"], errors="coerce").notna().any():
        summary["scatter_2d"]["confidence_vs_spectral_angle"] = scatter_2d_distribution(
            frame,
            x_field=confidence_field_resolved,
            y_field="spectral_angle",
            hover_fields=scatter_hover_fields,
        )
    if "spectral_angle" in frame.columns and "pearson_correlation" in frame.columns:
        if pd.to_numeric(frame["spectral_angle"], errors="coerce").notna().any() and pd.to_numeric(frame["pearson_correlation"], errors="coerce").notna().any():
            summary["scatter_2d"]["spectral_angle_vs_pearson"] = scatter_2d_distribution(
                frame,
                x_field="spectral_angle",
                y_field="pearson_correlation",
                hover_fields=scatter_hover_fields,
            )
    if "spectral_angle" in frame.columns and "log2_xcorr" in frame.columns:
        if pd.to_numeric(frame["spectral_angle"], errors="coerce").notna().any() and pd.to_numeric(frame["log2_xcorr"], errors="coerce").notna().any():
            summary["scatter_2d"]["spectral_angle_vs_log2_xcorr"] = scatter_2d_distribution(
                frame,
                x_field="spectral_angle",
                y_field="log2_xcorr",
                hover_fields=scatter_hover_fields,
            )
    if confidence_field_resolved and "matched_ion_count" in frame.columns:
        if pd.to_numeric(frame["matched_ion_count"], errors="coerce").notna().any():
            summary["scatter_2d"]["confidence_vs_matched_ion_count"] = scatter_2d_distribution(
                frame,
                x_field=confidence_field_resolved,
                y_field="matched_ion_count",
                hover_fields=scatter_hover_fields,
            )
    if "spectral_angle" in frame.columns and "precursor_mass_error_ppm" in frame.columns:
        if pd.to_numeric(frame["precursor_mass_error_ppm"], errors="coerce").notna().any():
            summary["scatter_2d"]["spectral_angle_vs_precursor_mass_error_ppm"] = scatter_2d_distribution(
                frame,
                x_field="spectral_angle",
                y_field="precursor_mass_error_ppm",
                hover_fields=scatter_hover_fields,
            )
    if "spectral_angle" in frame.columns and "fragment_mass_error_abs_mean_ppm" in frame.columns:
        if pd.to_numeric(frame["fragment_mass_error_abs_mean_ppm"], errors="coerce").notna().any():
            summary["scatter_2d"]["spectral_angle_vs_fragment_mass_error_abs_mean_ppm"] = scatter_2d_distribution(
                frame,
                x_field="spectral_angle",
                y_field="fragment_mass_error_abs_mean_ppm",
                hover_fields=scatter_hover_fields,
            )

    return _json_safe(summary)


def _entity_series(frame: pd.DataFrame, entity: str) -> pd.Series:
    if entity == "peptide":
        return frame["stripped_sequence"].astype(str)
    if entity == "peptidoform":
        return frame["peptidoform"].astype(str)
    if entity == "peptide_precursor":
        charge = pd.to_numeric(frame.get("precursor_charge", pd.Series(dtype=float)), errors="coerce")
        return frame["stripped_sequence"].astype(str) + "|" + charge.map(lambda value: str(int(value)) if pd.notna(value) else "NA")
    if entity == "peptidoform_precursor":
        charge = pd.to_numeric(frame.get("precursor_charge", pd.Series(dtype=float)), errors="coerce")
        return frame["peptidoform"].astype(str) + "|" + charge.map(lambda value: str(int(value)) if pd.notna(value) else "NA")
    raise ValueError(f"Unsupported overlap entity '{entity}'.")


def summarize_overlap(
    frame: pd.DataFrame,
    *,
    entity: str = "peptide",
    group_field: str = "run_key",
) -> dict[str, Any]:
    groups = _string_series(frame, group_field)
    entities = _entity_series(frame, entity)
    membership: dict[str, set[str]] = {}
    for label in sorted({value for value in groups.tolist() if value is not None}):
        mask = groups == label
        membership[label] = {str(value) for value in entities[mask].tolist() if _coerce_string(value)}

    group_rows: list[dict[str, Any]] = []
    all_groups = list(membership.keys())
    for label in all_groups:
        current = membership[label]
        others = set().union(*(membership[other] for other in all_groups if other != label)) if len(all_groups) > 1 else set()
        group_rows.append(
            {
                "group": label,
                "entity_count": int(len(current)),
                "unique_count": int(len(current - others)),
                "shared_count": int(len(current & others)),
            }
        )

    pairwise_rows: list[dict[str, Any]] = []
    for i, left in enumerate(all_groups):
        for right in all_groups[i + 1 :]:
            left_set = membership[left]
            right_set = membership[right]
            intersection = left_set & right_set
            union = left_set | right_set
            pairwise_rows.append(
                {
                    "group_a": left,
                    "group_b": right,
                    "intersection_count": int(len(intersection)),
                    "union_count": int(len(union)),
                    "jaccard": (len(intersection) / len(union)) if union else None,
                    "containment_a": (len(intersection) / len(left_set)) if left_set else None,
                    "containment_b": (len(intersection) / len(right_set)) if right_set else None,
                }
            )

    return {
        "entity": entity,
        "group_field": group_field,
        "group_count": int(len(all_groups)),
        "total_distinct_entities": int(len(set().union(*membership.values()))) if membership else 0,
        "groups": _json_safe(group_rows),
        "pairwise": _json_safe(pairwise_rows),
    }


def _modification_entries(peptidoform: str) -> list[dict[str, Any]]:
    text = str(peptidoform or "")
    entries: list[dict[str, Any]] = []
    residue_index = 0
    i = 0
    while i < len(text):
        char = text[i]
        if char.isalpha() and char.isupper():
            residue_index += 1
            residue = char
            i += 1
            while i < len(text) and text[i] == "[":
                end = text.find("]", i + 1)
                if end == -1:
                    break
                mod = text[i + 1 : end]
                entries.append({"residue": residue, "position": residue_index, "modification": mod, "site": f"{residue}{residue_index}"})
                i = end + 1
            continue
        if char == "[" and residue_index == 0:
            end = text.find("]", i + 1)
            if end == -1:
                break
            mod = text[i + 1 : end]
            entries.append({"residue": "N-term", "position": 0, "modification": mod, "site": "N-term"})
            i = end + 1
            continue
        i += 1
    return entries


def summarize_modifications(frame: pd.DataFrame, *, group_field: str = "run_key") -> dict[str, Any]:
    modification_rows: list[dict[str, Any]] = []
    for _, row in frame.iterrows():
        peptidoform = _coerce_string(row.get("peptidoform"))
        if not peptidoform:
            continue
        group = _coerce_string(_extract_field_series(frame.loc[[row.name]], group_field).iloc[0]) if group_field else None
        for entry in _modification_entries(peptidoform):
            modification_rows.append({**entry, "group": group})

    mod_frame = pd.DataFrame(modification_rows)
    run_groups = _string_series(frame, group_field) if group_field else pd.Series(dtype=object)
    modified_mask = frame["is_modified"] if "is_modified" in frame.columns else pd.Series([False] * len(frame))

    summary: dict[str, Any] = {
        "modified_psm_count": int(modified_mask.sum()) if not modified_mask.empty else 0,
        "unmodified_psm_count": int((~modified_mask).sum()) if not modified_mask.empty else int(len(frame)),
        "modified_fraction": (float(modified_mask.mean()) if not modified_mask.empty else None),
        "distinct_modified_peptidoform_count": int(frame.loc[modified_mask, "peptidoform"].nunique()) if "peptidoform" in frame.columns and not modified_mask.empty else 0,
        "modification_counts": {},
        "residue_counts": {},
        "group_modification_counts": {},
    }
    if not mod_frame.empty:
        summary["modification_counts"] = {str(k): int(v) for k, v in mod_frame["modification"].value_counts().items()}
        residue_frame = mod_frame[~mod_frame["residue"].isin(["N-term", "C-term"])]
        summary["residue_counts"] = {str(k): int(v) for k, v in residue_frame["residue"].value_counts().items()} if not residue_frame.empty else {}
        if "group" in mod_frame.columns:
            group_counts: dict[str, dict[str, int]] = {}
            for group, group_df in mod_frame.groupby("group", dropna=True, sort=True):
                label = _coerce_string(group)
                if not label:
                    continue
                group_counts[label] = {str(k): int(v) for k, v in group_df["modification"].value_counts().items()}
            summary["group_modification_counts"] = group_counts
    if group_field:
        group_modified: dict[str, int] = {}
        for label in sorted({value for value in run_groups.tolist() if value is not None}):
            mask = run_groups == label
            group_modified[label] = int(modified_mask[mask].sum())
        summary["group_modified_psm_counts"] = group_modified
    return _json_safe(summary)


_UNIMOD_TAG_RE = re.compile(r"\[UNIMOD:\d+\]|\[\+?-?\d+(?:\.\d+)?\]")
_NTERM_MOD_RE = re.compile(r"^(?:\[UNIMOD:\d+\]|\[\+?-?\d+(?:\.\d+)?\])-")


def _strip_modifications(sequence: str | None) -> str:
    if not isinstance(sequence, str):
        return ""
    return _UNIMOD_TAG_RE.sub("", _NTERM_MOD_RE.sub("", sequence))


def summarize_prediction_accuracy(
    frame: pd.DataFrame,
    *,
    target_field: str,
    confidence_field: str | None = "confidence",
    mass_error_field: str | None = None,
    confidence_thresholds: list[float] | None = None,
    reliability_bin_count: int = 10,
) -> dict[str, Any]:
    """Compare each PSM's predicted peptidoform to a ground-truth target column.

    Returns three correctness levels (exact, stripped sequence, I/L-equivalent),
    a reliability diagram, a coverage-vs-precision sweep, accuracy by charge /
    predicted-length / target-modification-state, and confidence/mass-error
    distributions split by correctness. Output is a JSON-safe dict suitable for
    rendering or programmatic consumption.
    """
    if "peptidoform" not in frame.columns:
        raise ValueError("Frame is missing 'peptidoform' column.")

    pred_raw = frame["peptidoform"].map(lambda v: v if isinstance(v, str) else None)
    try:
        target_raw = _extract_field_series(frame, target_field).map(
            lambda v: v if isinstance(v, str) else None
        )
    except KeyError as exc:
        raise ValueError(f"Target field '{target_field}' is not present on the frame.") from exc

    valid_mask = pred_raw.notna() & target_raw.notna()
    work = pd.DataFrame(
        {
            "prediction": pred_raw,
            "target": target_raw,
            "pred_stripped": pred_raw.map(_strip_modifications),
            "target_stripped": target_raw.map(_strip_modifications),
        }
    )[valid_mask].copy()
    work["pred_il"] = work["pred_stripped"].str.replace("I", "L", regex=False)
    work["target_il"] = work["target_stripped"].str.replace("I", "L", regex=False)
    work["correct_exact"] = work["prediction"] == work["target"]
    work["correct_seq"] = work["pred_stripped"] == work["target_stripped"]
    work["correct_il"] = work["pred_il"] == work["target_il"]
    work["pred_length"] = work["pred_stripped"].str.len()

    base_frame = frame.loc[work.index]
    if confidence_field:
        confidence = pd.to_numeric(_extract_field_series(base_frame, confidence_field), errors="coerce")
    else:
        confidence = pd.Series([float("nan")] * len(work), index=work.index)
    work["confidence"] = confidence

    work["target_has_mod"] = work["target"].str.contains(r"\[", regex=True, na=False)

    if "precursor_charge" in base_frame.columns:
        work["precursor_charge"] = pd.to_numeric(base_frame["precursor_charge"], errors="coerce")
    else:
        work["precursor_charge"] = pd.Series([None] * len(work), index=work.index)

    n = int(len(work))
    overview = {
        "psm_count": n,
        "skipped_missing_pred_or_target": int((~valid_mask).sum()),
        "exact_accuracy": float(work["correct_exact"].mean()) if n else None,
        "stripped_sequence_accuracy": float(work["correct_seq"].mean()) if n else None,
        "il_equivalent_accuracy": float(work["correct_il"].mean()) if n else None,
    }

    reliability_bins: list[dict[str, Any]] = []
    if n and confidence.notna().any():
        edges = np.linspace(0.0, 1.0, max(reliability_bin_count, 2) + 1)
        binned = work.assign(
            _bin=pd.cut(work["confidence"], bins=edges, include_lowest=True),
        )
        for interval, group in binned.groupby("_bin", observed=True):
            if interval is None or pd.isna(interval):
                continue
            reliability_bins.append(
                {
                    "left": float(interval.left),
                    "right": float(interval.right),
                    "count": int(len(group)),
                    "accuracy_il": float(group["correct_il"].mean()),
                    "accuracy_seq": float(group["correct_seq"].mean()),
                    "accuracy_exact": float(group["correct_exact"].mean()),
                    "mean_confidence": float(group["confidence"].mean()),
                }
            )

    thresholds = confidence_thresholds or [0.99, 0.95, 0.9, 0.8, 0.5, 0.25, 0.1, 0.0]
    threshold_rows: list[dict[str, Any]] = []
    if n:
        for thr in thresholds:
            mask = work["confidence"] >= thr if confidence.notna().any() else pd.Series([True] * n, index=work.index)
            kept = int(mask.sum())
            row: dict[str, Any] = {
                "threshold": float(thr),
                "kept": kept,
                "coverage": float(kept) / n if n else 0.0,
            }
            if kept:
                row["accuracy_il"] = float(work.loc[mask, "correct_il"].mean())
                row["accuracy_seq"] = float(work.loc[mask, "correct_seq"].mean())
                row["accuracy_exact"] = float(work.loc[mask, "correct_exact"].mean())
            else:
                row["accuracy_il"] = None
                row["accuracy_seq"] = None
                row["accuracy_exact"] = None
            threshold_rows.append(row)

    by_charge: list[dict[str, Any]] = []
    if "precursor_charge" in work.columns and work["precursor_charge"].notna().any():
        for charge_value, group in work.dropna(subset=["precursor_charge"]).groupby("precursor_charge"):
            by_charge.append(
                {
                    "charge": int(charge_value),
                    "count": int(len(group)),
                    "accuracy_il": float(group["correct_il"].mean()),
                    "accuracy_seq": float(group["correct_seq"].mean()),
                    "accuracy_exact": float(group["correct_exact"].mean()),
                }
            )
        by_charge.sort(key=lambda row: row["charge"])

    by_length: list[dict[str, Any]] = []
    if work["pred_length"].notna().any():
        for length_value, group in work.groupby("pred_length"):
            by_length.append(
                {
                    "length": int(length_value),
                    "count": int(len(group)),
                    "accuracy_il": float(group["correct_il"].mean()),
                }
            )
        by_length.sort(key=lambda row: row["length"])

    by_mod_state: list[dict[str, Any]] = []
    for label, group in work.groupby("target_has_mod"):
        by_mod_state.append(
            {
                "target_has_mod": bool(label),
                "label": "modified" if label else "unmodified",
                "count": int(len(group)),
                "accuracy_il": float(group["correct_il"].mean()),
                "accuracy_seq": float(group["correct_seq"].mean()),
                "accuracy_exact": float(group["correct_exact"].mean()),
            }
        )
    by_mod_state.sort(key=lambda row: row["label"])

    confidence_distribution: dict[str, Any] = {}
    if confidence.notna().any():
        edges = np.linspace(0.0, 1.0, 41)
        for label, mask in [("correct_il", work["correct_il"]), ("wrong_il", ~work["correct_il"])]:
            counts, _ = np.histogram(work.loc[mask, "confidence"].dropna(), bins=edges)
            confidence_distribution[label] = {
                "bins": [
                    {"left": float(edges[idx]), "right": float(edges[idx + 1]), "count": int(counts[idx])}
                    for idx in range(len(counts))
                ]
            }

    mass_error_distribution: dict[str, Any] = {}
    if mass_error_field:
        try:
            mass_error_series = pd.to_numeric(_extract_field_series(base_frame, mass_error_field), errors="coerce")
        except KeyError:
            mass_error_series = pd.Series([float("nan")] * n, index=work.index)
        if mass_error_series.notna().any():
            clipped = mass_error_series.clip(-5.0, 5.0)
            edges = np.linspace(-5.0, 5.0, 61)
            for label, mask in [("correct_il", work["correct_il"]), ("wrong_il", ~work["correct_il"])]:
                counts, _ = np.histogram(clipped.loc[mask].dropna(), bins=edges)
                mass_error_distribution[label] = {
                    "bins": [
                        {"left": float(edges[idx]), "right": float(edges[idx + 1]), "count": int(counts[idx])}
                        for idx in range(len(counts))
                    ]
                }

    return _json_safe(
        {
            "overview": overview,
            "reliability": {"bins": reliability_bins},
            "operating_points": threshold_rows,
            "by_charge": by_charge,
            "by_predicted_length": by_length,
            "by_target_modification_state": by_mod_state,
            "confidence_distribution_by_correctness": confidence_distribution,
            "mass_error_distribution_by_correctness": mass_error_distribution,
        }
    )


def summarize_quant_quality(
    frame: pd.DataFrame,
    *,
    quant_field: str,
    entity: str = "peptide",
    group_field: str = "run_key",
    min_groups_for_cv: int = 2,
) -> dict[str, Any]:
    entity_values = _entity_series(frame, entity)
    group_values = _string_series(frame, group_field)
    quant_values = _numeric_series(frame, quant_field)
    valid = pd.DataFrame({"entity": entity_values, "group": group_values, "quant": quant_values}).dropna(subset=["group"])
    if valid.empty:
        return {"entity": entity, "group_field": group_field, "quant_field": quant_field, "overview": {}, "entity_rows": []}

    grouped = valid.groupby(["entity", "group"], observed=False)["quant"].mean().unstack("group")
    total_groups = int(grouped.shape[1])
    entity_rows: list[dict[str, Any]] = []
    for entity_key, row in grouped.iterrows():
        non_null = row.dropna()
        mean = float(non_null.mean()) if not non_null.empty else None
        std = float(non_null.std(ddof=1)) if len(non_null) > 1 else None
        cv = (std / mean) if std is not None and mean not in (None, 0.0) and len(non_null) >= min_groups_for_cv else None
        entity_rows.append(
            {
                "entity": str(entity_key),
                "observed_group_count": int(non_null.shape[0]),
                "missing_group_count": int(total_groups - non_null.shape[0]),
                "completeness_fraction": (float(non_null.shape[0]) / total_groups) if total_groups else None,
                "mean_quant": mean,
                "std_quant": std,
                "cv": cv,
            }
        )
    entity_rows.sort(key=lambda row: (-row["observed_group_count"], row["entity"]))
    complete = [row for row in entity_rows if row["observed_group_count"] == total_groups]
    cv_rows = [row for row in entity_rows if row["cv"] is not None]
    top_variable = sorted(cv_rows, key=lambda row: (-float(row["cv"]), row["entity"]))[:10]
    top_stable = sorted(cv_rows, key=lambda row: (float(row["cv"]), row["entity"]))[:10]

    return _json_safe(
        {
            "entity": entity,
            "group_field": group_field,
            "quant_field": quant_field,
            "overview": {
                "entity_count": int(len(entity_rows)),
                "group_count": total_groups,
                "complete_entity_count": int(len(complete)),
                "complete_entity_fraction": (len(complete) / len(entity_rows)) if entity_rows else None,
                "entities_with_cv": int(len(cv_rows)),
            },
            "entity_rows": entity_rows,
            "top_variable": top_variable,
            "top_stable": top_stable,
        }
    )


def summarize_file_coverage(
    frame: pd.DataFrame,
    *,
    group_field: str = "run_key",
    mz_field: str = "precursor_mz",
    rt_field: str | None = None,
    confidence_field: str | None = None,
    quant_field: str | None = None,
) -> dict[str, Any]:
    catalog = distribution_field_catalog(frame)
    rt_resolved = _resolve_field_name(frame, rt_field) if rt_field else catalog["semantic_defaults"].get("retention_time")
    confidence_resolved = _resolve_field_name(frame, confidence_field) if confidence_field else catalog["semantic_defaults"].get("confidence")
    group_values = _string_series(frame, group_field)
    rows: dict[str, Any] = {}
    for label in sorted({value for value in group_values.tolist() if value is not None}):
        mask = group_values == label
        group_frame = frame[mask].copy()
        rows[label] = {
            "overview": {
                "total_psm_count": int(len(group_frame)),
                "distinct_peptide_count": int(group_frame["stripped_sequence"].nunique()),
                "distinct_peptidoform_count": int(group_frame["peptidoform"].nunique()),
                **_distinct_precursor_counts(group_frame),
                "modified_psm_count": int(group_frame["is_modified"].sum()),
            },
            "charge_counts": summarize_categorical_distribution(group_frame, field="precursor_charge")["counts"],
            "m/z_summary": summarize_numeric_distribution(group_frame, field=mz_field),
            "m/z_histogram": bin_numeric_distribution(group_frame, field=mz_field, bin_count=20),
            "peptide_length_histogram": bin_numeric_distribution(group_frame, field="peptide_length", bin_count=20),
            "missing_values": summarize_missing_values(group_frame, fields=[field for field in [mz_field, rt_resolved, confidence_resolved, quant_field] if field]),
        }
        if rt_resolved:
            rows[label]["rt_summary"] = summarize_numeric_distribution(group_frame, field=rt_resolved)
            rows[label]["rt_histogram"] = bin_numeric_distribution(group_frame, field=rt_resolved, bin_count=20)
        if confidence_resolved:
            rows[label]["confidence_summary"] = summarize_numeric_distribution(group_frame, field=confidence_resolved)
            rows[label]["confidence_histogram"] = bin_numeric_distribution(group_frame, field=confidence_resolved, bin_count=20)
        if rt_resolved:
            rows[label]["mz_vs_rt"] = bin_2d_numeric_distribution(group_frame, x_field=mz_field, y_field=rt_resolved, x_bin_count=20, y_bin_count=20, weight_field=quant_field)
    return _json_safe({"group_field": group_field, "files": rows})


def _stripped_sequence(peptidoform: str) -> str:
    parsed = parse_peptidoform(str(peptidoform))
    return "".join(parsed.residues)


def _peptide_length_from_sequence(sequence: str) -> int:
    return len(str(sequence or ""))


def _is_modified(peptidoform: str, stripped_sequence: str) -> bool:
    cleaned = str(peptidoform or "").strip()
    return bool("[" in cleaned or "-" in cleaned or cleaned != stripped_sequence)


def _coerce_string(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, float) and math.isnan(value):
        return None
    text = str(value).strip()
    return text or None


def _run_key(row: pd.Series) -> str | None:
    experiment_name = _coerce_string(row.get("experiment_name"))
    if experiment_name:
        return experiment_name
    spectrum_id = _coerce_string(row.get("spectrum_id"))
    if spectrum_id and ":" in spectrum_id:
        return spectrum_id.split(":", 1)[0]
    return spectrum_id


def load_scored_psm_table(path: str | Path) -> pd.DataFrame:
    frame = read_candidate_table(path).copy()
    frame["__row_order"] = range(len(frame))
    for column in NUMERIC_COLUMNS:
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    if "ion_type_counts" in frame.columns:
        frame["ion_type_counts"] = frame["ion_type_counts"].apply(_parse_mapping)
    else:
        frame["ion_type_counts"] = [{} for _ in range(len(frame))]
    if "peptidoform" not in frame.columns:
        raise ValueError("Scored PSM table must include a 'peptidoform' column.")
    frame["peptidoform"] = frame["peptidoform"].astype(str)
    frame["stripped_sequence"] = frame["peptidoform"].map(_stripped_sequence)
    frame["peptide_length"] = frame["stripped_sequence"].map(_peptide_length_from_sequence)
    frame["is_modified"] = [
        _is_modified(peptidoform, stripped)
        for peptidoform, stripped in zip(frame["peptidoform"], frame["stripped_sequence"], strict=False)
    ]
    frame["run_key"] = frame.apply(_run_key, axis=1)
    frame["source_engine"] = frame.get("source_engine", pd.Series(["unknown"] * len(frame))).fillna("unknown").astype(str)
    return frame


_FALLBACK_JOIN_KEYS = ("run_key", "scan_number", "peptidoform")


def merge_scored_enrichment(frame: pd.DataFrame, scored_frame: pd.DataFrame) -> pd.DataFrame:
    return join_scored_enrichment(frame, scored_frame)[0]


def _fallback_key(row: Any) -> tuple[str, int, str] | None:
    run_key, scan, peptidoform = (row.get(column) for column in _FALLBACK_JOIN_KEYS)
    if run_key is None or pd.isna(run_key) or scan is None or pd.isna(scan) or not peptidoform:
        return None
    return str(run_key), int(scan), str(peptidoform)


def join_scored_enrichment(frame: pd.DataFrame, scored_frame: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Merge scored metrics onto candidate rows and report how the rows were matched.

    Rows join on candidate_id. Most adapters derive candidate_id from the input file name, so a renamed or
    copied table changes every id; scored rows whose id is not found then join on run + scan + peptidoform,
    when that key identifies exactly one candidate row.
    """
    summary = {
        "scored_rows": int(len(scored_frame)),
        "candidate_rows": int(len(frame)),
        "matched_by_candidate_id": 0,
        "matched_by_run_scan_peptidoform": 0,
        "unmatched_scored_rows": int(len(scored_frame)),
    }
    if frame.empty or scored_frame.empty:
        return frame.copy(), summary
    if "candidate_id" not in frame.columns or "candidate_id" not in scored_frame.columns:
        raise ValueError("Scored enrichment requires 'candidate_id' in both frames.")

    enriched = frame.copy()
    scored = scored_frame.copy()
    scored["candidate_id"] = scored["candidate_id"].astype(str)
    enriched["candidate_id"] = enriched["candidate_id"].astype(str)

    candidate_ids = set(enriched["candidate_id"])
    by_id = scored["candidate_id"].isin(candidate_ids)
    summary["matched_by_candidate_id"] = int(by_id.sum())
    if not by_id.all() and all(column in enriched.columns and column in scored.columns for column in _FALLBACK_JOIN_KEYS):
        keys = pd.Series([_fallback_key(row) for row in enriched.to_dict("records")], index=enriched.index)
        unique = keys.dropna()
        unique = unique[~unique.duplicated(keep=False)]
        id_by_key = dict(zip(unique, enriched.loc[unique.index, "candidate_id"]))
        claimed = set(scored.loc[by_id, "candidate_id"])
        for index, row in zip(scored.index[~by_id], scored.loc[~by_id].to_dict("records")):
            candidate_id = id_by_key.get(_fallback_key(row))
            if candidate_id is not None and candidate_id not in claimed:
                scored.at[index, "candidate_id"] = candidate_id
                claimed.add(candidate_id)
                summary["matched_by_run_scan_peptidoform"] += 1
    summary["unmatched_scored_rows"] = summary["scored_rows"] - summary["matched_by_candidate_id"] - summary["matched_by_run_scan_peptidoform"]

    available_columns = [column for column in SCORING_ENRICHMENT_COLUMNS if column in scored.columns]
    if not available_columns:
        return enriched, summary

    scored_subset = scored[["candidate_id", *available_columns]].drop_duplicates(subset=["candidate_id"], keep="last")
    merged = enriched.merge(scored_subset, on="candidate_id", how="left", suffixes=("", "__scored"))

    for column in available_columns:
        scored_column = f"{column}__scored"
        if scored_column not in merged.columns:
            continue
        if column in enriched.columns:
            merged[column] = merged[column].where(merged[column].notna(), merged[scored_column])
            merged.drop(columns=[scored_column], inplace=True)
        else:
            merged.rename(columns={scored_column: column}, inplace=True)

    for column in NUMERIC_COLUMNS:
        if column in merged.columns:
            merged[column] = pd.to_numeric(merged[column], errors="coerce")
    if "ion_type_counts" in merged.columns:
        merged["ion_type_counts"] = merged["ion_type_counts"].apply(_parse_mapping)
    return merged, summary


def candidate_models_to_frame(candidates: list[CandidatePSM]) -> pd.DataFrame:
    frame = pd.DataFrame([model_dump(candidate) for candidate in candidates]).copy()
    if frame.empty:
        frame = pd.DataFrame(columns=["candidate_id", "source_engine", "peptidoform"])
    frame["__row_order"] = range(len(frame))
    for column in NUMERIC_COLUMNS:
        if column in frame.columns:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    if "peptidoform" not in frame.columns:
        raise ValueError("Candidate table must include a 'peptidoform' column.")
    frame["peptidoform"] = frame["peptidoform"].astype(str)
    frame["stripped_sequence"] = frame["peptidoform"].map(_stripped_sequence)
    frame["peptide_length"] = frame["stripped_sequence"].map(_peptide_length_from_sequence)
    frame["is_modified"] = [
        _is_modified(peptidoform, stripped)
        for peptidoform, stripped in zip(frame["peptidoform"], frame["stripped_sequence"], strict=False)
    ]
    frame["run_key"] = frame.apply(_run_key, axis=1)
    frame["source_engine"] = frame.get("source_engine", pd.Series(["unknown"] * len(frame))).fillna("unknown").astype(str)
    return frame


def _sorted_counts(values: list[Any]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for value in values:
        key = _coerce_string(value)
        if key is None:
            continue
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items(), key=lambda item: (item[0])))


def _best_row(group: pd.DataFrame) -> pd.Series:
    sortable = group.copy()
    sortable["__has_confidence"] = sortable["confidence"].notna().astype(int) if "confidence" in sortable.columns else 0
    sortable["__confidence"] = sortable["confidence"].fillna(float("-inf")) if "confidence" in sortable.columns else float("-inf")
    sortable["__has_raw_score"] = sortable["raw_score"].notna().astype(int) if "raw_score" in sortable.columns else 0
    sortable["__raw_score"] = sortable["raw_score"].fillna(float("-inf")) if "raw_score" in sortable.columns else float("-inf")
    sortable["__spectral_angle"] = sortable["spectral_angle"].fillna(float("-inf")) if "spectral_angle" in sortable.columns else float("-inf")
    sortable["__matched_ion_count"] = sortable["matched_ion_count"].fillna(float("-inf")) if "matched_ion_count" in sortable.columns else float("-inf")
    sortable["__candidate_id"] = sortable.get("candidate_id", pd.Series([""] * len(sortable))).fillna("").astype(str)
    ordered = sortable.sort_values(
        by=[
            "__has_confidence",
            "__confidence",
            "__has_raw_score",
            "__raw_score",
            "__spectral_angle",
            "__matched_ion_count",
            "__row_order",
            "__candidate_id",
        ],
        ascending=[False, False, False, False, False, False, True, True],
        kind="mergesort",
    )
    return ordered.iloc[0]


def _metric_summary(group: pd.DataFrame) -> dict[str, Any]:
    summary: dict[str, Any] = {}
    for column in METRIC_COLUMNS:
        if column not in group.columns:
            continue
        series = pd.to_numeric(group[column], errors="coerce").dropna()
        if series.empty:
            summary[f"{column}_mean"] = None
            summary[f"{column}_median"] = None
            summary[f"{column}_max"] = None
            continue
        summary[f"{column}_mean"] = float(series.mean())
        summary[f"{column}_median"] = float(series.median())
        summary[f"{column}_max"] = float(series.max())
    return summary


def _top_ion_counts(group: pd.DataFrame) -> dict[str, int]:
    counts: dict[str, int] = {}
    if "ion_type_counts" not in group.columns:
        return counts
    for mapping in group["ion_type_counts"].tolist():
        if not isinstance(mapping, dict):
            continue
        for key, value in mapping.items():
            try:
                counts[str(key)] = counts.get(str(key), 0) + int(value)
            except Exception:
                continue
    return dict(sorted(counts.items(), key=lambda item: (-item[1], item[0])))


def _group_sample_metadata(group: pd.DataFrame) -> dict[str, Any]:
    metadata_rows = [value for value in group.get("__sample_metadata", pd.Series(dtype=object)).tolist() if isinstance(value, dict)]
    if not metadata_rows:
        return {}
    first = metadata_rows[0]
    return {
        "assay_name": first.get("assay_name"),
        "source_name": first.get("source_name"),
        "label": first.get("label"),
        "instrument": first.get("instrument"),
        "organism": first.get("organism"),
        "organism_part": first.get("organism_part"),
        "disease": first.get("disease"),
        "factor_values": first.get("factor_values"),
    }


def _summarize_group(group: pd.DataFrame, *, include_peptide_context: bool) -> dict[str, Any]:
    best = _best_row(group)
    runs = [value for value in group["run_key"].tolist() if _coerce_string(value)]
    charges = [int(value) for value in pd.to_numeric(group.get("precursor_charge", pd.Series(dtype=float)), errors="coerce").dropna().tolist()]
    source_engines = _sorted_counts(group.get("source_engine", pd.Series(dtype=object)).tolist())
    fragmentation_methods = _sorted_counts(group.get("fragmentation_method", pd.Series(dtype=object)).tolist())
    instruments = _sorted_counts(group.get("instrument_model", pd.Series(dtype=object)).tolist())
    detector_types = _sorted_counts(group.get("detector_type", pd.Series(dtype=object)).tolist())
    collision_energies = sorted({round(float(value), 4) for value in pd.to_numeric(group.get("collision_energy", pd.Series(dtype=float)), errors="coerce").dropna().tolist()})

    summary = {
        "total_psm_count": int(len(group)),
        "distinct_run_count": int(len(set(runs))),
        "per_file_psm_counts": _sorted_counts(runs),
        "distinct_charge_count": int(len(set(charges))),
        "per_charge_counts": {str(key): value for key, value in sorted(_sorted_counts(charges).items(), key=lambda item: int(item[0]))},
        "best_candidate_id": _coerce_string(best.get("candidate_id")),
        "best_scan_number": int(best["scan_number"]) if pd.notna(best.get("scan_number")) else None,
        "best_spectrum_id": _coerce_string(best.get("spectrum_id")),
        "best_experiment_name": _coerce_string(best.get("experiment_name")) or _coerce_string(best.get("run_key")),
        "best_confidence": float(best["confidence"]) if pd.notna(best.get("confidence")) else None,
        "best_raw_score": float(best["raw_score"]) if pd.notna(best.get("raw_score")) else None,
        "best_next_beam": float(best["next_beam"]) if pd.notna(best.get("next_beam")) else None,
        "best_next_beam_margin": float(best["next_beam_margin"]) if pd.notna(best.get("next_beam_margin")) else None,
        "best_token_log_probability_mean": float(best["token_log_probability_mean"]) if pd.notna(best.get("token_log_probability_mean")) else None,
        "best_token_log_probability_min": float(best["token_log_probability_min"]) if pd.notna(best.get("token_log_probability_min")) else None,
        "best_token_confidence_mean": float(best["token_confidence_mean"]) if pd.notna(best.get("token_confidence_mean")) else None,
        "best_token_confidence_min": float(best["token_confidence_min"]) if pd.notna(best.get("token_confidence_min")) else None,
        "best_token_confidence_min_position": int(best["token_confidence_min_position"]) if pd.notna(best.get("token_confidence_min_position")) else None,
        "best_token_confidence_min_token": _coerce_string(best.get("token_confidence_min_token")),
        "best_spectral_angle": float(best["spectral_angle"]) if pd.notna(best.get("spectral_angle")) else None,
        "best_cosine_similarity": float(best["cosine_similarity"]) if pd.notna(best.get("cosine_similarity")) else None,
        "best_pearson_correlation": float(best["pearson_correlation"]) if pd.notna(best.get("pearson_correlation")) else None,
        "best_log2_xcorr": float(best["log2_xcorr"]) if pd.notna(best.get("log2_xcorr")) else None,
        "best_log2_hyperscore": float(best["log2_hyperscore"]) if pd.notna(best.get("log2_hyperscore")) else None,
        "best_matched_ion_count": int(best["matched_ion_count"]) if pd.notna(best.get("matched_ion_count")) else None,
        "best_matched_ion_fraction": float(best["matched_ion_fraction"]) if pd.notna(best.get("matched_ion_fraction")) else None,
        "best_explained_intensity_fraction": float(best["explained_intensity_fraction"]) if pd.notna(best.get("explained_intensity_fraction")) else None,
        "best_annotated_peak_count": int(best["annotated_peak_count"]) if pd.notna(best.get("annotated_peak_count")) else None,
        "best_annotated_intensity_fraction": float(best["annotated_intensity_fraction"]) if pd.notna(best.get("annotated_intensity_fraction")) else None,
        "source_engine_counts": source_engines,
        "fragmentation_methods": fragmentation_methods,
        "instrument_counts": instruments,
        "detector_type_counts": detector_types,
        "collision_energies": collision_energies,
        "observed_ion_series_counts": _top_ion_counts(group),
    }
    summary.update(_metric_summary(group))
    if include_peptide_context:
        peptidoforms = sorted(set(group["peptidoform"].tolist()))
        modified = group[group["is_modified"]]
        unmodified = group[~group["is_modified"]]
        summary.update(
            {
                "distinct_peptidoform_count": int(group["peptidoform"].nunique()),
                "peptidoforms": peptidoforms,
                "modified_psm_count": int(group["is_modified"].sum()),
                "unmodified_psm_count": int((~group["is_modified"]).sum()),
                "distinct_modified_peptidoform_count": int(modified["peptidoform"].nunique()) if not modified.empty else 0,
                "best_peptidoform": _coerce_string(best.get("peptidoform")),
                "best_modified_peptidoform": _coerce_string(_best_row(modified).get("peptidoform")) if not modified.empty else None,
                "best_unmodified_peptidoform": _coerce_string(_best_row(unmodified).get("peptidoform")) if not unmodified.empty else None,
            }
        )
    return summary


def _enrich_with_sdrf(frame: pd.DataFrame, sdrf_summary: dict[str, Any] | None) -> pd.DataFrame:
    enriched = frame.copy()
    file_metadata = (sdrf_summary or {}).get("file_metadata") or {}

    def sample_metadata(run_key: Any) -> dict[str, Any] | None:
        key = _coerce_string(run_key)
        if not key:
            return None
        return file_metadata.get(key) or file_metadata.get(Path(key).name) or file_metadata.get(Path(key).stem)

    enriched["__sample_metadata"] = enriched["run_key"].map(sample_metadata)
    enriched["__sample_key"] = enriched["__sample_metadata"].map(
        lambda metadata: metadata.get("assay_name") or metadata.get("source_name") if isinstance(metadata, dict) else None
    )
    return enriched


def _frame_from_records(records: list[dict[str, Any]], sort_by: str) -> pd.DataFrame:
    frame = pd.DataFrame(records)
    if sort_by in frame.columns:
        frame = frame.sort_values(by=[sort_by], ascending=[False], kind="mergesort").reset_index(drop=True)
    return frame


def _distinct_precursor_counts(frame: pd.DataFrame) -> dict[str, int]:
    valid_charge = pd.to_numeric(frame.get("precursor_charge", pd.Series(dtype=float)), errors="coerce")
    precursor_frame = frame.assign(__charge=valid_charge).dropna(subset=["__charge"])
    peptide_precursor_count = 0
    peptidoform_precursor_count = 0
    if not precursor_frame.empty:
        peptide_precursor_count = int(
            precursor_frame[["stripped_sequence", "__charge"]].drop_duplicates().shape[0]
        )
        peptidoform_precursor_count = int(
            precursor_frame[["peptidoform", "__charge"]].drop_duplicates().shape[0]
        )
    return {
        "distinct_peptide_precursor_count": peptide_precursor_count,
        "distinct_peptidoform_precursor_count": peptidoform_precursor_count,
    }


def _summarize_runs(frame: pd.DataFrame) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for run_key, group in frame.groupby("run_key", dropna=True, sort=True):
        if not _coerce_string(run_key):
            continue
        row = {
            "experiment_name": str(run_key),
            "total_psm_count": int(len(group)),
            "distinct_peptide_count": int(group["stripped_sequence"].nunique()),
            "distinct_peptidoform_count": int(group["peptidoform"].nunique()),
        }
        row.update(_distinct_precursor_counts(group))
        row.update(_metric_summary(group))
        row["source_engine_counts"] = _sorted_counts(group["source_engine"].tolist())
        row["per_charge_counts"] = _sorted_counts(pd.to_numeric(group.get("precursor_charge", pd.Series(dtype=float)), errors="coerce").dropna().astype(int).tolist())
        row.update(_group_sample_metadata(group))
        records.append(row)
    return _frame_from_records(records, "total_psm_count")


def _summarize_samples(frame: pd.DataFrame) -> pd.DataFrame | None:
    if "__sample_key" not in frame.columns:
        return None
    records: list[dict[str, Any]] = []
    for sample_key, group in frame.groupby("__sample_key", dropna=True, sort=True):
        if not _coerce_string(sample_key):
            continue
        row = {
            "sample_key": str(sample_key),
            "total_psm_count": int(len(group)),
            "distinct_run_count": int(group["run_key"].nunique()),
            "distinct_peptide_count": int(group["stripped_sequence"].nunique()),
            "distinct_peptidoform_count": int(group["peptidoform"].nunique()),
            "runs": sorted({_coerce_string(value) for value in group["run_key"].tolist() if _coerce_string(value)}),
        }
        row.update(_distinct_precursor_counts(group))
        row.update(_metric_summary(group))
        row.update(_group_sample_metadata(group))
        records.append(row)
    if not records:
        return None
    return _frame_from_records(records, "total_psm_count")


def _summarize_charge_long(frame: pd.DataFrame) -> pd.DataFrame:
    if "precursor_charge" not in frame.columns:
        return pd.DataFrame()
    rows: list[dict[str, Any]] = []
    for (sequence, charge), group in frame.groupby(["stripped_sequence", "precursor_charge"], dropna=True, sort=True):
        rows.append(
            {
                "stripped_sequence": sequence,
                "precursor_charge": int(charge),
                "total_psm_count": int(len(group)),
                "distinct_peptidoform_count": int(group["peptidoform"].nunique()),
                "per_file_psm_counts": _sorted_counts(group["run_key"].tolist()),
            }
        )
    return pd.DataFrame(rows).sort_values(by=["stripped_sequence", "precursor_charge"], kind="mergesort").reset_index(drop=True) if rows else pd.DataFrame()


def summarize_dataset_frame(frame: pd.DataFrame, *, sdrf_summary: dict[str, Any] | None = None) -> DatasetAggregation:
    enriched = _enrich_with_sdrf(frame, sdrf_summary)

    peptidoform_rows: list[dict[str, Any]] = []
    for peptidoform, group in enriched.groupby("peptidoform", sort=True):
        row = {"peptidoform": peptidoform, "stripped_sequence": _coerce_string(group["stripped_sequence"].iloc[0])}
        row.update(_summarize_group(group, include_peptide_context=False))
        peptidoform_rows.append(row)
    peptidoform_summary = _frame_from_records(peptidoform_rows, "total_psm_count")

    peptide_rows: list[dict[str, Any]] = []
    for sequence, group in enriched.groupby("stripped_sequence", sort=True):
        row = {"stripped_sequence": sequence}
        row.update(_summarize_group(group, include_peptide_context=True))
        peptide_rows.append(row)
    peptide_summary = _frame_from_records(peptide_rows, "total_psm_count")

    run_summary = _summarize_runs(enriched)
    sample_summary = _summarize_samples(enriched)
    charge_summary = _summarize_charge_long(enriched)

    confidence_series = pd.to_numeric(enriched.get("confidence"), errors="coerce").dropna()
    spectral_series = pd.to_numeric(enriched.get("spectral_angle"), errors="coerce").dropna()
    modified_count = int(enriched["is_modified"].sum())
    precursor_counts = _distinct_precursor_counts(enriched)
    dataset_summary = {
        "dataset_metadata": {
            "total_psm_count": int(len(enriched)),
            "distinct_peptide_count": int(enriched["stripped_sequence"].nunique()),
            "distinct_peptidoform_count": int(enriched["peptidoform"].nunique()),
            **precursor_counts,
            "distinct_run_count": int(enriched["run_key"].dropna().nunique()),
            "source_engines": sorted(set(enriched["source_engine"].tolist())),
            "metric_columns": [column for column in METRIC_COLUMNS if column in enriched.columns],
        },
        "identification_summary": {
            "psms_per_run_mean": float(run_summary["total_psm_count"].mean()) if not run_summary.empty else None,
            "psms_per_run_max": int(run_summary["total_psm_count"].max()) if not run_summary.empty else None,
            "psms_per_peptide_mean": float(peptide_summary["total_psm_count"].mean()) if not peptide_summary.empty else None,
            "psms_per_peptide_max": int(peptide_summary["total_psm_count"].max()) if not peptide_summary.empty else None,
        },
        "sequence_summary": {
            "peptides_with_multiple_peptidoforms": int((peptide_summary.get("distinct_peptidoform_count", pd.Series(dtype=float)) > 1).sum()) if not peptide_summary.empty else 0,
            "modified_psm_count": modified_count,
            "unmodified_psm_count": int(len(enriched) - modified_count),
        },
        "peptidoform_summary": {
            "top_confident_peptidoforms": peptidoform_summary.sort_values(by=["best_confidence", "total_psm_count"], ascending=[False, False], kind="mergesort").head(10).to_dict(orient="records") if not peptidoform_summary.empty else [],
        },
        "run_summary": {
            "runs": run_summary.to_dict(orient="records"),
        },
        "charge_summary": {
            "distinct_charge_values": sorted({int(value) for value in pd.to_numeric(enriched.get("precursor_charge", pd.Series(dtype=float)), errors="coerce").dropna().astype(int).tolist()}),
            "charge_summary_rows": charge_summary.head(20).to_dict(orient="records") if not charge_summary.empty else [],
        },
        "confidence_summary": {
            "best_confidence": float(confidence_series.max()) if not confidence_series.empty else None,
            "mean_confidence": float(confidence_series.mean()) if not confidence_series.empty else None,
            "median_confidence": float(confidence_series.median()) if not confidence_series.empty else None,
        },
        "spectrum_evidence_summary": {
            "best_spectral_angle": float(spectral_series.max()) if not spectral_series.empty else None,
            "mean_spectral_angle": float(spectral_series.mean()) if not spectral_series.empty else None,
            "top_supported_peptides": peptide_summary.sort_values(
                by=["spectral_angle_max", "matched_ion_count_max", "total_psm_count"],
                ascending=[False, False, False],
                kind="mergesort",
            ).head(10).to_dict(orient="records") if not peptide_summary.empty else [],
        },
        "modification_summary": {
            "modified_psm_fraction": (modified_count / len(enriched)) if len(enriched) else None,
            "top_modified_peptides": peptide_summary.sort_values(by=["modified_psm_count", "total_psm_count"], ascending=[False, False], kind="mergesort").head(10).to_dict(orient="records") if not peptide_summary.empty else [],
        },
        "sdrf_summary": {
            "available": sdrf_summary is not None,
            "assay_count": int(len((sdrf_summary or {}).get("assays") or [])),
            "data_file_count": int(len((sdrf_summary or {}).get("data_files") or [])),
            "sample_rows": sample_summary.to_dict(orient="records") if sample_summary is not None else [],
        },
    }

    return DatasetAggregation(
        frame=enriched,
        dataset_summary=_json_safe(dataset_summary),
        peptide_summary=peptide_summary,
        peptidoform_summary=peptidoform_summary,
        run_summary=run_summary,
        sample_summary=sample_summary,
        charge_summary=charge_summary,
        sdrf_summary=sdrf_summary,
    )


def summarize_dataset_path(path: str | Path, *, sdrf_summary: dict[str, Any] | None = None) -> DatasetAggregation:
    return summarize_dataset_frame(load_scored_psm_table(path), sdrf_summary=sdrf_summary)


def summarize_candidate_inspection_frame(frame: pd.DataFrame, *, sdrf_summary: dict[str, Any] | None = None) -> CandidateInspection:
    enriched = _enrich_with_sdrf(frame, sdrf_summary)
    peptide_rows: list[dict[str, Any]] = []
    for sequence, group in enriched.groupby("stripped_sequence", sort=True):
        row = {"stripped_sequence": sequence}
        row.update(_summarize_group(group, include_peptide_context=True))
        peptide_rows.append(row)
    peptide_summary = _frame_from_records(peptide_rows, "best_confidence")
    run_summary = _summarize_runs(enriched)
    sample_summary = _summarize_samples(enriched)

    confidence_series = pd.to_numeric(enriched.get("confidence", pd.Series(dtype=float)), errors="coerce").dropna()
    spectral_series = pd.to_numeric(enriched.get("spectral_angle", pd.Series(dtype=float)), errors="coerce").dropna()
    matched_series = pd.to_numeric(enriched.get("matched_ion_count", pd.Series(dtype=float)), errors="coerce").dropna()
    explained_series = pd.to_numeric(enriched.get("explained_intensity_fraction", pd.Series(dtype=float)), errors="coerce").dropna()
    precursor_counts = _distinct_precursor_counts(enriched)
    dataset_summary = {
        "dataset_metadata": {
            "total_psm_count": int(len(enriched)),
            "distinct_peptide_count": int(enriched["stripped_sequence"].nunique()),
            "distinct_peptidoform_count": int(enriched["peptidoform"].nunique()),
            **precursor_counts,
            "distinct_run_count": int(enriched["run_key"].dropna().nunique()),
            "source_engines": sorted(set(enriched["source_engine"].tolist())),
        },
        "confidence_summary": {
            "best_confidence": float(confidence_series.max()) if not confidence_series.empty else None,
            "mean_confidence": float(confidence_series.mean()) if not confidence_series.empty else None,
            "median_confidence": float(confidence_series.median()) if not confidence_series.empty else None,
        },
        "sequence_summary": {
            "peptides_with_multiple_psms": int((peptide_summary.get("total_psm_count", pd.Series(dtype=float)) > 1).sum()) if not peptide_summary.empty else 0,
            "peptides_with_multiple_peptidoforms": int((peptide_summary.get("distinct_peptidoform_count", pd.Series(dtype=float)) > 1).sum()) if not peptide_summary.empty else 0,
            "modified_psm_count": int(enriched["is_modified"].sum()),
            "unmodified_psm_count": int((~enriched["is_modified"]).sum()),
        },
        "scoring_summary": {
            "scored_psm_count": int(spectral_series.shape[0]),
            "best_spectral_angle": float(spectral_series.max()) if not spectral_series.empty else None,
            "mean_spectral_angle": float(spectral_series.mean()) if not spectral_series.empty else None,
            "mean_matched_ion_count": float(matched_series.mean()) if not matched_series.empty else None,
            "mean_explained_intensity_fraction": float(explained_series.mean()) if not explained_series.empty else None,
        },
        "run_summary": {
            "runs": run_summary.to_dict(orient="records"),
        },
        "sdrf_summary": {
            "available": sdrf_summary is not None,
            "assay_count": int(len((sdrf_summary or {}).get("assays") or [])),
            "data_file_count": int(len((sdrf_summary or {}).get("data_files") or [])),
            "sample_rows": sample_summary.to_dict(orient="records") if sample_summary is not None else [],
        },
    }
    return CandidateInspection(
        frame=enriched,
        peptide_summary=peptide_summary,
        dataset_summary=_json_safe(dataset_summary),
        run_summary=run_summary,
        sample_summary=sample_summary,
        sdrf_summary=sdrf_summary,
    )


def inspection_psm_table(inspection: CandidateInspection) -> pd.DataFrame:
    columns = [column for column in ID_PSM_COLUMNS if column in inspection.frame.columns]
    return inspection.frame[columns].copy()


def inspection_peptide_table(inspection: CandidateInspection) -> pd.DataFrame:
    columns: list[str] = []
    for column in ID_PEPTIDE_COLUMNS:
        if column not in inspection.peptide_summary.columns:
            continue
        if column in OPTIONAL_PEPTIDE_COLUMNS:
            series = inspection.peptide_summary[column]
            if not series.map(
                lambda value: (
                    value is not None
                    and not (isinstance(value, float) and math.isnan(value))
                    and value != {}
                    and value != []
                )
            ).any():
                continue
        columns.append(column)
    return inspection.peptide_summary[columns].copy()


def _write_frame(frame: pd.DataFrame, path: Path) -> None:
    serializable = frame.copy()
    for column in serializable.columns:
        serializable[column] = serializable[column].map(_serialize_cell)
    serializable.to_csv(path, index=False)


def render_table_browse_html(frame: pd.DataFrame, title: str, output_path: str | Path | None = None, csv_filename: str | None = None) -> str:
    serializable = frame.copy()
    for column in serializable.columns:
        serializable[column] = serializable[column].map(lambda value: _serialize_browse_table_cell(str(column), value))
    table_html = _browse_table_html(serializable)
    csv_payload = json.dumps(serializable.to_csv(index=False))
    csv_filename_json = json.dumps(csv_filename or f"{title.lower().replace(' ', '_')}.csv")
    column_options = "".join(
        f'''<label class="column-option"><input type="checkbox" data-column-visible="{index}" checked /> <span>{html_lib.escape(str(column))}</span><input type="checkbox" data-column-pinned="{index}" aria-label="Pin {html_lib.escape(str(column))}" /><span class="pin-label">Pin</span></label>'''
        for index, column in enumerate(serializable.columns)
    )
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>{title}</title>
  <style>
    :root {{
      color-scheme: light;
      --bg: #f6f8fb;
      --fg: #172033;
      --muted: #556179;
      --panel: #ffffff;
      --border: #d6ddea;
      --thead-bg: #eef3f9;
      --thead-fg: #38445b;
      --row-alt: #fafcff;
      --row-hover: #f1f6fc;
      --row-border: #e8edf5;
    }}
    body {{
      margin: 0;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: var(--bg);
      color: var(--fg);
    }}
    .page {{
      padding: 20px 24px 28px;
    }}
    h1 {{
      margin: 0 0 6px;
      font-size: 24px;
      font-weight: 650;
    }}
    .meta {{
      margin: 0 0 16px;
      color: var(--muted);
      font-size: 13px;
    }}
    .toolbar {{
      display: flex;
      gap: 8px;
      align-items: center;
      flex-wrap: wrap;
      margin: 0 0 14px;
    }}
    .toolbar-label {{
      font-size: 12px;
      color: var(--muted);
      text-transform: uppercase;
      letter-spacing: 0.02em;
    }}
    .toolbar button {{
      border: 1px solid var(--border);
      background: var(--panel);
      color: var(--fg);
      border-radius: 6px;
      padding: 6px 10px;
      font-size: 12px;
      cursor: pointer;
    }}
    .toolbar button.active {{
      border-color: #068D9D;
      box-shadow: inset 0 0 0 1px #068D9D;
    }}
    .column-control {{
      position: relative;
    }}
    .column-menu {{
      position: absolute;
      z-index: 10;
      top: calc(100% + 6px);
      left: 0;
      width: min(420px, calc(100vw - 48px));
      max-height: min(60vh, 540px);
      overflow: auto;
      padding: 10px;
      border: 1px solid var(--border);
      border-radius: 6px;
      background: var(--panel);
      box-shadow: 0 8px 24px rgba(23, 32, 51, 0.18);
    }}
    .column-menu[hidden] {{ display: none; }}
    .column-menu-actions {{
      display: flex;
      gap: 6px;
      padding-bottom: 8px;
      margin-bottom: 6px;
      border-bottom: 1px solid var(--row-border);
    }}
    .column-option {{
      display: grid;
      grid-template-columns: auto minmax(0, 1fr) auto auto;
      gap: 7px;
      align-items: center;
      padding: 6px 2px;
      font-size: 12px;
    }}
    .column-option span {{ overflow-wrap: anywhere; }}
    .column-option .pin-label {{ color: var(--muted); }}
    .table-wrap {{
      border: 1px solid var(--border);
      border-radius: 6px;
      background: var(--panel);
      overflow: auto;
      max-height: calc(100vh - 120px);
      box-shadow: 0 1px 3px rgba(23, 32, 51, 0.08);
    }}
    table {{
      border-collapse: collapse;
      width: max-content;
      min-width: 100%;
      font-size: 12px;
      line-height: 1.35;
    }}
    thead th {{
      position: sticky;
      top: 0;
      z-index: 1;
      background: var(--thead-bg);
      color: var(--thead-fg);
      text-align: left;
      padding: 9px 10px;
      border-bottom: 1px solid var(--border);
      white-space: nowrap;
    }}
    th.pinned-column {{
      z-index: 3;
      box-shadow: 2px 0 0 var(--border);
    }}
    tbody td {{
      padding: 8px 10px;
      border-bottom: 1px solid var(--row-border);
      vertical-align: top;
      white-space: nowrap;
    }}
    td.pinned-column {{
      position: sticky;
      z-index: 2;
      background: var(--panel);
      box-shadow: 2px 0 0 var(--border);
    }}
    tbody tr:nth-child(even) {{
      background: var(--row-alt);
    }}
    tbody tr:nth-child(even) td.pinned-column {{ background: var(--row-alt); }}
    tbody tr:hover {{
      background: var(--row-hover);
    }}
    tbody tr:hover td.pinned-column {{ background: var(--row-hover); }}
    body[data-theme="dark"] {{
      --bg: #0f1726;
      --fg: #edf2fb;
      --muted: #aab7cc;
      --panel: #141d2d;
      --border: #2b3b56;
      --thead-bg: #1b273a;
      --thead-fg: #dbe5f3;
      --row-alt: #111b2b;
      --row-hover: #1a2740;
      --row-border: #223149;
    }}
    @media (prefers-color-scheme: dark) {{
      body[data-theme="auto"] {{
        --bg: #0f1726;
        --fg: #edf2fb;
        --muted: #aab7cc;
        --panel: #141d2d;
        --border: #2b3b56;
        --thead-bg: #1b273a;
        --thead-fg: #dbe5f3;
        --row-alt: #111b2b;
        --row-hover: #1a2740;
        --row-border: #223149;
      }}
    }}
  </style>
</head>
<body data-theme="auto">
  <div class="page">
    <h1>{title}</h1>
    <p class="meta">{len(serializable)} rows, {len(serializable.columns)} columns</p>
    <div class="toolbar" aria-label="Theme">
      <button type="button" onclick="downloadCsv()">Download CSV</button>
      <div class="column-control">
        <button type="button" id="column-menu-toggle" aria-expanded="false" aria-controls="column-menu">Columns</button>
        <div class="column-menu" id="column-menu" hidden>
          <div class="column-menu-actions">
            <button type="button" id="show-all-columns">Show all</button>
            <button type="button" id="clear-pinned-columns">Clear pins</button>
          </div>
          {column_options}
        </div>
      </div>
      <span class="toolbar-label">Theme</span>
      <button type="button" data-theme-button="light" class="active">Light</button>
      <button type="button" data-theme-button="dark">Dark</button>
      <button type="button" data-theme-button="auto">Auto</button>
    </div>
    <div class="table-wrap">
      {table_html}
    </div>
  </div>
  <script>
    (function() {{
      const csvContent = {csv_payload};
      const csvFilename = {csv_filename_json};
      const body = document.body;
      const buttons = Array.from(document.querySelectorAll('[data-theme-button]'));
      const table = document.querySelector('.browse-table');
      const columnMenu = document.getElementById('column-menu');
      const columnMenuToggle = document.getElementById('column-menu-toggle');
      const visibilityControls = Array.from(document.querySelectorAll('[data-column-visible]'));
      const pinControls = Array.from(document.querySelectorAll('[data-column-pinned]'));
      window.downloadCsv = function downloadCsv() {{
        const blob = new Blob([csvContent], {{ type: 'text/csv;charset=utf-8;' }});
        const url = URL.createObjectURL(blob);
        const link = document.createElement('a');
        link.href = url;
        link.download = csvFilename;
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
        URL.revokeObjectURL(url);
      }};
      function applyTheme(theme) {{
        body.setAttribute('data-theme', theme);
        buttons.forEach((button) => button.classList.toggle('active', button.getAttribute('data-theme-button') === theme));
      }}
      buttons.forEach((button) => {{
        button.addEventListener('click', () => applyTheme(button.getAttribute('data-theme-button')));
      }});
      function cellsForColumn(index) {{
        return Array.from(table.querySelectorAll(`[data-column-index="${{index}}"]`));
      }}
      function updateColumns() {{
        let left = 0;
        visibilityControls.forEach((control) => {{
          const index = control.getAttribute('data-column-visible');
          const pinControl = document.querySelector(`[data-column-pinned="${{index}}"]`);
          const visible = control.checked;
          const pinned = visible && pinControl.checked;
          pinControl.disabled = !visible;
          cellsForColumn(index).forEach((cell) => {{
            cell.hidden = !visible;
            cell.classList.toggle('pinned-column', pinned);
            cell.style.left = '';
          }});
        }});
        visibilityControls.forEach((control) => {{
          const index = control.getAttribute('data-column-visible');
          const pinControl = document.querySelector(`[data-column-pinned="${{index}}"]`);
          if (!control.checked || !pinControl.checked) return;
          const header = table.querySelector(`th[data-column-index="${{index}}"]`);
          cellsForColumn(index).forEach((cell) => {{ cell.style.left = `${{left}}px`; }});
          left += header.getBoundingClientRect().width;
        }});
      }}
      columnMenuToggle.addEventListener('click', () => {{
        const isOpen = !columnMenu.hidden;
        columnMenu.hidden = isOpen;
        columnMenuToggle.setAttribute('aria-expanded', String(!isOpen));
      }});
      document.getElementById('show-all-columns').addEventListener('click', () => {{
        visibilityControls.forEach((control) => {{ control.checked = true; }});
        updateColumns();
      }});
      document.getElementById('clear-pinned-columns').addEventListener('click', () => {{
        pinControls.forEach((control) => {{ control.checked = false; }});
        updateColumns();
      }});
      visibilityControls.concat(pinControls).forEach((control) => control.addEventListener('change', updateColumns));
      window.addEventListener('resize', updateColumns);
      applyTheme('auto');
      updateColumns();
    }})();
  </script>
</body>
</html>
"""
    if output_path:
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text(html, encoding="utf-8")
    return html


def _summary_value_html(value: Any) -> str:
    if isinstance(value, dict):
        items = "".join(f"<tr><th>{key}</th><td>{_summary_value_html(val)}</td></tr>" for key, val in value.items())
        return f"<table class=\"summary-table\">{items}</table>"
    if isinstance(value, list):
        if value and all(isinstance(item, dict) for item in value):
            frame = pd.DataFrame(value)
            serializable = frame.copy()
            for column in serializable.columns:
                serializable[column] = serializable[column].map(_serialize_cell)
            return serializable.to_html(index=False, border=0, classes="summary-grid", escape=False)
        items = "".join(f"<li>{_summary_value_html(item)}</li>" for item in value)
        return f"<ul>{items}</ul>"
    if value is None:
        return "<span class=\"muted\">NA</span>"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def render_summary_html(summary: dict[str, Any], title: str, output_path: str | Path | None = None) -> str:
    sections = "".join(
        f"<section><h2>{key}</h2>{_summary_value_html(value)}</section>"
        for key, value in summary.items()
    )
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>{title}</title>
  <style>
    body {{
      margin: 0;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: #f6f8fb;
      color: #172033;
    }}
    .page {{
      padding: 20px 24px 32px;
      max-width: 1400px;
    }}
    h1 {{
      margin: 0 0 20px;
      font-size: 24px;
      font-weight: 650;
    }}
    h2 {{
      margin: 0 0 12px;
      font-size: 16px;
      font-weight: 650;
    }}
    section {{
      margin: 0 0 20px;
      padding: 14px 16px;
      background: #ffffff;
      border: 1px solid #d6ddea;
      border-radius: 6px;
    }}
    .summary-table, .summary-grid {{
      border-collapse: collapse;
      width: 100%;
      font-size: 12px;
    }}
    .summary-table th, .summary-table td, .summary-grid th, .summary-grid td {{
      padding: 8px 10px;
      border-bottom: 1px solid #e8edf5;
      text-align: left;
      vertical-align: top;
      white-space: nowrap;
    }}
    .summary-table th {{
      width: 280px;
      color: #38445b;
    }}
    ul {{
      margin: 0;
      padding-left: 18px;
    }}
    .muted {{
      color: #77839a;
    }}
  </style>
</head>
<body>
  <div class="page">
    <h1>{title}</h1>
    {sections}
  </div>
</body>
</html>
"""
    if output_path:
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text(html, encoding="utf-8")
    return html


def render_dashboard_bundle_html(summary: dict[str, Any], title: str, output_path: str | Path | None = None) -> str:
    # Accept either the bundle shape produced by build_id_dashboard_tool
    # ({"dashboard": ..., "overlap_peptide": ..., "modifications": ..., "file_coverage": ...})
    # or the flat output of summarize-dashboard (which has "overview"/"histograms" at the top level).
    if "dashboard" in summary:
        dashboard = summary.get("dashboard") or {}
    elif any(key in summary for key in ("overview", "histograms", "ids_per_file", "scatter_2d", "grids_2d", "charge_distribution")):
        dashboard = summary
    else:
        dashboard = {}
    overview = dashboard.get("overview") or {}
    ids_per_file = ((dashboard.get("ids_per_file") or {}).get("counts")) or {}
    histograms = dashboard.get("histograms") or {}
    scatter_2d = dashboard.get("scatter_2d") or {}
    overlap = summary.get("overlap_peptide") or {}
    modifications = summary.get("modifications") or {}
    file_coverage = (summary.get("file_coverage") or {}).get("files") or {}

    if go is None:  # pragma: no cover
        raise RuntimeError("Dashboard rendering requires plotly.")

    def _fmt_count(value: Any) -> str:
        if value is None:
            return "NA"
        try:
            return f"{int(value):,}"
        except Exception:
            return str(value)

    def _fmt_value(value: float, decimals: int = 0) -> str:
        if decimals == 0:
            return f"{int(round(value)):,}" if math.isfinite(value) else "NA"
        return f"{value:,.{decimals}f}" if math.isfinite(value) else "NA"

    def _fmt_bin(left: float, right: float, *, integer: bool = False, decimals: int = 2) -> str:
        if integer:
            return str(int(round(left)))
        return f"{left:,.{decimals}f}-{right:,.{decimals}f}"

    run_names = list(file_coverage.keys())
    common_prefix = ""
    if run_names:
        common_prefix = run_names[0]
        for run_name in run_names[1:]:
            while common_prefix and not run_name.startswith(common_prefix):
                common_prefix = common_prefix[:-1]
    common_prefix = common_prefix.rstrip("_:- ")

    def _run_label(run_name: str) -> str:
        trimmed = run_name[len(common_prefix) :].lstrip("_:- ") if common_prefix and run_name.startswith(common_prefix) else run_name
        match = re.search(r"(rep\d+)$", trimmed)
        return match.group(1) if match else trimmed or run_name

    run_labels = {_name: _run_label(_name) for _name in run_names}

    plots: list[tuple[str, str]] = []

    def histogram_plot(hist: dict[str, Any], title_text: str, div_id: str, color: str, *, integer: bool = False) -> str | None:
        bins = hist.get("bins") or []
        if not bins:
            return None
        x = [_fmt_bin(row["left"], row["right"], integer=integer, decimals=2) for row in bins]
        y = [row["count"] for row in bins]
        fig = go.Figure(
            data=[go.Bar(x=x, y=y, marker_color=color)],
            layout=go.Layout(title=title_text, template="plotly_white", height=320, margin=dict(l=40, r=20, t=48, b=120)),
        )
        fig.update_traces(hovertemplate=f"{title_text}: %{{x}}<br>Count: %{{y:,}}<extra></extra>")
        fig.update_yaxes(tickformat=",d")
        fig.update_xaxes(tickangle=-35)
        return _plotly_div(fig, div_id)

    def categorical_plot(counts: dict[str, Any], title_text: str, div_id: str, color: str) -> str | None:
        if not counts:
            return None
        items = list(counts.items())
        fig = go.Figure(
            data=[go.Bar(x=[str(key) for key, _ in items], y=[value for _, value in items], marker_color=color)],
            layout=go.Layout(title=title_text, template="plotly_white", height=320, margin=dict(l=40, r=20, t=48, b=120)),
        )
        fig.update_traces(hovertemplate=f"{title_text}: %{{x}}<br>Count: %{{y:,}}<extra></extra>")
        fig.update_yaxes(tickformat=",d")
        fig.update_xaxes(tickangle=-30)
        return _plotly_div(fig, div_id)

    def scatter_plot(
        payload: dict[str, Any],
        title_text: str,
        div_id: str,
        *,
        x_label: str,
        y_label: str,
    ) -> str | None:
        points = payload.get("points") or []
        if not points:
            return None
        palette = ["#068D9D", "#FA824C", "#B8336A", "#9FD356", "#9792E3", "#F3A712"]
        groups = sorted({str(point.get("group")) for point in points if point.get("group") is not None})
        if not groups:
            groups = ["all"]
        traces: list[go.Scattergl] = []
        buttons: list[dict[str, Any]] = []
        trace_groups: list[str] = []
        for idx, group in enumerate(groups):
            subset = [point for point in points if (str(point.get("group")) if point.get("group") is not None else "all") == group]
            if not subset:
                continue
            display_group = run_labels.get(group, group) if group != "all" else "all"
            customdata = [
                [
                    point.get("candidate_id"),
                    point.get("peptidoform"),
                    run_labels.get(str(point.get("run_key")), str(point.get("run_key"))) if point.get("run_key") is not None else None,
                    point.get("scan_number"),
                ]
                for point in subset
            ]
            traces.append(
                go.Scattergl(
                    x=[point["x"] for point in subset],
                    y=[point["y"] for point in subset],
                    mode="markers",
                    name=display_group,
                    visible=len(traces) == 0,
                    marker=dict(color=palette[idx % len(palette)], size=7, opacity=0.72),
                    customdata=customdata,
                    hovertemplate=(
                        f"{x_label}: %{{x:.4g}}<br>"
                        f"{y_label}: %{{y:.4g}}<br>"
                        "Run: %{customdata[2]}<br>"
                        "Scan: %{customdata[3]}<br>"
                        "Peptidoform: %{customdata[1]}<br>"
                        "Candidate: %{customdata[0]}<extra></extra>"
                    ),
                )
            )
            trace_groups.append(group)
        if not traces:
            return None
        for idx, group in enumerate(trace_groups):
            visible = [False] * len(traces)
            visible[idx] = True
            buttons.append(
                dict(
                    label=run_labels.get(group, group) if group != "all" else "all",
                    method="update",
                    args=[
                        {"visible": visible},
                        {"title": f"{title_text} ({run_labels.get(group, group) if group != 'all' else 'all'})"},
                    ],
                )
            )
        initial_group = trace_groups[0]
        fig = go.Figure(data=traces)
        fig.update_layout(
            title=f"{title_text} ({run_labels.get(initial_group, initial_group) if initial_group != 'all' else 'all'})",
            template="plotly_white",
            height=340,
            margin=dict(l=48, r=20, t=70, b=64),
            showlegend=False,
            updatemenus=[dict(type="dropdown", x=1.0, y=1.18, xanchor="right", buttons=buttons)] if len(trace_groups) > 1 else [],
        )
        fig.update_xaxes(title=x_label)
        fig.update_yaxes(title=y_label)
        return _plotly_div(fig, div_id)

    def add_plot(div: str | None, klass: str = "plot") -> None:
        if div:
            plots.append((div, klass))

    plot_mz = histogram_plot(histograms.get("precursor_mz") or {}, "Precursor m/z", "hist_mz", "#068D9D")
    plot_confidence = histogram_plot(histograms.get("confidence") or {}, "Confidence", "hist_confidence", "#9792E3")
    plot_rt = histogram_plot(histograms.get("retention_time") or {}, "RT distribution", "hist_rt", "#9FD356")
    plot_charge = categorical_plot((dashboard.get("charge_distribution") or {}).get("counts") or {}, "Charge distribution", "hist_charge", "#068D9D")
    plot_peptide_length = histogram_plot(histograms.get("peptide_length") or {}, "Peptide length", "hist_peptide_length", "#FA824C", integer=True)
    plot_spectral_angle = histogram_plot(histograms.get("spectral_angle") or {}, "Spectral angle", "hist_spectral_angle", "#B8336A")
    plot_matched_ion_fraction = histogram_plot(histograms.get("matched_ion_fraction") or {}, "Matched ion fraction", "hist_matched_fraction", "#9FD356")
    plot_explained_intensity_fraction = histogram_plot(histograms.get("explained_intensity_fraction") or {}, "Explained intensity fraction", "hist_explained_fraction", "#F3A712")
    plot_log2_xcorr = histogram_plot(histograms.get("log2_xcorr") or {}, "log2 XCorr", "hist_log2_xcorr", "#068D9D")
    plot_confidence_vs_spectral_angle = scatter_plot(
        scatter_2d.get("confidence_vs_spectral_angle") or {},
        "Conf. vs SA",
        "scatter_confidence_vs_spectral_angle",
        x_label="Confidence",
        y_label="Spectral angle",
    )
    plot_spectral_angle_vs_pearson = scatter_plot(
        scatter_2d.get("spectral_angle_vs_pearson") or {},
        "SA vs Pearson",
        "scatter_spectral_angle_vs_pearson",
        x_label="Spectral angle",
        y_label="Pearson correlation",
    )
    plot_spectral_angle_vs_log2_xcorr = scatter_plot(
        scatter_2d.get("spectral_angle_vs_log2_xcorr") or {},
        "SA vs XCorr",
        "scatter_spectral_angle_vs_log2_xcorr",
        x_label="Spectral angle",
        y_label="log2 XCorr",
    )
    plot_confidence_vs_matched_ion_count = scatter_plot(
        scatter_2d.get("confidence_vs_matched_ion_count") or {},
        "Conf. vs ions detected",
        "scatter_confidence_vs_matched_ion_count",
        x_label="Confidence",
        y_label="Matched ion count",
    )
    plot_spectral_angle_vs_precursor_mass_error = scatter_plot(
        scatter_2d.get("spectral_angle_vs_precursor_mass_error_ppm") or {},
        "SA vs precursor mass error",
        "scatter_spectral_angle_vs_precursor_mass_error",
        x_label="Spectral angle",
        y_label="Precursor mass error (ppm)",
    )
    plot_spectral_angle_vs_fragment_mass_error = scatter_plot(
        scatter_2d.get("spectral_angle_vs_fragment_mass_error_abs_mean_ppm") or {},
        "SA vs fragment avg mass error",
        "scatter_spectral_angle_vs_fragment_mass_error",
        x_label="Spectral angle",
        y_label="Fragment abs. avg mass error (ppm)",
    )
    plot_file_coverage: str | None = None

    if file_coverage:
        coverage_rows = []
        for run, run_data in file_coverage.items():
            coverage_rows.append(
                {
                    "run": run_labels.get(run, run),
                    "psms": run_data.get("overview", {}).get("total_psm_count"),
                    "peptides": run_data.get("overview", {}).get("distinct_peptide_count"),
                    "peptidoforms": run_data.get("overview", {}).get("distinct_peptidoform_count"),
                    "modified_psms": run_data.get("overview", {}).get("modified_psm_count"),
                }
            )
        coverage_frame = pd.DataFrame(coverage_rows)
        if not coverage_frame.empty:
            fig = go.Figure()
            fig.add_bar(name="PSMs", x=coverage_frame["run"], y=coverage_frame["psms"], marker_color="#068D9D")
            fig.add_bar(name="Peptides", x=coverage_frame["run"], y=coverage_frame["peptides"], marker_color="#9FD356")
            fig.add_bar(name="Peptidoforms", x=coverage_frame["run"], y=coverage_frame["peptidoforms"], marker_color="#F3A712")
            fig.update_layout(title="Per-file coverage", template="plotly_white", barmode="group", height=360, margin=dict(l=40, r=20, t=48, b=120))
            fig.update_traces(hovertemplate="%{fullData.name}<br>Run: %{x}<br>Count: %{y:,}<extra></extra>")
            fig.update_yaxes(tickformat=",d")
            fig.update_xaxes(tickangle=-35)
            plot_file_coverage = _plotly_div(fig, "file_coverage")

    pairwise = overlap.get("pairwise") or []
    if pairwise:
        labels = sorted({run_labels.get(row["group_a"], row["group_a"]) for row in pairwise} | {run_labels.get(row["group_b"], row["group_b"]) for row in pairwise})
        index = {label: idx for idx, label in enumerate(labels)}
        z = [[1.0 if i == j else 0.0 for j in range(len(labels))] for i in range(len(labels))]
        for row in pairwise:
            left = run_labels.get(row["group_a"], row["group_a"])
            right = run_labels.get(row["group_b"], row["group_b"])
            i = index[left]
            j = index[right]
            z[i][j] = row["jaccard"] or 0.0
            z[j][i] = row["jaccard"] or 0.0
        fig = go.Figure(
            data=[go.Heatmap(z=z, x=labels, y=labels, colorscale="Viridis", zmin=0.0, zmax=1.0)],
            layout=go.Layout(title="Peptide overlap (Jaccard)", template="plotly_white", height=360, margin=dict(l=80, r=20, t=48, b=80)),
        )
        fig.update_traces(hovertemplate="Run A: %{y}<br>Run B: %{x}<br>Jaccard: %{z:.3f}<extra></extra>")
        plot_overlap = _plotly_div(fig, "overlap_heatmap")
    else:
        plot_overlap = None

    modification_counts = modifications.get("modification_counts") or {}
    if modification_counts:
        items = list(modification_counts.items())[:12]
        fig = go.Figure(
            data=[go.Bar(x=[k for k, _ in items], y=[v for _, v in items], marker_color="#FA824C")],
            layout=go.Layout(title="Modification counts", template="plotly_white", height=320, margin=dict(l=40, r=20, t=48, b=120)),
        )
        fig.update_traces(hovertemplate="Modification: %{x}<br>Count: %{y:,}<extra></extra>")
        fig.update_yaxes(tickformat=",d")
        fig.update_xaxes(tickangle=-30)
        plot_modification_counts = _plotly_div(fig, "modification_counts")
    else:
        plot_modification_counts = None

    plot_confidence_by_run: str | None = None
    plot_charge_by_run: str | None = None
    plot_peptide_length_by_run: str | None = None
    plot_rt_by_run: str | None = None
    plot_modifications_by_run: str | None = None
    if file_coverage:
        def add_run_dropdown_plot(metric_title: str, div_id: str, extractor: Any, color: str, *, integer: bool = False) -> str | None:
            traces = []
            buttons = []
            trace_runs: list[str] = []
            for idx, run in enumerate(run_names):
                data = extractor(file_coverage[run])
                if not data:
                    continue
                x, y = data
                if not x or not y:
                    continue
                visible = len(traces) == 0
                traces.append(go.Bar(x=x, y=y, marker_color=color, visible=visible, name=run_labels.get(run, run)))
                trace_runs.append(run)
            if not traces:
                return
            for idx, run in enumerate(trace_runs):
                visible = [False] * len(traces)
                visible[idx] = True
                buttons.append(
                    dict(
                        label=run_labels.get(run, run),
                        method="update",
                        args=[{"visible": visible}, {"title": f"{metric_title} ({run_labels.get(run, run)})"}],
                    )
                )
            fig = go.Figure(data=traces)
            fig.update_layout(
                title=f"{metric_title} ({run_labels.get(trace_runs[0], trace_runs[0])})",
                template="plotly_white",
                height=340,
                margin=dict(l=40, r=20, t=70, b=120),
                updatemenus=[dict(type="dropdown", x=1.0, y=1.18, xanchor="right", buttons=buttons)],
            )
            fig.update_traces(hovertemplate="%{x}<br>Count: %{y:,}<extra></extra>")
            fig.update_yaxes(tickformat=",d")
            fig.update_xaxes(tickangle=-35)
            return _plotly_div(fig, div_id)

        plot_confidence_by_run = add_run_dropdown_plot(
            "Confidence distribution by run",
            "confidence_by_run",
            lambda run_data: (
                [f"{row['left']:.3f}-{row['right']:.3f}" for row in ((run_data.get("confidence_histogram") or {}).get("bins") or [])],
                [row["count"] for row in ((run_data.get("confidence_histogram") or {}).get("bins") or [])],
            ),
            "#9792E3",
        )
        plot_charge_by_run = add_run_dropdown_plot(
            "Charge distribution by run",
            "charge_by_run",
            lambda run_data: (
                list((run_data.get("charge_counts") or {}).keys()),
                list((run_data.get("charge_counts") or {}).values()),
            ),
            "#068D9D",
        )
        plot_peptide_length_by_run = add_run_dropdown_plot(
            "Peptide length by run",
            "peptide_length_by_run",
            lambda run_data: (
                [_fmt_bin(row["left"], row["right"], integer=True) for row in ((run_data.get("peptide_length_histogram") or {}).get("bins") or [])],
                [row["count"] for row in ((run_data.get("peptide_length_histogram") or {}).get("bins") or [])],
            ),
            "#FA824C",
            integer=True,
        )
        plot_rt_by_run = add_run_dropdown_plot(
            "RT distribution by run",
            "rt_by_run",
            lambda run_data: (
                [f"{row['left']:.2f}-{row['right']:.2f}" for row in ((run_data.get("rt_histogram") or {}).get("bins") or [])],
                [row["count"] for row in ((run_data.get("rt_histogram") or {}).get("bins") or [])],
            ),
            "#9FD356",
        )

        group_mod_counts = modifications.get("group_modification_counts") or {}
        if group_mod_counts:
            traces = []
            buttons = []
            available_runs = [run for run in run_names if run in group_mod_counts]
            for idx, run in enumerate(available_runs):
                items = list(group_mod_counts.get(run, {}).items())[:12]
                traces.append(
                    go.Bar(
                        x=[k for k, _ in items],
                        y=[v for _, v in items],
                        marker_color="#F3A712",
                        visible=idx == 0,
                        name=run_labels.get(run, run),
                    )
                )
            if traces:
                for idx, run in enumerate(available_runs):
                    visible = [False] * len(traces)
                    visible[idx] = True
                    buttons.append(dict(label=run_labels.get(run, run), method="update", args=[{"visible": visible}, {"title": f"Modification counts ({run_labels.get(run, run)})"}]))
                fig = go.Figure(data=traces)
                fig.update_layout(
                    title=f"Modification counts ({run_labels.get(available_runs[0], available_runs[0])})",
                    template="plotly_white",
                    height=340,
                    margin=dict(l=40, r=20, t=70, b=120),
                    updatemenus=[dict(type="dropdown", x=1.0, y=1.18, xanchor="right", buttons=buttons)],
                )
                fig.update_traces(hovertemplate="Modification: %{x}<br>Count: %{y:,}<extra></extra>")
                fig.update_yaxes(tickformat=",d")
                fig.update_xaxes(tickangle=-30)
                plot_modifications_by_run = _plotly_div(fig, "modifications_by_run")

    add_plot(plot_overlap)
    add_plot(plot_file_coverage, "plot plot-wide")
    add_plot(plot_mz)
    add_plot(plot_confidence, "plot plot-wide")
    add_plot(plot_confidence_by_run)
    add_plot(plot_spectral_angle, "plot plot-wide")
    add_plot(plot_matched_ion_fraction)
    add_plot(plot_explained_intensity_fraction)
    add_plot(plot_log2_xcorr)
    add_plot(plot_confidence_vs_spectral_angle, "plot plot-wide")
    add_plot(plot_spectral_angle_vs_pearson)
    add_plot(plot_spectral_angle_vs_log2_xcorr)
    add_plot(plot_confidence_vs_matched_ion_count)
    add_plot(plot_spectral_angle_vs_precursor_mass_error)
    add_plot(plot_spectral_angle_vs_fragment_mass_error)
    add_plot(plot_rt)
    add_plot(plot_rt_by_run)
    add_plot(plot_charge)
    add_plot(plot_charge_by_run)
    add_plot(plot_peptide_length)
    add_plot(plot_peptide_length_by_run)
    add_plot(plot_modification_counts)
    add_plot(plot_modifications_by_run)

    tiles = [
        ("Precursors", overview.get("distinct_peptidoform_precursor_count")),
        ("PSMs", overview.get("total_psm_count")),
        ("Peptides", overview.get("distinct_peptide_count")),
        ("Peptidoforms", overview.get("distinct_peptidoform_count")),
        ("Modified Peptidoforms", modifications.get("distinct_modified_peptidoform_count")),
        ("Modified PSM Fraction", f"{(modifications.get('modified_fraction') or 0.0):.1%}" if modifications else None),
        ("Runs", overview.get("distinct_run_count")),
    ]
    tile_html = "".join(
        f"<div class='tile'><div class='label'>{label}</div><div class='value'>{_fmt_count(value) if isinstance(value, int) or (isinstance(value, float) and value.is_integer()) else (value if value is not None else 'NA')}</div></div>"
        for label, value in tiles
    )

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>{title}</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    :root {{
      color-scheme: light;
      --bg: #f6f8fb;
      --fg: #172033;
      --muted: #556179;
      --panel: #ffffff;
      --border: #d6ddea;
      --tile-border: #068D9D;
      --button-bg: #ffffff;
      --button-active: #068D9D;
      --button-active-fg: #ffffff;
    }}
    body {{
      margin: 0;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      background: var(--bg);
      color: var(--fg);
    }}
    body[data-theme="dark"] {{
      color-scheme: dark;
      --bg: #0b1020;
      --fg: #e5edf6;
      --muted: #94a3b8;
      --panel: #111827;
      --border: #243044;
      --tile-border: #9FD356;
      --button-bg: #111827;
      --button-active: #9FD356;
      --button-active-fg: #0b1020;
    }}
    @media (prefers-color-scheme: dark) {{
      body[data-theme="auto"] {{
        color-scheme: dark;
        --bg: #0b1020;
        --fg: #e5edf6;
        --muted: #94a3b8;
        --panel: #111827;
        --border: #243044;
        --tile-border: #9FD356;
        --button-bg: #111827;
        --button-active: #9FD356;
        --button-active-fg: #0b1020;
      }}
    }}
    .page {{
      padding: 20px 24px 32px;
      max-width: 1500px;
    }}
    h1 {{
      margin: 0 0 18px;
      font-size: 24px;
      font-weight: 650;
    }}
    h2 {{
      margin: 28px 0 12px;
      font-size: 16px;
      font-weight: 650;
    }}
    .toolbar {{
      display: flex;
      gap: 8px;
      align-items: center;
      margin: 0 0 16px;
    }}
    .toolbar-label {{
      font-size: 12px;
      color: var(--muted);
      text-transform: uppercase;
      letter-spacing: 0.02em;
    }}
    .toolbar button {{
      border: 1px solid var(--border);
      background: var(--button-bg);
      color: var(--fg);
      border-radius: 6px;
      padding: 6px 10px;
      font-size: 12px;
      cursor: pointer;
    }}
    .toolbar button.active {{
      border-color: var(--button-active);
      background: var(--button-active);
      color: var(--button-active-fg);
      box-shadow: inset 0 0 0 1px var(--button-active);
    }}
    .tiles {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
      gap: 12px;
      margin-bottom: 16px;
    }}
    .tile {{
      background: var(--panel);
      border: 1px solid var(--border);
      border-left: 4px solid var(--tile-border);
      border-radius: 6px;
      padding: 12px 14px;
    }}
    .label {{
      font-size: 12px;
      color: var(--muted);
      margin-bottom: 6px;
      text-transform: uppercase;
    }}
    .value {{
      font-size: 22px;
      font-weight: 650;
    }}
    .plot-grid {{
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 16px;
    }}
    .plot {{
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 6px;
      padding: 4px;
    }}
    .plot-wide {{
      grid-column: span 2;
    }}
    @media (max-width: 1100px) {{
      .plot-grid {{
        grid-template-columns: 1fr;
      }}
      .plot-wide {{
        grid-column: span 1;
      }}
    }}
  </style>
</head>
<body>
  <div class="page">
    <h1>{title}</h1>
    <div class="toolbar">
      <span class="toolbar-label">Theme</span>
      <button type="button" data-theme-button="light">Light</button>
      <button type="button" data-theme-button="dark">Dark</button>
      <button type="button" data-theme-button="auto">Auto</button>
    </div>
    <div class="tiles">{tile_html}</div>
    <div class="plot-grid">
      {''.join(f"<div class='{klass}'>{html}</div>" for html, klass in plots)}
    </div>
  </div>
  <script>
    (function() {{
      const body = document.body;
      const buttons = Array.from(document.querySelectorAll('[data-theme-button]'));
      const media = window.matchMedia ? window.matchMedia('(prefers-color-scheme: dark)') : null;
      const themes = {{
        light: {{
          paper: '#ffffff',
          plot: '#ffffff',
          text: '#172033',
          muted: '#556179',
          grid: '#e8edf5',
          zero: '#b7c2d0',
          hoverBg: 'rgba(255,255,255,0.96)',
          hoverBorder: '#d6ddea',
          menuBg: '#ffffff',
          menuBorder: '#d6ddea',
        }},
        dark: {{
          paper: '#111827',
          plot: '#111827',
          text: '#e5edf6',
          muted: '#94a3b8',
          grid: '#223149',
          zero: '#334155',
          hoverBg: 'rgba(17,24,39,0.96)',
          hoverBorder: '#475569',
          menuBg: '#111827',
          menuBorder: '#243044',
        }},
      }};
      function storedTheme() {{
        return localStorage.getItem('proteomics-mcp-theme') || 'auto';
      }}
      function resolvedTheme(theme) {{
        if (theme === 'auto') {{
          return media && media.matches ? 'dark' : 'light';
        }}
        return theme;
      }}
      function relayoutPlot(gd, palette) {{
        const update = {{
          'paper_bgcolor': palette.paper,
          'plot_bgcolor': palette.plot,
          'font.color': palette.text,
          'title.font.color': palette.text,
          'legend.font.color': palette.text,
          'hoverlabel.bgcolor': palette.hoverBg,
          'hoverlabel.bordercolor': palette.hoverBorder,
          'hoverlabel.font.color': palette.text,
        }};
        Object.keys(gd.layout || {{}}).forEach((key) => {{
          if (key.startsWith('xaxis') || key.startsWith('yaxis')) {{
            update[`${{key}}.color`] = palette.text;
            update[`${{key}}.gridcolor`] = palette.grid;
            update[`${{key}}.linecolor`] = palette.grid;
            update[`${{key}}.zerolinecolor`] = palette.zero;
            update[`${{key}}.tickfont.color`] = palette.text;
            update[`${{key}}.title.font.color`] = palette.text;
          }}
        }});
        (gd.layout.updatemenus || []).forEach((_, index) => {{
          update[`updatemenus[${{index}}].bgcolor`] = palette.menuBg;
          update[`updatemenus[${{index}}].bordercolor`] = palette.menuBorder;
          update[`updatemenus[${{index}}].font.color`] = palette.text;
        }});
        Plotly.relayout(gd, update);
      }}
      function applyTheme(theme) {{
        body.setAttribute('data-theme', theme);
        buttons.forEach((button) => button.classList.toggle('active', button.getAttribute('data-theme-button') === theme));
        localStorage.setItem('proteomics-mcp-theme', theme);
        const palette = themes[resolvedTheme(theme)];
        document.querySelectorAll('.js-plotly-plot').forEach((gd) => relayoutPlot(gd, palette));
      }}
      buttons.forEach((button) => {{
        button.addEventListener('click', () => applyTheme(button.getAttribute('data-theme-button')));
      }});
      if (media) {{
        media.addEventListener('change', () => {{
          if (storedTheme() === 'auto') {{
            applyTheme('auto');
          }}
        }});
      }}
      window.addEventListener('load', () => applyTheme(storedTheme()));
    }})();
  </script>
</body>
</html>
"""
    if output_path:
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text(html, encoding="utf-8")
    return html


def render_prediction_accuracy_html(summary: dict[str, Any], title: str, output_path: str | Path | None = None) -> str:
    if go is None:  # pragma: no cover
        raise RuntimeError("Prediction accuracy rendering requires plotly.")

    overview = summary.get("overview") or {}
    reliability_bins = (summary.get("reliability") or {}).get("bins") or []
    operating_points = summary.get("operating_points") or []
    by_charge = summary.get("by_charge") or []
    by_length = summary.get("by_predicted_length") or []
    by_mod = summary.get("by_target_modification_state") or []
    confidence_dist = summary.get("confidence_distribution_by_correctness") or {}
    mass_error_dist = summary.get("mass_error_distribution_by_correctness") or {}

    def _fmt_pct(value: Any) -> str:
        if value is None or (isinstance(value, float) and not math.isfinite(value)):
            return "NA"
        return f"{float(value):.1%}"

    def _fmt_count(value: Any) -> str:
        if value is None:
            return "NA"
        try:
            return f"{int(value):,}"
        except Exception:
            return str(value)

    tiles = [
        ("PSMs scored", _fmt_count(overview.get("psm_count"))),
        ("Exact accuracy", _fmt_pct(overview.get("exact_accuracy"))),
        ("Sequence accuracy", _fmt_pct(overview.get("stripped_sequence_accuracy"))),
        ("I/L-equivalent accuracy", _fmt_pct(overview.get("il_equivalent_accuracy"))),
        ("Skipped (missing)", _fmt_count(overview.get("skipped_missing_pred_or_target"))),
    ]
    tile_html = "".join(
        f"<div class='tile'><div class='label'>{label}</div><div class='value'>{value}</div></div>"
        for label, value in tiles
    )

    plots: list[tuple[str, str]] = []

    def add_plot(html: str | None, klass: str = "plot") -> None:
        if html:
            plots.append((html, klass))

    if reliability_bins:
        labels = [f"[{row['left']:.1f},{row['right']:.1f}]" for row in reliability_bins]
        fig = go.Figure(
            data=[
                go.Bar(x=labels, y=[row["count"] for row in reliability_bins], name="PSM count", marker_color="#cbd5e1", opacity=0.6, yaxis="y2", hovertemplate="bin %{x}<br>n=%{y:,}<extra></extra>"),
                go.Scatter(x=labels, y=[row["accuracy_il"] for row in reliability_bins], mode="lines+markers", name="Accuracy (I/L-eq)", line=dict(color="#068D9D", width=2)),
                go.Scatter(x=labels, y=[row["mean_confidence"] for row in reliability_bins], mode="lines+markers", name="Mean confidence", line=dict(color="#9FD356", width=1, dash="dash")),
            ],
            layout=go.Layout(
                title="Reliability diagram (confidence vs empirical accuracy)",
                template="plotly_white",
                height=360,
                margin=dict(l=50, r=50, t=50, b=110),
                yaxis=dict(title="Accuracy / mean confidence", range=[0, 1.02]),
                yaxis2=dict(title="PSM count", overlaying="y", side="right", showgrid=False),
                legend=dict(orientation="h", y=1.12, x=0),
            ),
        )
        fig.update_xaxes(tickangle=-30)
        add_plot(_plotly_div(fig, "accuracy_reliability"), "plot plot-wide")

    if operating_points:
        thresholds_table = []
        thresholds_table.append(["conf ≥", "PSMs kept", "coverage", "acc (I/L-eq)", "acc (sequence)", "acc (exact)"])
        for row in operating_points:
            thresholds_table.append(
                [
                    f"{row['threshold']:.2f}",
                    _fmt_count(row.get("kept")),
                    _fmt_pct(row.get("coverage")),
                    _fmt_pct(row.get("accuracy_il")),
                    _fmt_pct(row.get("accuracy_seq")),
                    _fmt_pct(row.get("accuracy_exact")),
                ]
            )
        header_row = thresholds_table[0]
        body_rows = list(zip(*thresholds_table[1:]))
        fig = go.Figure(
            data=[
                go.Table(
                    header=dict(values=header_row, fill_color="#068D9D", font=dict(color="white", size=12), align="center"),
                    cells=dict(values=body_rows, align="center", height=24),
                )
            ],
            layout=go.Layout(title="Operating points (confidence thresholds)", height=320, margin=dict(l=20, r=20, t=50, b=20)),
        )
        add_plot(_plotly_div(fig, "accuracy_operating_points"), "plot plot-wide")

    if by_charge:
        x = [str(row["charge"]) for row in by_charge]
        fig = go.Figure(
            data=[
                go.Bar(x=x, y=[row["accuracy_exact"] for row in by_charge], name="Exact", marker_color="#0F4C81"),
                go.Bar(x=x, y=[row["accuracy_seq"] for row in by_charge], name="Sequence", marker_color="#068D9D"),
                go.Bar(x=x, y=[row["accuracy_il"] for row in by_charge], name="I/L-equiv", marker_color="#9FD356"),
            ],
            layout=go.Layout(
                title="Accuracy by precursor charge",
                template="plotly_white",
                height=320,
                barmode="group",
                margin=dict(l=50, r=20, t=50, b=60),
                yaxis=dict(title="Accuracy", range=[0, 1.0]),
            ),
        )
        add_plot(_plotly_div(fig, "accuracy_by_charge"))

    if by_length:
        fig = go.Figure(
            data=[
                go.Bar(
                    x=[row["length"] for row in by_length],
                    y=[row["accuracy_il"] for row in by_length],
                    customdata=[row["count"] for row in by_length],
                    marker_color="#068D9D",
                    hovertemplate="length %{x}<br>n=%{customdata:,}<br>acc=%{y:.1%}<extra></extra>",
                )
            ],
            layout=go.Layout(
                title="I/L-equivalent accuracy by predicted peptide length",
                template="plotly_white",
                height=320,
                margin=dict(l=50, r=20, t=50, b=60),
                yaxis=dict(title="Accuracy", range=[0, 1.02]),
                xaxis=dict(title="Predicted length (residues)"),
            ),
        )
        add_plot(_plotly_div(fig, "accuracy_by_length"))

    if by_mod:
        x = [row["label"] for row in by_mod]
        fig = go.Figure(
            data=[
                go.Bar(x=x, y=[row["accuracy_exact"] for row in by_mod], name="Exact", marker_color="#0F4C81"),
                go.Bar(x=x, y=[row["accuracy_seq"] for row in by_mod], name="Sequence", marker_color="#068D9D"),
                go.Bar(x=x, y=[row["accuracy_il"] for row in by_mod], name="I/L-equiv", marker_color="#9FD356"),
            ],
            layout=go.Layout(
                title="Accuracy by target modification state",
                template="plotly_white",
                height=320,
                barmode="group",
                margin=dict(l=50, r=20, t=50, b=60),
                yaxis=dict(title="Accuracy", range=[0, 1.0]),
            ),
        )
        add_plot(_plotly_div(fig, "accuracy_by_mod"))

    if confidence_dist.get("correct_il") and confidence_dist.get("wrong_il"):
        def _bin_centers(bins: list[dict[str, Any]]) -> list[float]:
            return [(row["left"] + row["right"]) / 2 for row in bins]

        correct_bins = confidence_dist["correct_il"]["bins"]
        wrong_bins = confidence_dist["wrong_il"]["bins"]
        fig = go.Figure(
            data=[
                go.Bar(x=_bin_centers(correct_bins), y=[row["count"] for row in correct_bins], name="Correct (I/L-eq)", marker_color="#9FD356", opacity=0.75),
                go.Bar(x=_bin_centers(wrong_bins), y=[row["count"] for row in wrong_bins], name="Wrong", marker_color="#E26D5C", opacity=0.6),
            ],
            layout=go.Layout(
                title="Confidence distribution by correctness",
                template="plotly_white",
                height=320,
                barmode="overlay",
                margin=dict(l=50, r=20, t=50, b=60),
                xaxis=dict(title="Confidence", range=[0, 1]),
                yaxis=dict(title="PSM count"),
                legend=dict(orientation="h", y=1.12, x=0),
            ),
        )
        add_plot(_plotly_div(fig, "accuracy_confidence_split"))

    if mass_error_dist.get("correct_il") and mass_error_dist.get("wrong_il"):
        correct_bins = mass_error_dist["correct_il"]["bins"]
        wrong_bins = mass_error_dist["wrong_il"]["bins"]
        centers_c = [(row["left"] + row["right"]) / 2 for row in correct_bins]
        centers_w = [(row["left"] + row["right"]) / 2 for row in wrong_bins]
        fig = go.Figure(
            data=[
                go.Bar(x=centers_c, y=[row["count"] for row in correct_bins], name="Correct (I/L-eq)", marker_color="#9FD356", opacity=0.75),
                go.Bar(x=centers_w, y=[row["count"] for row in wrong_bins], name="Wrong", marker_color="#E26D5C", opacity=0.6),
            ],
            layout=go.Layout(
                title="Mass error (Da, clipped ±5) by correctness",
                template="plotly_white",
                height=320,
                barmode="overlay",
                margin=dict(l=50, r=20, t=50, b=60),
                xaxis=dict(title="Mass error (Da)"),
                yaxis=dict(title="PSM count"),
                legend=dict(orientation="h", y=1.12, x=0),
            ),
        )
        add_plot(_plotly_div(fig, "accuracy_mass_error_split"))

    plot_html = "".join(f"<div class='{klass}'>{html}</div>" for html, klass in plots)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>{title}</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    body {{ margin: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; background: #f6f8fb; color: #172033; }}
    .page {{ padding: 20px 24px 32px; max-width: 1500px; }}
    h1 {{ margin: 0 0 18px; font-size: 24px; font-weight: 650; }}
    .tiles {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 12px; margin-bottom: 16px; }}
    .tile {{ background: #ffffff; border: 1px solid #d6ddea; border-left: 4px solid #068D9D; border-radius: 6px; padding: 12px 14px; }}
    .label {{ font-size: 12px; color: #556179; margin-bottom: 6px; text-transform: uppercase; }}
    .value {{ font-size: 22px; font-weight: 650; }}
    .plot-grid {{ display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 16px; }}
    .plot {{ background: #ffffff; border: 1px solid #d6ddea; border-radius: 6px; padding: 4px; }}
    .plot-wide {{ grid-column: span 2; }}
    @media (max-width: 1100px) {{ .plot-grid {{ grid-template-columns: 1fr; }} .plot-wide {{ grid-column: span 1; }} }}
  </style>
</head>
<body>
  <div class="page">
    <h1>{title}</h1>
    <div class="tiles">{tile_html}</div>
    <div class="plot-grid">{plot_html}</div>
  </div>
</body>
</html>
"""
    if output_path:
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text(html, encoding="utf-8")
    return html


def render_psm_table_html(inspection: CandidateInspection, output_path: str | Path) -> dict[str, str]:
    render_table_browse_html(
        inspection_psm_table(inspection),
        title="PSM table",
        output_path=output_path,
        csv_filename="psm_table.csv",
    )
    return {"format": "html", "output_path": str(output_path)}


def render_peptide_table_html(inspection: CandidateInspection, output_path: str | Path) -> dict[str, str]:
    render_table_browse_html(
        inspection_peptide_table(inspection),
        title="Peptide table",
        output_path=output_path,
        csv_filename="peptide_table.csv",
    )
    return {"format": "html", "output_path": str(output_path)}


def write_dataset_outputs(aggregation: DatasetAggregation, output_prefix: str | Path) -> dict[str, str]:
    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    dataset_summary_path = prefix.with_name(f"{prefix.name}_dataset_summary.json")
    peptide_path = prefix.with_name(f"{prefix.name}_peptide_summary.csv")
    peptidoform_path = prefix.with_name(f"{prefix.name}_peptidoform_summary.csv")
    run_path = prefix.with_name(f"{prefix.name}_run_summary.csv")
    charge_path = prefix.with_name(f"{prefix.name}_charge_summary.csv")

    dataset_summary_path.write_text(json.dumps(aggregation.dataset_summary, indent=2), encoding="utf-8")
    _write_frame(aggregation.peptide_summary, peptide_path)
    _write_frame(aggregation.peptidoform_summary, peptidoform_path)
    _write_frame(aggregation.run_summary, run_path)
    _write_frame(aggregation.charge_summary, charge_path)

    output_files = {
        "dataset_summary_json": str(dataset_summary_path),
        "peptide_summary_csv": str(peptide_path),
        "peptidoform_summary_csv": str(peptidoform_path),
        "run_summary_csv": str(run_path),
        "charge_summary_csv": str(charge_path),
    }
    if aggregation.sample_summary is not None:
        sample_path = prefix.with_name(f"{prefix.name}_sample_summary.csv")
        _write_frame(aggregation.sample_summary, sample_path)
        output_files["sample_summary_csv"] = str(sample_path)
    return output_files


def _plotly_div(fig: Any, div_id: str) -> str:
    fig_json = fig.to_json()
    return f"""
<div id="{div_id}" class="plot"></div>
<script>
  (function() {{
    const spec = {fig_json};
    Plotly.newPlot("{div_id}", spec.data, spec.layout, {{displayModeBar: false, responsive: true}});
  }})();
</script>
"""


def render_dataset_report_html(aggregation: DatasetAggregation, output_path: str | None = None) -> str:
    dataset = aggregation.dataset_summary["dataset_metadata"]
    confidence = aggregation.dataset_summary["confidence_summary"]
    spectrum = aggregation.dataset_summary["spectrum_evidence_summary"]
    run_summary = aggregation.run_summary.head(10)
    peptide_confident = aggregation.peptide_summary.sort_values(by=["best_confidence", "total_psm_count"], ascending=[False, False], kind="mergesort").head(10)
    peptide_supported = aggregation.peptide_summary.sort_values(by=["spectral_angle_max", "matched_ion_count_max", "total_psm_count"], ascending=[False, False, False], kind="mergesort").head(10)

    if go is None:  # pragma: no cover
        raise RuntimeError("Dataset report rendering requires plotly.")

    plots: list[str] = []
    if not run_summary.empty:
        run_fig = go.Figure(
            data=[go.Bar(x=run_summary["experiment_name"], y=run_summary["total_psm_count"], marker_color="#068D9D")],
            layout=go.Layout(title="PSMs per run", template="plotly_white", margin=dict(l=40, r=20, t=48, b=120), height=360),
        )
        run_fig.update_xaxes(tickangle=-35)
        plots.append(_plotly_div(run_fig, "dataset-run-plot"))
    if "best_confidence" in aggregation.peptide_summary.columns and not aggregation.peptide_summary["best_confidence"].dropna().empty:
        conf_fig = go.Figure(
            data=[go.Histogram(x=aggregation.peptide_summary["best_confidence"].dropna(), marker_color="#B8336A")],
            layout=go.Layout(title="Best peptide confidence distribution", template="plotly_white", margin=dict(l=40, r=20, t=48, b=40), height=320),
        )
        plots.append(_plotly_div(conf_fig, "dataset-confidence-plot"))
    if "distinct_charge_count" in aggregation.peptide_summary.columns and not aggregation.peptide_summary.empty:
        charge_fig = go.Figure(
            data=[go.Histogram(x=aggregation.peptide_summary["distinct_charge_count"], marker_color="#FA824C")],
            layout=go.Layout(title="Distinct precursor charges per peptide", template="plotly_white", margin=dict(l=40, r=20, t=48, b=40), height=320),
        )
        plots.append(_plotly_div(charge_fig, "dataset-charge-plot"))

    disagreement = aggregation.peptide_summary.copy()
    if {"best_confidence", "spectral_angle_mean"}.issubset(disagreement.columns):
        disagreement["confidence_spectral_gap"] = (pd.to_numeric(disagreement["best_confidence"], errors="coerce") - pd.to_numeric(disagreement["spectral_angle_mean"], errors="coerce")).abs()
        disagreement = disagreement.sort_values(by=["confidence_spectral_gap", "total_psm_count"], ascending=[False, False], kind="mergesort").head(10)
    else:
        disagreement = pd.DataFrame()

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>proteomics-mcp dataset report</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; margin: 24px; background: #ffffff; color: #172033; }}
    h1, h2 {{ margin: 0 0 12px; }}
    h2 {{ margin-top: 28px; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 12px; }}
    .tile {{ border: 1px solid #d6ddea; border-left: 4px solid #068D9D; border-radius: 6px; padding: 12px 14px; background: #fbfdff; }}
    .label {{ font-size: 12px; color: #556179; text-transform: uppercase; letter-spacing: 0.02em; }}
    .value {{ font-size: 24px; font-weight: 600; margin-top: 4px; }}
    .note {{ font-size: 12px; color: #556179; margin-top: 6px; }}
    .plot {{ width: 100%; }}
    table {{ border-collapse: collapse; width: 100%; margin-top: 10px; }}
    th, td {{ border-bottom: 1px solid #e5e9f2; text-align: left; padding: 8px 10px; vertical-align: top; }}
    th {{ font-size: 12px; color: #556179; text-transform: uppercase; }}
    td {{ font-size: 13px; }}
  </style>
</head>
<body>
  <h1>Dataset report</h1>
  <div class="grid">
    <div class="tile"><div class="label">Total PSMs</div><div class="value">{dataset["total_psm_count"]}</div></div>
    <div class="tile"><div class="label">Distinct peptides</div><div class="value">{dataset["distinct_peptide_count"]}</div></div>
    <div class="tile"><div class="label">Distinct peptidoforms</div><div class="value">{dataset["distinct_peptidoform_count"]}</div></div>
    <div class="tile"><div class="label">Distinct runs</div><div class="value">{dataset["distinct_run_count"]}</div></div>
    <div class="tile"><div class="label">Mean confidence</div><div class="value">{_format_percent(confidence.get("mean_confidence"))}</div></div>
    <div class="tile"><div class="label">Best spectral angle</div><div class="value">{_format_number(spectrum.get("best_spectral_angle"))}</div></div>
  </div>
  <h2>Distributions</h2>
  {''.join(plots)}
  <h2>Top confident peptides</h2>
  {_frame_to_html(peptide_confident[["stripped_sequence", "best_confidence", "total_psm_count", "distinct_peptidoform_count"]])}
  <h2>Top supported peptides</h2>
  {_frame_to_html(peptide_supported[["stripped_sequence", "spectral_angle_max", "matched_ion_count_max", "total_psm_count"]])}
  <h2>Confidence vs spectrum-evidence outliers</h2>
  {_frame_to_html(disagreement[["stripped_sequence", "best_confidence", "spectral_angle_mean", "confidence_spectral_gap"]]) if not disagreement.empty else "<p>No outliers available.</p>"}
  <h2>Run summary</h2>
  {_frame_to_html(run_summary[["experiment_name", "total_psm_count", "distinct_peptide_count", "distinct_peptidoform_count"]]) if not run_summary.empty else "<p>No run grouping available.</p>"}
  {"<h2>SDRF-linked sample summary</h2>" + _frame_to_html(aggregation.sample_summary[["sample_key", "total_psm_count", "distinct_run_count", "distinct_peptide_count"]].head(10)) if aggregation.sample_summary is not None and not aggregation.sample_summary.empty else ""}
</body>
</html>
"""
    if output_path:
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        output_file.write_text(html, encoding="utf-8")
    return html


def _format_percent(value: Any) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    return f"{float(value) * 100:.1f}%"


def _format_number(value: Any) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "NA"
    return f"{float(value):.3f}"


def _frame_to_html(frame: pd.DataFrame) -> str:
    serializable = frame.copy()
    for column in serializable.columns:
        serializable[column] = serializable[column].map(_serialize_cell)
    return serializable.to_html(index=False, border=0, classes="summary-table")


def dataset_groups(aggregation: DatasetAggregation) -> dict[str, Any]:
    return {
        "dataset": {
            "total_psm_count": aggregation.dataset_summary["dataset_metadata"]["total_psm_count"],
            "distinct_peptide_count": aggregation.dataset_summary["dataset_metadata"]["distinct_peptide_count"],
            "distinct_peptidoform_count": aggregation.dataset_summary["dataset_metadata"]["distinct_peptidoform_count"],
        },
        "runs": aggregation.run_summary["experiment_name"].tolist() if "experiment_name" in aggregation.run_summary.columns else [],
        "samples": aggregation.sample_summary["sample_key"].tolist() if aggregation.sample_summary is not None and "sample_key" in aggregation.sample_summary.columns else [],
    }


def get_peptide_summary_row(aggregation: DatasetAggregation, peptide: str) -> dict[str, Any] | None:
    if aggregation.peptide_summary.empty:
        return None
    matches = aggregation.peptide_summary[aggregation.peptide_summary["stripped_sequence"] == peptide]
    if matches.empty:
        return None
    return _json_safe(matches.iloc[0].to_dict())


def get_peptidoform_summary_row(aggregation: DatasetAggregation, peptidoform: str) -> dict[str, Any] | None:
    if aggregation.peptidoform_summary.empty:
        return None
    matches = aggregation.peptidoform_summary[aggregation.peptidoform_summary["peptidoform"] == peptidoform]
    if matches.empty:
        return None
    return _json_safe(matches.iloc[0].to_dict())
