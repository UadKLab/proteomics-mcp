from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import gzip
import inspect
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time
from typing import Any
import uuid
import hashlib

from proteomics_mcp.adapters.casanovo import load_casanovo_predictions
from proteomics_mcp.adapters.fragpipe import load_fragpipe_psms
from proteomics_mcp.adapters.generic import load_generic_candidates
from proteomics_mcp.adapters.instanovo import load_instanovo_predictions
from proteomics_mcp.adapters.maxquant import load_maxquant_psms
from proteomics_mcp.adapters.proteome_discoverer import load_proteome_discoverer_psms
from proteomics_mcp.config import load_config
from proteomics_mcp.core.annotation import annotate_observed_ions
from proteomics_mcp.core.dataset import (
    build_dashboard_summary,
    bin_2d_numeric_distribution,
    summarize_file_coverage,
    summarize_modifications,
    summarize_overlap,
    summarize_quant_quality,
    summarize_missing_values,
    summarize_cv_by_category,
    summarize_categorical_distribution,
    dataset_groups,
    distribution_field_catalog,
    bin_numeric_distribution,
    candidate_models_to_frame,
    get_peptide_summary_row,
    get_peptidoform_summary_row,
    inspection_peptide_table,
    inspection_psm_table,
    load_scored_psm_table,
    merge_scored_enrichment,
    render_peptide_table_html,
    render_psm_table_html,
    render_dashboard_bundle_html,
    render_prediction_accuracy_html,
    render_summary_html,
    render_dataset_report_html,
    summarize_candidate_inspection_frame,
    summarize_dataset_path,
    summarize_numeric_distribution,
    summarize_prediction_accuracy,
    write_dataset_outputs,
)
from proteomics_mcp.core.instanovo_export import export_table_to_instanovo
from proteomics_mcp.core.sdrf import load_sdrf_summary, retrieve_pride_sdrf
from proteomics_mcp.core.scoring import score_candidate_spectrum
from proteomics_mcp.core.spectrum_loader import bruker_d_to_feather, load_spectrum_from_table, mgf_to_feather, mzml_to_feather, retrieve_spectrum_from_usi
from proteomics_mcp.core.visualization import (
    build_spectrum_payload,
    render_annotated_observed_spectrum_html,
    render_observed_spectrum_html,
    render_spectrum_comparison_html,
    render_static_spectrum,
)
from proteomics_mcp.models.schemas import CandidatePSM, ObservedSpectrum, PredictedFragment, PredictedSpectrum, model_dump
from proteomics_mcp.predictors import mock_predicted_spectrum, predict_fragments
from proteomics_mcp.predictors.koina import describe_koina_model, predict_fragments_koina_batch
from proteomics_mcp.predictors.ms2pip import predict_fragments_ms2pip_batch
import pandas as pd
from pandas.errors import EmptyDataError
import pyarrow.feather as pyarrow_feather


def _preview_records(frame: pd.DataFrame, limit: int = 10) -> list[dict[str, Any]]:
    if frame.empty:
        return []
    return json.loads(frame.head(limit).to_json(orient="records"))


def summarize_score_metrics_tool(score_result: dict[str, Any]) -> dict[str, Any]:
    """Return an agent-friendly grouping of score, evidence, error, and candidate fields."""
    candidate = score_result.get("candidate") or {}
    acquisition = score_result.get("acquisition_provenance") or {}
    annotation_summary = score_result.get("annotation_summary") or {}
    matched_peaks = score_result.get("matched_peaks") or []
    matched_ion_type_counts: dict[str, int] = {}
    for peak in matched_peaks:
        if not isinstance(peak, dict):
            continue
        ion_type = peak.get("ion_type")
        if ion_type:
            matched_ion_type_counts[str(ion_type)] = matched_ion_type_counts.get(str(ion_type), 0) + 1
    return {
        "match_quality": {
            "spectral_angle": score_result.get("spectral_angle"),
            "cosine_similarity": score_result.get("cosine_similarity"),
            "pearson_correlation": score_result.get("pearson_correlation"),
            "xcorr": score_result.get("xcorr"),
            "log10_xcorr": score_result.get("log10_xcorr"),
            "log2_xcorr": score_result.get("log2_xcorr"),
            "hyperscore": score_result.get("hyperscore"),
            "log2_hyperscore": score_result.get("log2_hyperscore"),
        },
        "predicted_match_evidence": {
            "matched_ion_count": score_result.get("matched_ion_count"),
            "total_predicted_ions": score_result.get("total_predicted_ions"),
            "matched_ion_fraction": score_result.get("matched_ion_fraction"),
            "explained_intensity_fraction": score_result.get("explained_intensity_fraction"),
            "matched_predicted_ion_type_counts": matched_ion_type_counts,
            "matched_peaks": matched_peaks,
        },
        "observed_annotation_evidence": {
            "precursor_signal_fraction": score_result.get("precursor_signal_fraction"),
            "annotated_peak_count": annotation_summary.get("annotated_peak_count"),
            "ion_type_counts": annotation_summary.get("ion_type_counts"),
            "annotated_intensity_fraction": annotation_summary.get("annotated_intensity_fraction"),
            "annotation_summary": annotation_summary,
        },
        "mass_error": {
            "fragment_mass_error_mean_ppm": score_result.get("fragment_mass_error_mean_ppm"),
            "fragment_mass_error_abs_mean_ppm": score_result.get("fragment_mass_error_abs_mean_ppm"),
            "fragment_mass_error_mean_da": score_result.get("fragment_mass_error_mean_da"),
            "fragment_mass_error_abs_mean_da": score_result.get("fragment_mass_error_abs_mean_da"),
            "precursor_mass_error_ppm": score_result.get("precursor_mass_error_ppm"),
            "precursor_mass_error_da": score_result.get("precursor_mass_error_da"),
        },
        "candidate": {
            "candidate_id": candidate.get("candidate_id"),
            "source_engine": candidate.get("source_engine"),
            "peptidoform": candidate.get("peptidoform"),
            "confidence": candidate.get("confidence"),
            "confidence_percent": candidate.get("confidence") * 100.0 if candidate.get("confidence") is not None else None,
            "raw_score": candidate.get("raw_score"),
            "raw_score_name": candidate.get("raw_score_name"),
            "next_beam": candidate.get("next_beam"),
            "next_beam_margin": candidate.get("next_beam_margin"),
            "token_log_probability_mean": candidate.get("token_log_probability_mean"),
            "token_log_probability_min": candidate.get("token_log_probability_min"),
            "token_confidence_mean": candidate.get("token_confidence_mean"),
            "token_confidence_min": candidate.get("token_confidence_min"),
            "token_confidence_min_position": candidate.get("token_confidence_min_position"),
            "token_confidence_min_token": candidate.get("token_confidence_min_token"),
        },
        "acquisition": acquisition.get("summary") or {},
        "ion_evidence": {
            "matched_ion_count": score_result.get("matched_ion_count"),
            "total_predicted_ions": score_result.get("total_predicted_ions"),
            "matched_ion_fraction": score_result.get("matched_ion_fraction"),
            "explained_intensity_fraction": score_result.get("explained_intensity_fraction"),
            "precursor_signal_fraction": score_result.get("precursor_signal_fraction"),
            "annotated_peak_count": annotation_summary.get("annotated_peak_count"),
            "ion_type_counts": annotation_summary.get("ion_type_counts"),
            "annotated_intensity_fraction": annotation_summary.get("annotated_intensity_fraction"),
            "annotation_summary": annotation_summary,
        },
        "metric_profile": score_result.get("metric_profile"),
    }


def load_candidates_tool(path: str, adapter: str = "instanovo", **kwargs: Any) -> list[dict[str, Any]]:
    return [candidate.model_dump() for candidate in _load_candidate_models(path, adapter, **kwargs)]


# Keys of a generic-adapter column mapping file: the keyword arguments of load_generic_candidates.
GENERIC_MAPPING_KEYS = tuple(
    name
    for name, parameter in inspect.signature(load_generic_candidates).parameters.items()
    if parameter.kind is inspect.Parameter.KEYWORD_ONLY and name != "limit"
)


def load_column_mapping(column_mapping_path: str) -> dict[str, Any]:
    """Read a generic-adapter column mapping: a JSON object such as {"peptide_column": "Sequence", "scan_column": "Scan"}."""
    try:
        mapping = json.loads(Path(column_mapping_path).read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ValueError(f"Column mapping {column_mapping_path} is not valid JSON: {exc}") from exc
    if not isinstance(mapping, dict):
        raise ValueError(f"Column mapping {column_mapping_path} must be a JSON object.")
    unknown = sorted(set(mapping) - set(GENERIC_MAPPING_KEYS))
    if unknown:
        raise ValueError(f"Unknown column mapping keys: {', '.join(unknown)}. Allowed keys: {', '.join(GENERIC_MAPPING_KEYS)}.")
    not_text = sorted(key for key, value in mapping.items() if value is not None and not isinstance(value, str))
    if not_text:
        raise ValueError(f"Column mapping values must be strings (column names or score_type/source_engine): {', '.join(not_text)}.")
    return mapping


def _load_candidate_models(path: str, adapter: str, *, column_mapping_path: str | None = None, **kwargs: Any) -> list[CandidatePSM]:
    if column_mapping_path and adapter != "generic":
        raise ValueError(f"column_mapping_path applies only to adapter='generic', not '{adapter}'.")
    if adapter == "generic":
        mapping = load_column_mapping(column_mapping_path) if column_mapping_path else {}
        kwargs = {**mapping, **{key: value for key, value in kwargs.items() if value is not None}}
        if not kwargs.get("peptide_column"):
            raise ValueError(
                "adapter='generic' needs a column mapping: pass column_mapping_path, a JSON file with at least "
                f"\"peptide_column\". Allowed keys: {', '.join(GENERIC_MAPPING_KEYS)}."
            )
    if adapter == "instanovo":
        return load_instanovo_predictions(
            path,
            limit=kwargs.get("limit"),
            default_experiment_name=kwargs.get("default_experiment_name"),
        )
    if adapter == "casanovo":
        return load_casanovo_predictions(path, limit=kwargs.get("limit"))
    if adapter == "fragpipe":
        return load_fragpipe_psms(path, limit=kwargs.get("limit"))
    if adapter == "maxquant":
        return load_maxquant_psms(path, limit=kwargs.get("limit"))
    if adapter == "proteome_discoverer":
        return load_proteome_discoverer_psms(path, limit=kwargs.get("limit"), sheet_name=kwargs.get("sheet_name", 0))
    if adapter == "generic":
        return load_generic_candidates(path, **kwargs)
    raise ValueError(f"Unsupported candidate adapter: {adapter}")


def _model_from_argument(model: Any, value: Any, argument: str) -> Any:
    """Build a schema object from a tool argument; report bad input as ValueError naming the argument."""
    if not isinstance(value, dict):
        raise ValueError(f"'{argument}' must be a JSON object, got {type(value).__name__}.")
    try:
        return model(**value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Invalid '{argument}': {exc}") from exc


def _serialize_cell(value: Any) -> Any:
    if isinstance(value, (dict, list)):
        return json.dumps(value, sort_keys=True)
    return value


def _write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(json.dumps(row, sort_keys=True))
            handle.write("\n")


def _utc_timestamp() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _batch_status_path(prefix: Path, suffix: str) -> Path:
    return prefix.with_name(prefix.name + suffix)


def _status_update_interval(total: int) -> int:
    if total <= 10:
        return 1
    return max(1, min(25, total // 20))


def _read_batch_status(path: Path) -> dict[str, Any]:
    """The last status written by a batch. While unfinished, also how long ago it was written: a batch that
    stopped updating was most likely interrupted (outputs are written at the end, so it must be run again)."""
    status = json.loads(path.read_text(encoding="utf-8"))
    if status.get("phase") != "completed":
        status["seconds_since_update"] = max(time.time() - path.stat().st_mtime, 0.0)
    return status


def _write_batch_status(
    path: Path,
    *,
    batch_id: str,
    phase: str,
    completed: int,
    total: int,
    succeeded: int,
    failed: int,
    started: float,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    payload = {
        "batch_id": batch_id,
        "phase": phase,
        "completed": int(completed),
        "total": int(total),
        "remaining": max(int(total) - int(completed), 0),
        "succeeded": int(succeeded),
        "failed": int(failed),
        "progress_fraction": (float(completed) / float(total)) if total > 0 else None,
        "started_at": _utc_timestamp(),
        "elapsed_seconds": max(time.time() - started, 0.0),
        "updated_at": _utc_timestamp(),
    }
    if extra:
        payload.update(extra)
    _write_json(path, payload)
    return payload


def _read_csv_if_exists(path: str | None) -> pd.DataFrame:
    if not path or not Path(path).exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path)
    except EmptyDataError:
        return pd.DataFrame()


def _prediction_batch_key(
    *,
    peptidoform: str,
    precursor_charge: int,
    predictor_backend: str,
    model: str,
    max_fragment_charge: int | None,
) -> str:
    payload = {
        "peptidoform": peptidoform,
        "precursor_charge": precursor_charge,
        "predictor_backend": predictor_backend,
        "model": model,
        "max_fragment_charge": max_fragment_charge,
    }
    return hashlib.sha1(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()[:16]


def _predict_fragment_spectrum(config: Any, peptidoform: str, precursor_charge: int) -> PredictedSpectrum:
    return predict_fragments(
        peptidoform,
        precursor_charge,
        backend=config.predictor_backend,
        model=config.predictor_model,
        max_fragment_charge=config.max_fragment_charge,
        koina_server_url=config.koina_server_url,
        koina_collision_energy=config.koina_collision_energy,
        koina_fragmentation_type=config.koina_fragmentation_type,
        koina_instrument_type=config.koina_instrument_type,
        koina_timeout_seconds=config.koina_timeout_seconds,
        koina_cache_dir=config.koina_cache_dir,
    )


def _top_n_frame(frame: pd.DataFrame, top_n: int | None) -> pd.DataFrame:
    if top_n is None or top_n <= 0 or frame.empty:
        return frame
    if "confidence" in frame.columns and pd.to_numeric(frame["confidence"], errors="coerce").notna().any():
        return frame.sort_values(by="confidence", ascending=False, kind="mergesort").head(top_n).reset_index(drop=True)
    if "raw_score" in frame.columns and pd.to_numeric(frame["raw_score"], errors="coerce").notna().any():
        return frame.sort_values(by="raw_score", ascending=False, kind="mergesort").head(top_n).reset_index(drop=True)
    return frame.head(top_n).reset_index(drop=True)


def _load_candidate_frame(
    path: str,
    adapter: str,
    *,
    limit: int | None = None,
    scored_path: str | None = None,
    **adapter_options: Any,
) -> pd.DataFrame:
    candidates = _load_candidate_models(path, adapter, limit=limit, **adapter_options)
    frame = candidate_models_to_frame(candidates)
    if scored_path:
        frame = merge_scored_enrichment(frame, load_scored_psm_table(scored_path))
    return frame


def enrich_frame_with_spectra_metadata(
    frame: pd.DataFrame,
    *,
    spectra_path: str,
    rt_source_field: str = "rt",
) -> tuple[pd.DataFrame, dict[str, Any]]:
    available_columns = set(pyarrow_feather.read_table(spectra_path, columns=[]).schema.names)
    spectra = pd.read_feather(
        spectra_path,
        columns=[
            column
            for column in [
                "source_file",
                "source_name",
                "scan_number",
                rt_source_field,
                "filter_string",
                "scan_resolution",
                "collision_energy",
                "isolation_target_mz",
                "precursor_mz",
                "precursor_charge",
            ]
            if column in available_columns
        ],
    )
    if "scan_number" not in spectra.columns:
        raise ValueError("Spectrum Feather must contain scan_number for metadata enrichment.")

    run_source = None
    if "source_file" in spectra.columns:
        run_source = spectra["source_file"].map(lambda value: Path(str(value)).stem if pd.notna(value) else None)
    elif "source_name" in spectra.columns:
        run_source = spectra["source_name"].map(lambda value: Path(str(value)).stem if pd.notna(value) else None)
    else:
        raise ValueError("Spectrum Feather must contain source_file or source_name for run matching.")

    spectra = spectra.copy()
    spectra["run_key"] = run_source
    spectra = spectra.dropna(subset=["run_key", "scan_number"]).copy()
    spectra["scan_number"] = pd.to_numeric(spectra["scan_number"], errors="coerce")
    spectra = spectra.dropna(subset=["scan_number"]).copy()
    spectra["scan_number"] = spectra["scan_number"].astype(int)
    spectra = spectra.drop_duplicates(subset=["run_key", "scan_number"], keep="first")

    enriched = frame.copy()
    enriched["scan_number"] = pd.to_numeric(enriched.get("scan_number"), errors="coerce")
    enriched = enriched.dropna(subset=["run_key", "scan_number"]).copy()
    enriched["scan_number"] = enriched["scan_number"].astype(int)

    rt_column = rt_source_field if rt_source_field in spectra.columns else None
    join_columns = ["run_key", "scan_number"] + [
        column
        for column in [
            rt_column,
            "filter_string",
            "scan_resolution",
            "collision_energy",
            "isolation_target_mz",
            "precursor_mz",
            "precursor_charge",
        ]
        if column and column in spectra.columns
    ]
    joined = enriched.merge(
        spectra[join_columns].rename(
            columns={
                rt_source_field: "retention_time_seconds",
                "precursor_mz": "mzml_precursor_mz",
                "precursor_charge": "mzml_precursor_charge",
            }
        ),
        on=["run_key", "scan_number"],
        how="left",
        validate="many_to_one",
        suffixes=("", "__mzml"),
    )

    for column in ["retention_time_seconds", "scan_resolution", "collision_energy", "isolation_target_mz", "mzml_precursor_mz", "mzml_precursor_charge"]:
        alt = f"{column}__mzml"
        if alt not in joined.columns:
            continue
        if column in joined.columns:
            joined[column] = joined[column].where(joined[column].notna(), joined[alt])
            joined.drop(columns=[alt], inplace=True)
        else:
            joined.rename(columns={alt: column}, inplace=True)

    # Minimal de novo exports often omit precursor fields. Retain supplied values
    # where present, but recover them from the matched MS2 metadata otherwise.
    for candidate_column, mzml_column in [("precursor_mz", "mzml_precursor_mz"), ("precursor_charge", "mzml_precursor_charge")]:
        if mzml_column not in joined.columns:
            continue
        if candidate_column in joined.columns:
            joined[candidate_column] = joined[candidate_column].where(joined[candidate_column].notna(), joined[mzml_column])
        else:
            joined[candidate_column] = joined[mzml_column]

    rt_seconds = pd.to_numeric(joined["retention_time_seconds"], errors="coerce") if "retention_time_seconds" in joined.columns else pd.Series(dtype=float)
    if "retention_time_minutes" in joined.columns:
        joined["retention_time_minutes"] = pd.to_numeric(joined["retention_time_minutes"], errors="coerce").where(
            pd.to_numeric(joined["retention_time_minutes"], errors="coerce").notna(),
            rt_seconds / 60.0,
        )
    else:
        joined["retention_time_minutes"] = rt_seconds / 60.0
    matched_rows = int(rt_seconds.notna().sum())
    return joined, {
        "total_rows": int(len(joined)),
        "matched_rows": matched_rows,
        "unmatched_rows": int(len(joined) - matched_rows),
        "spectra_path": spectra_path,
        "rt_source_field": rt_source_field,
    }


def inspect_psm_table_tool(
    *,
    path: str,
    adapter: str = "instanovo",
    scored_path: str | None = None,
    sdrf_path: str | None = None,
    pride_project_accession: str | None = None,
    timeout_seconds: float = 20.0,
    limit: int | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    sdrf_summary = _resolve_sdrf_summary(sdrf_path=sdrf_path, pride_project_accession=pride_project_accession, timeout_seconds=timeout_seconds)
    inspection = summarize_candidate_inspection_frame(
        _load_candidate_frame(path, adapter, limit=limit, scored_path=scored_path, **adapter_options),
        sdrf_summary=sdrf_summary,
    )
    return {
        "dataset_metadata": inspection.dataset_summary["dataset_metadata"],
        "confidence_summary": inspection.dataset_summary["confidence_summary"],
        "scoring_summary": inspection.dataset_summary["scoring_summary"],
        "sequence_summary": inspection.dataset_summary["sequence_summary"],
        "run_summary": inspection.dataset_summary["run_summary"],
        "sdrf_summary": inspection.dataset_summary["sdrf_summary"],
        "psm_preview": _preview_records(inspection_psm_table(inspection)),
        "peptide_preview": _preview_records(inspection_peptide_table(inspection)),
    }


def inspect_peptide_table_tool(
    *,
    path: str,
    adapter: str = "instanovo",
    peptide: str | None = None,
    scored_path: str | None = None,
    sdrf_path: str | None = None,
    pride_project_accession: str | None = None,
    timeout_seconds: float = 20.0,
    limit: int | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    sdrf_summary = _resolve_sdrf_summary(sdrf_path=sdrf_path, pride_project_accession=pride_project_accession, timeout_seconds=timeout_seconds)
    inspection = summarize_candidate_inspection_frame(
        _load_candidate_frame(path, adapter, limit=limit, scored_path=scored_path, **adapter_options),
        sdrf_summary=sdrf_summary,
    )
    if peptide is not None:
        peptide_table = inspection_peptide_table(inspection)
        summary = peptide_table[peptide_table["stripped_sequence"] == peptide]
        if summary.empty:
            raise LookupError(f"Peptide '{peptide}' not found in candidate table.")
        return json.loads(summary.head(1).to_json(orient="records"))[0]
    return {
        "dataset_metadata": inspection.dataset_summary["dataset_metadata"],
        "scoring_summary": inspection.dataset_summary["scoring_summary"],
        "peptide_preview": _preview_records(inspection_peptide_table(inspection)),
    }


def list_distribution_fields_tool(
    *,
    path: str,
    adapter: str = "instanovo",
    scored_path: str | None = None,
    limit: int | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    frame = _load_candidate_frame(path, adapter, limit=limit, scored_path=scored_path, **adapter_options)
    return distribution_field_catalog(frame)


def summarize_numeric_field_tool(
    *,
    path: str,
    field: str,
    adapter: str = "instanovo",
    scored_path: str | None = None,
    group_by: str | None = None,
    percentiles: list[float] | None = None,
    weight_field: str | None = None,
    limit: int | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    frame = _load_candidate_frame(path, adapter, limit=limit, scored_path=scored_path, **adapter_options)
    return summarize_numeric_distribution(frame, field=field, group_by=group_by, percentiles=percentiles, weight_field=weight_field)


def bin_numeric_field_tool(
    *,
    path: str,
    field: str,
    adapter: str = "instanovo",
    scored_path: str | None = None,
    group_by: str | None = None,
    bin_edges: list[float] | None = None,
    bin_count: int | None = None,
    value_range: tuple[float, float] | None = None,
    weight_field: str | None = None,
    limit: int | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    frame = _load_candidate_frame(path, adapter, limit=limit, scored_path=scored_path, **adapter_options)
    return bin_numeric_distribution(
        frame,
        field=field,
        group_by=group_by,
        bin_edges=bin_edges,
        bin_count=bin_count,
        value_range=value_range,
        weight_field=weight_field,
    )


def summarize_categorical_field_tool(
    *,
    path: str,
    field: str,
    adapter: str = "instanovo",
    scored_path: str | None = None,
    group_by: str | None = None,
    top_n: int | None = None,
    limit: int | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    frame = _load_candidate_frame(path, adapter, limit=limit, scored_path=scored_path, **adapter_options)
    return summarize_categorical_distribution(frame, field=field, group_by=group_by, top_n=top_n)


def summarize_missing_values_tool(
    *,
    path: str,
    adapter: str = "instanovo",
    scored_path: str | None = None,
    fields: list[str] | None = None,
    limit: int | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    frame = _load_candidate_frame(path, adapter, limit=limit, scored_path=scored_path, **adapter_options)
    return summarize_missing_values(frame, fields=fields)


def summarize_cv_by_category_tool(
    *,
    path: str,
    value_field: str,
    category_field: str,
    adapter: str = "instanovo",
    scored_path: str | None = None,
    top_n: int | None = None,
    limit: int | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    frame = _load_candidate_frame(path, adapter, limit=limit, scored_path=scored_path, **adapter_options)
    return summarize_cv_by_category(frame, value_field=value_field, category_field=category_field, top_n=top_n)


def bin_2d_numeric_fields_tool(
    *,
    path: str,
    x_field: str,
    y_field: str,
    adapter: str = "instanovo",
    scored_path: str | None = None,
    x_bin_edges: list[float] | None = None,
    y_bin_edges: list[float] | None = None,
    x_bin_count: int | None = None,
    y_bin_count: int | None = None,
    weight_field: str | None = None,
    limit: int | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    frame = _load_candidate_frame(path, adapter, limit=limit, scored_path=scored_path, **adapter_options)
    return bin_2d_numeric_distribution(
        frame,
        x_field=x_field,
        y_field=y_field,
        x_bin_edges=x_bin_edges,
        y_bin_edges=y_bin_edges,
        x_bin_count=x_bin_count,
        y_bin_count=y_bin_count,
        weight_field=weight_field,
    )


def summarize_dashboard_tool(
    *,
    path: str,
    adapter: str = "instanovo",
    scored_path: str | None = None,
    rt_field: str | None = None,
    mz_field: str = "precursor_mz",
    confidence_field: str | None = "confidence",
    quant_field: str | None = None,
    category_field: str = "experiment_name",
    limit: int | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    frame = _load_candidate_frame(path, adapter, limit=limit, scored_path=scored_path, **adapter_options)
    return build_dashboard_summary(
        frame,
        rt_field=rt_field,
        mz_field=mz_field,
        confidence_field=confidence_field,
        quant_field=quant_field,
        category_field=category_field,
    )


def summarize_overlap_tool(
    *,
    path: str,
    adapter: str = "instanovo",
    scored_path: str | None = None,
    entity: str = "peptide",
    group_field: str = "run_key",
    limit: int | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    frame = _load_candidate_frame(path, adapter, limit=limit, scored_path=scored_path, **adapter_options)
    return summarize_overlap(frame, entity=entity, group_field=group_field)


def summarize_modifications_tool(
    *,
    path: str,
    adapter: str = "instanovo",
    scored_path: str | None = None,
    group_field: str = "run_key",
    limit: int | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    frame = _load_candidate_frame(path, adapter, limit=limit, scored_path=scored_path, **adapter_options)
    return summarize_modifications(frame, group_field=group_field)


def summarize_prediction_accuracy_tool(
    *,
    path: str,
    target_field: str,
    adapter: str = "instanovo",
    scored_path: str | None = None,
    confidence_field: str | None = "confidence",
    mass_error_field: str | None = None,
    confidence_thresholds: list[float] | None = None,
    reliability_bin_count: int = 10,
    limit: int | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    frame = _load_candidate_frame(path, adapter, limit=limit, scored_path=scored_path, **adapter_options)
    return summarize_prediction_accuracy(
        frame,
        target_field=target_field,
        confidence_field=confidence_field,
        mass_error_field=mass_error_field,
        confidence_thresholds=confidence_thresholds,
        reliability_bin_count=reliability_bin_count,
    )


def render_prediction_accuracy_html_tool(
    summary: dict[str, Any], output_path: str, title: str = "Prediction accuracy"
) -> dict[str, Any]:
    render_prediction_accuracy_html(summary, title=title, output_path=output_path)
    return {"format": "html", "output_path": str(output_path), "title": title}


def summarize_quant_quality_tool(
    *,
    path: str,
    quant_field: str,
    adapter: str = "instanovo",
    scored_path: str | None = None,
    entity: str = "peptide",
    group_field: str = "run_key",
    min_groups_for_cv: int = 2,
    limit: int | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    frame = _load_candidate_frame(path, adapter, limit=limit, scored_path=scored_path, **adapter_options)
    return summarize_quant_quality(frame, quant_field=quant_field, entity=entity, group_field=group_field, min_groups_for_cv=min_groups_for_cv)


def summarize_file_coverage_tool(
    *,
    path: str,
    adapter: str = "instanovo",
    scored_path: str | None = None,
    group_field: str = "run_key",
    mz_field: str = "precursor_mz",
    rt_field: str | None = None,
    confidence_field: str | None = None,
    quant_field: str | None = None,
    limit: int | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    frame = _load_candidate_frame(path, adapter, limit=limit, scored_path=scored_path, **adapter_options)
    return summarize_file_coverage(
        frame,
        group_field=group_field,
        mz_field=mz_field,
        rt_field=rt_field,
        confidence_field=confidence_field,
        quant_field=quant_field,
    )


def render_summary_html_tool(summary: dict[str, Any], output_path: str, title: str = "Summary") -> dict[str, Any]:
    render_summary_html(summary, title=title, output_path=output_path)
    return {"format": "html", "output_path": str(output_path), "title": title}


def render_dashboard_html_tool(summary: dict[str, Any], output_path: str, title: str = "Dashboard") -> dict[str, Any]:
    render_dashboard_bundle_html(summary, title=title, output_path=output_path)
    return {"format": "html", "output_path": str(output_path), "title": title}


def build_id_dashboard_tool(
    *,
    path: str,
    spectra_path: str | None = None,
    output_prefix: str,
    adapter: str = "instanovo",
    scored_path: str | None = None,
    top_n: int | None = None,
    rt_source_field: str = "rt",
    confidence_field: str | None = "confidence",
    mz_field: str = "precursor_mz",
    limit: int | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    frame = _top_n_frame(_load_candidate_frame(path, adapter, limit=limit, scored_path=scored_path, **adapter_options), top_n)
    if spectra_path:
        enriched, join_summary = enrich_frame_with_spectra_metadata(frame, spectra_path=spectra_path, rt_source_field=rt_source_field)
        rt_field_for_summary: str | None = "retention_time_minutes"
    else:
        enriched = frame
        join_summary = None
        rt_field_for_summary = None

    summary = {
        "dashboard": build_dashboard_summary(
            enriched,
            rt_field=rt_field_for_summary,
            mz_field=mz_field,
            confidence_field=confidence_field,
            category_field="run_key",
        ),
        "overlap_peptide": summarize_overlap(enriched, entity="peptide", group_field="run_key"),
        "modifications": summarize_modifications(enriched, group_field="run_key"),
        "file_coverage": summarize_file_coverage(
            enriched,
            group_field="run_key",
            mz_field=mz_field,
            rt_field=rt_field_for_summary,
            confidence_field=confidence_field,
        ),
        "metadata_join": join_summary,
    }

    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    enriched_csv = prefix.with_name(prefix.name + "_enriched.csv")
    serializable = enriched.copy()
    for column in serializable.columns:
        serializable[column] = serializable[column].map(_serialize_cell)
    serializable.to_csv(enriched_csv, index=False)

    summary_json = prefix.with_name(prefix.name + "_summary.json")
    summary_json.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    dashboard_html = prefix.with_name(prefix.name + "_dashboard.html")
    render_dashboard_bundle_html(summary, title=f"{prefix.name} dashboard", output_path=dashboard_html)

    return {
        "top_n": top_n,
        "output_files": {
            "enriched_csv": str(enriched_csv),
            "summary_json": str(summary_json),
            "dashboard_html": str(dashboard_html),
        },
        "metadata_join": join_summary,
        "dataset_metadata": summary["dashboard"]["overview"],
    }


def build_prediction_batch_tool(
    *,
    path: str,
    output_prefix: str,
    spectra_path: str | None = None,
    adapter: str = "instanovo",
    predictor_backend: str | None = None,
    predictor_model: str | None = None,
    ms2pip_model: str | None = None,
    max_fragment_charge: int | None = None,
    top_n: int | None = None,
    limit: int | None = None,
    config_path: str | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    started = time.time()
    started_at = _utc_timestamp()
    config = load_config(
        config_path,
        {
            "predictor_backend": predictor_backend,
            "predictor_model": predictor_model,
            "ms2pip_model": ms2pip_model,
            "max_fragment_charge": max_fragment_charge,
        },
    )
    candidates = _load_candidate_models(path, adapter, limit=limit, **adapter_options)
    frame = _top_n_frame(candidate_models_to_frame(candidates), top_n)
    if frame.empty:
        raise ValueError("No candidates available for prediction batch.")
    metadata_join = None
    if spectra_path:
        frame, metadata_join = enrich_frame_with_spectra_metadata(frame, spectra_path=spectra_path)

    run_id = f"pred_{time.strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)

    work_rows = (
        frame[["candidate_id", "run_key", "experiment_name", "scan_number", "spectrum_id", "source_engine", "peptidoform", "stripped_sequence", "precursor_charge"]]
        .copy()
    )
    work_rows["precursor_charge"] = pd.to_numeric(work_rows["precursor_charge"], errors="coerce")
    work_rows = work_rows.dropna(subset=["peptidoform", "precursor_charge"]).copy()
    if work_rows.empty:
        raise ValueError("No candidates have a usable precursor charge. Supply precursor_charge upstream or provide spectra_path for metadata enrichment.")
    work_rows["precursor_charge"] = work_rows["precursor_charge"].astype(int)
    work_rows["prediction_id"] = work_rows.apply(
        lambda row: "pred:" + _prediction_batch_key(
            peptidoform=str(row["peptidoform"]),
            precursor_charge=int(row["precursor_charge"]),
            predictor_backend=config.predictor_backend,
            model=config.predictor_model,
            max_fragment_charge=config.max_fragment_charge,
        ),
        axis=1,
    )

    unique_rows = work_rows.drop_duplicates(subset=["prediction_id"], keep="first").reset_index(drop=True)

    prediction_rows: list[dict[str, Any]] = []
    fragment_rows: list[dict[str, Any]] = []
    failure_rows: list[dict[str, Any]] = []
    log_rows: list[dict[str, Any]] = []
    status_json = _batch_status_path(prefix, "_prediction_status.json")
    update_interval = _status_update_interval(int(len(unique_rows)))
    predictions_csv = prefix.with_name(prefix.name + "_predictions.csv")
    fragments_csv = prefix.with_name(prefix.name + "_prediction_fragments.csv")
    map_csv = prefix.with_name(prefix.name + "_prediction_map.csv")
    failures_csv = prefix.with_name(prefix.name + "_prediction_failures.csv")
    log_jsonl = prefix.with_name(prefix.name + "_prediction_log.jsonl")
    manifest_json = prefix.with_name(prefix.name + "_prediction_manifest.json")
    output_files = {
        "predictions_csv": str(predictions_csv),
        "fragments_csv": str(fragments_csv),
        "prediction_map_csv": str(map_csv),
        "failures_csv": str(failures_csv),
        "log_jsonl": str(log_jsonl),
        "status_json": str(status_json),
        "manifest_json": str(manifest_json),
    }
    # Written now, so the status tools can follow the batch while it runs or after it was interrupted;
    # rewritten with the results (phase "completed") at the end.
    _write_json(manifest_json, {"batch_id": run_id, "phase": "running", "input_path": path, "spectra_path": spectra_path,
                                "adapter": adapter, "predictor_backend": config.predictor_backend,
                                "requested_model": config.predictor_model, "top_n": top_n,
                                "unique_prediction_requests": int(len(unique_rows)), "output_files": output_files})

    _write_batch_status(
        status_json,
        batch_id=run_id,
        phase="predicting",
        completed=0,
        total=int(len(unique_rows)),
        succeeded=0,
        failed=0,
        started=started,
        extra={
            "started_at": started_at,
            "input_path": path,
            "spectra_path": spectra_path,
            "metadata_join": metadata_join,
            "adapter": adapter,
            "predictor_backend": config.predictor_backend,
            "requested_model": config.predictor_model,
            "requested_max_fragment_charge": config.max_fragment_charge,
            "candidate_rows_considered": int(len(frame)),
            "candidate_rows_mappable": int(len(work_rows)),
            "unique_prediction_requests": int(len(unique_rows)),
        },
    )

    predicted_by_index: dict[int, PredictedSpectrum] = {}
    batch_prediction_errors: dict[int, Exception] = {}
    if config.predictor_backend in {"koina", "ms2pip"}:
        def predict_subset(indices: list[int]) -> None:
            if not indices:
                return
            try:
                requests = [
                    (str(unique_rows.iloc[index]["peptidoform"]), int(unique_rows.iloc[index]["precursor_charge"]))
                    for index in indices
                ]
                if config.predictor_backend == "koina":
                    predicted_batch = predict_fragments_koina_batch(
                        requests,
                        model=config.predictor_model,
                        max_fragment_charge=config.max_fragment_charge,
                        server_url=config.koina_server_url,
                        collision_energy=config.koina_collision_energy,
                        fragmentation_type=config.koina_fragmentation_type,
                        instrument_type=config.koina_instrument_type,
                        timeout_seconds=config.koina_timeout_seconds,
                        batch_size=config.koina_batch_size,
                        cache_dir=config.koina_cache_dir,
                    )
                else:
                    predicted_batch = predict_fragments_ms2pip_batch(
                        requests,
                        model=config.predictor_model,
                        max_fragment_charge=config.max_fragment_charge,
                    )
                predicted_by_index.update(dict(zip(indices, predicted_batch, strict=True)))
            except Exception as exc:
                if len(indices) == 1:
                    batch_prediction_errors[indices[0]] = exc
                    return
                midpoint = len(indices) // 2
                predict_subset(indices[:midpoint])
                predict_subset(indices[midpoint:])

        predict_subset(list(range(len(unique_rows))))

    for index, row in unique_rows.iterrows():
        event_base = {
            "batch_id": run_id,
            "prediction_id": row["prediction_id"],
            "peptidoform": str(row["peptidoform"]),
            "precursor_charge": int(row["precursor_charge"]),
            "predictor_backend": config.predictor_backend,
            "requested_model": config.predictor_model,
            "requested_max_fragment_charge": config.max_fragment_charge,
            "completed": int(index + 1),
            "total": int(len(unique_rows)),
            "timestamp": _utc_timestamp(),
        }
        attempt_started = time.time()
        try:
            if index in batch_prediction_errors:
                error = batch_prediction_errors[index]
                backend_name = "Koina" if config.predictor_backend == "koina" else "MS2PIP"
                raise RuntimeError(f"{backend_name} prediction failed: {error}") from error
            predicted = (
                predicted_by_index[index]
                if config.predictor_backend in {"koina", "ms2pip"}
                else _predict_fragment_spectrum(config, str(row["peptidoform"]), int(row["precursor_charge"]))
            )
            elapsed = time.time() - attempt_started
            fragment_count = len(predicted.fragments)
            metadata = predicted.metadata or {}
            prediction_rows.append(
                {
                    "prediction_id": row["prediction_id"],
                    "peptidoform": predicted.peptidoform,
                    "stripped_sequence": str(row["stripped_sequence"]) if pd.notna(row["stripped_sequence"]) else None,
                    "precursor_charge": predicted.precursor_charge,
                    "predictor_backend": config.predictor_backend,
                    "predictor_name": predicted.predictor_name,
                    "predictor_version": predicted.predictor_version,
                    "requested_model": config.predictor_model,
                    "resolved_model": predicted.model,
                    "max_fragment_charge": metadata.get("max_fragment_charge", config.max_fragment_charge),
                    "fragment_charge_rule": metadata.get("fragment_charge_rule"),
                    "available_fragment_charges": metadata.get("available_fragment_charges"),
                    "fragment_count": fragment_count,
                    "prediction_status": "ok",
                    "prediction_seconds": elapsed,
                }
            )
            for fragment in predicted.fragments:
                fragment_rows.append(
                    {
                        "prediction_id": row["prediction_id"],
                        "annotation": fragment.annotation,
                        "ion_type": fragment.ion_type,
                        "ordinal": fragment.ordinal,
                        "charge": fragment.charge,
                        "mz": fragment.mz,
                        "intensity": fragment.intensity,
                        "neutral_loss": fragment.neutral_loss,
                    }
                )
            log_rows.append(
                {
                    **event_base,
                    "event": "prediction_completed",
                    "status": "ok",
                    "predictor_name": predicted.predictor_name,
                    "predictor_version": predicted.predictor_version,
                    "resolved_model": predicted.model,
                    "fragment_count": fragment_count,
                    "elapsed_seconds": elapsed,
                }
            )
        except Exception as exc:
            elapsed = time.time() - attempt_started
            failure_rows.append(
                {
                    "prediction_id": row["prediction_id"],
                    "peptidoform": str(row["peptidoform"]),
                    "precursor_charge": int(row["precursor_charge"]),
                    "predictor_backend": config.predictor_backend,
                    "requested_model": config.predictor_model,
                    "requested_max_fragment_charge": config.max_fragment_charge,
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                    "prediction_seconds": elapsed,
                }
            )
            log_rows.append(
                {
                    **event_base,
                    "event": "prediction_failed",
                    "status": "failed",
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                    "elapsed_seconds": elapsed,
                }
            )

        completed = int(index + 1)
        if completed == int(len(unique_rows)) or completed % update_interval == 0:
            _write_batch_status(
                status_json,
                batch_id=run_id,
                phase="predicting",
                completed=completed,
                total=int(len(unique_rows)),
                succeeded=int(len(prediction_rows)),
                failed=int(len(failure_rows)),
                started=started,
                extra={
                    "started_at": started_at,
                    "input_path": path,
                    "spectra_path": spectra_path,
                    "metadata_join": metadata_join,
                    "adapter": adapter,
                    "predictor_backend": config.predictor_backend,
                    "requested_model": config.predictor_model,
                    "requested_max_fragment_charge": config.max_fragment_charge,
                    "candidate_rows_considered": int(len(frame)),
                    "candidate_rows_mappable": int(len(work_rows)),
                    "unique_prediction_requests": int(len(unique_rows)),
                    "latest_prediction_id": row["prediction_id"],
                    "latest_peptidoform": str(row["peptidoform"]),
                },
            )

    prediction_map = work_rows[
        [
            "candidate_id",
            "prediction_id",
            "source_engine",
            "run_key",
            "experiment_name",
            "scan_number",
            "spectrum_id",
            "peptidoform",
            "stripped_sequence",
            "precursor_charge",
        ]
    ].copy()

    pd.DataFrame(prediction_rows).to_csv(predictions_csv, index=False)
    pd.DataFrame(fragment_rows).to_csv(fragments_csv, index=False)
    prediction_map.to_csv(map_csv, index=False)
    pd.DataFrame(failure_rows).to_csv(failures_csv, index=False)
    _write_jsonl(log_jsonl, log_rows)

    manifest = {
        "batch_id": run_id,
        "phase": "completed",
        "input_path": path,
        "spectra_path": spectra_path,
        "metadata_join": metadata_join,
        "adapter": adapter,
        "column_mapping_path": adapter_options.get("column_mapping_path"),
        "predictor_backend": config.predictor_backend,
        "requested_model": config.predictor_model,
        "requested_max_fragment_charge": config.max_fragment_charge,
        "top_n": top_n,
        "candidates_read": int(len(candidates)),
        "candidate_rows_considered": int(len(frame)),
        "candidate_rows_mappable": int(len(work_rows)),
        "unique_prediction_requests": int(len(unique_rows)),
        "predictions_succeeded": int(len(prediction_rows)),
        "predictions_failed": int(len(failure_rows)),
        "batch_seconds": time.time() - started,
        "output_files": output_files,
    }
    _write_json(manifest_json, manifest)
    _write_batch_status(
        status_json,
        batch_id=run_id,
        phase="completed",
        completed=int(len(unique_rows)),
        total=int(len(unique_rows)),
        succeeded=int(len(prediction_rows)),
        failed=int(len(failure_rows)),
        started=started,
        extra={
            "started_at": started_at,
            "input_path": path,
            "spectra_path": spectra_path,
            "metadata_join": metadata_join,
            "adapter": adapter,
            "predictor_backend": config.predictor_backend,
            "requested_model": config.predictor_model,
            "requested_max_fragment_charge": config.max_fragment_charge,
            "candidate_rows_considered": int(len(frame)),
            "candidate_rows_mappable": int(len(work_rows)),
            "unique_prediction_requests": int(len(unique_rows)),
            "output_files": manifest["output_files"],
        },
    )
    return manifest


def inspect_prediction_batch_tool(manifest_path: str) -> dict[str, Any]:
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    output_files = manifest.get("output_files") or {}
    predictions = _read_csv_if_exists(output_files.get("predictions_csv"))
    fragments = _read_csv_if_exists(output_files.get("fragments_csv"))
    failures = _read_csv_if_exists(output_files.get("failures_csv"))
    prediction_map = _read_csv_if_exists(output_files.get("prediction_map_csv"))

    fragment_counts = pd.to_numeric(predictions.get("fragment_count", pd.Series(dtype=float)), errors="coerce").dropna()
    charge_counts = (
        pd.to_numeric(predictions.get("precursor_charge", pd.Series(dtype=float)), errors="coerce").dropna().astype(int).value_counts().sort_index()
    )
    model_counts = predictions.get("resolved_model", pd.Series(dtype=object)).dropna().astype(str).value_counts()
    ion_counts = fragments.get("ion_type", pd.Series(dtype=object)).dropna().astype(str).value_counts()
    error_counts = failures.get("error_type", pd.Series(dtype=object)).dropna().astype(str).value_counts()

    return {
        "manifest": manifest,
        "overview": {
            "unique_prediction_requests": int(manifest.get("unique_prediction_requests", 0)),
            "predictions_succeeded": int(manifest.get("predictions_succeeded", 0)),
            "predictions_failed": int(manifest.get("predictions_failed", 0)),
            "mapped_candidate_rows": int(len(prediction_map)),
            "distinct_prediction_ids_in_map": int(prediction_map.get("prediction_id", pd.Series(dtype=object)).nunique()) if not prediction_map.empty else 0,
        },
        "fragment_count_summary": {
            "count": int(len(fragment_counts)),
            "min": float(fragment_counts.min()) if not fragment_counts.empty else None,
            "max": float(fragment_counts.max()) if not fragment_counts.empty else None,
            "mean": float(fragment_counts.mean()) if not fragment_counts.empty else None,
            "median": float(fragment_counts.median()) if not fragment_counts.empty else None,
        },
        "precursor_charge_counts": {str(int(k)): int(v) for k, v in charge_counts.items()},
        "resolved_model_counts": {str(k): int(v) for k, v in model_counts.items()},
        "fragment_ion_type_counts": {str(k): int(v) for k, v in ion_counts.items()},
        "failure_type_counts": {str(k): int(v) for k, v in error_counts.items()},
        "preview": {
            "predictions": predictions.head(20).to_dict(orient="records") if not predictions.empty else [],
            "failures": failures.head(20).to_dict(orient="records") if not failures.empty else [],
            "prediction_map": prediction_map.head(20).to_dict(orient="records") if not prediction_map.empty else [],
        },
    }


def inspect_prediction_batch_status_tool(manifest_path: str) -> dict[str, Any]:
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    status_path = (manifest.get("output_files") or {}).get("status_json")
    if status_path and Path(status_path).exists():
        return _read_batch_status(Path(status_path))
    return {
        "batch_id": manifest.get("batch_id"),
        "phase": "completed",
        "completed": int(manifest.get("unique_prediction_requests", 0)),
        "total": int(manifest.get("unique_prediction_requests", 0)),
        "remaining": 0,
        "succeeded": int(manifest.get("predictions_succeeded", 0)),
        "failed": int(manifest.get("predictions_failed", 0)),
        "progress_fraction": 1.0 if int(manifest.get("unique_prediction_requests", 0)) > 0 else None,
        "batch_seconds": manifest.get("batch_seconds"),
    }


def _clean_optional_text(value: Any) -> str | None:
    if value is None:
        return None
    if isinstance(value, float) and pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() == "nan":
        return None
    return text


def _clean_optional_int(value: Any) -> int | None:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _clean_optional_float(value: Any) -> float | None:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _load_prediction_batch_artifacts(manifest_path: str) -> tuple[dict[str, Any], pd.DataFrame, dict[str, PredictedSpectrum]]:
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    output_files = manifest.get("output_files") or {}
    prediction_map = _read_csv_if_exists(output_files.get("prediction_map_csv"))
    predictions = _read_csv_if_exists(output_files.get("predictions_csv"))
    fragments = _read_csv_if_exists(output_files.get("fragments_csv"))

    spectra: dict[str, PredictedSpectrum] = {}
    prediction_rows = predictions.set_index("prediction_id").to_dict(orient="index") if not predictions.empty else {}
    for prediction_id, meta in prediction_rows.items():
        fragment_subset = fragments[fragments["prediction_id"] == prediction_id] if not fragments.empty else pd.DataFrame()
        predicted_fragments: list[PredictedFragment] = []
        for _, fragment_row in fragment_subset.iterrows():
            predicted_fragments.append(
                PredictedFragment(
                    mz=float(fragment_row["mz"]),
                    intensity=float(fragment_row["intensity"]),
                    annotation=str(fragment_row["annotation"]),
                    ion_type=_clean_optional_text(fragment_row.get("ion_type")),
                    ordinal=_clean_optional_int(fragment_row.get("ordinal")),
                    charge=_clean_optional_int(fragment_row.get("charge")) or 1,
                    neutral_loss=_clean_optional_text(fragment_row.get("neutral_loss")),
                )
            )
        spectra[str(prediction_id)] = PredictedSpectrum(
            peptidoform=str(meta["peptidoform"]),
            precursor_charge=int(meta["precursor_charge"]),
            fragments=predicted_fragments,
            predictor_name=_clean_optional_text(meta.get("predictor_name")) or str(meta.get("predictor_backend") or "unknown"),
            predictor_version=_clean_optional_text(meta.get("predictor_version")),
            model=_clean_optional_text(meta.get("resolved_model")),
            metadata={
                "requested_model": _clean_optional_text(meta.get("requested_model")),
                "max_fragment_charge": _clean_optional_int(meta.get("max_fragment_charge")),
                "fragment_charge_rule": _clean_optional_text(meta.get("fragment_charge_rule")),
            },
        )
    return manifest, prediction_map, spectra


def score_from_prediction_batch_tool(
    *,
    manifest_path: str,
    spectra_path: str,
    output_prefix: str,
    tolerance: float | None = None,
    tolerance_unit: str | None = None,
    annotation_ion_types: str | None = None,
    annotation_max_isotope: int | None = None,
    annotation_max_ion_charge: int | None = None,
    annotation_neutral_losses: bool | None = None,
    config_path: str | None = None,
) -> dict[str, Any]:
    manifest, prediction_map, predicted_lookup = _load_prediction_batch_artifacts(manifest_path)
    if manifest.get("phase", "completed") != "completed":
        raise ValueError(
            f"Prediction batch {manifest.get('batch_id')} has not finished (phase '{manifest.get('phase')}'); check it with "
            "inspect_prediction_batch_status, and build it again if the run was interrupted."
        )
    if prediction_map.empty:
        raise ValueError("Prediction batch contains no candidate-to-prediction map.")

    config = load_config(
        config_path,
        {
            "tolerance": tolerance,
            "tolerance_unit": tolerance_unit,
            "annotation_ion_types": annotation_ion_types,
            "annotation_max_isotope": annotation_max_isotope,
            "annotation_max_ion_charge": annotation_max_ion_charge,
            "annotation_neutral_losses": annotation_neutral_losses,
        },
    )

    input_path = manifest.get("input_path")
    adapter = manifest.get("adapter", "instanovo")
    candidates = _load_candidate_models(str(input_path), str(adapter), column_mapping_path=manifest.get("column_mapping_path"))
    candidate_by_id = {str(candidate.candidate_id): candidate for candidate in candidates}

    batch_id = f"score_{time.strftime('%Y%m%d_%H%M%S')}_{uuid.uuid4().hex[:8]}"
    started = time.time()
    started_at = _utc_timestamp()
    scored_rows: list[dict[str, Any]] = []
    failure_rows: list[dict[str, Any]] = []
    log_rows: list[dict[str, Any]] = []
    prefix = Path(output_prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    status_json = _batch_status_path(prefix, "_scoring_status.json")
    update_interval = _status_update_interval(int(len(prediction_map)))
    scored_csv = prefix.with_name(prefix.name + "_scored.csv")
    failures_csv = prefix.with_name(prefix.name + "_scoring_failures.csv")
    log_jsonl = prefix.with_name(prefix.name + "_scoring_log.jsonl")
    manifest_json = prefix.with_name(prefix.name + "_scoring_manifest.json")
    output_files = {
        "scored_csv": str(scored_csv),
        "failures_csv": str(failures_csv),
        "log_jsonl": str(log_jsonl),
        "status_json": str(status_json),
        "manifest_json": str(manifest_json),
    }
    _write_json(manifest_json, {"batch_id": batch_id, "phase": "running", "prediction_manifest_path": manifest_path,
                                "spectra_path": spectra_path, "candidate_rows_requested": int(len(prediction_map)),
                                "output_files": output_files})

    _write_batch_status(
        status_json,
        batch_id=batch_id,
        phase="scoring",
        completed=0,
        total=int(len(prediction_map)),
        succeeded=0,
        failed=0,
        started=started,
        extra={
            "started_at": started_at,
            "manifest_path": manifest_path,
            "spectra_path": spectra_path,
            "candidate_rows_requested": int(len(prediction_map)),
            "prediction_batch_id": manifest.get("batch_id"),
        },
    )

    for index, row in prediction_map.iterrows():
        candidate_id = str(row["candidate_id"])
        prediction_id = str(row["prediction_id"])
        candidate = candidate_by_id.get(candidate_id)
        event_base = {
            "batch_id": batch_id,
            "candidate_id": candidate_id,
            "prediction_id": prediction_id,
            "completed": int(index + 1),
            "total": int(len(prediction_map)),
            "timestamp": _utc_timestamp(),
        }
        attempt_started = time.time()
        try:
            if candidate is None:
                raise KeyError(f"Candidate '{candidate_id}' not found in original input.")
            predicted = predicted_lookup.get(prediction_id)
            if predicted is None:
                raise KeyError(f"Prediction '{prediction_id}' not found in batch artifacts.")
            observed = load_spectrum_from_table(
                spectra_path,
                scan_number=candidate.scan_number,
                spectrum_id=candidate.spectrum_id,
                experiment_name=candidate.experiment_name,
            )
            result = score_candidate_spectrum(
                candidate=candidate,
                observed=observed,
                predicted=predicted,
                tolerance=config.tolerance,
                tolerance_unit=config.tolerance_unit,
                metric_profile=config.metric_profile,
                annotation_ion_types=config.annotation_ion_types,
                annotation_max_isotope=config.annotation_max_isotope,
                annotation_max_ion_charge=config.annotation_max_ion_charge,
                annotation_neutral_losses=config.annotation_neutral_losses,
            )
            acquisition_summary = result.acquisition_provenance.get("summary", {})
            scored_rows.append(
                {
                    "candidate_id": candidate.candidate_id,
                    "prediction_id": prediction_id,
                    "source_engine": candidate.source_engine,
                    "source_engine_version": candidate.source_engine_version,
                    "experiment_name": candidate.experiment_name,
                    "scan_number": candidate.scan_number,
                    "spectrum_id": candidate.spectrum_id,
                    "usi": candidate.usi,
                    "peptidoform": candidate.peptidoform,
                    "precursor_mz": candidate.precursor_mz,
                    "precursor_charge": candidate.precursor_charge,
                    "rank": candidate.rank,
                    "confidence": candidate.confidence,
                    "raw_score": candidate.raw_score,
                    "raw_score_name": candidate.raw_score_name,
                    "next_beam": candidate.next_beam,
                    "next_beam_margin": candidate.next_beam_margin,
                    "token_log_probability_mean": candidate.token_log_probability_mean,
                    "token_log_probability_min": candidate.token_log_probability_min,
                    "token_confidence_mean": candidate.token_confidence_mean,
                    "token_confidence_min": candidate.token_confidence_min,
                    "token_confidence_min_position": candidate.token_confidence_min_position,
                    "token_confidence_min_token": candidate.token_confidence_min_token,
                    "spectral_angle": result.spectral_angle,
                    "cosine_similarity": result.cosine_similarity,
                    "pearson_correlation": result.pearson_correlation,
                    "xcorr": result.xcorr,
                    "log10_xcorr": result.log10_xcorr,
                    "log2_xcorr": result.log2_xcorr,
                    "hyperscore": result.hyperscore,
                    "log2_hyperscore": result.log2_hyperscore,
                    "matched_ion_count": result.matched_ion_count,
                    "total_predicted_ions": result.total_predicted_ions,
                    "matched_ion_fraction": result.matched_ion_fraction,
                    "explained_intensity_fraction": result.explained_intensity_fraction,
                    "precursor_signal_fraction": result.precursor_signal_fraction,
                    "fragment_mass_error_mean_ppm": result.fragment_mass_error_mean_ppm,
                    "fragment_mass_error_abs_mean_ppm": result.fragment_mass_error_abs_mean_ppm,
                    "fragment_mass_error_mean_da": result.fragment_mass_error_mean_da,
                    "fragment_mass_error_abs_mean_da": result.fragment_mass_error_abs_mean_da,
                    "precursor_mass_error_ppm": result.precursor_mass_error_ppm,
                    "precursor_mass_error_da": result.precursor_mass_error_da,
                    "fragmentation_method": acquisition_summary.get("fragmentation_method"),
                    "isolation_target_mz": acquisition_summary.get("isolation_target_mz"),
                    "retention_time_seconds": acquisition_summary.get("retention_time_seconds"),
                    "retention_time_minutes": acquisition_summary.get("retention_time_minutes"),
                    "collision_energy": acquisition_summary.get("collision_energy"),
                    "faims_compensation_voltage": acquisition_summary.get("faims_compensation_voltage"),
                    "injection_time": acquisition_summary.get("injection_time"),
                    "scan_header": acquisition_summary.get("scan_header"),
                    "scan_resolution": acquisition_summary.get("scan_resolution"),
                    "instrument_model": acquisition_summary.get("instrument_model"),
                    "detector_type": acquisition_summary.get("detector_type"),
                    "annotated_peak_count": result.annotation_summary.get("annotated_peak_count"),
                    "ion_type_counts": result.annotation_summary.get("ion_type_counts"),
                    "annotated_intensity_fraction": result.annotation_summary.get("annotated_intensity_fraction"),
                    "metric_profile": result.metric_profile,
                    "predictor_name": predicted.predictor_name,
                    "predictor_model": predicted.model,
                }
            )
            log_rows.append(
                {
                    **event_base,
                    "event": "scoring_completed",
                    "status": "ok",
                    "matched_ion_count": result.matched_ion_count,
                    "spectral_angle": result.spectral_angle,
                    "elapsed_seconds": time.time() - attempt_started,
                }
            )
        except Exception as exc:
            failure_rows.append(
                {
                    "candidate_id": candidate_id,
                    "prediction_id": prediction_id,
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                    "elapsed_seconds": time.time() - attempt_started,
                }
            )
            log_rows.append(
                {
                    **event_base,
                    "event": "scoring_failed",
                    "status": "failed",
                    "error_type": type(exc).__name__,
                    "error": str(exc),
                    "elapsed_seconds": time.time() - attempt_started,
                }
            )

        completed = int(index + 1)
        if completed == int(len(prediction_map)) or completed % update_interval == 0:
            _write_batch_status(
                status_json,
                batch_id=batch_id,
                phase="scoring",
                completed=completed,
                total=int(len(prediction_map)),
                succeeded=int(len(scored_rows)),
                failed=int(len(failure_rows)),
                started=started,
                extra={
                    "started_at": started_at,
                    "manifest_path": manifest_path,
                    "spectra_path": spectra_path,
                    "candidate_rows_requested": int(len(prediction_map)),
                    "prediction_batch_id": manifest.get("batch_id"),
                    "latest_candidate_id": candidate_id,
                    "latest_prediction_id": prediction_id,
                },
            )

    pd.DataFrame(scored_rows).to_csv(scored_csv, index=False)
    pd.DataFrame(failure_rows).to_csv(failures_csv, index=False)
    _write_jsonl(log_jsonl, log_rows)

    scoring_manifest = {
        "batch_id": batch_id,
        "phase": "completed",
        "prediction_manifest_path": manifest_path,
        "spectra_path": spectra_path,
        "candidate_rows_requested": int(len(prediction_map)),
        "rows_scored": int(len(scored_rows)),
        "rows_failed": int(len(failure_rows)),
        "batch_seconds": time.time() - started,
        "metric_groups": {
            "match_quality": ["spectral_angle", "cosine_similarity", "pearson_correlation", "xcorr", "log10_xcorr", "log2_xcorr", "hyperscore", "log2_hyperscore"],
            "predicted_match_evidence": ["matched_ion_count", "total_predicted_ions", "matched_ion_fraction", "explained_intensity_fraction"],
            "observed_annotation_evidence": ["precursor_signal_fraction", "annotated_peak_count", "ion_type_counts", "annotated_intensity_fraction"],
            "mass_error": ["fragment_mass_error_mean_ppm", "fragment_mass_error_abs_mean_ppm", "fragment_mass_error_mean_da", "fragment_mass_error_abs_mean_da", "precursor_mass_error_ppm", "precursor_mass_error_da"],
        },
        "output_files": output_files,
    }
    _write_json(manifest_json, scoring_manifest)
    _write_batch_status(
        status_json,
        batch_id=batch_id,
        phase="completed",
        completed=int(len(prediction_map)),
        total=int(len(prediction_map)),
        succeeded=int(len(scored_rows)),
        failed=int(len(failure_rows)),
        started=started,
        extra={
            "started_at": started_at,
            "manifest_path": manifest_path,
            "spectra_path": spectra_path,
            "candidate_rows_requested": int(len(prediction_map)),
            "prediction_batch_id": manifest.get("batch_id"),
            "output_files": scoring_manifest["output_files"],
        },
    )
    return scoring_manifest


def inspect_scored_batch_tool(manifest_path: str) -> dict[str, Any]:
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    output_files = manifest.get("output_files") or {}
    scored = _read_csv_if_exists(output_files.get("scored_csv"))
    failures = _read_csv_if_exists(output_files.get("failures_csv"))

    def _series_summary(frame: pd.DataFrame, column: str) -> dict[str, Any]:
        values = pd.to_numeric(frame.get(column, pd.Series(dtype=float)), errors="coerce").dropna()
        return {
            "count": int(len(values)),
            "min": float(values.min()) if not values.empty else None,
            "max": float(values.max()) if not values.empty else None,
            "mean": float(values.mean()) if not values.empty else None,
            "median": float(values.median()) if not values.empty else None,
        }

    failure_counts = failures.get("error_type", pd.Series(dtype=object)).dropna().astype(str).value_counts()
    run_counts = scored.get("experiment_name", pd.Series(dtype=object)).dropna().astype(str).value_counts()

    return {
        "manifest": manifest,
        "overview": {
            "rows_scored": int(manifest.get("rows_scored", 0)),
            "rows_failed": int(manifest.get("rows_failed", 0)),
            "distinct_runs": int(scored.get("experiment_name", pd.Series(dtype=object)).nunique()) if not scored.empty else 0,
        },
        "spectral_angle_summary": _series_summary(scored, "spectral_angle"),
        "matched_ion_count_summary": _series_summary(scored, "matched_ion_count"),
        "explained_intensity_summary": _series_summary(scored, "explained_intensity_fraction"),
        "log2_xcorr_summary": _series_summary(scored, "log2_xcorr"),
        "log2_hyperscore_summary": _series_summary(scored, "log2_hyperscore"),
        "run_counts": {str(k): int(v) for k, v in run_counts.items()},
        "failure_type_counts": {str(k): int(v) for k, v in failure_counts.items()},
        "preview": {
            "scored_rows": scored.head(20).to_dict(orient="records") if not scored.empty else [],
            "failures": failures.head(20).to_dict(orient="records") if not failures.empty else [],
        },
    }


def inspect_scored_batch_status_tool(manifest_path: str) -> dict[str, Any]:
    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    status_path = (manifest.get("output_files") or {}).get("status_json")
    if status_path and Path(status_path).exists():
        return _read_batch_status(Path(status_path))
    return {
        "batch_id": manifest.get("batch_id"),
        "phase": "completed",
        "completed": int(manifest.get("candidate_rows_requested", 0)),
        "total": int(manifest.get("candidate_rows_requested", 0)),
        "remaining": 0,
        "succeeded": int(manifest.get("rows_scored", 0)),
        "failed": int(manifest.get("rows_failed", 0)),
        "progress_fraction": 1.0 if int(manifest.get("candidate_rows_requested", 0)) > 0 else None,
        "batch_seconds": manifest.get("batch_seconds"),
    }


def render_psm_table_tool(
    *,
    path: str,
    adapter: str = "instanovo",
    output_path: str,
    scored_path: str | None = None,
    sdrf_path: str | None = None,
    pride_project_accession: str | None = None,
    timeout_seconds: float = 20.0,
    limit: int | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    sdrf_summary = _resolve_sdrf_summary(sdrf_path=sdrf_path, pride_project_accession=pride_project_accession, timeout_seconds=timeout_seconds)
    inspection = summarize_candidate_inspection_frame(
        _load_candidate_frame(path, adapter, limit=limit, scored_path=scored_path, **adapter_options),
        sdrf_summary=sdrf_summary,
    )
    return render_psm_table_html(inspection, output_path)


def render_peptide_table_tool(
    *,
    path: str,
    adapter: str = "instanovo",
    output_path: str,
    scored_path: str | None = None,
    sdrf_path: str | None = None,
    pride_project_accession: str | None = None,
    timeout_seconds: float = 20.0,
    limit: int | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    sdrf_summary = _resolve_sdrf_summary(sdrf_path=sdrf_path, pride_project_accession=pride_project_accession, timeout_seconds=timeout_seconds)
    inspection = summarize_candidate_inspection_frame(
        _load_candidate_frame(path, adapter, limit=limit, scored_path=scored_path, **adapter_options),
        sdrf_summary=sdrf_summary,
    )
    return render_peptide_table_html(inspection, output_path)


def predict_fragments_tool(
    peptidoform: str,
    charge: int,
    model: str | None = None,
    backend: str | None = None,
    max_fragment_charge: int | None = None,
    config_path: str | None = None,
) -> dict[str, Any]:
    config = load_config(
        config_path,
        {
            "predictor_backend": backend,
            "predictor_model": model,
            "ms2pip_model": model,
            "max_fragment_charge": max_fragment_charge,
        },
    )
    return _predict_fragment_spectrum(config, peptidoform, charge).model_dump()


def describe_predictor_model_tool(
    model: str,
    *,
    backend: str = "koina",
    server_url: str = "https://koina.wilhelmlab.org",
    timeout_seconds: float = 30.0,
) -> dict[str, Any]:
    """Describe a predictor model and validate remote Koina input/output requirements."""
    if backend == "koina":
        metadata = describe_koina_model(
            model=model,
            server_url=server_url,
            timeout_seconds=timeout_seconds,
        )
        return {
            "backend": "koina",
            "model": model,
            "server_url": server_url.rstrip("/"),
            "metadata": metadata,
            "supported_by_proteomics_mcp": all(
                str(item.get("name"))
                in {"peptide_sequences", "precursor_charges", "collision_energies", "fragmentation_types", "instrument_types"}
                for item in metadata.get("inputs", [])
                if isinstance(item, dict)
            ),
        }
    if backend in {"ms2pip", "mock"}:
        return {
            "backend": backend,
            "model": model,
            "supported_by_proteomics_mcp": True,
            "metadata": {"local": True},
        }
    raise ValueError(f"Unsupported predictor backend: {backend}")


def retrieve_usi_spectrum_tool(usi: str, backend: str = "pride", timeout_seconds: float = 15.0) -> dict[str, Any]:
    return retrieve_spectrum_from_usi(usi, backend=backend, timeout_seconds=timeout_seconds).model_dump()


def load_sdrf_tool(path: str) -> dict[str, Any]:
    return load_sdrf_summary(path)


def retrieve_pride_sdrf_tool(project_accession: str, timeout_seconds: float = 20.0) -> dict[str, Any]:
    return retrieve_pride_sdrf(project_accession, timeout_seconds=timeout_seconds)


def mzml_to_feather_tool(mzml_path: str, output_path: str, ms_level: int = 2) -> dict[str, Any]:
    return mzml_to_feather(mzml_path, output_path, ms_level=ms_level)


def mgf_to_feather_tool(mgf_path: str, output_path: str) -> dict[str, Any]:
    return mgf_to_feather(mgf_path, output_path)


def bruker_d_to_feather_tool(d_path: str, output_path: str) -> dict[str, Any]:
    return bruker_d_to_feather(d_path, output_path)


def export_to_instanovo_tool(
    input_path: str,
    output_path: str,
    *,
    output_format: str = "ipc",
    sequence_column: str | None = None,
) -> dict[str, Any]:
    return export_table_to_instanovo(
        input_path,
        output_path,
        output_format=output_format,
        sequence_column=sequence_column,
    )


def annotate_observed_spectrum_tool(
    observed: dict[str, Any],
    *,
    peptidoform: str,
    ion_types: str = "by",
    tolerance: float = 0.05,
    tolerance_unit: str = "Da",
    max_isotope: int = 2,
    max_ion_charge: int | None = None,
    neutral_losses: bool = False,
) -> dict[str, Any]:
    observed_model = _model_from_argument(ObservedSpectrum, observed, "observed")
    return annotate_observed_ions(
        observed_model,
        peptidoform=peptidoform,
        tolerance=tolerance,
        tolerance_unit=tolerance_unit,
        ion_types=ion_types,
        max_isotope=max_isotope,
        max_ion_charge=max_ion_charge,
        neutral_losses=neutral_losses,
    )


def score_spectrum_tool(
    candidate: dict[str, Any],
    observed: dict[str, Any],
    predicted: dict[str, Any] | None = None,
    *,
    predictor_backend: str | None = None,
    predictor_model: str | None = None,
    ms2pip_model: str | None = None,
    max_fragment_charge: int | None = None,
    tolerance: float | None = None,
    tolerance_unit: str | None = None,
    annotation_ion_types: str | None = None,
    annotation_max_isotope: int | None = None,
    annotation_max_ion_charge: int | None = None,
    annotation_neutral_losses: bool | None = None,
    config_path: str | None = None,
) -> dict[str, Any]:
    config = load_config(
        config_path,
        {
            "predictor_backend": predictor_backend,
            "predictor_model": predictor_model,
            "ms2pip_model": ms2pip_model,
            "max_fragment_charge": max_fragment_charge,
            "tolerance": tolerance,
            "tolerance_unit": tolerance_unit,
            "annotation_ion_types": annotation_ion_types,
            "annotation_max_isotope": annotation_max_isotope,
            "annotation_max_ion_charge": annotation_max_ion_charge,
            "annotation_neutral_losses": annotation_neutral_losses,
        },
    )
    candidate_model = _model_from_argument(CandidatePSM, candidate, "candidate")
    observed_model = _model_from_argument(ObservedSpectrum, observed, "observed")
    if predicted is None:
        charge = candidate_model.precursor_charge or observed_model.precursor_charge or 2
        predicted_model = _predict_fragment_spectrum(config, candidate_model.peptidoform, charge)
    else:
        fragments = [_model_from_argument(PredictedFragment, fragment, "predicted.fragments[]") for fragment in predicted.get("fragments", [])]
        predicted_model = PredictedSpectrum(
            peptidoform=predicted["peptidoform"],
            precursor_charge=int(predicted["precursor_charge"]),
            fragments=fragments,
            predictor_name=predicted.get("predictor_name", "unknown"),
            predictor_version=predicted.get("predictor_version"),
            model=predicted.get("model"),
            metadata=predicted.get("metadata", {}),
        )
    result = score_candidate_spectrum(
        candidate=candidate_model,
        observed=observed_model,
        predicted=predicted_model,
        tolerance=config.tolerance,
        tolerance_unit=config.tolerance_unit,
        metric_profile=config.metric_profile,
        annotation_ion_types=config.annotation_ion_types,
        annotation_max_isotope=config.annotation_max_isotope,
        annotation_max_ion_charge=config.annotation_max_ion_charge,
        annotation_neutral_losses=config.annotation_neutral_losses,
    ).model_dump()
    result["metric_summary"] = summarize_score_metrics_tool(result)
    return result


def score_candidates_tool(
    *,
    candidates_path: str,
    spectra_path: str | None = None,
    output_path: str | None = None,
    adapter: str = "instanovo",
    predictor_backend: str | None = None,
    predictor_model: str | None = None,
    ms2pip_model: str | None = None,
    max_fragment_charge: int | None = None,
    tolerance: float | None = None,
    tolerance_unit: str | None = None,
    annotation_ion_types: str | None = None,
    annotation_max_isotope: int | None = None,
    annotation_max_ion_charge: int | None = None,
    annotation_neutral_losses: bool | None = None,
    limit: int | None = None,
    status_path: str | None = None,
    config_path: str | None = None,
    **adapter_options: Any,
) -> dict[str, Any]:
    started = time.time()
    started_at = _utc_timestamp()
    batch_id = f"score_candidates_{uuid.uuid4().hex[:8]}"
    config = load_config(
        config_path,
        {
            "predictor_backend": predictor_backend,
            "predictor_model": predictor_model,
            "ms2pip_model": ms2pip_model,
            "max_fragment_charge": max_fragment_charge,
            "tolerance": tolerance,
            "tolerance_unit": tolerance_unit,
            "annotation_ion_types": annotation_ion_types,
            "annotation_max_isotope": annotation_max_isotope,
            "annotation_max_ion_charge": annotation_max_ion_charge,
            "annotation_neutral_losses": annotation_neutral_losses,
        },
    )
    candidates = _load_candidate_models(candidates_path, adapter, limit=limit, **adapter_options)
    rows: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    resolved_status_path = Path(status_path) if status_path else None
    update_interval = _status_update_interval(int(len(candidates)))
    if resolved_status_path is not None:
        _write_batch_status(
            resolved_status_path,
            batch_id=batch_id,
            phase="scoring",
            completed=0,
            total=int(len(candidates)),
            succeeded=0,
            failed=0,
            started=started,
            extra={
                "started_at": started_at,
                "candidates_path": candidates_path,
                "spectra_path": spectra_path,
                "adapter": adapter,
                "predictor_backend": config.predictor_backend,
                "requested_model": config.predictor_model,
            },
        )
    for index, candidate in enumerate(candidates):
        try:
            if candidate.usi:
                observed = retrieve_spectrum_from_usi(candidate.usi, backend=config.usi_backend)
            elif spectra_path:
                observed = load_spectrum_from_table(
                    spectra_path,
                    scan_number=candidate.scan_number,
                    spectrum_id=candidate.spectrum_id,
                    experiment_name=candidate.experiment_name,
                )
            else:
                raise ValueError("Candidate has no USI and no spectra_path was supplied.")
            charge = candidate.precursor_charge or observed.precursor_charge or 2
            predicted = _predict_fragment_spectrum(config, candidate.peptidoform, charge)
            result = score_candidate_spectrum(
                candidate=candidate,
                observed=observed,
                predicted=predicted,
                tolerance=config.tolerance,
                tolerance_unit=config.tolerance_unit,
                metric_profile=config.metric_profile,
                annotation_ion_types=config.annotation_ion_types,
                annotation_max_isotope=config.annotation_max_isotope,
                annotation_max_ion_charge=config.annotation_max_ion_charge,
                annotation_neutral_losses=config.annotation_neutral_losses,
            )
            acquisition_summary = result.acquisition_provenance.get("summary", {})
            rows.append(
                {
                    "candidate_id": candidate.candidate_id,
                    "source_engine": candidate.source_engine,
                    "source_engine_version": candidate.source_engine_version,
                    "experiment_name": candidate.experiment_name,
                    "scan_number": candidate.scan_number,
                    "spectrum_id": candidate.spectrum_id,
                    "usi": candidate.usi,
                    "peptidoform": candidate.peptidoform,
                    "precursor_mz": candidate.precursor_mz,
                    "precursor_charge": candidate.precursor_charge,
                    "rank": candidate.rank,
                    "confidence": candidate.confidence,
                    "raw_score": candidate.raw_score,
                    "raw_score_name": candidate.raw_score_name,
                    "next_beam": candidate.next_beam,
                    "next_beam_margin": candidate.next_beam_margin,
                    "token_log_probability_mean": candidate.token_log_probability_mean,
                    "token_log_probability_min": candidate.token_log_probability_min,
                    "token_confidence_mean": candidate.token_confidence_mean,
                    "token_confidence_min": candidate.token_confidence_min,
                    "token_confidence_min_position": candidate.token_confidence_min_position,
                    "token_confidence_min_token": candidate.token_confidence_min_token,
                    "spectral_angle": result.spectral_angle,
                    "cosine_similarity": result.cosine_similarity,
                    "pearson_correlation": result.pearson_correlation,
                    "xcorr": result.xcorr,
                    "log10_xcorr": result.log10_xcorr,
                    "log2_xcorr": result.log2_xcorr,
                    "hyperscore": result.hyperscore,
                    "log2_hyperscore": result.log2_hyperscore,
                    "matched_ion_count": result.matched_ion_count,
                    "total_predicted_ions": result.total_predicted_ions,
                    "matched_ion_fraction": result.matched_ion_fraction,
                    "explained_intensity_fraction": result.explained_intensity_fraction,
                    "precursor_signal_fraction": result.precursor_signal_fraction,
                    "fragment_mass_error_mean_ppm": result.fragment_mass_error_mean_ppm,
                    "fragment_mass_error_abs_mean_ppm": result.fragment_mass_error_abs_mean_ppm,
                    "fragment_mass_error_mean_da": result.fragment_mass_error_mean_da,
                    "fragment_mass_error_abs_mean_da": result.fragment_mass_error_abs_mean_da,
                    "precursor_mass_error_ppm": result.precursor_mass_error_ppm,
                    "precursor_mass_error_da": result.precursor_mass_error_da,
                    "fragmentation_method": acquisition_summary.get("fragmentation_method"),
                    "isolation_target_mz": acquisition_summary.get("isolation_target_mz"),
                    "retention_time_seconds": acquisition_summary.get("retention_time_seconds"),
                    "retention_time_minutes": acquisition_summary.get("retention_time_minutes"),
                    "collision_energy": acquisition_summary.get("collision_energy"),
                    "faims_compensation_voltage": acquisition_summary.get("faims_compensation_voltage"),
                    "injection_time": acquisition_summary.get("injection_time"),
                    "scan_header": acquisition_summary.get("scan_header"),
                    "scan_resolution": acquisition_summary.get("scan_resolution"),
                    "instrument_model": acquisition_summary.get("instrument_model"),
                    "detector_type": acquisition_summary.get("detector_type"),
                    "annotated_peak_count": result.annotation_summary.get("annotated_peak_count"),
                    "ion_type_counts": result.annotation_summary.get("ion_type_counts"),
                    "annotated_intensity_fraction": result.annotation_summary.get("annotated_intensity_fraction"),
                    "metric_profile": result.metric_profile,
                }
            )
        except Exception as exc:
            failures.append({"candidate_id": candidate.candidate_id, "error": str(exc)})
        if resolved_status_path is not None:
            completed = int(index + 1)
            if completed == int(len(candidates)) or completed % update_interval == 0:
                _write_batch_status(
                    resolved_status_path,
                    batch_id=batch_id,
                    phase="scoring",
                    completed=completed,
                    total=int(len(candidates)),
                    succeeded=int(len(rows)),
                    failed=int(len(failures)),
                    started=started,
                    extra={
                        "started_at": started_at,
                        "candidates_path": candidates_path,
                        "spectra_path": spectra_path,
                        "adapter": adapter,
                        "predictor_backend": config.predictor_backend,
                        "requested_model": config.predictor_model,
                        "latest_candidate_id": candidate.candidate_id,
                    },
                )
    if output_path:
        pd.DataFrame(rows).to_csv(output_path, index=False)
    if resolved_status_path is not None:
        _write_batch_status(
            resolved_status_path,
            batch_id=batch_id,
            phase="completed",
            completed=int(len(candidates)),
            total=int(len(candidates)),
            succeeded=int(len(rows)),
            failed=int(len(failures)),
            started=started,
            extra={
                "started_at": started_at,
                "candidates_path": candidates_path,
                "spectra_path": spectra_path,
                "adapter": adapter,
                "predictor_backend": config.predictor_backend,
                "requested_model": config.predictor_model,
                "output_path": output_path,
            },
        )
    return {
        "candidates_read": len(candidates),
        "rows_scored": len(rows),
        "rows_failed": len(failures),
        "output_path": output_path,
        "status_path": str(resolved_status_path) if resolved_status_path is not None else None,
        "rows_preview": rows[:20],
        "metric_groups": {
            "match_quality": ["spectral_angle", "cosine_similarity", "pearson_correlation", "xcorr", "log10_xcorr", "log2_xcorr", "hyperscore", "log2_hyperscore"],
            "predicted_match_evidence": ["matched_ion_count", "total_predicted_ions", "matched_ion_fraction", "explained_intensity_fraction"],
            "observed_annotation_evidence": ["precursor_signal_fraction", "annotated_peak_count", "ion_type_counts", "annotated_intensity_fraction"],
            "ion_evidence": ["matched_ion_count", "total_predicted_ions", "matched_ion_fraction", "explained_intensity_fraction", "precursor_signal_fraction", "annotated_peak_count", "ion_type_counts", "annotated_intensity_fraction"],
            "mass_error": [
                "fragment_mass_error_mean_ppm",
                "fragment_mass_error_abs_mean_ppm",
                "fragment_mass_error_mean_da",
                "fragment_mass_error_abs_mean_da",
                "precursor_mass_error_ppm",
                "precursor_mass_error_da",
            ],
            "candidate": [
                "confidence",
                "raw_score",
                "raw_score_name",
                "next_beam",
                "next_beam_margin",
                "token_log_probability_mean",
                "token_log_probability_min",
                "token_confidence_mean",
                "token_confidence_min",
            ],
            "acquisition": ["fragmentation_method", "isolation_target_mz", "retention_time_minutes", "collision_energy", "faims_compensation_voltage", "injection_time", "scan_header", "scan_resolution", "instrument_model", "detector_type"],
        },
        "failures": failures[:20],
    }


def _resolve_sdrf_summary(*, sdrf_path: str | None = None, pride_project_accession: str | None = None, timeout_seconds: float = 20.0) -> dict[str, Any] | None:
    if sdrf_path:
        return load_sdrf_summary(sdrf_path)
    if pride_project_accession:
        return retrieve_pride_sdrf(pride_project_accession, timeout_seconds=timeout_seconds)
    return None


def summarize_dataset_tool(
    *,
    scored_path: str,
    output_prefix: str | None = None,
    sdrf_path: str | None = None,
    pride_project_accession: str | None = None,
    timeout_seconds: float = 20.0,
) -> dict[str, Any]:
    sdrf_summary = _resolve_sdrf_summary(sdrf_path=sdrf_path, pride_project_accession=pride_project_accession, timeout_seconds=timeout_seconds)
    aggregation = summarize_dataset_path(scored_path, sdrf_summary=sdrf_summary)
    output_files = write_dataset_outputs(aggregation, output_prefix) if output_prefix else {}
    return {
        "dataset_metadata": aggregation.dataset_summary["dataset_metadata"],
        "identification_summary": aggregation.dataset_summary["identification_summary"],
        "sequence_summary": aggregation.dataset_summary["sequence_summary"],
        "peptidoform_summary": aggregation.dataset_summary["peptidoform_summary"],
        "run_summary": aggregation.dataset_summary["run_summary"],
        "charge_summary": aggregation.dataset_summary["charge_summary"],
        "confidence_summary": aggregation.dataset_summary["confidence_summary"],
        "spectrum_evidence_summary": aggregation.dataset_summary["spectrum_evidence_summary"],
        "modification_summary": aggregation.dataset_summary["modification_summary"],
        "sdrf_summary": aggregation.dataset_summary["sdrf_summary"],
        "output_files": output_files,
        "preview": {
            "peptide_summary": _preview_records(aggregation.peptide_summary),
            "peptidoform_summary": _preview_records(aggregation.peptidoform_summary),
            "run_summary": _preview_records(aggregation.run_summary),
            "sample_summary": _preview_records(aggregation.sample_summary) if aggregation.sample_summary is not None else [],
        },
    }


def render_dataset_report_tool(
    *,
    scored_path: str,
    output_path: str | None = None,
    sdrf_path: str | None = None,
    pride_project_accession: str | None = None,
    timeout_seconds: float = 20.0,
) -> dict[str, Any]:
    sdrf_summary = _resolve_sdrf_summary(sdrf_path=sdrf_path, pride_project_accession=pride_project_accession, timeout_seconds=timeout_seconds)
    aggregation = summarize_dataset_path(scored_path, sdrf_summary=sdrf_summary)
    html = render_dataset_report_html(aggregation, output_path)
    return {
        "format": "html",
        "output_path": output_path,
        "bytes": len(html),
        "dataset_metadata": aggregation.dataset_summary["dataset_metadata"],
    }


def list_dataset_groups_tool(
    *,
    scored_path: str,
    sdrf_path: str | None = None,
    pride_project_accession: str | None = None,
    timeout_seconds: float = 20.0,
) -> dict[str, Any]:
    sdrf_summary = _resolve_sdrf_summary(sdrf_path=sdrf_path, pride_project_accession=pride_project_accession, timeout_seconds=timeout_seconds)
    aggregation = summarize_dataset_path(scored_path, sdrf_summary=sdrf_summary)
    return dataset_groups(aggregation)


def get_peptide_summary_tool(
    *,
    scored_path: str,
    peptide: str,
    sdrf_path: str | None = None,
    pride_project_accession: str | None = None,
    timeout_seconds: float = 20.0,
) -> dict[str, Any]:
    sdrf_summary = _resolve_sdrf_summary(sdrf_path=sdrf_path, pride_project_accession=pride_project_accession, timeout_seconds=timeout_seconds)
    aggregation = summarize_dataset_path(scored_path, sdrf_summary=sdrf_summary)
    summary = get_peptide_summary_row(aggregation, peptide)
    if summary is None:
        raise LookupError(f"Peptide '{peptide}' not found in scored dataset.")
    return summary


def get_peptidoform_summary_tool(
    *,
    scored_path: str,
    peptidoform: str,
    sdrf_path: str | None = None,
    pride_project_accession: str | None = None,
    timeout_seconds: float = 20.0,
) -> dict[str, Any]:
    sdrf_summary = _resolve_sdrf_summary(sdrf_path=sdrf_path, pride_project_accession=pride_project_accession, timeout_seconds=timeout_seconds)
    aggregation = summarize_dataset_path(scored_path, sdrf_summary=sdrf_summary)
    summary = get_peptidoform_summary_row(aggregation, peptidoform)
    if summary is None:
        raise LookupError(f"Peptidoform '{peptidoform}' not found in scored dataset.")
    return summary


def render_spectrum_comparison_tool(payload: dict[str, Any], output_path: str | None = None) -> dict[str, Any]:
    html = render_spectrum_comparison_html(payload, output_path)
    return {"format": "html", "output_path": output_path, "bytes": len(html), "download_image": "plotly-modebar"}


def render_observed_spectrum_tool(payload_or_observed: dict[str, Any], output_path: str | None = None) -> dict[str, Any]:
    html = render_observed_spectrum_html(payload_or_observed, output_path)
    return {"format": "html", "output_path": output_path, "bytes": len(html), "view": "observed_spectrum", "download_image": "plotly-modebar"}


def render_annotated_observed_spectrum_tool(
    payload_or_observed: dict[str, Any],
    output_path: str | None = None,
    *,
    peptidoform: str | None = None,
    ion_types: str = "by",
    tolerance: float = 0.05,
    tolerance_unit: str = "Da",
    max_isotope: int = 2,
    max_ion_charge: int | None = None,
    neutral_losses: bool = False,
) -> dict[str, Any]:
    html = render_annotated_observed_spectrum_html(
        payload_or_observed,
        output_path,
        peptidoform=peptidoform,
        ion_types=ion_types,
        tolerance=tolerance,
        tolerance_unit=tolerance_unit,
        max_isotope=max_isotope,
        max_ion_charge=max_ion_charge,
        neutral_losses=neutral_losses,
    )
    return {"format": "html", "output_path": output_path, "bytes": len(html), "view": "annotated_observed_spectrum", "download_image": "plotly-modebar"}


def render_static_spectrum_tool(
    payload_or_observed: dict[str, Any],
    output_path: str,
    *,
    view: str = "mirror",
    peptidoform: str | None = None,
    ion_types: str = "by",
    tolerance: float = 0.05,
    tolerance_unit: str = "Da",
    max_isotope: int = 2,
    max_ion_charge: int | None = None,
    neutral_losses: bool = False,
    include_fragment_sequences: bool = False,
    label_top_n: int | None = None,
    title: str | None = None,
) -> dict[str, Any]:
    """Write a publication-ready static spectrum plot as SVG, PDF, or PNG."""
    rendered = render_static_spectrum(
        payload_or_observed,
        output_path,
        view=view,
        peptidoform=peptidoform,
        ion_types=ion_types,
        tolerance=tolerance,
        tolerance_unit=tolerance_unit,
        max_isotope=max_isotope,
        max_ion_charge=max_ion_charge,
        neutral_losses=neutral_losses,
        include_fragment_sequences=include_fragment_sequences,
        label_top_n=label_top_n,
        title=title,
    )
    return {"format": rendered.suffix.lstrip("."), "output_path": str(rendered), "view": view, "renderer": "matplotlib_seaborn"}


def _safe_output_stem(value: str | None, fallback: str) -> str:
    text = (value or "").strip()
    if not text:
        text = fallback
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text).strip("._-")
    return text or fallback


def _render_batch_item(
    *,
    view: str,
    item: dict[str, Any],
    output_dir: Path,
    item_index: int,
) -> dict[str, Any]:
    payload = item.get("payload")
    if not isinstance(payload, dict):
        raise ValueError("Each batch item must include a 'payload' object.")
    stem = _safe_output_stem(item.get("output_stem") or item.get("name") or item.get("candidate_id"), f"spectrum_{item_index + 1}")
    view_dir = output_dir / view
    view_dir.mkdir(parents=True, exist_ok=True)
    output_path = view_dir / f"{stem}.html"
    if view == "comparison":
        result = render_spectrum_comparison_tool(payload, str(output_path))
    elif view == "observed":
        result = render_observed_spectrum_tool(payload, str(output_path))
    elif view == "annotated":
        annotated_options = item.get("annotated_options") or {}
        if not isinstance(annotated_options, dict):
            raise ValueError("annotated_options must be a JSON object when provided.")
        result = render_annotated_observed_spectrum_tool(payload, str(output_path), **annotated_options)
    else:
        raise ValueError(f"Unsupported render view: {view}")
    return {
        "view": view,
        "item_index": item_index,
        "name": item.get("name"),
        "output_stem": stem,
        "output_path": str(output_path),
        "bytes": result.get("bytes"),
        "format": result.get("format"),
    }


def _load_render_batch_spec_from_path(input_path: str) -> tuple[list[dict[str, Any]], list[str] | None]:
    """Read a batch render spec: JSON (a list of items, or {"items": [...], "views": [...]}), or JSON Lines
    (.jsonl/.ndjson, one item per line); either may be gzip-compressed (.gz)."""
    path = Path(input_path)
    suffixes = [suffix.lower() for suffix in path.suffixes]
    compressed = bool(suffixes) and suffixes[-1] == ".gz"
    with (gzip.open(path, "rt", encoding="utf-8") if compressed else path.open(encoding="utf-8")) as handle:
        text = handle.read()
    if (suffixes[-2:-1] if compressed else suffixes[-1:]) in ([".jsonl"], [".ndjson"]):
        lines = [(number, line) for number, line in enumerate(text.splitlines(), start=1) if line.strip()]
        items = []
        for number, line in lines:
            try:
                items.append(json.loads(line))
            except json.JSONDecodeError as exc:
                raise ValueError(f"{input_path}, line {number}: not valid JSON ({exc.msg}).") from exc
        return items, None
    payload = json.loads(text)
    if isinstance(payload, dict) and "items" in payload:
        items = payload["items"]
        views = payload.get("views")
    else:
        items = payload
        views = None
    return items, views


def _validate_render_batch_item(item: Any, *, item_index: int) -> dict[str, Any] | None:
    if not isinstance(item, dict):
        return {"item_index": item_index, "error": "Each batch item must be an object.", "error_type": "ValueError"}
    if not isinstance(item.get("payload"), dict):
        return {
            "item_index": item_index,
            "name": item.get("name"),
            "error": "Each batch item must include a 'payload' object.",
            "error_type": "ValueError",
        }
    annotated_options = item.get("annotated_options")
    if annotated_options is not None and not isinstance(annotated_options, dict):
        return {
            "item_index": item_index,
            "name": item.get("name"),
            "error": "annotated_options must be a JSON object when provided.",
            "error_type": "ValueError",
        }
    return None


def render_spectra_batch_tool(
    items: list[dict[str, Any]] | None = None,
    output_dir: str | None = None,
    *,
    input_path: str | None = None,
    views: list[str] | None = None,
    max_workers: int | None = None,
    validate_only: bool = False,
    status_path: str | None = None,
) -> dict[str, Any]:
    if (items is None) == (input_path is None):
        raise ValueError("Exactly one of 'items' or 'input_path' must be provided.")
    if input_path is not None:
        loaded_items, file_views = _load_render_batch_spec_from_path(input_path)
        items = loaded_items
        if views is None:
            views = file_views
    if not isinstance(items, list) or not items:
        raise ValueError("items must be a non-empty list of render payload objects.")
    normalized_views = views or ["comparison", "annotated", "observed"]
    supported_views = {"comparison", "annotated", "observed"}
    invalid_views = [view for view in normalized_views if view not in supported_views]
    if invalid_views:
        raise ValueError(f"Unsupported render views: {', '.join(invalid_views)}")
    validation_failures = [failure for item_index, item in enumerate(items) if (failure := _validate_render_batch_item(item, item_index=item_index)) is not None]
    if validate_only:
        return {
            "format": "validation",
            "input_path": input_path,
            "views": normalized_views,
            "items_requested": len(items),
            "valid_items": len(items) - len(validation_failures),
            "invalid_items": len(validation_failures),
            "failures": validation_failures,
        }
    if validation_failures:
        return {
            "format": "html",
            "input_path": input_path,
            "output_dir": output_dir,
            "views": normalized_views,
            "items_requested": len(items),
            "items_rendered": 0,
            "items_per_view": {view: 0 for view in normalized_views},
            "max_workers": max(1, min(max_workers or (os.cpu_count() or 1), len(items))),
            "results": {view: [] for view in normalized_views},
            "failures": validation_failures,
        }
    if not output_dir:
        raise ValueError("output_dir is required unless validate_only=true.")
    batch_dir = Path(output_dir)
    batch_dir.mkdir(parents=True, exist_ok=True)
    worker_count = max(1, min(max_workers or (os.cpu_count() or 1), len(items)))
    started = time.time()
    started_at = _utc_timestamp()
    batch_id = f"render_{uuid.uuid4().hex[:8]}"
    resolved_status_path = Path(status_path) if status_path else batch_dir / "render_batch_status.json"
    results_by_view: dict[str, list[dict[str, Any]]] = {view: [] for view in normalized_views}
    failures: list[dict[str, Any]] = []
    total_tasks = int(len(items) * len(normalized_views))
    completed_tasks = 0
    _write_batch_status(
        resolved_status_path,
        batch_id=batch_id,
        phase="rendering",
        completed=0,
        total=total_tasks,
        succeeded=0,
        failed=0,
        started=started,
        extra={
            "started_at": started_at,
            "input_path": input_path,
            "output_dir": str(batch_dir),
            "views": normalized_views,
            "items_requested": len(items),
            "max_workers": worker_count,
        },
    )
    for view in normalized_views:
        with ThreadPoolExecutor(max_workers=worker_count, thread_name_prefix=f"render_{view}") as executor:
            future_map = {
                executor.submit(_render_batch_item, view=view, item=item, output_dir=batch_dir, item_index=item_index): (item_index, item)
                for item_index, item in enumerate(items)
            }
            for future in as_completed(future_map):
                item_index, item = future_map[future]
                try:
                    results_by_view[view].append(future.result())
                except Exception as exc:
                    failures.append(
                        {
                            "view": view,
                            "item_index": item_index,
                            "name": item.get("name"),
                            "error": str(exc),
                            "error_type": type(exc).__name__,
                        }
                    )
                completed_tasks += 1
                _write_batch_status(
                    resolved_status_path,
                    batch_id=batch_id,
                    phase="rendering",
                    completed=completed_tasks,
                    total=total_tasks,
                    succeeded=sum(len(rows) for rows in results_by_view.values()),
                    failed=int(len(failures)),
                    started=started,
                    extra={
                        "started_at": started_at,
                        "input_path": input_path,
                        "output_dir": str(batch_dir),
                        "views": normalized_views,
                        "items_requested": len(items),
                        "max_workers": worker_count,
                        "latest_view": view,
                        "latest_item_index": item_index,
                        "latest_item_name": item.get("name"),
                    },
                )
        results_by_view[view].sort(key=lambda row: row["item_index"])
    _write_batch_status(
        resolved_status_path,
        batch_id=batch_id,
        phase="completed",
        completed=total_tasks,
        total=total_tasks,
        succeeded=sum(len(rows) for rows in results_by_view.values()),
        failed=int(len(failures)),
        started=started,
        extra={
            "started_at": started_at,
            "input_path": input_path,
            "output_dir": str(batch_dir),
            "views": normalized_views,
            "items_requested": len(items),
            "max_workers": worker_count,
        },
    )
    return {
        "format": "html",
        "input_path": input_path,
        "output_dir": str(batch_dir),
        "status_path": str(resolved_status_path),
        "views": normalized_views,
        "items_requested": len(items),
        "items_rendered": sum(len(rows) for rows in results_by_view.values()),
        "items_per_view": {view: len(rows) for view, rows in results_by_view.items()},
        "max_workers": worker_count,
        "results": results_by_view,
        "failures": failures,
    }
