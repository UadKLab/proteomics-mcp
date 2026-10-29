from __future__ import annotations

import functools
import logging
import os
import urllib.error
from pathlib import Path
from typing import Any

from proteomics_mcp import __version__
from proteomics_mcp.resources.docs import SCHEMAS

# --- Logging --------------------------------------------------------------
# Honor PROTEOMICS_MCP_LOG_LEVEL (DEBUG / INFO / WARNING / ERROR). Default INFO.
# Configured at import time so the MCP runtime sees logs from the first request.
_log_level = os.environ.get("PROTEOMICS_MCP_LOG_LEVEL", "INFO").upper()
logging.basicConfig(
    level=getattr(logging, _log_level, logging.INFO),
    format="%(asctime)s %(levelname)s %(name)s: %(message)s",
)
logger = logging.getLogger("proteomics_mcp")

# --- Output directory fallback -------------------------------------------
# If PROTEOMICS_MCP_OUTPUT_DIR is set, any relative output_path / output_dir /
# output_prefix passed to a tool is resolved against it. Absolute paths and
# unset values are left untouched, preserving existing behavior.
_OUTPUT_BASE = os.environ.get("PROTEOMICS_MCP_OUTPUT_DIR")
if _OUTPUT_BASE:
    Path(_OUTPUT_BASE).expanduser().mkdir(parents=True, exist_ok=True)
    logger.info("Default output base: %s", _OUTPUT_BASE)


def _resolve(path: str | None) -> str | None:
    """Resolve a relative output path against PROTEOMICS_MCP_OUTPUT_DIR if set.

    No-op when path is None, already absolute, or the env var is unset.
    """
    if path is None or not _OUTPUT_BASE:
        return path
    p = Path(path).expanduser()
    if p.is_absolute():
        return str(p)
    return str(Path(_OUTPUT_BASE).expanduser() / p)
from proteomics_mcp.tools.workflows import (
    annotate_observed_spectrum_tool,
    bruker_d_to_feather_tool,
    export_to_instanovo_tool,
    get_peptide_summary_tool,
    get_peptidoform_summary_tool,
    inspect_peptide_table_tool,
    inspect_psm_table_tool,
    list_distribution_fields_tool,
    load_candidates_tool,
    mzml_to_feather_tool,
    thermo_raw_to_feather_tool,
    thermo_raw_to_mzml_tool,
    mgf_to_feather_tool,
    predict_fragments_tool,
    describe_predictor_model_tool,
    load_sdrf_tool,
    list_dataset_groups_tool,
    summarize_cv_by_category_tool,
    summarize_dashboard_tool,
    summarize_file_coverage_tool,
    summarize_missing_values_tool,
    summarize_modifications_tool,
    summarize_prediction_accuracy_tool,
    summarize_overlap_tool,
    summarize_quant_quality_tool,
    summarize_categorical_field_tool,
    bin_numeric_field_tool,
    bin_2d_numeric_fields_tool,
    build_id_dashboard_tool,
    build_prediction_batch_tool,
    inspect_prediction_batch_tool,
    inspect_prediction_batch_status_tool,
    score_from_prediction_batch_tool,
    inspect_scored_batch_tool,
    inspect_scored_batch_status_tool,
    render_peptide_table_tool,
    render_psm_table_tool,
    render_dashboard_html_tool,
    render_prediction_accuracy_html_tool,
    render_summary_html_tool,
    render_annotated_observed_spectrum_tool,
    render_static_spectrum_tool,
    render_dataset_report_tool,
    render_observed_spectrum_tool,
    render_spectra_batch_tool,
    render_spectrum_comparison_tool,
    retrieve_pride_sdrf_tool,
    retrieve_usi_spectrum_tool,
    score_candidates_tool,
    score_spectrum_tool,
    summarize_numeric_field_tool,
    summarize_dataset_tool,
    summarize_score_metrics_tool,
)

from proteomics_mcp.predictors.ms2pip import MissingOpenMPError

try:
    from mcp.server import MCPServer  # type: ignore
    from mcp.server.mcpserver.exceptions import ToolError
    from mcp.types import ToolAnnotations
except ImportError:  # pragma: no cover - exercised when optional dependency is absent
    MCPServer = None  # type: ignore
    ToolAnnotations = None  # type: ignore


if MCPServer is not None:
    _READ_ONLY_TOOLS = {
        "load_candidates",
        "inspect_psm_table",
        "inspect_peptide_table",
        "list_distribution_fields",
        "summarize_numeric_field",
        "summarize_categorical_field",
        "summarize_missing_values",
        "summarize_cv_by_category",
        "bin_numeric_field",
        "bin_2d_numeric_fields",
        "summarize_dashboard",
        "summarize_overlap",
        "summarize_modifications",
        "summarize_prediction_accuracy",
        "summarize_quant_quality",
        "summarize_file_coverage",
        "inspect_prediction_batch",
        "inspect_prediction_batch_status",
        "inspect_scored_batch",
        "inspect_scored_batch_status",
        "load_sdrf",
        "annotate_observed_spectrum",
        "summarize_score_metrics",
        "list_dataset_groups",
        "get_peptide_summary",
        "get_peptidoform_summary",
    }

    # The MCP SDK hides the message of any exception other than ToolError from the
    # client ("Error executing tool <name>"). Re-raise the errors the tools raise
    # deliberately for input the user or agent can correct (bad arguments or tables,
    # missing files, unknown fields or peptides, unreachable services) as ToolError,
    # so the message reaches the agent. Their messages name only the user's own inputs.
    # Anything else (RuntimeError, TypeError, ...) is a defect and stays masked and logged.
    _ACTIONABLE_ERRORS = (
        MissingOpenMPError,
        ValueError,  # includes pydantic ValidationError
        LookupError,  # KeyError, and LookupError for unknown peptides, spectra and SDRF files
        FileNotFoundError,
        IsADirectoryError,
        NotADirectoryError,
        PermissionError,
        ConnectionError,
        TimeoutError,
        urllib.error.URLError,
    )

    def _actionable_message(exc: BaseException) -> str:
        if isinstance(exc, MissingOpenMPError):
            return str(exc)
        # str(KeyError("x")) is "'x'"; show the message itself
        detail = exc.args[0] if isinstance(exc, KeyError) and len(exc.args) == 1 else exc
        return f"{type(exc).__name__}: {detail}"

    def _surface_actionable_errors(func: Any) -> Any:
        @functools.wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            try:
                return func(*args, **kwargs)
            except _ACTIONABLE_ERRORS as exc:
                raise ToolError(_actionable_message(exc)) from exc

        return wrapper

    class AnnotatedMCPServer(MCPServer):
        def tool(self, *args: Any, **kwargs: Any) -> Any:
            def register(func: Any) -> Any:
                read_only = func.__name__ in _READ_ONLY_TOOLS
                annotations = ToolAnnotations(
                    title=func.__name__.replace("_", " ").title(),
                    readOnlyHint=read_only,
                    destructiveHint=not read_only,
                    openWorldHint=func.__name__ in {
                        "predict_fragments",
                        "describe_predictor_model",
                        "retrieve_usi_spectrum",
                        "retrieve_pride_sdrf",
                    },
                )
                return super(AnnotatedMCPServer, self).tool(
                    *args,
                    title=annotations.title,
                    annotations=annotations,
                    **kwargs,
                )(_surface_actionable_errors(func))

            return register

    mcp = AnnotatedMCPServer("proteomics-mcp", version=__version__)

    @mcp.tool()
    def load_candidates(
        path: str,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        fixed_modifications: str | None = None,
        casanovo_version: str | None = None,
        limit: int | None = None,
    ) -> list[dict[str, Any]]:
        """Normalize an upstream candidate or identification table into the canonical PSM schema.

        For adapter="generic", column_mapping_path is a JSON file mapping the table's columns,
        e.g. {"peptide_column": "Sequence", "scan_column": "Scan"}; see method://generic_mapping_v1.
        For adapter="maxquant", fixed_modifications lists MaxQuant's fixed modifications, e.g. "Carbamidomethyl (C)"
        (";"-separated; "" for none). Without it they are read from the parameters.txt or mqpar.xml next to msms.txt;
        if neither exists, each candidate carries a warning.
        For adapter="casanovo", casanovo_version (e.g. "4.2.1") says what a peptide_score / search_engine_score[1]
        column means: a mean of amino-acid scores before 5.0, a sequence probability from 5.0. Without it that score
        gets no confidence. A log_probs column is read as a natural-log probability whatever the version.
        """
        return load_candidates_tool(path, adapter, column_mapping_path=column_mapping_path,
                                    fixed_modifications=fixed_modifications, casanovo_version=casanovo_version, limit=limit)

    @mcp.tool()
    def inspect_psm_table(
        path: str,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        scored_path: str | None = None,
        sdrf_path: str | None = None,
        pride_project_accession: str | None = None,
        timeout_seconds: float = 20.0,
        limit: int | None = None,
        max_rows: int | None = None,
    ) -> dict[str, Any]:
        """Inspect a normalized PSM table with dataset, confidence, run, and preview summaries.

        The summaries always cover the whole table; limit sets the number of preview rows (default 10). max_rows reads
        only the first rows of a very large table, and the result then reports it (row_limit, warnings).
        """
        return inspect_psm_table_tool(
            path=path,
            adapter=adapter,
            column_mapping_path=column_mapping_path,
            scored_path=scored_path,
            sdrf_path=sdrf_path,
            pride_project_accession=pride_project_accession,
            timeout_seconds=timeout_seconds,
            limit=limit,
            max_rows=max_rows,
        )

    @mcp.tool()
    def inspect_peptide_table(
        path: str,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        peptide: str | None = None,
        scored_path: str | None = None,
        sdrf_path: str | None = None,
        pride_project_accession: str | None = None,
        timeout_seconds: float = 20.0,
        limit: int | None = None,
        max_rows: int | None = None,
    ) -> dict[str, Any]:
        """Inspect the collapsed peptide table or return one peptide summary row by stripped sequence.

        The summaries always cover the whole table; limit sets the number of preview rows (default 10). max_rows reads
        only the first rows of a very large table, and the result then reports it (row_limit, warnings).
        """
        return inspect_peptide_table_tool(
            path=path,
            adapter=adapter,
            column_mapping_path=column_mapping_path,
            peptide=peptide,
            scored_path=scored_path,
            sdrf_path=sdrf_path,
            pride_project_accession=pride_project_accession,
            timeout_seconds=timeout_seconds,
            limit=limit,
            max_rows=max_rows,
        )

    @mcp.tool()
    def list_distribution_fields(
        path: str,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        scored_path: str | None = None,
        limit: int | None = None,
    ) -> dict[str, Any]:
        """List numeric and categorical fields available for summary and binning on a normalized table."""
        return list_distribution_fields_tool(path=path, adapter=adapter, column_mapping_path=column_mapping_path, scored_path=scored_path, limit=limit)

    @mcp.tool()
    def summarize_numeric_field(
        path: str,
        field: str,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        scored_path: str | None = None,
        group_by: str | None = None,
        percentiles: list[float] | None = None,
        weight_field: str | None = None,
        limit: int | None = None,
    ) -> dict[str, Any]:
        """Summarize one numeric field with optional grouping, weighting, and custom percentiles."""
        return summarize_numeric_field_tool(
            path=path,
            field=field,
            adapter=adapter,
            column_mapping_path=column_mapping_path,
            scored_path=scored_path,
            group_by=group_by,
            percentiles=percentiles,
            weight_field=weight_field,
            limit=limit,
        )

    @mcp.tool()
    def summarize_categorical_field(
        path: str,
        field: str,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        scored_path: str | None = None,
        group_by: str | None = None,
        top_n: int | None = None,
        limit: int | None = None,
    ) -> dict[str, Any]:
        """Count values for one categorical field with optional grouping and top-N truncation."""
        return summarize_categorical_field_tool(
            path=path,
            field=field,
            adapter=adapter,
            column_mapping_path=column_mapping_path,
            scored_path=scored_path,
            group_by=group_by,
            top_n=top_n,
            limit=limit,
        )

    @mcp.tool()
    def summarize_missing_values(
        path: str,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        scored_path: str | None = None,
        fields: list[str] | None = None,
        limit: int | None = None,
    ) -> dict[str, Any]:
        """Count missing values and missing fractions for selected or inferred fields."""
        return summarize_missing_values_tool(path=path, adapter=adapter, column_mapping_path=column_mapping_path, scored_path=scored_path, fields=fields, limit=limit)

    @mcp.tool()
    def summarize_cv_by_category(
        path: str,
        value_field: str,
        category_field: str,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        scored_path: str | None = None,
        top_n: int | None = None,
        limit: int | None = None,
    ) -> dict[str, Any]:
        """Summarize mean, standard deviation, and coefficient of variation across categories."""
        return summarize_cv_by_category_tool(
            path=path,
            value_field=value_field,
            category_field=category_field,
            adapter=adapter,
            column_mapping_path=column_mapping_path,
            scored_path=scored_path,
            top_n=top_n,
            limit=limit,
        )

    @mcp.tool()
    def bin_numeric_field(
        path: str,
        field: str,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        scored_path: str | None = None,
        group_by: str | None = None,
        bin_edges: list[float] | None = None,
        bin_count: int | None = None,
        range_min: float | None = None,
        range_max: float | None = None,
        weight_field: str | None = None,
        limit: int | None = None,
    ) -> dict[str, Any]:
        """Return histogram-style bins for one numeric field with optional grouping and weighting.

        bin_count equal-width bins span range_min to range_max; either may be given on its own, and the end not
        given is the data's minimum or maximum. bin_edges, when given, are used instead.
        """
        value_range = (range_min, range_max) if range_min is not None or range_max is not None else None
        return bin_numeric_field_tool(
            path=path,
            field=field,
            adapter=adapter,
            column_mapping_path=column_mapping_path,
            scored_path=scored_path,
            group_by=group_by,
            bin_edges=bin_edges,
            bin_count=bin_count,
            value_range=value_range,
            weight_field=weight_field,
            limit=limit,
        )

    @mcp.tool()
    def bin_2d_numeric_fields(
        path: str,
        x_field: str,
        y_field: str,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        scored_path: str | None = None,
        x_bin_edges: list[float] | None = None,
        y_bin_edges: list[float] | None = None,
        x_bin_count: int | None = None,
        y_bin_count: int | None = None,
        weight_field: str | None = None,
        limit: int | None = None,
    ) -> dict[str, Any]:
        """Return a 2D bin grid for two numeric fields, suitable for downstream plotting."""
        return bin_2d_numeric_fields_tool(
            path=path,
            x_field=x_field,
            y_field=y_field,
            adapter=adapter,
            column_mapping_path=column_mapping_path,
            scored_path=scored_path,
            x_bin_edges=x_bin_edges,
            y_bin_edges=y_bin_edges,
            x_bin_count=x_bin_count,
            y_bin_count=y_bin_count,
            weight_field=weight_field,
            limit=limit,
        )

    @mcp.tool()
    def summarize_dashboard(
        path: str,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        scored_path: str | None = None,
        rt_field: str | None = None,
        mz_field: str = "precursor_mz",
        confidence_field: str | None = "confidence",
        quant_field: str | None = None,
        category_field: str = "experiment_name",
        limit: int | None = None,
    ) -> dict[str, Any]:
        """Build a dashboard-oriented summary bundle from a normalized ID table with optional scored enrichment."""
        return summarize_dashboard_tool(
            path=path,
            adapter=adapter,
            column_mapping_path=column_mapping_path,
            scored_path=scored_path,
            rt_field=rt_field,
            mz_field=mz_field,
            confidence_field=confidence_field,
            quant_field=quant_field,
            category_field=category_field,
            limit=limit,
        )

    @mcp.tool()
    def summarize_overlap(
        path: str,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        scored_path: str | None = None,
        entity: str = "peptide",
        group_field: str = "run_key",
        limit: int | None = None,
    ) -> dict[str, Any]:
        """Summarize overlap across runs or groups for peptides, peptidoforms, or precursors."""
        return summarize_overlap_tool(
            path=path,
            adapter=adapter,
            column_mapping_path=column_mapping_path,
            scored_path=scored_path,
            entity=entity,
            group_field=group_field,
            limit=limit,
        )

    @mcp.tool()
    def summarize_modifications(
        path: str,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        scored_path: str | None = None,
        group_field: str = "run_key",
        limit: int | None = None,
    ) -> dict[str, Any]:
        """Summarize modification counts, sites, residues, and per-group modification support."""
        return summarize_modifications_tool(
            path=path,
            adapter=adapter,
            column_mapping_path=column_mapping_path,
            scored_path=scored_path,
            group_field=group_field,
            limit=limit,
        )

    @mcp.tool()
    def summarize_prediction_accuracy(
        path: str,
        target_field: str,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        scored_path: str | None = None,
        confidence_field: str | None = "confidence",
        mass_error_field: str | None = None,
        confidence_thresholds: list[float] | None = None,
        reliability_bin_count: int = 10,
        limit: int | None = None,
    ) -> dict[str, Any]:
        """Summarize prediction correctness, calibration, and operating points against a target sequence field.

        Prediction and target are both normalised to UNIMOD notation before they are compared, so a target written as
        ALPM[+15.9949]HIR equals the prediction ALPM[UNIMOD:35]HIR."""
        return summarize_prediction_accuracy_tool(
            path=path,
            target_field=target_field,
            adapter=adapter,
            column_mapping_path=column_mapping_path,
            scored_path=scored_path,
            confidence_field=confidence_field,
            mass_error_field=mass_error_field,
            confidence_thresholds=confidence_thresholds,
            reliability_bin_count=reliability_bin_count,
            limit=limit,
        )

    @mcp.tool()
    def render_prediction_accuracy_html(summary: dict[str, Any], output_path: str, title: str = "Prediction accuracy") -> dict[str, Any]:
        """Render a prediction-accuracy summary payload as standalone HTML."""
        return render_prediction_accuracy_html_tool(summary, output_path=_resolve(output_path), title=title)

    @mcp.tool()
    def summarize_quant_quality(
        path: str,
        quant_field: str,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        scored_path: str | None = None,
        entity: str = "peptide",
        group_field: str = "run_key",
        min_groups_for_cv: int = 2,
        limit: int | None = None,
    ) -> dict[str, Any]:
        """Summarize quant completeness and per-entity variability across runs or groups."""
        return summarize_quant_quality_tool(
            path=path,
            quant_field=quant_field,
            adapter=adapter,
            column_mapping_path=column_mapping_path,
            scored_path=scored_path,
            entity=entity,
            group_field=group_field,
            min_groups_for_cv=min_groups_for_cv,
            limit=limit,
        )

    @mcp.tool()
    def summarize_file_coverage(
        path: str,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        scored_path: str | None = None,
        group_field: str = "run_key",
        mz_field: str = "precursor_mz",
        rt_field: str | None = None,
        confidence_field: str | None = None,
        quant_field: str | None = None,
        limit: int | None = None,
    ) -> dict[str, Any]:
        """Summarize per-file or per-run ID coverage, histograms, and metadata-derived support."""
        return summarize_file_coverage_tool(
            path=path,
            adapter=adapter,
            column_mapping_path=column_mapping_path,
            scored_path=scored_path,
            group_field=group_field,
            mz_field=mz_field,
            rt_field=rt_field,
            confidence_field=confidence_field,
            quant_field=quant_field,
            limit=limit,
        )

    @mcp.tool()
    def render_psm_table(
        path: str,
        output_path: str,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        scored_path: str | None = None,
        sdrf_path: str | None = None,
        pride_project_accession: str | None = None,
        timeout_seconds: float = 20.0,
        limit: int | None = None,
    ) -> dict[str, Any]:
        """Render the normalized PSM table as browsable HTML with in-page CSV download."""
        return render_psm_table_tool(
            path=path,
            output_path=_resolve(output_path),
            adapter=adapter,
            column_mapping_path=column_mapping_path,
            scored_path=scored_path,
            sdrf_path=sdrf_path,
            pride_project_accession=pride_project_accession,
            timeout_seconds=timeout_seconds,
            limit=limit,
        )

    @mcp.tool()
    def render_peptide_table(
        path: str,
        output_path: str,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        scored_path: str | None = None,
        sdrf_path: str | None = None,
        pride_project_accession: str | None = None,
        timeout_seconds: float = 20.0,
        limit: int | None = None,
    ) -> dict[str, Any]:
        """Render the collapsed peptide table as browsable HTML with in-page CSV download."""
        return render_peptide_table_tool(
            path=path,
            output_path=_resolve(output_path),
            adapter=adapter,
            column_mapping_path=column_mapping_path,
            scored_path=scored_path,
            sdrf_path=sdrf_path,
            pride_project_accession=pride_project_accession,
            timeout_seconds=timeout_seconds,
            limit=limit,
        )

    @mcp.tool()
    def render_summary_html(summary: dict[str, Any], output_path: str, title: str = "Summary") -> dict[str, Any]:
        """Render any structured summary payload as simple inspection HTML."""
        return render_summary_html_tool(summary, output_path=_resolve(output_path), title=title)

    @mcp.tool()
    def render_dashboard_html(summary: dict[str, Any], output_path: str, title: str = "Dashboard") -> dict[str, Any]:
        """Render a dashboard summary bundle as standalone interactive HTML."""
        return render_dashboard_html_tool(summary, output_path=_resolve(output_path), title=title)

    @mcp.tool()
    def build_id_dashboard(
        path: str,
        output_prefix: str,
        spectra_path: str | None = None,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        scored_path: str | None = None,
        top_n: int | None = None,
        rt_source_field: str | None = None,
        rt_unit: str = "seconds",
        confidence_field: str | None = "confidence",
        mz_field: str = "precursor_mz",
        default_experiment_name: str | None = None,
        limit: int | None = None,
    ) -> dict[str, Any]:
        """Build a full HTML dashboard bundle from a normalized ID table, with optional spectra metadata and scored enrichment.

        With spectra_path, candidates are joined to spectra by run and scan number; metadata_join in the result counts the
        matches and the unmatched rows by reason. rt_source_field names the spectra table's retention-time column
        (default: rt, retention_time, retention_time_seconds or RT, whichever it has), in rt_unit ("seconds" or "minutes")."""
        return build_id_dashboard_tool(
            path=path,
            spectra_path=spectra_path,
            output_prefix=_resolve(output_prefix),
            adapter=adapter,
            column_mapping_path=column_mapping_path,
            scored_path=scored_path,
            top_n=top_n,
            rt_source_field=rt_source_field,
            rt_unit=rt_unit,
            confidence_field=confidence_field,
            mz_field=mz_field,
            default_experiment_name=default_experiment_name,
            limit=limit,
        )

    @mcp.tool()
    def build_prediction_batch(
        path: str,
        output_prefix: str,
        spectra_path: str | None = None,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
        predictor_backend: str | None = None,
        predictor_model: str | None = None,
        ms2pip_model: str | None = None,
        max_fragment_charge: int | None = None,
        top_n: int | None = None,
        default_experiment_name: str | None = None,
        limit: int | None = None,
        config_path: str | None = None,
        collision_energy: float | None = None,
    ) -> dict[str, Any]:
        """Build a reusable fragment-prediction batch with manifest, fragment, mapping, failure, and status artifacts.

        Koina predictions use collision_energy when given; otherwise the spectrum's own normalized collision energy when it states one (a Thermo scan filter such as @hcd28.00), otherwise koina_collision_energy from the config (30). eV values from other vendors are not used."""
        return build_prediction_batch_tool(
            path=path,
            output_prefix=_resolve(output_prefix),
            spectra_path=_resolve(spectra_path),
            adapter=adapter,
            column_mapping_path=column_mapping_path,
            predictor_backend=predictor_backend,
            predictor_model=predictor_model,
            ms2pip_model=ms2pip_model,
            max_fragment_charge=max_fragment_charge,
            top_n=top_n,
            default_experiment_name=default_experiment_name,
            limit=limit,
            config_path=config_path,
            collision_energy=collision_energy,
        )

    @mcp.tool()
    def inspect_prediction_batch(manifest_path: str) -> dict[str, Any]:
        """Inspect the static artifacts of a completed or partially completed prediction batch."""
        return inspect_prediction_batch_tool(manifest_path)

    @mcp.tool()
    def inspect_prediction_batch_status(manifest_path: str) -> dict[str, Any]:
        """Poll the persistent status file for a long-running prediction batch."""
        return inspect_prediction_batch_status_tool(manifest_path)

    @mcp.tool()
    def score_from_prediction_batch(
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
        """Score a candidate subset against local spectra by reusing a prediction batch instead of recomputing fragments."""
        return score_from_prediction_batch_tool(
            manifest_path=manifest_path,
            spectra_path=spectra_path,
            output_prefix=_resolve(output_prefix),
            tolerance=tolerance,
            tolerance_unit=tolerance_unit,
            annotation_ion_types=annotation_ion_types,
            annotation_max_isotope=annotation_max_isotope,
            annotation_max_ion_charge=annotation_max_ion_charge,
            annotation_neutral_losses=annotation_neutral_losses,
            config_path=config_path,
        )

    @mcp.tool()
    def inspect_scored_batch(manifest_path: str) -> dict[str, Any]:
        """Inspect the scored output table, failure rows, and summary statistics of a scored batch."""
        return inspect_scored_batch_tool(manifest_path)

    @mcp.tool()
    def inspect_scored_batch_status(manifest_path: str) -> dict[str, Any]:
        """Poll the persistent status file for a long-running scored batch."""
        return inspect_scored_batch_status_tool(manifest_path)

    @mcp.tool()
    def predict_fragments(
        peptidoform: str,
        charge: int,
        model: str | None = None,
        backend: str | None = None,
        max_fragment_charge: int | None = None,
        config_path: str | None = None,
        collision_energy: float | None = None,
    ) -> dict[str, Any]:
        """Predict fragment ions and intensities for one peptidoform and precursor charge.

        Koina predictions use collision_energy when given; otherwise the spectrum's own normalized collision energy when it states one (a Thermo scan filter such as @hcd28.00), otherwise koina_collision_energy from the config (30). eV values from other vendors are not used."""
        return predict_fragments_tool(peptidoform, charge, model, backend, max_fragment_charge=max_fragment_charge, config_path=config_path, collision_energy=collision_energy)

    @mcp.tool()
    def retrieve_usi_spectrum(usi: str, backend: str = "pride", timeout_seconds: float = 15.0) -> dict[str, Any]:
        """Retrieve one observed spectrum from a USI provider and normalize it into the local observed-spectrum schema."""
        return retrieve_usi_spectrum_tool(usi, backend, timeout_seconds)

    @mcp.tool()
    def load_sdrf(path: str) -> dict[str, Any]:
        """Load and summarize a local SDRF-Proteomics TSV file."""
        return load_sdrf_tool(path)

    @mcp.tool()
    def retrieve_pride_sdrf(project_accession: str, timeout_seconds: float = 20.0) -> dict[str, Any]:
        """Retrieve and summarize the SDRF metadata for one PRIDE project accession."""
        return retrieve_pride_sdrf_tool(project_accession, timeout_seconds)

    @mcp.tool()
    def mzml_to_feather(mzml_path: str, output_path: str, ms_level: int = 2) -> dict[str, Any]:
        """Convert an mzML file into the normalized Feather spectrum table used by local workflows."""
        return mzml_to_feather_tool(mzml_path, _resolve(output_path), ms_level)

    # The parser is the one on PATH or in $PROTEOMICS_MCP_THERMORAWFILEPARSER, never a tool argument: an MCP client
    # must not choose which program runs, or with which arguments (#172)
    @mcp.tool()
    def thermo_raw_to_mzml(raw_path: str, output_path: str, timeout_seconds: float = 3600.0) -> dict[str, Any]:
        """Convert a Thermo `.raw` file to mzML with ThermoRawFileParser, installed separately (no Thermo libraries are bundled)."""
        return thermo_raw_to_mzml_tool(raw_path, _resolve(output_path), timeout_seconds=timeout_seconds)

    @mcp.tool()
    def thermo_raw_to_feather(
        raw_path: str, output_path: str, mzml_output_path: str | None = None, ms_level: int = 2, timeout_seconds: float = 3600.0
    ) -> dict[str, Any]:
        """Convert a Thermo `.raw` file to the normalized Feather spectrum table, through mzML from ThermoRawFileParser (installed separately); the mzML is kept."""
        return thermo_raw_to_feather_tool(
            raw_path, _resolve(output_path), mzml_output_path=_resolve(mzml_output_path), ms_level=ms_level, timeout_seconds=timeout_seconds
        )

    @mcp.tool()
    def mgf_to_feather(mgf_path: str, output_path: str) -> dict[str, Any]:
        """Convert an MGF file into the normalized Feather spectrum table used by local workflows."""
        return mgf_to_feather_tool(mgf_path, _resolve(output_path))

    @mcp.tool()
    def bruker_d_to_feather(d_path: str, output_path: str) -> dict[str, Any]:
        """Convert a Bruker `.d` directory into the normalized Feather spectrum table, preferring AlphaTims direct extraction."""
        return bruker_d_to_feather_tool(d_path, _resolve(output_path))

    @mcp.tool()
    def export_to_instanovo(
        input_path: str,
        output_path: str,
        output_format: str = "ipc",
        sequence_column: str | None = None,
        drop_unlabelled: bool = False,
    ) -> dict[str, Any]:
        """Export a normalized spectrum table into an InstaNovo-ready IPC or Parquet layout.

        With sequence_column, spectra without a sequence are written with a null sequence (InstaNovo then refuses the
        table as annotated data); set drop_unlabelled=true to leave them out for a labelled training or evaluation set.
        """
        return export_to_instanovo_tool(
            input_path, _resolve(output_path), output_format=output_format, sequence_column=sequence_column, drop_unlabelled=drop_unlabelled
        )

    @mcp.tool()
    def annotate_observed_spectrum(
        observed: dict[str, Any],
        peptidoform: str,
        ion_types: str = "by",
        tolerance: float = 0.05,
        tolerance_unit: str = "Da",
        max_isotope: int = 2,
        max_ion_charge: int | None = None,
        neutral_losses: bool = False,
    ) -> dict[str, Any]:
        """Annotate an observed spectrum against a candidate peptidoform without predicting a mirror spectrum."""
        return annotate_observed_spectrum_tool(
            observed,
            peptidoform=peptidoform,
            ion_types=ion_types,
            tolerance=tolerance,
            tolerance_unit=tolerance_unit,
            max_isotope=max_isotope,
            max_ion_charge=max_ion_charge,
            neutral_losses=neutral_losses,
        )

    @mcp.tool()
    def score_spectrum(
        candidate: dict[str, Any],
        observed: dict[str, Any],
        predicted: dict[str, Any] | None = None,
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
        score_predictions_as_given: bool | None = None,
        config_path: str | None = None,
        collision_energy: float | None = None,
    ) -> dict[str, Any]:
        """Score one candidate against one observed spectrum, optionally using a supplied predicted spectrum.

        A supplied prediction goes through the same rules as one the package makes: fragments above the fragment-charge
        limit (max_fragment_charge, or precursor charge - 1) are left out, unless score_predictions_as_given. The result
        lists the scored ions, the ions left out with the reason, and the settings in effect.

        Koina predictions use collision_energy when given; otherwise the spectrum's own normalized collision energy when it states one (a Thermo scan filter such as @hcd28.00), otherwise koina_collision_energy from the config (30). eV values from other vendors are not used."""
        return score_spectrum_tool(
            candidate,
            observed,
            predicted,
            predictor_backend=predictor_backend,
            predictor_model=predictor_model,
            ms2pip_model=ms2pip_model,
            max_fragment_charge=max_fragment_charge,
            tolerance=tolerance,
            tolerance_unit=tolerance_unit,
            annotation_ion_types=annotation_ion_types,
            annotation_max_isotope=annotation_max_isotope,
            annotation_max_ion_charge=annotation_max_ion_charge,
            annotation_neutral_losses=annotation_neutral_losses,
            score_predictions_as_given=score_predictions_as_given,
            config_path=config_path,
            collision_energy=collision_energy,
        )

    @mcp.tool()
    def describe_predictor_model(
        model: str,
        backend: str = "koina",
        server_url: str = "https://koina.wilhelmlab.org",
        timeout_seconds: float = 30.0,
    ) -> dict[str, Any]:
        """Inspect a predictor model before scoring; for Koina this retrieves its live input/output schema."""
        return describe_predictor_model_tool(
            model,
            backend=backend,
            server_url=server_url,
            timeout_seconds=timeout_seconds,
        )

    @mcp.tool()
    def score_candidates(
        candidates_path: str,
        spectra_path: str | None = None,
        output_path: str | None = None,
        adapter: str = "instanovo",
        column_mapping_path: str | None = None,
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
        collision_energy: float | None = None,
    ) -> dict[str, Any]:
        """Score many candidates directly against local spectra or USI targets without a reusable prediction batch.

        Koina predictions use collision_energy when given; otherwise the spectrum's own normalized collision energy when it states one (a Thermo scan filter such as @hcd28.00), otherwise koina_collision_energy from the config (30). eV values from other vendors are not used."""
        return score_candidates_tool(
            candidates_path=candidates_path,
            spectra_path=spectra_path,
            output_path=_resolve(output_path),
            adapter=adapter,
            column_mapping_path=column_mapping_path,
            predictor_backend=predictor_backend,
            predictor_model=predictor_model,
            ms2pip_model=ms2pip_model,
            max_fragment_charge=max_fragment_charge,
            tolerance=tolerance,
            tolerance_unit=tolerance_unit,
            annotation_ion_types=annotation_ion_types,
            annotation_max_isotope=annotation_max_isotope,
            annotation_max_ion_charge=annotation_max_ion_charge,
            annotation_neutral_losses=annotation_neutral_losses,
            limit=limit,
            status_path=_resolve(status_path),
            config_path=config_path,
            collision_energy=collision_energy,
        )

    @mcp.tool()
    def summarize_score_metrics(score_result: dict[str, Any]) -> dict[str, Any]:
        """Collapse a raw spectrum score result into agent-friendly metric groups and evidence summaries."""
        return summarize_score_metrics_tool(score_result)

    @mcp.tool()
    def summarize_dataset(
        scored_path: str,
        output_prefix: str | None = None,
        sdrf_path: str | None = None,
        pride_project_accession: str | None = None,
        timeout_seconds: float = 20.0,
    ) -> dict[str, Any]:
        """Aggregate a scored PSM table into dataset, peptide, peptidoform, run, and optional sample summaries."""
        return summarize_dataset_tool(
            scored_path=scored_path,
            output_prefix=_resolve(output_prefix),
            sdrf_path=sdrf_path,
            pride_project_accession=pride_project_accession,
            timeout_seconds=timeout_seconds,
        )

    @mcp.tool()
    def render_dataset_report(
        scored_path: str,
        output_path: str | None = None,
        sdrf_path: str | None = None,
        pride_project_accession: str | None = None,
        timeout_seconds: float = 20.0,
    ) -> dict[str, Any]:
        """Render a scored dataset aggregation as standalone HTML."""
        return render_dataset_report_tool(
            scored_path=scored_path,
            output_path=_resolve(output_path),
            sdrf_path=sdrf_path,
            pride_project_accession=pride_project_accession,
            timeout_seconds=timeout_seconds,
        )

    @mcp.tool()
    def list_dataset_groups(
        scored_path: str,
        sdrf_path: str | None = None,
        pride_project_accession: str | None = None,
        timeout_seconds: float = 20.0,
    ) -> dict[str, Any]:
        """List dataset, run, and sample grouping keys available in a scored dataset summary."""
        return list_dataset_groups_tool(
            scored_path=scored_path,
            sdrf_path=sdrf_path,
            pride_project_accession=pride_project_accession,
            timeout_seconds=timeout_seconds,
        )

    @mcp.tool()
    def get_peptide_summary(
        scored_path: str,
        peptide: str,
        sdrf_path: str | None = None,
        pride_project_accession: str | None = None,
        timeout_seconds: float = 20.0,
    ) -> dict[str, Any]:
        """Return one peptide-level aggregate row from a scored dataset."""
        return get_peptide_summary_tool(
            scored_path=scored_path,
            peptide=peptide,
            sdrf_path=sdrf_path,
            pride_project_accession=pride_project_accession,
            timeout_seconds=timeout_seconds,
        )

    @mcp.tool()
    def get_peptidoform_summary(
        scored_path: str,
        peptidoform: str,
        sdrf_path: str | None = None,
        pride_project_accession: str | None = None,
        timeout_seconds: float = 20.0,
    ) -> dict[str, Any]:
        """Return one modification-aware peptidoform aggregate row from a scored dataset."""
        return get_peptidoform_summary_tool(
            scored_path=scored_path,
            peptidoform=peptidoform,
            sdrf_path=sdrf_path,
            pride_project_accession=pride_project_accession,
            timeout_seconds=timeout_seconds,
        )

    @mcp.tool()
    def render_spectrum_comparison(payload: dict[str, Any], output_path: str | None = None) -> dict[str, Any]:
        """Render the mirror plot view with observed and predicted spectra plus spectrum-level metrics."""
        return render_spectrum_comparison_tool(payload, _resolve(output_path))

    @mcp.tool()
    def render_observed_spectrum(payload_or_observed: dict[str, Any], output_path: str | None = None) -> dict[str, Any]:
        """Render the observed-only spectrum view without annotations or a predicted mirror spectrum."""
        return render_observed_spectrum_tool(payload_or_observed, _resolve(output_path))

    @mcp.tool()
    def render_annotated_observed_spectrum(
        payload_or_observed: dict[str, Any],
        output_path: str | None = None,
        peptidoform: str | None = None,
        ion_types: str = "by",
        tolerance: float = 0.05,
        tolerance_unit: str = "Da",
        max_isotope: int = 2,
        max_ion_charge: int | None = None,
        neutral_losses: bool = False,
    ) -> dict[str, Any]:
        """Render the observed-only annotated spectrum view using the internal ion-series annotator."""
        return render_annotated_observed_spectrum_tool(
            payload_or_observed,
            _resolve(output_path),
            peptidoform=peptidoform,
            ion_types=ion_types,
            tolerance=tolerance,
            tolerance_unit=tolerance_unit,
            max_isotope=max_isotope,
            max_ion_charge=max_ion_charge,
            neutral_losses=neutral_losses,
        )

    @mcp.tool()
    def render_static_spectrum(
        payload_or_observed: dict[str, Any],
        output_path: str,
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
        """Write a publication-ready static mirror, annotated, or plain spectrum as SVG, PDF, or PNG from a normalized payload."""
        return render_static_spectrum_tool(
            payload_or_observed,
            _resolve(output_path) or output_path,
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

    @mcp.tool()
    def render_spectra_batch(
        items: list[dict[str, Any]] | None = None,
        output_dir: str | None = None,
        input_path: str | None = None,
        views: list[str] | None = None,
        max_workers: int | None = None,
        validate_only: bool = False,
        status_path: str | None = None,
    ) -> dict[str, Any]:
        """Render many spectrum HTML artifacts at once from inline items or a file-backed batch manifest.

        input_path may be JSON (a list of items, or {"items": [...], "views": [...]}), JSON Lines
        (.jsonl or .ndjson, one item per line), or either of them gzip-compressed (.gz).
        """
        return render_spectra_batch_tool(
            items=items,
            output_dir=_resolve(output_dir),
            input_path=_resolve(input_path),
            views=views,
            max_workers=max_workers,
            validate_only=validate_only,
            status_path=_resolve(status_path),
        )

    @mcp.resource("method://metric_profile_v2")
    def metric_profile_v2() -> str:
        return SCHEMAS["method://metric_profile_v2"]

    @mcp.resource("method://metric_profile_legacy_v1")
    def metric_profile_legacy_v1() -> str:
        return SCHEMAS["method://metric_profile_legacy_v1"]

    @mcp.resource("method://tooling_surface_v1")
    def tooling_surface_v1() -> str:
        return SCHEMAS["method://tooling_surface_v1"]

    @mcp.resource("method://candidate_schema_v1")
    def candidate_schema_v1() -> str:
        return SCHEMAS["method://candidate_schema_v1"]

    @mcp.resource("method://generic_mapping_v1")
    def generic_mapping_v1() -> str:
        return SCHEMAS["method://generic_mapping_v1"]

    @mcp.resource("schema://candidate_psm_v1")
    def schema_candidate_psm_v1() -> str:
        return SCHEMAS["schema://candidate_psm_v1"]

    @mcp.resource("schema://observed_spectrum_v1")
    def schema_observed_spectrum_v1() -> str:
        return SCHEMAS["schema://observed_spectrum_v1"]

    @mcp.resource("schema://predicted_spectrum_v1")
    def schema_predicted_spectrum_v1() -> str:
        return SCHEMAS["schema://predicted_spectrum_v1"]

    @mcp.resource("schema://spectrum_metric_result_v1")
    def schema_spectrum_metric_result_v1() -> str:
        return SCHEMAS["schema://spectrum_metric_result_v1"]

    @mcp.resource("config://proteomics_mcp_v1")
    def config_proteomics_mcp_v1() -> str:
        return SCHEMAS["config://proteomics_mcp_v1"]

    @mcp.resource("method://supported_predictors")
    def supported_predictors() -> str:
        return SCHEMAS["method://supported_predictors"]

    @mcp.resource("method://terminal_ion_annotation_v1")
    def terminal_ion_annotation_v1() -> str:
        return SCHEMAS["method://terminal_ion_annotation_v1"]

    @mcp.resource("method://visualization_modes_v1")
    def visualization_modes_v1() -> str:
        return SCHEMAS["method://visualization_modes_v1"]

    @mcp.resource("method://sdrf_support_v1")
    def sdrf_support_v1() -> str:
        return SCHEMAS["method://sdrf_support_v1"]

    @mcp.resource("method://dataset_summary_v1")
    def dataset_summary_v1() -> str:
        return SCHEMAS["method://dataset_summary_v1"]

    @mcp.resource("method://id_table_inspection_v1")
    def id_table_inspection_v1() -> str:
        return SCHEMAS["method://id_table_inspection_v1"]
else:
    mcp = None


def run() -> None:
    if mcp is None:
        raise RuntimeError("MCP server requires the optional 'mcp' package to be installed.")
    mcp.run()
