from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from proteomics_mcp.server import run as run_server
from proteomics_mcp.tools.workflows import (
    annotate_observed_spectrum_tool,
    get_peptide_summary_tool,
    get_peptidoform_summary_tool,
    inspect_peptide_table_tool,
    inspect_psm_table_tool,
    list_distribution_fields_tool,
    load_candidates_tool,
    load_sdrf_tool,
    list_dataset_groups_tool,
    bruker_d_to_feather_tool,
    export_to_instanovo_tool,
    summarize_categorical_field_tool,
    summarize_cv_by_category_tool,
    summarize_dashboard_tool,
    summarize_file_coverage_tool,
    summarize_missing_values_tool,
    summarize_modifications_tool,
    summarize_overlap_tool,
    summarize_prediction_accuracy_tool,
    summarize_quant_quality_tool,
    mzml_to_feather_tool,
    mgf_to_feather_tool,
    describe_predictor_model_tool,
    predict_fragments_tool,
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
)


def _read_json(path: str) -> Any:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="proteomics-mcp")
    sub = parser.add_subparsers(dest="command", required=True)
    p = sub.add_parser("load-candidates")
    p.add_argument("--path", required=True)
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--fixed-modifications", help='MaxQuant fixed modifications for --adapter maxquant, e.g. "Carbamidomethyl (C)"')
    p.add_argument("--limit", type=int)
    p = sub.add_parser("inspect-psm-table")
    p.add_argument("--path", required=True)
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--scored-path")
    p.add_argument("--sdrf-path")
    p.add_argument("--project-accession")
    p.add_argument("--timeout-seconds", type=float, default=20.0)
    p.add_argument("--limit", type=int, help="number of preview rows (default 10); the summaries cover the whole table")
    p.add_argument("--max-rows", type=int, help="read only the first rows of the table; the result then reports it")
    p = sub.add_parser("inspect-peptide-table")
    p.add_argument("--path", required=True)
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--scored-path")
    p.add_argument("--peptide")
    p.add_argument("--sdrf-path")
    p.add_argument("--project-accession")
    p.add_argument("--timeout-seconds", type=float, default=20.0)
    p.add_argument("--limit", type=int, help="number of preview rows (default 10); the summaries cover the whole table")
    p.add_argument("--max-rows", type=int, help="read only the first rows of the table; the result then reports it")
    p = sub.add_parser("render-psm-table")
    p.add_argument("--path", required=True)
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--scored-path")
    p.add_argument("--output-path", required=True)
    p.add_argument("--sdrf-path")
    p.add_argument("--project-accession")
    p.add_argument("--timeout-seconds", type=float, default=20.0)
    p.add_argument("--limit", type=int)
    p = sub.add_parser("render-peptide-table")
    p.add_argument("--path", required=True)
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--scored-path")
    p.add_argument("--output-path", required=True)
    p.add_argument("--sdrf-path")
    p.add_argument("--project-accession")
    p.add_argument("--timeout-seconds", type=float, default=20.0)
    p.add_argument("--limit", type=int)
    p = sub.add_parser("render-summary-html")
    p.add_argument("--input-json", required=True)
    p.add_argument("--output-path", required=True)
    p.add_argument("--title", default="Summary")
    p = sub.add_parser("render-dashboard-html")
    p.add_argument("--input-json", required=True)
    p.add_argument("--output-path", required=True)
    p.add_argument("--title", default="Dashboard")
    p = sub.add_parser("build-id-dashboard")
    p.add_argument("--path", required=True)
    p.add_argument("--spectra-path")
    p.add_argument("--output-prefix", required=True)
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--scored-path")
    p.add_argument("--top-n", type=int)
    p.add_argument("--rt-source-field", default="rt")
    p.add_argument("--confidence-field", default="confidence")
    p.add_argument("--mz-field", default="precursor_mz")
    p.add_argument("--default-experiment-name", help="Raw-file stem to use when an InstaNovo export lacks run provenance.")
    p.add_argument("--limit", type=int)
    p = sub.add_parser("build-prediction-batch")
    p.add_argument("--path", required=True)
    p.add_argument("--spectra-path", help="Optional Feather/IPC spectra table used to fill missing precursor metadata before prediction.")
    p.add_argument("--output-prefix", required=True)
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--predictor-backend")
    p.add_argument("--predictor-model")
    p.add_argument("--ms2pip-model")
    p.add_argument("--max-fragment-charge", type=int)
    p.add_argument("--top-n", type=int)
    p.add_argument("--default-experiment-name", help="Raw-file stem to use when an InstaNovo export lacks run provenance.")
    p.add_argument("--limit", type=int)
    p.add_argument("--config")
    p.add_argument("--collision-energy", type=float, help="collision energy for Koina predictions (default: the spectrum's Thermo NCE, else the config's koina_collision_energy)")
    p = sub.add_parser("inspect-prediction-batch")
    p.add_argument("--manifest-path", required=True)
    p = sub.add_parser("inspect-prediction-batch-status")
    p.add_argument("--manifest-path", required=True)
    p = sub.add_parser("score-from-prediction-batch")
    p.add_argument("--manifest-path", required=True)
    p.add_argument("--spectra-path", required=True)
    p.add_argument("--output-prefix", required=True)
    p.add_argument("--tolerance", type=float)
    p.add_argument("--tolerance-unit")
    p.add_argument("--annotation-ion-types")
    p.add_argument("--annotation-max-isotope", type=int)
    p.add_argument("--annotation-max-ion-charge", type=int)
    p.add_argument("--annotation-neutral-losses", action="store_true")
    p.add_argument("--config")
    p = sub.add_parser("inspect-scored-batch")
    p.add_argument("--manifest-path", required=True)
    p = sub.add_parser("inspect-scored-batch-status")
    p.add_argument("--manifest-path", required=True)
    p = sub.add_parser("list-distribution-fields")
    p.add_argument("--path", required=True)
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--scored-path")
    p.add_argument("--limit", type=int)
    p = sub.add_parser("summarize-numeric-field")
    p.add_argument("--path", required=True)
    p.add_argument("--field", required=True)
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--scored-path")
    p.add_argument("--group-by")
    p.add_argument("--weight-field")
    p.add_argument("--percentiles-json")
    p.add_argument("--limit", type=int)
    p = sub.add_parser("summarize-categorical-field")
    p.add_argument("--path", required=True)
    p.add_argument("--field", required=True)
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--scored-path")
    p.add_argument("--group-by")
    p.add_argument("--top-n", type=int)
    p.add_argument("--limit", type=int)
    p = sub.add_parser("summarize-missing-values")
    p.add_argument("--path", required=True)
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--scored-path")
    p.add_argument("--fields-json")
    p.add_argument("--limit", type=int)
    p = sub.add_parser("summarize-cv-by-category")
    p.add_argument("--path", required=True)
    p.add_argument("--value-field", required=True)
    p.add_argument("--category-field", required=True)
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--scored-path")
    p.add_argument("--top-n", type=int)
    p.add_argument("--limit", type=int)
    p = sub.add_parser("bin-numeric-field")
    p.add_argument("--path", required=True)
    p.add_argument("--field", required=True)
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--scored-path")
    p.add_argument("--group-by")
    p.add_argument("--weight-field")
    p.add_argument("--bin-edges-json")
    p.add_argument("--bin-count", type=int)
    p.add_argument("--range-json")
    p.add_argument("--limit", type=int)
    p = sub.add_parser("bin-2d-numeric-fields")
    p.add_argument("--path", required=True)
    p.add_argument("--x-field", required=True)
    p.add_argument("--y-field", required=True)
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--scored-path")
    p.add_argument("--x-bin-edges-json")
    p.add_argument("--y-bin-edges-json")
    p.add_argument("--x-bin-count", type=int)
    p.add_argument("--y-bin-count", type=int)
    p.add_argument("--weight-field")
    p.add_argument("--limit", type=int)
    p = sub.add_parser("summarize-dashboard")
    p.add_argument("--path", required=True)
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--scored-path")
    p.add_argument("--rt-field")
    p.add_argument("--mz-field", default="precursor_mz")
    p.add_argument("--confidence-field", default="confidence")
    p.add_argument("--quant-field")
    p.add_argument("--category-field", default="experiment_name")
    p.add_argument("--limit", type=int)
    p = sub.add_parser("summarize-overlap")
    p.add_argument("--path", required=True)
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--scored-path")
    p.add_argument("--entity", default="peptide")
    p.add_argument("--group-field", default="run_key")
    p.add_argument("--limit", type=int)
    p = sub.add_parser("summarize-modifications")
    p.add_argument("--path", required=True)
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--scored-path")
    p.add_argument("--group-field", default="run_key")
    p.add_argument("--limit", type=int)
    p = sub.add_parser("summarize-prediction-accuracy")
    p.add_argument("--path", required=True)
    p.add_argument("--target-field", required=True, help="Frame column or 'source_columns.<name>' holding ground-truth peptidoforms.")
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--scored-path")
    p.add_argument("--confidence-field", default="confidence")
    p.add_argument("--mass-error-field")
    p.add_argument("--reliability-bin-count", type=int, default=10)
    p.add_argument("--confidence-thresholds-json", help="Optional JSON list of confidence thresholds, e.g. '[0.99,0.9,0.5]'.")
    p.add_argument("--limit", type=int)
    p = sub.add_parser("render-prediction-accuracy-html")
    p.add_argument("--input-json", required=True)
    p.add_argument("--output-path", required=True)
    p.add_argument("--title", default="Prediction accuracy")
    p = sub.add_parser("summarize-quant-quality")
    p.add_argument("--path", required=True)
    p.add_argument("--quant-field", required=True)
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--scored-path")
    p.add_argument("--entity", default="peptide")
    p.add_argument("--group-field", default="run_key")
    p.add_argument("--min-groups-for-cv", type=int, default=2)
    p.add_argument("--limit", type=int)
    p = sub.add_parser("summarize-file-coverage")
    p.add_argument("--path", required=True)
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--scored-path")
    p.add_argument("--group-field", default="run_key")
    p.add_argument("--mz-field", default="precursor_mz")
    p.add_argument("--rt-field")
    p.add_argument("--confidence-field")
    p.add_argument("--quant-field")
    p.add_argument("--limit", type=int)
    p = sub.add_parser("predict-fragments")
    p.add_argument("--peptidoform", required=True)
    p.add_argument("--charge", required=True, type=int)
    p.add_argument("--model")
    p.add_argument("--backend")
    p.add_argument("--max-fragment-charge", type=int)
    p.add_argument("--config")
    p.add_argument("--collision-energy", type=float, help="collision energy for Koina predictions (default: the spectrum's Thermo NCE, else the config's koina_collision_energy)")
    p = sub.add_parser("describe-predictor-model")
    p.add_argument("--model", required=True)
    p.add_argument("--backend", default="koina")
    p.add_argument("--server-url", default="https://koina.wilhelmlab.org")
    p.add_argument("--timeout-seconds", type=float, default=30.0)
    p = sub.add_parser("retrieve-usi-spectrum")
    p.add_argument("--usi", required=True)
    p.add_argument("--backend", default="pride")
    p = sub.add_parser("load-sdrf")
    p.add_argument("--path", required=True)
    p = sub.add_parser("retrieve-pride-sdrf")
    p.add_argument("--project-accession", required=True)
    p.add_argument("--timeout-seconds", type=float, default=20.0)
    p = sub.add_parser("mzml-to-feather")
    p.add_argument("--mzml-path", required=True)
    p.add_argument("--output-path", required=True)
    p.add_argument("--ms-level", type=int, default=2)
    p = sub.add_parser("mgf-to-feather")
    p.add_argument("--mgf-path", required=True)
    p.add_argument("--output-path", required=True)
    p = sub.add_parser("bruker-d-to-feather")
    p.add_argument("--d-path", required=True)
    p.add_argument("--output-path", required=True)
    p = sub.add_parser("export-to-instanovo")
    p.add_argument("--input-path", required=True)
    p.add_argument("--output-path", required=True)
    p.add_argument("--format", default="ipc", choices=["ipc", "parquet"])
    p.add_argument("--sequence-column")
    p.add_argument("--drop-unlabelled", action="store_true", help="leave out spectra without a sequence (needs --sequence-column)")
    p = sub.add_parser("annotate-observed-spectrum")
    p.add_argument("--spectrum-json", required=True)
    p.add_argument("--peptidoform", required=True)
    p.add_argument("--ion-types", default="by")
    p.add_argument("--tolerance", type=float, default=0.05)
    p.add_argument("--tolerance-unit", default="Da")
    p.add_argument("--max-isotope", type=int, default=2)
    p.add_argument("--max-ion-charge", type=int)
    p.add_argument("--neutral-losses", action="store_true")
    p = sub.add_parser("score-spectrum")
    p.add_argument("--candidate-json", required=True)
    p.add_argument("--spectrum-json", required=True)
    p.add_argument("--predicted-json")
    p.add_argument("--predictor-backend")
    p.add_argument("--predictor-model")
    p.add_argument("--ms2pip-model")
    p.add_argument("--max-fragment-charge", type=int)
    p.add_argument("--tolerance", type=float)
    p.add_argument("--tolerance-unit")
    p.add_argument("--annotation-ion-types")
    p.add_argument("--annotation-max-isotope", type=int)
    p.add_argument("--annotation-max-ion-charge", type=int)
    p.add_argument("--annotation-neutral-losses", action="store_true")
    p.add_argument("--config")
    p.add_argument("--collision-energy", type=float, help="collision energy for Koina predictions (default: the spectrum's Thermo NCE, else the config's koina_collision_energy)")
    p = sub.add_parser("score-candidates")
    p.add_argument("--candidates-path", required=True)
    p.add_argument("--spectra-path")
    p.add_argument("--output-path")
    p.add_argument("--adapter", default="instanovo")
    p.add_argument("--column-mapping", dest="column_mapping_path", help="JSON column mapping for --adapter generic")
    p.add_argument("--predictor-backend")
    p.add_argument("--predictor-model")
    p.add_argument("--ms2pip-model")
    p.add_argument("--max-fragment-charge", type=int)
    p.add_argument("--tolerance", type=float)
    p.add_argument("--tolerance-unit")
    p.add_argument("--annotation-ion-types")
    p.add_argument("--annotation-max-isotope", type=int)
    p.add_argument("--annotation-max-ion-charge", type=int)
    p.add_argument("--annotation-neutral-losses", action="store_true")
    p.add_argument("--config")
    p.add_argument("--limit", type=int)
    p.add_argument("--status-path")
    p.add_argument("--collision-energy", type=float, help="collision energy for Koina predictions (default: the spectrum's Thermo NCE, else the config's koina_collision_energy)")
    p = sub.add_parser("summarize-dataset")
    p.add_argument("--scored-path", required=True)
    p.add_argument("--output-prefix")
    p.add_argument("--sdrf-path")
    p.add_argument("--project-accession")
    p.add_argument("--timeout-seconds", type=float, default=20.0)
    p = sub.add_parser("render-dataset-report")
    p.add_argument("--scored-path", required=True)
    p.add_argument("--output-path")
    p.add_argument("--sdrf-path")
    p.add_argument("--project-accession")
    p.add_argument("--timeout-seconds", type=float, default=20.0)
    p = sub.add_parser("list-dataset-groups")
    p.add_argument("--scored-path", required=True)
    p.add_argument("--sdrf-path")
    p.add_argument("--project-accession")
    p.add_argument("--timeout-seconds", type=float, default=20.0)
    p = sub.add_parser("get-peptide-summary")
    p.add_argument("--scored-path", required=True)
    p.add_argument("--peptide", required=True)
    p.add_argument("--sdrf-path")
    p.add_argument("--project-accession")
    p.add_argument("--timeout-seconds", type=float, default=20.0)
    p = sub.add_parser("get-peptidoform-summary")
    p.add_argument("--scored-path", required=True)
    p.add_argument("--peptidoform", required=True)
    p.add_argument("--sdrf-path")
    p.add_argument("--project-accession")
    p.add_argument("--timeout-seconds", type=float, default=20.0)
    p = sub.add_parser("render-spectrum-comparison")
    p.add_argument("--input-json", required=True)
    p.add_argument("--output-path")
    p = sub.add_parser("render-observed-spectrum")
    p.add_argument("--input-json", required=True)
    p.add_argument("--output-path")
    p = sub.add_parser("render-annotated-observed-spectrum")
    p.add_argument("--input-json", required=True)
    p.add_argument("--output-path")
    p.add_argument("--peptidoform")
    p.add_argument("--ion-types", default="by")
    p.add_argument("--tolerance", type=float, default=0.05)
    p.add_argument("--tolerance-unit", default="Da")
    p.add_argument("--max-isotope", type=int, default=2)
    p.add_argument("--max-ion-charge", type=int)
    p.add_argument("--neutral-losses", action="store_true")
    p = sub.add_parser("render-static-spectrum")
    p.add_argument("--input-json", required=True)
    p.add_argument("--output-path", required=True)
    p.add_argument("--view", choices=["mirror", "annotated", "plain"], default="mirror")
    p.add_argument("--peptidoform")
    p.add_argument("--ion-types", default="by")
    p.add_argument("--tolerance", type=float, default=0.05)
    p.add_argument("--tolerance-unit", default="Da")
    p.add_argument("--max-isotope", type=int, default=2)
    p.add_argument("--max-ion-charge", type=int)
    p.add_argument("--neutral-losses", action="store_true")
    p.add_argument("--include-fragment-sequences", action="store_true")
    p.add_argument("--label-top-n", type=int, help="Optional maximum number of matched-ion labels; omit to label all matches.")
    p.add_argument("--title")
    p = sub.add_parser("render-spectra-batch")
    group = p.add_mutually_exclusive_group(required=True)
    group.add_argument("--input-json")
    group.add_argument("--input-path")
    p.add_argument("--output-dir")
    p.add_argument("--views-json")
    p.add_argument("--max-workers", type=int)
    p.add_argument("--validate-only", action="store_true")
    p.add_argument("--status-path")
    sub.add_parser("serve")
    return parser


def main(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)
    if args.command == "load-candidates":
        print(json.dumps(load_candidates_tool(args.path, args.adapter, column_mapping_path=args.column_mapping_path,
                                              fixed_modifications=args.fixed_modifications, limit=args.limit), indent=2))
    elif args.command == "inspect-psm-table":
        print(
            json.dumps(
                inspect_psm_table_tool(
                    path=args.path,
                    adapter=args.adapter,
                    column_mapping_path=args.column_mapping_path,
                    scored_path=args.scored_path,
                    sdrf_path=args.sdrf_path,
                    pride_project_accession=args.project_accession,
                    timeout_seconds=args.timeout_seconds,
                    limit=args.limit,
                    max_rows=args.max_rows,
                ),
                indent=2,
            )
        )
    elif args.command == "inspect-peptide-table":
        print(
            json.dumps(
                inspect_peptide_table_tool(
                    path=args.path,
                    adapter=args.adapter,
                    column_mapping_path=args.column_mapping_path,
                    peptide=args.peptide,
                    scored_path=args.scored_path,
                    sdrf_path=args.sdrf_path,
                    pride_project_accession=args.project_accession,
                    timeout_seconds=args.timeout_seconds,
                    limit=args.limit,
                    max_rows=args.max_rows,
                ),
                indent=2,
            )
        )
    elif args.command == "render-psm-table":
        print(
            json.dumps(
                render_psm_table_tool(
                    path=args.path,
                    adapter=args.adapter,
                    column_mapping_path=args.column_mapping_path,
                    output_path=args.output_path,
                    scored_path=args.scored_path,
                    sdrf_path=args.sdrf_path,
                    pride_project_accession=args.project_accession,
                    timeout_seconds=args.timeout_seconds,
                    limit=args.limit,
                ),
                indent=2,
            )
        )
    elif args.command == "render-peptide-table":
        print(
            json.dumps(
                render_peptide_table_tool(
                    path=args.path,
                    adapter=args.adapter,
                    column_mapping_path=args.column_mapping_path,
                    output_path=args.output_path,
                    scored_path=args.scored_path,
                    sdrf_path=args.sdrf_path,
                    pride_project_accession=args.project_accession,
                    timeout_seconds=args.timeout_seconds,
                    limit=args.limit,
                ),
                indent=2,
            )
        )
    elif args.command == "render-summary-html":
        print(json.dumps(render_summary_html_tool(_read_json(args.input_json), args.output_path, title=args.title), indent=2))
    elif args.command == "render-dashboard-html":
        print(json.dumps(render_dashboard_html_tool(_read_json(args.input_json), args.output_path, title=args.title), indent=2))
    elif args.command == "build-id-dashboard":
        print(
            json.dumps(
                build_id_dashboard_tool(
                    path=args.path,
                    spectra_path=args.spectra_path,
                    output_prefix=args.output_prefix,
                    adapter=args.adapter,
                    column_mapping_path=args.column_mapping_path,
                    scored_path=args.scored_path,
                    top_n=args.top_n,
                    rt_source_field=args.rt_source_field,
                    confidence_field=args.confidence_field,
                    mz_field=args.mz_field,
                    default_experiment_name=args.default_experiment_name,
                    limit=args.limit,
                ),
                indent=2,
            )
        )
    elif args.command == "build-prediction-batch":
        print(
            json.dumps(
                build_prediction_batch_tool(
                    path=args.path,
                    output_prefix=args.output_prefix,
                    spectra_path=args.spectra_path,
                    adapter=args.adapter,
                    column_mapping_path=args.column_mapping_path,
                    predictor_backend=args.predictor_backend,
                    predictor_model=args.predictor_model,
                    ms2pip_model=args.ms2pip_model,
                    max_fragment_charge=args.max_fragment_charge,
                    top_n=args.top_n,
                    default_experiment_name=args.default_experiment_name,
                    limit=args.limit,
                    config_path=args.config,
                    collision_energy=args.collision_energy,
                ),
                indent=2,
            )
        )
    elif args.command == "inspect-prediction-batch":
        print(json.dumps(inspect_prediction_batch_tool(args.manifest_path), indent=2))
    elif args.command == "inspect-prediction-batch-status":
        print(json.dumps(inspect_prediction_batch_status_tool(args.manifest_path), indent=2))
    elif args.command == "score-from-prediction-batch":
        print(
            json.dumps(
                score_from_prediction_batch_tool(
                    manifest_path=args.manifest_path,
                    spectra_path=args.spectra_path,
                    output_prefix=args.output_prefix,
                    tolerance=args.tolerance,
                    tolerance_unit=args.tolerance_unit,
                    annotation_ion_types=args.annotation_ion_types,
                    annotation_max_isotope=args.annotation_max_isotope,
                    annotation_max_ion_charge=args.annotation_max_ion_charge,
                    annotation_neutral_losses=args.annotation_neutral_losses,
                    config_path=args.config,
                ),
                indent=2,
            )
        )
    elif args.command == "inspect-scored-batch":
        print(json.dumps(inspect_scored_batch_tool(args.manifest_path), indent=2))
    elif args.command == "inspect-scored-batch-status":
        print(json.dumps(inspect_scored_batch_status_tool(args.manifest_path), indent=2))
    elif args.command == "list-distribution-fields":
        print(json.dumps(list_distribution_fields_tool(path=args.path, adapter=args.adapter, column_mapping_path=args.column_mapping_path, scored_path=args.scored_path, limit=args.limit), indent=2))
    elif args.command == "summarize-numeric-field":
        print(
            json.dumps(
                summarize_numeric_field_tool(
                    path=args.path,
                    field=args.field,
                    adapter=args.adapter,
                    column_mapping_path=args.column_mapping_path,
                    scored_path=args.scored_path,
                    group_by=args.group_by,
                    weight_field=args.weight_field,
                    percentiles=json.loads(args.percentiles_json) if args.percentiles_json else None,
                    limit=args.limit,
                ),
                indent=2,
            )
        )
    elif args.command == "summarize-categorical-field":
        print(
            json.dumps(
                summarize_categorical_field_tool(
                    path=args.path,
                    field=args.field,
                    adapter=args.adapter,
                    column_mapping_path=args.column_mapping_path,
                    scored_path=args.scored_path,
                    group_by=args.group_by,
                    top_n=args.top_n,
                    limit=args.limit,
                ),
                indent=2,
            )
        )
    elif args.command == "summarize-missing-values":
        print(
            json.dumps(
                summarize_missing_values_tool(
                    path=args.path,
                    adapter=args.adapter,
                    column_mapping_path=args.column_mapping_path,
                    scored_path=args.scored_path,
                    fields=json.loads(args.fields_json) if args.fields_json else None,
                    limit=args.limit,
                ),
                indent=2,
            )
        )
    elif args.command == "summarize-cv-by-category":
        print(
            json.dumps(
                summarize_cv_by_category_tool(
                    path=args.path,
                    value_field=args.value_field,
                    category_field=args.category_field,
                    adapter=args.adapter,
                    column_mapping_path=args.column_mapping_path,
                    scored_path=args.scored_path,
                    top_n=args.top_n,
                    limit=args.limit,
                ),
                indent=2,
            )
        )
    elif args.command == "bin-numeric-field":
        print(
            json.dumps(
                bin_numeric_field_tool(
                    path=args.path,
                    field=args.field,
                    adapter=args.adapter,
                    column_mapping_path=args.column_mapping_path,
                    scored_path=args.scored_path,
                    group_by=args.group_by,
                    weight_field=args.weight_field,
                    bin_edges=json.loads(args.bin_edges_json) if args.bin_edges_json else None,
                    bin_count=args.bin_count,
                    value_range=tuple(json.loads(args.range_json)) if args.range_json else None,
                    limit=args.limit,
                ),
                indent=2,
            )
        )
    elif args.command == "bin-2d-numeric-fields":
        print(
            json.dumps(
                bin_2d_numeric_fields_tool(
                    path=args.path,
                    x_field=args.x_field,
                    y_field=args.y_field,
                    adapter=args.adapter,
                    column_mapping_path=args.column_mapping_path,
                    scored_path=args.scored_path,
                    x_bin_edges=json.loads(args.x_bin_edges_json) if args.x_bin_edges_json else None,
                    y_bin_edges=json.loads(args.y_bin_edges_json) if args.y_bin_edges_json else None,
                    x_bin_count=args.x_bin_count,
                    y_bin_count=args.y_bin_count,
                    weight_field=args.weight_field,
                    limit=args.limit,
                ),
                indent=2,
            )
        )
    elif args.command == "summarize-dashboard":
        print(
            json.dumps(
                summarize_dashboard_tool(
                    path=args.path,
                    adapter=args.adapter,
                    column_mapping_path=args.column_mapping_path,
                    scored_path=args.scored_path,
                    rt_field=args.rt_field,
                    mz_field=args.mz_field,
                    confidence_field=args.confidence_field,
                    quant_field=args.quant_field,
                    category_field=args.category_field,
                    limit=args.limit,
                ),
                indent=2,
            )
        )
    elif args.command == "summarize-overlap":
        print(
            json.dumps(
                summarize_overlap_tool(
                    path=args.path,
                    adapter=args.adapter,
                    column_mapping_path=args.column_mapping_path,
                    scored_path=args.scored_path,
                    entity=args.entity,
                    group_field=args.group_field,
                    limit=args.limit,
                ),
                indent=2,
            )
        )
    elif args.command == "summarize-modifications":
        print(
            json.dumps(
                summarize_modifications_tool(
                    path=args.path,
                    adapter=args.adapter,
                    column_mapping_path=args.column_mapping_path,
                    scored_path=args.scored_path,
                    group_field=args.group_field,
                    limit=args.limit,
                ),
                indent=2,
            )
        )
    elif args.command == "summarize-prediction-accuracy":
        thresholds = json.loads(args.confidence_thresholds_json) if args.confidence_thresholds_json else None
        print(
            json.dumps(
                summarize_prediction_accuracy_tool(
                    path=args.path,
                    target_field=args.target_field,
                    adapter=args.adapter,
                    column_mapping_path=args.column_mapping_path,
                    scored_path=args.scored_path,
                    confidence_field=args.confidence_field,
                    mass_error_field=args.mass_error_field,
                    confidence_thresholds=thresholds,
                    reliability_bin_count=args.reliability_bin_count,
                    limit=args.limit,
                ),
                indent=2,
            )
        )
    elif args.command == "render-prediction-accuracy-html":
        print(
            json.dumps(
                render_prediction_accuracy_html_tool(_read_json(args.input_json), args.output_path, title=args.title),
                indent=2,
            )
        )
    elif args.command == "summarize-quant-quality":
        print(
            json.dumps(
                summarize_quant_quality_tool(
                    path=args.path,
                    quant_field=args.quant_field,
                    adapter=args.adapter,
                    column_mapping_path=args.column_mapping_path,
                    scored_path=args.scored_path,
                    entity=args.entity,
                    group_field=args.group_field,
                    min_groups_for_cv=args.min_groups_for_cv,
                    limit=args.limit,
                ),
                indent=2,
            )
        )
    elif args.command == "summarize-file-coverage":
        print(
            json.dumps(
                summarize_file_coverage_tool(
                    path=args.path,
                    adapter=args.adapter,
                    column_mapping_path=args.column_mapping_path,
                    scored_path=args.scored_path,
                    group_field=args.group_field,
                    mz_field=args.mz_field,
                    rt_field=args.rt_field,
                    confidence_field=args.confidence_field,
                    quant_field=args.quant_field,
                    limit=args.limit,
                ),
                indent=2,
            )
        )
    elif args.command == "predict-fragments":
        print(
            json.dumps(
                predict_fragments_tool(
                    args.peptidoform,
                    args.charge,
                    args.model,
                    args.backend,
                    max_fragment_charge=args.max_fragment_charge,
                    config_path=args.config,
                    collision_energy=args.collision_energy,
                ),
                indent=2,
            )
        )
    elif args.command == "describe-predictor-model":
        print(
            json.dumps(
                describe_predictor_model_tool(
                    args.model,
                    backend=args.backend,
                    server_url=args.server_url,
                    timeout_seconds=args.timeout_seconds,
                ),
                indent=2,
            )
        )
    elif args.command == "retrieve-usi-spectrum":
        print(json.dumps(retrieve_usi_spectrum_tool(args.usi, args.backend), indent=2))
    elif args.command == "load-sdrf":
        print(json.dumps(load_sdrf_tool(args.path), indent=2))
    elif args.command == "retrieve-pride-sdrf":
        print(json.dumps(retrieve_pride_sdrf_tool(args.project_accession, args.timeout_seconds), indent=2))
    elif args.command == "mzml-to-feather":
        print(json.dumps(mzml_to_feather_tool(args.mzml_path, args.output_path, args.ms_level), indent=2))
    elif args.command == "mgf-to-feather":
        print(json.dumps(mgf_to_feather_tool(args.mgf_path, args.output_path), indent=2))
    elif args.command == "bruker-d-to-feather":
        print(json.dumps(bruker_d_to_feather_tool(args.d_path, args.output_path), indent=2))
    elif args.command == "export-to-instanovo":
        print(
            json.dumps(
                export_to_instanovo_tool(
                    args.input_path,
                    args.output_path,
                    output_format=args.format,
                    sequence_column=args.sequence_column,
                    drop_unlabelled=args.drop_unlabelled,
                ),
                indent=2,
            )
        )
    elif args.command == "annotate-observed-spectrum":
        print(
            json.dumps(
                annotate_observed_spectrum_tool(
                    _read_json(args.spectrum_json),
                    peptidoform=args.peptidoform,
                    ion_types=args.ion_types,
                    tolerance=args.tolerance,
                    tolerance_unit=args.tolerance_unit,
                    max_isotope=args.max_isotope,
                    max_ion_charge=args.max_ion_charge,
                    neutral_losses=args.neutral_losses,
                ),
                indent=2,
            )
        )
    elif args.command == "score-spectrum":
        predicted = _read_json(args.predicted_json) if args.predicted_json else None
        print(
            json.dumps(
                score_spectrum_tool(
                    _read_json(args.candidate_json),
                    _read_json(args.spectrum_json),
                    predicted,
                    predictor_backend=args.predictor_backend,
                    predictor_model=args.predictor_model,
                    ms2pip_model=args.ms2pip_model,
                    max_fragment_charge=args.max_fragment_charge,
                    tolerance=args.tolerance,
                    tolerance_unit=args.tolerance_unit,
                    annotation_ion_types=args.annotation_ion_types,
                    annotation_max_isotope=args.annotation_max_isotope,
                    annotation_max_ion_charge=args.annotation_max_ion_charge,
                    annotation_neutral_losses=args.annotation_neutral_losses or None,
                    config_path=args.config,
                    collision_energy=args.collision_energy,
                ),
                indent=2,
            )
        )
    elif args.command == "score-candidates":
        print(
            json.dumps(
                score_candidates_tool(
                    candidates_path=args.candidates_path,
                    spectra_path=args.spectra_path,
                    output_path=args.output_path,
                    adapter=args.adapter,
                    column_mapping_path=args.column_mapping_path,
                    predictor_backend=args.predictor_backend,
                    predictor_model=args.predictor_model,
                    ms2pip_model=args.ms2pip_model,
                    max_fragment_charge=args.max_fragment_charge,
                    tolerance=args.tolerance,
                    tolerance_unit=args.tolerance_unit,
                    annotation_ion_types=args.annotation_ion_types,
                    annotation_max_isotope=args.annotation_max_isotope,
                    annotation_max_ion_charge=args.annotation_max_ion_charge,
                    annotation_neutral_losses=args.annotation_neutral_losses or None,
                    config_path=args.config,
                    collision_energy=args.collision_energy,
                    limit=args.limit,
                    status_path=args.status_path,
                ),
                indent=2,
            )
        )
    elif args.command == "summarize-dataset":
        print(
            json.dumps(
                summarize_dataset_tool(
                    scored_path=args.scored_path,
                    output_prefix=args.output_prefix,
                    sdrf_path=args.sdrf_path,
                    pride_project_accession=args.project_accession,
                    timeout_seconds=args.timeout_seconds,
                ),
                indent=2,
            )
        )
    elif args.command == "render-dataset-report":
        print(
            json.dumps(
                render_dataset_report_tool(
                    scored_path=args.scored_path,
                    output_path=args.output_path,
                    sdrf_path=args.sdrf_path,
                    pride_project_accession=args.project_accession,
                    timeout_seconds=args.timeout_seconds,
                ),
                indent=2,
            )
        )
    elif args.command == "list-dataset-groups":
        print(
            json.dumps(
                list_dataset_groups_tool(
                    scored_path=args.scored_path,
                    sdrf_path=args.sdrf_path,
                    pride_project_accession=args.project_accession,
                    timeout_seconds=args.timeout_seconds,
                ),
                indent=2,
            )
        )
    elif args.command == "get-peptide-summary":
        print(
            json.dumps(
                get_peptide_summary_tool(
                    scored_path=args.scored_path,
                    peptide=args.peptide,
                    sdrf_path=args.sdrf_path,
                    pride_project_accession=args.project_accession,
                    timeout_seconds=args.timeout_seconds,
                ),
                indent=2,
            )
        )
    elif args.command == "get-peptidoform-summary":
        print(
            json.dumps(
                get_peptidoform_summary_tool(
                    scored_path=args.scored_path,
                    peptidoform=args.peptidoform,
                    sdrf_path=args.sdrf_path,
                    pride_project_accession=args.project_accession,
                    timeout_seconds=args.timeout_seconds,
                ),
                indent=2,
            )
        )
    elif args.command == "render-spectrum-comparison":
        print(json.dumps(render_spectrum_comparison_tool(_read_json(args.input_json), args.output_path), indent=2))
    elif args.command == "render-observed-spectrum":
        print(json.dumps(render_observed_spectrum_tool(_read_json(args.input_json), args.output_path), indent=2))
    elif args.command == "render-annotated-observed-spectrum":
        print(
            json.dumps(
                render_annotated_observed_spectrum_tool(
                    _read_json(args.input_json),
                    args.output_path,
                    peptidoform=args.peptidoform,
                    ion_types=args.ion_types,
                    tolerance=args.tolerance,
                    tolerance_unit=args.tolerance_unit,
                    max_isotope=args.max_isotope,
                    max_ion_charge=args.max_ion_charge,
                    neutral_losses=args.neutral_losses,
                ),
                indent=2,
            )
        )
    elif args.command == "render-static-spectrum":
        print(
            json.dumps(
                render_static_spectrum_tool(
                    _read_json(args.input_json),
                    args.output_path,
                    view=args.view,
                    peptidoform=args.peptidoform,
                    ion_types=args.ion_types,
                    tolerance=args.tolerance,
                    tolerance_unit=args.tolerance_unit,
                    max_isotope=args.max_isotope,
                    max_ion_charge=args.max_ion_charge,
                    neutral_losses=args.neutral_losses,
                    include_fragment_sequences=args.include_fragment_sequences,
                    label_top_n=args.label_top_n,
                    title=args.title,
                ),
                indent=2,
            )
        )
    elif args.command == "render-spectra-batch":
        items = None
        views = json.loads(args.views_json) if args.views_json else None
        if args.input_json:
            batch_payload = _read_json(args.input_json)
            if isinstance(batch_payload, dict) and "items" in batch_payload:
                items = batch_payload["items"]
                if views is None:
                    views = batch_payload.get("views")
            else:
                items = batch_payload
        print(
            json.dumps(
                render_spectra_batch_tool(
                    items=items,
                    output_dir=args.output_dir,
                    input_path=args.input_path,
                    views=views,
                    max_workers=args.max_workers,
                    validate_only=args.validate_only,
                    status_path=args.status_path,
                ),
                indent=2,
            )
        )
    elif args.command == "serve":
        run_server()
