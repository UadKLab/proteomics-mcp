SCHEMAS = {
    "method://tooling_surface_v1": """proteomics-mcp MCP semantic surface.

Tool categories:
- inspect_*: normalize and inspect upstream ID tables without rescoring everything
- summarize_* / bin_*: return structured distributions, counts, and grouped summaries
- render_*: produce HTML inspection artifacts from tables, summaries, or spectrum payloads
- build_* / inspect_*_batch / score_from_prediction_batch: orchestration tools for reusable batch prediction, scoring, and status polling
- convert_* and export_*: move local spectral data into normalized Feather or InstaNovo-ready formats
- score_* and annotate_*: targeted spectrum-level evidence tools for one candidate or a selected subset

Design intent:
- keep heavy operations explicit
- expose typed tool arguments rather than opaque JSON option blobs
- return structured dictionaries that agents can inspect, summarize, or pass into downstream render steps

Recommended operating pattern:
1. inspect or summarize an upstream table
2. pick one peptide/PSM or a constrained subset
3. predict/score/render only that selected target or batch
4. use batch status tools for long-running prediction and scoring workflows""",
    "method://metric_profile_legacy_v1": """Metric profile legacy_v1.

Spectrum-level metrics:
- spectral_angle
- cosine_similarity
- pearson_correlation
- xcorr
- log10_xcorr
- log2_xcorr
- hyperscore
- log2_hyperscore
- matched_ion_count
- total_predicted_ions
- matched_ion_fraction
- explained_intensity_fraction
- precursor_signal_fraction
- fragment_mass_error_mean_ppm
- fragment_mass_error_abs_mean_ppm
- fragment_mass_error_mean_da
- fragment_mass_error_abs_mean_da
- precursor_mass_error_ppm
- precursor_mass_error_da

The profile preserves the sidecar-style formulas where practical while normalizing all outputs into SpectrumMetricResult.

spectral_angle, cosine_similarity and pearson_correlation compare observed and predicted intensities over the matched
predicted ions only, so a candidate matching a single peak scores 1.0. Each predicted ion takes its nearest peak, even one
another ion also takes, so a shared peak counts more than once in matched ions, explained intensity and the hyperscore.
Kept to reproduce scores from 0.1.10 and earlier; the default profile is v2 (method://metric_profile_v2).""",
    "method://metric_profile_v2": """Metric profile v2 (default).

The same metrics as legacy_v1 (method://metric_profile_legacy_v1). spectral_angle, cosine_similarity and
pearson_correlation compare observed and predicted intensities over every predicted ion, with zero observed intensity
for predicted ions that matched no peak, so unexplained predicted signal lowers the score. When no predicted ion
matches, spectral_angle and cosine_similarity are 0.0 and pearson_correlation is null.

Matching is one to one: each observed peak explains at most one predicted ion, assigned closest first, so a peak shared
by isobaric ions (e.g. b3/y2) counts once and explained_intensity_fraction is at most 1.

Predicted ions below min_predicted_intensity (default 0.001) of the most intense predicted ion are not scored: predictors
give ions they expect to be absent a tiny intensity, and those would count in total_predicted_ions and could be matched
to noise. matched_ion_fraction is therefore relative to the ions predicted above that level. Set
min_predicted_intensity = 0 to score every predicted ion.

Select a profile with metric_profile in a JSON or TOML config file passed as config_path (CLI: --config).""",
    "method://candidate_schema_v1": "Canonical candidate schema normalizes upstream engines into peptide/peptidoform, scan/spectrum identifiers, raw score, optional 0-1 confidence, and token-confidence summaries. Supported adapters include instanovo, casanovo, maxquant, fragpipe, proteome_discoverer, and generic mapped tables.",
    "schema://candidate_psm_v1": """CandidatePSM fields:
candidate_id, source_engine, source_engine_version, experiment_name, scan_number, spectrum_id, usi, peptidoform, precursor_mz, precursor_charge, rank, raw_score, raw_score_name, next_beam, next_beam_margin, confidence, token_log_probabilities, token_log_probability_mean, token_log_probability_min, token_confidences, token_confidence_mean, token_confidence_min, token_confidence_min_position, token_confidence_min_token, source_columns.

For InstaNovo beam outputs, `next_beam` is the beam-1 confidence on a 0-1 scale and `next_beam_margin` is `confidence - next_beam`. Search-engine score deltas remain in source_columns rather than being presented as beams.""",
    "schema://observed_spectrum_v1": """ObservedSpectrum fields:
spectrum_id, scan_number, usi, precursor_mz, precursor_charge, mz, intensity, retention_time, source, metadata.

mz and intensity are parallel numeric arrays. IPC and Feather inputs are normalized into this schema. Feather inputs can be produced from mzML, MGF, or Bruker `.d` sources.""",
    "schema://predicted_spectrum_v1": """PredictedSpectrum fields:
peptidoform, precursor_charge, fragments, predictor_name, predictor_version, model, metadata.

Each fragment contains annotation, ion_type, ordinal, charge, mz, intensity, neutral_loss, isotope, metadata.""",
    "schema://spectrum_metric_result_v1": """SpectrumMetricResult fields:
candidate, spectral_angle, cosine_similarity, pearson_correlation, xcorr, log10_xcorr, log2_xcorr, hyperscore, log2_hyperscore, matched_ion_count, total_predicted_ions, matched_ion_fraction, explained_intensity_fraction, precursor_signal_fraction, fragment and precursor mass-error summaries, matched_peaks, metric_profile, predictor_metadata, annotation_settings, annotation_summary, acquisition_provenance.

MCP score_spectrum also includes metric_summary, an agent-oriented grouping with match_quality, predicted_match_evidence, observed_annotation_evidence, mass_error, candidate, and acquisition fields. A backward-compatible ion_evidence group is retained as a merged view. Acquisition summary includes fragmentation_method, isolation_target_mz, retention_time_seconds/minutes, collision_energy, FAIMS compensation voltage, injection time, scan_header, scan_resolution, instrument_model, and detector_type when available from the scan header or fallbacks.""",
    "config://proteomics_mcp_v1": """Config file keys, JSON or TOML:
predictor_backend, predictor_model, ms2pip_model (legacy alias), max_fragment_charge, koina_server_url, koina_timeout_seconds, koina_collision_energy, koina_fragmentation_type, koina_instrument_type, koina_batch_size, koina_cache_dir, tolerance, tolerance_unit, usi_backend, metric_profile, min_predicted_intensity, annotation_ion_types, annotation_max_isotope, annotation_max_ion_charge, annotation_neutral_losses.

If max_fragment_charge is omitted, MS2PIP predictions use charge <= precursor_charge - 1 where the selected/resolved model supports those ion types. For HCD/CID precursor charge > 2, the wrapper can resolve to HCDch2/CIDch2 to include charge-2 fragments. Koina-backed Prosit models additionally depend on collision energy, fragmentation type, instrument type, and the selected remote model identifier.

CLI options and explicit MCP tool arguments override config-file values when explicitly provided.""",
    "method://generic_mapping_v1": """Generic candidate mapping contract.

Minimum required information:
- peptide or peptidoform sequence
- one spectrum locator: USI, spectrum_id, or experiment_name/raw-file stem plus scan_number

Recommended additional fields:
- precursor_mz
- precursor_charge
- raw_score and raw_score_name
- optional confidence on a 0-1 scale or a known convertible score type
- optional next_beam / next_beam_margin equivalents
- optional token-level log probabilities or confidences

The generic adapter accepts user-specified column mappings for peptide, scan, spectrum_id, usi, charge, precursor m/z, score, next-score equivalents, and token-log-probability arrays. Generic mapped tables can be .csv, .tsv, or .xlsx.

Pass the mapping as column_mapping_path (MCP) or --column-mapping (CLI): a JSON file whose keys name the table's columns, for example
{"peptide_column": "Sequence", "scan_column": "Scan", "charge_column": "Charge", "score_column": "Score", "score_type": "confidence_0_1"}.
Every tool and command that takes adapter accepts it; it is required for adapter="generic" and rejected for other adapters.
Keys: peptide_column (required), scan_column, spectrum_id_column, usi_column, charge_column, precursor_mz_column, score_column, next_beam_column, next_beam_margin_column, next_score_column, next_score_margin_column, token_log_probabilities_column, score_type (confidence_0_1, posterior_error_probability, evalue, score_unscaled, unknown; default unknown), source_engine (default "generic").

For multi-run local data, prefer joining by experiment_name/raw-file stem plus scan_number rather than scan_number alone.""",
    "method://supported_predictors": """Supported predictors:
- ms2pip: local MS2PIP fragment-intensity prediction.
- koina: remote prediction through a Koina server, including Prosit-family models selected by predictor_model.
- mock: deterministic lightweight predictor for tests and offline workflow checks.

InstaNovo is treated as an upstream candidate generator in v1, not as a local inference backend. The canonical public Koina/Prosit path is backend='koina' plus a model name such as Prosit_2020_intensity_HCD. Use describe_predictor_model before a new Koina model: it retrieves that model's live input/output schema and reports whether the currently supported input fields are sufficient.""",
    "method://terminal_ion_annotation_v1": """proteomics-mcp terminal-ion annotation is used as an independent observed-spectrum annotation summary in the scoring path and as a standalone MCP tool.

The primary match schema remains proteomics-mcp MatchedPeak, based on predicted-fragment to observed-peak matching. The internal annotator contributes annotation_summary fields including annotated peak count, annotated observed-intensity fraction, ion-type counts, per-peak labels, and detailed observed/theoretical m/z error evidence.

Supported annotation ion families include a, a+1, b, c, d, d-alt, x, y, z, z-dot, w, and w-alt plus auxiliary ion classes immonium, reporter, and precursor-related peaks. Neutral-loss variants can be generated when annotation_neutral_losses=true. Use ion_types='all' for standard, exotic, and auxiliary ions, or ion_types='exotic' for only a+1/d/d-alt/z-dot/w/w-alt.""",
    "method://sdrf_support_v1": """SDRF support in proteomics-mcp.

The MCP exposes:
- load_sdrf: parse a local SDRF-Proteomics TSV and return grouped sample/run metadata
- retrieve_pride_sdrf: retrieve an SDRF file from a PRIDE project accession through the PRIDE Archive API and summarize it

The returned summary includes:
- data-file columns and factor columns
- assay, organism, label, and instrument summaries
- file_metadata keyed by raw filename, basename, and filename stem for multi-run joins and metadata enrichment

This is intentionally optional. Users can score spectra without SDRF, but when SDRF exists it can provide sample and acquisition context through MCP-accessible metadata.""",
    "method://visualization_modes_v1": """Visualization modes:
- render_spectrum_comparison: mirror plot with observed/predicted b/y fragment evidence and score metrics. The HTML includes an optional checkbox to show or hide fragment-ion subsequences alongside ion labels.
- render_annotated_observed_spectrum: observed MS2 spectrum annotated with the internal terminal-ion annotator; default ion_types='by' and max_isotope=2, optionally pass ion_types='abcxyz', ion_types='all', or ion_types='exotic'. Auxiliary ion classes such as immonium, reporter, and precursor can also be requested. No predicted spectrum is shown. The HTML includes an optional checkbox to show or hide fragment-ion subsequences alongside ion labels.
- render_observed_spectrum: observed MS2 spectrum only with no annotations or predicted spectrum; includes metrics derivable from that spectrum, including peak count, TIC, base peak, m/z range, intensity summaries, precursor-region intensity fraction, and acquisition metadata.
- render_static_spectrum: write a publication-ready mirror, annotated, or plain spectrum from the same normalized payload. The output suffix selects SVG, PDF, PNG, or another Matplotlib-supported format; SVG and PDF are native vector files rather than browser screenshots.
- render_spectra_batch: batch rendering wrapper over the same three spectrum views. It accepts either inline `items` or a file-backed `input_path`, writes category-specific output directories (`comparison/`, `annotated/`, `observed/`) under one output directory, supports `validate_only=true` for schema checking without rendering, and runs in parallel within each selected category.""",
    "method://dataset_summary_v1": """Dataset-level aggregation over scored PSM tables.

The MCP exposes:
- summarize_dataset: aggregate a scored PSM table into dataset, peptide, peptidoform, run, and optional SDRF-linked sample summaries
- render_dataset_report: generate a static HTML dataset report
- list_dataset_groups: list dataset, run, and sample grouping keys
- get_peptide_summary / get_peptidoform_summary: return one aggregated row by sequence key

Dataset summaries are built on top of score_candidates output and preserve score, confidence, run, charge, modification, and spectrum-evidence fields. Peptide tables collapse modified peptidoforms while retaining distinct peptidoform counts and representative modified-form metadata.

This is a scored-table aggregation path. For the primary ID-table workflow with optional scored enrichment, use method://id_table_inspection_v1 instead.""",
    "method://id_table_inspection_v1": """ID-table inspection and optional scored enrichment.

The MCP exposes:
- inspect_psm_table: normalize a raw identification/candidate table into the canonical PSM schema and return an inspectable PSM preview plus grouped metadata
- inspect_peptide_table: return the collapsed peptide table or one peptide row by stripped sequence
- list_distribution_fields: list canonical and upstream source_columns fields that can be summarized numerically or categorically
- summarize_numeric_field: return numeric summary statistics for one field, optionally grouped by another field and optionally weighted by a numeric weight field
- summarize_categorical_field: return frequency counts and fractions for a categorical field, optionally grouped by another field
- summarize_missing_values: return missing counts and missing fractions for selected fields
- summarize_cv_by_category: return mean, standard deviation, and coefficient of variation for one numeric field across categories
- bin_numeric_field: return histogram-like bin counts for one numeric field, optionally grouped by another field and optionally weighted by a numeric weight field
- bin_2d_numeric_fields: return a 2D bin grid for two numeric fields, optionally weighted by a numeric weight field
- summarize_dashboard: return a dashboard-oriented bundle of overview counts, missingness, histograms, 2D grids, and score-aware scatter panels when scored columns are present
- summarize_overlap: return overlap counts and pairwise overlap metrics across runs/files or other grouping fields for peptides, peptidoforms, or precursor-level entities
- summarize_modifications: return modification counts, site counts, residue counts, and per-group modification summaries
- summarize_quant_quality: return quant completeness and entity-level CV summaries across groups for a quant field
- summarize_file_coverage: return per-file/run summaries including counts, missingness, charge counts, histograms, and m/z-vs-RT coverage grids
- render_psm_table: render the normalized PSM table as browsable HTML with in-page CSV download
- render_peptide_table: render the collapsed peptide table as browsable HTML with in-page CSV download
- render_summary_html: render any structured summary payload as simple HTML for inspection
- render_spectra_batch: render many spectrum HTML artifacts at once from precomputed payloads, optionally for only selected visualization categories; supports file-backed `input_path` and validation-only schema checks for large local batches

This path is ID-table-first. It can operate in two modes:
- ID-only: use only the fields already present upstream, such as peptidoform, run, scan, precursor charge, confidence, next beam, and token-confidence summaries
- ID + scored enrichment: merge a scored PSM table back onto the same canonical rows via `candidate_id`, so inspection tables, summaries, and dashboards also expose spectrum-evidence fields such as spectral angle, matched-ion fraction, explained intensity, log2_xcorr, mass-error summaries, and acquisition metadata

Supported input adapters for this path currently include instanovo, maxquant, fragpipe, proteome_discoverer, and generic mapped tables.

Dashboard behavior in this path:
- ID-only dashboards from normalized candidate tables
- metadata-enriched dashboards after Feather joins derived from mzML, MGF, or Bruker `.d`
- score-aware dashboards when a scored PSM table is merged with `scored_path`
- light/dark/auto theme support in rendered HTML
- run-specific dropdown panels for per-run histograms and scatter plots

Current score-aware dashboard panels include score distributions for spectral_angle, matched_ion_fraction, explained_intensity_fraction, and log2_xcorr, plus scatter relationships for confidence vs spectral_angle, spectral_angle vs pearson_correlation, spectral_angle vs log2_xcorr, confidence vs matched_ion_count, spectral_angle vs precursor_mass_error_ppm, and spectral_angle vs fragment_mass_error_abs_mean_ppm.

Score-aware dashboards are typically produced through build-prediction-batch -> score-from-prediction-batch -> build-id-dashboard --scored-path ... .""",
    "method://spectrum_conversion_v1": """Spectrum conversion utilities.

The MCP exposes:
- mzml_to_feather: convert mzML spectra to Feather through pyOpenMS
- mgf_to_feather: convert MGF spectra to Feather through pyteomics
- bruker_d_to_feather: convert a Bruker `.d` directory to Feather, preferring AlphaTims direct extraction and otherwise using an explicit `ProteinAnalysisResults.mgf` fallback when that export is present

All three converters produce a semi-wide Feather table with one row per spectrum and fields such as source_file/source_name, spectrum_id, scan_number, precursor_mz, precursor_charge, mz_array, intensity_array, retention time when available, and converter-specific metadata.""",
}
