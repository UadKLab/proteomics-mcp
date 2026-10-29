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
- xcorr_like, log10_xcorr_like, log2_xcorr_like
- hyperscore_like, log2_hyperscore_like
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
- precursor_isotope_error
- precursor_mass_error_corrected_ppm
- precursor_mass_error_corrected_da
- engine_precursor_isotope_error
- engine_precursor_mass_error_ppm

precursor_mass_error_* is the observed minus the theoretical monoisotopic precursor m/z. The theoretical mass includes
every residue and terminal modification; it is null when a modification cannot be resolved. precursor_isotope_error
(0 to 3) is the 13C isotope peak closest to the observed m/z, and precursor_mass_error_corrected_* is the error
relative to that peak. The observed m/z is the spectrum's precursor m/z, for every engine; the
candidate's is used only when the spectrum has none. engine_precursor_isotope_error and engine_precursor_mass_error_ppm
are the same from the candidate's precursor m/z as the identification table reports it, recalibrated by MaxQuant,
FragPipe and Proteome Discoverer: the engine's own error.

The profile preserves the sidecar-style formulas where practical while normalizing all outputs into SpectrumMetricResult.

spectral_angle, cosine_similarity and pearson_correlation compare observed and predicted intensities over the matched
predicted ions only, so a candidate matching a single peak scores 1.0. Each predicted ion takes its nearest peak, even one
another ion also takes, so a shared peak counts more than once in matched ions, explained intensity and the hyperscore.
xcorr_like and hyperscore_like are the 0.1.x variants (named xcorr, log2_xcorr, hyperscore, log2_hyperscore there):
xcorr_like bins observed and predicted intensities at 0.1 m/z and subtracts the mean over shifts of +-75 bins;
hyperscore_like is ln(n!) + ln(sum of matched observed) + ln(sum of matched predicted). They are not the engines'
scores; xcorr and hyperscore are null under this profile.
Kept to reproduce scores from 0.1.10 and earlier; the default profile is v2 (method://metric_profile_v2).""",
    "method://metric_profile_v2": """Metric profile v2 (default).

The metrics of legacy_v1 (method://metric_profile_legacy_v1), except that xcorr and hyperscore are the search
engines' scores (#51), computed on the candidate's theoretical b/y ions, independent of the predictor:
- xcorr: Comet 2026.01.1 XCorr (fast-XCorr preprocessing: sqrt intensities, 10 windows normalised to 50, the mean over
  +-75 bins subtracted; each theoretical bin counted once; x0.005, 3 decimals). xcorr_preset "high_res" (default:
  fragment_bin_tol 0.02, offset 0.0, flanking peaks) or "low_res" (1.0005, 0.4, M peak only) for ion-trap spectra.
- hyperscore: the X! Tandem 2015.12.15.2 hyperscore as X! Tandem reports it, 4*log10(H), with
  H = (sum of observed intensity x ion weight) x the two largest per-charge N_b! N_y!, on X! Tandem's conditioned
  spectrum (isotope clean-up, parent region removed, m/z > 150, dynamic range 100, top 50 peaks; null with fewer than
  15 peaks). hyperscore_fragment_error / _unit (default 0.4 Da, X! Tandem's default).
validation/engine_scores/ in the repository reproduces both against the engines themselves. They are not the
raw_score of a FragPipe (MSFragger hyperscore, a different definition) or Proteome Discoverer (Sequest HT XCorr) row.

spectral_angle, cosine_similarity and
pearson_correlation compare observed and predicted intensities over every predicted ion, with zero observed intensity
for predicted ions that matched no peak, so unexplained predicted signal lowers the score. When no predicted ion
matches, spectral_angle and cosine_similarity are 0.0 and pearson_correlation is null.

Matching is one to one: each observed peak explains at most one predicted ion, assigned closest first, so a peak shared
by isobaric ions (e.g. b3/y2) counts once and explained_intensity_fraction is at most 1.

Predicted-ion cut-off: a predicted ion is scored only if its intensity is positive and at least min_predicted_intensity
(default 0.001) times the most intense predicted ion. Predictors give ions they expect to be absent a tiny intensity or
0, and those would count in total_predicted_ions and could be matched to noise. total_predicted_ions,
matched_ion_count and matched_ion_fraction therefore refer to the ions above the cut-off. Set
min_predicted_intensity = 0 to score every ion with a positive intensity; legacy_v1 applies no cut-off.

Select a profile with metric_profile in a JSON or TOML config file passed as config_path (CLI: --config).""",
    "method://candidate_schema_v1": "Canonical candidate schema normalizes upstream engines into peptide/peptidoform, scan/spectrum identifiers, raw score, optional 0-1 confidence, and token-confidence summaries. Supported adapters include instanovo, casanovo, maxquant, fragpipe, proteome_discoverer, and generic mapped tables.",
    "schema://candidate_psm_v1": """CandidatePSM fields:
candidate_id, source_engine, source_engine_version, experiment_name, scan_number, spectrum_id, usi, peptidoform, precursor_mz, precursor_charge, rank, raw_score, raw_score_name, next_beam, next_beam_margin, confidence, token_log_probabilities, token_log_probability_mean, token_log_probability_min, token_confidences, token_confidence_mean, token_confidence_min, token_confidence_min_position, token_confidence_min_token, source_columns, warnings.

`warnings` lists problems found while reading the row that may make its peptidoform or scores wrong, e.g. MaxQuant fixed modifications that could not be applied.

For InstaNovo beam outputs, `next_beam` is the beam-1 confidence on a 0-1 scale and `next_beam_margin` is `confidence - next_beam`. Search-engine score deltas remain in source_columns rather than being presented as beams.""",
    "schema://observed_spectrum_v1": """ObservedSpectrum fields:
spectrum_id, scan_number, usi, precursor_mz, precursor_charge, mz, intensity, retention_time, source, metadata.

mz and intensity are parallel numeric arrays. IPC and Feather inputs are normalized into this schema. Feather inputs can be produced from mzML, MGF, or Bruker `.d` sources.""",
    "schema://predicted_spectrum_v1": """PredictedSpectrum fields:
peptidoform, precursor_charge, fragments, predictor_name, predictor_version, model, metadata.

Each fragment contains annotation, ion_type, ordinal, charge, mz, intensity, neutral_loss, isotope, metadata.""",
    "schema://spectrum_metric_result_v1": """SpectrumMetricResult fields:
candidate, spectral_angle, cosine_similarity, pearson_correlation, xcorr, log10_xcorr, log2_xcorr, hyperscore, log2_hyperscore, matched_ion_count, total_predicted_ions, matched_ion_fraction, explained_intensity_fraction, precursor_signal_fraction, fragment and precursor mass-error summaries, matched_peaks, metric_profile, predictor_metadata, scored_ions, excluded_predicted_ions, metric_n, scoring_settings, annotation_settings, annotation_summary, acquisition_provenance.

What was scored: scored_ions lists every scored predicted ion with its predicted intensity (divided by the most intense scored ion) and the intensity of the observed peak it matched, 0 if none; these are the vectors spectral angle, cosine and Pearson compare. excluded_predicted_ions lists the predicted ions left out, each with its reason (fragment_charge_above_limit, invalid_mz, non_finite_intensity, non_positive_intensity, below_min_predicted_intensity). metric_n gives the number of values behind each metric. scoring_settings records the predictor, model version, collision energy sent and where it came from, tolerance and unit, metric profile, min_predicted_intensity, the fragment-charge rule and whether the prediction was scored as given.

A prediction passed in (score_spectrum's predicted, or a prediction batch) goes through the same rules as one the package makes: under metric profile v2, fragments above max_fragment_charge, or above precursor charge - 1 when it is not set, are left out. Set score_predictions_as_given to score every fragment charge; the result records it. Batch scoring writes the per-ion detail to <prefix>_scored_ions.jsonl (score_candidates: <output>_scored_ions.jsonl) and the settings to the manifest.

MCP score_spectrum also includes metric_summary, an agent-oriented grouping with match_quality, predicted_match_evidence, observed_annotation_evidence, mass_error, candidate, and acquisition fields. A backward-compatible ion_evidence group is retained as a merged view. Acquisition summary includes fragmentation_method, isolation_target_mz, retention_time_seconds/minutes, collision_energy, FAIMS compensation voltage, injection time, scan_header, scan_resolution, instrument_model, and detector_type when available from the scan header or fallbacks.""",
    "config://proteomics_mcp_v1": """Config file keys, JSON or TOML:
predictor_backend, predictor_model, ms2pip_model (legacy alias), max_fragment_charge, koina_server_url, koina_timeout_seconds, koina_collision_energy, koina_fragmentation_type, koina_instrument_type, koina_batch_size, koina_cache_dir, tolerance, tolerance_unit, usi_backend, metric_profile, min_predicted_intensity, score_predictions_as_given, annotation_ion_types, annotation_max_isotope, annotation_max_ion_charge, annotation_neutral_losses.

Koina collision energy: a collision_energy tool argument, else the spectrum's own normalized collision energy when it states one (Thermo scan filter @hcd28.00 or a Thermo collision energy; never eV values from other vendors), else koina_collision_energy (30). The value used and its source are in predictor_metadata.

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
Keys: peptide_column (required), scan_column, spectrum_id_column, run_column (alias experiment_name_column: the PSM's run, for multi-run spectra tables), usi_column, charge_column, precursor_mz_column, score_column, next_beam_column, next_beam_margin_column, next_score_column, next_score_margin_column, token_log_probabilities_column, score_type (confidence_0_1, posterior_error_probability, evalue, score_unscaled, unknown; default unknown), source_engine (default "generic").

For multi-run local data, prefer joining by experiment_name/raw-file stem plus scan_number rather than scan_number alone.""",
    "method://supported_predictors": """Supported predictors:
- ms2pip: local MS2PIP fragment-intensity prediction.
- koina: remote prediction through a Koina server, including Prosit-family models selected by predictor_model.
- mock: deterministic lightweight predictor for tests and offline workflow checks.

These predict fragment intensities for a given peptidoform. The de novo sequencing and database search engines (InstaNovo, Casanovo, MaxQuant, FragPipe, Proteome Discoverer) are upstream inputs: their result tables are read through adapters, and none of them is run by this package. The canonical public Koina/Prosit path is backend='koina' plus a model name such as Prosit_2020_intensity_HCD. Use describe_predictor_model before a new Koina model: it retrieves that model's live input/output schema and reports whether the currently supported input fields are sufficient.""",
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
- multiplexed runs (TMT, iTRAQ, SILAC; one SDRF row per channel): file_metadata keeps every channel in "channels"
  (channel_count); fields that differ between channels hold their distinct values joined with "; ", e.g. label
  "TMT126; TMT127N", because a PSM of such a run belongs to all its channels

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
- thermo_raw_to_mzml: convert a Thermo `.raw` file to mzML with ThermoRawFileParser, which the user installs separately (on PATH or in $PROTEOMICS_MCP_THERMORAWFILEPARSER); no Thermo libraries are bundled. The result records the converter's version and command
- thermo_raw_to_feather: the same, then mzml_to_feather on the mzML, which is kept
- mgf_to_feather: convert MGF spectra to Feather through pyteomics
- bruker_d_to_feather: convert a Bruker `.d` directory to Feather, preferring AlphaTims direct extraction and otherwise using an explicit `ProteinAnalysisResults.mgf` fallback when that export is present

The Feather converters produce a semi-wide table with one row per spectrum and fields such as source_file/source_name, spectrum_id, scan_number, precursor_mz, precursor_charge, mz_array, intensity_array, retention time when available, and converter-specific metadata.""",
}
