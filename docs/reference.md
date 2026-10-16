# Technical Reference

## Operating Model

Use the package in two stages:

1. Normalize and inspect an upstream identification table.
2. Select one candidate or a constrained subset for spectrum retrieval, prediction, scoring, annotation, and rendering.

This keeps table exploration responsive and makes expensive scoring explicit. Batch workflows are available when many selected spectra need prediction or scoring.

## Installation

The project uses Python 3.11+ and `uv`.

```bash
uv sync
uv run pytest -q
```

For an editable install outside `uv`:

```bash
pip install -e ".[dev]"
```

## MCP Setup

`proteomics-mcp serve` is a stdio MCP server. An MCP client launches it, sends JSON-RPC messages, and makes its tools available to an LLM. It is not an interactive terminal application.

Use this configuration pattern, replacing `cwd` with your checkout:

```json
{
  "mcpServers": {
    "proteomics-mcp": {
      "command": "uv",
      "args": ["run", "proteomics-mcp", "serve"],
      "cwd": "/absolute/path/to/proteomics-mcp"
    }
  }
}
```

For HPC use, either run upstream inference/search on the cluster and inspect mounted outputs locally, or run the MCP process near the shared data and connect through a client that can launch that process.

## Input Tables

Supported adapters:

- `instanovo`: prediction CSVs of InstaNovo 1.0–1.2: transformer greedy and beam search (`preds`/`log_probs`, `predictions`/`log_probabilities`, `predictions`/`log_probs`), InstaNovo+ with and without refinement (`final_prediction`, `diffusion_predictions`, falling back to `transformer_predictions` where refinement produced nothing), verbose prediction columns, and combined runs. Notes:
  - For MGF input (and mzML input in InstaNovo 1.0), InstaNovo's `scan_number` column holds the 0-based position of the spectrum in the file. The adapter recognises this (the column contains 0) and stores it as `spectrum_index`, which spectrum lookup uses by position within the run, instead of as a scan number. Convert the same MGF or mzML file with `mgf_to_feather` / `mzml_to_feather`.
  - Confidence is the exponent of the sequence log-probability.
  - Token log-probabilities of legacy beam-search files are aligned to the sequence.
- `casanovo`: Casanovo prediction CSVs.
- `maxquant`: `msms.txt` exports.
- `fragpipe`: PSM tables including run/scan parsing from spectrum fields.
- `proteome_discoverer`: PSM exports, including `.xlsx`.
- `generic`: CSV, TSV, or XLSX using explicit column mappings. Pass a JSON mapping file with `--column-mapping` (CLI) or `column_mapping_path` (MCP), on any command or tool that takes an adapter:

  ```json
  {"peptide_column": "Sequence", "scan_column": "Scan", "charge_column": "Charge",
   "score_column": "Score", "score_type": "confidence_0_1"}
  ```

  `peptide_column` is required. The other keys are `scan_column`, `spectrum_id_column`, `usi_column`, `charge_column`, `precursor_mz_column`, `score_column`, `next_beam_column`, `next_beam_margin_column`, `next_score_column`, `next_score_margin_column`, `token_log_probabilities_column`, `score_type` (`confidence_0_1`, `posterior_error_probability`, `evalue`, `score_unscaled` or `unknown`, the default) and `source_engine` (default `generic`). The mapping is rejected for the named adapters.

Canonical PSM fields include identity/provenance, peptidoform, precursor values, native score and optional 0-1 confidence, de novo beam confidence when genuinely available, token evidence, and preserved upstream `source_columns`. Modifications are converted to inline UNIMOD notation. For multi-run data, match spectra using raw-file stem plus scan number rather than scan number alone.

Spectrum lookup (`load_spectrum_from_table`, and the scoring tools) never silently substitutes another spectrum:
- A run name that matches none of the runs in a multi-run spectrum table is an error. In a single-run table, a different name (for example a renamed file) still refers to that run.
- A scan number or spectrum id that matches several *different* spectra, such as the same scan in several runs when no run name is given, is an error. Identical copies, as written by MGF charge-state enumeration, are not.
- Scan numbers are read only from identifiers that state one: `scan=` in native IDs and USIs, Thermo-style MGF titles (`<run>.<scan>.<scan>.<charge>`, also `.dta`), `index=`, `<run>:<scan>`, or a plain number. SCIEX (`sample=… cycle=… experiment=…`) and Bruker (`merged=… frame=…`) native IDs have no single scan number. Their spectra are found by spectrum id instead of by an unrelated number from the identifier.

Recognised modification notations include:

| Source | Examples |
|---|---|
| InstaNovo | `M(+15.99)`, `N(+.98)`, `(+42.01)PEPTIDE`, `[UNIMOD:35]` |
| Casanovo | `+42.011PEPM+15.995K` |
| MaxQuant `Modified sequence` | `_(ac)PEPM(ox)K_`, `_(Acetyl (Protein N-term))PEPM(Oxidation (M))K_` |
| FragPipe | `Assigned Modifications` such as `N-term(42.0106), 5M(15.9949)`; `Modified Peptide` with total masses (`n[43]PEPTM[147]K`) or deltas |
| Proteome Discoverer `Modifications` | `M1(Oxidation); N-Term(Acetyl)`, `1xTMT6plex [N-Term]; 2xCarbamidomethyl [C5; C12]` |

Common modifications (acetyl, carbamidomethyl, carbamyl, deamidation, phospho, oxidation, pyro-Glu, ammonia loss, amidation, methyl, dimethyl, formyl, GG, iTRAQ4plex, TMT6plex/10plex/11plex, TMTpro) resolve without network access; other names and masses are looked up in the Unimod catalogue. A modification that cannot be resolved or placed, such as an unlocalised site (`1xPhospho [S/T]`), stays in the peptidoform as an explicit token (unlocalised ones as ProForma `[UNIMOD:21]?PEPSTK`). It is flagged as an unsupported modification rather than dropped, so the peptidoform is never silently scored as unmodified. MaxQuant `Modified sequence` omits fixed modifications (#101).

### Generic Mapping

Required: a peptide/peptidoform and a spectrum locator (`usi`, `spectrum_id`, or raw-file stem plus scan number). Strongly recommended: precursor m/z/charge, score semantics, and token-level evidence where available.

## Spectrum Sources And Conversion

- IPC/Feather: normalized local spectrum tables.
- mzML: `mzml-to-feather`, backed by pyOpenMS.
- MGF: `mgf-to-feather`, backed by pyteomics.
- Bruker `.d`: `bruker-d-to-feather`, preferring AlphaTims with an embedded-MGF fallback. For ddaPASEF, `scan_number` is the precursor `Id`, the identifier of each MS2 spectrum (also in `precursor_id` and `spectrum_id`, `<run>:precursor:<Id>`). `mobility` is the reduced ion mobility 1/K0 (V·s/cm²), and `scan_index` is the precursor's mean position on the TIMS scan axis. Precursors without an assigned charge are kept, with an empty `precursor_charge`.
- USI: `retrieve-usi-spectrum` through a configured provider.
- InstaNovo export: `export-to-instanovo` writes a normalized IPC or Parquet inference table.

USI spectra come from the provider's own processing of the raw file, so they can differ from a local mzML or MGF conversion of the same scan. For the 25 reference PSMs in `examples/` (PXD044934, PXD014877, PXD010595), PRIDE returned the same peaks as the local conversions for 2 PSMs. For the other 23 it returned 1–33 fewer peaks; the missing peaks were the only difference, and every remaining peak had the same m/z and intensity. Scores computed from a USI spectrum can therefore differ from scores computed from a local file. Record which source each result used.

## Fragment Prediction

Prediction backends:

- `ms2pip`: local MS2PIP; the default. It runs on your computer; its model files are downloaded once on first use (from `genesis.ugent.be`, with a Zenodo fallback) and cached in `~/.ms2pip`, after which predictions run offline. On macOS it needs the OpenMP runtime (`brew install libomp`). Single-spectrum requests use a bounded process cache, while `build_prediction_batch` groups compatible requests through MS2PIP's native batch API.
- `koina`: remote model inference, including `Prosit_2020_intensity_HCD`.
- `mock`: deterministic test backend; never use it for scientific interpretation.

For Koina, call `describe-predictor-model` before a new model or large batch. The wrapper retrieves live model metadata, submits only supported required fields, caches completed predictions locally, and rejects incompatible models before inference. By default, retained fragment charges are limited to `<= precursor charge - 1`.

```toml
predictor_backend = "koina"
predictor_model = "Prosit_2020_intensity_HCD"
koina_server_url = "https://koina.wilhelmlab.org"
koina_collision_energy = 30.0
koina_batch_size = 64
koina_cache_dir = ".proteomics_mcp_cache/koina"
```

JSON and TOML configuration files are accepted. Explicit tool arguments override configuration values.

## Spectrum-Level Evidence

`score_spectrum` and `score_candidates` provide candidate-to-observed-spectrum evidence. The default tolerance is `0.05 Da`.

Metrics include spectral angle, cosine similarity, Pearson correlation, XCorr/log2 XCorr, Hyperscore/log2 Hyperscore, matched ion count/fraction, explained intensity fraction, precursor signal fraction, observed annotation counts/intensity, and precursor/fragment mass error in Da and ppm.

The precursor mass error (`precursor_mass_error_ppm`, `precursor_mass_error_da`) is the observed minus the theoretical monoisotopic precursor m/z. The theoretical mass includes every residue and terminal modification; when a modification cannot be resolved, the error is left empty rather than computed without it. Instruments sometimes select a 13C isotope peak as the precursor, which shows up as an error of about 1.003/z m/z. `precursor_isotope_error` (0, 1 or 2) gives the isotope peak closest to the observed m/z, and `precursor_mass_error_corrected_ppm` / `precursor_mass_error_corrected_da` the error relative to that peak.

Observed-spectrum annotation is implemented internally. It supports standard terminal-ion families, exotic families, immonium ions, reporter ions, precursor-related peaks, isotope variants, and optional neutral losses. Use `ion_types=by` for b/y ions, `all` for the broad set, or `exotic` for exotic terminal-ion families.

## Rendering

- `render_spectrum_comparison`: observed-vs-predicted mirror plot with scoring evidence.
- `render_annotated_observed_spectrum`: annotated observed MS2 spectrum without a predicted spectrum.
- `render_observed_spectrum`: observed spectrum without annotations.
- `render_static_spectrum`: Matplotlib/Seaborn static mirror, annotated, or plain spectrum figure. The output suffix selects SVG, PDF, PNG, or another Matplotlib format. SVG and PDF are native vector output, suitable for figure assembly.

Outputs are Plotly HTML with hover data, zoom, image export, light/dark support, and the underlying structured payload. Interactive views load Plotly from a CDN when opened. The optional fragment-ion sequence control uses compact display notation for modifications, for example `PEPM(ox)`, while all internal scoring and prediction retain normalized UNIMOD notation. `render_spectra_batch` renders many precomputed payloads in parallel; use `input_path` and `validate_only=true` for realistic local batches. `input_path` may be JSON (a list of items, or an object with `items` and optional `views`), JSON Lines (`.jsonl` or `.ndjson`, one item per line), or either of them gzip-compressed (`.gz`).

## Tables, Summaries, And SDRF

`inspect_psm_table` and `inspect_peptide_table` normalize and inspect upstream results without rescoring the dataset. Peptide tables collapse modified forms while retaining peptidoform context, run support, charge species, and representative PSM evidence.

Structured tools support field discovery, numeric/categorical summaries, binning, missingness, overlap, modification summaries, quant quality, file coverage, and dashboard bundles. `load_sdrf` and `retrieve_pride_sdrf` add optional run/sample context. PRIDE retrieval recognizes `.sdrf`, SDRF-labelled TSV/TXT, and gzip-compressed variants; API failures report the project accession, request URL, and HTTP response detail.

## Batch Workflows

`build_prediction_batch` writes reusable prediction, fragment, mapping, failure, log, status, and manifest artifacts. MS2PIP and Koina batches isolate malformed requests recursively so valid rows still complete; `score_from_prediction_batch` reuses predictions against local spectra and creates analogous scored artifacts. Use the corresponding `inspect_*_status` tools to report durable progress during long operations: the manifest is written when a batch starts (`phase: running`) and completed at the end, so its status can be read while the batch runs and after an interruption (`seconds_since_update` shows when the batch last reported). Results are written when a batch finishes; an interrupted batch has to be run again.

The direct `score_candidates` path is suitable for smaller selected sets.

## Tool Groups

- Table inspection: `load_candidates`, `inspect_*`, `render_*_table`.
- Structured summaries: `summarize_*`, `bin_*`, `build_id_dashboard`.
- Spectrum conversion/retrieval: `mzml_to_feather`, `mgf_to_feather`, `bruker_d_to_feather`, `export_to_instanovo`, `retrieve_usi_spectrum`.
- Prediction/evidence: `describe_predictor_model`, `predict_fragments`, `annotate_observed_spectrum`, `score_spectrum`, `score_candidates`, `summarize_score_metrics`.
- Batch orchestration: `build_prediction_batch`, `inspect_prediction_batch*`, `score_from_prediction_batch`, `inspect_scored_batch*`.
- Visualization: `render_spectrum_comparison`, `render_annotated_observed_spectrum`, `render_observed_spectrum`, `render_static_spectrum`, `render_spectra_batch`.

## Validation

```bash
uv run pytest -q
uv run proteomics-mcp --help
```
