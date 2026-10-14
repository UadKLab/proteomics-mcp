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

- `instanovo`: legacy `preds`/`log_probs`, beam outputs, and current `instanovo_*`/`diffusion_*` variants.
- `casanovo`: Casanovo prediction CSVs.
- `maxquant`: `msms.txt` exports.
- `fragpipe`: PSM tables including run/scan parsing from spectrum fields.
- `proteome_discoverer`: PSM exports, including `.xlsx`.
- `generic`: CSV, TSV, or XLSX using explicit column mappings.

Canonical PSM fields include identity/provenance, peptidoform, precursor values, native score and optional 0-1 confidence, de novo beam confidence when genuinely available, token evidence, and preserved upstream `source_columns`. Modifications are converted to inline UNIMOD notation. For multi-run data, match spectra using raw-file stem plus scan number rather than scan number alone.

### Generic Mapping

Required: a peptide/peptidoform and a spectrum locator (`usi`, `spectrum_id`, or raw-file stem plus scan number). Strongly recommended: precursor m/z/charge, score semantics, and token-level evidence where available.

## Spectrum Sources And Conversion

- IPC/Feather: normalized local spectrum tables.
- mzML: `mzml-to-feather`, backed by pyOpenMS.
- MGF: `mgf-to-feather`, backed by pyteomics.
- Bruker `.d`: `bruker-d-to-feather`, preferring AlphaTims with an embedded-MGF fallback.
- USI: `retrieve-usi-spectrum` through a configured provider.
- InstaNovo export: `export-to-instanovo` writes a normalized IPC or Parquet inference table.

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

`build_prediction_batch` writes reusable prediction, fragment, mapping, failure, log, status, and manifest artifacts. MS2PIP and Koina batches isolate malformed requests recursively so valid rows still complete; `score_from_prediction_batch` reuses predictions against local spectra and creates analogous scored artifacts. Use the corresponding `inspect_*_status` tools to report durable progress during long operations.

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
