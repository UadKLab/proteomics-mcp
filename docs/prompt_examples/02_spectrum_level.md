# Spectrum-Level Prompts

## Example 1: Score one selected PSM

**User prompt**

Take the highest-confidence PSM for `PEPTIDESEQ`, load the matching spectrum, predict fragments with MS2PIP, and score it.

**What the MCP should do**

- inspect the table first to identify the target PSM
- load or retrieve the observed spectrum
- call `score_spectrum`
- optionally call `summarize_score_metrics`

**Expected output shape**

- one `SpectrumMetricResult`
- grouped metrics such as `match_quality`, `predicted_match_evidence`, `mass_error`, and `acquisition`

## Example 2: Render a mirror plot

**User prompt**

Render a mirror plot for scan 12345 from run `sample_A` using the candidate sequence from this table.

**What the MCP should do**

- identify the matching candidate row
- score the spectrum or build a payload from observed and predicted spectra
- call `render_spectrum_comparison`

**Expected output shape**

- standalone HTML artifact
- path to the rendered file

## Example 3: Annotate observed ions only

**User prompt**

Annotate the observed spectrum for this phosphopeptide using all ion types, including neutral losses, without showing the predicted spectrum.

**What the MCP should do**

- call `annotate_observed_spectrum` or `render_annotated_observed_spectrum`
- use `ion_types="all"`
- use `neutral_losses=true`

**Expected output shape**

- observed-spectrum annotation summary
- optionally an annotated observed-only HTML view

## Example 4: Compare several selected PSMs

**User prompt**

Compare the spectral angle, matched-ion fraction, explained intensity, and XCorr for these five selected PSMs.

**What the MCP should do**

- score each selected target
- summarize the returned metric groups across the subset

**Expected output shape**

- one short comparison table or structured list of per-PSM metrics
- no need to render HTML unless requested
