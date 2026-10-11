# Batch Prediction And Scoring Prompts

## Example 1: Build and score a constrained batch

**User prompt**

Build a prediction batch for the top 100 confident candidates from this table, then score them against the Feather spectra file.

**What the MCP should do**

- call `build_prediction_batch`
- call `inspect_prediction_batch_status` until complete if needed
- call `score_from_prediction_batch`
- call `inspect_scored_batch_status` until complete if needed

**Expected output shape**

- prediction manifest path
- scoring manifest path
- scored CSV path
- failure/status artifacts when present

## Example 2: Poll a running batch

**User prompt**

Poll the prediction batch and tell me whether it is still running, how many items are done, and whether anything has failed.

**What the MCP should do**

- call `inspect_prediction_batch_status`

**Expected output shape**

- `phase`
- `completed`
- `total`
- `succeeded`
- `failed`
- `progress_fraction`

## Example 3: Inspect batch failures

**User prompt**

Inspect the scored batch failures and group them by error type.

**What the MCP should do**

- call `inspect_scored_batch`
- summarize `failure_type_counts`

**Expected output shape**

- counts per error type
- preview rows for failures

## Example 4: Render the best scored subset

**User prompt**

Render comparison plots for the ten best-scoring PSMs from the scored batch.

**What the MCP should do**

- inspect the scored batch output
- select the top-scoring rows
- prepare payloads
- call `render_spectra_batch`

**Expected output shape**

- HTML files grouped by rendering mode
- batch render status path
