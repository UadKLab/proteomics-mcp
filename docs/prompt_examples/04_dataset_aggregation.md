# Dataset Aggregation Prompts

## Example 1: Summarize a scored table

**User prompt**

Summarize this scored PSM table into dataset, peptide, peptidoform, and run-level outputs.

**What the MCP should do**

- call `summarize_dataset`

**Expected output shape**

- `dataset_metadata`
- `identification_summary`
- `sequence_summary`
- `run_summary`
- `confidence_summary`
- preview rows for peptide and peptidoform tables

## Example 2: Retrieve one peptide summary

**User prompt**

Show me the peptide-level summary row for `LQLETELEALKEELLFMK`.

**What the MCP should do**

- call `get_peptide_summary`

**Expected output shape**

- one peptide-level aggregate row

## Example 3: Retrieve one peptidoform summary

**User prompt**

Return the best peptidoform-level row for `M[UNIMOD:35]PEPTIDE`.

**What the MCP should do**

- call `get_peptidoform_summary`

**Expected output shape**

- one modification-aware aggregate row

## Example 4: Inspect grouping structure

**User prompt**

List the dataset groups available in this scored table and tell me whether sample-level grouping is available.

**What the MCP should do**

- call `list_dataset_groups`

**Expected output shape**

- dataset/run/sample grouping keys
- indication of whether SDRF-linked sample summaries exist
