# Table Inspection Prompts

## Example 1: Inspect an upstream table

**User prompt**

Inspect this FragPipe PSM table and tell me how many PSMs, peptides, modified peptidoforms, and runs it contains.

**What the MCP should do**

- call `inspect_psm_table`
- use `adapter="fragpipe"`
- summarize `dataset_metadata`, `sequence_summary`, and `run_summary`

**Expected output shape**

- structured counts for PSMs, peptides, peptidoforms, and runs
- preview rows from the normalized PSM table
- no rescoring

## Example 2: Inspect one peptide

**User prompt**

Show me the peptide summary for `LQLETELEALKEELLFMK` and list the supporting PSMs.

**What the MCP should do**

- call `inspect_peptide_table` with `peptide="LQLETELEALKEELLFMK"`
- optionally call `inspect_psm_table` to retrieve preview rows or supporting candidate IDs

**Expected output shape**

- one collapsed peptide-level row
- supporting PSM counts
- modified/unmodified support if present

## Example 3: Summarize score-like upstream fields

**User prompt**

Load this MaxQuant `msms.txt` file and summarize the confidence, charge-state, and precursor m/z distributions.

**What the MCP should do**

- call `summarize_numeric_field` for `confidence`
- call `summarize_categorical_field` for `precursor_charge`
- call `bin_numeric_field` for `precursor_mz`

**Expected output shape**

- numeric summary stats for confidence
- counts/fractions for charge states
- histogram-style bins for precursor m/z

## Example 4: Inspect table fields before analysis

**User prompt**

List which fields in this Proteome Discoverer export can be summarized numerically or categorically.

**What the MCP should do**

- call `list_distribution_fields`

**Expected output shape**

- canonical fields
- detected upstream `source_columns.*`
- separation of numeric and categorical candidates
