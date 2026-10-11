# SDRF And Metadata Prompts

## Example 1: Load a local SDRF

**User prompt**

Load this SDRF file and tell me which raw files map to which assays and samples.

**What the MCP should do**

- call `load_sdrf`

**Expected output shape**

- grouped file metadata
- assay/sample summaries
- factor columns and instrument-related metadata when present

## Example 2: Retrieve SDRF from PRIDE

**User prompt**

Retrieve the SDRF for this PRIDE accession and summarize the sample layout.

**What the MCP should do**

- call `retrieve_pride_sdrf`

**Expected output shape**

- PRIDE-derived SDRF summary
- mapped raw file names, basenames, and stems

## Example 3: Combine ID summaries with SDRF

**User prompt**

Inspect this PSM table and include SDRF-linked sample metadata if it matches the raw files.

**What the MCP should do**

- call `inspect_psm_table` with `sdrf_path` or `pride_project_accession`

**Expected output shape**

- normal inspection summary
- additional `sdrf_summary`

## Example 4: Ask for acquisition context

**User prompt**

For this scored peptide, tell me which runs, samples, fragmentation methods, and instruments support it.

**What the MCP should do**

- if working from a scored table, call `get_peptide_summary`
- optionally enrich through SDRF if available

**Expected output shape**

- peptide-level support summary
- run/sample context where available
