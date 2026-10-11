# Conversion And Export Prompts

## Example 1: Convert mzML

**User prompt**

Convert this mzML file to Feather and tell me how many MS2 spectra were extracted.

**What the MCP should do**

- call `mzml_to_feather`

**Expected output shape**

- output file path
- source format
- row count or extracted spectrum count

## Example 2: Convert Bruker `.d`

**User prompt**

Convert this Bruker `.d` directory to Feather and tell me whether AlphaTims or the embedded MGF fallback was used.

**What the MCP should do**

- call `bruker_d_to_feather`

**Expected output shape**

- output file path
- source format
- fallback indicator when applicable

## Example 3: Convert and export to InstaNovo

**User prompt**

Convert this MGF-derived Feather file into an InstaNovo-ready IPC file.

**What the MCP should do**

- call `export_to_instanovo`

**Expected output shape**

- output path
- output format
- any export metadata returned by the tool

## Example 4: Ask whether a table is suitable downstream

**User prompt**

Check whether this converted Feather table contains the fields needed for downstream scoring and export.

**What the MCP should do**

- inspect the spectrum table externally if needed
- explain required fields based on the MCP schemas and converter contracts

**Expected output shape**

- concise checklist of required fields
- identification of missing fields, if any
