# Rendering And Dashboard Prompts

## Example 1: Render a peptide table

**User prompt**

Render the peptide table for this candidate file as HTML so I can browse it.

**What the MCP should do**

- call `render_peptide_table`

**Expected output shape**

- HTML output path
- no scoring unless a scored table is explicitly merged

## Example 2: Build an ID dashboard

**User prompt**

Build an ID dashboard for the top 1000 predictions and include retention time metadata from this Feather file.

**What the MCP should do**

- call `build_id_dashboard`
- use `top_n=1000`
- pass `spectra_path`

**Expected output shape**

- summary JSON path
- enriched CSV path
- dashboard HTML path

## Example 3: Render a batch of spectra

**User prompt**

Render observed, annotated, and mirror plots for these five selected PSMs.

**What the MCP should do**

- prepare a small batch manifest
- call `render_spectra_batch`

**Expected output shape**

- output directory with `comparison/`, `annotated/`, and `observed/`
- render status file

## Example 4: Ask for a rendering mode recommendation

**User prompt**

For this peptide, should I use the mirror plot, annotated observed spectrum, or plain observed spectrum?

**What the MCP should do**

- explain the differences between the three visualization modes
- recommend one based on the user’s goal

**Expected output shape**

- short rationale
- no rendering unless requested
