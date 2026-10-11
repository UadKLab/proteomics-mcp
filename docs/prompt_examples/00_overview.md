# Prompt Examples Overview

These examples are meant to teach the intended operating model of `proteomics-mcp`.

The main pattern is:

1. inspect or summarize an upstream ID/PSM table
2. select one peptide, PSM, or constrained subset
3. score, annotate, or render only the selected target
4. aggregate scored outputs only when needed

Each example includes:

- **User prompt**: what a user could ask an LLM or MCP-capable assistant
- **What the MCP should do**: the intended tool sequence
- **Expected output shape**: the kind of structured response or artifact the user should expect

## Files

- [01_table_inspection.md](01_table_inspection.md)
- [02_spectrum_level.md](02_spectrum_level.md)
- [03_batch_prediction_and_scoring.md](03_batch_prediction_and_scoring.md)
- [04_dataset_aggregation.md](04_dataset_aggregation.md)
- [05_sdrf_and_metadata.md](05_sdrf_and_metadata.md)
- [06_conversion_and_export.md](06_conversion_and_export.md)
- [07_rendering_and_dashboards.md](07_rendering_and_dashboards.md)
- [08_troubleshooting.md](08_troubleshooting.md)

## Guidance

Good prompts for this MCP usually:

- name the table or spectrum source
- indicate the adapter when obvious (`instanovo`, `fragpipe`, `maxquant`, `proteome_discoverer`)
- state whether the goal is inspection, scoring, rendering, or aggregation
- constrain the target set (`top 10`, one peptide, one run, one scan, one PRIDE accession)

Less useful prompts are vague requests like:

- "analyze this proteomics file"
- "score everything"
- "make all plots"

Those are underspecified and usually not the right operating pattern for this server.
