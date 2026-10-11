# Troubleshooting Prompts

## Example 1: Scoring failed

**User prompt**

Why did this candidate fail to score?

**What the MCP should do**

- inspect the scoring failure row or batch failure table
- identify whether the issue is spectrum lookup, missing prediction, malformed candidate metadata, or another error

**Expected output shape**

- concise root-cause explanation
- relevant failure fields

## Example 2: Batch is slow or silent

**User prompt**

Is this batch still running, and how far along is it?

**What the MCP should do**

- call `inspect_prediction_batch_status` or `inspect_scored_batch_status`

**Expected output shape**

- progress summary with `phase`, `completed`, `total`, and `failed`

## Example 3: Field semantics are unclear

**User prompt**

Explain what `next_beam` and `next_beam_margin` mean in this table.

**What the MCP should do**

- explain the normalized candidate schema
- describe how the values are derived for supported adapters

**Expected output shape**

- concise field-level explanation
- no tool call needed if the schema resource is sufficient

## Example 4: Unsure what to do next

**User prompt**

I have a candidate table and matching spectra. What is the right sequence of MCP steps if I want to inspect first and then validate a few targets?

**What the MCP should do**

- recommend the normal workflow:
  - inspect table
  - choose targets
  - score/render selected targets
  - batch only if needed

**Expected output shape**

- short stepwise plan
- optionally point to specific tools by name
