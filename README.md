# proteomics-mcp

[![tests](https://github.com/UadKLab/proteomics-mcp/actions/workflows/tests.yml/badge.svg)](https://github.com/UadKLab/proteomics-mcp/actions/workflows/tests.yml)
[![Python](https://img.shields.io/badge/python-3.11%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/UadKLab/proteomics-mcp/blob/main/LICENSE)
[![MCP](https://img.shields.io/badge/MCP-compatible-purple.svg)](https://modelcontextprotocol.io/)

<p align="center">
  <img src="https://raw.githubusercontent.com/UadKLab/proteomics-mcp/main/icon.svg" alt="Proteomics MCP" width="540">
</p>

`proteomics-mcp` is a local-first MCP server and Python package for inspecting proteomics identification results and interrogating selected MS2 spectra.

It turns heterogeneous de novo and search-engine outputs into one usable workflow: inspect IDs and peptides, retrieve the relevant spectrum, predict fragments, score the candidate, annotate the observed signal, and render an interactive report. Interactive HTML reports load Plotly from a CDN when opened; see the [privacy policy](https://github.com/UadKLab/proteomics-mcp/blob/main/PRIVACY.md).

## What It Supports

- Input adapters for InstaNovo, Casanovo, MaxQuant, FragPipe, Proteome Discoverer, and user-mapped tables.
- Modified peptidoforms normalized to inline UNIMOD notation, such as `C[UNIMOD:4]` and `M[UNIMOD:35]`.
- Local IPC/Feather spectra, mzML, MGF, Bruker `.d`, and USI retrieval.
- Local MS2PIP and remote Koina/Prosit fragment prediction.
- On-demand spectrum scoring, observed-spectrum annotation, three interactive HTML views, and publication-ready static SVG/PDF/PNG spectrum figures.
- PSM/peptide tables, structured dataset summaries, dashboards, SDRF context, and resumable batch prediction/scoring artifacts.

The package reads the output of de novo sequencing and database search engines (InstaNovo, Casanovo, MaxQuant, FragPipe, Proteome Discoverer, or a mapped table); it does not run them. It also does not perform FDR/q-value analysis or protein inference.

## Quick Start

Install from [PyPI](https://pypi.org/project/proteomics-mcp/):

```bash
pip install proteomics-mcp
proteomics-mcp --help
```

Or run it without installing it into the current environment, with [uv](https://docs.astral.sh/uv/):

```bash
uvx proteomics-mcp --help
```

For an MCP client, launch `proteomics-mcp serve` over stdio; see [Use With An Agent](https://github.com/UadKLab/proteomics-mcp#use-with-an-agent).

### Requirements

- Python 3.11, 3.12, 3.13 or 3.14 on Linux (x86-64), macOS or Windows (x86-64).
- Linux: glibc 2.28 or newer; Python 3.13 and 3.14 need glibc 2.34 or newer (e.g. Ubuntu 22.04+, Debian 12+, RHEL 9+), because `pyopenms` publishes no older wheels for them.
- macOS: Apple silicon needs macOS 14 or newer. Intel Macs need macOS 13 or newer and Python 3.11–3.13; `numba` publishes no Intel macOS wheels for Python 3.14.
- macOS, for local MS2PIP prediction: install the OpenMP runtime with `brew install libomp`. Without it, MS2PIP tools report an error asking you to install it; everything else works.

### Develop from source

```bash
git clone https://github.com/UadKLab/proteomics-mcp.git
cd proteomics-mcp
uv sync --extra dev
uv run pytest -q
```

In a source checkout, prefix the commands below with `uv run`.

Inspect an identification table before committing to prediction or scoring:

```bash
proteomics-mcp inspect-psm-table \
  --path results.csv \
  --adapter instanovo
```

Then render a browseable peptide table:

```bash
proteomics-mcp render-peptide-table \
  --path results.csv \
  --adapter instanovo \
  --output-path peptide_table.html
```

## Predict And Score

Use local MS2PIP by default. It runs on your computer, but downloads its model files once on first use (from `genesis.ugent.be`, with a Zenodo fallback) and caches them in `~/.ms2pip`; after that, predictions work offline. Repeated single-spectrum calls use a bounded in-process cache; reusable batch workflows use MS2PIP's native batch API. For a remote Prosit prediction through Koina:

```bash
proteomics-mcp predict-fragments \
  --backend koina \
  --model Prosit_2020_intensity_HCD \
  --peptidoform 'PEPM[UNIMOD:35]TIDE' \
  --charge 2
```

Inspect a new Koina model before a large batch:

```bash
proteomics-mcp describe-predictor-model \
  --model Prosit_2020_intensity_HCD
```

## Use With An Agent

An MCP-capable client launches `proteomics-mcp serve` over stdio and exposes the tools to an LLM. The server runs locally and reads your spectra and result tables from your computer. Optional tools contact external services; see the [privacy policy](https://github.com/UadKLab/proteomics-mcp/blob/main/PRIVACY.md).

For Claude Desktop, install the `.mcpb` file attached to the latest [GitHub release](https://github.com/UadKLab/proteomics-mcp/releases). The extension installs its Python dependencies through the managed `uv` runtime. Select an output directory in the extension settings; local input paths must be accessible to Claude Desktop.

For other MCP clients, run the package from PyPI with `uvx`. Pin a version, such as `proteomics-mcp@0.1.10`, for reproducible results. Step-by-step setup for Claude Desktop, Claude Code, Codex CLI, Cursor and Windsurf is in [Connect an MCP client](docs/clients.md):

```json
{
  "mcpServers": {
    "proteomics-mcp": {
      "command": "uvx",
      "args": ["proteomics-mcp@0.1.10", "serve"]
    }
  }
}
```

To run a source checkout instead, for development:

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

Do not run `proteomics-mcp serve` in a terminal by itself. It expects JSON-RPC messages from an MCP client.

## Example Requests

- "Inspect `/path/to/results.csv` as an InstaNovo PSM table. Show the number of identifications, confidence distribution, and top five rows." The `inspect_psm_table` tool returns a structured dataset summary and a constrained preview.
- "Score the first ten InstaNovo candidates in `/path/to/results.csv` against `/path/to/spectra.ipc` with local MS2PIP, save the scores to `/path/to/scored.csv`, and summarize matched ions." The `score_candidates` tool scores a constrained set and writes the scored table; `summarize_score_metrics` explains individual results.
- "Convert `/path/to/run.mgf` to `/path/to/run.feather`, then export that spectrum table as `/path/to/instanovo.ipc`." The `mgf_to_feather` and `export_to_instanovo` tools write local files and report spectrum counts and output paths.

For more task-specific examples, see the [prompt examples](https://github.com/UadKLab/proteomics-mcp/blob/main/docs/prompt_examples/00_overview.md).

## Privacy And Support

See the [privacy policy](https://github.com/UadKLab/proteomics-mcp/blob/main/PRIVACY.md) for local file handling, optional network requests, and generated report behavior. Report issues through [GitHub Issues](https://github.com/UadKLab/proteomics-mcp/issues).

For software citation, use the metadata in [CITATION.cff](https://github.com/UadKLab/proteomics-mcp/blob/main/CITATION.cff). The manuscript's [supplementary data and figure evidence](https://doi.org/10.5281/zenodo.23295902) are deposited at Zenodo; the original mass-spectrometry files remain in PRIDE under PXD044934, PXD014877, and PXD010595.

## Documentation

- [Documentation site](https://uadklab.github.io/proteomics-mcp/): installation, reference, and prompt examples.
- [Technical reference](https://github.com/UadKLab/proteomics-mcp/blob/main/docs/reference.md): installation, configuration, data formats, workflows, metrics, batch artifacts, and tool groups.
- [Prompt examples](https://github.com/UadKLab/proteomics-mcp/blob/main/docs/prompt_examples/00_overview.md): practical requests for an MCP-capable agent.
- [Examples](https://github.com/UadKLab/proteomics-mcp/blob/main/examples/README.md): rendered spectra and dashboard artifacts included with the repository.
- [GitHub releases](https://github.com/UadKLab/proteomics-mcp/releases): versioned source and build artifacts, and the Claude Desktop extension.
- [Changelog](https://github.com/UadKLab/proteomics-mcp/blob/main/CHANGELOG.md): changes per release.

## Testing

From a source checkout:

```bash
uv run pytest -q
uv run proteomics-mcp --help
```

Tests that download MS2PIP models are marked `network`; deselect them when offline with `uv run pytest -m "not network"`.
