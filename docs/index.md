# Proteomics MCP

Proteomics MCP is a local-first MCP server and Python package for inspecting proteomics identification results and selected MS2 spectra. It normalizes outputs from InstaNovo, Casanovo, MaxQuant, FragPipe, Proteome Discoverer, and user-mapped tables, then supports targeted spectrum retrieval, fragment prediction, scoring, annotation, and visualization.

It is an inspection and evidence layer, not a de novo sequencing model, search engine, FDR workflow, or protein-inference pipeline.

## Install

Python 3.11 to 3.14 is supported. Install from [PyPI](https://pypi.org/project/proteomics-mcp/):

```bash
pip install proteomics-mcp
proteomics-mcp --help
```

To run the server from an MCP-capable client without installing it in the current environment:

```bash
uvx proteomics-mcp serve
```

Platform requirements (glibc for Python 3.13+ on Linux, macOS versions, and `brew install libomp` for local MS2PIP on macOS) are listed in the [README](https://github.com/UadKLab/proteomics-mcp#requirements).

The server uses stdio and must be launched by an MCP client, not run interactively in a terminal. For Claude Desktop, the [GitHub release](https://github.com/UadKLab/proteomics-mcp/releases) also contains a `.mcpb` extension.

## Typical workflow

1. Inspect an upstream PSM table and its available fields.
2. Select a peptide, PSM, or constrained subset.
3. Retrieve the corresponding MS2 spectrum, predict fragments if needed, and score or annotate the candidate.
4. Render a spectrum or aggregate saved scores into a table or dashboard.

Start with the [technical reference](reference.md) for formats and commands, the [prompt examples](prompt_examples/00_overview.md) for agent workflows, or [an agent session](tutorials/agent_session.md) that follows one analysis from download to figure.

## Data and privacy

Local input and output paths remain on your computer unless you explicitly use a network-backed tool. Koina predictions, USI retrieval, and PRIDE metadata lookup contact external services. Interactive HTML reports load Plotly from a CDN when opened. Review the [privacy policy](https://github.com/UadKLab/proteomics-mcp/blob/main/PRIVACY.md) before using sensitive data.

Example outputs are available in the [source repository](https://github.com/UadKLab/proteomics-mcp/tree/main/examples). The manuscript's [supplementary data and figure evidence](https://doi.org/10.5281/zenodo.23295902) include processed identification and prediction tables, search settings, rendered figures, and an inventory of exact raw-file names. Original vendor RAW files are available from PRIDE under PXD044934, PXD014877, and PXD010595. Raw datasets and large converted spectra are not distributed with the Python package.
