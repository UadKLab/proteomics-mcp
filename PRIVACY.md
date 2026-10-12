# Proteomics MCP privacy policy

Proteomics MCP is an open-source, local MCP server maintained by UadKLab. It has no account system, telemetry, or first-party data collection service. It reads files whose paths you provide to its tools and writes generated files at the paths you request. By default, data stays on your computer, apart from data included in your conversation with the MCP client.

Some tools contact third-party services when you invoke them or select a remote backend:

- Koina predictions send peptide sequences or peptidoforms and prediction parameters, such as charge, collision energy, and model name, to the configured Koina endpoint. Koina may be configured to use a different server.
- USI retrieval sends the requested USI to the selected public spectrum provider. PRIDE SDRF retrieval sends a project accession to EMBL-EBI.
- Local MS2PIP downloads its model files once, on first use, from `genesis.ugent.be` (Ghent University), falling back to Zenodo. The request names only the model file; no peptide or spectrum data is sent. The files are cached in `~/.ms2pip`, and later predictions run offline.
- Generated HTML reports load the Plotly JavaScript library from a public CDN when opened in a browser. The HTML data and plots remain in the local report, but loading the library makes a normal network request to the CDN.

The server may cache Koina prediction responses in the operating system's temporary directory. Generated files and caches remain under your control and are not uploaded by Proteomics MCP. Remove them using your usual file-management tools when they are no longer needed.

External services operate under their own terms and privacy practices. See the [EMBL privacy policy](https://www.embl.org/info/privacy-policy/) and [Plotly privacy policy](https://plotly.com/privacy/). The public Koina service does not currently link to a privacy policy we can verify; use local MS2PIP or a Koina instance you trust for sensitive sequences. For privacy questions or issues, open an issue at [UadKLab/proteomics-mcp](https://github.com/UadKLab/proteomics-mcp/issues).
