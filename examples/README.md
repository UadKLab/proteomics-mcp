# Proteomics MCP Examples

This folder contains rendered example outputs from selected non-smORF datasets. They are the reference examples of the manuscript's Zenodo supplement ([record 23295902](https://doi.org/10.5281/zenodo.23295902)): the identification tables are deposited there, and the original RAW files are in PRIDE (PXD044934, PXD014877, PXD010595).

**Version:** the payloads, HTML views and `top5_summary.csv` files are scored with proteomics-mcp **0.2.0** (metric profile `v2`, Comet XCorr and X! Tandem hyperscore) and MS2PIP 4.1.2. `regeneration.json` records the package version, its source digest (equal to that of the `v0.2.0` tag), the commit, the command and what changed. Until 0.2.0 they held the 0.1.x scores. The file names are kept, so three MaxQuant names still show the peptide without its carbamidomethylated cysteines (`AQNSCPICQK`, `KAVEATQVCR`); the payloads carry `C[UNIMOD:4]` (#101).

The spectrum-rendering example sets are:

- `17_maxquant_top5`: top-5 MaxQuant PSMs rendered from `msms.txt` against the matched mzML-derived Feather table
- `18_fragpipe_top5`: top-5 FragPipe PSMs rendered from the Arabidopsis `psm.tsv` against the matched mzML-derived Feather table
- `19_proteome_discoverer_top5`: top-5 Proteome Discoverer PSMs rendered from the Sbrodae PD export
- `22_hela_instanovo_top5`: top-5 HeLa InstaNovo-style labelled predictions rendered from the HeLa mzML-derived Feather table
- `23_sbrodae_casanovo_top5`: top-5 Sbrodae Casanovo predictions with unique scan linkage recovered through the matched PD export

Each spectrum example folder contains:

- `*_mirror.html`: mirror plot with observed spectrum, MS2PIP-predicted spectrum, and score tiles
- `*_annotated.html`: observed-only annotated spectrum using the internal annotator with `ion_types=all`
- `*_plain.html`: observed-only unannotated spectrum
- `*_payload.json`: underlying render payload
- `top5_summary.csv`: compact metric table for the five selected spectra
- `README.txt`: dataset-specific provenance for that example set

The shared rendering settings in these example folders are:

- `tolerance = 0.05 Da`
- `max_isotope = 2`
- mirror predictions from MS2PIP `HCD`

## Regenerate

From a checkout whose package source is that of the release:

```bash
uv run python scripts/regenerate_examples.py --cache-dir path/to/cache --source-of v0.2.0
```

The script downloads the Zenodo supplement into `--cache-dir` once (md5-checked; `--zenodo-zip` takes a local copy). It resolves each of the 25 PSMs from the deposited identification table with the package's adapter, and scores it against its stored observed spectrum with the settings it was published with, under the package defaults. It then rewrites the payloads, the views and the summaries, and `regeneration.json`. It runs itself in a uv environment with MS2PIP 4.1.2, the release behind the manuscript figures, and uses the code of `scripts/regenerate_figures.py`. The 7 payloads that are also in the Zenodo deposit therefore come out identical to those of its 0.2.0 version. `--source-of v0.2.0` stops the run unless the package source equals the tag's. Running the script again gives the same files.

Other workflow examples kept here:

- `14_scored_batch_demo`: mock batch prediction and scoring artifacts for the batch workflow shape (mock predictor, not re-scored)

Machine-readable index: `manifest.json`.
