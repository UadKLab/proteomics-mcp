# Proteomics MCP Examples

This folder contains rendered example outputs from selected non-smORF datasets. The source identification tables, raw spectra, and converted spectrum tables referenced below are not distributed in this repository; their public accessions and reuse terms must be documented before these examples are presented as reproducible paper artifacts.

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

Other workflow examples kept here:

- `14_scored_batch_demo`: mock batch prediction and scoring artifacts for the batch workflow shape

Machine-readable index: `manifest.json`.
