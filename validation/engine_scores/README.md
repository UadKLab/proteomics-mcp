# Engine-score validation

proteomics-mcp reports two search-engine scores under metric profile `v2`:

- `xcorr`: Comet's XCorr
- `hyperscore`: X! Tandem's hyperscore

`src/proteomics_mcp/core/engine_scores.py` reimplements both from the engines' source code (#51). This directory checks the reimplementations against the engines themselves:
1. It runs each engine from a BioContainers image pinned by digest.
2. For every PSM the engine reports, at any rank, target or decoy, it scores the same spectrum and peptidoform with proteomics-mcp and compares the two.

| engine | image |
|---|---|
| Comet 2026.01.1 | `quay.io/biocontainers/comet-ms:2026011--h9ee0642_0@sha256:af486865…9c9c` |
| X! Tandem 2015.12.15.2 | `quay.io/biocontainers/xtandem:15.12.15.2--h4464bbb_11@sha256:712774bd…10ea` |

## Run it

From the repository root:

```bash
# against the committed engine output (no Docker needed); also part of the test suite
uv run python validation/engine_scores/validate.py

# run Comet and X! Tandem again and compare (needs Docker; about a minute)
docker pull quay.io/biocontainers/comet-ms:2026011--h9ee0642_0@sha256:af4868654b0652fb55e8b87e8d1da77bed87f92345e4e85895280b4a9e3e9c9c
docker pull quay.io/biocontainers/xtandem:15.12.15.2--h4464bbb_11@sha256:712774bdb3afbfd884d22de3cc2e5b93ccaf0074ba527a3b5b2cba748aac10ea
uv run --with matplotlib python validation/engine_scores/validate.py --run-engines
```

The script writes these files to `validation/engine_scores/out/`:
- `comparison.tsv`: one row per PSM, with the engine score, the proteomics-mcp score and the difference;
- `summary.json`;
- `agreement.png` / `.pdf`.

It exits 1 in either of two cases:
- a PSM differs by more than the tolerance: 0.001 for Comet, which reports 3 decimals, and 0.05 for X! Tandem, which reports 1 decimal;
- with `--run-engines`, the engines' fresh output differs from `golden/psms.tsv`.

The `engine validation` GitHub workflow runs the second command weekly, on release tags, and on pull requests that touch the scores.

## Result

| run | spectra | engine settings | PSMs (decoys) | equal to the engine |
|---|---|---|---|---|
| `comet_high_res` | high-res | `comet -p` defaults (bin 0.02, offset 0, flanking peaks), no enzyme, ±3 Da, top 10 | 1,250 (542) | 1,250 |
| `comet_low_res` | ion trap | bin 1.0005, offset 0.4, M peak only, no enzyme, ±3 Da, top 5 | 1,500 (726) | 1,500 |
| `xtandem_da` | high-res | `default_input.xml`, fragment error 0.4 Da | 117 | 117 |
| `xtandem_ppm` | high-res | fragment error 20 ppm | 106 | 106 |
| `xtandem_reversed` | high-res | reversed sequences only, so every hit is a wrong candidate | 63 (63) | 63 |
| `xtandem_low_res` | ion trap | `default_input.xml` | 92 | 92 |

All 3,128 PSMs are equal to the engines' scores, not merely within tolerance. The searches are deliberately broad (no enzyme, ±3 Da, internal decoys, reversed sequences), so most PSMs are wrong candidates. The agreement therefore holds across the score range, not only for good matches.

## Inputs

All inputs are committed, so the comparison is reproducible without PRIDE or UniProt:

| file | what | where from |
|---|---|---|
| high-res spectra | the 25 reference payloads of `examples/*_top5` (HeLa/InstaNovo, *S. brodae*/Casanovo and PD, *A. thaliana*/FragPipe, ProteomeTools/MaxQuant) | written by `validate.py` from the payloads |
| `data/arabidopsis_usi.mgf` | the 100 *A. thaliana* PSM spectra of the manuscript's Figure 4 subset | PRIDE, by USI (`mzspec:PXD014877:…`); `fetch_usi_spectra.py --check` fetches them again and compares the sha256 in `data/arabidopsis_usi.tsv` |
| `data/bsa1_ms2.mgf` | 300 ion-trap CID spectra (precursor charge 2–3) | OpenMS `BSA1.mzML`, tag Release3.1.0 (sha256-checked); `build_inputs.py` rebuilds it |
| `data/search.fasta` | 300 human Swiss-Prot entries, BSA, and the reference peptides as their own entries | `build_inputs.py` (the human part is a UniProt snapshot of 2026-10-22) |
| `params/` | `comet_high_res.params`, `comet_low_res.params`, `xtandem_default_input.xml` | the engines' own templates; each file's header lists what was changed |

The searches use these modifications:
- fixed carbamidomethyl C, +57.021464;
- variable oxidation M, +15.994915.

## Settings that matter for an exact match

- **X! Tandem "spectrum synthesis".** X! Tandem's default `refine, spectrum synthesis = yes` turns on residue-dependent ion weights in every pass, including the first: ×5 for y ions N-terminal to proline, ×2–5 for b ions after D, N, V, E, Q, I, L. It applies even with refinement off. proteomics-mcp applies it too (`XTandemHyperscoreParams.spectrum_synthesis`).
- **Refinement.** It is off in these searches (`refine = no`), so every reported score comes from the first pass. Refinement can rescore with other settings.
- **Comet's precursor tolerance.** It sets the size of the binned array, so the comparison uses the search's tolerance (±3 Da here). proteomics-mcp's default is Comet's 20 ppm.
