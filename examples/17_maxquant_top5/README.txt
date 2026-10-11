MaxQuant top-5 rendered spectra

Top-5 highest-confidence MaxQuant PSMs from `data/results/identifications/maxquant/txt/msms.txt`, rendered against `data/feather/maxquant_01974_r4_ms2.feather`.

Files in this folder:
- `*_mirror.html`: observed vs predicted mirror plot with score tiles
- `*_annotated.html`: observed-only annotated spectrum using the internal ion annotator with `ion_types=all`
- `*_plain.html`: observed-only unannotated spectrum
- `*_payload.json`: underlying render payload for the same spectrum
- `top5_summary.csv`: compact metric summary for the five examples

All renders use `0.05 Da` tolerance, `max_isotope=2`, and MS2PIP `HCD` prediction for the mirror comparison.
