FragPipe top-5 rendered spectra

Top-5 highest-confidence FragPipe PSMs from `data/results/identifications/fragpipe/20180905_QX4_JoMu_SA_EASY12-7_uPAC_500ng_ArabidopsisiRT-F5/psm.tsv`, rendered against `data/feather/fragpipe_20180905_ms2.feather`.

Files in this folder:
- `*_mirror.html`: observed vs predicted mirror plot with score tiles
- `*_annotated.html`: observed-only annotated spectrum using the internal ion annotator with `ion_types=all`
- `*_plain.html`: observed-only unannotated spectrum
- `*_payload.json`: underlying render payload for the same spectrum
- `top5_summary.csv`: compact metric summary for the five examples

All renders use `0.05 Da` tolerance, `max_isotope=2`, and MS2PIP `HCD` prediction for the mirror comparison.
