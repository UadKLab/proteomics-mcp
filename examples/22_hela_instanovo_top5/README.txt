HeLa InstaNovo-style top-5 rendered spectra

Top-5 highest-confidence labelled HeLa predictions from `data/results/identifications/instanovo/hela_qc_labelled_kpreds.csv`, rendered against `data/feather/20230606_MK_E1200_MVL_FAIMS_2CV_15cm_70min_0000_DDA_200ng_QC1_Hela_QC_ms2.feather`.

Files in this folder:
- `*_mirror.html`: observed vs predicted mirror plot with score tiles
- `*_annotated.html`: observed-only annotated spectrum using the internal ion annotator with `ion_types=all`
- `*_plain.html`: observed-only unannotated spectrum
- `*_payload.json`: underlying render payload for the same spectrum
- `top5_summary.csv`: compact metric summary for the five examples

All renders use `0.05 Da` tolerance, `max_isotope=2`, and MS2PIP `HCD` prediction for the mirror comparison.
