Sbrodae Casanovo top-5 rendered spectra

Top-5 highest-confidence Casanovo predictions from `data/results/identifications/casanovo/20230608_EV_EVO_KK_FAIMS_2CV_Endurance_30SPD_1475_Sbrodae.csv`, with unique scan linkage recovered through the matched Proteome Discoverer export from the same run.

Files in this folder:
- `*_mirror.html`: observed vs predicted mirror plot with score tiles
- `*_annotated.html`: observed-only annotated spectrum using the internal ion annotator with `ion_types=all`
- `*_plain.html`: observed-only unannotated spectrum
- `*_payload.json`: underlying render payload for the same spectrum
- `top5_summary.csv`: compact metric summary for the five examples

All renders use `0.05 Da` tolerance, `max_isotope=2`, and MS2PIP `HCD` prediction for the mirror comparison.
