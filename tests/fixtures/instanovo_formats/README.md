# InstaNovo output formats

Row excerpts of InstaNovo prediction CSVs from Zenodo record
[10.5281/zenodo.20756892](https://doi.org/10.5281/zenodo.20756892) ("InstaNovo prediction sample files for PDV import
validation", CC0 1.0). They were predicted from `SF_200217_U2OS_TiO2_HCD_OT_rep1` (PDV test data; 39,151 MS2 spectra)
with InstaNovo 1.0.0, 1.1.0, 1.1.3 and 1.2.2: transformer greedy and beam search, InstaNovo+ with and without
refinement, verbose prediction columns, and combined runs.

File names: `v<InstaNovo version>_<model>_<decoding>_<input format>.csv`. Each excerpt keeps the header and the first three
rows of the original file (row 0 shows whether `scan_number` holds positions), plus rows that cover a specific case:

- an empty InstaNovo+ prediction, where the transformer prediction applies (`v1.1.0_instanovoplus_alpha_refined_mgf.csv`)
- a `final_prediction` that differs from `diffusion_predictions` (`v1.1.3_instanovoplus_refined*`)
- a log-probability just above 0 (InstaNovo+ files)
- an empty prediction with log-probability `-inf` (`v1.0.0_transformer_beam5_mgf.csv`, row 0)
- the same spectrum and peptide in a greedy and a beam-search file, to check token log-probability order
  (scans 64, 12 and 12 in the 1.0.0, 1.1.0 and 1.1.3 transformer files)
