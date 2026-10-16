# Changelog

All notable changes to proteomics-mcp are listed here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/). Detailed notes for each GitHub release are in [`docs/releases/`](https://github.com/UadKLab/proteomics-mcp/tree/main/docs/releases).

## [Unreleased]

### Changed

- **Spectral angle, cosine and Pearson** compare every predicted ion, with zero observed intensity for ions that matched no peak (metric profile `v2`, now the default). Before, only matched ions were compared, so a wrong candidate that matched a single peak scored 1.0, and on real data wrong candidates scored higher than correct ones. Scores change for every PSM with unmatched predicted ions. Set `metric_profile = "legacy_v1"` in a configuration file to reproduce scores from 0.1.10 and earlier, including the manuscript figures ([#72](https://github.com/UadKLab/proteomics-mcp/issues/72)).
- **Definition change: predicted-ion cut-off.** Under metric profile `v2`, a predicted ion is scored only if its intensity is positive and at least `min_predicted_intensity` (new setting, default `0.001`) times the most intense predicted ion. Prosit and MS2PIP give the ions they expect to be absent a tiny intensity (Prosit down to ~1e-8 of the base peak) or 0, about a third of all predicted ions. Those ions counted in `total_predicted_ions` and could be matched to noise peaks. This changes what `total_predicted_ions`, `matched_ion_count` and `matched_ion_fraction` mean: they now refer to the ions predicted above the cut-off.
  - Figure 2 (`EFHLNESGDPSSK`): **22/24 → 18/18 matched ions**; explained intensity 0.480 → 0.471, spectral angle 0.833 → 0.834.
  - The 25 reference payloads: median `total_predicted_ions` 24 → 21, `matched_ion_fraction` 0.54 → 0.73. Spectral angle and cosine change by less than 0.01.
  - `min_predicted_intensity = 0` scores every ion with a positive intensity. `legacy_v1` scores every predicted ion, so it reproduces 0.1.10 for MS2PIP and Koina predictions ([#83](https://github.com/UadKLab/proteomics-mcp/issues/83)).
- **One peak, one ion:** under metric profile `v2`, each observed peak explains at most one predicted ion; pairs are assigned closest first. Before, isobaric ions (for example b3/y2) both claimed the same peak, which counted its intensity twice: explained intensity could exceed 1, and matched-ion counts and the hyperscore were inflated. Static figures labelled such a peak with one ion and coloured it as the other. Of the 25 reference payloads, 3 change (Proteome Discoverer ranks 1–3: 1–3 fewer matched ions, explained intensity 0.007–0.048 lower); `legacy_v1` reproduces all of them ([#74](https://github.com/UadKLab/proteomics-mcp/issues/74)).
- **Faster peptide summaries:** the best PSM of each peptide and peptidoform was found by sorting every group separately. The table is now ranked once and each group takes its lowest rank, with the same order and ties as before. On the manuscript's FragPipe table (52,530 PSMs, 27,603 peptides), `summarize_candidate_inspection_frame`, behind `inspect_psm_table` and `render_psm_table`, takes 10 s instead of 131 s, and `summarize_dataset_frame` 37 s instead of 228 s. The summaries are identical ([#68](https://github.com/UadKLab/proteomics-mcp/issues/68)).

### Fixed

- **mzML acquisition metadata:** `mzml_to_feather` left `isolation_target_mz` empty whenever the isolation target equalled the selected-ion m/z, and `isolation_lower_offset`, `isolation_upper_offset` and `ion_injection_time` always empty, because pyOpenMS keeps them elsewhere (Precursor attributes; the `<scan>` element). On the OpenMS BSA1 mzML, the isolation window is now filled for all 1,120 MS2 spectra. Scores are unchanged: the precursor m/z already came from the selected ion ([#77](https://github.com/UadKLab/proteomics-mcp/issues/77)).
- **InstaNovo export:** `export_to_instanovo` wrote missing sequences as the strings `"None"` / `"nan"`, which InstaNovo accepts as peptide labels for training and evaluation. They are now null: InstaNovo 1.3.0 then refuses the table as annotated data instead of using fake labels. `sequences_missing` and `spectra_dropped_without_precursor` report what was not written as a labelled spectrum, and the new `drop_unlabelled` option (CLI `--drop-unlabelled`) leaves out spectra without a sequence, for a labelled training or evaluation set. Missing run names and spectrum ids also became `"None"` / `"nan"`; they now fall back per row to the source file name and `<run>:<scan>` ([#80](https://github.com/UadKLab/proteomics-mcp/issues/80)).
- **Multiplexed SDRF:** for TMT, iTRAQ and SILAC designs (one SDRF row per channel of the same data file), the file's metadata kept only the last channel, so labels, samples, diseases and factor values of the other channels were lost. Every channel is now kept (`channels`, `channel_count`), and fields that differ between channels hold their distinct values joined with `; `, e.g. label `TMT126; TMT127N`. On PXD006430 (TMT10), 72 files × 10 channels: all 720 rows are now represented; before, 72 ([#81](https://github.com/UadKLab/proteomics-mcp/issues/81)).
- **Spectrum lookup never silently returns another spectrum** ([#76](https://github.com/UadKLab/proteomics-mcp/issues/76)):
  - A run name that matches no run of a multi-run spectrum table is an error; before, the first spectrum with that scan number from any run was returned.
  - A scan number or spectrum id that matches several different spectra is an error, instead of returning the first; identical copies from MGF charge-state enumeration are still accepted.
  - Spectra without a precursor charge (or m/z, or RT) load; a missing value read from Feather raised `cannot convert float NaN to integer`.
  - Scan numbers are no longer taken from any number in an identifier. Before, SCIEX native IDs gave the sample number (`sample=1 … cycle=812` → 1), `QE-HF-01.5000.5000.2` gave 1, `Run 7, scans: 1234` gave 7, and Bruker `merged=102 frame=5 …` gave 102.
- **timsTOF (`bruker_d_to_feather`)** ([#95](https://github.com/UadKLab/proteomics-mcp/issues/95)):
  - `scan_number` is now the precursor `Id`. Before, it was the precursor's truncated TIMS scan position, which repeats across frames: on a 21-min HeLa ddaPASEF run, 55,897 MS2 spectra shared 709 scan numbers, up to 164 per number, and lookups by scan number returned the first of them.
  - `mobility` is now 1/K0 (V·s/cm²) instead of the scan position; the scan position is in the new `scan_index` column.
  - Precursors without a charge (951 in that run) can be loaded.
- **Modification parsing:** several notations lost modifications, mostly N-terminal ones, and the peptides were scored as unmodified ([#75](https://github.com/UadKLab/proteomics-mcp/issues/75)):
  - MaxQuant: N-terminal modifications (`_(ac)PEPM(ox)K_`), and full names with nested parentheses (`(Oxidation (M))`).
  - FragPipe: `N-term(...)` / `C-term(...)` in `Assigned Modifications`. `Modified Peptide` masses are totals (`n[43]` is acetyl, not carbamyl; `M[147]` is oxidised Met), and both readings are tried. On the arabidopsis table of the manuscript supplement, 503 of 52,530 PSMs change: 406 carried carbamyl instead of acetyl, 97 had lost the acetyl. All 52,530 now match FragPipe's `Calculated Peptide Mass`.
  - InstaNovo and Casanovo: masses like `+.98`, parenthesised N- and C-terminal tokens (`(+42.01)PEPTIDE`) and unbracketed masses (`+42.011PEPM+15.995K`).
  - Proteome Discoverer: the `1xTMT6plex [N-Term]; 2xCarbamidomethyl [C5; C12]` notation, sites checked against the residue, and Met-loss no longer subtracted twice.
  - `summarize_prediction_accuracy` strips modifications in every notation.
- **Unplaceable modifications are flagged, not dropped:** a modification that cannot be resolved or placed stays in the peptidoform as an explicit token, and `parse_peptidoform` flags any text it cannot place as an unsupported modification ([#75](https://github.com/UadKLab/proteomics-mcp/issues/75)).
- **w-ion masses:** C-terminal fragments (x, y, z, w) took the residue at the wrong end of the peptide as their cleavage-site residue. So w and w-alt ions lost the side chain of the wrong residue, and the `residue` reported for x/y/z/w fragments was wrong in every peptide longer than two residues. Observed-spectrum annotation with `ion_types=all` or `exotic` changes: 20 of the 25 reference payloads get different w-ion assignments, annotated peak counts and annotated intensity fractions. Scores are not affected ([#79](https://github.com/UadKLab/proteomics-mcp/issues/79)).
- **Precursor mass error with modifications:** the theoretical precursor m/z came from pyOpenMS, which recognised only a few UNIMOD ids and silently left out other modifications, and read named ones such as `C[Carbamidomethyl]` as extra residues. TMT, N-terminal and many other modified PSMs were off by up to hundreds of Da. The mass now comes from the package's own peptidoform parser and includes every residue and terminal modification. The error is empty when a modification cannot be resolved. For unmodified peptides, values change by about 0.006 ppm because the residue masses are slightly more exact ([#73](https://github.com/UadKLab/proteomics-mcp/issues/73)).
- **InstaNovo prediction CSVs:** every output format of InstaNovo 1.0.0, 1.1.0, 1.1.3 and 1.2.2 in Zenodo [10.5281/zenodo.20756892](https://doi.org/10.5281/zenodo.20756892) (32 files) is now read correctly:
  - **Positions, not scan numbers:** for MGF input, and mzML input in 1.0.0, InstaNovo's `scan_number` is the spectrum's 0-based position in the file. It was used as a scan number, so spectrum lookup returned the wrong spectrum for 88% of rows and none for 12%. It is now stored as `spectrum_index`, and `load_spectrum_from_table` / scoring look spectra up by position within the run.
  - **1.1.0 transformer files** (`predictions` with `log_probabilities`) failed with "Could not detect an InstaNovo prediction schema".
  - **InstaNovo+ refined files:** `final_prediction` is used. Where refinement produced nothing, the transformer prediction is used instead of dropping the row (39% of rows in 1.1.0-alpha refined files).
  - **Confidence:** a log-probability just above 0 (e.g. 1.9e-08) gave a confidence of 1.9e-08 instead of 1.
  - **Beam files:** `next_beam` is read in every layout, when the beams rank alternatives to the reported prediction. Token log-probabilities of 1.0.x/1.1.x beam files are aligned to the sequence: they carry a trailing extra value, and 1.1.x stores them in reverse order.
- **MaxQuant fixed modifications:** `Modified sequence` lists only variable modifications, so every cysteine of a carbamidomethylated sample was predicted and scored unmodified. Fixed modifications are now read from the run's `parameters.txt` or `mqpar.xml`, or given with `fixed_modifications` / `--fixed-modifications`. In the manuscript supplement, 5,966 of 22,782 ProteomeTools PSMs contain cysteine, including the Figure 3B PSM `AQNSCPICQK` (now `AQNSC[UNIMOD:4]PIC[UNIMOD:4]QK`). Without a source, candidates carry a warning ([#101](https://github.com/UadKLab/proteomics-mcp/issues/101)).
- **HTML escaping:** run names, spectrum ids, peptides, file names, titles and summary keys and values were written into HTML reports as-is. A value containing `</script>` or `<img onerror=...>` broke the page or ran script when the report was opened. Text is now HTML-escaped. JSON in `<script>` blocks (Plotly figures, embedded CSV, spectrum payloads) escapes `<`, `>` and `&`, which keeps the data unchanged. This covers PSM and peptide tables, table browse, summary, dashboard, prediction accuracy, dataset report and spectrum pages ([#78](https://github.com/UadKLab/proteomics-mcp/issues/78)).

### Added

- Built-in masses for amidation, dimethyl, GG, iTRAQ4plex, ammonia loss, TMT6plex/10plex/11plex and TMTpro, so they resolve without downloading the Unimod catalogue.
- **Isotope-corrected precursor error:** `precursor_isotope_error` (0, 1 or 2) and `precursor_mass_error_corrected_ppm` / `precursor_mass_error_corrected_da` report the error relative to the 13C isotope peak closest to the observed m/z. Before, a precursor selected on its second isotope peak showed an error of hundreds of ppm ([#73](https://github.com/UadKLab/proteomics-mcp/issues/73)).
- **Per-candidate `warnings`:** problems found while reading a row that may make its peptidoform or scores wrong. MaxQuant candidates are checked against MaxQuant's `Mass` column; a difference above 0.01 Da is reported there ([#101](https://github.com/UadKLab/proteomics-mcp/issues/101)).

## [0.1.10] - 2026-10-14

Fixes from the end-to-end check of the published package on the manuscript's public data ([#30](https://github.com/UadKLab/proteomics-mcp/issues/30)). This is the version cited in the manuscript. How scores, annotations and metrics are computed is unchanged.

### Fixed

- **Static figures over MCP:** `render_static_spectrum` failed on macOS, and on Windows with a broken Tcl/Tk it aborted the server. Figures are now drawn without pyplot or a GUI backend and no longer change matplotlib's global settings; the output is unchanged ([#35](https://github.com/UadKLab/proteomics-mcp/issues/35), [#56](https://github.com/UadKLab/proteomics-mcp/pull/56)).
- **Error messages:** with mcp 2.3, MCP clients saw only "Error executing tool". Errors the agent can act on (wrong path, unsupported adapter, malformed table or argument, unknown field, unreachable service) now show their type and message ([#27](https://github.com/UadKLab/proteomics-mcp/issues/27), [#57](https://github.com/UadKLab/proteomics-mcp/pull/57)).
- **Generic adapter:** `adapter="generic"` could not be used over MCP or the CLI. Every tool and command that takes an adapter accepts a JSON column mapping (`column_mapping_path` / `--column-mapping`) ([#38](https://github.com/UadKLab/proteomics-mcp/issues/38), [#58](https://github.com/UadKLab/proteomics-mcp/pull/58)).
- **Scored tables:** a scored table joined nothing, silently, when the candidate table had been renamed. Unmatched rows now join on run + scan + peptidoform, one to one only. Tools report the join (`scored_join`) and warn when nothing matched or rows were ambiguous ([#39](https://github.com/UadKLab/proteomics-mcp/issues/39), [#59](https://github.com/UadKLab/proteomics-mcp/pull/59)).
- **PRIDE SDRF:** `retrieve_pride_sdrf` searched only the first 100 project files and could not download what it found. It now pages through the listing and downloads over HTTPS ([#40](https://github.com/UadKLab/proteomics-mcp/issues/40), [#60](https://github.com/UadKLab/proteomics-mcp/pull/60)).
- **USI spectra:** spectra retrieved by USI had no precursor m/z, so `precursor_signal_fraction` was missing ([#41](https://github.com/UadKLab/proteomics-mcp/issues/41), [#61](https://github.com/UadKLab/proteomics-mcp/pull/61)).
- **Batch rendering:** `render_spectra_batch` reads JSON Lines and gzip-compressed specs ([#45](https://github.com/UadKLab/proteomics-mcp/issues/45), [#63](https://github.com/UadKLab/proteomics-mcp/pull/63)).
- **Batch status:** prediction and scoring batches write their manifest at the start, so their status can be read while they run or after an interruption (`seconds_since_update`). Scoring refuses an unfinished prediction batch ([#46](https://github.com/UadKLab/proteomics-mcp/issues/46), [#64](https://github.com/UadKLab/proteomics-mcp/pull/64)).

### Added

- **Predictor version:** every prediction and score records it: the MS2PIP release, or the Koina model and version ([#22](https://github.com/UadKLab/proteomics-mcp/issues/22), [#62](https://github.com/UadKLab/proteomics-mcp/pull/62)).
- **Docs:** `docs/clients.md`, setup for Claude Desktop, Claude Code, Codex CLI, Cursor and Windsurf ([#53](https://github.com/UadKLab/proteomics-mcp/issues/53), [#65](https://github.com/UadKLab/proteomics-mcp/pull/65)).
- **End-to-end tests** (`tests/e2e`, `e2e.yml`): they check the published package against the Zenodo supplement, PRIDE spectra by USI, an OpenMS mzML, a PRIDE MGF and a timsTOF run ([#30](https://github.com/UadKLab/proteomics-mcp/issues/30), [#37](https://github.com/UadKLab/proteomics-mcp/pull/37), [#67](https://github.com/UadKLab/proteomics-mcp/pull/67)).

### Changed

- **Docs:** spectra retrieved by USI can differ from local conversions of the same scan ([#31](https://github.com/UadKLab/proteomics-mcp/issues/31), [#44](https://github.com/UadKLab/proteomics-mcp/pull/44)).

## [0.1.9] - 2026-10-13

Packaging, CI and documentation only. No changes to adapters, scoring, annotation, metrics or rendering.

### Fixed

- The Claude Desktop extension (`.mcpb`) failed to start when uv selected Python 3.14, because the locked `ms2rescore-rs` 0.4.3 has no Python 3.14 wheel ([#21](https://github.com/UadKLab/proteomics-mcp/issues/21)). The bundle now ships `.python-version` (3.12), and `uv.lock` pins `ms2pip` 4.2.0 / `ms2rescore-rs` 0.6.0 ([#20](https://github.com/UadKLab/proteomics-mcp/pull/20)). Workaround for v0.1.8: install Python 3.12, or set `UV_PYTHON=3.12`.
- MCP `initialize` reports the package version instead of an empty string ([#19](https://github.com/UadKLab/proteomics-mcp/pull/19)).
- On macOS without Homebrew `libomp`, MS2PIP prediction now fails with an actionable error (`brew install libomp`) instead of a raw `XGBoostError`. The message also reaches MCP clients ([#25](https://github.com/UadKLab/proteomics-mcp/pull/25)).
- The sdist ships the complete `tests/` directory, so its tests run against the installed wheel ([#18](https://github.com/UadKLab/proteomics-mcp/pull/18)).
- Removed the unused `JOBS_DIR` constant, which resolved outside the package when installed ([#18](https://github.com/UadKLab/proteomics-mcp/pull/18)).

### Changed

- **MS2PIP 4.2.0** in the locked environment used by the Claude Desktop extension (previously 4.1.2). Predictions are bit-identical for 2+ precursors. For 3+ precursors (`HCDch2`) they shift slightly, changing spectral angle, cosine and Pearson by less than 1e-7 on the example PSMs ([#20](https://github.com/UadKLab/proteomics-mcp/pull/20), [#23](https://github.com/UadKLab/proteomics-mcp/issues/23)).
- **Project metadata** ([#25](https://github.com/UadKLab/proteomics-mcp/pull/25)):
  - authors from `CITATION.cff`
  - classifiers for Python 3.11–3.14, `Development Status :: 4 - Beta` and `Intended Audience :: Science/Research`
  - `license-files`, a `setuptools>=77` build backend and a Changelog URL
- **Dependency:** dropped the unused `typer` ([#25](https://github.com/UadKLab/proteomics-mcp/pull/25)).
- **Docs** ([#25](https://github.com/UadKLab/proteomics-mcp/pull/25)):
  - README and install docs lead with `pip install proteomics-mcp` / `uvx proteomics-mcp serve`, and list platform requirements: glibc ≥ 2.34 for Python 3.13+ on Linux, macOS 14+ on Apple silicon, `brew install libomp` for MS2PIP on macOS.
  - The README, docs and `PRIVACY.md` state that MS2PIP downloads its models once, on first use.
- **Version:** `__version__` is read from the installed package metadata; `pyproject.toml` is the single source of truth ([#19](https://github.com/UadKLab/proteomics-mcp/pull/19)).

### Added

- **Release workflow:** tag-triggered with PyPI Trusted Publishing. It runs build → verify → TestPyPI → PyPI → GitHub release, plus a TestPyPI-only manual dry run ([#24](https://github.com/UadKLab/proteomics-mcp/pull/24)).
- **CI** tests the built wheel on Linux, macOS and Windows with Python 3.11–3.14, plus the Claude Desktop extension launch command ([#20](https://github.com/UadKLab/proteomics-mcp/pull/20)).
- **Scripts:**
  - `scripts/check_version.py`: version consistency across `pyproject.toml`, `manifest.json`, `CITATION.cff`, the tag and the release notes ([#19](https://github.com/UadKLab/proteomics-mcp/pull/19)).
  - `scripts/check_wheel.py`: wheel contents ([#18](https://github.com/UadKLab/proteomics-mcp/pull/18)).
  - `scripts/smoke_mcp.py`: MCP stdio smoke test ([#20](https://github.com/UadKLab/proteomics-mcp/pull/20)).
- **Docs:** `RELEASING.md` (release and contribution process) and this changelog ([#28](https://github.com/UadKLab/proteomics-mcp/pull/28)).

### Known issues

- With `mcp` ≥ 2.3, which fresh `pip`/`uvx` installs get, tool error messages are hidden from MCP clients ([#27](https://github.com/UadKLab/proteomics-mcp/issues/27)).
- On macOS, MS2PIP needs `brew install libomp`.

## [0.1.8] - 2026-10-11

This is the version cited in the manuscript, and the first release from the public repository snapshot. Scoring behaviour is unchanged from 0.1.7. The [release notes](https://github.com/UadKLab/proteomics-mcp/blob/main/docs/releases/v0.1.8.md) list the changes.

The GitHub release was published on 2026-10-11. Its wheel and sdist were uploaded unchanged to TestPyPI and PyPI on 2026-10-12; see [RELEASING.md](https://github.com/UadKLab/proteomics-mcp/blob/main/RELEASING.md#history-manual-upload-of-018) for the SHA-256 values.

## Earlier releases

Releases 0.1.1–0.1.7 predate the public repository and were not published to PyPI. Their notes:

- [0.1.7](https://github.com/UadKLab/proteomics-mcp/blob/main/docs/releases/v0.1.7.md)
- [0.1.6](https://github.com/UadKLab/proteomics-mcp/blob/main/docs/releases/v0.1.6.md)
- [0.1.5](https://github.com/UadKLab/proteomics-mcp/blob/main/docs/releases/v0.1.5.md)
- [0.1.4](https://github.com/UadKLab/proteomics-mcp/blob/main/docs/releases/v0.1.4.md)
- [0.1.3](https://github.com/UadKLab/proteomics-mcp/blob/main/docs/releases/v0.1.3.md)
- [0.1.2](https://github.com/UadKLab/proteomics-mcp/blob/main/docs/releases/v0.1.2.md)
- [0.1.1](https://github.com/UadKLab/proteomics-mcp/blob/main/docs/releases/v0.1.1.md)

[Unreleased]: https://github.com/UadKLab/proteomics-mcp/compare/v0.1.10...HEAD
[0.1.10]: https://github.com/UadKLab/proteomics-mcp/compare/v0.1.9...v0.1.10
[0.1.9]: https://github.com/UadKLab/proteomics-mcp/compare/v0.1.8...v0.1.9
[0.1.8]: https://github.com/UadKLab/proteomics-mcp/releases/tag/v0.1.8
