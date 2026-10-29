"""Assemble the new version of the manuscript's Zenodo supplement from a figure-regeneration run.

    uv run python scripts/build_zenodo_package.py --figures path/to/figures --cache-dir path/to/cache --out path/to/zenodo

The package has the structure of the current deposit (record 23295902): every deposited file is copied unchanged,
then the files the regeneration run staged in <figures>/deposit/ replace their deposited versions or are added. It
adds the run's records (regeneration/: its manifest, figure changes, environment and batch manifests, with the home
directory written as ~), the engine-score validation (validation/: validation/engine_scores/validate.py run against the
committed Comet and X! Tandem output, #154), CHANGES.md (what changed against the previous version and why),
PROVENANCE.md (versions, commands, the RAW -> mzML conversion) and MANIFEST.sha256 (every file). Nothing is uploaded, and
nothing is written inside the repository. Exit codes: 0 done, 2 the regeneration run did not pass its checks, 4 the
deposited supplement's checksum differs, 5 the engine-score validation failed.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
# the record, its checksum and how it is extracted are those the figures were regenerated from
from regenerate_figures import REPO_ROOT, ZENODO_MD5, ZENODO_RECORD, file_digest, git_state  # noqa: E402
from regenerate_figures import extract_supplement as extract_archive  # noqa: E402

ISSUES = "https://github.com/UadKLab/proteomics-mcp/issues"

# Why each kind of file changed (path prefix -> reason with issues), most specific prefix first.
REASONS = [
    ("figure_assets/pdfs/Figure2.pdf", "metric profile v2: spectral angle, cosine and Pearson over all scored predicted ions, one peak per ion, "
     "the 0.1% predicted-ion cut-off (#72, #74, #83); Comet XCorr and X! Tandem hyperscore instead of the 0.1.x variants (#51); MS2PIP 4.1.2 (#23)"),
    ("figure_assets/pdfs/Figure3.pdf", "3B: carbamidomethylated cysteines from the MaxQuant fixed modifications (#101) and corrected w-ion masses (#79); 3A unchanged"),
    ("figure_assets/pdfs/Figure4.pdf", "rebuilt from the first 100 PSMs in the deposited FragPipe PSM table, retained in "
     "source-file order (#42, #47), and their PRIDE spectra with retention times (#135); every PSM joined to its spectrum "
     "(#146); every histogram counts all 100 PSMs, where the published one left out the PSMs at each maximum (#169); "
     "scores as above"),
    ("figure_assets/pdfs/Suppl_figure2.pdf", "panel B: peaks of predicted ions below the 0.1% cut-off are no longer matched (#83)"),
    ("figure_assets/pdfs/Suppl_figure5.pdf", "panels B and C: fewer matched ions under metric profile v2 (#74, #83)"),
    ("figure_assets/interactive/", "rebuilt from the Figure 4 build (#42, #47): scores under metric profile v2 (#72, #74, #83), xcorr and hyperscore "
     "columns (#51), retention times from PRIDE (#135), histograms that count their maximum (#169), precursor errors from the "
     "spectrum with FragPipe's own as engine_precursor_* (#179)"),
    ("processed/arabidopsis/figure4_scored_psms_full.csv", "new: the Figure 4 batch scoring output with every column of the 0.2.0 scoring "
     "schema, next to the compact supplementary_scored_psms.csv (#137, #51)"),
    ("processed/arabidopsis/figure4_scored_ions.jsonl", "new: for each Figure 4 PSM, the predicted and observed intensities behind its "
     "scores and the predicted ions left out, with the reason (#137)"),
    ("processed/arabidopsis/figure4_usi_spectra.mgf", "new: the 100 PRIDE spectra of Figure 4 with their retention times, so the figure can be rebuilt "
     "from this deposit alone (`regenerate_figures.py --figure4-mgf`; #42, #135)"),
    ("processed/arabidopsis/", "rebuilt from the Figure 4 build (#42, #47): scores under metric profile v2 (#72, #74, #83), xcorr and hyperscore "
     "columns (#51), retention times from PRIDE (#135); the dashboard plot data's histograms count their maximum (#169); precursor "
     "errors from the spectrum's precursor m/z, i.e. FragPipe's uncalibrated error, with FragPipe's calibrated one as "
     "engine_precursor_mass_error_ppm (#179)"),
    ("figure_assets/examples/", "re-scored with 0.2.0: metric profile v2 (#72, #74, #83), Comet XCorr and X! Tandem hyperscore (#51), corrected "
     "C-terminal ions (#79), MaxQuant fixed modifications (#101), precursor-related ions from the peptidoform's own mass "
     "(#174), precursor errors from the spectrum for every engine, with the engine's own as engine_precursor_* (#179); "
     "candidate ids follow the deposited table names (#42)"),
    ("regeneration/", "new: the record of the regeneration run: package version and git SHA, MS2PIP 4.1.2, inputs, checks, every changed "
     "value, the environment (uv.lock, requirements, installed packages) and the Figure 4 prediction and scoring manifests"),
    ("validation/", "new: Comet's XCorr and X! Tandem's hyperscore as proteomics-mcp computes them, checked against the engines (#51, "
     "#154): every PSM's comparison, the summary, the agreement plots, the engines' images by digest and the sha256 of every input "
     "and parameter file (versions.json), and the commands to run the check again (README.md)"),
    ("PROVENANCE.md", "new: versions, commands and the RAW -> mzML conversion behind this version"),
]
# Deposited files the regeneration does not produce, with what a reader should know.
NOT_REGENERATED = {
    "figure_assets/pdfs/Figure1.pdf": "overview schematic; no data",
    "figure_assets/pdfs/Suppl_figure1.pdf": "unchanged in 0.2.0: panel A's a/b/y labels (Figure 2's spectrum, neutral losses, isotopes up to 2, "
                                            "0.05 Da) are identical to 0.1.10's; panel B has no annotation",
    "figure_assets/pdfs/Suppl_figure3.pdf": "picture of the PSM table's leftmost columns (identifiers), which do not change; the table itself "
                                            "(interactive HTML, CSV) is regenerated",
    "figure_assets/pdfs/Suppl_figure4.pdf": "picture of the peptide table's leftmost columns; the table itself (interactive HTML, CSV) is regenerated",
}


def linked(text: str) -> str:
    """'#51' -> a link to the issue."""
    import re

    return re.sub(r"#(\d+)", lambda m: f"[#{m[1]}]({ISSUES}/{m[1]})", text)


def reason(path: str) -> str:
    return next((text for prefix, text in REASONS if path.startswith(prefix)), "")


def extract_supplement(archive: Path, target: Path) -> Path:
    if file_digest(archive, "md5") != ZENODO_MD5:
        print(f"{archive}: md5 differs from the deposited supplement ({ZENODO_MD5})", file=sys.stderr)
        sys.exit(4)
    return extract_archive(archive, target)


RAW_CONVERSION = """\
## RAW -> mzML conversion (deposited spectra)

The spectra behind the deposited identifications were converted from Thermo .raw to mzML with:

- CompOmics ThermoRawFileParser 2.0.0-dev, built from the `arm64` branch at commit `e00b388`
  (<https://github.com/CompOmics/ThermoRawFileParser/tree/arm64>), on .NET 8, with Thermo CommonCore 8.0.42.
- Licences: ThermoRawFileParser is Apache-2.0. The Thermo RawFileReader libraries it uses have their own licence
  (<https://github.com/thermofisherlsms/RawFileReader/blob/main/License.doc>).
- Citation: Hulstaert N, Shofstahl J, Sachsenberg T, Walzer M, Barsnes H, Martens L, Perez-Riverol Y.
  ThermoRawFileParser: Modular, Scalable, and Cross-Platform RAW File Conversion. J Proteome Res. 2020;19(1):537-542.
  PMID 31755270.
"""


VALIDATION = REPO_ROOT / "validation" / "engine_scores"
# What a validate.py run writes that goes into the deposit; its high_res.mgf and low_res.mgf are rebuilt from the
# deposited reference payloads and the committed BSA1 spectra, so they are left out
VALIDATION_OUTPUTS = ("comparison.tsv", "summary.json", "agreement.png", "agreement.pdf")


def run_validate(out: Path) -> subprocess.CompletedProcess:
    """validate.py against the committed engine output (golden/psms.tsv): no Docker, the package's own scores."""
    return subprocess.run([sys.executable, str(VALIDATION / "validate.py"), "--out", str(out)], cwd=REPO_ROOT, capture_output=True, text=True)


def validation_readme(summary: dict, versions: dict, packaging: dict, release: str) -> str:
    rows = [f"| `{name}` | {item['engine']} | {item['spectra']} | {item['psms']:,} ({item['decoy_psms']:,}) | {item['exact']:,} | "
            f"{item['tolerance']} | {item['max_difference']} |" for name, item in summary.items()]
    total, exact = sum(item["psms"] for item in summary.values()), sum(item["exact"] for item in summary.values())
    return "\n".join([
        "# Engine-score validation", "",
        "Under metric profile v2, proteomics-mcp reports Comet's XCorr (`xcorr`) and X! Tandem's hyperscore (`hyperscore`), "
        "reimplemented from the engines' source code (#51). This folder checks them against the engines themselves "
        "(`validation/engine_scores/` in the repository, #154): for every PSM an engine reports, at any rank, target or decoy, "
        "proteomics-mcp scores the same spectrum and peptidoform and the two scores are compared.", "",
        f"Generated by `scripts/build_zenodo_package.py` at `{packaging.get('describe')}` (commit `{packaging.get('sha')}`), which runs "
        "`validation/engine_scores/validate.py` against the engines' committed output (`golden/psms.tsv`).", "",
        "| run | engine | spectra | PSMs (decoys) | equal to the engine | tolerance | max difference |",
        "|---|---|---|---|---|---|---|", *rows, "",
        f"{exact:,} of {total:,} PSMs are equal to the engines' scores.", "",
        "## Files", "",
        "- `comparison.tsv`: one row per PSM: the engine's score, the proteomics-mcp score (`ours`), the difference and whether it is within tolerance.",
        "- `summary.json`: the table above.",
        "- `agreement.png`, `agreement.pdf`: proteomics-mcp against the engine, per run.",
        "- `versions.json`: the engines' versions, their BioContainers images pinned by digest, the sha256 of every input and "
        "parameter file, and the PSMs per run.", "",
        "## Run it again", "",
        "From a checkout of the release:", "",
        "```bash",
        f"git clone https://github.com/UadKLab/proteomics-mcp.git && cd proteomics-mcp && git switch --detach v{release}",
        "# against the committed engine output (no Docker)",
        "uv run python validation/engine_scores/validate.py --out out",
        "# run Comet and X! Tandem again and compare (needs Docker)",
        f"docker pull {versions['comet_image']}",
        f"docker pull {versions['xtandem_image']}",
        "uv run --with matplotlib python validation/engine_scores/validate.py --run-engines --out out",
        "```", "",
        "The script exits 1 if a PSM differs by more than its tolerance or, with `--run-engines`, if the engines' fresh output "
        "differs from the committed one.", ""])


def engine_validation(target: Path, packaging: dict, release: str) -> int:
    """Run the engine-score validation and write its outputs, versions.json and a README into target; validate.py's exit code."""
    with tempfile.TemporaryDirectory() as work:
        completed = run_validate(Path(work))
        if completed.returncode != 0:
            print(completed.stdout + completed.stderr, file=sys.stderr)
            return completed.returncode
        missing = [name for name in VALIDATION_OUTPUTS if not (Path(work) / name).is_file()]
        if missing:
            print(f"validate.py wrote no {', '.join(missing)} (the plots need matplotlib)", file=sys.stderr)
            return 1
        target.mkdir(parents=True, exist_ok=True)
        for name in VALIDATION_OUTPUTS:
            shutil.copyfile(Path(work) / name, target / name)
    shutil.copyfile(VALIDATION / "golden" / "versions.json", target / "versions.json")
    summary, versions = json.loads((target / "summary.json").read_text()), json.loads((target / "versions.json").read_text())
    (target / "README.md").write_text(validation_readme(summary, versions, packaging, release))
    return 0


def sanitized(text: str) -> str:
    """Paths of the machine that ran the regeneration, with the home directory written as ~."""
    home = str(Path.home())
    for written in (json.dumps(home)[1:-1], home, Path.home().as_posix()):  # as JSON writes it (\\ on Windows), as is, with /
        text = text.replace(written, "~")
    return text


# MS2PIP's xgboost links @rpath/libomp.dylib; on macOS it resolves to the copy pyOpenMS ships, which the package loads first
OPENMP_MACOS = ("- **OpenMP (macOS):** MS2PIP's xgboost uses the OpenMP runtime bundled with pyOpenMS, so no separate libomp is "
                "installed. Installing a second one (e.g. `brew install libomp`) can abort with OpenMP Error #15 (\"libomp already initialized\").")


def provenance_markdown(figures: dict, release: str, build_command: str, packaging: dict | None = None, validation: bool = False) -> str:
    git, package = figures.get("git", {}), figures["package"]["proteomics-mcp"]
    released = git.get("describe") == f"v{release}" and package == release and not git.get("dirty")
    version_line = (f"proteomics-mcp **{release}**, built from the tag `v{release}` (commit `{git.get('sha')}`)." if released else
                    f"proteomics-mcp **{release} rc**: a release candidate, proteomics-mcp {package} at `{git.get('describe')}` "
                    f"(commit `{git.get('sha')}`{', with uncommitted changes' if git.get('dirty') else ''}). Regenerate from the "
                    f"`v{release}` tag before uploading.")
    fig4 = figures.get("figure4", {})
    return "\n".join([
        "# Provenance", "",
        f"- **Package:** {version_line}",
        f"- **MS2PIP:** {figures.get('ms2pip')} (model HCD), pinned in the figure environment only; the package requires `ms2pip>=4.1`.",
        f"- **Python:** {figures.get('python')} ({sanitized(figures.get('platform', ''))}).",
        "- **Environment:** `regeneration/uv.lock` (the checkout's lock file), `regeneration/requirements_figures.txt` (the lock with "
        "MS2PIP and the figure tools pinned) and `regeneration/environment.txt` (every installed package).",
        *([OPENMP_MACOS] if str(figures.get("platform", "")).lower().startswith(("macos", "darwin")) else []),
        *([f"- **Packaging:** `scripts/build_zenodo_package.py` at `{packaging.get('describe')}` (commit `{packaging.get('sha')}`"
           f"{', with uncommitted changes' if packaging.get('dirty') else ''})."] if packaging else []),
        *(["- **Engine-score validation:** `validation/`, from `validation/engine_scores/validate.py` run by the packaging step "
           "against the engines' committed output; `validation/README.md` has the result and the commands to run it again."] if validation else []),
        f"- **Figure 4 spectra:** PRIDE PXD014877 by USI, in `processed/arabidopsis/figure4_usi_spectra.mgf` "
        f"(sha256 `{(fig4.get('spectra') or {}).get('mgf_sha256')}`).",
        "- **Previous version:** Zenodo record " + ZENODO_RECORD + " reproduces proteomics-mcp 0.1.8 / 0.1.10; this version goes "
        f"with {release}.", "",
        "## Commands", "",
        "```bash", sanitized(figures.get("command", "")), build_command, "```", "",
        "Rebuild Figure 4 from this deposit alone: add `--figure4-mgf processed/arabidopsis/figure4_usi_spectra.mgf` to the first command.", "",
        RAW_CONVERSION])


def changes_markdown(changed: dict[str, list[str]], unchanged: list[str], figures: dict) -> str:
    git, fig4, examples = figures.get("git", {}), figures.get("figure4", {}), figures.get("examples", {})
    lines = [f"# Changes against Zenodo record {ZENODO_RECORD}", "",
             f"This version goes with proteomics-mcp {figures['package']['proteomics-mcp']} ({git.get('describe')}, commit {git.get('sha')}). "
             f"Figures and tables were regenerated by `scripts/regenerate_figures.py` with MS2PIP {figures.get('ms2pip')}, metric profile v2 "
             "and 0.05 Da; every input file of the previous version is kept unchanged. `regeneration/FIGURE_CHANGES.md` lists every changed "
             "value in the figures and the manuscript text to update; `regeneration/manifest.json` records the run. `MANIFEST.sha256` has "
             "the checksum of every file.", ""]
    for title, key in (("Regenerated (content changed)", "replaced"), ("Added", "added"), ("Regenerated, identical to the previous version", "same")):
        if changed[key]:
            lines += [f"## {title}", "", "| File | Why |", "|---|---|"]
            lines += [f"| `{path}` | {linked(reason(path))} |" for path in changed[key]]
            lines.append("")
    lines += ["## Not regenerated", "", "| File | Note |", "|---|---|"]
    lines += [f"| `{path}` | {linked(note)} |" for path, note in NOT_REGENERATED.items()]
    lines += ["", "## Unchanged inputs", "", "Copied byte for byte from the previous version:", ""]
    lines += [f"- `{path}`" for path in unchanged if path not in NOT_REGENERATED]
    lines += ["", "## Figure 4", "", f"Status: **{fig4.get('status')}**.", ""]
    lines += [f"- {linked(note)}" for note in fig4.get("notes", [])]
    lines += [f"- check: {name}: {'passed' if passed else 'FAILED'}" for name, passed in (fig4.get("checks") or {}).items()]
    for key in ("scored_psms_columns", "scored_peptides_columns"):
        columns = fig4.get(key) or {}
        lines.append(f"- `{key.removesuffix('_columns')}` columns removed: {columns.get('removed') or 'none'}; added: {columns.get('added') or 'none'}.")
    lines += ["", "## Reference examples", ""]
    for folder, report in examples.items():
        columns = report.get("top5_summary_columns") or {}
        peptides = "; ".join(f"rank {rank}: `{old}` → `{new}`" for rank, (old, new) in (report.get("peptidoform_changes") or {}).items())
        lines.append(f"- `{folder}`: {report.get('psms')} PSMs re-scored; `top5_summary.csv` columns {columns.get('removed')} → "
                     f"{columns.get('added')}" + (f"; peptidoforms: {peptides} ({linked('#101')})" if peptides else "") + ".")
    lines += ["", "## Pending for this version", "",
              f"- The figures' Illustrator touch-ups ({linked('#149')}): the PDFs here are the script's assembly.", ""]
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--figures", type=Path, required=True, help="output folder of scripts/regenerate_figures.py")
    parser.add_argument("--out", type=Path, required=True, help="new, empty folder outside the repository")
    parser.add_argument("--cache-dir", type=Path, help="the regeneration cache, which holds the deposited supplement.zip")
    parser.add_argument("--zenodo-zip", type=Path, help="the deposited supplement.zip (md5-checked)")
    parser.add_argument("--release", default="0.2.0", help="the release this version of the supplement goes with")
    parser.add_argument("--allow-provisional", action="store_true", help="package a run whose checks or Figure 4 did not pass")
    args = parser.parse_args(argv)

    figures_dir, out = args.figures.expanduser().resolve(), args.out.expanduser().resolve()
    if REPO_ROOT in out.parents or out == REPO_ROOT:
        parser.error("--out must be outside the repository")
    if out.exists() and any(out.iterdir()):
        parser.error(f"{out} is not empty")
    archive = args.zenodo_zip or (args.cache_dir and args.cache_dir / f"zenodo-{ZENODO_RECORD}-supplement.zip")
    if not archive or not archive.expanduser().is_file():
        parser.error("give --zenodo-zip, or the --cache-dir of the regeneration run")
    figures = json.loads((figures_dir / "manifest.json").read_text())
    passed = figures["verification"]["passed"] and figures.get("figure4", {}).get("status") == "REGENERATED"
    if not passed and not args.allow_provisional:
        print(f"the regeneration run did not pass (verification {figures['verification']}, Figure 4 {figures.get('figure4', {}).get('status')})",
              file=sys.stderr)
        return 2

    supplement = extract_supplement(archive.expanduser(), out)
    deposited = {p.relative_to(supplement).as_posix() for p in supplement.rglob("*") if p.is_file()}
    changed: dict[str, list[str]] = {"replaced": [], "added": [], "same": []}
    staged = figures_dir / "deposit"
    for source in sorted(p for p in staged.rglob("*") if p.is_file()):
        rel = source.relative_to(staged).as_posix()
        target = supplement / rel
        kind = "added" if rel not in deposited else "same" if file_digest(target) == file_digest(source) else "replaced"
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, target)
        changed[kind].append(rel)
    records = supplement / "regeneration"
    records.mkdir(exist_ok=True)
    shutil.copyfile(figures_dir / "manifest.json", records / "manifest.json")
    shutil.copyfile(figures_dir / "CHANGES.md", records / "FIGURE_CHANGES.md")
    changed["added"] += ["regeneration/manifest.json", "regeneration/FIGURE_CHANGES.md"]
    rewritten = []
    for path in sorted(p for p in records.rglob("*") if p.suffix in {".json", ".md", ".txt"}):
        text = path.read_text()
        if sanitized(text) != text:
            path.write_text(sanitized(text))
            rewritten.append(path.relative_to(supplement).as_posix())
    packaging = git_state()
    validation = supplement / "validation"
    status = engine_validation(validation, packaging, args.release)
    if status != 0:
        print(f"the engine-score validation failed (validate.py exit code {status})", file=sys.stderr)
        return 5
    changed["added"] += sorted(p.relative_to(supplement).as_posix() for p in validation.rglob("*") if p.is_file())
    # arguments are paths, flags and versions: paths are written with / on every platform
    build_command = "uv run python scripts/build_zenodo_package.py " + " ".join(
        sanitized(a.replace("\\", "/")) for a in (argv if argv is not None else sys.argv[1:]))
    (supplement / "PROVENANCE.md").write_text(provenance_markdown(figures, args.release, build_command, packaging, validation=True))
    changed["added"].append("PROVENANCE.md")
    unchanged = sorted(deposited - set(changed["replaced"]) - set(changed["same"]))
    (supplement / "CHANGES.md").write_text(changes_markdown(changed, unchanged, figures))
    files = sorted(p for p in supplement.rglob("*") if p.is_file() and p.name != "MANIFEST.sha256")
    (supplement / "MANIFEST.sha256").write_text("".join(f"{file_digest(p)}  {p.relative_to(supplement).as_posix()}\n" for p in files))
    if rewritten:
        print(f"home directory written as ~ in: {', '.join(rewritten)}")
    print(f"{out / 'supplement'}: {len(files) + 1} files; {len(changed['replaced'])} regenerated with changes, {len(changed['same'])} regenerated "
          f"identical, {len(changed['added'])} added, {len(unchanged)} kept")
    return 0


if __name__ == "__main__":
    sys.exit(main())
