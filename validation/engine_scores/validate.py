"""Reproduce Comet XCorr and the X! Tandem hyperscore with the engines themselves (#51).

proteomics-mcp computes `xcorr` and `hyperscore` (metric profile v2) by reimplementing Comet and X! Tandem
(src/proteomics_mcp/core/engine_scores.py). This script checks that claim against the engines, run from BioContainers
images pinned by digest. For every PSM an engine reports, at any rank and target or decoy, it scores the same
spectrum and peptidoform with proteomics-mcp and compares:

    Comet XCorr        |difference| <= 0.001 (Comet reports 3 decimals)
    X! Tandem score    |difference| <= 0.05  (X! Tandem reports 1 decimal)

Usage (from the repository root):

    uv run python validation/engine_scores/validate.py                 # against the committed engine outputs
    uv run python validation/engine_scores/validate.py --run-engines   # run both engines again (needs Docker)
    uv run python validation/engine_scores/validate.py --run-engines --update-golden

Writes comparison.tsv and summary.json (and agreement plots, with matplotlib) to --out, default
validation/engine_scores/out/. Exit code 1 if any PSM is outside its tolerance, or if the engines' fresh output
differs from the committed one.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

from pyteomics import mgf, pepxml

from proteomics_mcp.core.engine_scores import (
    COMET_VERSION,
    XTANDEM_VERSION,
    CometXcorrParams,
    XTandemHyperscoreParams,
    comet_xcorr,
    xtandem_hyperscore,
)
from proteomics_mcp.models.schemas import ObservedSpectrum

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
COMET_IMAGE = "quay.io/biocontainers/comet-ms:2026011--h9ee0642_0@sha256:af4868654b0652fb55e8b87e8d1da77bed87f92345e4e85895280b4a9e3e9c9c"
XTANDEM_IMAGE = "quay.io/biocontainers/xtandem:15.12.15.2--h4464bbb_11@sha256:712774bdb3afbfd884d22de3cc2e5b93ccaf0074ba527a3b5b2cba748aac10ea"
TOLERANCE = {"comet": 0.001, "xtandem": 0.05}
# The searches' modifications, as the engines were told them (comet_*.params, xtandem_default_input.xml)
MODIFICATIONS = {"C": 57.021464, "M": 15.994915}
GOLDEN = HERE / "golden" / "psms.tsv"
GOLDEN_FIELDS = ["run", "spectrum", "title", "charge", "rank", "peptidoform", "decoy", "engine_score"]


@dataclass(frozen=True)
class Run:
    name: str
    engine: str  # "comet" or "xtandem"
    spectra: str  # "high_res" or "low_res"
    comet_params: str | None = None
    comet_score: CometXcorrParams | None = None
    xtandem_input: dict[str, str] = field(default_factory=dict)
    xtandem_score: XTandemHyperscoreParams | None = None


_WIDE = {"peptide_mass_tolerance": 3.0, "peptide_mass_units": 0, "precursor_tolerance_type": 0}
RUNS = [
    Run("comet_high_res", "comet", "high_res", "comet_high_res.params", CometXcorrParams(**_WIDE)),
    Run("comet_low_res", "comet", "low_res", "comet_low_res.params",
        CometXcorrParams(fragment_bin_tol=1.0005, fragment_bin_offset=0.4, theoretical_fragment_ions=1, **_WIDE)),
    Run("xtandem_da", "xtandem", "high_res", xtandem_score=XTandemHyperscoreParams()),
    Run("xtandem_ppm", "xtandem", "high_res",
        xtandem_input={"spectrum, fragment monoisotopic mass error": "20", "spectrum, fragment monoisotopic mass error units": "ppm"},
        xtandem_score=XTandemHyperscoreParams(fragment_monoisotopic_mass_error=20, fragment_monoisotopic_mass_error_units="ppm")),
    Run("xtandem_reversed", "xtandem", "high_res", xtandem_input={"scoring, include reverse": "only"}, xtandem_score=XTandemHyperscoreParams()),
    Run("xtandem_low_res", "xtandem", "low_res", xtandem_score=XTandemHyperscoreParams()),
]


# --- inputs ----------------------------------------------------------------------------------------------------------


def _peaks(mz, intensity) -> str:
    return "\n".join(f"{m:.6f} {i:.4f}" for m, i in zip(mz, intensity) if i > 0)


def high_res_blocks() -> list[str]:
    """The 25 reference payloads of examples/*_top5, then the Arabidopsis spectra fetched by USI."""
    blocks = []
    for path in sorted((REPO / "examples").glob("*_top5/*_payload.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        observed, candidate = payload["observed"], payload["metrics"]["candidate"]
        charge = candidate.get("precursor_charge") or observed.get("precursor_charge") or payload["predicted"]["precursor_charge"]
        precursor = candidate.get("precursor_mz") or observed.get("precursor_mz")
        title = f"{path.parent.name}/{path.stem}"
        blocks.append(f"BEGIN IONS\nTITLE={title}\nPEPMASS={precursor}\nCHARGE={int(charge)}+\n{_peaks(observed['mz'], observed['intensity'])}\nEND IONS")
    text = (HERE / "data" / "arabidopsis_usi.mgf").read_text(encoding="utf-8")
    blocks += [block.strip() for block in text.split("\n\n") if block.strip()]
    return blocks


def write_spectra(kind: str, target: Path) -> Path:
    if kind == "high_res":
        blocks = high_res_blocks()
    else:
        text = (HERE / "data" / "bsa1_ms2.mgf").read_text(encoding="utf-8")
        blocks = [block.strip() for block in text.split("\n\n") if block.strip()]
    # SCANS numbers the spectra in file order; Comet reports it as start_scan
    numbered = [block.replace("\nPEPMASS=", f"\nSCANS={index}\nPEPMASS=", 1) for index, block in enumerate(blocks, start=1)]
    target.write_text("\n\n".join(numbered) + "\n", encoding="utf-8")
    return target


def read_spectra(path: Path) -> list[tuple[str, ObservedSpectrum]]:
    spectra = []
    for spectrum in mgf.read(str(path), use_index=False):
        params = spectrum["params"]
        observed = ObservedSpectrum(mz=list(spectrum["m/z array"]), intensity=list(spectrum["intensity array"]), precursor_mz=float(params["pepmass"][0]))
        spectra.append((params["title"], observed))
    return spectra


# --- engines ---------------------------------------------------------------------------------------------------------


def _docker(image: str, work: Path, *command: str) -> None:
    user = [] if platform.system() == "Windows" else ["-u", f"{os.getuid()}:{os.getgid()}"]
    # docker reads a relative -v source as the name of a volume, so the path must be absolute
    completed = subprocess.run(["docker", "run", "--rm", *user, "-v", f"{work.resolve()}:/work", "-w", "/work", image, *command],
                               capture_output=True, text=True)
    (work / "engine.log").write_text(completed.stdout + completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise SystemExit(f"{' '.join(command)} failed ({completed.returncode}):\n{completed.stderr[-2000:]}")


def run_comet(run: Run, work: Path) -> Path:
    shutil.copy(HERE / "params" / run.comet_params, work / "comet.params")
    shutil.copy(HERE / "data" / "search.fasta", work / "search.fasta")
    _docker(COMET_IMAGE, work, "comet", "-Pcomet.params", "spectra.mgf")
    return work / "spectra.pep.xml"


def run_xtandem(run: Run, work: Path) -> Path:
    shutil.copy(HERE / "params" / "xtandem_default_input.xml", work / "default_input.xml")
    shutil.copy(HERE / "data" / "search.fasta", work / "search.fasta")
    (work / "taxonomy.xml").write_text(
        '<?xml version="1.0"?>\n<bioml label="x! taxon-to-file matching list">\n  <taxon label="validation">\n'
        '    <file format="peptide" URL="/work/search.fasta" />\n  </taxon>\n</bioml>\n', encoding="utf-8")
    notes = {
        "list path, default parameters": "/work/default_input.xml",
        "list path, taxonomy information": "/work/taxonomy.xml",
        "protein, taxon": "validation",
        "spectrum, path": "/work/spectra.mgf",
        "output, path": "/work/xtandem.xml",
        "output, path hashing": "no",
        "residue, potential modification mass": "15.994915@M",
        "protein, quick acetyl": "no",
        "protein, quick pyrolidone": "no",
        **run.xtandem_input,
    }
    body = "\n".join(f'  <note type="input" label="{label}">{value}</note>' for label, value in notes.items())
    (work / "input.xml").write_text(f'<?xml version="1.0"?>\n<bioml>\n{body}\n</bioml>\n', encoding="utf-8")
    _docker(XTANDEM_IMAGE, work, "xtandem", "/work/input.xml")
    return work / "xtandem.xml"


def _peptidoform(sequence: str, modified: set[int]) -> str:
    """Residues with the searches' modification masses; `modified` holds 0-based positions."""
    out = []
    for index, residue in enumerate(sequence):
        out.append(residue)
        if index in modified:
            if residue not in MODIFICATIONS:
                raise ValueError(f"unexpected modification on {residue} in {sequence}")
            out.append(f"[+{MODIFICATIONS[residue]}]")
    return "".join(out)


def parse_comet(run: Run, path: Path, titles: list[str]) -> list[dict]:
    rows = []
    for query in pepxml.read(str(path)):
        index = int(query["start_scan"])
        for hit in query.get("search_hit", []):
            modified = {mod["position"] - 1 for mod in hit.get("modifications", []) if mod["position"] >= 1}
            rows.append({"run": run.name, "spectrum": index, "title": titles[index - 1], "charge": int(query["assumed_charge"]),
                         "rank": int(hit["hit_rank"]), "peptidoform": _peptidoform(hit["peptide"], modified),
                         "decoy": hit["proteins"][0]["protein"].startswith("DECOY_"), "engine_score": float(hit["search_score"]["xcorr"])})
    return rows


def parse_xtandem(run: Run, path: Path, titles: list[str]) -> list[dict]:
    rows = []
    reversed_only = run.xtandem_input.get("scoring, include reverse") == "only"
    for group in ET.parse(path).getroot().findall("group"):
        if group.get("type") != "model":
            continue
        index = int(group.get("id"))
        seen = set()
        for protein in group.findall("protein"):
            for domain in protein.iter("domain"):
                start = int(domain.get("start"))
                peptidoform = _peptidoform(domain.get("seq"), {int(aa.get("at")) - start for aa in domain.findall("aa")})
                if peptidoform in seen:
                    continue
                seen.add(peptidoform)
                rows.append({"run": run.name, "spectrum": index, "title": titles[index - 1], "charge": int(group.get("z")), "rank": 1,
                             # with "include reverse = only" every sequence is reversed; X! Tandem truncates long labels, ":reversed" too
                             "peptidoform": peptidoform, "decoy": reversed_only or ":reversed" in (protein.get("label") or ""),
                             "engine_score": float(domain.get("hyperscore"))})
    return rows


def run_engines(work_root: Path) -> list[dict]:
    rows = []
    for run in RUNS:
        work = work_root / run.name
        shutil.rmtree(work, ignore_errors=True)
        work.mkdir(parents=True)
        spectra = write_spectra(run.spectra, work / "spectra.mgf")
        titles = [title for title, _ in read_spectra(spectra)]
        print(f"running {run.name} ({len(titles)} spectra)", flush=True)
        if run.engine == "comet":
            rows += parse_comet(run, run_comet(run, work), titles)
        else:
            rows += parse_xtandem(run, run_xtandem(run, work), titles)
    return rows


# --- comparison ------------------------------------------------------------------------------------------------------


def read_golden(path: Path = GOLDEN) -> list[dict]:
    rows = []
    with path.open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append({**row, "spectrum": int(row["spectrum"]), "charge": int(row["charge"]), "rank": int(row["rank"]),
                         "decoy": row["decoy"] == "True", "engine_score": float(row["engine_score"])})
    return rows


def write_golden(rows: list[dict], path: Path = GOLDEN) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=GOLDEN_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in sorted(rows, key=lambda r: (r["run"], r["spectrum"], r["rank"], r["peptidoform"])):
            writer.writerow({key: row[key] for key in GOLDEN_FIELDS})


def rescore(rows: list[dict], work: Path) -> list[dict]:
    """Score every engine PSM with proteomics-mcp; adds `ours` and `difference`."""
    runs = {run.name: run for run in RUNS}
    spectra = {kind: read_spectra(write_spectra(kind, work / f"{kind}.mgf")) for kind in ("high_res", "low_res")}
    out = []
    for row in rows:
        run = runs[row["run"]]
        title, observed = spectra[run.spectra][row["spectrum"] - 1]
        if title != row["title"]:
            raise SystemExit(f"{row['run']}: spectrum {row['spectrum']} is {title!r}, expected {row['title']!r}; the inputs changed")
        if run.engine == "comet":
            ours = comet_xcorr(observed, row["peptidoform"], charge=row["charge"], params=run.comet_score)
        else:
            ours = xtandem_hyperscore(observed, row["peptidoform"], charge=row["charge"], params=run.xtandem_score)
        difference = None if ours is None else abs(ours - row["engine_score"])
        out.append({**row, "ours": ours, "difference": difference,
                    "ok": difference is not None and difference <= TOLERANCE[run.engine] + 1e-9})
    return out


def summarize(rows: list[dict]) -> dict:
    summary = {}
    for run in RUNS:
        mine = [row for row in rows if row["run"] == run.name]
        differences = [row["difference"] for row in mine if row["difference"] is not None]
        summary[run.name] = {
            "engine": f"Comet {COMET_VERSION}" if run.engine == "comet" else f"X! Tandem {XTANDEM_VERSION}",
            "spectra": run.spectra,
            "psms": len(mine),
            "decoy_psms": sum(row["decoy"] for row in mine),
            "exact": sum(1 for d in differences if d < 1e-9),
            "within_tolerance": sum(row["ok"] for row in mine),
            "tolerance": TOLERANCE[run.engine],
            "max_difference": max(differences) if differences else None,
            "not_scored": sum(1 for row in mine if row["ours"] is None),
        }
    return summary


def plot(rows: list[dict], target: Path) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not installed; no plots", file=sys.stderr)
        return
    figure, axes = plt.subplots(1, 2, figsize=(9, 4.2))
    for axis, engine, label in ((axes[0], "comet", "Comet XCorr"), (axes[1], "xtandem", "X! Tandem hyperscore")):
        for run in [r for r in RUNS if r.engine == engine]:
            points = [(row["engine_score"], row["ours"]) for row in rows if row["run"] == run.name and row["ours"] is not None]
            axis.scatter([p[0] for p in points], [p[1] for p in points], s=6, alpha=0.6, label=f"{run.name} (n={len(points)})")
        low, high = axis.get_xlim()
        axis.plot([low, high], [low, high], color="black", linewidth=0.6)
        axis.set_xlabel(f"{label}, engine")
        axis.set_ylabel(f"{label}, proteomics-mcp")
        axis.legend(fontsize=7)
    figure.tight_layout()
    figure.savefig(target, dpi=200)
    figure.savefig(target.with_suffix(".pdf"))


def versions(rows: list[dict]) -> dict:
    return {
        "comet_image": COMET_IMAGE,
        "xtandem_image": XTANDEM_IMAGE,
        "comet_version": COMET_VERSION,
        "xtandem_version": XTANDEM_VERSION,
        "generated": datetime.now(timezone.utc).strftime("%Y-%m-%d"),
        "inputs_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted((HERE / "data").glob("*")) + sorted((HERE / "params").glob("*"))},
        "psms": {run.name: sum(1 for row in rows if row["run"] == run.name) for run in RUNS},
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Check proteomics-mcp's XCorr and hyperscore against Comet and X! Tandem.")
    parser.add_argument("--run-engines", action="store_true", help="run Comet and X! Tandem (Docker) instead of reading golden/psms.tsv")
    parser.add_argument("--update-golden", action="store_true", help="with --run-engines: replace golden/ with this run")
    parser.add_argument("--out", type=Path, default=HERE / "out")
    args = parser.parse_args(argv)
    args.out = args.out.resolve()
    args.out.mkdir(parents=True, exist_ok=True)
    engine_rows = run_engines(args.out / "engines") if args.run_engines else read_golden()
    status = 0
    if args.run_engines:
        fresh = args.out / "psms.tsv"
        write_golden(engine_rows, fresh)
        if args.update_golden:
            shutil.copy(fresh, GOLDEN)
            (HERE / "golden" / "versions.json").write_text(json.dumps(versions(engine_rows), indent=1) + "\n", encoding="utf-8")
        elif fresh.read_bytes() != GOLDEN.read_bytes():
            print("The engines' output differs from golden/psms.tsv", file=sys.stderr)
            status = 1
    rows = rescore(engine_rows, args.out)
    with (args.out / "comparison.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=GOLDEN_FIELDS + ["ours", "difference", "ok"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = summarize(rows)
    (args.out / "summary.json").write_text(json.dumps(summary, indent=1) + "\n", encoding="utf-8")
    plot(rows, args.out / "agreement.png")
    print(f"{'run':18s} {'PSMs':>6s} {'decoys':>6s} {'exact':>6s} {'within':>6s} {'max |diff|':>10s}")
    for name, item in summary.items():
        print(f"{name:18s} {item['psms']:6d} {item['decoy_psms']:6d} {item['exact']:6d} {item['within_tolerance']:6d} {item['max_difference']!s:>10}")
    failed = [row for row in rows if not row["ok"]]
    for row in failed[:20]:
        print(f"FAIL {row['run']} {row['title']} {row['peptidoform']} engine={row['engine_score']} ours={row['ours']}", file=sys.stderr)
    return 1 if failed else status


if __name__ == "__main__":
    raise SystemExit(main())
