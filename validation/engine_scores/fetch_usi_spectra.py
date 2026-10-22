"""Fetch the Arabidopsis PSM spectra of the manuscript's Figure 4 subset from PRIDE by USI, once.

The result, data/arabidopsis_usi.mgf, is committed so the validation does not depend on PRIDE being up; the USI of each
spectrum is its TITLE, and data/arabidopsis_usi.tsv records a sha256 of every peak list, so anyone can fetch the
spectra again and check they are the same:

    uv run python validation/engine_scores/fetch_usi_spectra.py --check

Input: processed/arabidopsis/supplementary_scored_psms.csv of the Zenodo supplement (10.5281/zenodo.23295902).
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
import time
from pathlib import Path

from proteomics_mcp.core.spectrum_loader import retrieve_spectrum_from_usi

HERE = Path(__file__).resolve().parent
ACCESSION = "PXD014877"


def peak_text(mz: list[float], intensity: list[float]) -> str:
    return "\n".join(f"{m:.6f} {i:.4f}" for m, i in zip(mz, intensity) if i > 0)


def fetch(usi: str, attempts: int = 4):
    for attempt in range(attempts):
        try:
            return retrieve_spectrum_from_usi(usi, timeout_seconds=120)
        except Exception as exc:  # PRIDE PROXI is occasionally slow or unavailable
            if attempt == attempts - 1:
                raise
            print(f"retry {usi}: {exc}", file=sys.stderr)
            time.sleep(10 * (attempt + 1))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--psms", type=Path, default=Path.home() / ".cache/proteomics-mcp-e2e/zenodo-23295902/supplement/processed/arabidopsis/supplementary_scored_psms.csv")
    parser.add_argument("--check", action="store_true", help="fetch again and compare with the committed sha256 values")
    args = parser.parse_args()
    rows = list(csv.DictReader(args.psms.open(encoding="utf-8")))
    manifest = HERE / "data" / "arabidopsis_usi.tsv"
    expected = {row["usi"]: row["peaks_sha256"] for row in csv.DictReader(manifest.open(encoding="utf-8"), delimiter="\t")} if args.check else {}
    blocks, records, seen = [], [], set()
    for row in rows:
        usi = f"mzspec:{ACCESSION}:{row['experiment_name']}:scan:{row['scan_number']}"
        if usi in seen:
            continue
        seen.add(usi)
        spectrum = fetch(usi)
        peaks = peak_text(spectrum.mz, spectrum.intensity)
        digest = hashlib.sha256(peaks.encode()).hexdigest()
        if args.check and expected.get(usi) != digest:
            print(f"MISMATCH {usi}", file=sys.stderr)
            return 1
        charge = int(row["precursor_charge"])
        precursor = spectrum.precursor_mz or float(row["precursor_mz"])
        blocks.append(f"BEGIN IONS\nTITLE={usi}\nPEPMASS={precursor:.6f}\nCHARGE={charge}+\n{peaks}\nEND IONS")
        records.append({"usi": usi, "precursor_mz": f"{precursor:.6f}", "charge": charge, "peaks": peaks.count("\n") + 1, "peaks_sha256": digest})
        print(len(records), usi, flush=True)
    if args.check:
        print(f"all {len(records)} spectra match")
        return 0
    (HERE / "data" / "arabidopsis_usi.mgf").write_text("\n\n".join(blocks) + "\n", encoding="utf-8")
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(records)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
