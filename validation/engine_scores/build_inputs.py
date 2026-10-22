"""Build the committed inputs of the engine-score validation: data/bsa1_ms2.mgf and data/search.fasta.

    uv run python validation/engine_scores/build_inputs.py            # check the committed files can be rebuilt
    uv run python validation/engine_scores/build_inputs.py --write    # rebuild them

- data/bsa1_ms2.mgf: the first 300 MS2 spectra with precursor charge 2 or 3 of BSA1.mzML, the OpenMS example BSA run
  (tag Release3.1.0; ion-trap CID fragment spectra), downloaded and checked against its sha256.
- data/search.fasta: the FASTA both engines search. Its human part is the first 300 reviewed human entries UniProt's
  REST API returned on 2026-10-22, which cannot be fetched again exactly, so it is kept as committed; to it are added
  BSA (P02769) and, as their own entries, the peptides of the reference payloads and of the Arabidopsis PSMs, so the
  true peptides are among the candidates.
- data/arabidopsis_usi.mgf is built by fetch_usi_spectra.py.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
import urllib.request
from pathlib import Path

from pyteomics import mzml

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
BSA_MZML = ("https://raw.githubusercontent.com/OpenMS/OpenMS/Release3.1.0/share/OpenMS/examples/BSA/BSA1.mzML",
            "dc9ed61d595328d4ef2f1de47d21f41b83e2eae7c9145e1d9b88e910c8cec2f7")
BSA_SPECTRA = 300
BSA_FASTA = "https://rest.uniprot.org/uniprotkb/P02769.fasta"
CACHE = Path.home() / ".cache" / "proteomics-mcp-e2e"


def download(url: str, sha256: str, name: str) -> Path:
    target = CACHE / name
    if not target.is_file() or hashlib.sha256(target.read_bytes()).hexdigest() != sha256:
        CACHE.mkdir(parents=True, exist_ok=True)
        with urllib.request.urlopen(url, timeout=300) as response:
            target.write_bytes(response.read())
    digest = hashlib.sha256(target.read_bytes()).hexdigest()
    if digest != sha256:
        raise SystemExit(f"{name}: sha256 {digest} != {sha256}")
    return target


def bsa1_mgf() -> str:
    path = download(*BSA_MZML, "openms-Release3.1.0-BSA1.mzML")
    blocks = []
    for spectrum in mzml.read(str(path)):
        if spectrum.get("ms level") != 2:
            continue
        ion = spectrum["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]
        charge = int(ion.get("charge state", 0) or 0)
        if charge not in (2, 3):
            continue
        peaks = "\n".join(f"{m:.6f} {i:.4f}" for m, i in zip(spectrum["m/z array"], spectrum["intensity array"]) if i > 0)
        blocks.append(f"BEGIN IONS\nTITLE=BSA1:{spectrum['id']}\nPEPMASS={float(ion['selected ion m/z']):.6f}\nCHARGE={charge}+\n{peaks}\nEND IONS")
        if len(blocks) == BSA_SPECTRA:
            break
    return "\n\n".join(blocks) + "\n"


def _stripped(peptidoform: str) -> str:
    return re.sub(r"[^A-Z]", "", re.sub(r"\[[^\]]*\]", "", peptidoform))


def search_fasta(current: str) -> str:
    """Keep the committed human entries; rebuild BSA and the peptide entries."""
    human = re.split(r"(?=^>)", current, flags=re.M)
    human = "".join(entry for entry in human if entry.startswith(">sp|") and "OX=9606" in entry)
    with urllib.request.urlopen(BSA_FASTA, timeout=120) as response:
        bsa = response.read().decode()
    peptides = set()
    for path in sorted((REPO / "examples").glob("*_top5/*_payload.json")):
        peptides.add(_stripped(json.loads(path.read_text(encoding="utf-8"))["metrics"]["candidate"]["peptidoform"]))
    tsv = (HERE / "data" / "arabidopsis_usi.tsv").read_text(encoding="utf-8").splitlines()[1:]
    psms = Path.home() / ".cache/proteomics-mcp-e2e/zenodo-23295902/supplement/processed/arabidopsis/supplementary_scored_psms.csv"
    if psms.is_file():
        import csv

        peptides |= {_stripped(row["peptidoform"]) for row in csv.DictReader(psms.open(encoding="utf-8"))}
    elif tsv:
        print("Zenodo supplement not in the cache: Arabidopsis peptides left out", file=sys.stderr)
    entries = "".join(f">pep{index}\n{sequence}\n" for index, sequence in enumerate(sorted(peptides)))
    return human + bsa + entries


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    fasta_path = HERE / "data" / "search.fasta"
    built = {HERE / "data" / "bsa1_ms2.mgf": bsa1_mgf(), fasta_path: search_fasta(fasta_path.read_text(encoding="utf-8") if fasta_path.is_file() else "")}
    status = 0
    for path, text in built.items():
        if args.write:
            path.write_text(text, encoding="utf-8")
            print(f"wrote {path.relative_to(REPO)}")
        elif not path.is_file() or path.read_text(encoding="utf-8") != text:
            print(f"{path.relative_to(REPO)} differs from a rebuild", file=sys.stderr)
            status = 1
        else:
            print(f"{path.relative_to(REPO)}: rebuilt identically")
    return status


if __name__ == "__main__":
    raise SystemExit(main())
