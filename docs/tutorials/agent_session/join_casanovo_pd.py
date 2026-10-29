"""Attach Proteome Discoverer scan number and retention time to each Casanovo prediction.

Casanovo's `targets` column holds the Proteome Discoverer peptide for the same spectrum, so each
prediction is matched on (target peptidoform, precursor charge, precursor m/z). When a key still
matches several PSMs (the same precursor was fragmented in more than one scan), the scan cannot
be determined from the data; those rows keep all candidate scans and are flagged instead.

Usage: python join_casanovo_pd.py CASANOVO_CSV PD_XLSX OUTPUT_CSV
"""

import re
import sys

import pandas as pd

UNIMOD = {"Oxidation": "UNIMOD:35", "Carbamidomethyl": "UNIMOD:4", "Acetyl": "UNIMOD:1"}


def pd_peptidoform(sequence: str, modifications) -> str:
    """Proteome Discoverer sequence + modification string -> UNIMOD peptidoform.

    Uses the same convention as Casanovo's `targets`: the sequence is kept as Proteome
    Discoverer reports it (including an initial Met flagged as Met-loss) and residue
    modifications stay at Proteome Discoverer's positions.
    """
    residues = list(sequence)
    nterm = ""
    if isinstance(modifications, str):
        for mod in modifications.split("; "):
            site = re.fullmatch(r"([A-Z])(\d+)\((.+)\)", mod)
            if site:
                residues[int(site[2]) - 1] += f"[{UNIMOD[site[3]]}]"
            elif mod == "N-Term(Prot)(Acetyl)":
                nterm = "[UNIMOD:1]-"
            elif mod != "N-Term(Prot)(Met-loss)":
                raise ValueError(f"unhandled modification {mod!r} on {sequence}")
    return nterm + "".join(residues)


def main(casanovo_csv: str, pd_xlsx: str, output_csv: str) -> None:
    cas = pd.read_csv(casanovo_csv)
    psms = pd.read_excel(pd_xlsx)
    psms["pd_peptidoform"] = [
        pd_peptidoform(s, m) for s, m in zip(psms["Sequence"], psms["Modifications"])
    ]
    pd_cols = psms[
        ["pd_peptidoform", "Charge", "m/z [Da]", "First Scan", "RT [min]", "Modifications",
         "Master Protein Accessions", "Percolator q-Value", "Percolator PEP", "XCorr"]
    ].rename(columns={
        "First Scan": "scan_number", "RT [min]": "rt_min", "Modifications": "pd_modifications",
        "Master Protein Accessions": "pd_master_proteins", "Percolator q-Value": "pd_q_value",
        "Percolator PEP": "pd_pep", "XCorr": "pd_xcorr",
    })

    merged = cas.merge(
        pd_cols,
        left_on=["targets", "precursor_charge", "precursor_mz"],
        right_on=["pd_peptidoform", "Charge", "m/z [Da]"],
        how="left",
    ).drop(columns=["pd_peptidoform", "Charge", "m/z [Da]"])

    n_candidates = merged.groupby("prediction_id")["scan_number"].transform("count")
    merged["match_status"] = n_candidates.map(
        lambda n: "unmatched" if n == 0 else ("unique" if n == 1 else "ambiguous")
    )
    candidates = (
        merged.groupby("prediction_id")["scan_number"]
        .agg(lambda s: ";".join(str(int(x)) for x in sorted(s.dropna())))
        .rename("candidate_scans")
    )
    pd_fields = ["scan_number", "rt_min", "pd_modifications", "pd_master_proteins",
                 "pd_q_value", "pd_pep", "pd_xcorr"]
    ambiguous = merged["match_status"] == "ambiguous"
    merged.loc[ambiguous, pd_fields] = pd.NA

    out = (
        merged.drop_duplicates("prediction_id")
        .merge(candidates, on="prediction_id")
        .sort_values("prediction_id")
    )
    out["scan_number"] = out["scan_number"].astype("Int64")
    assert len(out) == len(cas), (len(out), len(cas))
    out.to_csv(output_csv, index=False)

    print(out["match_status"].value_counts().to_string())
    unused = set(psms["First Scan"]) - {
        int(s) for c in out["candidate_scans"] if c for s in c.split(";")
    }
    print(f"Proteome Discoverer PSMs with no Casanovo row: {len(unused)}")


if __name__ == "__main__":
    main(*sys.argv[1:4])
