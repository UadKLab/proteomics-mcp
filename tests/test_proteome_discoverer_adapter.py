from __future__ import annotations

import pandas as pd

from proteomics_mcp.adapters.proteome_discoverer import load_proteome_discoverer_psms


def test_load_proteome_discoverer_psms_extracts_join_and_score_fields(tmp_path):
    path = tmp_path / "pd.xlsx"
    frame = pd.DataFrame(
        [
            {
                "Annotated Sequence": "[R].KTQTVcNFTDGALVQHQEWDGKESTITR.[K]",
                "Modifications": "C6(Carbamidomethyl)",
                "Spectrum File": "20231127_sample_run.raw",
                "First Scan": 11770,
                "Charge": 4,
                "m/z [Da]": 813.15192,
                "Rank": 1,
                "XCorr": 9.55,
                "DeltaScore": 0.7717,
                "Percolator PEP": 3.739192e-28,
            }
        ]
    )
    frame.to_excel(path, index=False)

    [candidate] = load_proteome_discoverer_psms(path)

    assert candidate.source_engine == "proteome_discoverer"
    assert candidate.experiment_name == "20231127_sample_run"
    assert candidate.scan_number == 11770
    assert candidate.spectrum_id == "20231127_sample_run:11770"
    assert candidate.precursor_charge == 4
    assert candidate.precursor_mz == 813.15192
    assert candidate.peptidoform == "KTQTVC[UNIMOD:4]NFTDGALVQHQEWDGKESTITR"
    assert candidate.raw_score == 9.55
    assert candidate.raw_score_name == "XCorr"
    assert candidate.next_beam is None
    assert candidate.next_beam_margin is None
    assert candidate.source_columns["DeltaScore"] == 0.7717
    assert candidate.confidence is not None and candidate.confidence > 0.999999


def test_load_proteome_discoverer_psms_preserves_multiple_modifications_and_met_loss(tmp_path):
    path = tmp_path / "pd.xlsx"
    frame = pd.DataFrame(
        [
            {
                "Annotated Sequence": "[-].mVHLTPEEKSAVTALWGK.[V]",
                "Modifications": "N-Term(Prot)(Met-loss)",
                "Spectrum File": "run.raw",
                "First Scan": 42,
                "Charge": 2,
                "m/z [Da]": 622.67810,
            },
            {
                "Annotated Sequence": "[K].KIYPTVNcQPLGmISLmK.[R]",
                "Modifications": "C8(Carbamidomethyl); M13(Oxidation); M17(Oxidation)",
                "Spectrum File": "run.raw",
                "First Scan": 43,
                "Charge": 3,
                "m/z [Da]": 709.03687,
            },
        ]
    )
    frame.to_excel(path, index=False)

    candidates = load_proteome_discoverer_psms(path)

    assert candidates[0].peptidoform == "VHLTPEEKSAVTALWGK"
    assert candidates[1].peptidoform == "KIYPTVNC[UNIMOD:4]QPLGM[UNIMOD:35]ISLM[UNIMOD:35]K"


def test_met_loss_sites_are_shifted_and_mismatched_sites_are_reported(tmp_path):
    # PD numbers modification sites from the methionine it reports as lost (#125).
    path = tmp_path / "pd.xlsx"
    rows = [
        ("[K].mSYELPAPFKEmVK.[A]", "N-Term(Prot)(Met-loss); M12(Oxidation)"),
        ("[K].PEPTIDEK.[A]", "M3(Oxidation)"),
    ]
    pd.DataFrame(
        [
            {"Annotated Sequence": seq, "Modifications": mods, "Spectrum File": "run.raw", "First Scan": scan, "Charge": 2, "m/z [Da]": 500.0}
            for scan, (seq, mods) in enumerate(rows, start=1)
        ]
    ).to_excel(path, index=False)

    met_loss, mismatched = load_proteome_discoverer_psms(path)

    assert met_loss.peptidoform == "SYELPAPFKEM[UNIMOD:35]VK"
    assert met_loss.warnings == []
    assert mismatched.peptidoform == "[UNIMOD:35]?PEPTIDEK"
    assert mismatched.warnings == ["Modification 'Oxidation' at site M3 does not match the residue at that position; kept unlocalised."]
