"""Modification notations of each supported engine, resolved without the online Unimod catalogue (#75)."""

from __future__ import annotations

import pandas as pd
import pytest

from proteomics_mcp.adapters import modifications
from proteomics_mcp.adapters.fragpipe import _peptidoform_from_row
from proteomics_mcp.adapters.maxquant import _peptidoform_from_modified_sequence
from proteomics_mcp.adapters.modifications import normalize_peptidoform_to_unimod_style
from proteomics_mcp.adapters.proteome_discoverer import _peptidoform_from_annotated_sequence
from proteomics_mcp.core.annotation import parse_peptidoform
from proteomics_mcp.core.dataset import _strip_modifications


@pytest.fixture(autouse=True)
def offline_unimod(monkeypatch):
    """Only the built-in modification table resolves: no catalogue download, no cached online results."""
    monkeypatch.setattr(modifications, "_load_unimod_catalog", lambda: None)
    modifications.unimod_mass.cache_clear()
    yield
    modifications.unimod_mass.cache_clear()


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        # InstaNovo
        ("M(+15.99)PEPTIDE", "M[UNIMOD:35]PEPTIDE"),
        ("PEPN(+.98)TIDE", "PEPN[UNIMOD:7]TIDE"),
        ("(+42.01)PEPTIDE", "[UNIMOD:1]-PEPTIDE"),
        ("(+25.98)PEPTIDE", "[+25.98]-PEPTIDE"),
        ("(-17.03)QPEPTIDE", "[UNIMOD:385]-QPEPTIDE"),
        ("C(+57.02)PEPTIDE", "C[UNIMOD:4]PEPTIDE"),
        # Casanovo
        ("+42.011PEPM+15.995K", "[UNIMOD:1]-PEPM[UNIMOD:35]K"),
        ("+43.006-17.027PEPTIDE", "[UNIMOD:5][UNIMOD:385]-PEPTIDE"),
        ("PEPC+57.021K", "PEPC[UNIMOD:4]K"),
        # Bracketed names, masses and UNIMOD ids, unchanged behaviour
        ("C[Carbamidomethyl]M[Oxidation]", "C[UNIMOD:4]M[UNIMOD:35]"),
        ("[Acetyl]PEPTIDE", "[UNIMOD:1]-PEPTIDE"),
        ("M[15.994915]PEPTIDE", "M[UNIMOD:35]PEPTIDE"),
        ("[UNIMOD:1]-PEPTIDE", "[UNIMOD:1]-PEPTIDE"),
        ("PEPTIDE-(-.98)", "PEPTIDE-[UNIMOD:2]"),
    ],
)
def test_normalize_handles_every_notation(text, expected):
    assert normalize_peptidoform_to_unimod_style(text) == expected
    assert not parse_peptidoform(expected).unsupported_modification


@pytest.mark.parametrize(
    ("modified_sequence", "expected"),
    [
        ("_(ac)PEPM(ox)K_", "[UNIMOD:1]-PEPM[UNIMOD:35]K"),
        ("_(Acetyl (Protein N-term))PEPM(Oxidation (M))K_", "[UNIMOD:1]-PEPM[UNIMOD:35]K"),
        ("_PEPS(Phospho (STY))K_", "PEPS[UNIMOD:21]K"),
        ("_PEPTIDEK(Amidated (Peptide C-term))_", "PEPTIDEK-[UNIMOD:2]"),
        ("_PEPK(GlyGly (K))R_", "PEPK[UNIMOD:121]R"),
        ("_AETLM(ox)QARKR_", "AETLM[UNIMOD:35]QARKR"),
        ("_EPEPTIDE_", "EPEPTIDE"),
    ],
)
def test_maxquant_modified_sequence(modified_sequence, expected):
    assert _peptidoform_from_modified_sequence(modified_sequence, None) == expected


@pytest.mark.parametrize(
    ("row", "expected"),
    [
        ({"Peptide": "PEPTMK", "Assigned Modifications": "N-term(42.0106), 5M(15.9949)", "Modified Peptide": "n[43]PEPTM[147]K"},
         "[UNIMOD:1]-PEPTM[UNIMOD:35]K"),
        ({"Peptide": "PEPTMK", "Assigned Modifications": "N-term(42.0106)", "Modified Peptide": "n[43]PEPTMK"}, "[UNIMOD:1]-PEPTMK"),
        ({"Peptide": "PEPTMK", "Assigned Modifications": float("nan"), "Modified Peptide": float("nan")}, "PEPTMK"),
        # Modified Peptide only: total masses, integer-rounded or exact, and delta masses
        ({"Peptide": "PEPTMK", "Modified Peptide": "n[43]PEPTM[147]K"}, "[UNIMOD:1]-PEPTM[UNIMOD:35]K"),
        ({"Peptide": "PEPTIDEK", "Modified Peptide": "n[230]PEPTIDEK[357]"}, "[UNIMOD:737]-PEPTIDEK[UNIMOD:737]"),
        ({"Peptide": "PEPSK", "Modified Peptide": "PEPS[166.998359]K"}, "PEPS[UNIMOD:21]K"),
        ({"Peptide": "MPEPTIDE", "Modified Peptide": "M[15.994915]PEPTIDE"}, "M[UNIMOD:35]PEPTIDE"),
        # Pyro-Glu from Q as a total mass; UNIMOD:385 (ammonia loss) and UNIMOD:28 have the same composition and mass
        ({"Peptide": "QPEPTIDE", "Modified Peptide": "Q[111.032]PEPTIDE"}, "Q[UNIMOD:385]PEPTIDE"),
    ],
)
def test_fragpipe_modifications(row, expected):
    assert _peptidoform_from_row(pd.Series(row)) == expected


@pytest.mark.parametrize(
    ("annotated", "mods", "expected"),
    [
        ("PEPTIDEK", "1xTMT6plex [N-Term]; 1xTMT6plex [K8]", "[UNIMOD:737]-PEPTIDEK[UNIMOD:737]"),
        ("PEPTCDECK", "2xCarbamidomethyl [C5; C8]", "PEPTC[UNIMOD:4]DEC[UNIMOD:4]K"),
        ("PEPSTK", "1xPhospho [S4(99.2)]", "PEPS[UNIMOD:21]TK"),
        ("PEPTIDEK", "1xAcetyl [N-Term(Prot)]", "[UNIMOD:1]-PEPTIDEK"),
        ("PEPTIDEK", "K8(TMTpro); N-Term(TMTpro)", "[UNIMOD:2016]-PEPTIDEK[UNIMOD:2016]"),
        ("[K].mPEPTIDEk.[L]", "N-Term(Acetyl); M1(Oxidation)", "[UNIMOD:1]-M[UNIMOD:35]PEPTIDEK"),
        # Met-loss: the methionine is removed from the sequence once, not again through the modification mass
        ("MPEPTIDEK", "N-Term(Prot)(Met-loss+Acetyl)", "[UNIMOD:1]-PEPTIDEK"),
        ("MPEPTIDEK", "1xMet-loss+Acetyl [N-Term]", "[UNIMOD:1]-PEPTIDEK"),
        # PD numbers sites from the removed methionine, so they shift by one after a Met-loss (#125)
        ("[K].mSYELPAPFKEmVK.[A]", "N-Term(Prot)(Met-loss); M12(Oxidation)", "SYELPAPFKEM[UNIMOD:35]VK"),
        ("[-].mIDKELLDILAcPLcK.[A]", "N-Term(Prot)(Met-loss); C12(Carbamidomethyl); C15(Carbamidomethyl)", "IDKELLDILAC[UNIMOD:4]PLC[UNIMOD:4]K"),
        ("[-].mIDKELLDILAcPLcK.[A]", "1xMet-loss [N-Term]; 2xCarbamidomethyl [C12; C15]", "IDKELLDILAC[UNIMOD:4]PLC[UNIMOD:4]K"),
        ("[-].mSDKELLDILAcPLcK.[A]", "N-Term(Prot)(Met-loss+Acetyl); C12(Carbamidomethyl)", "[UNIMOD:1]-SDKELLDILAC[UNIMOD:4]PLCK"),
        # without a Met-loss, a leading methionine keeps its own site
        ("[K].MSYELPAPFKEMVK.[A]", "M1(Oxidation); M12(Oxidation)", "M[UNIMOD:35]SYELPAPFKEM[UNIMOD:35]VK"),
    ],
)
def test_proteome_discoverer_modifications(annotated, mods, expected):
    assert _peptidoform_from_annotated_sequence(annotated, mods) == expected


@pytest.mark.parametrize(
    "peptidoform",
    [
        _peptidoform_from_annotated_sequence("PEPSTK", "1xPhospho [S/T]"),  # unlocalised
        _peptidoform_from_annotated_sequence("PEPTCDECK", "1xCarbamidomethyl [C4]"),  # site does not match the residue
        _peptidoform_from_annotated_sequence("PEPTIDEK", "M3(Oxidation)"),  # site does not match the residue (#125)
        _peptidoform_from_modified_sequence("_PEPK(NotAModification)R_", None),
        _peptidoform_from_row(pd.Series({"Peptide": "PEPTMK", "Assigned Modifications": "5M(15.9949), garbled"})),
        "(+42.01)PEPTIDE",  # not normalised
        "n[UNIMOD:5]PEPTMK",
        "PEPXIDE",
    ],
)
def test_modifications_that_cannot_be_placed_are_flagged_not_dropped(peptidoform):
    assert parse_peptidoform(peptidoform).unsupported_modification, peptidoform


def test_unlocalised_modification_is_kept_as_unknown_position_token():
    assert _peptidoform_from_annotated_sequence("PEPSTK", "1xPhospho [S/T]") == "[UNIMOD:21]?PEPSTK"


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("PEPN[+.98]TIDE", "PEPNTIDE"),
        ("PEPN(+.98)TIDE", "PEPNTIDE"),
        ("PEPM[Oxidation]K", "PEPMK"),
        ("_PEPM(ox)K_", "PEPMK"),
        ("+42.011PEPM+15.995K", "PEPMK"),
        ("(+42.01)PEPTIDE", "PEPTIDE"),
        ("[UNIMOD:1]-PEPTIDE", "PEPTIDE"),
    ],
)
def test_strip_modifications_handles_every_notation(text, expected):
    assert _strip_modifications(text) == expected
