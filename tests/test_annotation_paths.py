"""Precursor-related ions, peptidoform parsing edge cases and ion-type selection (#166)."""

from __future__ import annotations

import pytest

from proteomics_mcp.core import annotation
from proteomics_mcp.core.annotation import (
    NEUTRON_MASS,
    PROTON_MASS,
    aa_sequence_from_proforma,
    generate_precursor_related_ions,
    generate_terminal_fragments,
    parse_peptidoform,
    peptidoform_neutral_mass,
)

pytest.importorskip("pyopenms")

# Every modification aa_sequence_from_proforma maps to OpenMS: N-terminal acetyl, carbamidomethyl, deamidation and
# phospho (site-specific), oxidation
MAPPED = "[UNIMOD:1]-PEPM[UNIMOD:35]S[UNIMOD:21]C[UNIMOD:4]N[UNIMOD:7]K"


@pytest.mark.parametrize("peptidoform", ["PEPTIDEK", MAPPED, "PEPT[UNIMOD:21]Y[UNIMOD:21]Q[UNIMOD:7]K/3"])
def test_precursor_ion_mass_agrees_with_the_peptidoform_mass(peptidoform):
    [ion] = generate_precursor_related_ions(peptidoform, precursor_charge=1)

    assert ion["annotation"] == "p+"
    # OpenMS and our own element masses differ in the sixth decimal
    assert ion["mz"] == pytest.approx(peptidoform_neutral_mass(peptidoform) + PROTON_MASS, abs=1e-4)


def test_precursor_ions_with_charges_losses_and_isotopes():
    ions = generate_precursor_related_ions("PEPS[UNIMOD:21]K", precursor_charge=2, max_isotope=1, neutral_losses=True)
    by_annotation = {ion["annotation"]: ion for ion in ions}

    # 2 charges x (no loss, H2O, NH3, and H3PO4 for the phospho peptide) x 2 isotopes
    assert len(ions) == 16
    assert {"p+", "p++", "p+-H3PO4", "p++-H2O[M+1]"} <= set(by_annotation)
    assert by_annotation["p+-H3PO4"]["mz"] == pytest.approx(by_annotation["p+"]["mz"] - 97.976895573)
    assert by_annotation["p++[M+1]"]["mz"] == pytest.approx(by_annotation["p++"]["mz"] + NEUTRON_MASS / 2)
    assert by_annotation["p++-H2O[M+1]"]["neutral_loss"] == "H2O" and by_annotation["p++-H2O[M+1]"]["isotope"] == 1


def test_explicit_neutral_losses_replace_the_standard_ones():
    ions = generate_precursor_related_ions("PEPTIDEK", precursor_charge=1, neutral_losses={"CO2": 43.98983})

    assert [ion["annotation"] for ion in ions] == ["p+", "p+-CO2"]


@pytest.mark.parametrize("peptidoform", ["", "[UNIMOD:1]-", "PEP[UNIMOD:35"])
def test_no_precursor_ions_without_a_sequence(peptidoform):
    assert aa_sequence_from_proforma(peptidoform) is None
    assert generate_precursor_related_ions(peptidoform, precursor_charge=2) == []


def test_without_pyopenms_there_are_no_precursor_ions(monkeypatch):
    import builtins

    real_import = builtins.__import__

    def no_pyopenms(name, *args, **kwargs):
        if name == "pyopenms":
            raise ImportError(name)
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_pyopenms)

    assert aa_sequence_from_proforma("PEPTIDEK") is None


def test_a_sequence_openms_rejects_gives_no_precursor_ions(monkeypatch):
    class Rejecting:
        @staticmethod
        def fromString(_):
            raise RuntimeError("not a sequence OpenMS knows")

    import pyopenms

    monkeypatch.setattr(pyopenms, "AASequence", Rejecting)

    assert aa_sequence_from_proforma("PEPTIDEK") is None


@pytest.mark.parametrize(
    ("peptidoform", "unsupported"),
    [
        ("[Phospho]?PEPSK", True),  # unknown position
        ("PEPS[Phospho", True),  # unclosed bracket
        ("PEPK[NotAModification]", True),
        ("[NotAModification]-PEPK", True),
        ("PEPK-[NotAModification]", True),
        ("K.PEPTIDEK.A", False),  # flanking residues are dropped
        ("PEPTIDEK/2", False),
    ],
)
def test_parse_peptidoform_flags_what_it_cannot_place(monkeypatch, peptidoform, unsupported):
    monkeypatch.setattr(annotation, "unimod_mass", lambda token: None if "NotA" in token or token == "Phospho" else 79.966331)

    parsed = parse_peptidoform(peptidoform)

    assert parsed.unsupported_modification is unsupported
    if not unsupported:
        assert "".join(parsed.residues) == "PEPTIDEK"


@pytest.mark.parametrize(
    ("token", "delta"),
    [("+15.9949", 15.9949), ("-17.0265", -17.0265), ("15.9949", 15.9949), ("", None), ("+abc", None)],
)
def test_mod_delta_reads_masses(monkeypatch, token, delta):
    monkeypatch.setattr(annotation, "unimod_mass", lambda token: None)

    assert annotation._mod_delta(token) == (pytest.approx(delta) if delta is not None else None)


def test_a_single_residue_has_no_terminal_fragments():
    assert generate_terminal_fragments("K", ion_types="by") == []


@pytest.mark.parametrize(
    ("ion_types", "expected"),
    [("b, y", ["b", "y"]), ("b y z", ["b", "y", "z"]), ("by", ["b", "y"]), ("w,d", ["w", "d"]), ("b,q", ["b"])],
)
def test_ion_type_lists_accept_commas_spaces_and_letters(ion_types, expected):
    assert annotation._parse_ion_types(ion_types)[: len(expected)] == expected
