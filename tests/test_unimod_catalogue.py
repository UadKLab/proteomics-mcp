"""Modification lookups that go to the Unimod catalogue, run against a small stand-in catalogue (#166).

The real catalogue is downloaded by pyteomics on first use; these tests never touch the network.
"""

from __future__ import annotations

import sys

import pytest

from proteomics_mcp.adapters import modifications
from proteomics_mcp.adapters.modifications import (
    bracket_modification_tokens,
    build_peptidoform,
    format_peptidoform_for_display,
    matching_parenthesis,
    normalize_peptidoform_to_unimod_style,
    parse_fragpipe_assigned_modifications,
    resolve_mod_token,
    resolve_unimod_id,
    unimod_mass,
    unresolved_modification_tokens,
)

RARE = {"record_id": 1234, "title": "Rare", "full_name": "A rare modification", "mono_mass": 123.4567, "approved": True}
TWIN = {"record_id": 1235, "title": "Rare2", "full_name": "A rare modification, twin", "mono_mass": 223.4567, "approved": True}
UNAPPROVED = {"record_id": 5555, "title": "Unapproved", "full_name": "Not approved", "mono_mass": 222.2222, "approved": False}
NO_MASS = {"record_id": 7777, "title": "NoMass", "full_name": "No usable mass", "mono_mass": "n/a", "approved": True}
BROKEN = {"record_id": "not a number", "title": "Broken", "mono_mass": 1.0}


class FakeCatalogue:
    """The parts of pyteomics.mass.Unimod the lookups use: mods, by_title and by_name."""

    mods = [RARE, TWIN, UNAPPROVED, NO_MASS, BROKEN]

    def by_title(self, title: str):
        return next((entry for entry in self.mods if entry["title"] == title), None)

    def by_name(self, name: str):
        # pyteomics returns every entry whose name matches; the lookup picks the exact one
        return [entry for entry in self.mods if name.lower() in str(entry.get("full_name", "")).lower()]


class TitleLookupFails(FakeCatalogue):
    def by_title(self, title: str):
        raise KeyError(title)


@pytest.fixture
def catalogue(monkeypatch):
    def use(instance):
        monkeypatch.setattr(modifications, "_load_unimod_catalog", lambda: instance)
        modifications.unimod_mass.cache_clear()
        return instance

    yield use
    modifications.unimod_mass.cache_clear()


def test_a_name_resolves_through_the_catalogue_title(catalogue):
    catalogue(FakeCatalogue())

    assert resolve_unimod_id("Rare") == "UNIMOD:1234"
    assert unimod_mass("Rare") == pytest.approx(123.4567)


def test_a_full_name_resolves_through_by_name_when_the_title_lookup_fails(catalogue):
    catalogue(TitleLookupFails())

    # by_name("A rare modification") matches both entries; the exact full name wins
    assert resolve_unimod_id("A rare modification") == "UNIMOD:1234"


def test_a_maxquant_title_with_sites_falls_back_to_the_title(catalogue):
    catalogue(FakeCatalogue())

    assert resolve_unimod_id("Rare (K)") == "UNIMOD:1234"


def test_an_accession_outside_the_built_in_table_gets_its_mass_from_the_catalogue(catalogue):
    catalogue(FakeCatalogue())

    assert unimod_mass("UNIMOD:1234") == pytest.approx(123.4567)
    assert unimod_mass("UNIMOD:7777") is None  # the entry has no usable mass
    assert unimod_mass("UNIMOD:424242") is None  # not in the catalogue


def test_a_mass_matches_approved_catalogue_entries_only(catalogue):
    catalogue(FakeCatalogue())

    assert resolve_mod_token("+123.4567") == "UNIMOD:1234"
    assert resolve_mod_token("+123.46") == "UNIMOD:1234"  # within MASS_MATCH_TOLERANCE
    assert resolve_mod_token("+222.2222") == "+222.2222"  # unapproved: kept as a mass


def test_a_catalogue_that_cannot_be_listed_leaves_the_built_in_table(catalogue):
    class Unlistable(FakeCatalogue):
        @property
        def mods(self):
            raise RuntimeError("catalogue unavailable")

    catalogue(Unlistable())

    assert resolve_mod_token("+15.9949") == "UNIMOD:35"
    assert resolve_mod_token("+123.4567") == "+123.4567"
    assert resolve_unimod_id("UNIMOD:1234") == "UNIMOD:1234"  # accessions resolve without the catalogue
    assert unimod_mass("UNIMOD:1234") is None


def test_without_a_catalogue_only_built_in_names_resolve(catalogue):
    catalogue(None)

    assert resolve_unimod_id("Oxidation") == "UNIMOD:35"
    assert resolve_unimod_id("Rare") is None
    assert resolve_mod_token("Rare") == "Rare"  # kept as written, and listed as unresolved
    assert unresolved_modification_tokens("PEPK[Rare]") == ["Rare"]


@pytest.mark.parametrize("failure", ["missing pyteomics", "download fails"])
def test_the_catalogue_is_none_when_it_cannot_be_loaded(monkeypatch, failure):
    monkeypatch.setattr(modifications, "_UNIMOD_CATALOG", None)
    if failure == "missing pyteomics":
        monkeypatch.setitem(sys.modules, "pyteomics.mass", None)
    else:
        import pyteomics.mass

        def unreachable():
            raise OSError("unimod.org unreachable")

        monkeypatch.setattr(pyteomics.mass, "Unimod", unreachable)

    assert modifications._load_unimod_catalog() is None


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("PEPM[UNIMOD:35]K", "PEPM(ox)K"),
        ("[UNIMOD:1]-PEPTIDEK", "(acetyl)-PEPTIDEK"),
        ("PEPK[Foo]", "PEPK(FOO)"),
        (None, None),
        ("  ", None),
    ],
)
def test_format_peptidoform_for_display(catalogue, text, expected):
    catalogue(None)
    assert format_peptidoform_for_display(text) == expected


def test_empty_inputs(catalogue):
    catalogue(None)

    assert resolve_mod_token(None) is None and resolve_mod_token(" ") is None
    assert resolve_unimod_id(None) is None and resolve_unimod_id(" ") is None and resolve_unimod_id("UNIMOD:") is None
    assert normalize_peptidoform_to_unimod_style(None) is None and normalize_peptidoform_to_unimod_style(" ") is None
    assert unresolved_modification_tokens("") == []
    assert build_peptidoform(None) is None and build_peptidoform("123") is None


def test_unclosed_brackets_are_kept_as_written():
    assert matching_parenthesis("M(ox", 1) is None
    assert bracket_modification_tokens("PEPM(oxK") == "PEPM(oxK"
    assert bracket_modification_tokens("PEPM[UNIMOD:35K") == "PEPM[UNIMOD:35K"


def test_assigned_modifications_on_the_wrong_residue_are_not_placed():
    # 3 is P, not M: a residue mismatch is not placed on the peptide
    assert parse_fragpipe_assigned_modifications("3M(15.9949), 4M(15.9949)", "PEPMK") == {4: ["15.9949"]}
