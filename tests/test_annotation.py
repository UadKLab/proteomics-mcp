from __future__ import annotations

import csv
from pathlib import Path

from proteomics_mcp.core.annotation import annotate_observed_ions, generate_terminal_fragments, generate_theoretical_ions
from proteomics_mcp.models.schemas import ObservedSpectrum


def test_internal_annotation_counts_terminal_ions():
    fragments = generate_terminal_fragments("PEPTIDE", ion_types="by", max_ion_charge=1)
    b1 = next(fragment for fragment in fragments if fragment["annotation"] == "b1+")
    y1 = next(fragment for fragment in fragments if fragment["annotation"] == "y1+")
    observed = ObservedSpectrum(
        mz=[b1["mz"] + 0.01, y1["mz"] - 0.01, 900.0],
        intensity=[10.0, 20.0, 5.0],
        precursor_charge=2,
    )

    result = annotate_observed_ions(observed, peptidoform="PEPTIDE", tolerance=0.05, ion_types="by")

    assert result["engine"] == "proteomics_mcp_terminal_ion_annotator"
    assert result["annotated_peak_count"] == 2
    assert result["ion_type_counts"] == {"b": 1, "y": 1}
    assert result["annotations"][:2] == ["b1+", "y1+"]
    assert result["matches"][0]["mz_error_unit"] == "Da"


def test_terminal_and_exotic_ions_match_reference_table():
    table = Path(__file__).parent / "fixtures" / "residue_ion_reference.csv"
    rows = list(csv.DictReader(table.read_text(encoding="utf-8").splitlines()))

    for row in rows:
        residue = row["residue"]
        n_terminal = {fragment["annotation"]: fragment for fragment in generate_terminal_fragments(f"{residue}G", ion_types="all", max_ion_charge=1)}
        c_terminal = {fragment["annotation"]: fragment for fragment in generate_terminal_fragments(f"G{residue}", ion_types="all", max_ion_charge=1)}
        expected = {
            "a1+": row["a_1plus"],
            "a1•+": row["a+1_1plus"],
            "b1+": row["b_1plus"],
            "c1+": row["c_1plus"],
            "d1+": row["d_primary_mz_1plus"],
            "d-alt1+": row["d_alt_mz_1plus"],
        }
        for label, value in expected.items():
            if value:
                assert abs(n_terminal[label]["mz"] - float(value)) < 1e-5
        expected = {
            "x1+": row["x_1plus"],
            "y1+": row["y_1plus"],
            "z1+": row["z_1plus"],
            "z1•+": row["z-dot_1plus"],
            "w1+": row["w_primary_mz_1plus"],
            "w-alt1+": row["w_alt_mz_1plus"],
        }
        for label, value in expected.items():
            if value:
                assert abs(c_terminal[label]["mz"] - float(value)) < 1e-5


def test_exotic_annotations_contribute_to_observed_signal():
    d_fragment = next(fragment for fragment in generate_terminal_fragments("IG", ion_types="exotic", max_ion_charge=1) if fragment["annotation"] == "d1+")
    observed = ObservedSpectrum(mz=[d_fragment["mz"], 900.0], intensity=[25.0, 75.0], precursor_charge=2)

    result = annotate_observed_ions(observed, peptidoform="IG", tolerance=0.05, ion_types="exotic")

    assert result["ion_type_counts"] == {"d": 1}
    assert result["annotated_intensity"] == 25.0
    assert result["annotated_intensity_fraction"] == 0.25


def test_modified_terminal_fragments_include_mass_shift():
    unmodified = next(fragment for fragment in generate_terminal_fragments("MPEPTIDE", ion_types="b", max_ion_charge=1) if fragment["annotation"] == "b1+")
    oxidized = next(
        fragment
        for fragment in generate_terminal_fragments("M[UNIMOD:35]PEPTIDE", ion_types="b", max_ion_charge=1)
        if fragment["annotation"] == "b1+"
    )

    assert abs((oxidized["mz"] - unmodified["mz"]) - 15.99491461957) < 1e-6


def test_numeric_bracket_mass_shift_is_supported():
    unmodified = next(fragment for fragment in generate_terminal_fragments("MPEPTIDE", ion_types="b", max_ion_charge=1) if fragment["annotation"] == "b1+")
    shifted = next(fragment for fragment in generate_terminal_fragments("M[15.994915]PEPTIDE", ion_types="b", max_ion_charge=1) if fragment["annotation"] == "b1+")

    assert abs((shifted["mz"] - unmodified["mz"]) - 15.994915) < 1e-6


def test_unimod_catalog_mass_shift_is_supported():
    unmodified = next(fragment for fragment in generate_terminal_fragments("KPEPTIDE", ion_types="b", max_ion_charge=1) if fragment["annotation"] == "b1+")
    shifted = next(fragment for fragment in generate_terminal_fragments("K[UNIMOD:34]PEPTIDE", ion_types="b", max_ion_charge=1) if fragment["annotation"] == "b1+")

    assert abs((shifted["mz"] - unmodified["mz"]) - 14.01565) < 1e-5


def test_c_terminal_ordinals_increase_with_mass():
    fragments = {fragment["annotation"]: fragment for fragment in generate_terminal_fragments("PEPTIDE", ion_types="y", max_ion_charge=1, max_isotope=0)}

    assert fragments["y1+"]["mz"] < fragments["y2+"]["mz"] < fragments["y3+"]["mz"]
    assert fragments["y6+"]["mz"] > fragments["y1+"]["mz"]


def test_auxiliary_ion_types_are_generated_and_matched():
    theoretical = {fragment["annotation"]: fragment for fragment in generate_theoretical_ions("PEPTIDE", ion_types="immonium,reporter,precursor", precursor_charge=2, max_isotope=0)}
    observed = ObservedSpectrum(
        mz=[theoretical["Immonium(P)"]["mz"], theoretical["TMT126"]["mz"], theoretical["p+"]["mz"]],
        intensity=[20.0, 30.0, 50.0],
        precursor_charge=2,
    )

    result = annotate_observed_ions(observed, peptidoform="PEPTIDE", tolerance=0.05, ion_types="immonium,reporter,precursor")

    assert result["annotated_peak_count"] == 3
    assert result["ion_type_counts"] == {"immonium": 1, "reporter": 1, "precursor": 1}
    assert result["annotations"] == ["Immonium(P)", "TMT126", "p+"]
