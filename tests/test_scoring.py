from __future__ import annotations

from proteomics_mcp.core.acquisition import acquisition_metadata_summary
from proteomics_mcp.core.scoring import match_predicted_to_observed, precursor_mass_errors, score_candidate_spectrum
from proteomics_mcp.models.schemas import CandidatePSM, ObservedSpectrum, PredictedFragment, PredictedSpectrum


def test_match_uses_closest_peak_within_da_tolerance():
    observed = ObservedSpectrum(mz=[99.99, 100.01, 100.015], intensity=[10, 50, 100])
    fragments = [PredictedFragment(mz=100.0, intensity=1.0, annotation="b1+", ion_type="b", ordinal=1)]

    [match] = match_predicted_to_observed(observed, fragments, tolerance=0.02, tolerance_unit="Da")

    assert match.observed_mz == 100.01
    assert match.observed_intensity == 50


def test_score_candidate_spectrum_returns_broad_metrics():
    candidate = CandidatePSM(candidate_id="c1", source_engine="test", peptidoform="PEPTIDE", precursor_charge=2)
    observed = ObservedSpectrum(mz=[100.0, 200.0, 300.0], intensity=[10.0, 20.0, 30.0], precursor_mz=300.0, precursor_charge=2)
    predicted = PredictedSpectrum(
        peptidoform="PEPTIDE",
        precursor_charge=2,
        predictor_name="mock",
        fragments=[
            PredictedFragment(mz=100.0, intensity=0.5, annotation="b1+", ion_type="b", ordinal=1),
            PredictedFragment(mz=200.0, intensity=1.0, annotation="y1+", ion_type="y", ordinal=1),
        ],
    )

    result = score_candidate_spectrum(candidate=candidate, observed=observed, predicted=predicted)

    assert result.metric_profile == "legacy_v1"
    assert result.matched_ion_count == 2
    assert result.total_predicted_ions == 2
    assert result.matched_ion_fraction == 1.0
    assert result.spectral_angle is not None
    assert result.cosine_similarity is not None
    assert result.xcorr is not None
    assert result.log2_xcorr is not None
    assert result.hyperscore is not None
    assert result.log2_hyperscore is not None
    assert result.precursor_signal_fraction == 0.5
    assert result.precursor_mass_error_ppm is not None
    assert result.annotation_summary["engine"] == "proteomics_mcp_terminal_ion_annotator"
    assert "ion_type_counts" in result.annotation_summary


def test_precursor_mass_error_is_zero_for_theoretical_openms_mz():
    import pyopenms  # type: ignore

    theoretical_mz = float(pyopenms.AASequence.fromString("PEPTIDE").getMZ(2))

    ppm, da = precursor_mass_errors("PEPTIDE", 2, theoretical_mz)

    assert ppm is not None
    assert da is not None
    assert abs(ppm) < 1e-9
    assert abs(da) < 1e-12


def test_score_candidate_uses_precursor_target_metadata_for_precursor_error():
    import pyopenms  # type: ignore

    theoretical_mz = float(pyopenms.AASequence.fromString("PEPTIDE").getMZ(2))
    candidate = CandidatePSM(candidate_id="c1", source_engine="test", peptidoform="PEPTIDE", precursor_charge=2)
    observed = ObservedSpectrum(
        mz=[100.0],
        intensity=[10.0],
        precursor_charge=2,
        metadata={"attributes": [{"accession": "MS:1000827", "name": "isolation window target m/z", "value": theoretical_mz}]},
    )
    predicted = PredictedSpectrum(
        peptidoform="PEPTIDE",
        precursor_charge=2,
        predictor_name="mock",
        fragments=[PredictedFragment(mz=100.0, intensity=1.0, annotation="b1+", ion_type="b", ordinal=1)],
    )

    result = score_candidate_spectrum(candidate=candidate, observed=observed, predicted=predicted)

    assert result.precursor_mass_error_ppm is not None
    assert result.precursor_mass_error_da is not None
    assert abs(result.precursor_mass_error_ppm) < 1e-9
    assert abs(result.precursor_mass_error_da) < 1e-12


def test_acquisition_separates_instrument_model_from_detector_type():
    observed = ObservedSpectrum(
        mz=[100.0],
        intensity=[10.0],
        metadata={
            "attributes": [
                {"accession": "MS:1000463", "name": "instrument", "value": "Orbitrap Fusion Lumos"},
                {"accession": "MS:1000484", "name": "orbitrap"},
            ]
        },
    )

    summary = acquisition_metadata_summary(observed)

    assert summary["instrument_model"] == "Orbitrap Fusion Lumos"
    assert summary["detector_type"] == "orbitrap"


def test_acquisition_does_not_treat_instrument_model_as_detector_type():
    observed = ObservedSpectrum(
        mz=[100.0],
        intensity=[10.0],
        metadata={"attributes": [{"accession": "MS:1002732", "name": "Orbitrap Fusion Lumos"}]},
    )

    summary = acquisition_metadata_summary(observed)

    assert summary["instrument_model"] == "Orbitrap Fusion Lumos"
    assert summary["detector_type"] is None
