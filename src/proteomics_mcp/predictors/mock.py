from __future__ import annotations

import re

from proteomics_mcp.models.schemas import PredictedFragment, PredictedSpectrum


def mock_predicted_spectrum(peptidoform: str, precursor_charge: int) -> PredictedSpectrum:
    """Deterministic test and demonstration predictor; not scientific inference."""
    fragments: list[PredictedFragment] = []
    clean = re.sub(r"\[[^\]]+\]", "", peptidoform)
    for ion in ["b", "y"]:
        for ordinal in range(1, max(len(clean), 2)):
            mz = 100.0 + ordinal * 100.0 + (0.2 if ion == "y" else 0.0)
            fragments.append(
                PredictedFragment(
                    mz=mz,
                    intensity=1.0 / ordinal,
                    annotation=f"{ion}{ordinal}+",
                    ion_type=ion,
                    ordinal=ordinal,
                    charge=1,
                )
            )
    return PredictedSpectrum(
        peptidoform=peptidoform,
        precursor_charge=precursor_charge,
        fragments=fragments,
        predictor_name="mock",
        predictor_version="1.0",
        model="mock",
    )
