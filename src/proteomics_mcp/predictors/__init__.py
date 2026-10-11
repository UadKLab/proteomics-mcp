"""Spectrum prediction backends."""

from proteomics_mcp.models.schemas import PredictedSpectrum
from proteomics_mcp.predictors.koina import describe_koina_model, predict_fragments_koina, predict_fragments_koina_batch
from proteomics_mcp.predictors.mock import mock_predicted_spectrum
from proteomics_mcp.predictors.ms2pip import predict_fragments_ms2pip, predict_fragments_ms2pip_batch


def predict_fragments(
    peptidoform: str,
    precursor_charge: int,
    *,
    backend: str,
    model: str,
    max_fragment_charge: int | None = None,
    koina_server_url: str = "https://koina.wilhelmlab.org",
    koina_collision_energy: float = 30.0,
    koina_fragmentation_type: str = "HCD",
    koina_instrument_type: str = "QE",
    koina_timeout_seconds: float = 30.0,
    koina_cache_dir: str | None = None,
) -> PredictedSpectrum:
    if backend == "mock":
        return mock_predicted_spectrum(peptidoform, precursor_charge)
    if backend == "ms2pip":
        return predict_fragments_ms2pip(
            peptidoform,
            precursor_charge,
            model=model,
            max_fragment_charge=max_fragment_charge,
        )
    if backend == "koina":
        return predict_fragments_koina(
            peptidoform,
            precursor_charge,
            model=model,
            max_fragment_charge=max_fragment_charge,
            server_url=koina_server_url,
            collision_energy=koina_collision_energy,
            fragmentation_type=koina_fragmentation_type,
            instrument_type=koina_instrument_type,
            timeout_seconds=koina_timeout_seconds,
            cache_dir=koina_cache_dir,
        )
    raise ValueError(f"Unsupported predictor backend: {backend}")


__all__ = [
    "mock_predicted_spectrum",
    "describe_koina_model",
    "predict_fragments",
    "predict_fragments_koina",
    "predict_fragments_koina_batch",
    "predict_fragments_ms2pip",
    "predict_fragments_ms2pip_batch",
]
