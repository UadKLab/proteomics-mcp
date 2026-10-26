from __future__ import annotations

from dataclasses import asdict, dataclass, field, is_dataclass
from typing import Any, Literal

ToleranceUnit = Literal["ppm", "Da"]
NormalizationMode = Literal["basepeak", "tic", "l2", "none"]


def model_dump(value: Any) -> Any:
    """Recursively convert dataclass model objects to JSON-serializable values."""
    if is_dataclass(value):
        return {k: model_dump(v) for k, v in asdict(value).items()}
    if isinstance(value, list):
        return [model_dump(v) for v in value]
    if isinstance(value, dict):
        return {str(k): model_dump(v) for k, v in value.items()}
    return value


@dataclass(slots=True)
class CandidatePSM:
    candidate_id: str
    source_engine: str
    peptidoform: str
    scan_number: int | None = None
    source_engine_version: str | None = None
    experiment_name: str | None = None
    spectrum_id: str | None = None
    # 0-based position of the spectrum among the run's spectra in the source file, when the engine reports that instead
    # of a scan number (InstaNovo on MGF input)
    spectrum_index: int | None = None
    usi: str | None = None
    precursor_mz: float | None = None
    precursor_charge: int | None = None
    rank: int | None = None
    raw_score: float | None = None
    raw_score_name: str | None = None
    next_beam: float | None = None
    next_beam_margin: float | None = None
    confidence: float | None = None
    token_log_probabilities: list[float] | None = None
    token_log_probability_mean: float | None = None
    token_log_probability_min: float | None = None
    token_confidences: list[float] | None = None
    token_confidence_mean: float | None = None
    token_confidence_min: float | None = None
    token_confidence_min_position: int | None = None
    token_confidence_min_token: str | None = None
    confidence_type: str | None = None
    source_columns: dict[str, Any] = field(default_factory=dict)
    # Problems found while reading this row that may make the peptidoform or its scores wrong, e.g. MaxQuant fixed
    # modifications that could not be applied (#101)
    warnings: list[str] = field(default_factory=list)

    def model_dump(self) -> dict[str, Any]:
        return model_dump(self)


@dataclass(slots=True)
class ObservedSpectrum:
    mz: list[float]
    intensity: list[float]
    spectrum_id: str | None = None
    scan_number: int | None = None
    precursor_mz: float | None = None
    precursor_charge: int | None = None
    rt_seconds: float | None = None
    source: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def model_dump(self) -> dict[str, Any]:
        return model_dump(self)


@dataclass(slots=True)
class PredictedFragment:
    mz: float
    intensity: float
    annotation: str
    ion_type: str | None = None
    ordinal: int | None = None
    charge: int = 1
    neutral_loss: str | None = None

    def model_dump(self) -> dict[str, Any]:
        return model_dump(self)


@dataclass(slots=True)
class PredictedSpectrum:
    peptidoform: str
    precursor_charge: int
    fragments: list[PredictedFragment]
    predictor_name: str
    predictor_version: str | None = None
    model: str | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def model_dump(self) -> dict[str, Any]:
        return model_dump(self)


@dataclass(slots=True)
class MatchedPeak:
    annotation: str
    predicted_mz: float
    observed_mz: float
    mz_error: float
    mz_error_unit: str
    predicted_intensity: float
    observed_intensity: float
    ion_type: str | None = None
    ordinal: int | None = None
    fragment_charge: int = 1
    observed_index: int | None = None
    predicted_index: int | None = None

    def model_dump(self) -> dict[str, Any]:
        return model_dump(self)


@dataclass(slots=True)
class SpectrumMetricResult:
    candidate: CandidatePSM
    spectral_angle: float | None
    cosine_similarity: float | None
    pearson_correlation: float | None
    # Comet XCorr and the X! Tandem hyperscore (core/engine_scores.py); None under legacy_v1, which fills the *_like fields
    xcorr: float | None
    hyperscore: float | None
    matched_ion_count: int
    total_predicted_ions: int
    matched_ion_fraction: float | None
    explained_intensity_fraction: float | None
    precursor_signal_fraction: float | None
    fragment_mass_error_mean_ppm: float | None
    fragment_mass_error_abs_mean_ppm: float | None
    fragment_mass_error_mean_da: float | None
    fragment_mass_error_abs_mean_da: float | None
    precursor_mass_error_ppm: float | None
    precursor_mass_error_da: float | None
    matched_peaks: list[MatchedPeak]
    precursor_isotope_error: int | None = None
    precursor_mass_error_corrected_ppm: float | None = None
    precursor_mass_error_corrected_da: float | None = None
    # legacy_v1 only: the 0.1.x variants (0.1 m/z bins against the predicted spectrum; ln n! + ln sums), kept so earlier
    # scores reproduce. They are not the engines' scores.
    xcorr_like: float | None = None
    log10_xcorr_like: float | None = None
    log2_xcorr_like: float | None = None
    hyperscore_like: float | None = None
    log2_hyperscore_like: float | None = None
    metric_profile: str = "v2"
    predictor_metadata: dict[str, Any] = field(default_factory=dict)
    # What was scored (#137): the aligned predicted/observed intensities behind spectral angle, cosine and Pearson, the
    # predicted ions left out with the reason, the number of values each metric was computed on, and the settings
    scored_ions: list[dict[str, Any]] = field(default_factory=list)
    excluded_predicted_ions: list[dict[str, Any]] = field(default_factory=list)
    metric_n: dict[str, int] = field(default_factory=dict)
    scoring_settings: dict[str, Any] = field(default_factory=dict)
    annotation_settings: dict[str, Any] = field(default_factory=dict)
    annotation_summary: dict[str, Any] = field(default_factory=dict)
    acquisition_provenance: dict[str, Any] = field(default_factory=dict)

    def model_dump(self) -> dict[str, Any]:
        return model_dump(self)
