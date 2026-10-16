from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

# Metric profiles: how spectral angle, cosine and Pearson are computed (see core.scoring.similarity_vectors).
METRIC_PROFILES = ("v2", "legacy_v1")
DEFAULT_METRIC_PROFILE = "v2"
# Under metric profile v2, predicted fragments below this fraction of the most intense one are not scored (#83).
DEFAULT_MIN_PREDICTED_INTENSITY = 0.001


@dataclass(frozen=True)
class ProteomicsMCPConfig:
    predictor_backend: str = "ms2pip"
    predictor_model: str = "HCD"
    ms2pip_model: str = "HCD"
    max_fragment_charge: int | None = None
    koina_server_url: str = "https://koina.wilhelmlab.org"
    koina_timeout_seconds: float = 30.0
    koina_collision_energy: float = 30.0
    koina_fragmentation_type: str = "HCD"
    koina_instrument_type: str = "QE"
    koina_batch_size: int = 64
    koina_cache_dir: str | None = None
    tolerance: float = 0.05
    tolerance_unit: str = "Da"
    usi_backend: str = "pride"
    metric_profile: str = DEFAULT_METRIC_PROFILE
    min_predicted_intensity: float = DEFAULT_MIN_PREDICTED_INTENSITY
    annotation_ion_types: str = "by"
    annotation_max_isotope: int = 1
    annotation_max_ion_charge: int | None = None
    annotation_neutral_losses: bool = False

    @classmethod
    def from_mapping(cls, values: dict[str, Any]) -> "ProteomicsMCPConfig":
        normalized = dict(values)
        predictor_backend = str(normalized.get("predictor_backend") or "ms2pip")
        predictor_model = normalized.get("predictor_model")
        ms2pip_model = normalized.get("ms2pip_model")
        if predictor_backend == "koina" and predictor_model in {None, "HCD"} and ms2pip_model in {None, "HCD"}:
            normalized["predictor_model"] = "Prosit_2020_intensity_HCD"
            normalized["ms2pip_model"] = "Prosit_2020_intensity_HCD"
            predictor_model = normalized["predictor_model"]
            ms2pip_model = normalized["ms2pip_model"]
        if predictor_model is None and ms2pip_model is not None:
            normalized["predictor_model"] = ms2pip_model
        elif predictor_model is not None and ms2pip_model is None:
            normalized["ms2pip_model"] = predictor_model
        metric_profile = normalized.get("metric_profile")
        if metric_profile is not None and metric_profile not in METRIC_PROFILES:
            raise ValueError(f"Unknown metric_profile {metric_profile!r}; expected one of {', '.join(METRIC_PROFILES)}.")
        min_intensity = normalized.get("min_predicted_intensity")
        if min_intensity is not None and not 0 <= float(min_intensity) < 1:
            raise ValueError(f"min_predicted_intensity must be in [0, 1), relative to the most intense fragment; got {min_intensity!r}.")
        known = {field: normalized[field] for field in cls.__dataclass_fields__ if field in normalized}
        return cls(**known)


def load_config(path: str | Path | None = None, overrides: dict[str, Any] | None = None) -> ProteomicsMCPConfig:
    values: dict[str, Any] = {}
    if path:
        config_path = Path(path)
        if config_path.suffix.lower() == ".json":
            values.update(json.loads(config_path.read_text(encoding="utf-8")))
        elif config_path.suffix.lower() in {".toml", ".tml"}:
            import tomllib

            values.update(tomllib.loads(config_path.read_text(encoding="utf-8")))
        else:
            raise ValueError("Config files must be JSON or TOML.")
    if overrides:
        values.update({k: v for k, v in overrides.items() if v is not None})
    return ProteomicsMCPConfig.from_mapping(values)
