from __future__ import annotations

import math
from typing import Iterable

import numpy as np


def cosine_similarity(a: list[float], b: list[float]) -> float | None:
    if not a or not b or len(a) != len(b):
        return None
    dot = sum(x * y for x, y in zip(a, b))
    norm_a = math.sqrt(sum(x * x for x in a))
    norm_b = math.sqrt(sum(y * y for y in b))
    if norm_a == 0 or norm_b == 0:
        return None
    return dot / (norm_a * norm_b)


def spectral_angle(a: list[float], b: list[float]) -> float | None:
    cos_val = cosine_similarity(a, b)
    if cos_val is None:
        return None
    clipped = max(-1.0, min(1.0, cos_val))
    return 1.0 - (2.0 * math.acos(clipped) / math.pi)


def pearson_correlation(a: list[float], b: list[float]) -> float | None:
    if len(a) < 2 or len(a) != len(b):
        return None
    arr_a = np.asarray(a, dtype=float)
    arr_b = np.asarray(b, dtype=float)
    if arr_a.std() == 0 or arr_b.std() == 0:
        return None
    return float(np.corrcoef(arr_a, arr_b)[0, 1])


def xcorr_like(
    observed_mz: list[float],
    observed_intensity: list[float],
    predicted_mz: list[float],
    predicted_intensity: list[float],
    *,
    bin_width: float = 0.1,
    max_shift: int = 75,
) -> float | None:
    if not observed_mz or not predicted_mz:
        return None
    max_mz = max(max(observed_mz), max(predicted_mz)) + 1.0
    bins = np.arange(0.0, max_mz + bin_width, bin_width)
    obs_hist = np.zeros(len(bins) - 1, dtype=float)
    pred_hist = np.zeros(len(bins) - 1, dtype=float)
    for mz, intensity in zip(observed_mz, observed_intensity):
        idx = np.searchsorted(bins, mz, side="right") - 1
        if 0 <= idx < len(obs_hist):
            obs_hist[idx] += intensity
    for mz, intensity in zip(predicted_mz, predicted_intensity):
        idx = np.searchsorted(bins, mz, side="right") - 1
        if 0 <= idx < len(pred_hist):
            pred_hist[idx] += intensity
    main = float(np.dot(obs_hist, pred_hist))
    shifted_scores = [float(np.dot(obs_hist, np.roll(pred_hist, shift))) for shift in range(-max_shift, max_shift + 1) if shift != 0]
    return main - (sum(shifted_scores) / len(shifted_scores)) if shifted_scores else main


def hyperscore(matched_count: int, matched_observed: list[float], matched_predicted: list[float]) -> float | None:
    if matched_count == 0:
        return None
    sum_obs = sum(v for v in matched_observed if v > 0)
    sum_pred = sum(v for v in matched_predicted if v > 0)
    if sum_obs <= 0 or sum_pred <= 0:
        return None
    log_factorial = sum(math.log(i) for i in range(1, matched_count + 1))
    return log_factorial + math.log(sum_obs) + math.log(sum_pred)


def mean(values: Iterable[float]) -> float | None:
    values = list(values)
    return sum(values) / len(values) if values else None
