from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from proteomics_mcp.adapters.instanovo import _as_float, _as_int, _is_missing, _probability_from_log_probability, normalize_peptidoform_to_instanovo_style
from proteomics_mcp.adapters.table_io import describe_columns, read_candidate_table
from proteomics_mcp.models.schemas import CandidatePSM

# Columns holding Casanovo's own peptide score: its PSM field and its mzTab column
_PEPTIDE_SCORE_COLUMNS = ("peptide_score", "search_engine_score[1]")
# InstaNovo writes log-probabilities a hair above 0 (e.g. 1.9e-08); anything clearly above 0 is not a log-probability
_LOG_PROBABILITY_NOISE = 1e-6


def _get_str(row: pd.Series, column: str) -> str | None:
    if column not in row.index or _is_missing(row[column]):
        return None
    value = str(row[column]).strip()
    return value or None


def _parse_version(casanovo_version: str | None) -> tuple[str, tuple[int, int]] | None:
    if casanovo_version is None:
        return None
    match = re.fullmatch(r"v?(?P<major>\d+)(?:\.(?P<minor>\d+))?(?:\.\d+)*", str(casanovo_version).strip())
    if not match:
        raise ValueError(f"casanovo_version must look like '4.2.1' or '5.2.1', not {casanovo_version!r}.")
    return str(casanovo_version).strip().removeprefix("v"), (int(match.group("major")), int(match.group("minor") or 0))


def _read_peptide_score(score: float, version: tuple[int, int] | None) -> tuple[float | None, str | None, str, list[str]]:
    """Confidence, confidence type, score name and warnings for Casanovo's own peptide score (#157).

    Casanovo 4.x scores a peptide by the mean of its amino-acid scores, 5.x by their product (a sequence probability).
    Up to 5.1 it subtracts 1 when the peptide does not fit the precursor m/z, so a score below 0 marks a mismatch;
    from 5.2 the score stays in [0, 1].
    """
    if version is None:
        name = "Casanovo peptide score"
    elif version < (5, 0):
        name = "Casanovo peptide score (mean amino-acid score)"
    else:
        name = "Casanovo peptide score (amino-acid score product)"
    if score < 0:
        return None, None, name, [
            f"Casanovo peptide score {score:.4g} is below 0, Casanovo's mark (v4.0-5.1) for a peptide that does not fit the "
            "precursor m/z; no confidence given."
        ]
    if score > 1:
        return None, None, name, [f"Casanovo peptide score {score:.4g} is outside Casanovo's range [-1, 1]; no confidence given."]
    if version is None:
        return None, None, name, [
            "Casanovo version unknown: its peptide score is a mean of amino-acid scores before v5.0 and a sequence "
            "probability from v5.0, so no confidence is given. Pass casanovo_version."
        ]
    if version < (5, 0):
        # a mean of amino-acid scores is not the probability of the sequence
        return None, None, name, []
    return score, "casanovo_sequence_probability", name, []


def _score_column(frame: pd.DataFrame) -> str | None:
    """The score column, by what its name declares, never by the sign of its values (#157)."""
    present = [column for column in ("log_probs", *_PEPTIDE_SCORE_COLUMNS) if column in frame.columns]
    if len(present) > 1:
        raise ValueError(f"Casanovo table has several score columns ({', '.join(present)}); keep the one to use.")
    if not present:
        return None
    column = present[0]
    if column == "log_probs":
        values = pd.to_numeric(frame[column], errors="coerce")
        if (values > _LOG_PROBABILITY_NOISE).any():
            raise ValueError(
                f"Column 'log_probs' has values above 0 (largest {values.max():.4g}), so it does not hold natural-log "
                "probabilities. If it holds Casanovo's own peptide score (search_engine_score[1] in its mzTab), name the "
                "column 'peptide_score' and pass casanovo_version."
            )
    return column


def load_casanovo_predictions(path: str | Path, *, limit: int | None = None, casanovo_version: str | None = None) -> list[CandidatePSM]:
    """Load Casanovo prediction exports into canonical candidates.

    The score is read from 'log_probs' (a natural-log sequence probability) or from 'peptide_score' /
    'search_engine_score[1]' (Casanovo's own peptide score, whose meaning depends on casanovo_version).
    """
    frame = read_candidate_table(path, limit=limit)
    if not {"predictions", "prediction"} & set(frame.columns):
        # Without it every row is skipped and the table would load as empty, silently (#140).
        raise ValueError(f"Missing required Casanovo column 'predictions'. Columns found: {describe_columns(frame)}.")
    score_column = _score_column(frame)
    parsed_version = _parse_version(casanovo_version)
    candidates: list[CandidatePSM] = []
    for row_idx, row in frame.iterrows():
        peptidoform = _get_str(row, "predictions") or _get_str(row, "prediction")
        if not peptidoform:
            continue
        peptidoform = normalize_peptidoform_to_instanovo_style(peptidoform)
        raw_score = _as_float(row.get(score_column)) if score_column else None
        confidence, confidence_type, raw_score_name, warnings = None, None, None, []
        if raw_score is not None and score_column == "log_probs":
            confidence = _probability_from_log_probability(raw_score)
            confidence_type, raw_score_name = "sequence_probability_from_log_probability", "log_probs"
        elif raw_score is not None:
            confidence, confidence_type, raw_score_name, warnings = _read_peptide_score(raw_score, parsed_version[1] if parsed_version else None)
        candidate_id = str(row.get("prediction_id")) if "prediction_id" in row.index and not _is_missing(row.get("prediction_id")) else f"{Path(path).stem}:{row_idx}"
        candidates.append(
            CandidatePSM(
                candidate_id=candidate_id,
                source_engine="casanovo",
                source_engine_version=parsed_version[0] if parsed_version else None,
                experiment_name=_get_str(row, "group") or _get_str(row, "experiment_name"),
                scan_number=_as_int(row.get("scan_number")),
                spectrum_id=_get_str(row, "spectrum_id"),
                usi=_get_str(row, "usi"),
                peptidoform=peptidoform,
                precursor_mz=_as_float(row.get("precursor_mz")),
                precursor_charge=_as_int(row.get("precursor_charge")),
                rank=0,
                raw_score=raw_score,
                raw_score_name=raw_score_name,
                next_beam=None,
                next_beam_margin=None,
                confidence=confidence,
                confidence_type=confidence_type,
                source_columns={"schema": "casanovo_predictions", **{str(k): row[k] for k in row.index}},
                warnings=warnings,
            )
        )
    return candidates
