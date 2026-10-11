from __future__ import annotations

import html as html_lib
import json
import math
import re
from pathlib import Path
from typing import Any

from proteomics_mcp.adapters.modifications import format_peptidoform_for_display
from proteomics_mcp.core.acquisition import acquisition_metadata_summary
from proteomics_mcp.core.annotation import annotate_observed_ions, parse_peptidoform
from proteomics_mcp.core.scoring import precursor_mass_errors
from proteomics_mcp.models.schemas import ObservedSpectrum, PredictedSpectrum, SpectrumMetricResult, model_dump


def build_spectrum_payload(observed: ObservedSpectrum, predicted: PredictedSpectrum, result: SpectrumMetricResult | None = None) -> dict[str, Any]:
    return {
        "observed": observed.model_dump(),
        "predicted": predicted.model_dump(),
        "metrics": result.model_dump() if result else None,
    }


def _peptide_residues(peptidoform: str | None) -> list[str]:
    if not peptidoform:
        return []
    try:
        return parse_peptidoform(peptidoform).residues
    except Exception:
        return []


def _display_peptide_residues(peptidoform: str | None) -> list[str]:
    display = format_peptidoform_for_display(peptidoform)
    if not display:
        return []
    residues = re.findall(r"[A-Z](?:\([^()]+\))*", display)
    return residues or []


def _fragment_sequence_from_annotation(peptidoform: str | None, annotation: Any) -> str | None:
    if not peptidoform or not annotation:
        return None
    residues = _display_peptide_residues(peptidoform)
    if not residues:
        return None
    label = str(annotation)
    if label.startswith("Immonium("):
        match = re.match(r"Immonium\(([A-Z])\)", label)
        return match.group(1) if match else None
    if label.startswith("p+"):
        return "".join(residues)
    if label.startswith("TMT") or label.startswith("iTRAQ"):
        return None
    match = re.match(r"^([A-Za-z][A-Za-z+-]*)(\d+)", label)
    if not match:
        return None
    ion_type = match.group(1).lower()
    ordinal = int(match.group(2))
    if ordinal <= 0 or ordinal > len(residues):
        return None
    if ion_type in {"a", "b", "c", "d", "d-alt"}:
        return "".join(residues[:ordinal])
    if ion_type in {"x", "y", "z", "w", "w-alt"}:
        return "".join(residues[-ordinal:])
    return None


def _plot_controls_markup(sequence_control_id: str | None = None) -> str:
    if not sequence_control_id:
        return ""
    return (
        '<div class="plot-controls">'
        + f'<label class="plot-checkbox"><input type="checkbox" id="{html_lib.escape(sequence_control_id)}" />'
        + "<span>Show fragment ion sequences</span></label>"
        + "</div>"
    )


def _ion_ladder_controls_markup(ion_types: list[str], colors: dict[str, str]) -> str:
    if not ion_types:
        return ""
    controls = []
    for ion_type in ion_types:
        color = colors.get(ion_type, colors["other"])
        escaped_ion_type = html_lib.escape(ion_type)
        controls.append(
            '<label class="plot-checkbox ion-checkbox">'
            f'<input type="checkbox" data-ion-type="{escaped_ion_type}" />'
            f'<span class="ion-swatch" style="--ion-color:{html_lib.escape(color)}"></span>'
            f"<span>{escaped_ion_type}</span>"
            "</label>"
        )
    return '<div class="ion-ladder-controls" aria-label="Annotated ion series">' + "".join(controls) + "</div>"


def _sequence_toggle_script(control_id: str) -> str:
    return f"""
    function formatSequenceLabel(baseLabel, fragmentSequence, showSequence) {{
      if (!showSequence || !fragmentSequence) {{
        return baseLabel || '';
      }}
      if (!baseLabel) {{
        return fragmentSequence;
      }}
      return `${{baseLabel}} | ${{fragmentSequence}}`;
    }}
    function applySequenceToggle(gd, showSequence) {{
      traces.forEach((trace, index) => {{
        const meta = trace.meta || {{}};
        if (!Array.isArray(meta.baseText) || !Array.isArray(meta.fragmentSequenceText)) {{
          return;
        }}
        const nextText = meta.baseText.map((label, textIndex) =>
          formatSequenceLabel(label, meta.fragmentSequenceText[textIndex], showSequence)
        );
        trace.text = nextText;
        Plotly.restyle(gd, {{text: [nextText]}}, [index]);
      }});
    }}
    const sequenceToggle = document.getElementById({json.dumps(control_id)});
    if (sequenceToggle) {{
      sequenceToggle.addEventListener('change', function(event) {{
        applySequenceToggle(document.getElementById('plot'), Boolean(event.target.checked));
      }});
    }}
    """.strip()


def _ion_ladder_toggle_script() -> str:
    return """
    document.querySelectorAll('.ion-ladder-controls input[data-ion-type]').forEach((control) => {
      control.addEventListener('change', function(event) {
        const ionType = event.target.dataset.ionType;
        const visible = Boolean(event.target.checked);
        const indices = traces
          .map((trace, index) => ((trace.meta || {}).ionType === ionType ? index : null))
          .filter((index) => index !== null);
        if (indices.length) {
          Plotly.restyle(document.getElementById('plot'), {visible: visible}, indices);
        }
      });
    });
    """.strip()


def render_spectrum_comparison_html(payload: dict[str, Any], output_path: str | Path | None = None) -> str:
    """Render a self-contained interactive mirror plot with Plotly image download enabled."""
    observed = payload["observed"]
    predicted = payload["predicted"]
    metrics = payload.get("metrics") or {}
    matched_peaks = metrics.get("matched_peaks") or []
    obs_x = observed.get("mz", [])
    obs_y = observed.get("intensity", [])
    fragments = predicted.get("fragments", [])
    peptidoform = (metrics.get("candidate") or {}).get("peptidoform") or predicted.get("peptidoform")
    max_observed = max((abs(float(i)) for i in obs_y), default=1.0) or 1.0
    max_predicted = max((abs(float(f.get("intensity", 0.0))) for f in fragments), default=1.0) or 1.0

    def normalized_observed(value: float) -> float:
        return abs(float(value)) / max_observed

    def normalized_predicted(value: float) -> float:
        return -abs(float(value)) / max_predicted

    matched_observed_by_index = {m.get("observed_index"): m for m in matched_peaks if m.get("observed_index") is not None}
    matched_predicted_by_index = {m.get("predicted_index"): m for m in matched_peaks if m.get("predicted_index") is not None}

    obs_unmatched_x: list[float] = []
    obs_unmatched_y: list[float] = []
    obs_unmatched_raw_y: list[float] = []
    matched_observed_groups: dict[str, dict[str, list[Any]]] = {}
    for idx, (mz, intensity) in enumerate(zip(obs_x, obs_y)):
        raw_intensity = float(intensity)
        y = normalized_observed(raw_intensity)
        match = matched_observed_by_index.get(idx)
        if match is None:
            obs_unmatched_x.append(mz)
            obs_unmatched_y.append(y)
            obs_unmatched_raw_y.append(raw_intensity)
            continue
        ion_type = str(match.get("ion_type") or "other")
        group = matched_observed_groups.setdefault(ion_type, {"x": [], "y": [], "raw_intensity": [], "text": [], "sequence": []})
        group["x"].append(mz)
        group["y"].append(y)
        group["raw_intensity"].append(raw_intensity)
        annotation = match.get("annotation", "")
        group["text"].append(annotation)
        group["sequence"].append(_fragment_sequence_from_annotation(peptidoform, annotation))

    pred_unmatched_x: list[float] = []
    pred_unmatched_y: list[float] = []
    pred_unmatched_raw_y: list[float] = []
    pred_unmatched_text: list[str] = []
    pred_unmatched_sequence: list[str | None] = []
    matched_predicted_groups: dict[str, dict[str, list[Any]]] = {}
    for idx, fragment in enumerate(fragments):
        mz = float(fragment.get("mz", 0.0))
        raw_intensity = float(fragment.get("intensity", 0.0))
        intensity = normalized_predicted(raw_intensity)
        annotation = fragment.get("annotation", "")
        match = matched_predicted_by_index.get(idx)
        if match is None:
            pred_unmatched_x.append(mz)
            pred_unmatched_y.append(intensity)
            pred_unmatched_raw_y.append(raw_intensity)
            pred_unmatched_text.append(annotation)
            pred_unmatched_sequence.append(_fragment_sequence_from_annotation(peptidoform, annotation))
            continue
        ion_type = str(match.get("ion_type") or fragment.get("ion_type") or "other")
        group = matched_predicted_groups.setdefault(ion_type, {"x": [], "y": [], "raw_intensity": [], "text": [], "sequence": []})
        group["x"].append(mz)
        group["y"].append(intensity)
        group["raw_intensity"].append(raw_intensity)
        group["text"].append(annotation)
        group["sequence"].append(_fragment_sequence_from_annotation(peptidoform, annotation))

    observed_palette = {
        "b": "#0072B2",
        "y": "#D55E00",
        "a": "#009E73",
        "c": "#CC79A7",
        "x": "#56B4E9",
        "z": "#E69F00",
        "precursor": "#F0E442",
        "other": "#666666",
    }
    predicted_palette = {
        "b": "#B8336A",
        "y": "#B8336A",
        "a": "#B8336A",
        "c": "#B8336A",
        "x": "#B8336A",
        "z": "#B8336A",
        "precursor": "#B8336A",
        "other": "#B8336A",
    }

    def make_sticks(x_values: list[Any], y_values: list[Any]) -> tuple[list[Any], list[Any]]:
        xs: list[Any] = []
        ys: list[Any] = []
        for x, y in zip(x_values, y_values):
            xs.extend([x, x, None])
            ys.extend([0, y, None])
        return xs, ys

    def repeat_stick_values(values: list[Any]) -> list[Any]:
        return [value for value in values for _ in (0, 1, 2)]

    def hover_scaffold_values(
        x_values: list[Any],
        y_values: list[Any],
        raw_intensities: list[Any],
        labels: list[Any] | None = None,
        fragment_sequences: list[Any] | None = None,
        *,
        steps: int = 9,
    ) -> tuple[list[Any], list[Any], list[Any], list[Any], list[Any]]:
        xs: list[Any] = []
        ys: list[Any] = []
        customdata: list[Any] = []
        text: list[Any] = []
        sequence_text: list[Any] = []
        labels = labels or [""] * len(x_values)
        fragment_sequences = fragment_sequences or [""] * len(x_values)
        for x, y, raw, label, fragment_sequence in zip(x_values, y_values, raw_intensities, labels, fragment_sequences):
            for step in range(1, steps + 1):
                fraction = step / steps
                xs.append(x)
                ys.append(float(y) * fraction)
                customdata.append(raw)
                text.append(label)
                sequence_text.append(fragment_sequence)
        return xs, ys, customdata, text, sequence_text

    def hover_scaffold_trace(
        *,
        name: str,
        x_values: list[Any],
        y_values: list[Any],
        raw_intensities: list[Any],
        labels: list[Any] | None,
        fragment_sequences: list[Any] | None,
        intensity_label: str,
    ) -> str:
        hover_x, hover_y, hover_customdata, hover_text, hover_sequence_text = hover_scaffold_values(x_values, y_values, raw_intensities, labels, fragment_sequences)
        return json.dumps(
            {
                "x": hover_x,
                "y": hover_y,
                "customdata": hover_customdata,
                "text": hover_text,
                "meta": {
                    "baseText": hover_text,
                    "fragmentSequenceText": hover_sequence_text,
                },
                "type": "scatter",
                "mode": "markers",
                "name": f"{name} hover area",
                "marker": {"color": "rgba(0,0,0,0)", "size": 12},
                "showlegend": False,
                "hovertemplate": "m/z=%{x:.5f}<br>" + intensity_label + "=%{customdata:.3e}<br>%{text}<extra>" + name + "</extra>",
            }
        )

    def matched_traces(name: str, group: dict[str, list[Any]], sign: str) -> list[str]:
        palette = observed_palette if sign == "observed" else predicted_palette
        color = palette.get(name, palette["other"])
        trace_name = f"{name}-ions observed" if sign == "observed" else f"{name}-ions predicted"
        stick_x, stick_y = make_sticks(group["x"], group["y"])
        text_position = "top center" if sign == "observed" else "bottom center"
        intensity_label = "observed intensity" if sign == "observed" else "predicted intensity"
        line_trace = json.dumps(
            {
                "x": stick_x,
                "y": stick_y,
                "customdata": repeat_stick_values(group["raw_intensity"]),
                "text": repeat_stick_values(group["text"]),
                "meta": {
                    "baseText": repeat_stick_values(group["text"]),
                    "fragmentSequenceText": repeat_stick_values(group["sequence"]),
                },
                "type": "scatter",
                "mode": "lines",
                "name": trace_name,
                "line": {"color": color, "width": 1.55},
                "showlegend": True,
                "hovertemplate": "m/z=%{x:.5f}<br>" + intensity_label + "=%{customdata:.3e}<br>%{text}<extra>" + trace_name + "</extra>",
            }
        )
        hover_trace = hover_scaffold_trace(
            name=trace_name,
            x_values=group["x"],
            y_values=group["y"],
            raw_intensities=group["raw_intensity"],
            labels=group["text"],
            fragment_sequences=group["sequence"],
            intensity_label=intensity_label,
        )
        label_trace = json.dumps(
            {
                "x": group["x"],
                "y": group["y"],
                "text": group["text"],
                "meta": {
                    "baseText": group["text"],
                    "fragmentSequenceText": group["sequence"],
                },
                "type": "scatter",
                "mode": "text",
                "name": f"{trace_name} labels",
                "textfont": {"color": "#27313f", "size": 11},
                "textposition": text_position,
                "cliponaxis": False,
                "showlegend": False,
                "hoverinfo": "skip",
            }
        )
        if sign == "predicted":
            return [line_trace, hover_trace]
        return [line_trace, hover_trace, label_trace]

    obs_unmatched_stick_x, obs_unmatched_stick_y = make_sticks(obs_unmatched_x, obs_unmatched_y)
    pred_unmatched_stick_x, pred_unmatched_stick_y = make_sticks(pred_unmatched_x, pred_unmatched_y)
    traces = [
        json.dumps(
            {
                "x": obs_unmatched_stick_x,
                "y": obs_unmatched_stick_y,
                "customdata": repeat_stick_values(obs_unmatched_raw_y),
                "type": "scatter",
                "mode": "lines",
                "name": "Observed unmatched",
                "line": {"color": "#9AA6B2", "width": 1.05},
                "hovertemplate": "m/z=%{x:.5f}<br>observed intensity=%{customdata:.3e}<extra>observed unmatched</extra>",
            }
        )
    ]
    if obs_unmatched_x:
        traces.append(
            hover_scaffold_trace(
                name="observed unmatched",
                x_values=obs_unmatched_x,
                y_values=obs_unmatched_y,
                raw_intensities=obs_unmatched_raw_y,
                labels=None,
                fragment_sequences=None,
                intensity_label="observed intensity",
            )
        )
    for ion_type, group in sorted(matched_observed_groups.items()):
        traces.extend(matched_traces(ion_type, group, "observed"))
    traces.append(
        json.dumps(
            {
                "x": pred_unmatched_stick_x,
                "y": pred_unmatched_stick_y,
                "customdata": repeat_stick_values(pred_unmatched_raw_y),
                "text": repeat_stick_values(pred_unmatched_text),
                "meta": {
                    "baseText": repeat_stick_values(pred_unmatched_text),
                    "fragmentSequenceText": repeat_stick_values(pred_unmatched_sequence),
                },
                "type": "scatter",
                "mode": "lines",
                "name": "Predicted unmatched",
                "line": {"color": "#B8336A", "width": 1.05},
                "hovertemplate": "m/z=%{x:.5f}<br>predicted intensity=%{customdata:.3e}<br>%{text}<extra>predicted unmatched</extra>",
            }
        )
    )
    if pred_unmatched_x:
        traces.append(
            hover_scaffold_trace(
                name="predicted unmatched",
                x_values=pred_unmatched_x,
                y_values=pred_unmatched_y,
                raw_intensities=pred_unmatched_raw_y,
                labels=pred_unmatched_text,
                fragment_sequences=pred_unmatched_sequence,
                intensity_label="predicted intensity",
            )
        )
    for ion_type, group in sorted(matched_predicted_groups.items()):
        traces.extend(matched_traces(ion_type, group, "predicted"))

    def metric_value(value: Any, *, percent: bool = False, decimals: int = 3) -> str:
        if value is None:
            return "NA"
        try:
            number = float(value)
        except (TypeError, ValueError):
            return html_lib.escape(str(value))
        if not math.isfinite(number):
            return "NA"
        if percent:
            return f"{number * 100:.1f}%"
        if abs(number) >= 1000 or (0 < abs(number) < 0.001):
            return f"{number:.2e}"
        return f"{number:.{decimals}f}"

    def metric_tile(label: str, value: str, note: str = "") -> str:
        escaped_label = html_lib.escape(label)
        escaped_value = html_lib.escape(value)
        escaped_note = html_lib.escape(note)
        note_html = f'<div class="metric-note">{escaped_note}</div>' if note else ""
        return f'<div class="metric-tile"><div class="metric-label">{escaped_label}</div><div class="metric-value">{escaped_value}</div>{note_html}</div>'

    def numeric_list(values: Any) -> list[float]:
        if not isinstance(values, list):
            return []
        numbers: list[float] = []
        for value in values:
            try:
                number = float(value)
            except (TypeError, ValueError):
                continue
            if math.isfinite(number):
                numbers.append(number)
        return numbers

    def metric_group(title: str, tiles: list[str]) -> str:
        if not tiles:
            return ""
        escaped_title = html_lib.escape(title)
        return f'<section class="metric-group"><h2>{escaped_title}</h2><div class="metric-grid">{"".join(tiles)}</div></section>'

    def ion_count_summary(counts: Any) -> str:
        if not isinstance(counts, dict):
            return "NA"
        parts = []
        for ion, count in sorted(counts.items()):
            try:
                number = int(count)
            except (TypeError, ValueError):
                continue
            if number > 0:
                parts.append(f"{ion}:{number}")
        return ", ".join(parts) if parts else "NA"

    metric_groups: list[str] = []
    metadata_panel = ""
    if metrics:
        candidate = metrics.get("candidate") or {}
        acquisition = _normalized_acquisition_summary(
            observed,
            predicted,
            (metrics.get("acquisition_provenance") or {}).get("summary") or {},
        )
        metrics.setdefault("acquisition_provenance", {})["summary"] = acquisition
        precursor_target_mz = _first_number(
            acquisition.get("isolation_target_mz"),
            candidate.get("precursor_mz"),
            predicted.get("precursor_mz"),
        )
        precursor_charge = _first_int(candidate.get("precursor_charge"), observed.get("precursor_charge"), predicted.get("precursor_charge"))
        precursor_error_ppm = metrics.get("precursor_mass_error_ppm")
        precursor_error_da = metrics.get("precursor_mass_error_da")
        if (precursor_error_ppm is None or precursor_error_da is None) and candidate.get("peptidoform"):
            fallback_ppm, fallback_da = precursor_mass_errors(candidate["peptidoform"], precursor_charge, precursor_target_mz)
            precursor_error_ppm = precursor_error_ppm if precursor_error_ppm is not None else fallback_ppm
            precursor_error_da = precursor_error_da if precursor_error_da is not None else fallback_da
        metadata_panel = _prediction_metadata_panel(
            candidate=candidate,
            observed=observed,
            predicted=predicted,
            acquisition=acquisition,
            precursor_target_mz=precursor_target_mz,
            precursor_charge=precursor_charge,
        )
        matched = metrics.get("matched_ion_count")
        total = metrics.get("total_predicted_ions")
        matched_value = "NA" if matched is None or total is None else f"{matched} / {total}"
        metric_groups.append(
            metric_group(
                "Match Quality",
                [
                    metric_tile("Spectral angle", metric_value(metrics.get("spectral_angle")), "spectrum similarity"),
                    metric_tile("Cosine", metric_value(metrics.get("cosine_similarity")), "matched-ion intensities"),
                    metric_tile("Pearson", metric_value(metrics.get("pearson_correlation")), "matched-ion intensities"),
                    metric_tile("log2 XCorr", metric_value(metrics.get("log2_xcorr")), "0.1 m/z bins"),
                    metric_tile("log2 Hyperscore", metric_value(metrics.get("log2_hyperscore")), "matched-ion evidence"),
                ],
            )
        )
        metric_groups.append(
            metric_group(
                "Ion Evidence",
                [
                    metric_tile("Matched ions", matched_value, "observed / predicted"),
                    metric_tile("Matched fraction", metric_value(metrics.get("matched_ion_fraction"), percent=True), "of predicted ions"),
                    metric_tile("Explained intensity", metric_value(metrics.get("explained_intensity_fraction"), percent=True), "of observed intensity"),
                    metric_tile("Precursor signal", metric_value(metrics.get("precursor_signal_fraction"), percent=True), "within +/-1.5 m/z"),
                    metric_tile("Observed ion series", ion_count_summary((metrics.get("annotation_summary") or {}).get("ion_type_counts")), "annotated peaks"),
                ],
            )
        )
        metric_groups.append(
            metric_group(
                "Mass Error",
                [
                    metric_tile("Fragment abs ppm", metric_value(metrics.get("fragment_mass_error_abs_mean_ppm"))),
                    metric_tile("Fragment abs Da", metric_value(metrics.get("fragment_mass_error_abs_mean_da"))),
                    metric_tile("Precursor ppm", metric_value(precursor_error_ppm), "target minus theoretical"),
                    metric_tile("Precursor Da", metric_value(precursor_error_da), "target minus theoretical"),
                ],
            )
        )
        candidate_tiles = [
            metric_tile("Confidence", metric_value(candidate.get("confidence"), percent=True), candidate.get("source_engine") or "candidate"),
            metric_tile("Raw score", metric_value(candidate.get("raw_score")), candidate.get("raw_score_name") or ""),
            metric_tile("Next beam", metric_value(candidate.get("next_beam")), "second-ranked sequence"),
            metric_tile("Next beam margin", metric_value(candidate.get("next_beam_margin")), "best minus second-ranked"),
        ]
        if candidate.get("token_confidence_mean") is not None:
            candidate_tiles.append(metric_tile("Token confidence mean", metric_value(candidate.get("token_confidence_mean"), percent=True), "exp(logp)"))
        if candidate.get("token_confidence_min") is not None:
            token = candidate.get("token_confidence_min_token")
            position = candidate.get("token_confidence_min_position")
            note = f"{token} at {position}" if token is not None and position is not None else "lowest token"
            candidate_tiles.append(metric_tile("Token confidence min", metric_value(candidate.get("token_confidence_min"), percent=True), note))
        metric_groups.append(metric_group("Candidate", candidate_tiles))
    metrics_panel = f'<div id="metrics-panel" class="metrics-panel">{metadata_panel}{"".join(metric_groups)}</div>' if metric_groups else ""

    html = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8" />
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <title>Spectrum comparison</title>
  <style>
    body {{ margin: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; color: #27313f; background: #FFFFFF; }}
    body.dark {{ color: #E5EDF6; background: #0B1020; }}
    .plot-shell {{ position: relative; width: 100%; height: 720px; }}
    .plot-controls {{ position: absolute; left: 80px; bottom: 92px; z-index: 20; display: flex; align-items: center; gap: 12px; padding: 6px 10px; border: 1px solid #D8E0EA; border-radius: 8px; background: rgba(248, 250, 252, 0.92); backdrop-filter: blur(4px); box-sizing: border-box; }}
    .ion-ladder-controls {{ position: absolute; left: 80px; right: 30px; bottom: 132px; z-index: 20; display: flex; flex-wrap: wrap; align-items: center; gap: 7px 12px; padding: 6px 10px; border: 1px solid #D8E0EA; border-radius: 8px; background: rgba(248, 250, 252, 0.92); backdrop-filter: blur(4px); box-sizing: border-box; }}
    .plot-checkbox {{ display: inline-flex; align-items: center; gap: 8px; font-size: 13px; color: #475569; user-select: none; }}
    .plot-checkbox input {{ width: 14px; height: 14px; margin: 0; accent-color: #068D9D; }}
    .ion-checkbox {{ gap: 5px; font-size: 12px; }}
    .ion-swatch {{ display: inline-block; width: 9px; height: 9px; border-radius: 50%; background: var(--ion-color); }}
    .metrics-panel {{ max-width: 1180px; margin: 0 auto 34px auto; padding: 0 22px 10px 22px; }}
    .metadata-summary {{ border: 1px solid #D8E0EA; border-radius: 8px; padding: 14px; margin: 0 0 18px 0; background: #F8FAFC; }}
    .metadata-summary h2 {{ font-size: 14px; font-weight: 650; letter-spacing: 0; margin: 0 0 10px 0; color: #334155; }}
    .metadata-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); gap: 10px 14px; }}
    .metadata-item {{ min-width: 0; }}
    .metadata-label {{ font-size: 11px; color: #64748B; line-height: 1.2; margin-bottom: 3px; }}
    .metadata-value {{ font-size: 14px; font-weight: 650; color: #111827; line-height: 1.25; overflow-wrap: anywhere; word-break: break-word; }}
    .metric-group {{ margin: 0 0 18px 0; }}
    .metric-group h2 {{ font-size: 14px; font-weight: 650; letter-spacing: 0; margin: 0 0 8px 0; color: #334155; }}
    .metric-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 10px; }}
    .metric-tile {{ border: 1px solid #D8E0EA; border-left: 4px solid #068D9D; border-radius: 8px; padding: 10px 12px; background: #F8FAFC; min-height: 70px; box-sizing: border-box; }}
    .metric-label {{ font-size: 12px; color: #64748B; line-height: 1.2; margin-bottom: 5px; }}
    .metric-value {{ font-size: 22px; font-weight: 720; color: #111827; line-height: 1.15; font-variant-numeric: tabular-nums; overflow-wrap: anywhere; word-break: break-word; }}
    .metric-note {{ font-size: 11px; color: #64748B; line-height: 1.2; margin-top: 5px; overflow-wrap: anywhere; }}
    body.dark .metric-group h2 {{ color: #CBD5E1; }}
    body.dark .plot-checkbox {{ color: #CBD5E1; }}
    body.dark .plot-checkbox input {{ accent-color: #9FD356; }}
    body.dark .plot-controls, body.dark .ion-ladder-controls {{ border-color: #243044; background: rgba(17, 24, 39, 0.92); }}
    body.dark .metadata-summary {{ border-color: #243044; background: #111827; }}
    body.dark .metadata-summary h2 {{ color: #CBD5E1; }}
    body.dark .metadata-label {{ color: #94A3B8; }}
    body.dark .metadata-value {{ color: #F8FAFC; }}
    body.dark .metric-tile {{ border-color: #243044; border-left-color: #9FD356; background: #111827; }}
    body.dark .metric-label, body.dark .metric-note {{ color: #94A3B8; }}
    body.dark .metric-value {{ color: #F8FAFC; }}
    @media print {{
      #plot {{ height: 92vh !important; }}
      body {{ margin: 0; }}
      .modebar, .plot-controls, .ion-ladder-controls {{ display: none !important; }}
    }}
  </style>
</head>
<body>
  <div class="plot-shell">
    <div id="plot" style="width:100%;height:720px;"></div>
    {_plot_controls_markup("fragment-sequence-toggle")}
  </div>
  {metrics_panel}
  <script>
    const traces = [{", ".join(traces)}];
    const colorThemes = {{
      light: {{
        observedUnmatched: '#9AA6B2',
        predictedUnmatched: '#B8336A',
        observed: {{b: '#0072B2', y: '#D55E00', a: '#009E73', c: '#CC79A7', x: '#56B4E9', z: '#E69F00', precursor: '#F0E442', other: '#666666'}},
        predicted: {{b: '#B8336A', y: '#B8336A', a: '#B8336A', c: '#B8336A', x: '#B8336A', z: '#B8336A', precursor: '#B8336A', other: '#B8336A'}},
        label: '#27313f',
        text: '#27313f',
        zero: '#4D6A7F',
        hoverBg: 'rgba(255, 255, 255, 0.96)',
        hoverBorder: '#9AA6B2',
        plot: '#FFFFFF',
        paper: '#FFFFFF'
      }},
      dark: {{
        observedUnmatched: '#64748B',
        predictedUnmatched: '#9792E3',
        observed: {{b: '#2DD4BF', y: '#F59E0B', a: '#34D399', c: '#F472B6', x: '#38BDF8', z: '#FBBF24', precursor: '#FDE047', other: '#94A3B8'}},
        predicted: {{b: '#9792E3', y: '#9792E3', a: '#9792E3', c: '#9792E3', x: '#9792E3', z: '#9792E3', precursor: '#9792E3', other: '#9792E3'}},
        label: '#E5EDF6',
        text: '#E5EDF6',
        zero: '#6B7F99',
        hoverBg: 'rgba(17, 24, 39, 0.96)',
        hoverBorder: '#94A3B8',
        plot: '#111827',
        paper: '#0B1020'
      }}
    }};
    function ionFromTraceName(name) {{
      const match = name.match(/^([a-z]+|precursor)-ions (observed|predicted)/);
      return match ? {{ion: match[1], side: match[2]}} : null;
    }}
    function applyTraceTheme(themeName) {{
      const theme = colorThemes[themeName];
      traces.forEach((trace) => {{
        if (trace.name === 'Observed unmatched') {{
          trace.line.color = theme.observedUnmatched;
        }} else if (trace.name === 'Predicted unmatched') {{
          trace.line.color = theme.predictedUnmatched;
        }} else {{
          const ion = ionFromTraceName(trace.name || '');
          if (ion && trace.line) {{
            trace.line.color = theme[ion.side][ion.ion] || theme[ion.side].other;
          }}
        }}
        if (trace.mode === 'text' && trace.textfont) {{
          trace.textfont.color = theme.label;
        }}
      }});
    }}
    const storedTheme = localStorage.getItem('proteomics-mcp-theme');
    let currentTheme = storedTheme || (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
    applyTraceTheme(currentTheme);
    document.body.classList.toggle('dark', currentTheme === 'dark');
    const layout = {{
      title: 'Observed vs predicted spectrum',
      xaxis: {{title: 'm/z', zeroline: false, showgrid: false}},
      yaxis: {{
        title: 'Relative intensity; observed positive, predicted mirrored below zero',
        exponentformat: 'e',
        tickformat: '.2f',
        separatethousands: false,
        showgrid: false,
        zeroline: true,
        zerolinecolor: colorThemes[currentTheme].zero,
        zerolinewidth: 1
      }},
      bargap: 0.02,
      legend: {{orientation: 'h', y: -0.38}},
      margin: {{l: 80, r: 30, t: 60, b: 210}},
      plot_bgcolor: colorThemes[currentTheme].plot,
      paper_bgcolor: colorThemes[currentTheme].paper,
      hovermode: 'closest',
      hoverlabel: {{
        bgcolor: colorThemes[currentTheme].hoverBg,
        bordercolor: colorThemes[currentTheme].hoverBorder,
        font: {{family: '-apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif', color: colorThemes[currentTheme].text, size: 12}},
        namelength: -1
      }},
      font: {{family: '-apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif', color: colorThemes[currentTheme].text}}
    }};
    function updateTheme(gd, themeName) {{
      currentTheme = themeName;
      localStorage.setItem('proteomics-mcp-theme', themeName);
      document.body.classList.toggle('dark', themeName === 'dark');
      const theme = colorThemes[themeName];
      const lineUpdates = {{}};
      const textUpdates = {{}};
      traces.forEach((trace, index) => {{
        let color = null;
        if (trace.name === 'Observed unmatched') {{
          color = theme.observedUnmatched;
        }} else if (trace.name === 'Predicted unmatched') {{
          color = theme.predictedUnmatched;
        }} else {{
          const ion = ionFromTraceName(trace.name || '');
          if (ion && trace.line) {{
            color = theme[ion.side][ion.ion] || theme[ion.side].other;
          }}
        }}
        if (color !== null) {{
          lineUpdates[index] = color;
        }}
        if (trace.mode === 'text') {{
          textUpdates[index] = theme.label;
        }}
      }});
      Object.entries(lineUpdates).forEach(([index, color]) => Plotly.restyle(gd, {{'line.color': color}}, [Number(index)]));
      Object.entries(textUpdates).forEach(([index, color]) => Plotly.restyle(gd, {{'textfont.color': color}}, [Number(index)]));
      Plotly.relayout(gd, {{
        'plot_bgcolor': theme.plot,
        'paper_bgcolor': theme.paper,
        'font.color': theme.text,
        'hoverlabel.bgcolor': theme.hoverBg,
        'hoverlabel.bordercolor': theme.hoverBorder,
        'hoverlabel.font.color': theme.text,
        'xaxis.showgrid': false,
        'yaxis.showgrid': false,
        'yaxis.zerolinecolor': theme.zero
      }});
    }}
    Plotly.newPlot('plot', traces, layout, {{
      responsive: true,
      scrollZoom: true,
      displaylogo: false,
      modeBarButtons: [
        [
          'toImage',
          {{
            name: 'Print page / Save PDF',
            icon: Plotly.Icons.disk,
            click: function() {{ window.print(); }}
          }}
        ],
        ['zoom2d', 'pan2d'],
        ['zoomIn2d', 'zoomOut2d', 'autoScale2d'],
        [
        {{
          name: 'Toggle dark mode',
          icon: Plotly.Icons.pencil,
          click: function(gd) {{
            updateTheme(gd, currentTheme === 'dark' ? 'light' : 'dark');
          }}
        }}
        ]
      ],
      toImageButtonOptions: {{format: 'png', filename: 'spectrum_comparison', scale: 2}}
    }});
    {_sequence_toggle_script("fragment-sequence-toggle")}
  </script>
  <script type="application/json" id="proteomics-mcp-payload">{json.dumps(payload)}</script>
</body>
</html>"""
    if output_path is not None:
        Path(output_path).write_text(html, encoding="utf-8")
    return html


def _as_float_list(values: Any) -> list[float]:
    numbers: list[float] = []
    if not isinstance(values, list):
        return numbers
    for value in values:
        try:
            number = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(number):
            numbers.append(number)
    return numbers


def _stick_values(x_values: list[float], y_values: list[float]) -> tuple[list[Any], list[Any]]:
    xs: list[Any] = []
    ys: list[Any] = []
    for x, y in zip(x_values, y_values):
        xs.extend([x, x, None])
        ys.extend([0, y, None])
    return xs, ys


def _repeat_stick_values(values: list[Any]) -> list[Any]:
    return [value for value in values for _ in (0, 1, 2)]


def _hover_scaffold_trace(
    *,
    name: str,
    x_values: list[Any],
    y_values: list[Any],
    raw_intensities: list[Any],
    labels: list[Any] | None = None,
    fragment_sequences: list[Any] | None = None,
    intensity_label: str = "observed intensity",
    ion_type: str | None = None,
    steps: int = 9,
) -> dict[str, Any]:
    hover_x: list[Any] = []
    hover_y: list[Any] = []
    hover_customdata: list[Any] = []
    hover_text: list[Any] = []
    hover_sequence_text: list[Any] = []
    labels = labels or [""] * len(x_values)
    fragment_sequences = fragment_sequences or [""] * len(x_values)
    for x, y, raw, label, fragment_sequence in zip(x_values, y_values, raw_intensities, labels, fragment_sequences):
        for step in range(1, steps + 1):
            hover_x.append(x)
            hover_y.append(float(y) * (step / steps))
            hover_customdata.append(raw)
            hover_text.append(label)
            hover_sequence_text.append(fragment_sequence)
    metadata = {
        "baseText": hover_text,
        "fragmentSequenceText": hover_sequence_text,
    }
    if ion_type is not None:
        metadata["ionType"] = ion_type
    return {
        "x": hover_x,
        "y": hover_y,
        "customdata": hover_customdata,
        "text": hover_text,
        "meta": metadata,
        "type": "scatter",
        "mode": "markers",
        "name": f"{name} hover area",
        "marker": {"color": "rgba(0,0,0,0)", "size": 12},
        "showlegend": False,
        "hovertemplate": "m/z=%{x:.5f}<br>" + intensity_label + "=%{customdata:.3e}<br>%{text}<extra>" + name + "</extra>",
    }


def _metric_value(value: Any, *, percent: bool = False, decimals: int = 3) -> str:
    if value is None:
        return "NA"
    try:
        number = float(value)
    except (TypeError, ValueError):
        return html_lib.escape(str(value))
    if not math.isfinite(number):
        return "NA"
    if percent:
        return f"{number * 100:.1f}%"
    if decimals == 0:
        return f"{number:.0f}"
    if abs(number) >= 1000 or (0 < abs(number) < 0.001):
        return f"{number:.2e}"
    return f"{number:.{decimals}f}"


def _metric_tile(label: str, value: str, note: str = "") -> str:
    escaped_label = html_lib.escape(label)
    escaped_value = html_lib.escape(value)
    escaped_note = html_lib.escape(note)
    note_html = f'<div class="metric-note">{escaped_note}</div>' if note else ""
    return f'<div class="metric-tile"><div class="metric-label">{escaped_label}</div><div class="metric-value">{escaped_value}</div>{note_html}</div>'


def _metric_group(title: str, tiles: list[str]) -> str:
    if not tiles:
        return ""
    return f'<section class="metric-group"><h2>{html_lib.escape(title)}</h2><div class="metric-grid">{"".join(tiles)}</div></section>'


def _ion_count_summary(counts: Any) -> str | None:
    if not isinstance(counts, dict):
        return None
    parts: list[str] = []
    for ion, count in sorted(counts.items()):
        try:
            number = int(count)
        except (TypeError, ValueError):
            continue
        if number > 0:
            parts.append(f"{ion}:{number}")
    return ", ".join(parts) if parts else None


def _scan_from_usi(usi: Any) -> str | None:
    if not usi:
        return None
    match = re.search(r":scan:(\d+)", str(usi))
    return match.group(1) if match else None


def _first_number(*values: Any) -> float | None:
    for value in values:
        try:
            number = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(number):
            return number
    return None


def _first_int(*values: Any) -> int | None:
    for value in values:
        try:
            number = int(value)
        except (TypeError, ValueError):
            continue
        if number > 0:
            return number
    return None


def _metadata_item(label: str, value: Any) -> str:
    if value is None or value == "":
        value = "NA"
    return f'<div class="metadata-item"><div class="metadata-label">{html_lib.escape(label)}</div><div class="metadata-value">{html_lib.escape(str(value))}</div></div>'


def _retention_minutes_display(seconds: Any) -> str | None:
    value = _first_number(seconds)
    if value is None:
        return None
    return f"{value / 60.0:.2f} min"


def _prediction_metadata_panel(
    *,
    candidate: dict[str, Any],
    observed: dict[str, Any],
    predicted: dict[str, Any],
    acquisition: dict[str, Any],
    precursor_target_mz: float | None,
    precursor_charge: int | None,
    include_sequence: bool = True,
    include_predictor: bool = True,
) -> str:
    peptidoform = candidate.get("peptidoform") or predicted.get("peptidoform")
    source_engine = candidate.get("source_engine") or "candidate"
    source_version = candidate.get("source_engine_version")
    source = f"{source_engine} {source_version}".strip() if source_version else source_engine
    predictor = predicted.get("predictor_name")
    model = predicted.get("model")
    predictor_label = " / ".join(str(part) for part in (predictor, model) if part)
    scan_number = candidate.get("scan_number") or observed.get("scan_number") or _scan_from_usi(candidate.get("usi"))
    spectrum_id = candidate.get("spectrum_id") or observed.get("spectrum_id")
    isolation_target = f"{precursor_target_mz:.5f}" if precursor_target_mz is not None else None
    fragmentation_method = acquisition.get("fragmentation_method") or predictor_label or "NA"
    items: list[str] = []
    if include_sequence:
        items.append(_metadata_item("Predicted sequence", peptidoform))
    items.extend(
        [
            _metadata_item("Source", source),
            _metadata_item("Fragmentation method", fragmentation_method),
            _metadata_item("Isolation target m/z", isolation_target),
            _metadata_item("Precursor charge", precursor_charge),
            _metadata_item("Retention time", _retention_minutes_display(acquisition.get("retention_time_seconds"))),
            _metadata_item("Collision energy", acquisition.get("collision_energy")),
        ]
    )
    if acquisition.get("scan_header"):
        items.append(_metadata_item("Scan header", acquisition.get("scan_header")))
    if acquisition.get("scan_resolution"):
        items.append(_metadata_item("Scan resolution", acquisition.get("scan_resolution")))
    if acquisition.get("instrument_model"):
        items.append(_metadata_item("Instrument", acquisition.get("instrument_model")))
    if acquisition.get("detector_type"):
        items.append(_metadata_item("Detector", acquisition.get("detector_type")))
    items.extend(
        [
            _metadata_item("Scan", scan_number),
            _metadata_item("Spectrum USI", spectrum_id),
        ]
    )
    if include_predictor:
        items.append(_metadata_item("Predictor", predictor_label or "NA"))
    return (
        '<section class="metadata-summary">'
        '<h2>Prediction Metadata</h2>'
        f'<div class="metadata-grid">{"".join(items)}</div>'
        "</section>"
    )


def _normalized_acquisition_summary(
    observed: dict[str, Any],
    predicted: dict[str, Any] | None = None,
    existing: dict[str, Any] | None = None,
) -> dict[str, Any]:
    extracted = acquisition_metadata_summary(observed, predicted)
    summary = {**extracted, **(existing or {})}
    if not summary.get("instrument_model") and extracted.get("instrument_model"):
        summary["instrument_model"] = extracted["instrument_model"]
    if extracted.get("detector_type"):
        summary["detector_type"] = extracted["detector_type"]
    elif summary.get("detector_type") and summary.get("detector_type") == summary.get("instrument_model"):
        summary["detector_type"] = None
    return summary


def _basic_plot_html(
    *,
    traces: list[dict[str, Any]],
    layout_title: str,
    yaxis_title: str,
    payload: dict[str, Any],
    output_path: str | Path | None,
    metrics_panel: str = "",
    controls_html: str | None = None,
    extra_script: str = "",
) -> str:
    html = f"""<!doctype html>
<html>
<head>
  <meta charset="utf-8" />
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <title>{html_lib.escape(layout_title)}</title>
  <style>
    body {{ margin: 0; font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif; color: #27313f; background: #FFFFFF; }}
    body.dark {{ color: #E5EDF6; background: #0B1020; }}
    .plot-shell {{ position: relative; width: 100%; height: 720px; }}
    .plot-controls {{ position: absolute; left: 80px; bottom: 92px; z-index: 20; display: flex; align-items: center; gap: 12px; padding: 6px 10px; border: 1px solid #D8E0EA; border-radius: 8px; background: rgba(248, 250, 252, 0.92); backdrop-filter: blur(4px); box-sizing: border-box; }}
    .ion-ladder-controls {{ position: absolute; left: 80px; right: 30px; bottom: 132px; z-index: 20; display: flex; flex-wrap: wrap; align-items: center; gap: 7px 12px; padding: 6px 10px; border: 1px solid #D8E0EA; border-radius: 8px; background: rgba(248, 250, 252, 0.92); backdrop-filter: blur(4px); box-sizing: border-box; }}
    .plot-checkbox {{ display: inline-flex; align-items: center; gap: 8px; font-size: 13px; color: #475569; user-select: none; }}
    .plot-checkbox input {{ width: 14px; height: 14px; margin: 0; accent-color: #068D9D; }}
    .ion-checkbox {{ gap: 5px; font-size: 12px; }}
    .ion-swatch {{ display: inline-block; width: 9px; height: 9px; border-radius: 50%; background: var(--ion-color); }}
    .metrics-panel {{ max-width: 1180px; margin: 0 auto 34px auto; padding: 0 22px 10px 22px; }}
    .metadata-summary {{ border: 1px solid #D8E0EA; border-radius: 8px; padding: 14px; margin: 0 0 18px 0; background: #F8FAFC; }}
    .metadata-summary h2 {{ font-size: 14px; font-weight: 650; letter-spacing: 0; margin: 0 0 10px 0; color: #334155; }}
    .metadata-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(190px, 1fr)); gap: 10px 14px; }}
    .metadata-item {{ min-width: 0; }}
    .metadata-label {{ font-size: 11px; color: #64748B; line-height: 1.2; margin-bottom: 3px; }}
    .metadata-value {{ font-size: 14px; font-weight: 650; color: #111827; line-height: 1.25; overflow-wrap: anywhere; word-break: break-word; }}
    .metric-group {{ margin: 0 0 18px 0; }}
    .metric-group h2 {{ font-size: 14px; font-weight: 650; letter-spacing: 0; margin: 0 0 8px 0; color: #334155; }}
    .metric-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 10px; }}
    .metric-tile {{ border: 1px solid #D8E0EA; border-left: 4px solid #068D9D; border-radius: 8px; padding: 10px 12px; background: #F8FAFC; min-height: 70px; box-sizing: border-box; }}
    .metric-label {{ font-size: 12px; color: #64748B; line-height: 1.2; margin-bottom: 5px; }}
    .metric-value {{ font-size: 22px; font-weight: 720; color: #111827; line-height: 1.15; font-variant-numeric: tabular-nums; overflow-wrap: anywhere; word-break: break-word; }}
    .metric-note {{ font-size: 11px; color: #64748B; line-height: 1.2; margin-top: 5px; overflow-wrap: anywhere; }}
    body.dark .metric-group h2 {{ color: #CBD5E1; }}
    body.dark .plot-checkbox {{ color: #CBD5E1; }}
    body.dark .plot-checkbox input {{ accent-color: #9FD356; }}
    body.dark .plot-controls, body.dark .ion-ladder-controls {{ border-color: #243044; background: rgba(17, 24, 39, 0.92); }}
    body.dark .metadata-summary {{ border-color: #243044; background: #111827; }}
    body.dark .metadata-summary h2 {{ color: #CBD5E1; }}
    body.dark .metadata-label {{ color: #94A3B8; }}
    body.dark .metadata-value {{ color: #F8FAFC; }}
    body.dark .metric-tile {{ border-color: #243044; border-left-color: #9FD356; background: #111827; }}
    body.dark .metric-label, body.dark .metric-note {{ color: #94A3B8; }}
    body.dark .metric-value {{ color: #F8FAFC; }}
    @media print {{
      #plot {{ height: 92vh !important; }}
      body {{ margin: 0; }}
      .modebar, .plot-controls, .ion-ladder-controls {{ display: none !important; }}
    }}
  </style>
</head>
<body>
  <div class="plot-shell">
    <div id="plot" style="width:100%;height:720px;"></div>
    {controls_html if controls_html is not None else _plot_controls_markup()}
  </div>
  {metrics_panel}
  <script>
    const traces = {json.dumps(traces)};
    const light = {{
      text: '#27313f',
      zero: '#4D6A7F',
      plot: '#FFFFFF',
      paper: '#FFFFFF',
      hoverBg: 'rgba(255, 255, 255, 0.96)',
      hoverBorder: '#A6B3C2'
    }};
    const dark = {{
      text: '#E5EDF6',
      zero: '#6B7F99',
      plot: '#111827',
      paper: '#0B1020',
      hoverBg: 'rgba(17, 24, 39, 0.96)',
      hoverBorder: '#94A3B8'
    }};
    const storedTheme = localStorage.getItem('proteomics-mcp-theme');
    let currentTheme = storedTheme || (window.matchMedia && window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
    document.body.classList.toggle('dark', currentTheme === 'dark');
    function theme() {{ return currentTheme === 'dark' ? dark : light; }}
    function applyLabelTheme() {{
      traces.forEach((trace) => {{
        if (trace.mode === 'text') {{
          trace.textfont = trace.textfont || {{}};
          trace.textfont.color = theme().text;
        }}
      }});
    }}
    applyLabelTheme();
    const layout = {{
      title: {json.dumps(layout_title)},
      xaxis: {{title: 'm/z', zeroline: false, showgrid: false}},
      yaxis: {{title: {json.dumps(yaxis_title)}, tickformat: '.2f', showgrid: false, zeroline: true, zerolinecolor: theme().zero, zerolinewidth: 1}},
      legend: {{orientation: 'h', y: -0.38}},
      margin: {{l: 80, r: 30, t: 60, b: 210}},
      plot_bgcolor: theme().plot,
      paper_bgcolor: theme().paper,
      hovermode: 'closest',
      hoverlabel: {{
        bgcolor: theme().hoverBg,
        bordercolor: theme().hoverBorder,
        font: {{family: '-apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif', color: theme().text, size: 12}},
        namelength: -1
      }},
      font: {{family: '-apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif', color: theme().text}}
    }};
    function updateTheme(gd, themeName) {{
      currentTheme = themeName;
      localStorage.setItem('proteomics-mcp-theme', themeName);
      document.body.classList.toggle('dark', themeName === 'dark');
      const textUpdates = {{}};
      traces.forEach((trace, index) => {{
        if (trace.mode === 'text') {{
          textUpdates[index] = theme().text;
        }}
      }});
      Object.entries(textUpdates).forEach(([index, color]) => Plotly.restyle(gd, {{'textfont.color': color}}, [Number(index)]));
      Plotly.relayout(gd, {{
        'plot_bgcolor': theme().plot,
        'paper_bgcolor': theme().paper,
        'font.color': theme().text,
        'hoverlabel.bgcolor': theme().hoverBg,
        'hoverlabel.bordercolor': theme().hoverBorder,
        'hoverlabel.font.color': theme().text,
        'yaxis.zerolinecolor': theme().zero
      }});
    }}
    Plotly.newPlot('plot', traces, layout, {{
      responsive: true,
      scrollZoom: true,
      displaylogo: false,
      modeBarButtons: [
        [
          'toImage',
          {{
            name: 'Print page / Save PDF',
            icon: Plotly.Icons.disk,
            click: function() {{ window.print(); }}
          }}
        ],
        ['zoom2d', 'pan2d'],
        ['zoomIn2d', 'zoomOut2d', 'autoScale2d'],
        [{{
          name: 'Toggle dark mode',
          icon: Plotly.Icons.pencil,
          click: function(gd) {{ updateTheme(gd, currentTheme === 'dark' ? 'light' : 'dark'); }}
        }}]
      ],
      toImageButtonOptions: {{format: 'png', filename: 'spectrum', scale: 2}}
    }});
    {extra_script}
  </script>
  <script type="application/json" id="proteomics-mcp-payload">{json.dumps(payload)}</script>
</body>
</html>"""
    if output_path is not None:
        Path(output_path).write_text(html, encoding="utf-8")
    return html


def observed_spectrum_metrics(observed: dict[str, Any]) -> dict[str, Any]:
    mz = _as_float_list(observed.get("mz", []))
    intensity = _as_float_list(observed.get("intensity", []))
    paired = list(zip(mz, intensity))
    total = sum(intensity)
    base_mz = None
    base_intensity = None
    if paired:
        base_mz, base_intensity = max(paired, key=lambda pair: pair[1])
    precursor_mz = observed.get("precursor_mz")
    precursor_region = None
    try:
        precursor = float(precursor_mz)
    except (TypeError, ValueError):
        precursor = None
    if precursor is not None and total > 0:
        precursor_region = sum(i for m, i in paired if abs(m - precursor) <= 1.5) / total
    sorted_intensity = sorted(intensity)
    median_intensity = None
    if sorted_intensity:
        mid = len(sorted_intensity) // 2
        if len(sorted_intensity) % 2:
            median_intensity = sorted_intensity[mid]
        else:
            median_intensity = (sorted_intensity[mid - 1] + sorted_intensity[mid]) / 2.0
    acquisition = acquisition_metadata_summary(observed)
    return {
        "peak_count": len(paired),
        "tic": total,
        "base_peak_mz": base_mz,
        "base_peak_intensity": base_intensity,
        "mz_min": min(mz) if mz else None,
        "mz_max": max(mz) if mz else None,
        "intensity_mean": total / len(intensity) if intensity else None,
        "intensity_median": median_intensity,
        "precursor_region_intensity_fraction": precursor_region,
        "acquisition": acquisition,
    }


def render_observed_spectrum_html(payload_or_observed: dict[str, Any], output_path: str | Path | None = None) -> str:
    """Render an observed-spectrum-only QC view with metrics derivable from the MS2 spectrum itself."""
    observed = payload_or_observed.get("observed", payload_or_observed)
    metrics_payload = payload_or_observed.get("metrics") or {}
    candidate = metrics_payload.get("candidate") or {}
    predicted = payload_or_observed.get("predicted") or {}
    mz = _as_float_list(observed.get("mz", []))
    intensity = _as_float_list(observed.get("intensity", []))
    max_intensity = max((abs(i) for i in intensity), default=1.0) or 1.0
    rel_intensity = [abs(i) / max_intensity for i in intensity]
    stick_x, stick_y = _stick_values(mz, rel_intensity)
    metrics = observed_spectrum_metrics(observed)
    traces = [
        {
            "x": stick_x,
            "y": stick_y,
            "customdata": _repeat_stick_values(intensity),
            "type": "scatter",
            "mode": "lines",
            "name": "Observed spectrum",
            "line": {"color": "#068D9D", "width": 1.15},
            "hovertemplate": "m/z=%{x:.5f}<br>observed intensity=%{customdata:.3e}<extra>observed</extra>",
        }
    ]
    if mz:
        traces.append(
            _hover_scaffold_trace(
                name="observed",
                x_values=mz,
                y_values=rel_intensity,
                raw_intensities=intensity,
                intensity_label="observed intensity",
            )
        )
    spectrum_tiles = [
        _metric_tile("Peaks", _metric_value(metrics["peak_count"], decimals=0)),
        _metric_tile("TIC", _metric_value(metrics["tic"])),
        _metric_tile("Base peak m/z", _metric_value(metrics["base_peak_mz"])),
        _metric_tile("Base peak intensity", _metric_value(metrics["base_peak_intensity"])),
        _metric_tile("Median intensity", _metric_value(metrics["intensity_median"])),
    ]
    acquisition = _normalized_acquisition_summary(observed, predicted, (metrics_payload.get("acquisition_provenance") or {}).get("summary") or metrics["acquisition"])
    precursor_target_mz = _first_number(acquisition.get("isolation_target_mz"), candidate.get("precursor_mz"), predicted.get("precursor_mz"), observed.get("precursor_mz"))
    precursor_charge = _first_int(candidate.get("precursor_charge"), observed.get("precursor_charge"), predicted.get("precursor_charge"))
    metadata_panel = _prediction_metadata_panel(
        candidate=candidate,
        observed=observed,
        predicted=predicted,
        acquisition=acquisition,
        precursor_target_mz=precursor_target_mz,
        precursor_charge=precursor_charge,
        include_sequence=False,
        include_predictor=False,
    )
    metrics_panel = f'<div id="metrics-panel" class="metrics-panel">{metadata_panel}{_metric_group("Observed Spectrum", spectrum_tiles)}</div>'
    payload = {"observed": observed, "observed_metrics": metrics}
    return _basic_plot_html(
        traces=traces,
        layout_title="Observed spectrum",
        yaxis_title="Relative intensity",
        payload=payload,
        output_path=output_path,
        metrics_panel=metrics_panel,
    )


def _peptidoform_from_payload(payload_or_observed: dict[str, Any], peptidoform: str | None) -> str | None:
    if peptidoform:
        return peptidoform
    metrics = payload_or_observed.get("metrics") or {}
    candidate = metrics.get("candidate") or {}
    if candidate.get("peptidoform"):
        return str(candidate["peptidoform"])
    predicted = payload_or_observed.get("predicted") or {}
    if predicted.get("peptidoform"):
        return str(predicted["peptidoform"])
    if payload_or_observed.get("peptidoform"):
        return str(payload_or_observed["peptidoform"])
    return None


def render_annotated_observed_spectrum_html(
    payload_or_observed: dict[str, Any],
    output_path: str | Path | None = None,
    *,
    peptidoform: str | None = None,
    ion_types: str = "by",
    tolerance: float = 0.05,
    tolerance_unit: str = "Da",
    max_isotope: int = 2,
    max_ion_charge: int | None = None,
    neutral_losses: bool = False,
) -> str:
    """Render an observed-spectrum-only annotation view using the internal terminal-ion annotator."""
    observed = payload_or_observed.get("observed", payload_or_observed)
    peptide = _peptidoform_from_payload(payload_or_observed, peptidoform)
    mz = _as_float_list(observed.get("mz", []))
    intensity = _as_float_list(observed.get("intensity", []))
    max_intensity = max((abs(i) for i in intensity), default=1.0) or 1.0
    rel_intensity = [abs(i) / max_intensity for i in intensity]
    labels = [""] * len(mz)
    annotation_summary: dict[str, Any] = {"engine": "proteomics_mcp_terminal_ion_annotator", "annotated_peak_count": 0, "ion_type_counts": {}, "ion_types": ion_types}
    if peptide:
        try:
            observed_model = ObservedSpectrum(**observed)
            annotation_summary = annotate_observed_ions(
                observed_model,
                peptidoform=peptide,
                tolerance=tolerance,
                tolerance_unit=tolerance_unit,
                ion_types=ion_types,
                max_isotope=max_isotope,
                max_ion_charge=max_ion_charge,
                neutral_losses=neutral_losses,
            )
            labels = list(annotation_summary.get("annotations") or labels)
        except Exception as exc:
            annotation_summary["annotation_error"] = str(exc)
    ion_colors = {
        "a": "#009E73",
        "a+1": "#5AB4AC",
        "b": "#0072B2",
        "c": "#CC79A7",
        "d": "#8E6C8A",
        "d-alt": "#B07AA1",
        "x": "#56B4E9",
        "y": "#D55E00",
        "z": "#E69F00",
        "z-dot": "#F0B44C",
        "w": "#7A5195",
        "w-alt": "#A77BC2",
        "other": "#6B7280",
        "unannotated": "#A6B3C2",
    }
    groups: dict[str, dict[str, list[Any]]] = {}
    unannotated = {"x": [], "y": [], "raw": [], "text": []}
    match_by_index = {match.get("observed_index"): match for match in annotation_summary.get("matches", []) if isinstance(match, dict)}
    for index, (peak_mz, rel_y, raw_intensity, label) in enumerate(zip(mz, rel_intensity, intensity, labels)):
        if not label or str(label).startswith("?"):
            unannotated["x"].append(peak_mz)
            unannotated["y"].append(rel_y)
            unannotated["raw"].append(raw_intensity)
            unannotated["text"].append("")
            continue
        ion_type = str((match_by_index.get(index) or {}).get("ion_type") or str(label)[0].lower())
        group = groups.setdefault(ion_type, {"x": [], "y": [], "raw": [], "text": [], "sequence": []})
        group["x"].append(peak_mz)
        group["y"].append(rel_y)
        group["raw"].append(raw_intensity)
        group["text"].append(str(label))
        group["sequence"].append(_fragment_sequence_from_annotation(peptide, label))
    traces: list[dict[str, Any]] = []
    if unannotated["x"]:
        stick_x, stick_y = _stick_values(unannotated["x"], unannotated["y"])
        traces.append(
            {
                "x": stick_x,
                "y": stick_y,
                "customdata": _repeat_stick_values(unannotated["raw"]),
                "type": "scatter",
                "mode": "lines",
                "name": "Unannotated peaks",
                "line": {"color": ion_colors["unannotated"], "width": 1.0},
                "hovertemplate": "m/z=%{x:.5f}<br>observed intensity=%{customdata:.3e}<extra>unannotated</extra>",
            }
        )
    for ion_type, group in sorted(groups.items()):
        color = ion_colors.get(ion_type, ion_colors["other"])
        primary_ion_series = ion_type in {"b", "y"}
        stick_x, stick_y = _stick_values(group["x"], group["y"])
        traces.append(
            {
                "x": stick_x,
                "y": stick_y,
                "customdata": _repeat_stick_values(group["raw"]),
                "text": _repeat_stick_values(group["text"]),
                "meta": {
                    "baseText": _repeat_stick_values(group["text"]),
                    "fragmentSequenceText": _repeat_stick_values(group["sequence"]),
                    "ionType": ion_type,
                },
                "type": "scatter",
                "mode": "lines",
                "name": f"{ion_type}-ions",
                "line": {"color": color, "width": 1.55},
                "showlegend": primary_ion_series,
                "visible": primary_ion_series,
                "hovertemplate": "m/z=%{x:.5f}<br>observed intensity=%{customdata:.3e}<br>%{text}<extra>" + ion_type + "-ions</extra>",
            }
        )
        traces.append(
            _hover_scaffold_trace(
                name=f"{ion_type}-ions",
                x_values=group["x"],
                y_values=group["y"],
                raw_intensities=group["raw"],
                labels=group["text"],
                fragment_sequences=group["sequence"],
                intensity_label="observed intensity",
                ion_type=ion_type,
            )
        )
        traces[-1]["visible"] = primary_ion_series
        traces.append(
            {
                "x": group["x"],
                "y": group["y"],
                "text": group["text"],
                "meta": {
                    "baseText": group["text"],
                    "fragmentSequenceText": group["sequence"],
                    "ionType": ion_type,
                },
                "type": "scatter",
                "mode": "text",
                "name": f"{ion_type}-ion labels",
                "textfont": {"color": "#27313f", "size": 11},
                "textposition": "top center",
                "cliponaxis": False,
                "showlegend": False,
                "hoverinfo": "skip",
                "visible": primary_ion_series,
            }
        )
    counts = annotation_summary.get("ion_type_counts") or {}
    annotated_peak_count = annotation_summary.get("annotated_peak_count")
    total_theoretical_ions = annotation_summary.get("total_theoretical_ions")
    observed_metrics = observed_spectrum_metrics(observed)
    metrics_payload = payload_or_observed.get("metrics") or {}
    acquisition = _normalized_acquisition_summary(observed, payload_or_observed.get("predicted") or {}, (metrics_payload.get("acquisition_provenance") or {}).get("summary") or observed_metrics["acquisition"])
    annotation_tiles = [
        _metric_tile("Ion types", _metric_value(ion_types)),
        _metric_tile(
            "Matched peaks",
            (
                f"{_metric_value(annotated_peak_count, decimals=0)} / {_metric_value(total_theoretical_ions, decimals=0)}"
                if annotated_peak_count is not None and total_theoretical_ions is not None
                else _metric_value(annotated_peak_count, decimals=0)
            ),
            "of theoretical ions",
        ),
        _metric_tile("Annotated signal", _metric_value(annotation_summary.get("annotated_intensity_fraction"), percent=True), "of observed intensity"),
        _metric_tile("Observed ion series", _metric_value(_ion_count_summary(counts)), "annotated peaks"),
        _metric_tile("Tolerance", _metric_value(tolerance), tolerance_unit),
    ]
    if annotation_summary.get("annotation_error"):
        annotation_tiles.append(_metric_tile("Annotation error", _metric_value(annotation_summary["annotation_error"])))
    spectrum_tiles = [
        _metric_tile("Peaks", _metric_value(observed_metrics["peak_count"], decimals=0)),
        _metric_tile("TIC", _metric_value(observed_metrics["tic"])),
        _metric_tile("Base peak m/z", _metric_value(observed_metrics["base_peak_mz"])),
        _metric_tile("Base peak intensity", _metric_value(observed_metrics["base_peak_intensity"])),
    ]
    candidate = metrics_payload.get("candidate") or {}
    predicted = dict(payload_or_observed.get("predicted") or {})
    if peptide and not predicted.get("peptidoform"):
        predicted["peptidoform"] = peptide
    precursor_target_mz = _first_number(acquisition.get("isolation_target_mz"), candidate.get("precursor_mz"), predicted.get("precursor_mz"), observed.get("precursor_mz"))
    precursor_charge = _first_int(candidate.get("precursor_charge"), observed.get("precursor_charge"), predicted.get("precursor_charge"))
    metadata_panel = _prediction_metadata_panel(
        candidate=candidate,
        observed=observed,
        predicted=predicted,
        acquisition=acquisition,
        precursor_target_mz=precursor_target_mz,
        precursor_charge=precursor_charge,
        include_predictor=False,
    )
    candidate_tiles = [
        _metric_tile("Confidence", _metric_value(candidate.get("confidence"), percent=True), candidate.get("source_engine") or "candidate"),
        _metric_tile("Raw score", _metric_value(candidate.get("raw_score")), candidate.get("raw_score_name") or ""),
        _metric_tile("Next beam", _metric_value(candidate.get("next_beam")), "second-ranked sequence"),
        _metric_tile("Next beam margin", _metric_value(candidate.get("next_beam_margin")), "best minus second-ranked"),
    ]
    if candidate.get("token_confidence_mean") is not None:
        candidate_tiles.append(_metric_tile("Token confidence mean", _metric_value(candidate.get("token_confidence_mean"), percent=True), "exp(logp)"))
    if candidate.get("token_confidence_min") is not None:
        token = candidate.get("token_confidence_min_token")
        position = candidate.get("token_confidence_min_position")
        note = f"{token} at {position}" if token is not None and position is not None else "lowest token"
        candidate_tiles.append(_metric_tile("Token confidence min", _metric_value(candidate.get("token_confidence_min"), percent=True), note))
    metrics_panel = (
        f'<div id="metrics-panel" class="metrics-panel">'
        f"{metadata_panel}"
        f'{_metric_group("Observed Spectrum", spectrum_tiles)}'
        f'{_metric_group("Observed Annotation", annotation_tiles)}'
        f'{_metric_group("Candidate", candidate_tiles)}'
        "</div>"
    )
    payload = {
        "observed": observed,
        "peptidoform": peptide,
        "annotation_summary": annotation_summary,
        "annotation_settings": {
            "ion_types": ion_types,
            "tolerance": tolerance,
            "tolerance_unit": tolerance_unit,
            "max_isotope": max_isotope,
            "max_ion_charge": max_ion_charge,
            "neutral_losses": neutral_losses,
        },
    }
    return _basic_plot_html(
        traces=traces,
        layout_title="Annotated observed spectrum",
        yaxis_title="Relative intensity",
        payload=payload,
        output_path=output_path,
        metrics_panel=metrics_panel,
        controls_html=_ion_ladder_controls_markup(sorted(ion for ion in groups if ion not in {"b", "y"}), ion_colors) + _plot_controls_markup("fragment-sequence-toggle"),
        extra_script=_sequence_toggle_script("fragment-sequence-toggle") + _ion_ladder_toggle_script(),
    )


_STATIC_ION_COLORS = {
    "a": "#009E73",
    "a+1": "#5AB4AC",
    "b": "#0072B2",
    "c": "#CC79A7",
    "d": "#8E6C8A",
    "d-alt": "#B07AA1",
    "x": "#56B4E9",
    "y": "#D55E00",
    "z": "#E69F00",
    "z-dot": "#F0B44C",
    "w": "#7A5195",
    "w-alt": "#A77BC2",
    "immonium": "#5AB4AC",
    "reporter": "#C38D1A",
    "precursor": "#6B7280",
    "other": "#6B7280",
}


def _static_ion_type(label: str, match: dict[str, Any] | None) -> str:
    if match and match.get("ion_type"):
        return str(match["ion_type"])
    lowered = label.lower()
    if lowered.startswith("immonium"):
        return "immonium"
    if lowered.startswith("tmt") or lowered.startswith("itraq"):
        return "reporter"
    if lowered.startswith("p+"):
        return "precursor"
    parsed = re.match(r"^([a-z]+(?:[+-]\d+|-dot|-alt)?)", lowered)
    return parsed.group(1) if parsed else "other"


def _static_label(ax: Any, mz: float, intensity: float, label: str, color: str) -> None:
    ax.annotate(
        label,
        xy=(mz, intensity),
        xytext=(0, 4 if intensity >= 0 else -4),
        textcoords="offset points",
        ha="center",
        va="bottom" if intensity >= 0 else "top",
        fontsize=7,
        color=color,
        clip_on=True,
    )


def render_static_spectrum(
    payload_or_observed: dict[str, Any],
    output_path: str | Path,
    *,
    view: str = "mirror",
    peptidoform: str | None = None,
    ion_types: str = "by",
    tolerance: float = 0.05,
    tolerance_unit: str = "Da",
    max_isotope: int = 2,
    max_ion_charge: int | None = None,
    neutral_losses: bool = False,
    include_fragment_sequences: bool = False,
    label_top_n: int | None = None,
    title: str | None = None,
) -> Path:
    """Render a publication-ready static spectrum figure from a normalized payload.

    ``view`` is one of ``mirror``, ``annotated``, or ``plain``. The output format is
    determined by the suffix of ``output_path`` and can be SVG, PDF, PNG, or any
    Matplotlib-supported static format. When ``include_fragment_sequences`` is
    enabled for annotated spectra, sequence subsequences are added only to b/y
    fragment labels.
    """
    try:
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D
        import seaborn as sns
    except ImportError as exc:  # pragma: no cover - dependency contract is tested through installation
        raise RuntimeError("Static spectrum rendering requires matplotlib and seaborn. Install proteomics-mcp with its declared dependencies.") from exc

    if view not in {"mirror", "annotated", "plain"}:
        raise ValueError("view must be one of: mirror, annotated, plain")
    path = Path(output_path)
    if not path.suffix:
        raise ValueError("output_path must include a static image suffix such as .svg, .pdf, or .png")
    observed = payload_or_observed.get("observed", payload_or_observed)
    mz = _as_float_list(observed.get("mz", []))
    intensity = _as_float_list(observed.get("intensity", []))
    if not mz or len(mz) != len(intensity):
        raise ValueError("Observed spectrum must contain equally sized mz and intensity arrays.")
    maximum = max((abs(value) for value in intensity), default=1.0) or 1.0
    relative = [abs(value) / maximum for value in intensity]
    sns.set_theme(style="white", context="paper")
    fig, ax = plt.subplots(figsize=(11.5, 5.8), constrained_layout=True)
    fig.patch.set_facecolor("#FFFFFF")
    ax.set_facecolor("#FFFFFF")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.grid(False)
    ax.set_xlabel("m/z")

    if view == "plain":
        ax.vlines(mz, 0, relative, color="#068D9D", linewidth=0.8)
        ax.axhline(0, color="#4D6A7F", linewidth=0.7)
        ax.set_ylabel("Relative intensity")
        ax.set_ylim(0, 1.08)
        ax.set_title(title or "Observed spectrum")
    elif view == "annotated":
        peptide = _peptidoform_from_payload(payload_or_observed, peptidoform)
        if not peptide:
            raise ValueError("Annotated static rendering requires a peptidoform in the payload or peptidoform argument.")
        summary = annotate_observed_ions(
            ObservedSpectrum(**observed),
            peptidoform=peptide,
            tolerance=tolerance,
            tolerance_unit=tolerance_unit,
            ion_types=ion_types,
            max_isotope=max_isotope,
            max_ion_charge=max_ion_charge,
            neutral_losses=neutral_losses,
        )
        annotations = list(summary.get("annotations") or [""] * len(mz))
        match_by_index = {match.get("observed_index"): match for match in summary.get("matches", []) if isinstance(match, dict)}
        annotated: list[tuple[float, float, str, str, str]] = []
        for index, (peak_mz, rel_intensity, annotation) in enumerate(zip(mz, relative, annotations)):
            if not annotation or str(annotation).startswith("?"):
                ax.vlines(peak_mz, 0, rel_intensity, color="#A6B3C2", linewidth=0.65)
                continue
            label = str(annotation)
            ion_type = _static_ion_type(label, match_by_index.get(index))
            color = _STATIC_ION_COLORS.get(ion_type, _STATIC_ION_COLORS["other"])
            ax.vlines(peak_mz, 0, rel_intensity, color=color, linewidth=1.15)
            annotated.append((peak_mz, rel_intensity, label, color, ion_type))
        labels_to_render = sorted(annotated, key=lambda item: item[1], reverse=True)
        if label_top_n is not None:
            labels_to_render = labels_to_render[: max(0, label_top_n)]
        for peak_mz, rel_intensity, label, color, ion_type in labels_to_render:
            if include_fragment_sequences and ion_type in {"b", "y"}:
                sequence = _fragment_sequence_from_annotation(peptide, label)
                label = f"{label} | {sequence}" if sequence else label
            _static_label(ax, peak_mz, rel_intensity, label, color)
        ion_types_present = sorted({_static_ion_type(label, match_by_index.get(index)) for index, label in enumerate(annotations) if label and not str(label).startswith("?")})
        handles = [Line2D([0], [0], color="#A6B3C2", lw=1, label="Unannotated")] + [
            Line2D([0], [0], color=_STATIC_ION_COLORS.get(ion_type, _STATIC_ION_COLORS["other"]), lw=1.5, label=f"{ion_type}-ions")
            for ion_type in ion_types_present
        ]
        ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=min(5, len(handles)), frameon=False)
        ax.axhline(0, color="#4D6A7F", linewidth=0.7)
        ax.set_ylabel("Relative intensity")
        ax.set_ylim(0, 1.12)
        ax.set_title(title or "Annotated observed spectrum")
    else:
        predicted = payload_or_observed.get("predicted") or {}
        fragments = predicted.get("fragments") or []
        if not fragments:
            raise ValueError("Mirror static rendering requires predicted fragments in the payload.")
        metrics = payload_or_observed.get("metrics") or {}
        matched = metrics.get("matched_peaks") or []
        observed_matches = {match.get("observed_index"): str(match.get("ion_type") or "other") for match in matched if isinstance(match, dict)}
        predicted_matches = {match.get("predicted_index"): str(match.get("ion_type") or "other") for match in matched if isinstance(match, dict)}
        ax.vlines(mz, 0, relative, color="#9AA6B2", linewidth=0.65)
        labelled_observed: list[tuple[float, float, str, str]] = []
        for index, (peak_mz, rel_intensity) in enumerate(zip(mz, relative)):
            ion_type = observed_matches.get(index)
            if ion_type is None:
                continue
            color = _STATIC_ION_COLORS.get(ion_type, _STATIC_ION_COLORS["other"])
            ax.vlines(peak_mz, 0, rel_intensity, color=color, linewidth=1.15)
            match = next((item for item in matched if isinstance(item, dict) and item.get("observed_index") == index), {})
            labelled_observed.append((peak_mz, rel_intensity, str(match.get("annotation") or ""), color))
        predicted_mz = [float(fragment["mz"]) for fragment in fragments if fragment.get("mz") is not None]
        predicted_intensity = [float(fragment.get("intensity") or 0.0) for fragment in fragments if fragment.get("mz") is not None]
        predicted_maximum = max(predicted_intensity, default=1.0) or 1.0
        for index, (fragment_mz, fragment_intensity) in enumerate(zip(predicted_mz, predicted_intensity)):
            color = "#B8336A" if index in predicted_matches else "#D8A4BC"
            ax.vlines(fragment_mz, 0, -(fragment_intensity / predicted_maximum), color=color, linewidth=0.85)
        labels_to_render = sorted(labelled_observed, key=lambda item: item[1], reverse=True)
        if label_top_n is not None:
            labels_to_render = labels_to_render[: max(0, label_top_n)]
        for peak_mz, rel_intensity, label, color in labels_to_render:
            _static_label(ax, peak_mz, rel_intensity, label, color)
        ax.axhline(0, color="#4D6A7F", linewidth=0.7)
        ax.set_ylabel("Relative intensity\nobserved above, predicted below")
        ax.set_ylim(-1.12, 1.12)
        ax.set_title(title or "Observed vs predicted spectrum")
        handles = [
            Line2D([0], [0], color="#9AA6B2", lw=1, label="Observed unmatched"),
            Line2D([0], [0], color="#0072B2", lw=1.5, label="b-ions observed"),
            Line2D([0], [0], color="#D55E00", lw=1.5, label="y-ions observed"),
            Line2D([0], [0], color="#B8336A", lw=1.5, label="Predicted"),
        ]
        ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=4, frameon=False)

    ax.set_xlim(min(mz) - 15, max(mz) + 15)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=300, facecolor=fig.get_facecolor())
    plt.close(fig)
    return path
