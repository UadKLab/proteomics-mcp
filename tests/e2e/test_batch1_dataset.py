"""Batch 1, item 2 of #30: distribution summaries, dashboards and reports on the Zenodo supplement."""

from __future__ import annotations

import json
import re
import shutil

import pytest

pytestmark = pytest.mark.e2e

FRAGPIPE = "processed/arabidopsis/fragpipe_psm.tsv"
SCORED = "processed/arabidopsis/supplementary_scored_psms.csv"
SDRF = "processed/arabidopsis/sdrf.tsv"
FIGURE4 = "processed/arabidopsis/figure4_dashboard_plot_data.json"
INSTANOVO_LABELLED = "processed/hela/instanovo_labelled_kpreds.csv"


def _ok(tool, name, arguments):
    is_error, result = tool(name, arguments)
    assert not is_error, f"{name}: {result}"
    return result


def test_distribution_tools(zenodo, tool, report):
    base = {"path": str(zenodo / FRAGPIPE), "adapter": "fragpipe"}
    calls = {
        "list_distribution_fields": {},
        "summarize_numeric_field": {"field": "precursor_mz"},
        "summarize_categorical_field": {"field": "precursor_charge"},
        "summarize_missing_values": {},
        "summarize_modifications": {},
        "summarize_file_coverage": {},
        "summarize_overlap": {},
        "summarize_cv_by_category": {"value_field": "source_columns.Intensity", "category_field": "precursor_charge"},
        "summarize_quant_quality": {"quant_field": "source_columns.Intensity"},
        "bin_numeric_field": {"field": "precursor_mz", "bin_count": 20},
        "bin_2d_numeric_fields": {"x_field": "precursor_mz", "y_field": "source_columns.Retention", "x_bin_count": 10, "y_bin_count": 10},
        "summarize_dashboard": {},
    }
    results = {}
    for name, extra in calls.items():
        is_error, result = tool(name, {**base, **extra})
        results[name] = "error: " + str(result)[:200] if is_error else "ok"
    report("dataset", "distribution tools on fragpipe_psm.tsv", "pass" if all(v == "ok" for v in results.values()) else "fail", results)
    assert all(v == "ok" for v in results.values()), results


def test_prediction_accuracy(zenodo, tool, report, tmp_path):
    summary = _ok(tool, "summarize_prediction_accuracy", {"path": str(zenodo / INSTANOVO_LABELLED), "adapter": "instanovo",
                                                          "target_field": "source_columns.targets"})
    target = tmp_path / "accuracy.html"
    _ok(tool, "render_prediction_accuracy_html", {"summary": summary, "output_path": str(target)})
    assert target.is_file() and target.stat().st_size > 1000
    report("dataset", "summarize_prediction_accuracy / render_prediction_accuracy_html (InstaNovo labelled)", "pass",
           {"summary_keys": sorted(summary)[:12], "html_bytes": target.stat().st_size})


def test_reports_and_summary_html(zenodo, tool, report, tmp_path):
    scored, sdrf = str(zenodo / SCORED), str(zenodo / SDRF)
    report_html = tmp_path / "dataset_report.html"
    _ok(tool, "render_dataset_report", {"scored_path": scored, "sdrf_path": sdrf, "output_path": str(report_html)})
    summary = _ok(tool, "summarize_dashboard", {"path": str(zenodo / FRAGPIPE), "adapter": "fragpipe"})
    dashboard_html = tmp_path / "dashboard.html"
    _ok(tool, "render_dashboard_html", {"summary": summary, "output_path": str(dashboard_html)})
    dataset_summary = _ok(tool, "summarize_dataset", {"scored_path": scored, "sdrf_path": sdrf})
    summary_html = tmp_path / "summary.html"
    _ok(tool, "render_summary_html", {"summary": dataset_summary, "output_path": str(summary_html)})
    sizes = {p.name: p.stat().st_size for p in (report_html, dashboard_html, summary_html)}
    assert all(size > 1000 for size in sizes.values()), sizes
    report("dataset", "render_dataset_report / render_dashboard_html / render_summary_html", "pass", {"html_bytes": sizes})


def _plot_ids(html_path) -> set[str]:
    return set(re.findall(r"Plotly\.newPlot\(\s*[\"']([^\"']+)[\"']", html_path.read_text(encoding="utf-8")))


@pytest.mark.parametrize("table_name", [
    pytest.param("fragpipe_psm.tsv", id="zenodo-name",
                 marks=pytest.mark.xfail(strict=True, reason="#39: candidate_id depends on the file name; the join matches nothing")),
    pytest.param("psm.tsv", id="original-name"),
])
def test_figure4_score_aware_dashboard(zenodo, tool, report, tmp_path, table_name):
    """Rebuild the Figure 4 FragPipe dashboard and compare its plots with the published plot data.

    The scored table joins on candidate_id, which the adapters derive from the input file name
    ("<stem>:<row>"). The Zenodo table is named fragpipe_psm.tsv, the scored table refers to "psm:<row>".
    """
    table = tmp_path / table_name
    shutil.copy(zenodo / FRAGPIPE, table)
    prefix = tmp_path / "fig4"
    result = _ok(tool, "build_id_dashboard", {"path": str(table), "adapter": "fragpipe", "scored_path": str(zenodo / SCORED),
                                              "output_prefix": str(prefix)})
    html = tmp_path / "fig4_dashboard.html"
    reference = set(json.loads((zenodo / FIGURE4).read_text()))
    rebuilt = _plot_ids(html)
    import pandas as pd

    enriched = pd.read_csv(result["output_files"]["enriched_csv"], usecols=["spectral_angle"], low_memory=False)
    joined = int(enriched["spectral_angle"].notna().sum())
    report("dataset", f"build_id_dashboard Figure 4 ({table_name})", "pass" if rebuilt == reference else "differs", {
        "psms": result["dataset_metadata"]["total_psm_count"], "rows_with_scores_after_join": joined,
        "plots_rebuilt": len(rebuilt), "plots_reference": len(reference),
        "missing_vs_reference": sorted(reference - rebuilt), "extra_vs_reference": sorted(rebuilt - reference)})
    assert rebuilt <= reference, sorted(rebuilt - reference)
    assert joined > 0, "scored table matched no candidate (join on candidate_id)"
