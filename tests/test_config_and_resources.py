from __future__ import annotations

from proteomics_mcp.config import load_config
from proteomics_mcp.resources.docs import SCHEMAS


def test_load_config_from_toml_with_explicit_overrides(tmp_path):
    config_path = tmp_path / "proteomics.toml"
    config_path.write_text(
        'predictor_backend = "mock"\n'
        'ms2pip_model = "CID"\n'
        "max_fragment_charge = 2\n"
        "tolerance = 20.0\n"
        'tolerance_unit = "ppm"\n'
        "annotation_max_isotope = 2\n",
        encoding="utf-8",
    )

    config = load_config(config_path, {"tolerance": 0.02, "tolerance_unit": "Da"})

    assert config.predictor_backend == "mock"
    assert config.ms2pip_model == "CID"
    assert config.max_fragment_charge == 2
    assert config.tolerance == 0.02
    assert config.tolerance_unit == "Da"
    assert config.annotation_max_isotope == 2


def test_load_config_maps_predictor_model_aliases_and_koina_defaults(tmp_path):
    config_path = tmp_path / "proteomics.toml"
    config_path.write_text(
        'predictor_backend = "koina"\n',
        encoding="utf-8",
    )

    config = load_config(config_path)

    assert config.predictor_backend == "koina"
    assert config.predictor_model == "Prosit_2020_intensity_HCD"
    assert config.ms2pip_model == "Prosit_2020_intensity_HCD"


def test_explicit_config_predictor_model_is_not_replaced_by_cli_defaults(tmp_path):
    config_path = tmp_path / "proteomics.toml"
    config_path.write_text(
        'predictor_backend = "koina"\n'
        'predictor_model = "Prosit_2023_intensity_timsTOF"\n',
        encoding="utf-8",
    )

    config = load_config(config_path, {"predictor_backend": None, "predictor_model": None, "ms2pip_model": None})

    assert config.predictor_model == "Prosit_2023_intensity_timsTOF"


def test_mcp_resources_document_metrics_and_config():
    assert "precursor_mass_error_ppm" in SCHEMAS["method://metric_profile_legacy_v1"]
    assert "predictor_model" in SCHEMAS["config://proteomics_mcp_v1"]
    assert "koina_server_url" in SCHEMAS["config://proteomics_mcp_v1"]
    assert "max_fragment_charge" in SCHEMAS["config://proteomics_mcp_v1"]
    assert "annotation_max_isotope" in SCHEMAS["config://proteomics_mcp_v1"]
    assert "koina" in SCHEMAS["method://supported_predictors"]
    assert "CandidatePSM fields" in SCHEMAS["schema://candidate_psm_v1"]
    assert "Generic candidate mapping contract" in SCHEMAS["method://generic_mapping_v1"]
    assert "SDRF" in SCHEMAS["method://sdrf_support_v1"]
    assert "render_observed_spectrum" in SCHEMAS["method://visualization_modes_v1"]
    assert "summarize_dataset" in SCHEMAS["method://dataset_summary_v1"]
    assert "inspect_psm_table" in SCHEMAS["method://id_table_inspection_v1"]
