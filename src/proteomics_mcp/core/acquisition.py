from __future__ import annotations

import math
import re
from typing import Any

from proteomics_mcp.models.schemas import ObservedSpectrum, PredictedSpectrum


def _get(mapping: Any, key: str, default: Any = None) -> Any:
    if isinstance(mapping, dict):
        return mapping.get(key, default)
    return getattr(mapping, key, default)


def _metadata(spectrum: ObservedSpectrum | dict[str, Any]) -> dict[str, Any]:
    metadata = _get(spectrum, "metadata", {}) or {}
    return metadata if isinstance(metadata, dict) else {}


def _attributes(spectrum: ObservedSpectrum | dict[str, Any]) -> list[dict[str, Any]]:
    metadata = _metadata(spectrum)
    attrs = metadata.get("attributes") or metadata.get("params") or []
    return [attr for attr in attrs if isinstance(attr, dict)] if isinstance(attrs, list) else []


def _finite_float(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _first_metadata_number(spectrum: ObservedSpectrum | dict[str, Any], *keys: str) -> float | None:
    metadata = _metadata(spectrum)
    for key in keys:
        value = _finite_float(metadata.get(key))
        if value is not None:
            return value
    return None


def _format_value_with_unit(value: Any, unit_name: Any = None, unit_accession: Any = None) -> str | None:
    if value is None or value == "":
        return None
    suffix = str(unit_name or unit_accession or "").strip()
    return f"{value} {suffix}".strip()


def precursor_target_mz(spectrum: ObservedSpectrum | dict[str, Any]) -> float | None:
    direct = _finite_float(_get(spectrum, "precursor_mz"))
    if direct is not None:
        return direct
    metadata_value = _first_metadata_number(spectrum, "precursor_mz", "isolation_target_mz", "isolation window target m/z")
    if metadata_value is not None:
        return metadata_value
    target_accessions = {"MS:1000827", "MS:1000744"}
    target_names = ("isolation window target m/z", "selected ion m/z", "precursor m/z")
    for attr in _attributes(spectrum):
        accession = str(attr.get("accession") or "")
        name = str(attr.get("name") or "").lower()
        if accession in target_accessions or any(target in name for target in target_names):
            target = _finite_float(attr.get("value"))
            if target is not None:
                return target
    return None


def retention_time_seconds(spectrum: ObservedSpectrum | dict[str, Any]) -> float | None:
    direct = _finite_float(_get(spectrum, "rt_seconds"))
    if direct is not None:
        return direct
    metadata_value = _first_metadata_number(spectrum, "rt_seconds", "retention_time_seconds", "retention_time", "rt")
    if metadata_value is not None:
        return metadata_value
    for attr in _attributes(spectrum):
        name = str(attr.get("name") or "").lower()
        accession = str(attr.get("accession") or "")
        if accession != "MS:1000016" and "scan start time" not in name and "retention time" not in name:
            continue
        value = _finite_float(attr.get("value"))
        if value is None:
            continue
        unit = str(attr.get("unitName") or attr.get("unit_name") or attr.get("unitAccession") or "").lower()
        if "minute" in unit or unit == "uo:0000031":
            return value * 60.0
        return value
    return None


def collision_energy(spectrum: ObservedSpectrum | dict[str, Any]) -> str | None:
    metadata = _metadata(spectrum)
    for key in ("collision_energy", "normalized_collision_energy", "hcd_collision_energy"):
        formatted = _format_value_with_unit(metadata.get(key), metadata.get(f"{key}_unit"))
        if formatted:
            return formatted
    for attr in _attributes(spectrum):
        name = str(attr.get("name") or "").lower()
        accession = str(attr.get("accession") or "")
        if accession not in {"MS:1000045", "MS:1000138"} and "collision energy" not in name:
            continue
        formatted = _format_value_with_unit(attr.get("value"), attr.get("unitName") or attr.get("unit_name"), attr.get("unitAccession"))
        if formatted:
            return formatted
    return None


def faims_compensation_voltage(spectrum: ObservedSpectrum | dict[str, Any]) -> str | None:
    metadata = _metadata(spectrum)
    for key in ("compensation_voltage", "faims_cv", "cv", "comp_voltage"):
        formatted = _format_value_with_unit(metadata.get(key), metadata.get(f"{key}_unit"))
        if formatted:
            return formatted
        value = metadata.get(key)
        if value not in (None, ""):
            return str(value)
    for attr in _attributes(spectrum):
        name = str(attr.get("name") or "").lower()
        if "compensation voltage" not in name:
            continue
        formatted = _format_value_with_unit(attr.get("value"), attr.get("unitName") or attr.get("unit_name"), attr.get("unitAccession"))
        if formatted:
            return formatted
    return None


def injection_time(spectrum: ObservedSpectrum | dict[str, Any]) -> str | None:
    metadata = _metadata(spectrum)
    for key in ("ion_inject_time", "ion_injection_time", "injection_time", "inject_time"):
        formatted = _format_value_with_unit(metadata.get(key), metadata.get(f"{key}_unit"))
        if formatted:
            return formatted
        value = metadata.get(key)
        if value not in (None, ""):
            return str(value)
    for attr in _attributes(spectrum):
        name = str(attr.get("name") or "").lower()
        if "ion injection time" not in name and "inject time" not in name:
            continue
        formatted = _format_value_with_unit(attr.get("value"), attr.get("unitName") or attr.get("unit_name"), attr.get("unitAccession"))
        if formatted:
            return formatted
    return None


def fragmentation_method(spectrum: ObservedSpectrum | dict[str, Any], predicted: PredictedSpectrum | dict[str, Any] | None = None) -> str | None:
    metadata = _metadata(spectrum)
    for key in ("fragmentation_method", "activation_method", "dissociation_method"):
        value = metadata.get(key)
        if value:
            return str(value)
    method_accessions = {
        "MS:1000133": "CID",
        "MS:1000422": "HCD",
        "MS:1002481": "HCD",
        "MS:1000598": "ETD",
        "MS:1002631": "EThcD",
        "MS:1002679": "UVPD",
    }
    method_keywords = {
        "higher energy": "HCD",
        "beam-type collision": "HCD",
        "collision-induced": "CID",
        "electron transfer": "ETD",
        "ethcd": "EThcD",
        "ultraviolet photodissociation": "UVPD",
    }
    for attr in _attributes(spectrum):
        accession = str(attr.get("accession") or "")
        name = str(attr.get("name") or "")
        lower_name = name.lower()
        if accession in method_accessions:
            return method_accessions[accession]
        for keyword, method in method_keywords.items():
            if keyword in lower_name:
                return method
    if predicted is not None:
        model = _get(predicted, "model") or _get(predicted, "predictor_name")
        if model:
            return str(model)
    return None


def instrument_model(spectrum: ObservedSpectrum | dict[str, Any]) -> str | None:
    metadata = _metadata(spectrum)
    for key in ("instrument_model", "instrument_name", "instrument"):
        value = metadata.get(key)
        if value:
            return str(value)
    instrument_accessions = {
        "MS:1000031",  # instrument model
        "MS:1000463",  # instrument
        "MS:1000483",  # Thermo Fisher Scientific instrument model
        "MS:1000494",  # LCQ Deca
        "MS:1000495",  # LCQ Deca XP
        "MS:1000554",  # LTQ
        "MS:1000638",  # LTQ Orbitrap
        "MS:1001742",  # LTQ Orbitrap Velos
        "MS:1002416",  # Q Exactive
        "MS:1002732",  # Orbitrap Fusion Lumos
        "MS:1003028",  # Orbitrap Eclipse
    }
    instrument_keywords = (
        "instrument model",
        "mass spectrometer",
        "orbitrap fusion",
        "orbitrap eclipse",
        "orbitrap exploris",
        "q exactive",
        "ltq orbitrap",
        "timsTOF".lower(),
    )
    for attr in _attributes(spectrum):
        accession = str(attr.get("accession") or "")
        name = str(attr.get("name") or "")
        lower_name = name.lower()
        value = attr.get("value")
        if accession in instrument_accessions:
            return str(value or name) if (value or name) else None
        if any(keyword in lower_name for keyword in instrument_keywords):
            return str(value or name) if (value or name) else None
    return None


def detector_type(spectrum: ObservedSpectrum | dict[str, Any]) -> str | None:
    metadata = _metadata(spectrum)
    for key in ("detector_type", "detector", "mass_analyzer", "instrument_analyzer"):
        value = metadata.get(key)
        if value:
            return str(value)
    detector_accessions = {
        "MS:1000082": "electron multiplier",
        "MS:1000253": "faraday cup",
        "MS:1000111": "microchannel plate detector",
        "MS:1000624": "inductive detector",
        "MS:1000484": "orbitrap",
        "MS:1000291": "linear ion trap",
        "MS:1000084": "time-of-flight",
        "MS:1000079": "fourier transform ion cyclotron resonance",
    }
    instrument_accessions = {"MS:1000031", "MS:1000463", "MS:1000483", "MS:1000638", "MS:1001742", "MS:1002416", "MS:1002732", "MS:1003028"}
    detector_keywords = ("detector", "orbitrap", "ion trap", "tof", "time-of-flight", "quadrupole", "ft-icr", "mass analyzer")
    instrument_keywords = ("fusion", "lumos", "eclipse", "exploris", "q exactive", "ltq orbitrap", "mass spectrometer", "instrument model")
    for attr in _attributes(spectrum):
        accession = str(attr.get("accession") or "")
        name = str(attr.get("name") or "")
        lower_name = name.lower()
        if accession in instrument_accessions or any(keyword in lower_name for keyword in instrument_keywords):
            continue
        if accession in detector_accessions:
            return detector_accessions[accession]
        if any(keyword in lower_name for keyword in detector_keywords):
            return name or None
    return None


def scan_resolution(spectrum: ObservedSpectrum | dict[str, Any]) -> str | None:
    metadata = _metadata(spectrum)
    for key in ("scan_resolution", "resolution", "mass_resolution", "resolving_power"):
        formatted = _format_value_with_unit(metadata.get(key), metadata.get(f"{key}_unit"))
        if formatted:
            return formatted
        value = metadata.get(key)
        if value not in (None, ""):
            return str(value)
    resolution_accessions = {
        "MS:1000011",  # mass resolution
        "MS:1000024",  # resolving power
    }
    resolution_keywords = ("resolution", "resolving power")
    for attr in _attributes(spectrum):
        accession = str(attr.get("accession") or "")
        name = str(attr.get("name") or "")
        lower_name = name.lower()
        if accession not in resolution_accessions and not any(keyword in lower_name for keyword in resolution_keywords):
            continue
        formatted = _format_value_with_unit(attr.get("value"), attr.get("unitName") or attr.get("unit_name"), attr.get("unitAccession"))
        if formatted:
            return formatted
        if name:
            return name
    return None


def scan_header(spectrum: ObservedSpectrum | dict[str, Any]) -> str | None:
    metadata = _metadata(spectrum)
    for key in ("filter_string", "scan_header", "spectrum_title", "title"):
        value = metadata.get(key)
        if value:
            return str(value)
    spectrum_meta = metadata.get("spectrum_meta")
    if isinstance(spectrum_meta, dict):
        for key in ("filter string", "spectrum title"):
            value = spectrum_meta.get(key)
            if value:
                return str(value)
    for attr in _attributes(spectrum):
        accession = str(attr.get("accession") or "")
        name = str(attr.get("name") or "").lower()
        value = attr.get("value")
        if accession == "MS:1000512" or "filter string" in name or "spectrum title" in name:
            if value:
                return str(value)
    return None


# The activation and its normalized collision energy in a Thermo scan filter: "... Full ms2 498.7812@hcd28.00 [...]"
_FILTER_NCE = re.compile(r"@(hcd|cid|etd|ethcd|etcid)(\d+(?:\.\d+)?)", re.IGNORECASE)
_THERMO_FILTER_PREFIXES = ("FTMS", "ITMS", "ASTMS")
_THERMO_INSTRUMENTS = ("thermo", "orbitrap", "exactive", "exploris", "eclipse", "lumos", "fusion", "ltq", "velos", "elite", "astral", "tribrid")


def normalized_collision_energy(spectrum: ObservedSpectrum | dict[str, Any]) -> tuple[float | None, str]:
    """The spectrum's normalized collision energy (NCE), when it is known to be one, and where it came from (#128).

    Prosit and the other collision-energy-dependent Koina models take an NCE. Thermo instruments state it in the scan
    filter ("@hcd28.00"), and Thermo mzML files carry the same number as the collision energy. Bruker timsTOF and SCIEX
    report an absolute energy in eV, which is not an NCE, so it is never returned.
    """
    header = scan_header(spectrum)
    if header:
        match = _FILTER_NCE.search(header)
        if match:
            return float(match.group(2)), f"scan filter {match.group(0)}"
    energy = _first_metadata_number(spectrum, "collision_energy", "normalized_collision_energy", "hcd_collision_energy")
    if energy is None:
        for attr in _attributes(spectrum):
            if str(attr.get("accession") or "") in {"MS:1000045", "MS:1000138"} or "collision energy" in str(attr.get("name") or "").lower():
                energy = _finite_float(attr.get("value"))
                if energy is not None:
                    break
    if energy is None:
        return None, "the spectrum states no collision energy"
    model = (instrument_model(spectrum) or "").lower()
    if (header and header.strip().upper().startswith(_THERMO_FILTER_PREFIXES)) or any(name in model for name in _THERMO_INSTRUMENTS):
        return energy, "Thermo collision energy (an NCE)"
    return None, f"collision energy {energy:g} is not known to be an NCE (not a Thermo spectrum); not used"


def acquisition_metadata_summary(
    spectrum: ObservedSpectrum | dict[str, Any],
    predicted: PredictedSpectrum | dict[str, Any] | None = None,
) -> dict[str, Any]:
    rt = retention_time_seconds(spectrum)
    return {
        "fragmentation_method": fragmentation_method(spectrum, predicted),
        "isolation_target_mz": precursor_target_mz(spectrum),
        "retention_time_seconds": rt,
        "retention_time_minutes": rt / 60.0 if rt is not None else None,
        "collision_energy": collision_energy(spectrum),
        "faims_compensation_voltage": faims_compensation_voltage(spectrum),
        "injection_time": injection_time(spectrum),
        "scan_header": scan_header(spectrum),
        "scan_resolution": scan_resolution(spectrum),
        "instrument_model": instrument_model(spectrum),
        "detector_type": detector_type(spectrum),
    }
