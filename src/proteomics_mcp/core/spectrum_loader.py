from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

import pandas as pd

from proteomics_mcp.models.schemas import ObservedSpectrum


def parse_scan_number(value: Any) -> int | None:
    if value is None:
        return None
    if isinstance(value, int):
        return value
    try:
        if not isinstance(value, str):
            return int(value)
    except (TypeError, ValueError):
        pass
    text = str(value)
    for pattern in [r"scan[=:](\d+)", r":scan:(\d+)", r":(\d+)$", r"\b(\d+)\b"]:
        match = re.search(pattern, text)
        if match:
            return int(match.group(1))
    return None


def to_float_list(value: Any) -> list[float]:
    if value is None:
        return []
    if isinstance(value, list):
        return [float(v) for v in value]
    if hasattr(value, "tolist"):
        return [float(v) for v in value.tolist()]
    if isinstance(value, str):
        stripped = value.strip()
        if not stripped:
            return []
        try:
            parsed = json.loads(stripped)
            if isinstance(parsed, list):
                return [float(v) for v in parsed]
        except Exception:
            pass
        return [float(v) for v in re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", stripped)]
    return [float(v) for v in value]


def _normalize_meta_value(value: Any) -> Any:
    if isinstance(value, bytes):
        try:
            return value.decode()
        except Exception:
            return str(value)
    if hasattr(value, "decode") and not isinstance(value, str):
        try:
            return value.decode()
        except Exception:
            pass
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    try:
        return float(value)
    except (TypeError, ValueError):
        return str(value)


def _is_missing(value: Any) -> bool:
    return value is None or pd.isna(value)


def _meta_dict(obj: Any) -> dict[str, Any]:
    keys: list[Any] = []
    try:
        obj.getKeys(keys)
    except Exception:
        return {}
    out: dict[str, Any] = {}
    for key in keys:
        try:
            name = key.decode() if isinstance(key, bytes) else str(key)
            out[name] = _normalize_meta_value(obj.getMetaValue(key))
        except Exception:
            continue
    return out


def _read_arrow_table(path: str | Path) -> pd.DataFrame:
    path = Path(path)
    try:
        if path.suffix.lower() in {".ipc", ".arrow", ".feather"}:
            return pd.read_feather(path)
    except ImportError as exc:
        raise RuntimeError("Reading IPC/Feather requires pyarrow to be installed.") from exc
    except Exception:
        try:
            import pyarrow.feather as feather  # type: ignore

            return feather.read_table(path).to_pandas()
        except ImportError as exc:
            raise RuntimeError("Reading IPC/Feather requires pyarrow to be installed.") from exc
    raise ValueError(f"Unsupported spectrum table format: {path.suffix}")


def load_spectrum_table(path: str | Path) -> pd.DataFrame:
    """Read an IPC/Feather spectrum table."""
    return _read_arrow_table(path)


def _write_spectrum_frame(frame: pd.DataFrame, output_path: str | Path) -> dict[str, Any]:
    try:
        frame.to_feather(output_path)
    except ImportError as exc:
        raise RuntimeError("Writing Feather requires pyarrow to be installed.") from exc
    return {"output_path": str(output_path), "spectra_written": int(len(frame))}


def _source_stem(path_value: Any) -> str | None:
    if path_value is None:
        return None
    text = str(path_value).strip()
    if not text:
        return None
    return Path(text).stem


def load_spectrum_from_table(
    path: str | Path,
    *,
    scan_number: int | None = None,
    spectrum_id: str | None = None,
    experiment_name: str | None = None,
) -> ObservedSpectrum:
    frame = load_spectrum_table(path)
    candidates = frame
    if experiment_name is not None and "source_file" in frame.columns:
        source_stems = frame["source_file"].map(_source_stem)
        filtered = frame[source_stems == str(experiment_name)]
        if not filtered.empty:
            candidates = filtered
    if scan_number is not None:
        if "scan_number" in candidates.columns:
            candidates = candidates[candidates["scan_number"].astype(int) == int(scan_number)]
        elif "scan" in candidates.columns:
            candidates = candidates[candidates["scan"].map(parse_scan_number) == int(scan_number)]
        else:
            raise ValueError("Spectrum table has no scan_number or scan column.")
    elif spectrum_id is not None:
        if "spectrum_id" in candidates.columns:
            candidates = candidates[candidates["spectrum_id"].astype(str) == str(spectrum_id)]
        elif "scan" in candidates.columns:
            candidates = candidates[candidates["scan"].astype(str) == str(spectrum_id)]
        else:
            raise ValueError("Spectrum table has no spectrum_id or scan column.")
    else:
        raise ValueError("Either scan_number or spectrum_id must be provided.")
    if candidates.empty:
        raise LookupError("No matching spectrum found in table.")
    row = candidates.iloc[0].to_dict()
    mz_col = "mz_array" if "mz_array" in row else "mz"
    int_col = "intensity_array" if "intensity_array" in row else "intensity"
    scan = row.get("scan_number") if row.get("scan_number") is not None else parse_scan_number(row.get("scan"))
    nested_metadata = row.get("metadata")
    metadata = dict(nested_metadata) if isinstance(nested_metadata, dict) else {}
    metadata.update({k: v for k, v in row.items() if k not in {mz_col, int_col, "metadata"}})
    return ObservedSpectrum(
        mz=to_float_list(row.get(mz_col)),
        intensity=to_float_list(row.get(int_col)),
        spectrum_id=str(row.get("spectrum_id") or row.get("scan") or scan),
        scan_number=int(scan) if scan is not None else None,
        precursor_mz=float(row["precursor_mz"]) if row.get("precursor_mz") is not None else None,
        precursor_charge=int(row["precursor_charge"]) if row.get("precursor_charge") is not None else None,
        rt_seconds=float(row["rt"]) if row.get("rt") is not None else None,
        source=str(path),
        metadata=metadata,
    )


# PROXI attributes carrying the precursor m/z, most specific first: selected ion m/z,
# isolation window target m/z, then a generic "precursor m/z" attribute.
_PROXI_PRECURSOR_MZ = (("MS:1000744", "selected ion m/z"), ("MS:1000827", "isolation window target m/z"), (None, "precursor m/z"))


def _proxi_precursor_mz(attributes: list[dict[str, Any]]) -> float | None:
    for accession, name in _PROXI_PRECURSOR_MZ:
        for attr in attributes:
            matches = (accession is not None and attr.get("accession") == accession) or name in str(attr.get("name", "")).lower()
            if matches and attr.get("value") is not None:
                try:
                    return float(attr["value"])
                except (TypeError, ValueError):
                    continue
    return None


def retrieve_spectrum_from_usi(usi: str, *, backend: str = "pride", timeout_seconds: float = 15.0) -> ObservedSpectrum:
    """Retrieve and normalize a USI spectrum through pyteomics PROXI."""
    try:
        from pyteomics import usi as pyteomics_usi  # type: ignore
    except ImportError as exc:
        raise RuntimeError("USI retrieval requires pyteomics to be installed.") from exc
    response = pyteomics_usi.proxi(usi, backend=backend, timeout=timeout_seconds)
    mz = response.get("m/z array")
    if mz is None:
        mz = response.get("mz")
    if mz is None:
        mz = response.get("mzs")
    if mz is None:
        mz = []
    intensity = response.get("intensity array")
    if intensity is None:
        intensity = response.get("intensity")
    if intensity is None:
        intensity = response.get("intensities")
    if intensity is None:
        intensity = []
    attributes = response.get("attributes") or []
    metadata = {
        "usi": usi,
        "backend": backend,
        "status": response.get("status"),
        "attributes": attributes,
    }
    precursor_mz = _proxi_precursor_mz(attributes if isinstance(attributes, list) else [])
    precursor_charge = None
    for attr in attributes if isinstance(attributes, list) else []:
        name = str(attr.get("name", "")).lower()
        value = attr.get("value")
        if "charge state" in name:
            precursor_charge = int(value)
    parsed_scan = parse_scan_number(usi)
    return ObservedSpectrum(
        mz=to_float_list(mz),
        intensity=to_float_list(intensity),
        spectrum_id=usi,
        scan_number=parsed_scan,
        precursor_mz=precursor_mz,
        precursor_charge=precursor_charge,
        source=f"usi:{backend}",
        metadata=metadata,
    )


def mzml_to_feather(mzml_path: str | Path, output_path: str | Path, *, ms_level: int = 2) -> dict[str, Any]:
    """Convert mzML spectra to a semi-wide Feather table via pyOpenMS."""
    try:
        import pyopenms  # type: ignore
    except ImportError as exc:
        raise RuntimeError("mzML conversion requires pyOpenMS and pyarrow to be installed.") from exc
    spectra = []
    exp = pyopenms.MSExperiment()
    pyopenms.MzMLFile().load(str(mzml_path), exp)
    source_name = Path(mzml_path).stem
    for spec in exp.getSpectra():
        if int(spec.getMSLevel()) != int(ms_level):
            continue
        mz, intensity = spec.get_peaks()
        spectrum_meta = _meta_dict(spec)
        precursors = spec.getPrecursors()
        precursor_mz = float(precursors[0].getMZ()) if precursors else None
        precursor_charge = int(precursors[0].getCharge()) if precursors and precursors[0].getCharge() else None
        precursor_meta = _meta_dict(precursors[0]) if precursors else {}
        native_id = spec.getNativeID()
        spectra.append(
            {
                "source_file": str(mzml_path),
                "source_name": source_name,
                "spectrum_id": native_id,
                "scan_number": parse_scan_number(native_id),
                "ms_level": int(spec.getMSLevel()),
                "rt": float(spec.getRT()),
                "precursor_mz": precursor_mz,
                "precursor_charge": precursor_charge,
                "spectrum_title": spectrum_meta.get("spectrum title"),
                "filter_string": spectrum_meta.get("filter string"),
                "scan_resolution": spectrum_meta.get("mass resolving power") or spectrum_meta.get("resolution"),
                "compensation_voltage": spectrum_meta.get("FAIMS compensation voltage") or spectrum_meta.get("compensation voltage"),
                "ion_injection_time": spectrum_meta.get("ion injection time"),
                "preset_scan_configuration": spectrum_meta.get("preset scan configuration"),
                "base_peak_mz": spectrum_meta.get("base peak m/z"),
                "base_peak_intensity": spectrum_meta.get("base peak intensity"),
                "total_ion_current": spectrum_meta.get("total ion current"),
                "lowest_observed_mz": spectrum_meta.get("lowest observed m/z"),
                "highest_observed_mz": spectrum_meta.get("highest observed m/z"),
                "precursor_spectrum_ref": precursor_meta.get("spectrum_ref"),
                "precursor_ms_level": precursor_meta.get("ms level"),
                "isolation_target_mz": precursor_meta.get("isolation window target m/z"),
                "isolation_lower_offset": precursor_meta.get("isolation window lower offset"),
                "isolation_upper_offset": precursor_meta.get("isolation window upper offset"),
                "collision_energy": precursor_meta.get("collision energy"),
                "metadata": {
                    "spectrum_header": {
                        "native_id": native_id,
                        "ms_level": int(spec.getMSLevel()),
                        "rt_seconds": float(spec.getRT()),
                        "source_file": str(mzml_path),
                        "source_name": source_name,
                    },
                    "spectrum_meta": spectrum_meta,
                    "precursor_header": {
                        "precursor_mz": precursor_mz,
                        "precursor_charge": precursor_charge,
                    },
                    "precursor_meta": precursor_meta,
                },
                "mz_array": [float(v) for v in mz],
                "intensity_array": [float(v) for v in intensity],
            }
        )
    frame = pd.DataFrame(spectra)
    result = _write_spectrum_frame(frame, output_path)
    result["ms_level"] = ms_level
    result["source_format"] = "mzml"
    return result


def mgf_to_feather(mgf_path: str | Path, output_path: str | Path) -> dict[str, Any]:
    """Convert MGF spectra to a semi-wide Feather table via pyteomics."""
    try:
        from pyteomics import mgf as pyteomics_mgf  # type: ignore
    except ImportError as exc:
        raise RuntimeError("MGF conversion requires pyteomics and pyarrow to be installed.") from exc

    mgf_path = Path(mgf_path)
    source_name = mgf_path.stem
    spectra: list[dict[str, Any]] = []
    with pyteomics_mgf.MGF(str(mgf_path), convert_arrays=1, read_charges=False, read_ions=False) as reader:
        for index, spectrum in enumerate(reader, start=1):
            params = spectrum.get("params") or {}
            title = str(params.get("title") or params.get("name") or f"{source_name}:index:{index}")
            scans = params.get("scans") or params.get("rawscans")
            scan_number = parse_scan_number(scans) or parse_scan_number(title) or index
            pepmass = params.get("pepmass")
            precursor_mz = None
            precursor_intensity = None
            if isinstance(pepmass, (list, tuple)) and pepmass:
                precursor_mz = float(pepmass[0])
                if len(pepmass) > 1 and pepmass[1] is not None:
                    precursor_intensity = float(pepmass[1])
            elif pepmass is not None:
                precursor_mz = float(pepmass)
            charge = params.get("charge")
            precursor_charge = None
            if isinstance(charge, (list, tuple)) and charge:
                precursor_charge = parse_scan_number(str(charge[0]).replace("+", ""))
            elif charge is not None:
                precursor_charge = parse_scan_number(str(charge).replace("+", ""))
            rt_seconds = params.get("rtinseconds")
            mobility = params.get("ion_mobility")
            mobility_value = None
            if isinstance(mobility, (list, tuple)) and len(mobility) >= 3:
                mobility_value = float(mobility[2])
            elif mobility is not None:
                try:
                    mobility_value = float(mobility)
                except (TypeError, ValueError):
                    mobility_value = None
            mz_array = [float(v) for v in spectrum.get("m/z array", [])]
            intensity_array = [float(v) for v in spectrum.get("intensity array", [])]
            spectra.append(
                {
                    "source_file": str(mgf_path),
                    "source_name": source_name,
                    "spectrum_id": title,
                    "scan_number": scan_number,
                    "ms_level": 2,
                    "rt": float(rt_seconds) if rt_seconds is not None else None,
                    "precursor_mz": precursor_mz,
                    "precursor_charge": precursor_charge,
                    "precursor_intensity": precursor_intensity,
                    "mobility": mobility_value,
                    "raw_scans": params.get("rawscans"),
                    "scan_label": params.get("scans"),
                    "metadata": {
                        "spectrum_header": {
                            "title": title,
                            "source_file": str(mgf_path),
                            "source_name": source_name,
                            "scan_number": scan_number,
                        },
                        "mgf_params": {str(k): _normalize_meta_value(v) for k, v in params.items()},
                    },
                    "mz_array": mz_array,
                    "intensity_array": intensity_array,
                }
            )
    result = _write_spectrum_frame(pd.DataFrame(spectra), output_path)
    result["source_format"] = "mgf"
    return result


def bruker_d_to_feather(d_path: str | Path, output_path: str | Path) -> dict[str, Any]:
    """Convert a Bruker .d directory to a semi-wide Feather table.

    Primary path: direct AlphaTims extraction from the .d directory.
    Fallback path: convert an embedded ProteinAnalysisResults.mgf when present.
    """
    d_path = Path(d_path)
    if not d_path.exists():
        raise FileNotFoundError(f"Bruker .d path not found: {d_path}")
    if not d_path.is_dir():
        raise ValueError(f"Bruker .d path must be a directory: {d_path}")

    try:
        from alphatims.bruker import TimsTOF  # type: ignore
    except ImportError:
        mgf_fallback = d_path / "ProteinAnalysisResults.mgf"
        if mgf_fallback.exists():
            result = mgf_to_feather(mgf_fallback, output_path)
            frame = pd.read_feather(output_path)
            frame["source_file"] = str(d_path)
            frame["source_name"] = d_path.stem
            frame["metadata"] = frame["metadata"].map(
                lambda meta: {
                    **(meta if isinstance(meta, dict) else {}),
                    "bruker_d_fallback": True,
                    "bruker_d_path": str(d_path),
                    "embedded_mgf_path": str(mgf_fallback),
                }
            )
            frame.to_feather(output_path)
            result["source_format"] = "bruker_d_embedded_mgf"
            result["fallback_used"] = "ProteinAnalysisResults.mgf"
            return result
        raise RuntimeError("Bruker .d conversion requires AlphaTims to be installed, or an embedded ProteinAnalysisResults.mgf fallback.")

    mgf_fallback = d_path / "ProteinAnalysisResults.mgf"
    try:
        raw = TimsTOF(str(d_path))
    except Exception:
        if mgf_fallback.exists():
            result = mgf_to_feather(mgf_fallback, output_path)
            frame = pd.read_feather(output_path)
            frame["source_file"] = str(d_path)
            frame["source_name"] = d_path.stem
            frame["metadata"] = frame["metadata"].map(
                lambda meta: {
                    **(meta if isinstance(meta, dict) else {}),
                    "bruker_d_fallback": True,
                    "bruker_d_path": str(d_path),
                    "embedded_mgf_path": str(mgf_fallback),
                }
            )
            frame.to_feather(output_path)
            result["source_format"] = "bruker_d_embedded_mgf"
            result["fallback_used"] = "ProteinAnalysisResults.mgf"
            return result
        raise
    source_name = d_path.stem
    spectra: list[dict[str, Any]] = []

    frames = getattr(raw, "frames", None)
    if frames is None:
        raise RuntimeError("AlphaTims reader did not expose a frames table.")
    fragment_frames = getattr(raw, "fragment_frames", None)
    precursors = getattr(raw, "precursors", None)
    frames_df = pd.DataFrame(frames).copy()
    fragment_df = pd.DataFrame(fragment_frames).copy() if fragment_frames is not None else pd.DataFrame()
    precursors_df = pd.DataFrame(precursors).copy() if precursors is not None else pd.DataFrame()

    if precursors_df.empty:
        raise RuntimeError("AlphaTims loaded the .d path but no precursors table was available.")

    if fragment_df.empty:
        if mgf_fallback.exists():
            result = mgf_to_feather(mgf_fallback, output_path)
            frame = pd.read_feather(output_path)
            frame["source_file"] = str(d_path)
            frame["source_name"] = source_name
            frame["metadata"] = frame["metadata"].map(
                lambda meta: {
                    **(meta if isinstance(meta, dict) else {}),
                    "bruker_d_fallback": True,
                    "bruker_d_path": str(d_path),
                    "embedded_mgf_path": str(mgf_fallback),
                }
            )
            frame.to_feather(output_path)
            result["source_format"] = "bruker_d_embedded_mgf"
            result["fallback_used"] = "ProteinAnalysisResults.mgf"
            return result
        raise RuntimeError("AlphaTims loaded the .d path but no fragment_frames table was available.")

    frame_id_column = next((col for col in ("Id", "frame", "Frame", "frame_id", "FrameID") if col in frames_df.columns), None)
    precursor_key = next((col for col in ("Id", "precursor", "Precursor", "precursor_id") if col in precursors_df.columns), None)
    fragment_precursor_column = next((col for col in ("Precursor", "precursor", "precursor_id", "PrecursorID") if col in fragment_df.columns), None)

    frames_by_id = frames_df.set_index(frame_id_column).to_dict(orient="index") if frame_id_column else {}
    fragment_by_precursor: dict[Any, dict[str, Any]] = {}
    if fragment_precursor_column:
        for _, fragment_row in fragment_df.iterrows():
            precursor_id = fragment_row.get(fragment_precursor_column)
            if precursor_id is None or precursor_id in fragment_by_precursor:
                continue
            fragment_by_precursor[precursor_id] = fragment_row.to_dict()

    spectrum_indptr, spectrum_tof_indices, spectrum_intensity_values = raw.index_precursors()

    precursor_rows = (
        precursors_df.sort_values(precursor_key).to_dict(orient="records")
        if precursor_key
        else precursors_df.to_dict(orient="records")
    )

    for precursor_meta in precursor_rows:
        precursor_id = precursor_meta.get(precursor_key) if precursor_key else None
        if precursor_id in (None, 0):
            continue
        mz_array, intensity_array = _extract_alphatims_precursor_peaks(
            raw,
            int(precursor_id),
            spectrum_indptr=spectrum_indptr,
            spectrum_tof_indices=spectrum_tof_indices,
            spectrum_intensity_values=spectrum_intensity_values,
        )
        if not mz_array or not intensity_array:
            continue
        frame_id = _first_present(precursor_meta, ["Parent", "parent", "Frame", "frame", "frame_id"])
        frame_meta = frames_by_id.get(frame_id, {}) if frame_id is not None else {}
        fragment_meta = fragment_by_precursor.get(precursor_id, {})
        precursor_mz = _first_present(precursor_meta, ["MonoisotopicMz", "mz", "Mz", "AverageMz", "LargestPeakMz"])
        precursor_charge = _first_present(precursor_meta, ["Charge", "charge"])
        rt_seconds = _first_present(frame_meta, ["Time", "time", "rt", "Rt"])
        mobility = _first_present(precursor_meta, ["ScanNumber", "scan_number", "Mobility", "mobility", "OneOverK0", "one_over_k0"])
        scan_number = _first_present(precursor_meta, ["ScanNumber", "scan_number"]) or frame_id or precursor_id
        spectrum_id = f"{source_name}:precursor:{int(precursor_id)}"
        spectra.append(
            {
                "source_file": str(d_path),
                "source_name": source_name,
                "spectrum_id": spectrum_id,
                "scan_number": int(float(scan_number)) if scan_number is not None else int(precursor_id),
                "ms_level": 2,
                "rt": float(rt_seconds) if rt_seconds is not None else None,
                "precursor_mz": float(precursor_mz) if precursor_mz is not None else None,
                "precursor_charge": int(precursor_charge) if precursor_charge is not None else None,
                "mobility": float(mobility) if mobility is not None else None,
                "frame_id": int(frame_id) if frame_id is not None else None,
                "precursor_id": int(precursor_id),
                "isolation_target_mz": _first_present(fragment_meta, ["IsolationMz", "isolation_mz"]),
                "isolation_width": _first_present(fragment_meta, ["IsolationWidth", "isolation_width"]),
                "collision_energy": _first_present(fragment_meta, ["CollisionEnergy", "collision_energy"]),
                "metadata": {
                    "spectrum_header": {
                        "source_file": str(d_path),
                        "source_name": source_name,
                        "frame_id": frame_id,
                        "precursor_id": precursor_id,
                        "scan_number": scan_number,
                    },
                    "alphatims_frame_meta": {str(k): _normalize_meta_value(v) for k, v in frame_meta.items()},
                    "alphatims_precursor_meta": {str(k): _normalize_meta_value(v) for k, v in precursor_meta.items()},
                    "alphatims_fragment_meta": {str(k): _normalize_meta_value(v) for k, v in fragment_meta.items()},
                },
                "mz_array": mz_array,
                "intensity_array": intensity_array,
            }
        )
    result = _write_spectrum_frame(pd.DataFrame(spectra), output_path)
    result["source_format"] = "bruker_d"
    return result


def _first_present(row: dict[str, Any], candidates: list[str]) -> Any:
    for column in candidates:
        if column in row and not _is_missing(row[column]):
            return row[column]
    return None


def _extract_alphatims_precursor_peaks(
    raw: Any,
    precursor_id: int,
    *,
    spectrum_indptr: Any | None = None,
    spectrum_tof_indices: Any | None = None,
    spectrum_intensity_values: Any | None = None,
) -> tuple[list[float], list[float]]:
    if spectrum_indptr is not None and spectrum_tof_indices is not None and spectrum_intensity_values is not None:
        start = int(spectrum_indptr[int(precursor_id)])
        stop = int(spectrum_indptr[int(precursor_id) + 1])
        if stop <= start:
            return [], []
        mz_values = raw.mz_values[spectrum_tof_indices[start:stop]]
        intensity_values = spectrum_intensity_values[start:stop]
        grouped = (
            pd.DataFrame({"mz": mz_values, "intensity": intensity_values})
            .groupby("mz", as_index=False)["intensity"]
            .sum()
            .sort_values("mz", kind="stable")
        )
        return (
            [float(v) for v in grouped["mz"].tolist()],
            [float(v) for v in grouped["intensity"].tolist()],
        )
    try:
        result = raw[:, :, int(precursor_id), :, "df"]
    except Exception as exc:
        raise RuntimeError(
            "AlphaTims precursor extraction failed in this environment; "
            "use the embedded MGF fallback or update the adapter."
        ) from exc
    if not isinstance(result, pd.DataFrame) or result.empty:
        return [], []
    mz_column = "mz_values" if "mz_values" in result.columns else "mz"
    intensity_column = "intensity_values" if "intensity_values" in result.columns else "intensity"
    if mz_column not in result.columns or intensity_column not in result.columns:
        return [], []
    grouped = (
        result[[mz_column, intensity_column]]
        .groupby(mz_column, as_index=False)[intensity_column]
        .sum()
        .sort_values(mz_column, kind="stable")
    )
    return (
        [float(v) for v in grouped[mz_column].tolist()],
        [float(v) for v in grouped[intensity_column].tolist()],
    )
