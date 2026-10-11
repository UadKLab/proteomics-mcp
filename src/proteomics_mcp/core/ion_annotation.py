from __future__ import annotations

import re


ION_ANNOTATION_RE = re.compile(r"^([A-Za-z][A-Za-z0-9-]*?)(\d+)(?:\+(\d+)|(\+*))?")
MS2PIP_ION_TYPE_RE = re.compile(r"^([A-Za-z]+)(\d+)?$")


def parse_ion_annotation(annotation: str) -> tuple[str | None, int | None, int]:
    """Return ion series, ordinal, and charge from common fragment labels."""
    match = ION_ANNOTATION_RE.match(str(annotation).strip())
    if not match:
        return None, None, 1
    if match.group(3):
        charge = int(match.group(3))
    else:
        charge = len(match.group(4)) or 1
    return match.group(1), int(match.group(2)), charge


def parse_ms2pip_ion_type(ion_type: str) -> tuple[str, int]:
    """Return MS2PIP ion-series and encoded charge from keys such as ``b2``."""
    match = MS2PIP_ION_TYPE_RE.fullmatch(str(ion_type))
    if not match:
        return str(ion_type), 1
    return match.group(1), int(match.group(2) or 1)
