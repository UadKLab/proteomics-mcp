"""Keep untrusted strings from becoming markup or script in the generated HTML reports (#78).

Run names, spectrum ids, peptides, protein descriptions and file names come from input tables and spectra. They are
written into HTML text with ``text`` and into ``<script>`` blocks as JSON with ``to_script_json`` / ``script_json``.
"""

from __future__ import annotations

import html
import json
from typing import Any

# json.dumps leaves "</script>" and "<!--" intact, which would end or alter the surrounding <script> element.
# Escaping these characters keeps the JSON valid: they can only occur inside JSON strings, where \uXXXX means the same.
_SCRIPT_ESCAPES = str.maketrans({"<": "\\u003c", ">": "\\u003e", "&": "\\u0026", " ": "\\u2028", " ": "\\u2029"})


def text(value: Any) -> str:
    """A value as HTML text, safe in element content and in quoted attribute values."""
    return html.escape("" if value is None else str(value), quote=True)


def script_json(json_text: str) -> str:
    """JSON text made safe to embed in a <script> element."""
    return json_text.translate(_SCRIPT_ESCAPES)


def to_script_json(value: Any, **dumps_kwargs: Any) -> str:
    """json.dumps for embedding in a <script> element."""
    return script_json(json.dumps(value, **dumps_kwargs))
