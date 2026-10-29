"""Re-score the repository's example PSMs (examples/*_top5) with this checkout's package.

    uv run python scripts/regenerate_examples.py --cache-dir path/to/cache --source-of v0.2.0

The 25 PSMs (five per engine folder) are those of the manuscript's Zenodo supplement (record 23295902, md5-checked):
each is resolved from the supplement's identification table with the package's adapter and scored against its stored
observed spectrum with the settings it was published with (0.05 Da; annotation of all ion types, isotopes up to 2,
neutral losses) and MS2PIP HCD, under the package defaults (metric profile v2). The code is regenerate_figures.py's
(regenerate_examples), so the examples are scored exactly as the figures and the deposit's reference examples are, with
the same MS2PIP (FIGURE_MS2PIP, 4.1.2): the script runs itself again in a uv environment with this checkout, every other
dependency at its uv.lock version and that MS2PIP.

Writes, in examples/: each PSM's payload and its three HTML views (mirror, annotated, plain), under the file names they
already have, each folder's top5_summary.csv, and regeneration.json (package version and source digest, git commit,
MS2PIP and Python versions, the command and what changed). Exit codes: 0 done, 4 the supplement's checksum differs, 5 the
installed package is not this checkout's code, or with --source-of its source differs from that tag's.
"""

from __future__ import annotations

import argparse
import json
import platform
import shlex
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import regenerate_figures as figures  # noqa: E402

EXAMPLES = figures.REPO_ROOT / "examples"


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cache-dir", type=Path, required=True, help="where the Zenodo supplement is downloaded (or already is)")
    parser.add_argument("--zenodo-zip", type=Path, help="a local copy of the Zenodo supplement.zip (md5-checked); downloaded otherwise")
    parser.add_argument("--source-of", metavar="TAG", help="stop unless the package source equals this tag's (e.g. v0.2.0)")
    parser.add_argument("--python", default="3.12", help="Python for the MS2PIP environment (4.1.2 has no 3.14 wheels)")
    parser.add_argument("--in-figure-environment", action="store_true", help=argparse.SUPPRESS)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    args = parse_args(argv)
    if not args.in_figure_environment:
        with tempfile.TemporaryDirectory() as work:
            command = [*figures.pinned_environment([f"ms2pip=={figures.FIGURE_MS2PIP}"], args.python, Path(work), "examples"),
                       str(Path(__file__).resolve()), *argv, "--in-figure-environment"]
            figures.log(f"running with MS2PIP {figures.FIGURE_MS2PIP}: {shlex.join(command)}")
            return subprocess.run(command).returncode
    git = figures.git_state()  # before anything is written
    source_sha256 = figures.check_installed_package()
    source_of = figures.check_source_of(args.source_of, source_sha256) if args.source_of else None
    cache_dir = args.cache_dir.expanduser().resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)
    archive, supplement = figures.zenodo_supplement(args.zenodo_zip and args.zenodo_zip.expanduser().resolve(), cache_dir)
    figures.log(f"re-scoring examples/ with proteomics-mcp {figures.version('proteomics-mcp')} ({git['describe']}), "
                f"MS2PIP {figures.version('ms2pip')}, source {source_sha256[:16]}" + (f" = {source_of['tag']}" if source_of else ""))
    with tempfile.TemporaryDirectory() as unused:
        report = figures.regenerate_examples(supplement, Path(unused), {}, target_root=EXAMPLES, every_payload=True)
    record = {
        "generated": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "command": shlex.join(["uv", "run", "python", "scripts/regenerate_examples.py",
                               *(a.replace("\\", "/") for a in argv if a != "--in-figure-environment")]),
        "package": {"proteomics-mcp": figures.version("proteomics-mcp"), "source_sha256": source_sha256}, "source_of": source_of,
        "git": {key: git[key] for key in ("sha", "describe", "dirty")}, "ms2pip": figures.version("ms2pip"), "python": platform.python_version(),
        "inputs": {"zenodo_supplement": {"record": figures.ZENODO_RECORD, "md5": figures.file_digest(archive, "md5")}},
        "settings": {"metric_profile": "v2 (package default)", "predictor": figures.PREDICTOR,
                     "scoring_and_annotation": "as each payload was published (0.05 Da, ion types all, max isotope 2, neutral losses)",
                     "annotated_view": figures.EXAMPLE_HTML["annotated"]},
        "examples": report,
    }
    (EXAMPLES / "regeneration.json").write_text(json.dumps(record, indent=1, default=str).replace(str(Path.home()), "~") + "\n")
    figures.log(f"done: {sum(item['psms'] for item in report.values())} PSMs in {len(report)} folders")
    return 0


if __name__ == "__main__":
    sys.exit(main())
