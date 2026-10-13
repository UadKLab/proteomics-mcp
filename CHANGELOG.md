# Changelog

All notable changes to proteomics-mcp are listed here. The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow [Semantic Versioning](https://semver.org/). Detailed notes for each GitHub release are in [`docs/releases/`](https://github.com/UadKLab/proteomics-mcp/tree/main/docs/releases).

## [Unreleased]

Packaging, CI and documentation only. No changes to adapters, scoring, annotation, metrics or rendering.

### Fixed

- The Claude Desktop extension (`.mcpb`) failed to start when uv selected Python 3.14, because the locked `ms2rescore-rs` 0.4.3 has no Python 3.14 wheel ([#21](https://github.com/UadKLab/proteomics-mcp/issues/21)). The bundle now ships `.python-version` (3.12), and `uv.lock` pins `ms2pip` 4.2.0 / `ms2rescore-rs` 0.6.0 ([#20](https://github.com/UadKLab/proteomics-mcp/pull/20)). Workaround for v0.1.8: install Python 3.12, or set `UV_PYTHON=3.12`.
- MCP `initialize` reports the package version instead of an empty string ([#19](https://github.com/UadKLab/proteomics-mcp/pull/19)).
- On macOS without Homebrew `libomp`, MS2PIP prediction now fails with an actionable error (`brew install libomp`) instead of a raw `XGBoostError`. The message also reaches MCP clients ([#25](https://github.com/UadKLab/proteomics-mcp/pull/25)).
- The sdist ships the complete `tests/` directory, so its tests run against the installed wheel ([#18](https://github.com/UadKLab/proteomics-mcp/pull/18)).
- Removed the unused `JOBS_DIR` constant, which resolved outside the package when installed ([#18](https://github.com/UadKLab/proteomics-mcp/pull/18)).

### Changed

- **MS2PIP 4.2.0** in the locked environment used by the Claude Desktop extension (previously 4.1.2). Predictions are bit-identical for 2+ precursors. For 3+ precursors (`HCDch2`) they shift slightly, changing spectral angle, cosine and Pearson by less than 1e-7 on the example PSMs ([#20](https://github.com/UadKLab/proteomics-mcp/pull/20), [#23](https://github.com/UadKLab/proteomics-mcp/issues/23)).
- **Project metadata** ([#25](https://github.com/UadKLab/proteomics-mcp/pull/25)):
  - authors from `CITATION.cff`
  - classifiers for Python 3.11–3.14, `Development Status :: 4 - Beta` and `Intended Audience :: Science/Research`
  - `license-files`, a `setuptools>=77` build backend and a Changelog URL
- **Dependency:** dropped the unused `typer` ([#25](https://github.com/UadKLab/proteomics-mcp/pull/25)).
- **Docs** ([#25](https://github.com/UadKLab/proteomics-mcp/pull/25)):
  - README and install docs lead with `pip install proteomics-mcp` / `uvx proteomics-mcp serve`, and list platform requirements: glibc ≥ 2.34 for Python 3.13+ on Linux, macOS 14+ on Apple silicon, `brew install libomp` for MS2PIP on macOS.
  - The README, docs and `PRIVACY.md` state that MS2PIP downloads its models once, on first use.
- **Version:** `__version__` is read from the installed package metadata; `pyproject.toml` is the single source of truth ([#19](https://github.com/UadKLab/proteomics-mcp/pull/19)).

### Added

- **Release workflow:** tag-triggered with PyPI Trusted Publishing. It runs build → verify → TestPyPI → PyPI → GitHub release, plus a TestPyPI-only manual dry run ([#24](https://github.com/UadKLab/proteomics-mcp/pull/24)).
- **CI** tests the built wheel on Linux, macOS and Windows with Python 3.11–3.14, plus the Claude Desktop extension launch command ([#20](https://github.com/UadKLab/proteomics-mcp/pull/20)).
- **Scripts:**
  - `scripts/check_version.py`: version consistency across `pyproject.toml`, `manifest.json`, `CITATION.cff`, the tag and the release notes ([#19](https://github.com/UadKLab/proteomics-mcp/pull/19)).
  - `scripts/check_wheel.py`: wheel contents ([#18](https://github.com/UadKLab/proteomics-mcp/pull/18)).
  - `scripts/smoke_mcp.py`: MCP stdio smoke test ([#20](https://github.com/UadKLab/proteomics-mcp/pull/20)).
- **Docs:** `RELEASING.md` (release and contribution process) and this changelog ([#7](https://github.com/UadKLab/proteomics-mcp/issues/7)).

## [0.1.8] - 2026-10-11

This is the version cited in the manuscript, and the first release from the public repository snapshot. Scoring behaviour is unchanged from 0.1.7. The [release notes](https://github.com/UadKLab/proteomics-mcp/blob/main/docs/releases/v0.1.8.md) list the changes.

The GitHub release was published on 2026-10-11. Its wheel and sdist were uploaded unchanged to TestPyPI and PyPI on 2026-10-12; see [RELEASING.md](https://github.com/UadKLab/proteomics-mcp/blob/main/RELEASING.md#history-manual-upload-of-018) for the SHA-256 values.

## Earlier releases

Releases 0.1.1–0.1.7 predate the public repository and were not published to PyPI. Their notes:

- [0.1.7](https://github.com/UadKLab/proteomics-mcp/blob/main/docs/releases/v0.1.7.md)
- [0.1.6](https://github.com/UadKLab/proteomics-mcp/blob/main/docs/releases/v0.1.6.md)
- [0.1.5](https://github.com/UadKLab/proteomics-mcp/blob/main/docs/releases/v0.1.5.md)
- [0.1.4](https://github.com/UadKLab/proteomics-mcp/blob/main/docs/releases/v0.1.4.md)
- [0.1.3](https://github.com/UadKLab/proteomics-mcp/blob/main/docs/releases/v0.1.3.md)
- [0.1.2](https://github.com/UadKLab/proteomics-mcp/blob/main/docs/releases/v0.1.2.md)
- [0.1.1](https://github.com/UadKLab/proteomics-mcp/blob/main/docs/releases/v0.1.1.md)

[Unreleased]: https://github.com/UadKLab/proteomics-mcp/compare/v0.1.8...HEAD
[0.1.8]: https://github.com/UadKLab/proteomics-mcp/releases/tag/v0.1.8
