# Publishing

## Releases, PyPI and the Claude Desktop extension

Releases are built, verified and published by the `release` workflow when a version tag is pushed. The workflow publishes to TestPyPI, then to PyPI using Trusted Publishing, and finally creates the GitHub release with the wheel, sdist and `.mcpb` extension attached. The step-by-step process is in [RELEASING.md](https://github.com/UadKLab/proteomics-mcp/blob/main/RELEASING.md): the one-time setup, the release checklist, TestPyPI verification and yanking.

Install the `.mcpb` attachment in Claude Desktop to test the packaged extension. A GitHub release does not automatically submit it to Anthropic's directory. Before directory submission, test installation on clean macOS and Windows systems, and review the third-party privacy notices listed in `PRIVACY.md`.

## Documentation site

MkDocs builds the documentation from `docs/` with `mkdocs.yml`. Run `uvx --from mkdocs==1.6.1 mkdocs build --strict` locally. The `docs` workflow checks the build on pushes and pull requests, and deploys pushes to `main` to GitHub Pages (publishing source: **GitHub Actions**). The site is live at `https://uadklab.github.io/proteomics-mcp/`.

## Paper and public-repository checklist

- Keep public history limited to the sanitized snapshot. Do not reintroduce old local paths or generated files from private archives.
- Confirm that each example dataset can be shared, document its public accession or DOI and license, and provide the commands and input versions needed to reproduce the paper figures and benchmark results. The large raw and converted datasets are not in this repository.
- Add accurate paper author and citation metadata once the manuscript title and author list are final. Archive the exact released code version and record its DOI in the paper.
- Verify a clean install, local server startup, the Claude extension, the live documentation site, and all public privacy/support links before directory submission.
