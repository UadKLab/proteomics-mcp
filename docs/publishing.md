# Publishing

## GitHub release and Claude Desktop extension

The `release` workflow runs when a `v*` tag is pushed. It builds the Python wheel and source archive, builds the `.mcpb` extension, and creates a GitHub release. The release description comes from `docs/releases/<tag>.md`; create this file before pushing a new tag.

Install the `.mcpb` attachment in Claude Desktop to test the packaged extension. A GitHub release does not automatically submit it to Anthropic's directory. Before directory submission, also test installation on clean macOS and Windows systems and review the third-party privacy notices listed in `PRIVACY.md`.

## PyPI

PyPI is optional for the Claude Desktop extension. It makes `pip install proteomics-mcp` and `uvx proteomics-mcp serve` available to Python and MCP users.

The `publish-pypi` workflow is manual so a GitHub release does not automatically upload to PyPI. It accepts only an existing, non-draft release tag and checks that the tag matches the package version. To enable it:

1. Make the source repository and its privacy documentation public after reviewing tracked examples for local paths and data-sharing permissions. This is especially important for a paper and the Claude directory.
2. In GitHub repository settings, create a `pypi` environment, restrict it to release tags, and require manual approval by a maintainer. Required reviewers may not be available while the repository is private, depending on the GitHub plan.
3. In [PyPI's publishing settings](https://pypi.org/manage/account/publishing/), configure a pending Trusted Publisher for project `proteomics-mcp`, owner `UadKLab`, repository `proteomics-mcp`, workflow `publish-pypi.yml`, and environment `pypi`. A pending publisher does not reserve the package name until the first upload.
4. After the new GitHub release is built and tested, run the `publish-pypi` workflow from GitHub Actions with its tag, such as `v0.1.8`. Approve the `pypi` deployment when prompted.
5. Verify the resulting PyPI page and install the exact version in a clean environment with `pip install proteomics-mcp==0.1.8`.

Before the first upload, verify the name is still available and review the build and `twine check` output. PyPI versions cannot be replaced after upload. The GitHub release workflow tests the tag and validates both Python distributions before creating the release.

## Documentation site

MkDocs builds the documentation from `docs/` with `mkdocs.yml`. Run `uvx --from mkdocs==1.6.1 mkdocs build --strict` locally. The `docs` workflow checks the build on pushes and pull requests; it deploys to GitHub Pages only after the repository is public. At that point, set the Pages publishing source to **GitHub Actions** in repository settings and verify `https://uadklab.github.io/proteomics-mcp/` before linking it from PyPI or the Claude directory.

## Paper and public-repository checklist

- Keep public history limited to the sanitized snapshot. Do not reintroduce old local paths or generated files from private archives.
- Confirm that each example dataset can be shared, document its public accession or DOI and license, and provide the commands and input versions needed to reproduce the paper figures and benchmark results. The large raw and converted datasets are not in this repository.
- Add accurate paper author and citation metadata once the manuscript title and author list are final. Archive the exact released code version and record its DOI in the paper.
- Verify a clean install, local server startup, the Claude extension, the live documentation site, and all public privacy/support links before directory submission.
