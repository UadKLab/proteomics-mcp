# Releasing proteomics-mcp

Releases are published to [PyPI](https://pypi.org/project/proteomics-mcp/) by the [`release`](.github/workflows/release.yml) GitHub Actions workflow using [Trusted Publishing](https://docs.pypi.org/trusted-publishers/). No API tokens are stored anywhere. A release is started by pushing a version tag. The workflow then runs these jobs in order:

1. **build**: version check against the tag, sdist + wheel built once, `.mcpb` packed.
2. **verify**: the built wheel is tested on Linux, macOS and Windows with Python 3.11–3.14.
3. **publish to TestPyPI**
4. **publish to PyPI**: waits for a maintainer to approve the `pypi` deployment.
5. **GitHub release**: notes from `docs/releases/<tag>.md`, with sdist, wheel and `.mcpb` attached.

The GitHub release is created only after the PyPI upload succeeds.

## One-time setup

These settings live outside the repository and can only be changed by maintainers. Tracked in [#8](https://github.com/UadKLab/proteomics-mcp/issues/8).

### Trusted publishers (configured)

| Index | Owner | Repository | Workflow | Environment |
|---|---|---|---|---|
| [PyPI](https://pypi.org/manage/project/proteomics-mcp/settings/publishing/) | `UadKLab` | `proteomics-mcp` | `release.yml` | `pypi` |
| [TestPyPI](https://test.pypi.org/manage/project/proteomics-mcp/settings/publishing/) | `UadKLab` | `proteomics-mcp` | `release.yml` | `testpypi` |

The workflow file name must stay `release.yml`, and the environment names must match. Otherwise the upload is rejected.

### GitHub environments

Create both under *Settings → Environments* **before the first workflow run**. If an environment is missing, GitHub creates it on first use without any protection rules.

- **`pypi`**:
  - Required reviewers: at least one maintainer.
  - Deployment branches and tags: *Selected*, tag pattern `v*` only.
- **`testpypi`**:
  - Deployment branches and tags: *Selected*, tag pattern `v*` and branch `main`. The branch is needed for manual TestPyPI dry runs.
  - Reviewers are optional.

### Branch protection on `main`

Use a ruleset under *Settings → Rules*, or classic branch protection:

- Require a pull request before merging.
- Require status checks to pass. The jobs from the `tests` workflow are `build`, `test (<os>, <python>)` for `ubuntu-latest`, `macos-latest` and `windows-latest` with Python 3.11, 3.12, 3.13 and 3.14, and `extension launch (macos-latest)` / `extension launch (windows-latest)`. Optionally, also require the `docs` workflow's `build` job. Both workflows have a job named `build`, so pick the source workflow when you add it.
- Block force pushes and deletions.

### Merge strategy

Under *Settings → General → Pull Requests*, allow **merge commits only** and disable squash and rebase merging. This preserves every commit from an issue branch on `main`, and each merge commit links its pull request. Branches can be deleted after merge.

### PyPI owners

@marcoreverenna owns the PyPI and TestPyPI projects. To add a second owner (planned: @sviglios), on each index go to *Your projects → proteomics-mcp → Collaborators*, invite their username with role **Owner**, and ask them to accept the invitation. Every owner needs two-factor authentication. A lab PyPI organization is optional and can be set up later.

## Release checklist

Replace `X.Y.Z` with the new version.

1. **Prepare.** Open an issue "Prepare X.Y.Z release" and a branch `chore/release-X.Y.Z`.
2. **Bump the version** in every file that carries it:
   - `pyproject.toml` (`[project].version`, the single source of truth)
   - `manifest.json` (`version`)
   - `CITATION.cff` (`version`)
3. **Release notes.** Add `docs/releases/vX.Y.Z.md`. The workflow uses it as the GitHub release text.
4. **Changelog.** In `CHANGELOG.md`, move the *Unreleased* entries under a new `X.Y.Z` heading with the release date.
5. **Check:**
   ```bash
   python3 scripts/check_version.py --tag vX.Y.Z
   ```
6. **Merge.** Open a pull request, wait for CI, and merge it with a merge commit.
7. **Optional TestPyPI dry run.** Run *Actions → release → Run workflow* on `main`. This builds, verifies and publishes to TestPyPI only, and fails early if tag `vX.Y.Z` already exists. Then check the TestPyPI install (see below).
8. **Tag** the merge commit on `main`:
   ```bash
   git switch main && git pull --ff-only
   git tag -a vX.Y.Z -m "proteomics-mcp X.Y.Z"
   git push origin vX.Y.Z
   ```
9. **Approve.** When *verify* and *publish to TestPyPI* have passed, approve the `pypi` deployment in the workflow run. The *build* job summary lists the SHA-256 of the wheel, sdist and `.mcpb`.
10. **Check the result.** The workflow creates the GitHub release. Check that it lists the wheel, sdist and `.mcpb`, and that its text matches `docs/releases/vX.Y.Z.md`.
11. **Verify the PyPI install** in a clean environment:
    ```bash
    uvx --refresh proteomics-mcp@X.Y.Z --help
    python3 -m venv /tmp/pmcp && /tmp/pmcp/bin/pip install proteomics-mcp==X.Y.Z && /tmp/pmcp/bin/proteomics-mcp --help
    ```
12. **Track.** Tick the release items in the tracking issue.

### TestPyPI verification

TestPyPI hosts only `proteomics-mcp`; its dependencies come from PyPI:

```bash
uvx --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ --from proteomics-mcp proteomics-mcp --help
```

Add `==X.Y.Z` after `--from proteomics-mcp` to test a specific version.

### If a release run fails

- **Before *publish to PyPI*:** nothing has reached PyPI. Fix the problem on `main` through a pull request. Then delete the tag (`git push --delete origin vX.Y.Z`, then `git tag -d vX.Y.Z`) and tag the fixed commit again. TestPyPI may already hold files for `X.Y.Z`; the TestPyPI step skips existing files (`skip-existing`).
- **After *publish to PyPI*:** the version is on PyPI and cannot be re-uploaded. If only *GitHub release* failed, re-run that job from the workflow run page.
- **A dry run followed by a tag** leaves the dry-run files on TestPyPI, because the tag run skips existing files there. PyPI gets the tag-built files; their SHA-256 values are in the tag run's summary.

## Yanking a release

PyPI never accepts a re-upload of a version, even after deletion. To withdraw a broken release:

1. On PyPI, go to *Your projects → proteomics-mcp → Manage → Releases → X.Y.Z → Options → Yank*, and give a short reason. Do the same on TestPyPI if needed.
2. A yanked version stays downloadable for anyone who pins it exactly (`==X.Y.Z`), but installers skip it otherwise. Prefer yanking over deleting, which breaks existing pins and lock files.
3. Release the fix as a new version (e.g. `X.Y.Z+1`) following the checklist above.
4. Edit the GitHub release text to point to the fixed version, and add a note under that version in `CHANGELOG.md`.

## History: manual upload of 0.1.8

0.1.8 is the version cited in the manuscript. It was published once, by hand, before this workflow existed (decision 6 in [#1](https://github.com/UadKLab/proteomics-mcp/issues/1)). Details are in [#9](https://github.com/UadKLab/proteomics-mcp/issues/9).

- **What and when:** on 2026-10-12, @marcoreverenna uploaded the wheel and sdist attached to the GitHub release [`v0.1.8`](https://github.com/UadKLab/proteomics-mcp/releases/tag/v0.1.8) (tag → `ef7be33`) unchanged with `twine`, first to TestPyPI and then to PyPI.
- **Credentials:** short-lived, account-scoped API tokens, because the projects did not exist yet. They were revoked after the upload.
- **Integrity:** the SHA-256 values match on the GitHub release, TestPyPI and PyPI:

| File | SHA-256 |
|---|---|
| `proteomics_mcp-0.1.8-py3-none-any.whl` | `5f9fe2c04690b68d387eb93134f81c6b5555baddaaa890f62aaa8f3ef9747493` |
| `proteomics_mcp-0.1.8.tar.gz` | `f887e67a54bb3f29a453cb7fbd70a78724c40d53057a294c8e896264793e1ecb` |

All later releases go through the workflow above. Do not upload by hand.

## Contributing

1. **Issue:** every change starts with an issue that states the problem, the scope and the acceptance criteria.
2. **Branch:** one branch per issue, named by type, e.g. `fix/…`, `chore/…`, `ci/…` or `docs/…`.
3. **Commits:**
   - Small, one concern per commit.
   - [Conventional Commits](https://www.conventionalcommits.org/) with a scope and an imperative subject of at most 72 characters.
   - A body that explains *why*, and a `Refs #N` footer.
   - Each commit leaves the tests passing.
4. **Pull request:** open a draft PR early, push as you go, and describe the verification (commands and output). Link the issue with `Closes #N`.
5. **Merge:** a maintainer reviews and merges with a merge commit.
6. **No history rewriting after push:** no amend, rebase or force push. Address review feedback with new commits, and update a branch by merging `main` into it.
7. **Scope:** changes to scientific behaviour (adapters, scoring, annotation, metrics, rendering) need their own issue and maintainer review. Keep packaging, CI and documentation changes separate from them.
