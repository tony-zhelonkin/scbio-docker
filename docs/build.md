# Build Guide

How to build the `scdock-r-dev` base image (and the optional ArchR wrapper) with `scripts/build.sh`. The image tag is read from the [`VERSION`](../VERSION) file, so nothing here hardcodes a version.

## Entry point

`scripts/build.sh` is the preferred way to build. It wraps a multi-stage `docker build` of [`docker/base/Dockerfile`](../docker/base/Dockerfile), tags the result `scdock-r-dev:$(cat VERSION)`, streams output to `build.log`, and prints post-build next steps.

```bash
scripts/build.sh                       # generic build, tag from VERSION
scripts/build.sh --github-pat ghp_...  # with a GitHub token (see below)
```

Direct build (equivalent, no wrapper conveniences):

```bash
DOCKER_BUILDKIT=1 docker build . \
  -f docker/base/Dockerfile \
  --build-arg USER_ID=1000 --build-arg GROUP_ID=1000 \
  --build-arg USER=devuser --build-arg GROUP=devgroup \
  -t "scdock-r-dev:$(cat VERSION)"
```

## Build modes

The script bakes a user/group into the image. Three modes control which identity that is.

| Mode | Invocation | Baked identity | Shareable |
|------|-----------|----------------|-----------|
| **Generic** (default) | `scripts/build.sh` | `devuser:1000` / `devgroup:1000` | Yes — publish to a registry, everyone remaps at runtime |
| **Personal** | `scripts/build.sh --personal` | Your `id -u`/`id -g`/`$USER`/`id -gn` | No — only useful to you |
| **Custom** | `scripts/build.sh --user-id 2000 ...` | Whatever IDs/names you pass | Depends on the values |

Passing `--user-id` switches the mode to `custom`; `--personal` fills all four values from your current login.

**Prefer the generic build.** A single `devuser:1000` image works for the whole team and is remapped to each person's real UID at runtime (see below). Only build `--personal` when you need files owned by your literal UID inside a standalone `docker run` with no remapping.

## Flags

| Flag | Effect |
|------|--------|
| `--github-pat TOKEN` | GitHub token, passed to BuildKit as `--secret id=github_pat` (kept out of image layers). Avoids GitHub API rate limits during package installs. |
| `--tag TAG` | Override the output tag (default `scdock-r-dev:$(cat VERSION)`). |
| `--user-id UID` | Custom UID; sets mode to `custom` (default `1000`). |
| `--group-id GID` | Custom GID (default `1000`). |
| `--user NAME` | Custom username (default `devuser`). |
| `--group NAME` | Custom group name (default `devgroup`). |
| `--personal` | Build with your current UID/GID/user/group. |
| `-y`, `--yes` | Skip the confirmation prompt (for unattended/agentic builds). |
| `--help` | Print usage and exit. |

## GitHub PAT and rate limits

Some R/GitHub package installs pull from GitHub during the build. Without a token you may hit anonymous API rate limits and the build can fail. Provide a token either way:

```bash
scripts/build.sh --github-pat ghp_xxxxx
# or via environment (the script honors $GITHUB_PAT):
export GITHUB_PAT=ghp_xxxxx
scripts/build.sh
```

The script writes the token to a temp file and passes it as a BuildKit `--secret` (`id=github_pat`), then deletes the temp file. The secret is mounted only for the RUN steps that need it, so the token never lands in an image layer or in `build.log`.

## UID/GID remapping strategy

The generic image ships as `devuser:1000`. That identity is remapped to your real user at runtime — the image itself stays untouched and shareable:

| Runtime | How remapping happens |
|---------|----------------------|
| **VS Code Dev Containers** | `updateRemoteUserUID: true` in `devcontainer.json` remaps `1000` → your host UID automatically. |
| **Docker Compose** | The compose service runs as `${LOCAL_UID:-1000}:${LOCAL_GID:-1000}`; `init-container.sh` writes your UID/GID into `.env`. |
| **`docker run`** | Pass `-u $(id -u):$(id -g)` explicitly. |

This is why generic is the default: build once, everyone runs it as their own user with correctly owned files. See [devcontainer.md](devcontainer.md) and [operations.md](operations.md) for the runtime side.

## R package resolution and build records

R package resolution is not deterministic. `install_renv_project.R` would
restore `/opt/settings/renv.lock` if that file existed at the start of the
build, but the Dockerfile does not copy the repository lockfile into the image.
Every build therefore installs the R stack from `install_core.R` and snapshots
the result afterward. The repository lockfile contains only `renv` and is not a
build input.

Package sources have mixed version guarantees:

| Source | Current constraint |
|--------|--------------------|
| CRAN | `install_core.R` configures the RSPM `2026-04-15` snapshot as its default, but many calls explicitly use `cloud.r-project.org` and bypass that snapshot. |
| Bioconductor | Release `3.22` is selected; individual package versions are not locked. |
| r-universe | Package versions float. |
| GitHub | Repositories generally float; `bulkiRNA` is the explicit commit-pinned exception. |

The generated lockfile and manifests record what a particular build resolved;
they do not make the next build reproduce it. They can be extracted for audit:

```bash
TAG="scdock-r-dev:$(cat VERSION)"
CID=$(docker create "$TAG")
mkdir -p build-artifacts
docker cp "$CID":/opt/settings/renv.lock build-artifacts/renv.lock
docker cp "$CID":/opt/settings/R-packages-manifest.csv build-artifacts/R-packages-manifest.csv
docker cp "$CID":/opt/settings/install_failures.csv build-artifacts/install_failures.csv
docker rm "$CID"
```

The required-package contract in `install_core.R` fails the build when its
small required floor is absent or misidentified. It checks presence after
resolution; it does not pin versions or guarantee that the wider package set
matches an earlier image.

## Sanity check

Validate a freshly built image before using it:

`poststart_sanity.sh` is **not baked into the image** — it is mounted at
runtime (the devcontainer copies it into the project). To run it against a bare
image, mount the repo:

```bash
TAG="scdock-r-dev:$(cat VERSION)"
docker run --rm -v "$PWD:/repo" -w /repo "$TAG" bash scripts/poststart_sanity.sh
# interactive smoke test:
docker run --rm -it "$TAG" bash
```

Docker may report a large image `Size` because of layer accounting; the real on-disk filesystem is ~20GB. To inspect the largest paths:

```bash
docker run --rm "$TAG" du -hsx /* 2>/dev/null | sort -h | tail -n 10
```

## ArchR wrapper (legacy)

ArchR is a deprecated sidecar. The official upstream image `greenleaflab/archr:1.0.3-base-r4.4` (R 4.4 + ArchR 1.0.3) uses an `rstudio` user, which clashes with the `devuser` convention. `scripts/build-archr-wrapper.sh` builds a thin wrapper ([`.devcontainer/Dockerfile.archr-wrapper`](../.devcontainer/Dockerfile.archr-wrapper)) that removes `rstudio` and recreates `devuser`, so UID handling matches the base image.

```bash
scripts/build-archr-wrapper.sh              # generic (devuser:1000)
scripts/build-archr-wrapper.sh --personal   # personal build
scripts/build-archr-wrapper.sh --tag scdock-r-archr:custom
```

It accepts the same identity flags (`--user-id`, `--group-id`, `--user`, `--group`, `--personal`, `-y`, `--tag`), produces `scdock-r-archr:$(cat VERSION)`, and logs to `build-archr-wrapper.log`. In compose the `dev-archr` service is gated behind the `archr` profile (`docker compose --profile archr up`).

ArchR upstream has been unmaintained for over two years; the field has moved to per-language stacks (Seurat/Signac in R, scanpy/snapATAC in Python). Treat the wrapper as legacy and on the path to removal — see [roadmap.md](roadmap.md).

## Related docs

- [environments.md](environments.md) — what R/Python packages land in the image
- [architecture.md](architecture.md) — multi-stage layering
- [devcontainer.md](devcontainer.md) — running the image via VS Code / compose
- [changelog.md](changelog.md) — version history and release notes
