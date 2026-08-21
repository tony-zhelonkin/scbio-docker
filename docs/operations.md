# Operations

Runbook for building, running, and operating the scbio-docker environment outside of the VS Code Dev Container flow. For the container-first workflow see [devcontainer.md](devcontainer.md); for the image internals see [environments.md](environments.md).

Throughout, the base image tag is written as `scdock-r-dev:$(cat VERSION)`. Substitute the current tag when copying commands.

## Building

Building is covered in full in [build.md](build.md). In short:

```bash
scripts/build.sh                       # generic build (devuser:1000, shareable)
scripts/build.sh --github-pat ghp_...  # avoids GitHub rate limits during build
scripts/build.sh --personal            # bakes your UID/GID/USER/GROUP
```

The optional ArchR wrapper image (`scdock-r-archr`) is built with `scripts/build-archr-wrapper.sh`. ArchR is legacy — see [ArchR (legacy)](#archr-legacy).

## Running the image manually

For a one-off shell without VS Code, run the image directly as your own UID/GID so files land with the right ownership.

```bash
docker run --rm -it \
  -u $(id -u):$(id -g) \
  -v /path/to/project:/workspaces/project \
  -v /path/to/data:/workspaces/project/00_data/raw:ro \
  --memory=450g --cpus=50 \
  scdock-r-dev:$(cat VERSION) bash
```

Notes:
- `-u $(id -u):$(id -g)` runs as you; the generic `devuser:1000` in the image is transparently overridden.
- `WORKDIR` is `/workspaces/project`; mount your project there.
- Mount read-only data with `:ro`.
- `--memory` / `--cpus` cap resource use on shared machines.

## Docker Compose

The rendered `.devcontainer/docker-compose.yml` defines two services: `dev-core` (default) and `dev-archr` (gated behind the `archr` profile). Compose reads UID/GID and workspace path from the environment (or the `.env` written by `init-container.sh`).

```bash
export LOCAL_UID=$(id -u)
export LOCAL_GID=$(id -g)
export WORKSPACE_FOLDER=$PWD

# Start the base service in the background
docker compose -f .devcontainer/docker-compose.yml up -d dev-core

# Open a shell inside it
docker compose -f .devcontainer/docker-compose.yml exec dev-core bash

# Tear down
docker compose -f .devcontainer/docker-compose.yml down
```

On a shared machine, different users override UID/GID per invocation:

```bash
LOCAL_UID=1001 LOCAL_GID=1001 docker compose -f .devcontainer/docker-compose.yml up -d dev-core
```

## ArchR (legacy)

> ArchR is **deprecated / on the path to removal**. Upstream has been unmaintained for over two years; the field has moved to per-language stacks (Seurat/Signac in R, scanpy/snapATAC in Python). Only use `dev-archr` for reproducing older ArchR analyses.

The `dev-archr` service uses the wrapper image `scdock-r-archr` (R 4.4 + ArchR 1.0.3) and is gated behind a Compose profile so it never starts by default.

```bash
# Compose: bring up the ArchR service explicitly
docker compose -f .devcontainer/docker-compose.yml --profile archr up -d dev-archr
docker compose -f .devcontainer/docker-compose.yml --profile archr exec dev-archr bash
```

Or run the official upstream image directly:

```bash
docker run --rm -it \
  -u $(id -u):$(id -g) \
  -v /path/to/project:/workspaces/project \
  greenleaflab/archr:1.0.3-base-r4.4 bash
```

In VS Code, switch by setting `"service": "dev-archr"` in `.devcontainer/devcontainer.json`, then **Dev Containers: Rebuild and Reopen in Container**.

Verify inside the container:

```r
library(ArchR)
packageVersion("ArchR")   # 1.0.3
```

## Python environment switching

Environment layout and package pins live in [environments.md](environments.md). Quick reference:

```bash
usepy base    # core single-cell stack (default, already active)
usepy squid   # spatial (squidpy, spatialdata) — created on first use
usepy atac    # scATAC (snapatac2, episcanpy)
usepy comms   # cell communication / GRN

which python && python -V   # confirm the active venv
```

Layered venvs are created on demand from `docker/requirements/{name}.txt` with `--system-site-packages`, so they inherit the base stack.

## Persistent R sessions (radian + tmux)

Radian is the R terminal (`radian` or the `r-base` wrapper). For SSH work, run it inside `tmux` so the session survives disconnects.

```bash
tmux new-session -s analysis radian   # start radian in a named session
# detach with: Ctrl+B then D

tmux attach -t analysis               # re-attach later
tmux ls                               # list sessions
tmux kill-session -t analysis         # stop it
```

In VS Code, start radian in a tmux-backed terminal and send code to it from `.R` files; plots render via httpgd. R terminal settings are documented in [devcontainer.md](devcontainer.md).

## Troubleshooting

### Permission issues on the workspace

Files written inside the container should be owned by you. If writes fail, confirm the runtime UID and mount ownership:

```bash
id -u; id -g; id -gn
ls -ld /workspaces/project
```

Fix ownership on the host if needed:

```bash
sudo chown -R $(id -u):$(id -g) /path/to/project
```

Always pass `-u $(id -u):$(id -g)` (docker run) or set `LOCAL_UID`/`LOCAL_GID` (compose). In VS Code, `updateRemoteUserUID: true` handles this automatically.

### R user library not writable

Runtime R installs go to the writable user library, which takes precedence over the read-only system library. Check:

```r
.libPaths()
# [1] "/home/devuser/R/x86_64-pc-linux-gnu-library/4.5"   # writable, first
# [2] "/usr/local/lib/R/library"                          # system, read-only
file.access(.libPaths()[1], 2) == 0                        # TRUE = writable
```

The `Installation paths not writeable, unable to update packages` warning from `BiocManager::install()` is **expected and harmless** — your package installed to the user library; only the read-only system library can't be touched. Suppress with `update = FALSE`.

### Missing ArchR

`library(ArchR)` fails in `dev-core` because ArchR is not in the base image. Switch to the `dev-archr` service or the official image — see [ArchR (legacy)](#archr-legacy).

### GitHub rate limits during build

Package installs that pull from GitHub can hit anonymous rate limits. Supply a PAT:

```bash
scripts/build.sh --github-pat ghp_your_token_here
```

The PAT is passed via BuildKit `--secret` and is not baked into image layers.

### Sanity check

The post-start sanity script validates the environment (Python, R, httpgd, scanpy import, default venv, radian):

```bash
scripts/poststart_sanity.sh
```

## See also

- [build.md](build.md) — image build modes and flags
- [environments.md](environments.md) — R/Python stacks and package pins
- [devcontainer.md](devcontainer.md) — VS Code Dev Container workflow and settings
- [isolation.md](isolation.md) — filesystem/resource hardening for Compose
- [changelog.md](changelog.md) — version history
