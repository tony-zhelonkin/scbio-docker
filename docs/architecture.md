# Architecture

System design of the scbio-docker container substrate: how the image is layered, how dev containers wire it up, and how the Python and R environments are organized. This is the overview — see [build.md](build.md), [environments.md](environments.md), and [devcontainer.md](devcontainer.md) for detail.

Throughout, the base image tag is written version-agnostically as `scdock-r-dev:$(cat VERSION)`.

## Image layering

The runtime image is a single multi-stage build from `docker/base/Dockerfile`.

```
ubuntu:22.04
  └─ builder stage        compile R from source, build CLI bio tools, resolve venvs
       └─ scdock-r-dev     runtime: only final artifacts copied in (true ~20GB)

greenleaflab/archr:1.0.3-base-r4.4   (official, legacy sidecar — R 4.4 + ArchR 1.0.3)
       └─ scdock-r-archr   optional wrapper (removes rstudio user, adds devuser)
```

- **`scdock-r-dev`** is the default and only actively developed image: R 4.5.3 (built from source with Cairo/BLAS/LAPACK/R-shlib) + Bioconductor 3.22, Python 3.11 (deadsnakes), and compiled CLI bio tools (samtools/bcftools/htslib 1.21, bedtools 2.31.1, MACS3, scIBD).
- **ArchR is deprecated.** The upstream image has been unmaintained for over two years; the field has moved to per-language stacks (Seurat/Signac in R, scanpy/snapATAC in Python). The `scdock-r-archr` wrapper (built from `.devcontainer/Dockerfile.archr-wrapper` via `scripts/build-archr-wrapper.sh`) exists only so the legacy service uses the same `devuser` layout. Treat it as on-path-to-removal.

Build modes and flags are covered in [build.md](build.md).

## Dev containers

`scripts/init-container.sh` (symlinked as `init-project.sh`) renders `templates/devcontainer/` into a project directory. The rendered `docker-compose.yml` defines two services:

| Service | Image | State | How to start |
|---------|-------|-------|--------------|
| `dev-core` | `scdock-r-dev` | default | `docker compose up` |
| `dev-archr` | `scdock-r-archr` | legacy, gated behind `profiles: ["archr"]` | `docker compose --profile archr up` |

Both run as `${LOCAL_UID:-1000}:${LOCAL_GID:-1000}`, share the project `.env`, mount a `tmpfs` `/tmp`, and apply CPU/memory/pids limits. VS Code remaps `devuser` (1000) to the invoking user via `updateRemoteUserUID`. See [devcontainer.md](devcontainer.md) for the full template, placeholders, mount model, and postStart hooks.

## Python environments

A single fully-resolved base venv plus on-demand layered venvs, rather than several full environments.

- **Base venv** `/opt/venvs/base` — preinstalled during the build from `docker/requirements/base.txt` (scanpy, anndata, scvi-tools, cellrank, scvelo, muon, MACS3, radian, jupyter/ipykernel, …).
- **Layered venvs** `squid` / `atac` / `comms` — created on first use with `python3.11 -m venv --system-site-packages` from `docker/requirements/{name}.txt`. `--system-site-packages` means each layer inherits the base stack and adds only its specialization, so it stays a few GB instead of tens.

Switching helpers (`usepy`, `py-base`, etc.) and per-layer package lists live in [environments.md](environments.md).

## R environment

**Two-tier library model:**

| Tier | Path | Writable | Contents |
|------|------|----------|----------|
| System | `/usr/local/lib/R/library` | read-only (root) | ~80 core packages resolved at build time |
| User | `~/R/x86_64-pc-linux-gnu-library/4.5` | writable (devuser) | runtime installs, takes precedence |

Core packages are installed at build time by `docker/base/R/install_core.R`.
Although `install_renv_project.R` can restore `/opt/settings/renv.lock`, the
Dockerfile does not copy the repository lockfile, so builds resolve packages
and snapshot the result afterward. The manifest at
`/opt/settings/R-packages-manifest.csv` records the resolved environment but
does not pin the next build. Runtime `install.packages()` /
`BiocManager::install()` calls land in the writable user library with no sudo
needed — the expected "installation paths not writeable" notice about the
read-only system tree is harmless. See [build.md](build.md) for source-specific
version constraints.

**R startup:** `.devcontainer/.Rprofile` is interactive- and VS Code-aware. It enables httpgd for in-editor plotting only when running inside VS Code, keeping non-interactive scripts clean.

The core package set and runtime-install workflow are detailed in [environments.md](environments.md).

## Sanity check

`scripts/poststart_sanity.sh` runs on container start (wired into the postStart chain) and validates the environment: R and Python resolve, httpgd is available, venv paths exist, and the user R library is writable. It fails loudly if the container is misconfigured before you start work.

## Image size

The build is optimized so Docker reports and disk usage stay close to the true ~20GB footprint.

- **Multi-stage build.** Compilation (R from source, CLI bio tools, venv resolution) happens in a builder stage; only final artifacts are copied into a fresh runtime image, so build caches and intermediate layers never bloat the shipped image. Single-stage builds inflated Docker's reported size via layer accounting even when real filesystem use was ~20GB.
- **TinyTeX, not full TeX.** A minimal LaTeX install covers Quarto/rmarkdown PDF rendering without the multi-GB TeX Live distribution.
- **Base venv + on-demand layers.** One resolved base venv instead of several full environments; `squid`/`atac`/`comms` are layered with `--system-site-packages` at runtime (see above).
- **~80 core R packages + runtime installs.** The image ships the common ~95% of workflows; heavy annotation packages (BSgenome.*, EnsDb.*, org.*.eg.db) and anything niche install at runtime into the user library.
- **Build tools retained intentionally.** Compilers stay in the runtime image so R/Python packages can be built at runtime.

Verify actual on-disk usage from inside the container rather than trusting `docker images`:

```bash
docker run --rm scdock-r-dev:$(cat VERSION) bash -lc \
  'du -hsx /* 2>/dev/null | sort -h | tail -n 10'
```
