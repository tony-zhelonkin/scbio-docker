# Single-Cell Docker Dev Environment

![Docker Image Version](https://img.shields.io/badge/Docker-v0.5.1-blue?style=flat-square)
![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)

A Docker-based development environment for bioinformatics, particularly single-cell RNA-seq analyses. This repository is structured for use with Visual Studio Code's [**Dev Containers**](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) extension, enabling seamless local or remote development against powerful server resources.

## What's New in v0.5.1 (Multi-Stage Build)

** Size Reduction: 500GB → 20GB (Docker-reported)**

- **Multi-stage build**: Completely discards build artifacts, no layer bloat
- **Build-essential preserved**: Should still compile R/Python packages at runtime
- **Same functionality**: All features from v0.5.0, but with size efficiency

**Previous Optimizations (v0.5.0):**
- Aggressive cache cleanup: renv, pip, and build artifacts removed
- TinyTeX: Lightweight TeX distribution instead of full texlive
- Layered Python venvs: Base venv + runtime-created specialized venvs
- Core R packages approach: ~80 essential packages pre-installed, others at runtime on-demand
- R 4.5 + Bioconductor 3.21: Updated for anndataR and modern package support
- Official ArchR image: Use `greenleaflab/archr:1.0.3-base-r4.4` instead of custom build, custom build discarded

**Key Benefits:**
- ✅ **~20GB final image** (no layer accounting issues)
- ✅ Same functionality - core packages cover many standard use cases
- ✅ **Runtime package installation fully supported** (build tools preserved)
- ✅ Improved reproducibility with layered venv approach
- ✅ No R version conflicts (dev-core: R 4.5, ArchR: R 4.4)

**Project Templates:**
- ✅ `init-project.sh` script for quick project scaffolding
- ✅ 4 templates: `basic-rna`, `multimodal`, `archr-focused`, `example-DMATAC`
- ✅ Universal `.vscode/settings.json` with Python REPL + R configuration

**Build:** Use `scripts/build.sh` for the multi-stage build. See [docs/build.md](docs/build.md). Note: build assets moved to `docker/base/Dockerfile`.

> **Inspiration**: This setup is **strongly inspired by** [Rami Krispin’s vscode-r repository](https://github.com/RamiKrispin/vscode-r). Special thanks for the excellent reference on R + VS Code configurations!

---

## AI Integration (Optional)

AI features are provided via the external SciAgent‑toolkit. The `dev-claude-integration` and `dev-gpt-codex-integration` branches are deprecated.

- MCP/agent configuration: https://github.com/tony-zhelonkin/SciAgent-toolkit

### Archived branches

The deprecated AI branches have been archived as tags for reference:

- `archived/dev-claude-integration`
- `archived/dev-gpt-codex-integration`

Please use `dev` for active work and open PRs into `main` for releases.

---

## Table of Contents
1. [Overview](#overview)
2. [Features](#features)
3. [Usage Guide](#usage-guide)
4. [Repository Structure](#repository-structure)
5. [Getting Started](#getting-started)
    - [Prerequisites](#prerequisites)
    - [Build the Image](#build-the-image)
    - [Run the Container](#run-the-container)
    - [Using with VS Code Remote Containers](#using-with-vs-code-remote-containers)
6. [Typical Workflow for New Projects](#typical-workflow-for-new-projects)
7. [Included Tools & Packages](#included-tools--packages)
8. [Troubleshooting](#troubleshooting)
9. [Contributing](#contributing)
10. [License](#license)
11. [Acknowledgments](#acknowledgments)

---

## Overview

This repository contains Dockerfiles and configuration files for creating a **versatile bioinformatics environment** that focuses on single-cell RNA-seq and epigenomics workflows. It packages:

- **R** (v4.5.0) and comprehensive single-cell and genomic libraries (Seurat, Signac, Bioconductor ecosystem, etc.)
  - ~80 core packages pre-installed (Seurat, edgeR, limma, clusterProfiler, GSVA, etc.)
  - Runtime installation supported for additional packages
  - **ArchR**: Use official image on-demand (`greenleaflab/archr:1.0.3-base-r4.4`)
- **Python** (v3.10) with scientific libraries (NumPy, SciPy, Pandas, scVI, Scanpy, etc.)
  - Base venv with core single-cell stack
  - Layered venvs for specialized tools (spatial, ATAC, cell communication)
- **Command-line tools** focused on epigenomics-friendly workflows: **samtools**, **bcftools**, **bedtools**, and selected extras (e.g., **scIBD** for scATAC doublet detection). Bulk aligners and pre-processing tools (STAR, BWA, Salmon, kallisto, Picard, FastQC/Trimmomatic, Trim Galore, featureCounts, etc.) are intentionally excluded to keep the image slim. Use pipeline-specific containers if you need them.

The Docker setup is integrated with **VS Code Remote Containers** to streamline development on remote compute nodes—ideal for large-scale single-cell data processing.

### Quick Start

```bash
# 1. Build base image (multi-stage, ~20GB true size)
scripts/build.sh

# 2. Pull official ArchR image
docker pull greenleaflab/archr:1.0.3-base-r4.4

# 3. Initialize new project
./init-project.sh ~/projects/my-analysis basic-rna

# 4. Open in VS Code and reopen in container
code ~/projects/my-analysis
# Cmd/Ctrl+Shift+P → "Dev Containers: Reopen in Container"
```

See [DEVOPS.md](DEVOPS.md) for complete build and operational instructions.

**Runtime Package Installation:** See [RUNTIME_INSTALL.md](RUNTIME_INSTALL.md) for installing additional R/Python packages at runtime.

---

## Features

- **R** with a rich set of CRAN and Bioconductor packages for single-cell and bulk RNA-seq analysis.
- **Python** environment (virtualenv) preloaded with popular data science packages, single-cell toolkits, and epigenomics libraries.
- **CLI Tools** focused on downstream and epigenomics analysis (samtools, bcftools, bedtools) plus **scIBD**. Heavy aligners and trimming/QC tools are excluded to reduce image size.
- **Quarto** for R Markdown, Jupyter, and scientific documentation workflows.
- **Automated Build** process (via .devcontainer/Dockerfile) including compilation of R, specialized R packages, and Python requirements.
- **VS Code** integration with recommended extensions for R, Python, Docker, Markdown, and more.

---

## Usage Guide

### Building the image (with renv snapshot)

From the repository root (non-root runtime, nf-core style UID/GID passthrough baked in):

```bash
# First build: installs R packages and snapshots with renv
docker build . \
  -f docker/base/Dockerfile \
  --build-arg GITHUB_PAT=$GITHUB_PAT \
  --build-arg USER_ID=$(id -u) \
  --build-arg GROUP_ID=$(id -g) \
  --build-arg USER=$USER \
  --build-arg GROUP=$(id -gn) \
  -t scdock-r-dev:v0.4.0
```

- GITHUB_PAT is optional but recommended to avoid GitHub API rate limits during R package installs.
- USER_ID/GROUP_ID/USER/GROUP align the container user with your host user for mounted volume permissions.

### R package pinning with renv

- Purpose: renv captures exact R/CRAN/Bioconductor/GitHub package versions for deterministic rebuilds.
- First build (no lock present):
  - `install_R_packages.R` installs packages (CRAN snapshot + Bioc 3.20).
  - renv snapshots to `/opt/settings/renv.lock`.
  - A human-readable CSV is written to `/opt/settings/R-packages-manifest.csv`.
- After the first successful build, extract and commit the lock and manifest:

```bash
CID=$(docker create scdock-r-dev:v0.4.0)
docker cp $CID:/opt/settings/renv.lock ./renv.lock
docker cp $CID:/opt/settings/R-packages-manifest.csv ./R-packages-manifest.csv
docker rm $CID

git add renv.lock R-packages-manifest.csv
git commit -m "Pin R via renv; add manifest"
```

- Deterministic builds thereafter: uncomment the lock copy in `docker/base/Dockerfile`:

```Dockerfile
# COPY renv.lock /opt/settings/renv.lock
```

Rebuild and the image will restore via `renv::restore()` from the committed lockfile.

Tip: You can set a CRAN snapshot date with `RSPM_SNAPSHOT` (see `install_R_packages.R`) to reproduce historical CRAN states.

#### Handling GitHub PAT (optional)
- Export for the build session to avoid rate limits during GitHub installs:
```bash
export GITHUB_PAT=ghp_your_token_here
```
- Unset after builds if desired:
```bash
unset GITHUB_PAT
```

### Python environments (v0.5.0: Layered venvs)

**Base venv** (pre-installed):
- `/opt/venvs/base` - Core single-cell stack (scanpy, scvi-tools, muon, cellrank, scvelo, radian)
- Installed during build from `/opt/environments/base_requirements.txt`

**Layered venvs** (created at runtime with `--system-site-packages`):
- `/opt/venvs/squid` - Spatial transcriptomics (squidpy, spatialdata) - inherits base packages
- `/opt/venvs/atac` - scATAC-seq (snapatac2, episcanpy) - inherits base packages
- `/opt/venvs/comms` - Cell communication (liana, cellphonedb, scglue) - inherits base packages

**Creating and switching venvs:**

```bash
# First use automatically creates the venv
usepy squid      # Creates /opt/venvs/squid if needed, then activates
usepy atac       # Creates /opt/venvs/atac if needed, then activates
usepy comms      # Creates /opt/venvs/comms if needed, then activates
usepy base       # Switch back to base env
```

**Manual venv creation:**

```bash
# Create specific layered venv
create_layered_venv.sh squid squid_requirements.txt
create_layered_venv.sh atac atac_requirements.txt
create_layered_venv.sh comms comms_requirements.txt
```

**One-off commands:**

```bash
py-base python -V    # Uses base venv
```

**Verify active environment:**

```bash
which python && python -V && pip list | head
```

**Freeze for reproducibility:**

```bash
# After installing additional packages
pip freeze > requirements-project-frozen.txt
```

### Version pinning strategy

- R: renv.lock (commit it; enable restore via `COPY renv.lock /opt/settings/renv.lock`).
- Python: `.environments/*.txt` (pin versions there and commit).
- CLI tools: pinned by explicit versions/tags in `.devcontainer/Dockerfile`.
- Image: tag builds explicitly (e.g., `-t scdock-r-dev:v0.3`).

---

## Repository Structure

```bash
.
├── .devcontainer
│   ├── devcontainer.json             # VS Code Remote Containers configuration (docker-compose)
│   ├── docker-compose.yml            # Two services: base and archr
│   ├── Dockerfile                    # Primary Dockerfile (builds R, Python envs, tools)
│   ├── Dockerfile.archr              # ArchR variant (FROM base image)
│   ├── install_R_archr.R             # Dedicated ArchR installer (isolated lib)
│   ├── install_renv_project.R        # Clean renv init/restore runner
│   ├── install_httpgd.R              # httpgd install (CRAN-first, GitHub-fallback)
│   ├── install_quarto.sh             # Quarto installer
│   ├── install_R_packages.R          # R packages (CRAN/Bioc/GitHub); used with renv
│   └── .Rprofile                     # R profile (VS Code/httpgd; CRAN mirror; ArchR toggle)
├── .environments
│   ├── base_requirements.txt         # Base Python environment
│   └── squid_requirements.txt        # Alternative Python environment (squid)
├── README.md                         # This README
└── LICENCE.md
```

- **`docker/base/Dockerfile`** is the primary Dockerfile for building R, Python envs, and CLI tools.
- **`.devcontainer/install_R_packages.R`** installs R/Bioconductor/GitHub packages; paired with `renv` for pinning.
- **`docker/requirements/*.txt`** are the Python requirement sets installed into separate venvs inside the image.

---

## Architecture: Compose vs Dev Containers and the Image Stack

This repo supports two ways to run your development environment:

- Dev Containers with a single image (simple): `devcontainer.json` uses an `image:` directly.
- Dev Containers with Docker Compose (multi-service): `devcontainer.json` references `docker-compose.yml` and selects a `service:`.

Role of `.devcontainer/docker-compose.yml`:
- Defines two services that both mount the same workspace and run as your UID/GID:
  - `dev-core`: uses image `scdock-r-dev` (R+Python baseline)
  - `dev-archr`: uses image `scdock-r-archr` (adds ArchR in isolated R library)
- Bring up either or both containers, attach VS Code to one, and keep mounts/permissions consistent.

How this aligns with a single-image `devcontainer.json`:
- Instead of specifying `"image": "scdock-r-dev:…"`, you set:
  - `"dockerComposeFile": ".devcontainer/docker-compose.yml"`
  - `"service": "dev-archr"` (recommended default) or `"dev-core"`
- VS Code then starts the chosen service container with the same mounts. Switch by changing `service` or attach to the running container.

Image lineage (Docker layering; not an image-inside-image):

```text
ubuntu:22.04
  └─ scdock-r-dev:v0.4.1 (docker/base/Dockerfile)
       • System deps, R 4.4.2 built from source
       • R packages via install_R_packages.R, pinned with renv by install_renv_project.R
       • Python venvs: /opt/venvs/{base,squid,atac,comms} from docker/requirements/*.txt
       • CLI tools (samtools/bcftools/bedtools), Quarto
       • Wrappers: usepy, py-*, r-base, r-archr
       • httpgd installed via install_httpgd.R (CRAN-first, GitHub-fallback)

  └─ scdock-r-archr:v0.4.1 (.devcontainer/Dockerfile.archr; FROM scdock-r-dev:v0.4.1)
       • install_R_archr.R installs ArchR into ARCHR_LIB (separate R lib path)
       • Ensures macs2 compatibility (symlink if only macs3 exists)
       • Uses same r-base/r-archr wrappers (r-archr sets USE_ARCHR=1)
```

Runtime hooks and interop:
- `.Rprofile` (interactive-only): enables httpgd in VS Code; prepends `ARCHR_LIB` to `.libPaths()` when `USE_ARCHR=1`.
- `poststart_sanity.sh`: OK/NOT OK for Python, R, httpgd, Scanpy, default venv.
- Radian wrappers: `r-base` (base libs) and `r-archr` (sets `USE_ARCHR=1`).

Compose + Dev Containers flow (when using compose):

```text
Host VS Code
  └─ devcontainer.json
       └─ dockerComposeFile -> .devcontainer/docker-compose.yml
            ├─ service: dev-core  -> container from scdock-r-dev:v0.4.1
            └─ service: dev-archr -> container from scdock-r-archr:v0.4.1

Inside scdock-r-dev:
  - install_renv_project.R  (renv init/restore or snapshot)
  - install_httpgd.R        (httpgd install with fallback)
  - install_R_packages.R    (package set used by renv on first build)
  - Python venv setup       (base/squid/atac/comms)
  - poststart_sanity.sh     (run at container start)
```

When to prefer compose:
- You want both `base` and `archr` containers available with identical mounts/UIDs.
- You want to switch VS Code between them by toggling `service` or attaching to a running container.

If you prefer the simpler flow, keep using a single-image `devcontainer.json`. Compose is optional but convenient for multi-variant setups.

---

## Getting Started

### Prerequisites

- [Docker](https://docs.docker.com/get-docker/) installed on your local machine or remote server.
- (Optional) [Visual Studio Code](https://code.visualstudio.com/) with the [Remote Development](https://code.visualstudio.com/docs/remote/remote-overview) extensions installed if you plan to develop within VS Code.

> **Note**: If you plan to install packages from private GitHub repositories or use GitHub-limited APIs, you may need a **GitHub Personal Access Token**. Set it as a build argument `GITHUB_PAT` when you build the image.

### Build the Image

From the root of this repository, run:

```bash
docker build . \
    -f docker/base/Dockerfile \
    --build-arg GITHUB_PAT=<your_github_pat> \
    --build-arg R_VERSION_MAJOR=4 \
    --build-arg R_VERSION_MINOR=4 \
    --build-arg R_VERSION_PATCH=2 \
    -t scdock-r-dev:v0.2
```

- **`GITHUB_PAT`** (optional) is used to authenticate GitHub requests during R package installation (especially if you hit rate limits).
- Adjust **R_VERSION_*** build arguments as needed to override the default R version. But beware the build was only tested by myself for the R 4.4.2, and some of the R libraries depend on R version no older than 4.3-4.4

### Build the ArchR variant image

```bash
docker build /data1/users/antonz/pipeline/scbio-docker \
  -f /data1/users/antonz/pipeline/scbio-docker/.devcontainer/Dockerfile.archr \
  --build-arg BASE_IMAGE=scdock-r-dev:v0.4.0 \
  -t scdock-r-archr:v0.4.0
```

Sanity check:

```bash
docker run --rm scdock-r-archr:v0.4.0 bash -lc 'whoami && R -q -e "packageVersion(\"ArchR\")"'
```

## Typical Workflow for New Projects

After you’ve successfully built and tested your Docker image (e.g., `scdock-r-dev:v0.2`), you can quickly spin up a new analysis project by following these steps:

1. **Create a folder** for the analysis. For example:
   ```bash
   mkdir Yasmine-retroT && cd Yasmine-retroT
   ```

2. **Create `.devcontainer` and `.vscode` subfolders** inside your project directory:
   ```bash
   mkdir .devcontainer .vscode
   ```

3. Use docker-compose with two services (base and archr). Use the same mounts on both services.

```yaml
# .devcontainer/docker-compose.yml
version: "3.8"
services:
  base:
    image: scdock-r-dev:v0.4.0
    user: "${LOCAL_UID:-1000}:${LOCAL_GID:-1000}"
    working_dir: /workspaces/project
    volumes:
      - ${WORKSPACE_FOLDER}:/workspaces/project
      - /path/to/data:/workspaces/project/00_Data:ro
  archr:
    image: scdock-r-archr:v0.4.0
    user: "${LOCAL_UID:-1000}:${LOCAL_GID:-1000}"
    working_dir: /workspaces/project
    volumes:
      - ${WORKSPACE_FOLDER}:/workspaces/project
      - /path/to/data:/workspaces/project/00_Data:ro
```

```jsonc
// .devcontainer/devcontainer.json
{
  "name": "project",
  "dockerComposeFile": "docker-compose.yml",
  "service": "dev-archr",
  "workspaceFolder": "/workspaces/project",
  "customizations": { "vscode": { "extensions": [
    "reditorsupport.r", "ms-vscode-remote.remote-containers", "ms-python.python"
  ]}}
}
```

4. **Add a `settings.json`** to `.vscode` to configure R settings, radian, etc.:
   ```jsonc
   // .vscode/settings.json
   {
       "r.alwaysUseActiveTerminal": true,
       "r.bracketedPaste": true,
       "r.rterm.linux": "/opt/venvs/base/bin/radian",
       "r.sessionWatcher": true,
       "r.plot.useHttpgd": true,
       "r.lsp.diagnostics": false
   }
   ```

5. **Open the project in VS Code**. Then:
   - Press `Cmd+P` (Mac) or `Ctrl+P` (Windows/Linux) and type:
     ```
     >Dev Containers: Build & Attach
     ```
   - Select the **Dev Container** you just defined (service `base`).
   - To switch to ArchR: change `service` to `archr` and Reopen in Container.
   - Or run both services, then use “Dev Containers: Attach to Running Container” to attach a second VS Code window to `archr`.

7. **Long-running tasks**: For tasks that must continue even if you disconnect (like heavy alignment or processing), you can use `tmux` inside the container:
   ```bash
   tmux new -s your_session_name
   # run your alignment/analysis steps
   # detach with Ctrl+B, then D
   ```
   This ensures your process remains running if your SSH or VS Code session is interrupted.

---

## Operating the images and containers

### Manual runs (without VS Code)

Base image:

```bash
docker run --rm -it \
  -u $(id -u):$(id -g) \
  -v /path/to/project:/workspaces/project \
  -v /path/to/data:/workspaces/project/00_Data:ro \
  --memory=450g --cpus=50 \
  scdock-r-dev:v0.4.0 bash
```

ArchR image:

```bash
docker run --rm -it \
  -u $(id -u):$(id -g) \
  -v /path/to/project:/workspaces/project \
  -v /path/to/data:/workspaces/project/00_Data:ro \
  --memory=450g --cpus=50 \
  scdock-r-archr:v0.4.0 bash
```

Inside the container:

```bash
# regular R (radian)
r-base    # wrapper to start radian without ArchR toggle

# ArchR session (uses separate ArchR library path)
r-archr   # wrapper sets USE_ARCHR=1 then starts radian

# Manual toggle (if wrappers unavailable)
export USE_ARCHR=1
radian

Note: ArchR and its extras are installed into a separate R library path (ARCHR_LIB, default ~/R/archr-lib). With USE_ARCHR=1, this path is prepended to .libPaths(), so ArchR’s versions take precedence without altering the base R library.
```

If wrappers are not in PATH, toggle manually:

```bash
USE_ARCHR=1 radian
```

### Compose runs (without VS Code)

From your project directory containing `.devcontainer/docker-compose.yml`:

```bash
export LOCAL_UID=$(id -u); export LOCAL_GID=$(id -g); export WORKSPACE_FOLDER=$PWD
docker compose -f .devcontainer/docker-compose.yml up -d base
# or
docker compose -f .devcontainer/docker-compose.yml up -d archr
```

Enter a shell:

```bash
docker compose -f .devcontainer/docker-compose.yml exec base bash
docker compose -f .devcontainer/docker-compose.yml exec archr bash
```

Stop:

```bash
docker compose -f .devcontainer/docker-compose.yml down
```

### Nextflow/nf-core style runtime

```groovy
profiles {
  docker {
    docker.enabled = true
    docker.runOptions   = "-u ${System.env.NXF_UID ?: '1000'}:${System.env.NXF_GID ?: '1000'} -v /data2/nxf_tmp:/tmp"
    docker.fixOwnership = true
  }
}
```

### Tips
- VS Code: exclude heavy I/O paths in `.vscode/settings.json` to speed up indexing.
- Reproducibility: commit `renv.lock` and enable the `COPY renv.lock` line in the Dockerfile.
- MACS: MACS3 comes from Python; the ArchR image symlinks `macs2` if only `macs3` exists.
- Cairo: R Cairo package optional; ArchR works without it (plots vectorized).
- GSL: libgsl is installed in the base image.

### Engineering notes (Dev Container robustness)
- `.Rprofile` is safe and interactive-only; it checks for `~/.vscode-R/init.R` and enables httpgd only inside VS Code. This prevents crashes during non-interactive installs.
- A site-wide CRAN mirror is set at `/etc/R/Rprofile.site`, so scripted installs never see `@CRAN@`.
- When invoking R in automated steps, the build uses `R --vanilla` and explicit repos to avoid reading user/site profiles unexpectedly.
- For any postCreate scripts you add later, prefer:
  ```bash
  R_PROFILE_USER=/dev/null R_ENVIRON_USER=/dev/null Rscript --vanilla your_script.R
  ```

---

## R Library Architecture: Two-Tier Design

### Understanding the System Library vs User Library

The container uses a **two-tier R library architecture** for reproducibility, shareability, and efficient runtime package installation:

```
System Library (read-only)          User Library (writable)
/usr/local/lib/R/library            ~/R/x86_64-pc-linux-gnu-library/4.5
├─ Core packages (~80)              ├─ Runtime installs
├─ Pinned via renv.lock            ├─ Project-specific packages
├─ Same across all containers      ├─ Can differ per user/project
├─ Built into image (~10GB)        ├─ Persists in home directory
└─ Owned by root (read-only)       └─ Owned by devuser (writable)
```

**Check your library paths in R:**
```r
.libPaths()
# [1] "/home/devuser/R/x86_64-pc-linux-gnu-library/4.5"  # User library (writable)
# [2] "/usr/local/lib/R/library"                          # System library (read-only)
```

### Why Read-Only System Library?

**1. Reproducibility**
- Core packages are pinned via `renv.lock` at build time
- Same package versions across all container instances
- Updates to system packages are tested before baking into new image versions

**2. Image Shareability**
- Generic image with `devuser:1000` can be pushed to registry and shared with team
- System packages identical for all users
- Only user-installed packages differ (stored in home directory)

**3. Disk Space Efficiency**
- Core 80 packages (~10GB) shared across all container instances
- Each user's runtime installs (~1-5GB) stored separately in home directory
- Alternative (writable system lib) → every update creates new Docker layers

**4. Separation of Concerns**
- **System library**: Stable, tested baseline for all users
- **User library**: Experimental, project-specific, can be wiped without rebuilding image

### Installing R Packages at Runtime

**Runtime R installs automatically go to the writable user library** (`~/R/...`). No sudo needed:

```r
# Install from CRAN
install.packages("AnnotationHub")

# Install from Bioconductor
BiocManager::install("AnnotationHub")

# Install from GitHub
remotes::install_github("user/package")
```

Packages install to `/home/devuser/R/x86_64-pc-linux-gnu-library/4.5` and take **precedence** over system library versions.

### Expected Warning: "Installation paths not writeable"

**This warning is NORMAL and EXPECTED:**

```r
r$> BiocManager::install("AnnotationHub")
Bioconductor version 3.21 (BiocManager 1.30.26), R 4.5.0 (2025-04-11)
Installing package(s) 'AnnotationHub'
...
* DONE (AnnotationHub)

Installation paths not writeable, unable to update packages
  path: /usr/local/lib/R/library
  packages:
    aplot, BiocGenerics, boot, colorspace, Matrix, Seurat, ...
```

**What this means:**
- ✅ **AnnotationHub installed successfully** to user library
- ⚠️ BiocManager checked if any system packages need updates (they do)
- ❌ Cannot update system packages because you're not root (by design)

**This is intentional because:**
- Updating system packages would break reproducibility
- Updates should be tested before baking into new image versions
- You can install newer versions in user library (they take precedence)

**To suppress these warnings:**
```r
# Disable update checks (installs new packages only)
BiocManager::install("AnnotationHub", update = FALSE)
```

### Verifying Package Installation

**Check where a package is installed:**
```r
find.package("AnnotationHub")
# [1] "/home/devuser/R/x86_64-pc-linux-gnu-library/4.5/AnnotationHub"

# Check if package loads
library(AnnotationHub)
```

**List all packages and their locations:**
```r
ip <- installed.packages()[, c("Package", "LibPath")]
head(ip[ip[, "LibPath"] == .libPaths()[1], ])  # User library
head(ip[ip[, "LibPath"] == .libPaths()[2], ])  # System library
```

### ArchR Library Path (dev-archr service only)

When using the `dev-archr` service with `r-archr` wrapper:
```bash
r-archr  # Sets USE_ARCHR=1
```

The library path becomes:
```r
.libPaths()
# [1] "/home/devuser/R/archr-lib"                         # ArchR library (highest priority)
# [2] "/home/devuser/R/x86_64-pc-linux-gnu-library/4.5"  # User library
# [3] "/usr/local/lib/R/library"                          # System library
```

This allows ArchR's specific package versions to override system versions without conflicts.

### Creating Output Folders

R won't create intermediate directories for `saveRDS`. Ensure the path exists before writing:
```r
dir.create("03_Results/DEG", recursive = TRUE, showWarnings = FALSE)
saveRDS(dge, "03_Results/DEG/preprocessed_DGEList.rds")
```

If you bind-mount the workspace, ensure the host path is writable by your UID/GID (the container runs as your host IDs). Adjust permissions on the host if needed.

---

## Included Tools & Packages

1. **R** (v4.4.2) and key single-cell packages
   - CRAN: `tidyverse`, `ragg`, `devtools`, etc.
   - Bioconductor: `Seurat`, `SingleCellExperiment`, `DESeq2`, `ArchR`, etc.
   - GitHub: Additional packages like `Signac`, `DoubletFinder`, and more.

2. **Python** (v3.10) environment
   - Data science libraries: `numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`, ...
   - Single-cell: `scanpy`, `scvi-tools`, `squidpy`, `scvelo`, etc.
   - Epigenomics: `MACS3`, `episcanpy`, ...

3. **System Tools** 
   - Compiler toolchain, Java, TeX packages, etc. (for building R from source and generating PDFs).
   - Epigenomics-oriented CLI: **samtools**, **bcftools**, **bedtools**, and **scIBD**; heavy aligners and trimming/QC tools are intentionally excluded for a slim image.

4. **Quarto** for rendering R Markdown, Jupyter Notebooks, and other scientific documents.

See the Dockerfile at `docker/base/Dockerfile` and core installer at `docker/base/R/install_core.R` for full details on what gets installed.

---

## Troubleshooting

1. **Permission Errors**:  
   If you run into file permission issues when mounting a host directory, ensure the container user matches your host user `UID`/`GID`. You can do this by passing additional build/run arguments (e.g., `--build-arg USER_ID=1001 --build-arg GROUP_ID=1001`). Or, adjust folder permissions on the host.

2. **Out-of-Memory Errors**:  
   Large single-cell data analyses can require significant RAM. Assign more memory/CPU to Docker if you’re on Docker Desktop, or scale up your server if you’re remote.

3. **GitHub Rate Limits**:  
   If you see error messages about rate limits or missing Git credentials, try setting a GitHub PAT:  
   ```bash
   docker build . \
       -f ./.devcontainer/Dockerfile \
       --build-arg GITHUB_PAT=<your_github_pat> \
       -t scdock-r-dev:latest
   ```

4. **R Package Version Mismatch**:  
   Libraries compiled under a different R version might cause “object not found” or `.so` loading errors. Rebuild the container with consistent R version arguments or reinstall packages inside the container to match R’s version.

---

## Backlog and roadmap

This section reflects the current design (dev-core vs dev-archr) and what remains to validate or improve.

### Now
- Build v0.4.1 images and validate
  - dev-core: build and run `.devcontainer/scripts/poststart_sanity.sh`
  - dev-archr: build and validate ArchR (`library(ArchR)`) and macs2 symlink
- Validate httpgd in VS Code radian; ensure fallback works as documented
- Smoke-test Python venvs: base/squid/atac/comms core imports

### Next
- Decide scGLUE TensorFlow wheels strategy (pin `tensorflow-cpu` vs runtime optional)
- Document `scenicplus` runtime install and known dependency pins in COMMS env
- Add round-trip QA script for `.h5mu`/`.h5ad` conversions (embeddings/graphs/layers checks)
- Make `dev-archr` the recommended default in templates/examples
- Add a training example (RNA+ATAC) in `02_Analysis` showing R↔Python interop

### Later
- Add CI smoke tests for sanity script and env imports
- Consider slim (no toolchain) and GPU variants; document supported scopes
- Evaluate packaging some shared R functionality as internal packages

---

## Working with shared code inside the container

For shared code I used to mount its external Git repo inside the container, and tracked changes of shared code separately. This works, but there are more maintainable options. Below are three patterns, from most recommended to situational.

### 1) Git submodule (recommended)
Keep shared code as a submodule within your project tree. Benefits: pins an exact commit per project, single Git workflow from VS Code, no separate mounts.

Add the submodule:
```bash
cd /scratch/current/antonz/projects/13036-DM_DMlab_summer_2025
git submodule add git@github.com:tony-zhelonkin/R_GSEA_visualisations.git 01_Scripts/shared
```

Daily usage:
```bash
# Parent repo (your analysis project)
git status
# If the submodule commit changed, parent sees it as modified

# Work inside the submodule itself
cd 01_Scripts/shared
git checkout main
git pull                 # or edit, commit, push here
# After changes:
cd ../..
git add 01_Scripts/shared
git commit -m "Update shared scripts to latest"
git push
```

Cloning a project with submodules:
```bash
git clone --recurse-submodules <your-project-repo>
# or
git clone <your-project-repo>
cd <repo>
git submodule init && git submodule update
```

VS Code devcontainer: nothing special required; the submodule lives inside the workspace and is tracked by the parent repo. Use it like any folder in `01_Scripts/shared`.

### 2) R package for shared code
If the shared code becomes stable and broadly reused, convert it to an R package for cleaner dependency management and versioning.

Install in the container:
```r
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
devtools::install_github("tony-zhelonkin/R_GSEA_visualisations")
```
Pros: semantic versioning, DESCRIPTION-based deps, easy to pin via renv. Cons: more overhead to structure as a package.

### 3) External mount (works, but less ideal)
Mount an external shared repo into your container, and manage its Git lifecycle in another terminal outside the container.

Example devcontainer.json optional mount:
```jsonc
"mounts": [
  "source=/data1/users/antonz/pipeline/R_GSEA_visualisations,\
   target=/workspaces/project/01_Scripts/shared,\
   type=bind"
]
```
Caveats: two Git flows; easier to drift. Prefer submodules or packages when possible.

Recommendation: in general, submodules are fine; migrate to an R package if/when the shared code matures.

---

## Contributing

Contributions are welcome! To propose a change:

1. Fork this repository.
2. Create a new branch for your feature/fix.
3. Commit and push your changes.
4. Submit a Pull Request (PR) describing your modifications.

Please open an issue for any bugs, installation problems, or improvement ideas.

---

## License

This project is distributed under the [MIT License](https://github.com/tony-zhelonkin/scbio-docker/blob/main/LICENCE.md). Feel free to use, modify, and distribute it as permitted.

---

## Acknowledgments

- **Maintainer**: [Anton Zhelonkin](mailto:anton.bioinf.md@gmail.com)
- **Huge thanks** to [Rami Krispin’s vscode-r repo](https://github.com/RamiKrispin/vscode-r) for serving as a fantastic inspiration.
- Thanks to all authors of open-source bioinformatics tools included here.

---

## Beginner-friendly vignette: My standard project layout and workflows

### Typical analysis repository layout

```text
├── .devcontainer/devcontainer.json     # Container setup
├── .vscode/                            # VS Code settings
├── 00_Data/                            # Raw data and reference files
├── 01_Scripts/                         # Analysis scripts and custom functions
│   ├── module/                         # External git submodules for shared code
│   ├── Py_scripts/                     # Custom project-specific Python scripts
│   └── R_scripts/                      # Custom project-specific R scripts
├── 02_Analysis/                        # Main analysis pipeline
├── 03_Results/                         # Output files and results
│   ├── 01_Preprocessing/               # Data preprocessing results
│   └── 02_Analysis/                    # Analysis results and plots
└── README.md                           # Project-level README
```

### Dev Container workflows

You can run your environment either as a single-image devcontainer (simple) or via Docker Compose (two services: base and archr). Both approaches bind your project into the container and let you use VS Code Remote - Containers.

#### A) Single-image devcontainer.json (simple)

Example:
```jsonc
{
  "name": "your-project",
  "image": "scdock-r-dev:v0.4.1",
  "workspaceMount": "source=/abs/path/to/your/project,target=/workspaces/project,type=bind",
  "workspaceFolder": "/workspaces/project",
  "runArgs": ["--memory=128g","--cpus=60"],
  "customizations": {
    "vscode": { "extensions": [
      "rdebugger.r-debugger","reditorsupport.r","quarto.quarto",
      "ms-azuretools.vscode-docker","ms-vscode-remote.remote-containers",
      "ms-python.python"
    ]}
  },
  "postStartCommand": "bash -lc 'chmod +x .devcontainer/scripts/poststart_sanity.sh && .devcontainer/scripts/poststart_sanity.sh'"
}
```
- VS Code: Dev Containers: Reopen in Container
- R: `r-base` for radian (base lib), `r-archr` for ArchR toggle
- Python venv switching: `usepy base|squid|atac|comms`

#### B) Compose-based devcontainer.json (multi-variant)

`.devcontainer/docker-compose.yml` defines `dev-core` and `dev-archr` services. Use this `devcontainer.json` (recommended default is `dev-archr`):
```jsonc
{
  "name": "your-project",
  "dockerComposeFile": ".devcontainer/docker-compose.yml",
  "service": "dev-archr",      # or "dev-core"
  "workspaceFolder": "/workspaces/project",
  "customizations": { "vscode": { "extensions": [
    "reditorsupport.r","ms-python.python","ms-vscode-remote.remote-containers"
  ]}}
}
```
- Switch services by changing `service` and Reopen in Container, or attach to the running container.
- Same commands inside as in the single-image flow (radian wrappers, `usepy` switchers).

Important: r-base/r-archr vs service choice

- Service `dev-core` (image: scdock-r-dev):
  - `r-base` works (base R libs).
  - `r-archr` only sets `USE_ARCHR=1`; ArchR is NOT installed here → `library(ArchR)` fails. Use the `dev-archr` service when you need ArchR.
- Service `dev-archr` (image: scdock-r-archr):
  - `r-base` works and does NOT prepend ArchR lib path.
  - `r-archr` sets `USE_ARCHR=1` so `.Rprofile` prepends `ARCHR_LIB`; `library(ArchR)` succeeds.

Recommendation: default to `dev-archr` in development unless you specifically want the leaner `dev-core` image.

### Switching environments quickly

- R (radian front-end):
  - Base: `r-base`
  - ArchR: `r-archr` (sets `USE_ARCHR=1`, prepends `ARCHR_LIB` in `.Rprofile`)
- Python venvs:
  - `usepy base|squid|atac|comms`
  - One-off: `py-base`, `py-squid`, `py-atac`, `py-comms`

### Worked example: RNA + ATAC interoperability and integration

Scenario: You have parallel scRNA-seq and scATAC-seq datasets preprocessed/UMAPed and annotated in Seurat. You want to integrate in Python (scVI/PeakVI or SCGLUE) and run GRN inference (SCENIC+), then round-trip back to Seurat.

1) Start in R for Seurat preprocessing
```bash
r-base   # radian (base libs)
```
- Ensure assays and reductions are in good shape (Seurat v5: counts/data/scale.data).
- Check reductions and graphs:
```r
names(s@reductions); names(s@graphs); Layers(s[["RNA"]])
```

2) Export to Python-friendly formats
- Single-modality RNA or ATAC: write `.h5ad` (faithful with loadings/graphs):
```r
library(MuDataSeurat)
WriteH5AD(object = s, file = "rna.h5ad", assay = "RNA")
```
- True multiome (paired) or to keep modalities together: write `.h5mu`:
```r
WriteH5MU(s, "multiome.h5mu")
```
- Alternative R-only path: `anndataR` reads/writes H5AD natively:
```r
library(anndataR)
s <- read_h5ad("in.h5ad", as = "Seurat")
# or write_h5ad(s, "out.h5ad")
```

3) Switch to Python (base) for scVI/SCGLUE
```bash
usepy base
python -c "import scanpy as sc, muon as mu; print('OK')"
```
- Load data:
```python
import scanpy as sc, muon as mu
m = mu.read_h5mu("multiome.h5mu")  # if multiome
rna = m.mod["RNA"]; atac = m.mod["ATAC"]
# or: rna = sc.read_h5ad("rna.h5ad"); atac = sc.read_h5ad("atac.h5ad")
```
- Compute embeddings (examples):
  - scVI on RNA → store in `obsm['X_scvi']`
  - PeakVI on ATAC → `obsm['X_peakvi']`
  - SCGLUE → `obsm['X_scglue']`
- Ensure standard keys so converters recognize them.

4) Optional: use `usepy atac` if you keep ATAC-only stack separate
```bash
usepy atac   # snapatac2 and ATAC-focused deps
```

5) Write back to disk
- Keep multiome fidelity: `.h5mu`
- Or per-modality `.h5ad` files with embeddings in `obsm` and graphs in `obsp`.

6) Return to R and import
```bash
r-base
```
- Via MuDataSeurat (recommended for fidelity):
```r
library(MuDataSeurat)
# Single assay back from H5AD
s <- Seurat::ReadH5AD("rna.h5ad")
# Multiome
s_multi <- ReadH5MU("multiome.h5mu")
```
- Check round-trip:
```r
names(s@reductions)    # expect pca, umap, tsne, scvi, peakvi, scglue ...
head(Embeddings(s, "pca"))
names(s@graphs)
```

7) Conversion checklists (don’t lose artifacts)
- Before writing H5AD (Python):
  - `adata.X` vs `adata.raw.X` consistent with expectations
  - Embeddings in `obsm`: `X_pca`, `X_umap`, `X_scvi`, `X_peakvi`, `X_scglue`
  - Graphs in `obsp`: connectivities/distances
- After import in R:
  - `names(s@reductions)` contains your embeddings
  - `names(s@graphs)` populated
  - For Seurat v5, `Layers(s[["RNA"]])` includes counts/data/scale

8) GRN inference (SCENIC+)
- Use the COMMS venv for LR/GRN tools:
```bash
usepy comms
python -c "import pycistarget, pycistopic; print('OK')"  # scenicplus via source if needed
```
- If `scenicplus` is required and not preinstalled, install at runtime:
```bash
usepy comms
pip install 'scenicplus @ git+https://github.com/aertslab/SCENICplus.git'
```

### Tips and pitfalls
- Prefer `.h5mu` for true multiome; use `.h5ad` per-assay when needed.
- Use predictable embedding keys (X_scvi, X_peakvi, X_scglue) to preserve reductions on round-trip.
- Large datasets: prefer backed I/O (.h5mu with muon; Seurat v5 disk-aware flows).
- If a reduction/graph is missing after import, re-emit using `MuDataSeurat::WriteH5AD`.

---

## Permissions, UID/GID, and AD users (group names with spaces)

- Containers run as your host UID/GID (passed at build and used at runtime). This avoids root-owned files on your mounts.
- For AD-controlled users (e.g., primary group name contains a space like `domain users`):
  - The Dockerfile already normalizes/handles group creation via numeric GID and falls back to safe names if needed.
  - When binding a workspace path, ensure the host directory is writable by your UID/GID:
    ```bash
    id -u; id -g; id -gn
    ls -ld /path/to/workspace
    # Adjust on host if needed:
    sudo chown -R $(id -u):$(id -g) /path/to/workspace
    # or grant group write:
    sudo chmod -R g+rwX /path/to/workspace
    ```
- R user library: the container sets `R_LIBS_USER` and `.Rprofile` prepends it. Startup sanity now checks:
  - R user lib writable
  - Workspace write (creates and removes a temp dir)
  If either fails, the sanity script prints IDs and path info to guide fixes.

- VS Code Dev Containers: When using compose, both services inherit the same mount and user IDs; choose `dev-archr` by default if you want ArchR availability.

---

## Image slimming (what we removed and how to keep it small)

- Heavy R annotation/data packages moved to optional (install on-demand):
  - BSgenome.* (hg19, hg38, mm10, mm39), EnsDb.* (multiple versions), org.*.eg.db, reactome.db
  - Reason: multi-GB footprint; not always needed for development
  - To install when required:
    ```bash
    # Build with heavy data included
    docker build -f .devcontainer/Dockerfile \
      --build-arg INCLUDE_HEAVY_R_DATA=1 \
      -t scdock-r-dev:heavy .
    ```
    Or inside R at runtime via BiocManager installs.

- Python venv trims:
  - Moved cross R/Py bridging packages to on-demand: pyreadr, rpy2, anndata2ri (commented in base_requirements.txt)
  - Keep heavy stacks only where necessary; avoid duplicating scanpy/scvi across multiple venvs

- Cleanup and build hygiene:
  - Remove build sources after install (R sources are removed; keep iterating where applicable)
  - Consider TinyTeX or minimal TeX set if PDF output isn’t central
  - Optionally purge renv cache in final stage if you don’t need rebuild acceleration

- Quick size audit (inside a container):
  ```bash
  du -hsx /* 2>/dev/null | sort -h | tail -n 20
  du -hs /usr/local/lib/R/library/* | sort -h | tail -n 30
  du -hs /opt/venvs/* | sort -h
  ```

---
