# DevOps Guide - scbio-docker v0.5.1+

Quick reference for building, running, and operating the single-cell Docker environment.

---

## Build Instructions

### Recommended: Multi-Stage Build (v0.5.1+)

**Use the optimized build script for best size efficiency:**

```bash
# Simple build with defaults
./build-optimized.sh

# With GitHub PAT (recommended to avoid rate limits)
./build-optimized.sh --github-pat ghp_your_token_here

# Custom tag
./build-optimized.sh --tag scdock-r-dev:v0.5.1
```

**Benefits of multi-stage build:**
- ✅ **True ~20GB final image** (no layer bloat)
- ✅ Discards all build artifacts and caches
- ✅ Preserves build-essential and dev libraries for runtime package installation
- ✅ Same functionality as v0.5.0

**Manual multi-stage build:**

```bash
export GITHUB_PAT=ghp_your_token_here

docker build . \
  -f .devcontainer/Dockerfile.optimized \
  --build-arg GITHUB_PAT=$GITHUB_PAT \
  --build-arg USER_ID=$(id -u) \
  --build-arg GROUP_ID=$(id -g) \
  --build-arg USER=$USER \
  --build-arg GROUP=$(id -gn) \
  -t scdock-r-dev:v0.5.1

unset GITHUB_PAT
```

### Legacy: Single-Stage Build (v0.5.0)

**Note:** This build method results in 533GB reported size due to Docker layer accounting. Use multi-stage build instead.

```bash
docker build . \
  -f .devcontainer/Dockerfile \
  --build-arg GITHUB_PAT=$GITHUB_PAT \
  --build-arg USER_ID=$(id -u) \
  --build-arg GROUP_ID=$(id -g) \
  --build-arg USER=$USER \
  --build-arg GROUP=$(id -gn) \
  -t scdock-r-dev:v0.5.0
```

### Extract renv.lock (Do Once After First Build)

```bash
# Extract lockfile and manifest
CID=$(docker create scdock-r-dev:v0.5.1)
docker cp $CID:/opt/settings/renv.lock ./renv.lock
docker cp $CID:/opt/settings/R-packages-manifest.csv ./R-packages-manifest.csv
docker rm $CID

# Commit for deterministic builds
git add renv.lock R-packages-manifest.csv
git commit -m "Pin R packages via renv.lock (v0.5.1)"
```

### Pull Official ArchR Image (Recommended)

```bash
# Pull the official ArchR Docker image
docker pull greenleaflab/archr:1.0.3-base-r4.4

# Tag for consistency with compose files
docker tag greenleaflab/archr:1.0.3-base-r4.4 archr-official:1.0.3
```

**Why use the official image?**
- Maintained by ArchR developers
- Stable R 4.4 (tested with ArchR)
- Smaller size (~15GB vs custom build)
- No version conflicts with dev-core (R 4.5)

---

## Running Containers

### Option 1: Docker Run (Manual)

**Base R image:**
```bash
docker run --rm -it \
  -u $(id -u):$(id -g) \
  -v /path/to/project:/workspaces/project \
  -v /path/to/data:/workspaces/project/00_Data:ro \
  --memory=450g --cpus=50 \
  scdock-r-dev:v0.5.0 bash
```

**ArchR image (official):**
```bash
docker run --rm -it \
  -u $(id -u):$(id -g) \
  -v /path/to/project:/workspaces/project \
  --memory=450g --cpus=50 \
  greenleaflab/archr:1.0.3-base-r4.4 bash
```

### Option 2: Docker Compose (Multi-Service)

**Setup project with compose:**
```bash
cd /path/to/your/project

# Create .devcontainer directory
mkdir -p .devcontainer

# Copy compose file from repo
cp /path/to/scbio-docker/.devcontainer/docker-compose.yml .devcontainer/

# Edit to add your data mounts
vim .devcontainer/docker-compose.yml
```

**Start services:**
```bash
export LOCAL_UID=$(id -u)
export LOCAL_GID=$(id -g)
export WORKSPACE_FOLDER=$PWD

# Start base R service
docker compose -f .devcontainer/docker-compose.yml up -d dev-core

# Or start ArchR service
docker compose -f .devcontainer/docker-compose.yml up -d dev-archr

# Enter container
docker compose -f .devcontainer/docker-compose.yml exec dev-core bash

# Stop services
docker compose -f .devcontainer/docker-compose.yml down
```

### Option 3: VS Code Dev Containers (Recommended)

**Quick start with init-project.sh:**
```bash
# Initialize new project from template
./init-project.sh ~/projects/my-analysis basic-rna

# Open in VS Code
code ~/projects/my-analysis

# In VS Code: Cmd/Ctrl+Shift+P → "Dev Containers: Reopen in Container"
```

**Manual setup from existing project:**
1. Copy `.devcontainer/` folder to your project root
2. Edit `docker-compose.yml` to add data mounts
3. Open project in VS Code
4. `Cmd/Ctrl+Shift+P` → "Dev Containers: Reopen in Container"

**Switch between dev-core and dev-archr:**
- Edit `.devcontainer/devcontainer.json`
- Change `"service": "dev-core"` to `"service": "dev-archr"`
- Reopen in container

See **"Starting from Project Directory"** section below for full workflow.

---

## Environment Switching

### Python Virtual Environments

**Base venv (pre-installed):**
```bash
# Default environment (already active)
which python  # /opt/venvs/base/bin/python
```

**Layered venvs (created on-demand):**
```bash
# Create and activate squidpy venv (inherits base packages)
usepy squid

# Create and activate ATAC venv
usepy atac

# Create and activate communication venv
usepy comms

# Switch back to base
usepy base
```

**Manual venv creation:**
```bash
# Create specific venv
create_layered_venv.sh squid squid_requirements.txt

# Activate
source /opt/venvs/squid/bin/activate
```

**Check active environment:**
```bash
which python && python -V
pip list | head
```

### R Sessions

**Standard R (radian terminal):**
```bash
# Launch radian
radian

# Or via wrapper
r-base
```

**ArchR-enabled R (using official ArchR image):**
```bash
# Switch to dev-archr service (see "Switching to ArchR" section below)
# Then start R normally
radian

# ArchR is available directly
library(ArchR)
```

**Persistent R sessions with tmux:**
```bash
# Start persistent session
tmux new-session -s my-analysis radian

# Detach: Ctrl+B, then D
# Session persists even if SSH disconnects

# Re-attach later
tmux attach -t my-analysis

# List active sessions
tmux ls

# Kill session
tmux kill-session -t my-analysis
```

---

## VS Code Integration

### R Settings (`.vscode/settings.json`)

```jsonc
{
  "r.rterm.linux": "/opt/venvs/base/bin/radian",
  "r.alwaysUseActiveTerminal": true,
  "r.bracketedPaste": true,
  "r.sessionWatcher": true,
  "r.plot.useHttpgd": true,
  "r.lsp.diagnostics": false
}
```

### Workflow

1. Start container (VS Code Dev Container or manual)
2. Open VS Code terminal
3. Start radian in tmux: `tmux new-session -s analysis radian`
4. Write R code in `.R` files
5. Send to radian terminal: `Cmd/Ctrl+Enter`
6. Plots appear via httpgd in VS Code

---

## Starting from Project Directory

### Workflow 1: Initialize New Project

**Step 1: Create project from template**
```bash
cd /path/to/scbio-docker

# Choose template: basic-rna, multimodal, archr-focused, example-DMATAC
./init-project.sh ~/projects/my-scrna-analysis basic-rna
```

**Step 2: Configure data mounts**
```bash
cd ~/projects/my-scrna-analysis

# Edit docker-compose.yml to add your data mounts
nano .devcontainer/docker-compose.yml
```

Example data mount:
```yaml
services:
  dev-core:
    volumes:
      - ${WORKSPACE_FOLDER:-.}:/workspaces/project
      - /scratch/sequencing/dataset1:/workspaces/project/data/raw:ro
```

**Step 3: Open in VS Code**
```bash
code ~/projects/my-scrna-analysis
```

**Step 4: Reopen in container**
- `Cmd/Ctrl+Shift+P` → "Dev Containers: Reopen in Container"
- Wait for container to start (check post-start output)

**Step 5: Start analysis**
```bash
# Terminal opens inside container
whoami  # Should show your username
pwd     # /workspaces/project

# Start R session
tmux new-session -s analysis radian

# Or Python
python
```

### Workflow 2: Existing Project Setup

**If you already have a project directory:**

```bash
cd /path/to/existing-project

# Copy devcontainer setup from scbio-docker repo
cp -r /path/to/scbio-docker/.devcontainer .

# Copy VS Code settings
mkdir -p .vscode
cp /path/to/scbio-docker/templates/.vscode/settings.json .vscode/

# Edit docker-compose.yml to add data mounts
nano .devcontainer/docker-compose.yml

# Open in VS Code
code .

# Reopen in container (Cmd/Ctrl+Shift+P)
```

### Switching to ArchR

**When you need ArchR for scATAC-seq analysis:**

**Method 1: Change devcontainer.json (persistent)**
```bash
# Edit .devcontainer/devcontainer.json
nano .devcontainer/devcontainer.json

# Change line:
"service": "dev-core"  →  "service": "dev-archr"

# In VS Code: Cmd/Ctrl+Shift+P → "Dev Containers: Rebuild and Reopen in Container"
```

**Method 2: Attach to running container (temporary)**
```bash
# Start dev-archr service manually
cd /path/to/project
docker compose -f .devcontainer/docker-compose.yml up -d dev-archr

# In VS Code:
# Cmd/Ctrl+Shift+P → "Dev Containers: Attach to Running Container"
# Select container with "dev-archr" in name
```

**Method 3: Manual docker run**
```bash
docker run --rm -it \
  -u $(id -u):$(id -g) \
  -v /path/to/project:/workspaces/project \
  -v /path/to/data:/workspaces/project/data/raw:ro \
  greenleaflab/archr:1.0.3-base-r4.4 bash
```

**Verify ArchR environment:**
```bash
# Inside container
R --version  # Should show R 4.4.x

radian
# In R:
library(ArchR)
packageVersion("ArchR")  # Should be 1.0.3
```

---

## Package Management

### R Package Installation

**Two-tier library design:**
- **System library** (`/usr/local/lib/R/library`): Read-only, ~80 core packages, pinned via renv
- **User library** (`~/R/...`): Writable, runtime installs, project-specific packages

**Runtime installation (to user library, no sudo needed):**
```r
# Automatically installs to ~/R/...
if (!require("GSVA")) BiocManager::install("GSVA")
install.packages("ggExtra")

# Check library paths
.libPaths()
# [1] "/home/devuser/R/x86_64-pc-linux-gnu-library/4.5"  # writable (runtime installs)
# [2] "/usr/local/lib/R/library"                          # system (core packages)
```

**Expected warning (NORMAL and harmless):**
```r
BiocManager::install("AnnotationHub")
...
* DONE (AnnotationHub)

Installation paths not writeable, unable to update packages
  path: /usr/local/lib/R/library
  packages:
    aplot, BiocGenerics, Matrix, Seurat, ...
```

**What this means:**
- ✅ Package installed successfully to user library
- ⚠️ BiocManager checked for updates to system packages (default behavior)
- ❌ Cannot update system packages (by design for reproducibility)

**To suppress warnings:**
```r
BiocManager::install("PACKAGE", update = FALSE)
```

**Pre-installed core packages (~80):**
- Seurat ecosystem: Seurat, Signac, BPCells, sctransform
- RNA-seq: edgeR, limma, DESeq2, scran
- GSEA: clusterProfiler, GSVA, fgsea, msigdbr
- Multi-factorial: muscat, MOFA2, mixOmics, liger
- Visualization: ggplot2, patchwork, ComplexHeatmap
- Interop: anndataR, MuDataSeurat, reticulate

See `.devcontainer/install_R_core.R` for full list.

**Project reproducibility with renv:**
```r
# Initialize renv for project
renv::init()

# After installing packages
renv::snapshot()  # Creates project renv.lock

# Restore project packages
renv::restore()
```

### Python Package Installation

**Install to active venv:**
```bash
# In base venv
pip install package-name

# In layered venv
usepy squid
pip install additional-package

# Freeze for reproducibility
pip freeze > requirements-frozen.txt
```

**Install from requirements:**
```bash
pip install -r requirements-extra.txt
```

---

## Troubleshooting

### Permission Issues

**Check container user:**
```bash
id -u && id -g && id -gn
ls -ld /workspaces/project
```

**Fix on host:**
```bash
sudo chown -R $(id -u):$(id -g) /path/to/project
```

### R Library Not Writable

**Check:**
```r
.libPaths()
file.access(.libPaths()[1], 2) == 0  # Should be TRUE
```

**Fix:** User library should be created automatically at `~/R/...`

### Python Venv Issues

**Recreate layered venv:**
```bash
rm -rf /opt/venvs/squid
usepy squid  # Recreates
```

### Sanity Check

**Run post-start checks:**
```bash
.devcontainer/scripts/poststart_sanity.sh
```

**Expected output:**
```
==== Devcontainer sanity checks ====
Python: OK
R:      OK
httpgd: OK
scanpy import: OK
default venv: /opt/venvs/base OK
radian: OK
...
```

---

## Multi-User Setup (Shared Machine)

**Each user builds with their UID/GID:**
```bash
# User A (UID 1001)
docker build ... \
  --build-arg USER_ID=1001 \
  --build-arg GROUP_ID=1001 \
  --build-arg USER=usera \
  -t scdock-r-dev:v0.5.0-usera

# User B (UID 1002)
docker build ... \
  --build-arg USER_ID=1002 \
  --build-arg GROUP_ID=1002 \
  --build-arg USER=userb \
  -t scdock-r-dev:v0.5.0-userb
```

**Or use compose with LOCAL_UID:**
```bash
# User A
LOCAL_UID=1001 LOCAL_GID=1001 docker compose up -d

# User B
LOCAL_UID=1002 LOCAL_GID=1002 docker compose up -d
```

---

## Common Operations

### Check Image Size

```bash
docker images scdock-r-dev:v0.5.0
```

### Size Audit Inside Container

```bash
# Top-level directories
du -hsx /* 2>/dev/null | sort -h | tail -n 20

# R libraries
du -hs /usr/local/lib/R/library/* | sort -h | tail -n 30

# Python venvs
du -hs /opt/venvs/* | sort -h
```

### Clean Up Old Images

```bash
# Remove dangling images
docker image prune

# Remove specific old version
docker rmi scdock-r-dev:v0.4.1

# Remove all stopped containers
docker container prune
```

### Export/Import Images

```bash
# Export
docker save scdock-r-dev:v0.5.0 | gzip > scdock-r-dev-v0.5.0.tar.gz

# Import
docker load < scdock-r-dev-v0.5.0.tar.gz
```

---

## Quick Reference

| Task | Command |
|------|---------|
| Build base image | `docker build -f .devcontainer/Dockerfile -t scdock-r-dev:v0.5.0 .` |
| Pull ArchR image | `docker pull greenleaflab/archr:1.0.3-base-r4.4` |
| Initialize project | `./init-project.sh ~/projects/my-analysis basic-rna` |
| Run container | `docker run -it -u $(id -u):$(id -g) -v $PWD:/workspaces/project scdock-r-dev:v0.5.0` |
| Start radian | `radian` or `r-base` |
| Start radian in tmux | `tmux new-session -s analysis radian` |
| Switch Python venv | `usepy squid\|atac\|comms\|base` |
| Create layered venv | `create_layered_venv.sh <name> <requirements.txt>` |
| Install R package | `BiocManager::install("PACKAGE")` |
| Install Python package | `pip install package-name` |
| Snapshot R packages | `renv::snapshot()` |
| Freeze Python packages | `pip freeze > requirements-frozen.txt` |
| Sanity check | `.devcontainer/scripts/poststart_sanity.sh` |

---

## Project Templates

### Typical Project Structure

```
my-project/
├── .devcontainer/
│   ├── devcontainer.json
│   └── docker-compose.yml
├── .vscode/
│   └── settings.json
├── renv.lock                    # R package snapshot
├── requirements-frozen.txt       # Python package snapshot
├── .environments/
│   └── squidpy_requirements.txt  # Extra Python tools
├── 00_Data/                     # Raw data (read-only mount)
├── 01_Scripts/
│   ├── R/
│   └── Py/
├── 02_Analysis/                 # Main analysis scripts
└── 03_Results/                  # Output files
```

### devcontainer.json Template

```jsonc
{
  "name": "My Project",
  "dockerComposeFile": "docker-compose.yml",
  "service": "dev-core",
  "workspaceFolder": "/workspaces/project",
  "customizations": {
    "vscode": {
      "extensions": [
        "reditorsupport.r",
        "ms-python.python",
        "quarto.quarto"
      ],
      "settings": {
        "r.rterm.linux": "/opt/venvs/base/bin/radian",
        "r.alwaysUseActiveTerminal": true,
        "r.plot.useHttpgd": true
      }
    }
  },
  "postCreateCommand": "bash .devcontainer/post-create.sh"
}
```

### docker-compose.yml Template

```yaml
version: "3.8"
services:
  dev-core:
    image: scdock-r-dev:v0.5.0
    user: "${LOCAL_UID:-1000}:${LOCAL_GID:-1000}"
    working_dir: /workspaces/project
    volumes:
      - ${WORKSPACE_FOLDER:-.}:/workspaces/project
      - /path/to/data:/workspaces/project/00_Data:ro
    stdin_open: true
    tty: true
```

---

## Known Issues

### Image Size Discrepancy (v0.5.0 - FIXED in v0.5.1)

**v0.5.0 Issue:**
- Docker reported 533GB size due to layer accounting
- Actual filesystem was only ~20GB
- Caused by Docker BuildKit storing intermediate layer states

**Solution (v0.5.1+):**
Use the **multi-stage build** (`Dockerfile.optimized`):
```bash
./build-optimized.sh
```

**Result:** True ~20GB image size with no layer bloat.

**For v0.5.0 users:** Migrate to v0.5.1 multi-stage build, or use export/import workaround:
```bash
docker export $(docker create scdock-r-dev:v0.5.0) | docker import - scdock-r-dev:v0.5.0-flat
```

---

## Version History

- **v0.5.1** (Current): Multi-stage build, **true ~20GB image**, preserves runtime package installation
- **v0.5.0**: Size-optimized attempt, R 4.5 + Bioc 3.21, layered venvs (533GB Docker size issue)
- **v0.4.1**: R 4.4.2 + Bioc 3.20, 4 full Python venvs (500GB)
- **v0.4.0**: Full TeX, no cache cleanup

## Runtime Package Installation

See **[RUNTIME_INSTALL.md](RUNTIME_INSTALL.md)** for comprehensive guide on:
- Installing R packages (CRAN, Bioconductor, GitHub)
- Installing Python packages to venvs
- Installing system dependencies with `sudo apt-get`
- Handling compilation failures
- Best practices for reproducibility

---

## Support

- **Documentation**: See `README.md` and `CLAUDE.md`
- **Issues**: Check `.devcontainer/scripts/poststart_sanity.sh` output
- **Logs**: `docker logs <container_id>`
- **GitHub**: [Repository Issues](https://github.com/tony-zhelonkin/scbio-docker/issues)
