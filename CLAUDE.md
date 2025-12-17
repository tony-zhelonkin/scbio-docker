# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a Docker-based development environment for single-cell RNA-seq and epigenomics analyses, integrated with VS Code Remote Containers. The repository provides optimized Docker images with a focus on reproducibility and size efficiency.

**Current Version:** v0.5.3 (AI Tools Ready)

**Image Variants:**
- **scdock-r-dev:v0.5.3** (base): R 4.5 + Bioc 3.21 + Python 3.10 with core bioinformatics packages (**true ~20GB image**)
- **greenleaflab/archr:1.0.3-base-r4.4** (official ArchR): R 4.4 + ArchR 1.0.3, maintained by ArchR developers

**Key Changes in v0.5.3:**
- **AI Tools Ready**: Pre-installed dependencies for SciAgent-toolkit MCP servers
  - Node.js 20 LTS (npm, npx) for Sequential Thinking, Context7 MCP servers
  - uv/uvx for ToolUniverse, PAL, Serena MCP servers
  - `toml` Python package for Codex CLI configuration
- No more permission errors or installation hangs during `setup-ai.sh`
- Maintains ~20GB target size (Node.js and uv add ~100MB)

**v0.5.2 Changes (carried forward):**
- **Additional R packages pre-installed**: chromVAR, motifmatchr, TFBSTools, JASPAR2022, SingleR, celldex, AnnotationHub, EnsDb.Mmusculus.v79, crescendo (GitHub)
- **Ubuntu libraries added**: libhdf5-dev, libgsl-dev (enables hdf5r, DirichletMultinomial compilation)
- **Bug fix**: safe_install() now handles meta-packages correctly (tidyverse installs properly)

**v0.5.1 Changes (carried forward):**
- **Multi-stage build** - completely discards build artifacts, no layer bloat
- **True ~20GB Docker image** (not just filesystem, actual reported size)
- **Build tools preserved** - can compile R/Python packages at runtime

**v0.5.0 Optimizations (carried forward):**
- Aggressive cache cleanup (renv, pip, build artifacts)
- TinyTeX instead of full TeX distribution
- Single Python base venv with layered venvs for specialized tools
- Core R packages pre-installed (~80), additional packages installed at runtime
- Official ArchR image instead of custom build
- Project templates + init-project.sh for quick scaffolding

Build: prefer `scripts/build.sh` or `docker build -f docker/base/Dockerfile`

## Build Commands

### Building the base image (v0.5.3 - Multi-Stage)

**IMPORTANT: Build Strategy (Shareable vs Personal)**

The build system supports TWO modes:
- **GENERIC (default)**: Creates `devuser:1000` - shareable with team/registry
- **PERSONAL**: Creates image with YOUR UID - only for your use

**Recommended: Generic build (shareable)**
```bash
scripts/build.sh                       # Generic build (devuser:1000)
scripts/build.sh --github-pat ghp_...  # With GitHub PAT
```

**Personal build (your UID only):**
```bash
scripts/build.sh --personal            # Bakes your UID into image
scripts/build.sh --personal --github-pat ghp_...
```

**How UID remapping works:**
- **Generic image (devuser:1000)**: VS Code remaps 1000 → your UID automatically via `updateRemoteUserUID: true`
- **Personal image (your UID)**: No remapping needed, but NOT shareable with others
- **Best practice**: Use generic build + push to registry for team sharing

**Manual generic build:**
```bash
docker build . \
  -f docker/base/Dockerfile \
  --build-arg GITHUB_PAT=$GITHUB_PAT \
  -t scdock-r-dev:v0.5.3
# Note: USER_ID defaults to 1000 (no need to specify)
```

**Manual personal build:**
```bash
docker build . \
  -f docker/base/Dockerfile \
  --build-arg GITHUB_PAT=$GITHUB_PAT \
  --build-arg USER_ID=$(id -u) \
  --build-arg GROUP_ID=$(id -g) \
  --build-arg USER=$USER \
  --build-arg GROUP=$(id -gn) \
  -t scdock-r-dev:v0.5.3-personal
```

### Building ArchR Wrapper Image

The ArchR wrapper provides UID-compatible layer over official ArchR image:

```bash
./build-archr-wrapper.sh                   # Generic build (devuser:1000)
./build-archr-wrapper.sh --personal        # Personal build (your UID)
```

**What it does:**
- Pulls official `greenleaflab/archr:1.0.3-base-r4.4`
- Removes `rstudio` user, creates `devuser` with same pattern as base image
- Result: Consistent UID handling across dev-core ↔ dev-archr switching

**Why use wrapper instead of official image directly:**
- Official image has `rstudio:1000` - different user from `devuser:1000`
- Wrapper ensures both images use identical user setup
- Seamless switching between dev-core and dev-archr services
- VS Code UID remapping works consistently

**Note:** Custom ArchR Dockerfile (`.devcontainer/Dockerfile.archr`) is deprecated. Use wrapper instead.

### Extracting renv lockfile after first build

After the first successful build (when renv snapshots the packages):

```bash
CID=$(docker create scdock-r-dev:v0.5.3)
docker cp $CID:/opt/settings/renv.lock ./renv.lock
docker cp $CID:/opt/settings/R-packages-manifest.csv ./R-packages-manifest.csv
docker rm $CID
git add renv.lock R-packages-manifest.csv
git commit -m "Pin R via renv; add manifest"
```

Then add the following to `docker/base/Dockerfile` (after line ~95) to enable deterministic builds:
```dockerfile
COPY renv.lock /opt/settings/renv.lock
```

### Running sanity checks

```bash
docker run --rm scdock-r-dev:v0.5.3 bash -lc 'scripts/poststart_sanity.sh'
docker run --rm scdock-r-archr:v0.5.3 bash -lc 'scripts/poststart_sanity.sh'
```

## Architecture

### Image Layering

```
ubuntu:22.04
  └─ scdock-r-dev:v0.5.0 (docker/base/Dockerfile)
       • R 4.5.0 built from source with Cairo, BLAS, LAPACK
       • ~80 core R packages pre-installed (install_R_core.R)
       • Python base venv: /opt/venvs/base
       • Layered venvs created on-demand: {squid,atac,comms}
       • CLI tools: samtools, bcftools, bedtools, scIBD
       • TinyTeX (lightweight)
       • Quarto for scientific documentation
       • httpgd for VS Code R graphics

greenleaflab/archr:1.0.3-base-r4.4 (official, standalone)
       • R 4.4.x + ArchR 1.0.3
       • Maintained by ArchR developers
       • Use via docker-compose "dev-archr" service
```

### Docker Compose vs Single Image

The repo supports two devcontainer approaches:

1. **Docker Compose** (current setup in `.devcontainer/devcontainer.json`):
   - Defines two services: `dev-core` and `dev-archr`
   - Switch by changing `service` field in `devcontainer.json`
   - Both services use identical mounts and UIDs

2. **Single Image**: Directly specify `image` instead of `dockerComposeFile` in `devcontainer.json`

### Key Scripts

**Build-time:**
- `.devcontainer/install_R_core.R`: Installs ~80 core R packages (Seurat, edgeR, limma, etc.)
- `.devcontainer/install_R_packages.R`: Full package installer (deprecated in favor of core + runtime)
- `.devcontainer/install_renv_project.R`: Handles renv init/restore/snapshot
- `.devcontainer/install_httpgd.R`: Installs httpgd with CRAN-first, GitHub-fallback
- `.devcontainer/install_quarto.sh`: Installs Quarto
- `.devcontainer/create_layered_venv.sh`: Helper to create Python venvs with `--system-site-packages`

**Runtime:**
- `.devcontainer/scripts/poststart_sanity.sh`: Validates container environment on startup
- `init-project.sh`: Project scaffolding script (creates from templates)

**Deprecated:**
- `.devcontainer/install_R_archr.R`: ArchR installer (use official image instead)
- `.devcontainer/Dockerfile.archr`: Custom ArchR build (use official image instead)

### Python Virtual Environments (v0.5.0: Layered Approach)

**Base venv (pre-installed in image):**
- **/opt/venvs/base**: Core single-cell stack (scanpy, scvi-tools, muon, cellrank, scvelo, radian)
  - Pre-installed during build (~25GB)
  - Requirements: `docker/requirements/base.txt`

**Layered venvs (created at runtime with `--system-site-packages`):**
- **/opt/venvs/squid**: Spatial transcriptomics (squidpy, spatialdata)
  - Inherits from base, adds only squidpy-specific packages (~3-5GB additional)
  - Requirements: `docker/requirements/squid.txt`

- **/opt/venvs/atac**: scATAC-seq tools (snapatac2, episcanpy)
  - Inherits from base, adds only ATAC-specific packages (~2-3GB additional)
  - Requirements: `docker/requirements/atac.txt`

- **/opt/venvs/comms**: Cell communication and GRN (liana, cellphonedb, scglue)
  - Inherits from base, adds only communication tools (~3-4GB additional)
  - Requirements: `docker/requirements/comms.txt`

**Creating layered venvs:**
```bash
# Automatic on first use:
usepy squidpy  # Creates if doesn't exist, then activates

# Manual creation:
create_layered_venv.sh squidpy squidpy_requirements.txt
create_layered_venv.sh atac atac_requirements.txt
create_layered_venv.sh comms comms_requirements.txt
```

**Size savings:** 4 full venvs (100GB) → 1 base + 3 layered (~40GB total)

### R Environment (v0.5.0: Core + Runtime)

**R Version:** 4.5.0 + Bioconductor 3.21

**Pre-installed Core Packages (~80):**
- **Seurat ecosystem**: Seurat, Signac, BPCells, presto, glmGamPoi, sctransform
- **RNA-seq foundations**: edgeR, limma, DESeq2, scran, scater
- **GSEA & pathways**: clusterProfiler, GSVA, fgsea, msigdbr, decoupleR
- **Multi-factorial**: muscat, MOFA2, mixOmics, lemur, liger, harmony
- **Visualization**: ggplot2, patchwork, ggpubr, ComplexHeatmap, pheatmap
- **Statistics**: lme4, brms, Matrix, future
- **Interoperability**: anndataR, MuDataSeurat, sceasy, reticulate
- **VS Code support**: languageserver, httpgd, Cairo

Full list in `.devcontainer/install_R_core.R`

**Runtime Package Installation:**
Additional packages can be installed at runtime into user library (`~/R/...`):
```r
# Automatically installs to writable user library
if (!require("PACKAGE")) BiocManager::install("PACKAGE")
```

**Reproducibility:**
- **Image**: `renv.lock` at `/opt/settings/renv.lock` (extract after first build)
- **Project**: Project-level `renv.lock` for per-project snapshots
- **ArchR library** (ArchR image only): Separate library path at `$ARCHR_LIB` (default: `~/R/archr-lib`)

**Key Changes:**
- Heavy annotation packages (BSgenome.*, EnsDb.*, org.*.eg.db) NOT pre-installed - install at runtime if needed
- Core packages cover 95% of workflows, specialized packages installed on-demand
- User library is writable for runtime installs (no sudo needed)

### AI Tools Integration (Runtime Setup)

AI tools (Claude Code, MCP servers) are installed **at runtime** via **SciAgent-toolkit**, not during image build. This provides:
- Per-project configuration with correct absolute paths
- Clean separation between container infrastructure and AI tooling
- User control over when to spend time on optional installations

**Why runtime instead of build-time:**
1. `.mcp.json` requires absolute paths - must be generated with actual project path
2. AI context files (CLAUDE.md, etc.) are project-specific
3. ToolUniverse creates `tooluniverse-env/` per-project for isolation
4. Avoids maintaining separate AI-enabled image

**Setup workflow:**
1. Run `init-project.sh --with-submodules` to create project structure with SciAgent-toolkit submodule
2. Open project in VS Code Dev Container
3. Run: `./01_modules/SciAgent-toolkit/scripts/setup-ai.sh`
4. Fill in `context.md` with your scientific question

**Note:** SciAgent-toolkit is maintained at `toolkits/SciAgent-toolkit/` in the scbio-docker repo and added as a submodule to `01_modules/SciAgent-toolkit/` in new projects.

**Configuration files created:**
- `.mcp.json` - MCP server configuration (project-local, gitignored)
- `.claude/agents/` - Symlinks to active agents from SciAgent-toolkit
- `.claude/skills/` - Symlinks to active skills from SciAgent-toolkit
- `CLAUDE.md`, `GEMINI.md`, `AGENTS.md` - AI context files
- `context.md` - Scientific context (user fills in)
- `02_analysis/config/analysis_config.yaml` - Project parameters

**First-run timing:**
- Full setup: 5-15 minutes (optional Serena Rust compilation)
- Minimal setup (`--minimal`): 2-3 minutes
- Subsequent runs: <1 minute (checks if already installed)

**Available MCP servers:**
- **Sequential Thinking**: Structured reasoning for complex decisions
- **PAL**: Multi-model AI collaboration
- **Context7**: Up-to-date library documentation
- **ToolUniverse**: 600+ scientific tools (ChEMBL, UniProt, PubMed, etc.)
- **Serena**: Code intelligence (optional, requires compilation)

**Methodology Guidelines:**
Analysis guidelines are maintained in `SciAgent-toolkit/docs/guidelines/` (single source of truth).
These are referenced from AI context files, not copied to each project.

## Common Commands

### Python Environment Switching

```bash
# Interactive shell switching
usepy base       # switch to base env
usepy squid      # switch to squid env
usepy atac       # switch to ATAC env
usepy comms      # switch to COMMS env

# One-off commands
py-base python -V
py-squid python -c "import squidpy"
py-atac python -c "import snapatac2"
py-comms python -c "import pyscenic"

# Check active environment
which python && python -V
```

### R Session Management (Radian + Tmux Workflow)

**Standard R terminal:** Radian (pre-installed in base venv)

Inside **scdock-r-dev** container:
```bash
r-base           # Launch radian with base R libraries only
radian           # Direct invocation (same as r-base)
```

Inside **greenleaflab/archr** container (official ArchR image):
```bash
radian           # Launch radian (ArchR available directly)

# In R:
library(ArchR)
packageVersion("ArchR")  # 1.0.3
```

**Persistent R Sessions with Tmux (recommended for SSH):**
```bash
# Start persistent R session in tmux
tmux new-session -s my-analysis radian

# Detach: Ctrl+B, then D
# Session continues running even if SSH disconnects

# Re-attach later:
tmux attach -t my-analysis

# List sessions:
tmux ls
```

**VS Code Integration:**
- R code sent to active terminal (radian) via VS Code R extension
- Set in `.vscode/settings.json`:
  ```jsonc
  "r.rterm.linux": "/opt/venvs/base/bin/radian",
  "r.alwaysUseActiveTerminal": true,
  "r.bracketedPaste": true
  ```
- Start radian in tmux → select that terminal in VS Code → send R code
- Graphics display via httpgd (interactive plots in VS Code)

### Switching to ArchR Container

**Method 1: Edit devcontainer.json (persistent)**
```bash
# In .devcontainer/devcontainer.json, change:
"service": "dev-core"  →  "service": "dev-archr"

# In VS Code: Cmd/Ctrl+Shift+P → "Dev Containers: Rebuild and Reopen in Container"
```

**Method 2: Attach to running container (temporary)**
```bash
docker compose -f .devcontainer/docker-compose.yml up -d dev-archr
# In VS Code: Cmd/Ctrl+Shift+P → "Dev Containers: Attach to Running Container"
```

See DEVOPS.md for complete workflow.

### Installing Additional R Packages

**Two-Tier R Library Architecture:**

The container uses a read-only system library + writable user library design:

```
System Library (read-only)          User Library (writable)
/usr/local/lib/R/library            ~/R/x86_64-pc-linux-gnu-library/4.5
├─ Core packages (~80)              ├─ Runtime installs
├─ Pinned via renv.lock            ├─ Project-specific packages
├─ Owned by root                   ├─ Owned by devuser
└─ Same across all containers      └─ Takes precedence over system
```

**Installing packages (no sudo needed):**

```r
install.packages("PACKAGE")
BiocManager::install("PACKAGE")

# Check library paths:
.libPaths()
# [1] "/home/devuser/R/x86_64-pc-linux-gnu-library/4.5"  # writable (runtime installs)
# [2] "/usr/local/lib/R/library"                          # system (core packages)
```

**Expected Warning (NORMAL):**

When installing with BiocManager, you'll see:
```r
BiocManager::install("AnnotationHub")
...
* DONE (AnnotationHub)

Installation paths not writeable, unable to update packages
  path: /usr/local/lib/R/library
  packages:
    aplot, BiocGenerics, Matrix, Seurat, ...
```

**This is expected and harmless:**
- ✅ Your package installed successfully to user library
- ⚠️ BiocManager checked if system packages need updates (by design with `update=TRUE`)
- ❌ Cannot update system packages because they're read-only (intentional for reproducibility)

**To suppress warnings:**
```r
BiocManager::install("PACKAGE", update = FALSE)
```

**Why read-only system library?**
1. Reproducibility: Core packages pinned via renv.lock
2. Shareability: Same baseline for all users
3. Disk efficiency: ~10GB core packages shared across containers
4. Separation: System packages stable, user packages experimental

### Container Operations

**Run base image manually:**
```bash
docker run --rm -it \
  -u $(id -u):$(id -g) \
  -v /path/to/project:/workspaces/project \
  --memory=450g --cpus=50 \
  scdock-r-dev:v0.5.3 bash
```

**Run ArchR image manually:**
```bash
docker run --rm -it \
  -u $(id -u):$(id -g) \
  -v /path/to/project:/workspaces/project \
  --memory=450g --cpus=50 \
  greenleaflab/archr:1.0.3-base-r4.4 bash
```

**Using Docker Compose:**
```bash
export LOCAL_UID=$(id -u); export LOCAL_GID=$(id -g); export WORKSPACE_FOLDER=$PWD
docker compose -f .devcontainer/docker-compose.yml up -d dev-archr
docker compose -f .devcontainer/docker-compose.yml exec dev-archr bash
docker compose -f .devcontainer/docker-compose.yml down
```

## Code Architecture

### Package Management Strategy

**R packages:**
- Deterministic builds via renv lockfile at `/opt/settings/renv.lock`
- CRAN snapshot controlled by `RSPM_SNAPSHOT` env var (default: 2025-02-15)
- Bioconductor version pinned to 3.21 (for R 4.5)
- Core packages (~80) pre-installed via `install_R_core.R`
- Heavy annotation packages (BSgenome.*, EnsDb.*, org.*.eg.db) are optional; install with `--build-arg INCLUDE_HEAVY_R_DATA=1` or at runtime
- **ArchR**: Use official image instead of custom build

**Python packages:**
- Pinned versions in `docker/requirements/*.txt` files
- Base venv fully resolved during image build
- Layered venvs created on-demand with `--system-site-packages`
- No conflicting dependencies between venvs

**CLI tools:**
- Explicit version pinning in Dockerfile (samtools 1.21, bcftools 1.21, bedtools 2.31.1)

### UID/GID Handling (Generic + Runtime Remapping)

**New Strategy (v0.5.3): Generic Images with Runtime UID Remapping**

Images are built with **generic user (devuser:1000)** by default, then mapped to actual user at runtime:

**How it works:**
1. **Build**: Create generic `devuser:1000` (shareable image)
2. **Runtime**: Remap 1000 → your actual UID
   - **VS Code**: Automatic via `updateRemoteUserUID: true` in devcontainer.json
   - **Docker Compose**: Manual via `LOCAL_UID=${LOCAL_UID:-1000}` env vars
   - **Docker run**: Manual via `-u $(id -u):$(id -g)`

**Benefits:**
- ✅ **Single image works for everyone** (shareable via registry)
- ✅ **Runs as YOUR UID** when you use it (trackable in htop)
- ✅ **Runs as Bob's UID** when Bob uses it (auto-remapped)
- ✅ **No permission conflicts** (files owned by actual user)
- ✅ **Team-friendly** (build once, everyone can use)

**VS Code Setup (Automatic Remapping):**
```json
// .devcontainer/devcontainer.json
{
  "remoteUser": "devuser",
  "updateRemoteUserUID": true  // ← Remaps 1000 to your UID
}
```

**Docker Compose Setup (Env Var Override):**
```yaml
# docker-compose.yml
services:
  dev-core:
    user: "${LOCAL_UID:-1000}:${LOCAL_GID:-1000}"  # Override via env
```

```bash
# .env file
LOCAL_UID=788715489  # Your UID (run: id -u)
LOCAL_GID=788600513  # Your GID (run: id -g)
```

**Manual Docker Run:**
```bash
docker run -u $(id -u):$(id -g) scdock-r-dev:v0.5.3
```

**Active Directory / Special Characters:**
The Dockerfile still handles AD users with group names containing spaces:
- Normalizes group names by replacing spaces with underscores
- Falls back to `grp_<gid>` if group creation fails
- Uses numeric IDs for `chown` operations to avoid name parsing issues

**Legacy (Personal Build):**
If you need a personal-only image (not shareable):
```bash
scripts/build.sh --personal  # Bakes YOUR UID into image
```

### R Profile and Startup

- `.Rprofile` is interactive-only and VS Code-aware
- Enables httpgd only when running inside VS Code
- Site-wide CRAN mirror set at `/etc/R/Rprofile.site` for non-interactive scripts
- **Note:** `USE_ARCHR` mechanism deprecated; use official ArchR image instead

### Image Slimming Decisions (v0.5.0)

**Removed from base image** (to reduce size):
- Bulk aligners: STAR, BWA, Bowtie2, Salmon, kallisto
- Pre-processing tools: FastQC, Trimmomatic, Trim Galore, featureCounts, Picard
- Heavy R annotation packages: moved to optional install
- R/Python bridging: pyreadr, rpy2, anndata2ri (commented in base_requirements.txt)
- **Full R package installs**: Replaced with ~80 core packages
- **Multiple full Python venvs**: Replaced with base + layered venvs
- **Full TeX distribution**: Replaced with TinyTeX
- **Custom ArchR build**: Use official image instead

**Retained** (essential for downstream analysis):
- samtools, bcftools, bedtools (epigenomics-friendly workflows)
- scIBD (scATAC doublet detection)
- MACS3 (Python)

**Size Impact:**
- Before: 500GB reported, ~200GB actual filesystem
- After: 533GB reported (Docker layer issue), ~20GB actual filesystem

## Interoperability Workflows

### R ↔ Python Round-Trip

For multi-modal or integrated analyses (e.g., scRNA + scATAC):

**Export from R (Seurat):**
```r
library(MuDataSeurat)
WriteH5AD(object = s, file = "rna.h5ad", assay = "RNA")  # single assay
WriteH5MU(s, "multiome.h5mu")                            # multiome
```

**Process in Python:**
```bash
usepy base
```
```python
import scanpy as sc, muon as mu
m = mu.read_h5mu("multiome.h5mu")
rna = m.mod["RNA"]; atac = m.mod["ATAC"]
# Run scVI/PeakVI/SCGLUE; store embeddings in obsm['X_scvi'], etc.
mu.write_h5mu("multiome.h5mu", m)
```

**Import back to R:**
```r
library(MuDataSeurat)
s_multi <- ReadH5MU("multiome.h5mu")
names(s_multi@reductions)  # Check for pca, umap, scvi, peakvi, etc.
```

### GRN Inference (SCENIC+)

Use the COMMS venv for cell communication and GRN tools:

```bash
usepy comms
python -c "import pycistarget, pycistopic"
```

If `scenicplus` is not preinstalled:
```bash
pip install 'scenicplus @ git+https://github.com/aertslab/SCENICplus.git'
```

## Project Setup Pattern

### Separation of Concerns

This repository follows a strict separation between **container infrastructure** and **AI tooling**:

| Repository | Responsibility |
|------------|----------------|
| **scbio-docker** | Docker images, container setup, project directory structure |
| **SciAgent-toolkit** | AI tools installation, MCP servers, agents/skills, methodology guidelines |

### Quick Start

```bash
# From scbio-docker repository
./init-project.sh ~/projects/my-analysis              # Uses "base" template (default)
./init-project.sh ~/projects/my-analysis base --git-init --with-submodules
```

**What init-project.sh creates:**
```
my-project/
├── .devcontainer/
│   ├── devcontainer.json      # Pre-configured VS Code container setup
│   ├── docker-compose.yml     # Multi-service setup (dev-core + dev-archr)
│   ├── .env                   # LOCAL_UID, LOCAL_GID, API keys
│   └── scripts/
│       └── poststart_sanity.sh
├── .vscode/
│   └── settings.json          # Universal Python + R configuration
├── .gitignore                 # Excludes data/, results/, caches
├── 00_data/
│   ├── raw/                   # Raw data (add mount in docker-compose.yml)
│   ├── processed/             # Processed objects
│   └── references/            # Reference files
├── 01_modules/                # Reusable toolkits (submodules)
│   └── RNAseq-toolkit/        # R analysis functions (optional)
├── 02_analysis/
│   ├── config/                # Configuration files
│   │   ├── config.R
│   │   ├── pipeline.yaml
│   │   └── color_config.R
│   └── helpers/               # Project-specific scripts
├── 03_results/
│   ├── checkpoints/           # Cached intermediate objects
│   ├── plots/                 # Generated figures
│   └── tables/                # Output tables
├── logs/
├── README.md                  # Project instructions
├── tasks.md                   # Task tracker
└── notes.md                   # Research notes
```

### After Container Setup (AI Tooling)

Run inside the container to set up AI tools:
```bash
./01_modules/SciAgent-toolkit/scripts/setup-ai.sh
```

**What setup-ai.sh creates:**
```
my-project/
├── .claude/                   # Claude Code configuration
│   ├── agents/                # Symlinks to active agents
│   └── skills/                # Symlinks to active skills
├── .mcp.json                  # MCP server configuration
├── CLAUDE.md                  # AI context for Claude
├── GEMINI.md                  # AI context for Gemini
├── AGENTS.md                  # Universal AI rules
├── context.md                 # Scientific context (user fills in)
└── 02_analysis/config/
    └── analysis_config.yaml   # Project parameters
```

### Manual Setup (Without Templates)

1. Copy `.devcontainer/` folder from scbio-docker repo
2. Copy `.vscode/settings.json` from `templates/.vscode/settings.json`
3. Edit `docker-compose.yml` to add data mounts
4. Open in VS Code: `Dev Containers: Reopen in Container`

See DEVOPS.md for complete instructions.

## Troubleshooting

### Permission Issues

Check mount ownership if workspace writes fail:
```bash
id -u; id -g; id -gn
ls -ld /workspaces/project
# On host if needed:
sudo chown -R $(id -u):$(id -g) /path/to/project
```

### R Library Not Writable

The sanity script checks this at startup. Ensure user library path is writable:
```r
.libPaths()  # First entry should be writable
```

### Missing ArchR

If `library(ArchR)` fails, you're using the `dev-core` service which doesn't include ArchR.

**Solution:** Switch to `dev-archr` service (uses official ArchR image)
1. Edit `.devcontainer/devcontainer.json`: change `"service": "dev-core"` to `"service": "dev-archr"`
2. VS Code: `Cmd/Ctrl+Shift+P` → "Dev Containers: Rebuild and Reopen in Container"

See "Switching to ArchR" section in DEVOPS.md.

### GitHub Rate Limits

Set `GITHUB_PAT` during build:
```bash
export GITHUB_PAT=ghp_your_token_here
# build commands...
unset GITHUB_PAT
```

## Notes

- Current version: v0.5.3
- Default branch for PRs: `main`
- This repository uses git; current branch is `dev`
- Heavy annotation packages are excluded by default; enable with `--build-arg INCLUDE_HEAVY_R_DATA=1` or install at runtime
- For long-running tasks, use `tmux` inside the container to survive disconnections
- VS Code should exclude heavy I/O paths in settings to speed up indexing
