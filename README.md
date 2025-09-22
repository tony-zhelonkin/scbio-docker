# Single-Cell Docker Dev Environment

![Docker Image Version](https://img.shields.io/badge/Docker-v0.4.0-blue?style=flat-square)
![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)

A Docker-based development environment for bioinformatics, particularly single-cell RNA-seq analyses. This repository is structured for use with Visual Studio Code’s **Remote - Containers** extension, enabling seamless local or remote development against powerful server resources.

> **Inspiration**: This setup is **strongly inspired by** [Rami Krispin’s vscode-r repository](https://github.com/RamiKrispin/vscode-r). Special thanks for the excellent reference on R + VS Code configurations!

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

- **R** (v4.4.2) and comprehensive single-cell and genomic libraries (Seurat, ArchR, Bioconductor ecosystem, etc.).
- **Python** (v3.10) with scientific libraries (NumPy, SciPy, Pandas, scVI, Scanpy, etc.).
- Command-line tools focused on epigenomics-friendly workflows: **samtools**, **bcftools**, **bedtools**, and selected extras (e.g., **scIBD** for scATAC doublet detection). Bulk aligners and pre-processing tools (STAR, BWA, Salmon, kallisto, Picard, FastQC/Trimmomatic, Trim Galore, featureCounts, etc.) are intentionally excluded to keep the image slim. Use pipeline-specific containers if you need them.

The Docker setup is integrated with **VS Code Remote Containers** to streamline development on remote compute nodes—ideal for large-scale single-cell data processing.

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
  -f .devcontainer/Dockerfile \
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

- Deterministic builds thereafter: uncomment the lock copy in `.devcontainer/Dockerfile`:

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

### Python environments (base vs squid)

Two venvs are preinstalled and fully resolved during the image build:
- `/opt/venvs/base` from `/opt/environments/base_requirements.txt`
- `/opt/venvs/squid` from `/opt/environments/squid_requirements.txt`

Default PATH uses the base venv. Switch interactively inside the container:

```bash
usepy squid      # switch shell to squid env
usepy atac       # switch shell to ATAC env (snapatac2)
usepy base       # switch back to base env
```

One-off commands under a specific env:

```bash
py-squid python -V
py-atac python -V
py-base python -V
```

Verify which env is active:

```bash
which python && python -V && pip list | head
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
│   ├── install_quarto.sh             # Quarto installer
│   ├── install_R_packages.R          # R packages (CRAN/Bioc/GitHub); used with renv
│   └── .Rprofile                     # R profile (VS Code/httpgd; CRAN mirror; ArchR toggle)
├── .environments
│   ├── base_requirements.txt         # Base Python environment
│   └── squid_requirements.txt        # Alternative Python environment (squid)
├── README.md                         # This README
└── LICENCE.md
```

- **`.devcontainer/Dockerfile`** is the primary Dockerfile for building R, Python envs, and CLI tools.
- **`.devcontainer/install_R_packages.R`** installs R/Bioconductor/GitHub packages; paired with `renv` for pinning.
- **`.environments/*.txt`** are the two Python requirement sets installed into separate venvs inside the image.

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
    -f ./.devcontainer/Dockerfile \
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
  "service": "base",
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

See the [Dockerfile](./.devcontainer/Dockerfile) and [install_R_packages.R](./.devcontainer/install_R_packages.R) for full details on what gets installed.

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

## Future directions

* R packages to be added 
```
install.packages('msigdbdf', repos = 'https://igordot.r-universe.dev')
```

* Python env update 
The python pip modules are installed not into the venv right now, and so python doesn`t see imported modules. 

I activated pip3 directly. And this lead to incorrect path of installation.
>[!error]
>```
># Set up virtual environment
>RUN python3 -m venv /opt/$VENV_NAME \
>&& export PATH=/opt/$VENV_NAME/bin:$PATH \
>&& echo "source /opt/$VENV_NAME/bin/activate" >> ~/.bashrc
>
># Install dependencies first
>RUN pip3 install setuptools wheel
>RUN pip3 install numpy scipy pandas
>RUN pip3 install anndata==0.10.9
>RUN pip3 install "python-igraph==0.10.4"
>
># Install requirements
>RUN python3 -m pip install --upgrade pip setuptools wheel && \
>pip3 install -r ./settings/requirements.txt
>```

Instead of using pip3 directly, I should first activate the venv, and then install inside it
```
# Instead of using pip3 directly, activate the virtual environment first
RUN . /opt/$VENV_NAME/bin/activate && \
    python -m pip install setuptools wheel && \
    python -m pip install numpy scipy pandas && \
    python -m pip install anndata==0.10.9 && \
    python -m pip install "python-igraph==0.10.4" && \
    python -m pip install --upgrade pip setuptools wheel && \
    python -m pip install -r ./settings/requirements.txt
```

Right now got away with the 
```
python -m pip install --upgrade -r ./settings/requirements.txt
```
Where the requirements.txt is 
```
setuptools
wheel
numpy
scipy
pandas
anndata==0.10.9
python-igraph==0.10.4

# Core Scientific Computing
#numpy (pre-install for better dependency handling)
#scipy
#pandas
numba
statsmodels
umap-learn
scikit-learn

# Data Storage & Processing
h5py
cooler
tables
htseq
PySam
SWIG
Cython==0.29.32
pybedtools==0.10.0
pyBigWig==0.3.23

# Feature counters
TEtranscripts==2.2.3

# Visualization Libraries
matplotlib
seaborn
plotly
bokeh

# Single-Cell Analysis
scanpy==1.10.4
scvi-tools==1.2.0
cellrank==2.0.6
scvelo==0.3.2
snapatac2==2.7.1
scrublet==0.2.3
harmony-pytorch
scib==1.1.5
scirpy==0.19.0
squidpy==1.6.2
muon==0.1.7

# Epigenomics 
MACS3==3.0.2
episcanpy==0.4.0
deeptools

# Network Analysis
python-louvain==0.16
networkx==3.4.2
leidenalg==0.10.1
pyscenic==0.12.1


# Quality Control & Processing
RSeQC
multiqc==1.25.2
cutadapt

# R and Python Integration
pyreadr==0.5.2
rpy2==3.5.17
anndata2ri==1.3.2
```




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