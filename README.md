# Single-Cell Docker Dev Environment

![Docker Image Version](https://img.shields.io/badge/Docker-v0.2-blue?style=flat-square)
![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)

A Docker-based development environment for bioinformatics, particularly single-cell RNA-seq analyses. This repository is structured for use with Visual Studio Code’s **Remote - Containers** extension, enabling seamless local or remote development against powerful server resources.

> **Inspiration**: This setup is **strongly inspired by** [Rami Krispin’s vscode-r repository](https://github.com/RamiKrispin/vscode-r). Special thanks for the excellent reference on R + VS Code configurations!

---

## Table of Contents
1. [Overview](#overview)
2. [Features](#features)
3. [Repository Structure](#repository-structure)
4. [Getting Started](#getting-started)
    - [Prerequisites](#prerequisites)
    - [Build the Image](#build-the-image)
    - [Run the Container](#run-the-container)
    - [Using with VS Code Remote Containers](#using-with-vs-code-remote-containers)
5. [Typical Workflow for New Projects](#typical-workflow-for-new-projects)
6. [Included Tools & Packages](#included-tools--packages)
7. [Troubleshooting](#troubleshooting)
8. [Contributing](#contributing)
9. [License](#license)
10. [Acknowledgments](#acknowledgments)

---

## Overview

This repository contains Dockerfiles and configuration files for creating a **versatile bioinformatics environment** that focuses on single-cell RNA-seq and epigenomics workflows. It packages:

- **R** (v4.4.2) and comprehensive single-cell and genomic libraries (Seurat, ArchR, Bioconductor ecosystem, etc.).
- **Python** (v3.10) with scientific libraries (NumPy, SciPy, Pandas, scVI, Scanpy, etc.).
- Command-line tools (e.g., **bowtie2**, **chromap**, **samtools**, **bcftools**, **STAR**, **kallisto**, **salmon**, **bedtools**, **FastQC**, **Trimmomatic**, **Picard**, **Trim Galore**, **featureCounts**, **BWA**, and many more).

The Docker setup is integrated with **VS Code Remote Containers** to streamline development on remote compute nodes—ideal for large-scale single-cell data processing.

---

## Features

- **R** with a rich set of CRAN and Bioconductor packages for single-cell and bulk RNA-seq analysis.
- **Python** environment (virtualenv) preloaded with popular data science packages, single-cell toolkits, and epigenomics libraries.
- **CLI Tools** for alignment, quality control, and post-processing (STAR, bowtie2, bcftools, samtools, salmon, etc.).
- **Quarto** for R Markdown, Jupyter, and scientific documentation workflows.
- **Automated Build** process (via Dockerfile.dev) including compilation of R, specialized R packages, and Python requirements.
- **VS Code** integration with recommended extensions for R, Python, Docker, Markdown, and more.

---

## Repository Structure

```bash
.
├── .devcontainer
│   ├── devcontainer.back.json        # Backup devcontainer config (to be deleted in the next update)
│   ├── devcontainer.json             # VS Code Remote Containers template configuration
│   ├── Dockerfile                    # Backup base Dockerfile (to be deleted in the next update)
│   ├── Dockerfile.dev               # Latest successfully built Dockerfile
│   ├── install_quarto.sh            # script installing Quarto
│   ├── install_R_packages.R          # script installing R packages (CRAN/Bioc/GitHub)
│   ├── requirements.txt             # Python packages to install
│   └── .Rprofile                     # R profile tweaks/aliases
├── .git                              # Git versioning directory
├── .gitignore
├── README.md                         # This README
└── .vscode
    └── settings.json                 # Example VS Code settings (R, radian, etc.)
```

- **`.devcontainer/Dockerfile.dev`** is the primary Dockerfile used to build the environment with R and Python tools.
- **`.devcontainer/install_R_packages.R`** outlines and installs a comprehensive set of R/Bioconductor packages.
- **`.devcontainer/requirements.txt`** details Python packages installed via `pip`.

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
    -f ./.devcontainer/Dockerfile.dev \
    --build-arg GITHUB_PAT=<your_github_pat> \
    --build-arg R_VERSION_MAJOR=4 \
    --build-arg R_VERSION_MINOR=4 \
    --build-arg R_VERSION_PATCH=2 \
    -t scdock-r-dev:v0.2
```

- **`GITHUB_PAT`** (optional) is used to authenticate GitHub requests during R package installation (especially if you hit rate limits).
- Adjust **R_VERSION_*** build arguments as needed to override the default R version. But beware the build was only tested by myself for the R 4.4.2, and some of the R libraries depend on R version no older than 4.3-4.4


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

3. **Add a `devcontainer.json`** to `.devcontainer`. Below is a sample template that references the `scdock-r-dev:v0.2` image and mounts data directories read-only:
   ```jsonc
   // .devcontainer/devcontainer.json
   {
       "name": "Yasmine-RetroT",
       "image": "scdock-r-dev:v0.2", 
       "runArgs": [
           "--memory=58g",
           "--cpus=9"
       ],
       "customizations": {
           "vscode": {
               "extensions": [
                   // R Extensions
                   "reditorsupport.r",
                   //"rdebugger.r-debugger", // if nedeed
                   
                   // Documentation Extensions
                   "purocean.drawio-preview",
                   "yzhang.markdown-all-in-one",
                   //"redhat.vscode-yaml", // if needed
                   // //"quarto.quarto", // if needed
                   
                   // Docker Supporting Extensions
                   "ms-azuretools.vscode-docker",
                   "ms-vscode-remote.remote-containers",
                   
                   // Python Extensions
                   "ms-python.python",
                   //"ms-toolsai.jupyter" // Disabled if it causes lag
               ]
           }
       },
       "mounts": [
           "source=/path/to/data/reads, target=/workspaces/project/0_Data/data, type=bind, readonly",
           "source=/path/to/ref_genome/, target=/workspaces/project/0_Data/ref_genome,type=bind,readonly",
       ],
       "postStartCommand": "echo 'Hello World!' "
   }
   ```

4. **Add a `settings.json`** to `.vscode` to configure R settings, radian, etc.:
   ```jsonc
   // .vscode/settings.json
   {
       "r.alwaysUseActiveTerminal": true,
       "r.bracketedPaste": true,
       "r.rterm.linux": "/usr/local/bin/radian",
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
   - Select the **Dev Container** you just defined.

6. **Long-running tasks**: For tasks that must continue even if you disconnect (like heavy alignment or processing), you can use `tmux` inside the container:
   ```bash
   tmux new -s your_session_name
   # run your alignment/analysis steps
   # detach with Ctrl+B, then D
   ```
   This ensures your process remains running if your SSH or VS Code session is interrupted.

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
   - Common NGS tools: **bowtie2**, **STAR**, **samtools**, **bcftools**, **salmon**, **bedtools**, **bwa**, **FastQC**, **Trimmomatic**, **Picard**, **featureCounts**, **Trim Galore**, etc.

4. **Quarto** for rendering R Markdown, Jupyter Notebooks, and other scientific documents.

See the [Dockerfile.dev](./.devcontainer/Dockerfile.dev) and [install_R_packages.R](./.devcontainer/install_R_packages.R) for full details on what gets installed.

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
       -f ./.devcontainer/Dockerfile.dev \
       --build-arg GITHUB_PAT=<your_github_pat> \
       -t scdock-r-dev:latest
   ```

4. **R Package Version Mismatch**:  
   Libraries compiled under a different R version might cause “object not found” or `.so` loading errors. Rebuild the container with consistent R version arguments or reinstall packages inside the container to match R’s version.

---

## Future directions

* R packages to be added 
```
install.packages("msigdbr")
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

This project is distributed under the [MIT License](./LICENSE). Feel free to use, modify, and distribute it as permitted.

---

## Acknowledgments

- **Maintainer**: [Anton Zhelonkin](mailto:anton.bioinf.md@gmail.com)
- **Huge thanks** to [Rami Krispin’s vscode-r repo](https://github.com/RamiKrispin/vscode-r) for serving as a fantastic inspiration.
- Thanks to all authors of open-source bioinformatics tools included here.