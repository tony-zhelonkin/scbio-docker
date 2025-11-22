# Runtime Package Installation Guide

This guide explains how to install additional R packages, Python packages, and system tools inside the running container.

## Philosophy

The v0.5.1+ image uses a **multi-stage build** to dramatically reduce size (~20GB actual filesystem) while preserving the ability to install packages at runtime:

- **Build-essential and dev libraries are KEPT** in the final image
- R, Python, and system packages can be compiled from source if needed
- User has `sudo` access (NOPASSWD) for system-level installs

---

## R Package Installation

### Standard Installation (User Library)

R packages install to your **user library** by default (`~/R/x86_64-pc-linux-gnu-library/4.5`), which is writable without sudo:

```r
# CRAN packages
install.packages("ggExtra")
install.packages("patchwork")

# Bioconductor packages
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("EnsDb.Hsapiens.v86")

# GitHub packages
remotes::install_github("satijalab/seurat-data")

# Check where packages are installed
.libPaths()
# [1] "/home/user/R/x86_64-pc-linux-gnu-library/4.5"  # User library (writable)
# [2] "/usr/local/lib/R/library"                       # System library (core packages)
```

### Heavy Annotation Packages

Large annotation packages (BSgenome.*, EnsDb.*, org.*.eg.db) are **not** pre-installed to save space. Install them on-demand:

```r
# Install heavy annotation packages as needed
BiocManager::install(c(
  "BSgenome.Hsapiens.UCSC.hg38",
  "BSgenome.Mmusculus.UCSC.mm10",
  "EnsDb.Hsapiens.v86",
  "org.Hs.eg.db",
  "org.Mm.eg.db"
))
```

### Troubleshooting R Package Installation

**Problem: Package compilation fails due to missing system dependencies**

Example error:
```
configuration failed for package 'XML'
ERROR: configuration failed for package 'XML'
```

**Solution:** Install system dependencies with sudo:

```bash
# Example: Install libxml2 dev libraries
sudo apt-get update
sudo apt-get install -y libxml2-dev

# Then retry R installation
R -e 'install.packages("XML")'
```

**Common system dependencies:**

| R Package | System Dependency | Install Command |
|-----------|------------------|-----------------|
| XML | libxml2-dev | `sudo apt-get install -y libxml2-dev` |
| RCurl | libcurl4-openssl-dev | `sudo apt-get install -y libcurl4-openssl-dev` |
| rJava | openjdk-11-jdk | `sudo apt-get install -y openjdk-11-jdk` |
| sf | libgdal-dev libproj-dev | `sudo apt-get install -y libgdal-dev libproj-dev libgeos-dev` |
| rgdal | libgdal-dev | `sudo apt-get install -y libgdal-dev` |
| Cairo | libcairo2-dev | Already installed |
| png | libpng-dev | Already installed |
| jpeg | libjpeg-dev | Already installed |

**Problem: Out of memory during installation**

Large packages (e.g., BSgenome, Seurat) may fail with memory errors.

**Solution:** Install one at a time instead of in a vector:

```r
# Don't do this:
BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Mmusculus.UCSC.mm10"))

# Do this instead:
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
```

### Project-Level Reproducibility with renv

Use `renv` to snapshot your project's R package environment:

```r
# Initialize renv for your project
renv::init()

# After installing additional packages
renv::snapshot()

# Commit renv.lock to git
# Others can restore your exact environment:
renv::restore()
```

---

## Python Package Installation

### Standard Installation (Active venv)

Python packages install to the **active virtual environment**:

```bash
# Install to base venv
usepy base
pip install scikit-learn umap-learn

# Install to specialized venv (e.g., squid)
usepy squid
pip install napari  # Adds to squid venv, inherits base packages

# Freeze for reproducibility
pip freeze > requirements-project.txt
```

### Installing Packages with Compilation

Most Python packages install via wheel (pre-compiled), but some require compilation:

```bash
# Example: Install package requiring C++ compilation
usepy base
pip install pybind11  # May require g++, already installed

# If compilation fails, check error message for missing dependencies
```

**Common Python package dependencies:**

| Python Package | System Dependency | Install Command |
|---------------|------------------|-----------------|
| numpy/scipy | build-essential gfortran | Already installed |
| scikit-image | libjpeg-dev libpng-dev | Already installed |
| cartopy | libgeos-dev libproj-dev | `sudo apt-get install -y libgeos-dev libproj-dev` |
| gdal | libgdal-dev | `sudo apt-get install -y libgdal-dev` |
| h5py | libhdf5-dev | `sudo apt-get install -y libhdf5-dev` |

### Creating Additional Layered venvs

Create custom layered venvs for project-specific tools:

```bash
# Create a new layered venv manually
python3 -m venv --system-site-packages /opt/venvs/myproject

# Activate and install packages
source /opt/venvs/myproject/bin/activate
pip install custom-package

# Add to usepy function for convenience
# (Edit /etc/bash.bashrc or add to ~/.bashrc)
```

---

## System Tools Installation

### Installing Additional CLI Tools

You have `sudo` access to install system packages:

```bash
# Update package lists
sudo apt-get update

# Install tools
sudo apt-get install -y tree jq parallel

# Example: Install STAR aligner (excluded from base image for size)
sudo apt-get install -y rna-star

# Clean up after installation
sudo apt-get clean
sudo rm -rf /var/lib/apt/lists/*
```

### Installing from Source

The image includes `build-essential`, so you can compile tools:

```bash
# Example: Install a tool from source
cd /tmp
wget https://example.com/tool.tar.gz
tar -xf tool.tar.gz
cd tool
./configure --prefix=/usr/local
make -j$(nproc)
sudo make install

# Clean up
cd /
rm -rf /tmp/tool*
```

### Persistent vs. Ephemeral Installs

**Important:** System-level installs (`sudo apt-get install`, compiled tools) are **ephemeral** if you're using `docker run --rm` or rebuilding containers.

**For persistent tools:**
1. Install them in your project directory (if possible)
2. OR add them to the Dockerfile and rebuild
3. OR create a custom layer on top of the base image

**For project-specific needs:**
- Install to project directory when possible
- Document in project README for reproducibility

---

## Example Workflows

### Workflow 1: Install Heavy Annotation Packages for Analysis

```r
# Inside R (radian)
BiocManager::install(c(
  "BSgenome.Hsapiens.UCSC.hg38",
  "EnsDb.Hsapiens.v86",
  "org.Hs.eg.db"
))

# Snapshot for project reproducibility
renv::snapshot()
```

### Workflow 2: Add Spatial Analysis Tools (Python)

```bash
# Create squid venv if not exists
usepy squid

# Install additional spatial tools
pip install napari scikit-image

# Freeze
pip freeze > requirements-squid-frozen.txt
```

### Workflow 3: Install Missing System Dependency

```bash
# R package 'sf' fails to install
# Error message mentions "libgdal"

# Install system dependencies
sudo apt-get update
sudo apt-get install -y libgdal-dev libproj-dev libgeos-dev

# Retry R installation
R -e 'install.packages("sf")'
```

### Workflow 4: Install Tool from GitHub

```r
# Install development version from GitHub
remotes::install_github("satijalab/seurat", ref = "develop")
```

---

## Size Management

### Check Disk Usage

Monitor your container's disk usage:

```bash
# Check total disk usage
du -hsx /* 2>/dev/null | sort -h | tail -n 10

# Check R library size
du -hs ~/R
du -hs /usr/local/lib/R/library

# Check Python venv sizes
du -hs /opt/venvs/*
```

### Clean Up After Installation

```bash
# Clean R package caches
rm -rf ~/.cache/R /tmp/Rtmp* /tmp/downloaded_packages

# Clean pip caches
rm -rf ~/.cache/pip

# Clean apt caches (after system installs)
sudo apt-get clean
sudo rm -rf /var/lib/apt/lists/*
```

---

## Pre-installed Build Tools

The following development tools are **already installed** and available for package compilation:

### Compilers and Build Tools
- `build-essential` (gcc, g++, make)
- `gfortran` (Fortran compiler for R packages)
- `git` (version control)
- `wget`, `curl` (download tools)

### R Development Libraries
- `libcurl4-openssl-dev` (for RCurl, httr)
- `libssl-dev` (for openssl)
- `libxml2-dev` (for XML, xml2)
- `libcairo2-dev` (for Cairo graphics)
- `libfreetype6-dev`, `libpng-dev`, `libjpeg-dev` (for image processing)
- `libharfbuzz-dev`, `libfribidi-dev` (for text rendering)
- `libblas-dev`, `liblapack-dev` (for linear algebra)

### Python Development Libraries
- `python3-dev` (Python headers)
- `python3-pip`, `python3-venv` (package management)

---

## FAQ

### Q: Can I install packages without rebuilding the image?
**A:** Yes! That's the whole point. Install R packages to user library, Python packages to venvs, and system packages with `sudo apt-get`.

### Q: Will my installations persist?
**A:**
- **User library R packages:** Yes, if your home directory is mounted
- **Python venv packages:** Yes, if `/opt/venvs` is mounted (usually not)
- **System packages:** No, unless you commit the container or rebuild the image

For project reproducibility, use `renv::snapshot()` (R) and `pip freeze` (Python), then commit lockfiles to git.

### Q: How do I know if I need a system dependency?
**A:** The error message will usually tell you. Look for lines like:
```
checking for libxml2... no
configure: error: libxml2 not found
```
Then install with: `sudo apt-get install -y libxml2-dev`

### Q: Can I install STAR, BWA, or other aligners?
**A:** Yes, but they're intentionally excluded from the base image to save size. Install with:
```bash
sudo apt-get update
sudo apt-get install -y rna-star bwa
```

Or compile from source if you need specific versions.

### Q: How do I make system installs permanent?
**A:**
1. **Option 1:** Add to Dockerfile and rebuild
2. **Option 2:** Create a custom Dockerfile that extends the base image:
   ```dockerfile
   FROM scdock-r-dev:v0.5.1
   RUN sudo apt-get update && sudo apt-get install -y my-tool
   ```
3. **Option 3:** Commit running container to new image:
   ```bash
   docker commit <container_id> my-custom-image:v1
   ```

---

## Best Practices

1. **Use project-level renv/pip freeze** for reproducibility
2. **Install heavy packages on-demand** rather than pre-installing everything
3. **Clean up after installation** to save space
4. **Document system dependencies** in project README
5. **Use layered venvs** for specialized Python tools
6. **Check image/container size** periodically with `du -hsx`

---

## Support

If you encounter issues with package installation:
1. Check error messages for missing system dependencies
2. Search for the package name + "ubuntu 22.04 installation"
3. Consult package documentation for build requirements
4. Open an issue in the scbio-docker repository with error details
