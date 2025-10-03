#!/usr/bin/env Rscript

# CORE R PACKAGES for scbio-docker
# These are pre-installed in the image. Additional packages can be installed at runtime.


# Reproducible CRAN snapshot via Posit Package Manager
snapshot <- Sys.getenv("RSPM_SNAPSHOT", unset = "2025-02-15")
options(repos = c(CRAN = paste0("https://packagemanager.posit.co/cran/__linux__/jammy/", snapshot)))

# GitHub token passthrough (optional)
github_pat <- Sys.getenv("GITHUB_PAT")
if (nzchar(github_pat)) Sys.setenv(GITHUB_PAT = github_pat)

# Helper function
safe_install <- function(pkgs, installer, ...) {
  for (pkg in pkgs) {
    tryCatch({
      if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        message(sprintf("Installing %s ...", pkg))
        installer(pkg, ...)
      }
    }, error = function(e) {
      message(sprintf("Failed to install %s: %s", pkg, e$message))
    })
  }
}

# Base tooling
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes", repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

# Bioconductor 3.21 (for R 4.5 + anndataR support)
BiocManager::install(version = "3.21", ask = FALSE)

# ---- CRAN Core Packages ----

# Essential tidyverse ecosystem
essential_tidyverse <- c(
  "tidyverse",      # dplyr, ggplot2, tidyr, etc.
  "Matrix",         # Sparse matrices (Seurat dependency)
  "future",         # Parallelization
  "future.apply",
  "parallelly"
)

# Visualization & reporting
visualization <- c(
  "ggpubr",         # Publication-ready plots
  "ggridges",       # Ridge plots
  "patchwork",      # Combine plots
  "pheatmap",       # Heatmaps
  "RColorBrewer",   # Color palettes
  "Cairo",          # Graphics device
  "ragg",           # Anti-aliased graphics
  "kableExtra",     # Tables
  "plotly",         # Interactive plots
  "ggrastr"         # Rasterize ggplot layers
)

# Statistical modeling
stats_modeling <- c(
  "lme4",           # Mixed models
  "brms",           # Bayesian regression
  "broom.mixed",    # Tidy model outputs
  "DescTools"       # Descriptive statistics
)

# Data manipulation
data_manip <- c(
  "reshape2",       # Data reshaping
  "Rcpp",           # C++ integration
  "RcppEigen"       # Linear algebra
)

# VS Code integration
vscode_support <- c(
  "languageserver", # LSP for VS Code
  "httpgd"          # Graphics device (installed separately but declaring here)
)

# R-Python interoperability
r_python <- c(
  "reticulate"      # R <-> Python bridge
)

# Combine CRAN packages
cran_core <- c(
  essential_tidyverse,
  visualization,
  stats_modeling,
  data_manip,
  vscode_support,
  r_python,
  "devtools",       # Development tools
  "hdf5r"           # HDF5 file format
)

safe_install(cran_core, install.packages, repos = "https://cloud.r-project.org")

# ---- Bioconductor Core Packages ----

bioc_core <- c(
  # Single-cell infrastructure
  "SingleCellExperiment",
  "scran",
  "scater",
  "scuttle",

  # RNA-seq foundations
  "edgeR",
  "limma",
  "DESeq2",

  # Genomics infrastructure
  "GenomicRanges",
  "GenomeInfoDb",
  "IRanges",
  "S4Vectors",
  "SummarizedExperiment",
  "BiocParallel",

  # QC and utilities
  "SoupX",
  "DropletUtils",
  "scDblFinder",

  # Visualization
  "ComplexHeatmap",
  "dittoSeq",

  # Annotation
  "biomaRt",

  # HDF5 and data structures
  "HDF5Array",
  "DelayedArray",
  "DelayedMatrixStats"
)

safe_install(bioc_core, BiocManager::install, ask = FALSE, update = FALSE)

# ---- Seurat Ecosystem ----

# Ensure r-universe repositories
setRepositories(ind = 1:3, addURLs = c(
  satijalab = "https://satijalab.r-universe.dev",
  bnprks    = "https://bnprks.r-universe.dev/"
))

seurat_packages <- c(
  "Seurat",
  "BPCells",        # Efficient on-disk matrices
  "presto",         # Fast Wilcoxon
  "glmGamPoi",      # GLM for Seurat
  "Signac",         # scATAC-seq
  "sctransform"     # Variance stabilization
)

safe_install(seurat_packages, install.packages)

# ---- GSEA & Pathway Analysis ----

gsea_packages <- c(
  "clusterProfiler",
  "GSVA",
  "fgsea",
  "msigdbr",
  "enrichplot"
)

safe_install(gsea_packages, BiocManager::install, ask = FALSE, update = FALSE)

# decoupleR and PROGENy/DoRothEA (Bioc + GitHub)
safe_install("decoupleR", BiocManager::install, ask = FALSE, update = FALSE)

# ---- Multi-factorial & Latent Embeddings ----

multifactorial <- c(
  "muscat",         # Multi-sample multi-group scRNA-seq
  "harmony",        # Batch correction
  "mbkmeans"        # Mini-batch k-means
)

safe_install(multifactorial, BiocManager::install, ask = FALSE, update = FALSE)

# WGCNA from CRAN
safe_install("WGCNA", install.packages, repos = "https://cloud.r-project.org")

# ---- GitHub Packages (Seurat ecosystem & key tools) ----

github_packages <- c(
  # Seurat ecosystem
  "satijalab/seurat-data",
  "satijalab/azimuth",
  "mojaveazure/seurat-disk",

  # Interoperability
  "pmbio/MuDataSeurat",    # R <-> Python muon/anndata
  "cellgeni/sceasy",       # Format conversion
  "zellkonverter/zellkonverter",  # Also try Bioc version first

  # NMF and factorization
  "carmonalab/GeneNMF",

  # LIGER (multi-modal)
  "welch-lab/liger"
)

for (pkg in github_packages) {
  try({
    if (pkg == "zellkonverter/zellkonverter") {
      # Try Bioc first
      if (!requireNamespace("zellkonverter", quietly = TRUE)) {
        BiocManager::install("zellkonverter", ask = FALSE, update = FALSE)
      }
    } else {
      remotes::install_github(pkg, quiet = TRUE, upgrade = "never")
    }
  }, silent = TRUE)
}

# ---- MOFA2, mixOmics, lemur (multi-modal) ----
safe_install("MOFA2", BiocManager::install, ask = FALSE, update = FALSE)
safe_install("mixOmics", BiocManager::install, ask = FALSE, update = FALSE)
safe_install("lemur", BiocManager::install, ask = FALSE, update = FALSE)

# ---- anndataR (requires R 4.5!) ----
# Try CRAN/Bioc first, then GitHub
if (!requireNamespace("anndataR", quietly = TRUE)) {
  try({
    BiocManager::install("anndataR", ask = FALSE, update = FALSE)
  }, silent = TRUE)

  if (!requireNamespace("anndataR", quietly = TRUE)) {
    try(remotes::install_github("scverse/anndataR"), silent = TRUE)
  }
}

# ---- PROGENy and DoRothEA (if not in decoupleR) ----
try({
  if (!requireNamespace("progeny", quietly = TRUE)) {
    BiocManager::install("progeny", ask = FALSE, update = FALSE)
  }
}, silent = TRUE)

try({
  if (!requireNamespace("dorothea", quietly = TRUE)) {
    BiocManager::install("dorothea", ask = FALSE, update = FALSE)
  }
}, silent = TRUE)

# ---- Log installed packages ----
ip <- as.data.frame(installed.packages()[, c("Package", "Version", "Built")], stringsAsFactors = FALSE)
log_dir <- "/opt/settings"
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
write.csv(ip, file.path(log_dir, "installed_R_core_packages.csv"), row.names = FALSE)

message("Core R package installation completed.")
message(sprintf("Total packages installed: %d", nrow(ip)))
