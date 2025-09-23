#!/usr/bin/env Rscript

# Dedicated ArchR installation into ARCHR_LIB

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Determine ArchR library path
archr_lib <- Sys.getenv("ARCHR_LIB", unset = file.path(Sys.getenv("HOME"), "R", "archr-lib"))
if (!dir.exists(archr_lib)) dir.create(archr_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(archr_lib, .libPaths()))

# Base installers
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Pin a Bioconductor release for reproducibility
BiocManager::install(version = "3.20", ask = FALSE)

# Install ArchR and extras
devtools::install_github("GreenleafLab/ArchR", ref = "master", repos = BiocManager::repositories())
devtools::install_github("GreenleafLab/chromVARmotifs", ref = "master", repos = BiocManager::repositories())

# Optional extras; do not fail build if these fail
if (requireNamespace("ArchR", quietly = TRUE)) {
  try(ArchR::installExtraPackages(), silent = TRUE)
}

cat("ArchR installation completed in:", archr_lib, "\n")
