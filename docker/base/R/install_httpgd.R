#!/usr/bin/env Rscript

# Install httpgd with CRAN-first, GitHub fallback
if (!requireNamespace("httpgd", quietly = TRUE)) {
  try(install.packages("httpgd", repos = "https://cloud.r-project.org"), silent = TRUE)
}

if (!requireNamespace("httpgd", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org")
  }
  try(remotes::install_github("nx10/httpgd"), silent = TRUE)
}

invisible(TRUE)

