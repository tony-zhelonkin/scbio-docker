#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

installed <- function(pkg) requireNamespace(pkg, quietly = TRUE)

if (!installed("httpgd")) {
  try(install.packages("httpgd"), silent = TRUE)
}

if (!installed("httpgd")) {
  if (!installed("remotes")) install.packages("remotes")
  Sys.unsetenv("GITHUB_PAT")
  remotes::install_github("nx10/httpgd", upgrade = "never", quiet = TRUE)
}

if (installed("httpgd")) {
  cat("httpgd version:", as.character(utils::packageVersion("httpgd")), "\n")
} else {
  cat("httpgd not installed. You can install via CRAN or GitHub.\n")
}


