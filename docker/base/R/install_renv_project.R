#!/usr/bin/env Rscript

options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

project_dir <- "/opt/settings"
lockfile <- file.path(project_dir, "renv.lock")

if (file.exists(lockfile)) {
  renv::init(project = project_dir, bare = TRUE)
  renv::restore(project = project_dir, lockfile = lockfile, prompt = FALSE)
} else {
  # Expect install script to be present as install_R_packages.R
  source(file.path(project_dir, "install_R_packages.R"))
  renv::init(project = project_dir, bare = TRUE)
  renv::snapshot(project = project_dir, lockfile = lockfile, prompt = FALSE)
}

pk <- as.data.frame(installed.packages()[, c("Package","Version","LibPath")])
write.csv(pk, file.path(project_dir, "R-packages-manifest.csv"), row.names = FALSE)

