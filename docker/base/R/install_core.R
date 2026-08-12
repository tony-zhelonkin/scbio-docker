#!/usr/bin/env Rscript

# CORE R PACKAGES for scbio-docker (moved under docker/base/R)
# These are pre-installed in the image. Additional packages can be installed at runtime.

snapshot <- Sys.getenv("RSPM_SNAPSHOT", unset = "2026-04-15")
options(repos = c(CRAN = paste0("https://packagemanager.posit.co/cran/__linux__/jammy/", snapshot)))

github_pat <- Sys.getenv("GITHUB_PAT")
if (nzchar(github_pat)) Sys.setenv(GITHUB_PAT = github_pat)

# Failures are recorded to a file, not just warned about: R warnings are
# deferred and were being lost in the build log, which hid chromVAR and six
# other absent packages across two releases.
.failures <- new.env(parent = emptyenv())
.failures$rows <- list()

record_failure <- function(pkg, msg) {
  .failures$rows[[length(.failures$rows) + 1L]] <-
    data.frame(package = pkg, error = gsub("[\r\n]+", " ", msg),
               stringsAsFactors = FALSE)
  message(sprintf("FAILED: %s: %s", pkg, msg))
}

write_failure_report <- function(dir = "/opt/settings") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  path <- file.path(dir, "install_failures.csv")
  if (length(.failures$rows)) {
    df <- do.call(rbind, .failures$rows)
    write.csv(df, path, row.names = FALSE)
    message(sprintf("\n=== %d PACKAGE(S) FAILED TO INSTALL ===", nrow(df)))
    for (i in seq_len(nrow(df))) message(sprintf("  %s: %s", df$package[i], df$error[i]))
    message(sprintf("=== report: %s ===\n", path))
  } else {
    write.csv(data.frame(package = character(), error = character()), path, row.names = FALSE)
    message("All requested packages installed; no failures.")
  }
}

safe_install <- function(pkgs, installer, ...) {
  installed <- installed.packages()[, "Package"]
  for (pkg in pkgs) {
    if (!(pkg %in% installed)) {
      message(sprintf("Installing %s ...", pkg))
      tryCatch({
        installer(pkg, ...)
      }, error = function(e) {
        record_failure(pkg, conditionMessage(e))
      })
      # install.packages() signals a warning rather than an error when a build
      # fails, so the tryCatch above never fires; confirm the package landed.
      if (!requireNamespace(pkg, quietly = TRUE)) {
        record_failure(pkg, "install returned without error but package is not loadable")
      }
    } else {
      message(sprintf("Package %s already installed, skipping", pkg))
    }
  }
}

if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes", repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install(version = "3.22", ask = FALSE)

essential_tidyverse <- c(
  "tidyverse","Matrix","data.table","future","future.apply","parallelly"
)

visualization <- c(
  "ggpubr","ggridges","patchwork","pheatmap","RColorBrewer","Cairo","textshaping","ragg","kableExtra","plotly","ggrastr"
)

stats_modeling <- c("lme4","brms","broom.mixed","DescTools","Rfast")
data_manip <- c("reshape2","Rcpp","RcppEigen")
vscode_support <- c("languageserver","httpgd")
r_python <- c("reticulate")
notebooks <- c("IRkernel")

cran_core <- c(
  essential_tidyverse, visualization, stats_modeling, data_manip, vscode_support, r_python, notebooks,
  "devtools","hdf5r","pandoc",
  "tictoc"  # PACS dependency
)
safe_install(cran_core, install.packages, repos = "https://cloud.r-project.org")

bioc_core <- c(
  "SingleCellExperiment","scran","scater","scuttle",
  "edgeR","limma","DESeq2",
  "GenomicRanges","GenomeInfoDb","IRanges","S4Vectors","SummarizedExperiment","BiocParallel",
  "Rsamtools",
  "SoupX","DropletUtils","scDblFinder",
  "ComplexHeatmap","dittoSeq",
  "biomaRt","AnnotationHub","AnnotationDbi",
  "HDF5Array","DelayedArray","DelayedMatrixStats"
)
safe_install(bioc_core, BiocManager::install, ask = FALSE, update = FALSE)

setRepositories(ind = 1:3, addURLs = c(
  satijalab = "https://satijalab.r-universe.dev",
  bnprks    = "https://bnprks.r-universe.dev/"
))
seurat_packages <- c("Seurat","BPCells","presto","glmGamPoi","Signac","sctransform")
safe_install(seurat_packages, install.packages)

gsea_packages <- c("clusterProfiler","GSVA","fgsea","msigdbr","enrichplot")
safe_install(gsea_packages, BiocManager::install, ask = FALSE, update = FALSE)
safe_install("decoupleR", BiocManager::install, ask = FALSE, update = FALSE)

chromatin_packages <- c("chromVAR","motifmatchr","TFBSTools","JASPAR2022","JASPAR2024","SingleR","celldex")
safe_install(chromatin_packages, BiocManager::install, ask = FALSE, update = FALSE)

# Ortholog conversion (offline mappings; biomaRt above is the online path)
ortholog_packages <- c("homologene","babelgene")
safe_install(ortholog_packages, install.packages, repos = "https://cloud.r-project.org")
safe_install("orthogene", BiocManager::install, ask = FALSE, update = FALSE)

# Multivariate exploratory analysis (PCA/MCA/MFA + ggplot2 viz)
mva_packages <- c("FactoMineR","factoextra")
safe_install(mva_packages, install.packages, repos = "https://cloud.r-project.org")

organism_packages <- c("EnsDb.Mmusculus.v79","org.Hs.eg.db","org.Mm.eg.db")
safe_install(organism_packages, BiocManager::install, ask = FALSE, update = FALSE)

# Azimuth's dependency chain, declared explicitly rather than pulled in as an
# invisible side effect. ~700MB of human annotation; see AGENTS.md, which
# documents this as the exception to the "no heavy annotation packages" rule.
azimuth_deps_bioc <- c("BSgenome.Hsapiens.UCSC.hg38","EnsDb.Hsapiens.v86","JASPAR2020")
safe_install(azimuth_deps_bioc, BiocManager::install, ask = FALSE, update = FALSE)
azimuth_deps_cran <- c("shinyBS","shinydashboard","shinyjs")
safe_install(azimuth_deps_cran, install.packages, repos = "https://cloud.r-project.org")

multifactorial <- c("muscat","harmony","mbkmeans")
safe_install(multifactorial, BiocManager::install, ask = FALSE, update = FALSE)
# Via BiocManager: deps impute/preprocessCore are Bioconductor-only.
safe_install("WGCNA", BiocManager::install, ask = FALSE, update = FALSE)

github_packages <- c(
  # seurat-disk precedes azimuth: it is an Azimuth dependency.
  "satijalab/seurat-data","mojaveazure/seurat-disk","satijalab/azimuth",
  "pmbio/MuDataSeurat","cellgeni/sceasy","zellkonverter/zellkonverter",
  "carmonalab/GeneNMF","immunogenomics/crescendo",
  "Zhen-Miao/PICsnATAC","Zhen-Miao/PACS",
  "GreenleafLab/chromVARmotifs"
)
# The repo name is not always the package name (satijalab/seurat-data ->
# SeuratData), so a derived name silently breaks both the "already installed"
# check and the failure report.
GH_PKG_NAME <- c(
  "satijalab/seurat-data"   = "SeuratData",
  "satijalab/azimuth"       = "Azimuth",
  "mojaveazure/seurat-disk" = "SeuratDisk"
)

install_gh_pkg <- function(slug) {
  pkg_name <- if (slug %in% names(GH_PKG_NAME)) GH_PKG_NAME[[slug]] else sub(".*/", "", slug)
  if (requireNamespace(pkg_name, quietly = TRUE)) {
    message(sprintf("Package %s already installed, skipping", pkg_name))
    return(invisible(TRUE))
  }
  # api.github.com intermittently drops HTTP/2 streams mid-transfer; a plain
  # single attempt loses several packages per build. Retry before falling back.
  ok <- FALSE
  for (attempt in 1:3) {
    ok <- tryCatch({
      remotes::install_github(slug, quiet = TRUE, upgrade = "never")
      TRUE
    }, error = function(e) {
      message(sprintf("install_github attempt %d/3 failed for %s: %s",
                      attempt, slug, e$message))
      FALSE
    })
    if (ok) break
    Sys.sleep(5 * attempt)
  }
  if (!ok) {
    # Fallback: manual git clone + install_local (avoids api.github.com HTTP/2 flakiness).
    message(sprintf("Falling back to git clone + install_local for %s", slug))
    tryCatch({
      tmp <- tempfile()
      dir.create(tmp)
      url <- sprintf("https://github.com/%s.git", slug)
      Sys.setenv(GIT_HTTP_VERSION = "HTTP/1.1")
      # install_local still reaches api.github.com for remote deps; force
      # libcurl to HTTP/1.1 there as well.
      Sys.setenv(R_LIBCURL_HTTP_VERSION = "1.1")
      system2("git", c("-c", "http.version=HTTP/1.1", "clone", "--depth", "1",
                       url, file.path(tmp, pkg_name)),
              stdout = TRUE, stderr = TRUE)
      remotes::install_local(file.path(tmp, pkg_name),
                             quiet = TRUE, upgrade = "never", dependencies = TRUE)
      unlink(tmp, recursive = TRUE)
    }, error = function(e) message(sprintf("clone+install_local failed for %s: %s",
                                           pkg_name, e$message)))
  }
  # Last resort: git transport with Remotes resolution off. api.github.com
  # drops HTTP/2 streams from this host, and both install_github and
  # install_local go through it to resolve DESCRIPTION `Remotes:`.
  if (!requireNamespace(pkg_name, quietly = TRUE)) {
    message(sprintf("Trying install_git(dependencies=FALSE) for %s", slug))
    try(remotes::install_git(sprintf("https://github.com/%s", slug),
                             quiet = TRUE, upgrade = "never",
                             dependencies = FALSE), silent = TRUE)
  }
  if (!requireNamespace(pkg_name, quietly = TRUE))
    record_failure(pkg_name, sprintf("github install of %s produced no loadable package", slug))
}

for (pkg in github_packages) {
  if (pkg == "zellkonverter/zellkonverter") {
    if (!requireNamespace("zellkonverter", quietly = TRUE)) {
      try(BiocManager::install("zellkonverter", ask = FALSE, update = FALSE), silent = TRUE)
    }
  } else {
    install_gh_pkg(pkg)
  }
}

# rliger (NOT "liger": the GitHub slug welch-lab/liger yields the wrong package
# name, so it was never detected as installed). RcppPlanc lives on r-universe.
safe_install("rliger", install.packages,
             repos = c("https://welch-lab.r-universe.dev",
                       "https://cloud.r-project.org"))

safe_install("MOFA2", BiocManager::install, ask = FALSE, update = FALSE)
safe_install("mixOmics", BiocManager::install, ask = FALSE, update = FALSE)
safe_install("lemur", BiocManager::install, ask = FALSE, update = FALSE)

if (!requireNamespace("anndataR", quietly = TRUE)) {
  try({ BiocManager::install("anndataR", ask = FALSE, update = FALSE) }, silent = TRUE)
  if (!requireNamespace("anndataR", quietly = TRUE)) {
    try(remotes::install_github("scverse/anndataR"), silent = TRUE)
  }
}

try({ if (!requireNamespace("progeny", quietly = TRUE)) BiocManager::install("progeny", ask = FALSE, update = FALSE) }, silent = TRUE)
try({ if (!requireNamespace("dorothea", quietly = TRUE)) BiocManager::install("dorothea", ask = FALSE, update = FALSE) }, silent = TRUE)

ip <- as.data.frame(installed.packages()[, c("Package", "Version", "Built")], stringsAsFactors = FALSE)
log_dir <- "/opt/settings"
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
write.csv(ip, file.path(log_dir, "installed_R_core_packages.csv"), row.names = FALSE)

write_failure_report()

message("Core R package installation completed.")
message(sprintf("Total packages installed: %d", nrow(ip)))

