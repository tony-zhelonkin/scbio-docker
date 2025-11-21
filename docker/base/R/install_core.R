#!/usr/bin/env Rscript

# CORE R PACKAGES for scbio-docker (moved under docker/base/R)
# These are pre-installed in the image. Additional packages can be installed at runtime.

snapshot <- Sys.getenv("RSPM_SNAPSHOT", unset = "2025-02-15")
options(repos = c(CRAN = paste0("https://packagemanager.posit.co/cran/__linux__/jammy/", snapshot)))

github_pat <- Sys.getenv("GITHUB_PAT")
if (nzchar(github_pat)) Sys.setenv(GITHUB_PAT = github_pat)

safe_install <- function(pkgs, installer, ...) {
  installed <- installed.packages()[, "Package"]
  for (pkg in pkgs) {
    if (!(pkg %in% installed)) {
      message(sprintf("Installing %s ...", pkg))
      tryCatch({
        installer(pkg, ...)
      }, error = function(e) {
        warning(sprintf("Failed to install %s: %s", pkg, e$message))
      })
    } else {
      message(sprintf("Package %s already installed, skipping", pkg))
    }
  }
}

if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes", repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

BiocManager::install(version = "3.21", ask = FALSE)

essential_tidyverse <- c(
  "tidyverse","Matrix","future","future.apply","parallelly"
)

visualization <- c(
  "ggpubr","ggridges","patchwork","pheatmap","RColorBrewer","Cairo","ragg","kableExtra","plotly","ggrastr"
)

stats_modeling <- c("lme4","brms","broom.mixed","DescTools")
data_manip <- c("reshape2","Rcpp","RcppEigen")
vscode_support <- c("languageserver","httpgd")
r_python <- c("reticulate")

cran_core <- c(
  essential_tidyverse, visualization, stats_modeling, data_manip, vscode_support, r_python,
  "devtools","hdf5r","pandoc"
)
safe_install(cran_core, install.packages, repos = "https://cloud.r-project.org")

bioc_core <- c(
  "SingleCellExperiment","scran","scater","scuttle",
  "edgeR","limma","DESeq2",
  "GenomicRanges","GenomeInfoDb","IRanges","S4Vectors","SummarizedExperiment","BiocParallel",
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

chromatin_packages <- c("chromVAR","motifmatchr","TFBSTools","JASPAR2022","SingleR","celldex")
safe_install(chromatin_packages, BiocManager::install, ask = FALSE, update = FALSE)

organism_packages <- c("EnsDb.Mmusculus.v79")
safe_install(organism_packages, BiocManager::install, ask = FALSE, update = FALSE)

multifactorial <- c("muscat","harmony","mbkmeans")
safe_install(multifactorial, BiocManager::install, ask = FALSE, update = FALSE)
safe_install("WGCNA", install.packages, repos = "https://cloud.r-project.org")

github_packages <- c(
  "satijalab/seurat-data","satijalab/azimuth","mojaveazure/seurat-disk",
  "pmbio/MuDataSeurat","cellgeni/sceasy","zellkonverter/zellkonverter",
  "carmonalab/GeneNMF","welch-lab/liger","immunogenomics/crescendo"
)
for (pkg in github_packages) {
  try({
    if (pkg == "zellkonverter/zellkonverter") {
      if (!requireNamespace("zellkonverter", quietly = TRUE)) {
        BiocManager::install("zellkonverter", ask = FALSE, update = FALSE)
      }
    } else {
      remotes::install_github(pkg, quiet = TRUE, upgrade = "never")
    }
  }, silent = TRUE)
}

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

message("Core R package installation completed.")
message(sprintf("Total packages installed: %d", nrow(ip)))

