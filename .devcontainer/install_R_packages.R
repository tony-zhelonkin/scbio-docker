#!/usr/bin/env Rscript

# Reproducible CRAN snapshot via Posit Package Manager
snapshot <- Sys.getenv("RSPM_SNAPSHOT", unset = "2025-02-15")
options(repos = c(CRAN = paste0("https://packagemanager.posit.co/cran/__linux__/jammy/", snapshot)))

# ---- GitHub token passthrough (optional) ----
github_pat <- Sys.getenv("GITHUB_PAT")
if (nzchar(github_pat)) Sys.setenv(GITHUB_PAT = github_pat)

# ---- helpers ----
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

install_cran_version <- function(pkg, version, repos = "https://cloud.r-project.org") {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    remotes::install_version(pkg, version = version, repos = repos, upgrade = "never")
  }
}

# ---- base tooling ----
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes", repos = "https://cloud.r-project.org")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://cloud.r-project.org")

# Stick to a Bioc release (pins major/minor versions reasonably)
BiocManager::install(version = "3.20", ask = FALSE)

# ---- CRAN packages ----
cran_packages <- c(
  "future","future.apply","textshaping","pandoc","lattice","ragg","NMF","IRkernel","Rfast",
  "harmony","kableExtra","plotly","shiny","quarto","pheatmap","devtools","RColorBrewer",
  "caTools","parallel","RSpectra","irlba","SCINA","DescTools","lme4","reshape2","ggridges",
  "mbkmeans","mice","broom.mixed","remotes","Rcpp","RcppEigen","languageserver","hdf5r",
  "httpgd","ggrastr","networkD3","r2d3","Matrix","tidyverse","ggpubr","Cairo","imager",
  "lightgbm","rliger","splines","sleepwalk","singleCellHaystack","ClusterR","DDRTree",
  "densityClust","stringi","WGCNA","msigdbr",
  # ---- your requested additions ----
  "RcppML","GeneNMF","aricode","cluster","FNN"
)
safe_install(cran_packages, install.packages, repos = "https://cloud.r-project.org")

# Some packages you wanted from specific repos
try(remotes::install_github("immunogenomics/lisi"), silent = TRUE)
try(install.packages("RcppPlanc", repos = "https://welch-lab.r-universe.dev"), silent = TRUE)

# Locfit source (as in your original)
try(install.packages('https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz',
                     repos = NULL, type = 'source'), silent = TRUE)

# ---- Bioconductor packages ----
bioc_packages <- c(
  "scran","txdbmaker","dittoSeq","impute","SingleR","celldex","preprocessCore",
  "GenomicRanges","GenomeInfoDb","DESeq2","Rsamtools","S4Vectors","IRanges",
  "BiocParallel","DelayedArray","biovizBase","SoupX","scater","scDblFinder","scry",
  "zellkonverter","SingleCellExperiment","ComplexHeatmap","tidySingleCellExperiment",
  "BiocGenerics","DelayedMatrixStats","limma","SummarizedExperiment","batchelor",
  "HDF5Array","terra","Gviz","rtracklayer","chromVAR","scmap","DOSE","pathview",
  "clusterProfiler","AnnotationHub","biomaRt","ensembldb",
  "BSgenome.Hsapiens.UCSC.hg19","BSgenome.Hsapiens.UCSC.hg38",
  "BSgenome.Mmusculus.UCSC.mm10","BSgenome.Mmusculus.UCSC.mm39",
  "clusterExperiment","msigdb","EnsDb.Hsapiens.v75","EnsDb.Hsapiens.v79","EnsDb.Hsapiens.v86",
  "EnsDb.Mmusculus.v75","EnsDb.Mmusculus.v79","org.Hs.eg.db","org.Mm.eg.db",
  "DropletUtils",
  "JASPAR2022","JASPAR2024","TFBSTools","motifmatchr","scTensor",
  "SingleCellSignalR","slingshot","sctransform","splatter","sva",
  "UCell","mixOmics","MOFA2","lemur","SingleCellMultiModal"
)
safe_install(bioc_packages, BiocManager::install, ask = FALSE, update = FALSE)

# ---- Seurat ecosystem ----
setRepositories(ind = 1:3, addURLs = c("https://satijalab.r-universe.dev", "https://bnprks.r-universe.dev/"))
seurat_packages <- c("Seurat","BPCells","presto","glmGamPoi","Signac")
safe_install(seurat_packages, install.packages)

# ---- GitHub extras (non-ArchR) ----
github_packages <- c(
  "renozao/xbioc",

  # Seurat ecosystem
  "satijalab/seurat-data", "satijalab/azimuth", "satijalab/azimuth",
  "mojaveazure/seurat-disk",

  # Epigenomics tools (ArchR is installed in ArchR image variant)

  # QC tools
  "powellgenomicslab/DropletQC",
  "chris-mcginnis-ucsf/DoubletFinder",

  "settylab/convert2anndata",
  "cole-trapnell-lab/monocle3",
  "meichendong/SCDC",
  "cole-trapnell-lab/cicero-release",
  "aertslab/RcisTarget","aertslab/AUCell","aertslab/cisTopic",
  "SydneyBioX/scClustBench","SydneyBioX/scDC",
  "immunogenomics/SCENT",
  "jokergoo/circlize",
  "pcahan1/singleCellNet",
  "jinworks/CellChat",
  "carmonalab/SignatuR",
  "Zhen-Miao/PICsnATAC",
  "Zhen-Miao/PACS",
  "quadbio/Pando",
  "buenrostrolab/FigR"
)

for (pkg in github_packages) {
  try({
    if (pkg == "powellgenomicslab/DropletQC") {
      devtools::install_github("powellgenomicslab/DropletQC", build_vignettes = FALSE)
    } else if (pkg == "cole-trapnell-lab/cicero-release") {
      devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3", quiet = TRUE)
    } else if (pkg == "mojaveazure/loomR@develop") {
      devtools::install_github("mojaveazure/loomR", ref = "develop", quiet = TRUE)
    } else {
      remotes::install_github(pkg, quiet = TRUE)
    }
  }, silent = TRUE)
}

# ArchR extras are handled in the ArchR image variant

# Log all installed packages (exact versions)
ip <- as.data.frame(installed.packages()[, c("Package","Version","Built")], stringsAsFactors = FALSE)
log_dir <- "/opt/settings"
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
write.csv(ip, file.path(log_dir, "installed_R_packages.csv"), row.names = FALSE)

# ---- Log GitHub remotes (commit SHAs if available) ----
# Some packages have no Remote* fields; never error on that.
log_dir <- "/opt/settings"
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

ip <- as.data.frame(installed.packages()[, c("Package","Version","Built")],
                    stringsAsFactors = FALSE)
write.csv(ip, file.path(log_dir, "installed_R_packages.csv"), row.names = FALSE)

fields <- c("RemoteType","RemoteRepo","RemoteUsername","RemoteRef","RemoteSha")

get_remote <- function(p) {
  d <- tryCatch(utils::packageDescription(p), error = function(e) NULL)
  if (is.null(d)) return(NULL)
  # always produce all columns; fill with NA when missing
  vals <- vapply(fields,
                 function(f) {
                   x <- d[[f]]
                   if (is.null(x)) NA_character_ else as.character(x)
                 },
                 FUN.VALUE = character(1),
                 USE.NAMES = TRUE)
  # one-row data.frame, ordered columns
  out <- as.data.frame(as.list(vals), stringsAsFactors = FALSE, check.names = FALSE)
  out$Package <- p
  out[c("Package", fields)]
}

gh_list <- lapply(ip$Package, get_remote)
gh_list <- Filter(Negate(is.null), gh_list)

# Write an empty CSV with headers if nothing has Remote* metadata
gh_df <- if (length(gh_list)) {
  do.call(rbind, gh_list)
} else {
  data.frame(Package = character(),
             RemoteType = character(),
             RemoteRepo = character(),
             RemoteUsername = character(),
             RemoteRef = character(),
             RemoteSha = character(),
             check.names = FALSE)
}

# Never let logging break the build
try(write.csv(gh_df, file.path(log_dir, "installed_R_github_remotes.csv"), row.names = FALSE),
    silent = TRUE)

message("R package installation completed and logs written.")