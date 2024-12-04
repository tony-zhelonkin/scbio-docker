#!/usr/bin/env Rscript

# GitHub Api rate limit is likely to hit. 
# usethis::create_github_token() # to create a token on GitHub, follow prompts 
# gitcreds::gitcreds_set() # to add the token to session, follow prompts 
github_pat <- Sys.getenv("GITHUB_PAT")
if (github_pat != "") {
  Sys.setenv(GITHUB_PAT = github_pat)
}
# docker build --build-arg GITHUB_PAT=${GITHUB_PAT} . # command arg for build

# Function to safely install packages
safe_install <- function(packages, installer, ...) {
  for(package in packages) {
    tryCatch({
      if (!require(package, character.only = TRUE, quietly = TRUE)) {
        message(sprintf("Installing %s...", package))
        installer(package, ...)
      }
    }, error = function(e) {
      message(sprintf("Failed to install %s: %s", package, e$message))
    })
  }
}

# # Check R version (Deprecate version check)
# if (getRversion() < "4.4.2") {
#   stop("This script requires R version 4.4.2 or higher")
# }

# Install BiocManager if not present
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}

# Set Bioconductor version
BiocManager::install(version = "3.20", ask = FALSE)
message(sprintf("Using Bioconductor version: %s", BiocManager::version()))

# Base CRAN packages
message("Installing CRAN packages...")
cran_packages <- c(
  "future",
  "future.apply",
  "textshaping",
  "pandoc",
  "lattice",
  "ragg",
  "NMF",
  "harmony",
  "kableExtra",
  "plotly",
  "shiny",
  "quarto",
  "pheatmap",
  "devtools",
  "RColorBrewer",
  "caTools",
  "parallel",
  "RSpectra",
  "irlba",
  "DescTools",
  "lme4",
  "reshape2",
  "ggridges",
  "mice",
  "broom.mixed",
  "remotes",
  "Rcpp",
  "RcppEigen",
  "languageserver",
  "hdf5r",
  "httpgd",
  "ggrastr",
  "networkD3",
  "r2d3",
  "Matrix",
  "tidyverse",
  "ggpubr",
  "Cairo",
  "imager",
  "lightgbm",
  "rliger",
  "splines",
  "reshape2",
  "sleepwalk",
  "singleCellHaystack",
  "ClusterR",
  "DDRTree",
  "densityClust",
  "stringi",
  "WGCNA"
)

# Install CRAN packages
safe_install(cran_packages, install.packages, repos = "https://cloud.r-project.org/")


# Install source
install.packages('https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz', repos=NULL, type='source')


# Bioconductor packages
message("Installing Bioconductor packages...")
bioc_packages <- c(
  # Core Bioconductor packages
  "scran",
  "txdbmaker",
  "impute",
  "preprocessCore",
  "GenomicRanges",
  "GenomeInfoDb",
  "DESeq2",
  "Rsamtools",
  "S4Vectors",
  "IRanges",
  "BiocParallel",
  "DelayedArray",
  "biovizBase",
  "SoupX",
  "scater",
  "scDblFinder",
  "scry",
  "zellkonverter",
  "SingleCellExperiment",
  "ComplexHeatmap",
  "tidySingleCellExperiment",
  "BiocGenerics",
  "DelayedMatrixStats",
  "limma",
  "SummarizedExperiment",
  "batchelor",
  "HDF5Array",
  "terra",
  "Gviz",
  "rtracklayer",
  "chromVAR",
  "scmap",
  "DOSE",
  "pathview",
  "clusterProfiler",
  # Additional packages
  "AnnotationHub",
  "biomaRt",
  "ensembldb",
  "BSgenome.Hsapiens.UCSC.hg19",
  "BSgenome.Hsapiens.UCSC.hg38",
  "BSgenome.Mmusculus.UCSC.mm10",
  "BSgenome.Mmusculus.UCSC.mm39",
  "clusterExperiment",
  "msigdb",
  "EnsDb.Hsapiens.v75",
  "EnsDb.Hsapiens.v79",
  "EnsDb.Hsapiens.v86",
  "EnsDb.Mmusculus.v75",
  "EnsDb.Mmusculus.v79",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "DropletUtils",
  "JASPAR2022",
  "TFBSTools",
  "motifmatchr",
  "scTensor",
  "SingleCellSignalR",
  "slingshot",
  "sctransform",
  "splatter",
  "sva"
)

# Install Bioconductor packages
safe_install(bioc_packages, BiocManager::install, ask = FALSE)

# Set up Seurat repositories
message("Installing Seurat and related packages...")
setRepositories(ind = 1:3, addURLs = c("https://satijalab.r-universe.dev", "https://bnprks.r-universe.dev/"))

# Seurat packages
seurat_packages <- c(
  "Seurat",
  "BPCells",
  "presto",
  "glmGamPoi",
  "Signac"
)
safe_install(seurat_packages, install.packages)

# Install RcppPlanc from specific repository
install.packages("RcppPlanc", repos = "https://welch-lab.r-universe.dev")

# GitHub packages
message("Installing GitHub packages...")
github_packages <- c(
  "renozao/xbioc",
  "satijalab/seurat-data",
  "satijalab/azimuth",
  "satijalab/seurat-wrappers",
  "GreenleafLab/ArchR",
  "powellgenomicslab/DropletQC",
  "chris-mcginnis-ucsf/DoubletFinder",
  "cole-trapnell-lab/monocle3",
  #"velocyto-team/velocyto.R", # doesn`t compile, should velocity analysis be performed I`ll stick to scvelo python
  "meichendong/SCDC",
  "GreenleafLab/chromVARmotifs",
  "cole-trapnell-lab/cicero-release",
  "aertslab/RcisTarget",
  "aertslab/AUCell",
  "SydneyBioX/scClustBench",
  "SydneyBioX/scDC",
  "aertslab/cisTopic",
  "immunogenomics/SCENT",
  "mojaveazure/loomR@develop",
  "jokergoo/circlize",
  "pcahan1/singleCellNet",
  "jinworks/CellChat",
  "carmonalab/SignatuR"
)

# Install GitHub packages with special handling for certain packages
for (package in github_packages) {
  tryCatch({
    message(sprintf("Installing %s from GitHub...", package))
    if (package == "powellgenomicslab/DropletQC") {
      devtools::install_github("powellgenomicslab/DropletQC", build_vignettes = FALSE)
    } else if (package == "cole-trapnell-lab/cicero-release") {
      devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3", quiet = TRUE)
    } else if (package == "mojaveazure/loomR@develop") {
      devtools::install_github("mojaveazure/loomR", ref = "develop", quiet = TRUE)
    #} else if (package == "velocyto-team/velocyto.R") {
    #  devtools::install_github("velocyto-team/velocyto.R", quiet = TRUE)
    } else {
      remotes::install_github(package, quiet = TRUE)
    }
  }, error = function(e) {
    message(sprintf("Failed to install %s: %s", package, e$message))
  })
}

# Install ArchR extra packages if ArchR was successfully installed
if (require("ArchR")) {
  message("Installing ArchR extra packages...")
  tryCatch({
    ArchR::installExtraPackages()
  }, error = function(e) {
    message("Failed to install ArchR extra packages: ", e$message)
  })
}

message("R package installation completed!")