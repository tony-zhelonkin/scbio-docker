# RNA-seq Analysis Guidelines for AI Agents

**Version:** 1.0.0
**Maintained in:** scbio-docker repository
**Last Updated:** 2025-12-12
**Purpose:** Guiding principles for AI agents performing bulk RNA-seq analysis

---

## Table of Contents

1. [Core Philosophy](#1-core-philosophy)
2. [Data Flow Principles](#2-data-flow-principles)
3. [Directory Structure](#3-directory-structure)
4. [Script Organization](#4-script-organization)
5. [Configuration Management](#5-configuration-management)
6. [Checkpoint Caching](#6-checkpoint-caching)
7. [Master Tables Pattern](#7-master-tables-pattern)
8. [Separation of Concerns](#8-separation-of-concerns)
9. [DRY Principles](#9-dry-principles)
10. [Code Syntax Standards](#10-code-syntax-standards)
11. [Visualization Standards](#11-visualization-standards)
12. [Documentation Requirements](#12-documentation-requirements)

---

## 1. Core Philosophy

### 1.1 Fundamental Principles

1. **Normalize Once, Visualize Many** - All data processing happens in Phase 1. Visualization scripts only read pre-computed results.

2. **Single Source of Truth** - Configuration lives in one place (YAML), colors defined once, schemas validated everywhere.

3. **Checkpoint Everything Expensive** - Any computation taking >1 minute must be cached. First run may take 45-60 min; subsequent runs take 5-10 min.

4. **Master Tables as Bridges** - CSV exports connect R computation to Python visualization. No language lock-in.

5. **Separation of Concerns** - Data processing scripts never plot. Visualization scripts never compute. Config files never contain logic.

6. **Reproducibility First** - Every analysis must be reproducible from checkpoints. Document what each checkpoint contains.

### 1.2 Guiding Questions Before Writing Code

- Does a checkpoint already exist for this computation?
- Is this the right phase for this code (processing vs visualization)?
- Is there a single source of truth for this parameter/color/threshold?
- Can this be done with existing toolkit functions?
- Will the output be reusable by both R and Python?

---

## 2. Data Flow Principles

### 2.1 Standard RNA-seq Pipeline Stages

```
Stage 0: Raw Data -> Counts Matrix
|-- Input: FASTQ files (external)
|-- Process: QC, alignment, quantification
+-- Output: counts matrix (genes x samples)

Stage 1: Counts -> Normalized Data
|-- Input: Raw counts + sample metadata
|-- Process: Filtering + TMM/VST normalization
+-- Output: DGEList object (checkpoint)

Stage 2: Normalized -> Differential Expression
|-- Input: Normalized counts + design formula
|-- Process: limma-voom/DESeq2 modeling
+-- Output: DE results per contrast (checkpoint)

Stage 3: DE Results -> Pathway Enrichment
|-- Input: Ranked gene list from DE
|-- Process: GSEA across multiple databases
+-- Output: Enrichment tables (checkpoint)

Stage 4: Results -> Master Tables
|-- Input: All checkpoints
|-- Process: Normalize schemas, combine
+-- Output: CSV files (language-agnostic)

Stage 5: Master Tables -> Visualization
|-- Input: CSV tables (Python) or checkpoints (R)
|-- Process: Plotting only
+-- Output: Publication figures
```

### 2.2 Data Flow Rules

1. **Immutable Inputs** - Raw data in `00_data/` is read-only. Never modify source files.

2. **Checkpoints are Intermediate** - RDS files in `checkpoints/` can be regenerated. They're cache, not archive.

3. **Master Tables are the Contract** - CSVs in `tables/` define the interface between R and Python.

4. **Outputs are Disposable** - Everything in `plots/` can be regenerated from checkpoints.

---

## 3. Directory Structure

### 3.1 Standard Layout

```
project-root/
|-- 00_data/                    # READ-ONLY input data
|   |-- processed/              # Count matrices
|   +-- references/             # Gene sets, annotations
|
|-- 01_scripts/                 # Shared tools
|   +-- RNAseq-toolkit/         # Git submodule (reusable)
|
|-- 02_analysis/                # Project-specific code
|   |-- config/                 # Configuration files
|   |   |-- pipeline.yaml       # Single source of truth
|   |   |-- config.R            # R-side config loader
|   |   +-- color_config.R      # Color definitions
|   |-- helpers/                # Project-specific utilities
|   |-- 1.x.*.R                 # Phase 1: Data processing
|   |-- 2.x.*.R                 # Phase 3: R visualizations
|   +-- 3.x.*.py                # Phase 4-5: Python/interactive
|
|-- 03_results/                 # Generated outputs
|   |-- checkpoints/            # Cached RDS objects
|   |-- tables/                 # Master CSV exports
|   |-- plots/                  # Publication figures
|   |   |-- QC/                 # Quality control
|   |   |-- Volcano/            # Volcano plots
|   |   +-- GSEA/               # Pathway visualizations
|   +-- interactive/            # HTML dashboards
|
|-- docs/
|   +-- guidelines/             # Analysis guidelines (this doc)
|
|-- CLAUDE.md                   # AI assistant context
|-- plan.md                     # Analysis strategy
|-- tasks.md                    # Execution tracker
+-- notes.md                    # Research findings
```

### 3.2 Naming Conventions

**Directories:** lowercase with underscores (`00_data`, `checkpoints`)

**Scripts:** `{phase}.{order}.{description}.{ext}`
- `1.1.core_pipeline.R` - Phase 1, step 1
- `1.5.create_master_tables.R` - Phase 1, step 5
- `2.1.visualizations.R` - Phase 3 (viz), step 1
- `3.1.pathway_explorer.py` - Phase 4-5, step 1

**Checkpoints:** `{phase}.{step}_{description}.rds`
- `1.1_dge_normalized.rds`
- `1.1_de_results.rds`
- `1.4_gsea_custom.rds`

**Master Tables:** `master_{type}.csv`
- `master_gsea_table.csv`
- `master_de_table.csv`
- `master_tf_activities.csv`

---

## 4. Script Organization

### 4.1 Phase-Based Execution Model

| Phase | Scripts | Purpose | Output |
|-------|---------|---------|--------|
| 1 | `1.x.*.R` | Core analysis | Checkpoints (RDS) |
| 2 | `1.x.*.R` | Master tables | CSV exports |
| 3 | `2.x.*.R` | R visualizations | PDF/PNG plots |
| 4-5 | `3.x.*.py` | Python/Interactive | HTML, publication figs |

### 4.2 Script Dependencies

```
1.1.core_pipeline.R ----------------------+
    |                                     |
1.2.annotate_de.R                         |
    |                                     |
1.3.custom_gsea.R                         |
    |                                     |
1.5.create_master_tables.R <--------------+
    |
+-------+
|       |
2.x.R   3.x.py
(R viz) (Python viz)
```

### 4.3 Script Template (R)

```r
#!/usr/bin/env Rscript
# {script_name} - {Brief description}
# Project: {PROJECT_ID}
# Phase: {1|2|3}
# Description: {What this script does}
# Dependencies: {List required checkpoints}
# Outputs: {List generated files}
# Updated: {YYYY-MM-DD}

# ============================================================================
# SETUP
# ============================================================================

message("=================================================================")
message("{PROJECT_ID}: {Script Title}")
message("=================================================================")

# Source configuration (ALWAYS FIRST)
source("02_analysis/config/config.R")
load_packages()
source_toolkit()

# ============================================================================
# LOAD DEPENDENCIES
# ============================================================================

# Load required checkpoints
dge <- readRDS(file.path(DIR_CHECKPOINTS, "1.1_dge_normalized.rds"))

# ============================================================================
# MAIN ANALYSIS
# ============================================================================

result <- load_or_compute(
  checkpoint_file = "1.x_result.rds",
  description = "Result description",
  compute_fn = function() {
    # Expensive computation here
  }
)

# ============================================================================
# COMPLETE
# ============================================================================

message("Script complete. Outputs saved to: ", DIR_RESULTS)
```

---

## 5. Configuration Management

### 5.1 YAML Configuration (Single Source of Truth)

**File:** `02_analysis/config/pipeline.yaml`

```yaml
# Project metadata
project:
  id: "PROJECT-ID"
  title: "Project Title"
  species: "Mus musculus"
  genome_build: "mm10"

# Colors (Okabe-Ito colorblind-safe)
colors:
  diverging:
    negative: "#2166AC"   # Blue
    neutral: "#F7F7F7"    # White
    positive: "#B35806"   # Orange
  databases:
    Hallmark: "#E69F00"
    KEGG: "#56B4E9"
    Reactome: "#009E73"
    GO_BP: "#0072B2"

# Analysis parameters
analysis:
  de_fdr_cutoff: 0.05
  de_logfc_cutoff: 2.0
  gsea_nperm: 100000
  gsea_seed: 123

# Visualization parameters
visualization:
  plot_width: 10
  plot_height: 8
  plot_dpi: 300
  n_top_pathways: 20

# Schema definitions (for validation)
schemas:
  master_gsea_table:
    required_columns:
      - pathway_id
      - pathway_name
      - database
      - nes
      - pvalue
      - padj
      - core_enrichment
```

### 5.2 R Configuration Loader

**File:** `02_analysis/config/config.R`

```r
# Load YAML config
YAML_CONFIG <- yaml::read_yaml("02_analysis/config/pipeline.yaml")

# Project metadata
PROJECT_ID <- YAML_CONFIG$project$id
SPECIES <- YAML_CONFIG$project$species

# Directories (relative to project root)
DIR_CHECKPOINTS <- "03_results/checkpoints"
DIR_TABLES <- "03_results/tables"
DIR_PLOTS <- "03_results/plots"

# Analysis parameters
DE_FDR_CUTOFF <- YAML_CONFIG$analysis$de_fdr_cutoff
GSEA_NPERM <- YAML_CONFIG$analysis$gsea_nperm

# Helper function: load_or_compute
load_or_compute <- function(checkpoint_file, compute_fn,
                            force_recompute = FALSE, description = "Result") {
  checkpoint_path <- file.path(DIR_CHECKPOINTS, checkpoint_file)

  if (file.exists(checkpoint_path) && !force_recompute) {
    message("[CACHE] Loading ", description, " from: ", checkpoint_file)
    return(readRDS(checkpoint_path))
  }

  message("[COMPUTE] Computing ", description, "...")
  start_time <- Sys.time()
  result <- compute_fn()
  elapsed <- round(difftime(Sys.time(), start_time, units = "mins"), 2)

  message("[SAVE] Saving ", description, " (took ", elapsed, " min)")
  saveRDS(result, checkpoint_path)
  return(result)
}
```

### 5.3 Configuration Rules

1. **Never hardcode** thresholds, paths, or colors in analysis scripts
2. **YAML is canonical** - R and Python both read from it
3. **Fallback defaults** - Handle missing YAML gracefully
4. **Version schemas** - Track schema versions for data validation

---

## 6. Checkpoint Caching

### 6.1 The `load_or_compute()` Pattern

This is the **core caching mechanism**. Every expensive computation uses it:

```r
result <- load_or_compute(
  checkpoint_file = "1.1_expensive_result.rds",
  description = "Expensive computation",
  force_recompute = FALSE,  # Set TRUE during development
  compute_fn = function() {
    # Code that takes >1 minute
    expensive_function(data)
  }
)
```

### 6.2 What to Cache

| Always Cache | Sometimes Cache | Never Cache |
|--------------|-----------------|-------------|
| DGEList (raw, normalized) | Intermediate DE steps | Plot objects |
| limma fit objects | Gene annotations | Temporary variables |
| DE results (all genes) | Custom gene sets | |
| GSEA results per database | | |
| TF/pathway activities | | |

### 6.3 Checkpoint Naming Convention

```
{phase}.{step}_{description}.rds

Examples:
1.1_dge_raw.rds           # Initial DGEList
1.1_dge_normalized.rds    # After TMM
1.1_fit_object.rds        # limma fit
1.1_de_results.rds        # topTable output
1.1_gsea_H_C2.rds         # GSEA for Hallmark + C2
1.1_gsea_C5.rds           # GSEA for GO terms
1.4_gsea_custom.rds       # Custom gene set GSEA
1.7_tf_activities.rds     # TF analysis results
```

### 6.4 Force Recompute Toggle

During development, set `force_recompute = TRUE` to regenerate:

```r
# Development mode
result <- load_or_compute(..., force_recompute = TRUE)

# Production mode (default)
result <- load_or_compute(..., force_recompute = FALSE)
```

---

## 7. Master Tables Pattern

### 7.1 Purpose

Master tables serve as the **contract between R and Python**:
- R computes complex objects (gseaResult, limma fit)
- Master tables normalize to simple DataFrames
- Python reads CSVs without needing R objects

### 7.2 Schema Standardization

All GSEA results normalized to this schema:

```r
normalize_gsea_results <- function(gsea_obj, database, contrast) {
  tibble(
    pathway_id = gsea_obj@result$ID,
    pathway_name = clean_pathway_name(gsea_obj@result$Description),
    database = database,
    nes = gsea_obj@result$NES,
    pvalue = gsea_obj@result$pvalue,
    padj = gsea_obj@result$p.adjust,
    set_size = gsea_obj@result$setSize,
    core_enrichment = gsea_obj@result$core_enrichment,
    contrast = contrast,
    direction = ifelse(gsea_obj@result$NES > 0, "Up", "Down")
  )
}
```

### 7.3 Master Table Types

| Table | Content | Use |
|-------|---------|-----|
| `master_gsea_table.csv` | All GSEA results | Pathway visualization |
| `master_de_table.csv` | All DE results | Volcano, heatmaps |
| `master_tf_activities.csv` | TF scores | TF analysis plots |
| `master_progeny_activities.csv` | Pathway activities | PROGENy viz |

### 7.4 Validation Pattern

```python
# Python validation
def validate_schema(df, table_name, schemas):
    required = schemas[table_name]['required_columns']
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns in {table_name}: {missing}")
```

---

## 8. Separation of Concerns

### 8.1 Layer Responsibilities

| Layer | Location | Does | Does NOT |
|-------|----------|------|----------|
| **Config** | `config/` | Define parameters | Contain logic |
| **Processing** | `1.x.*.R` | Compute results | Generate plots |
| **Aggregation** | `1.x.*.R` | Create master tables | Visualize data |
| **R Visualization** | `2.x.*.R` | Create plots | Recompute results |
| **Python Viz** | `3.x.*.py` | Interactive dashboards | Read RDS files |

### 8.2 Anti-Patterns to Avoid

```r
# BAD: Mixing computation and visualization
volcano_plot <- function(dge, contrast) {
  # Computing DE in viz function - WRONG
  fit <- lmFit(dge, design)
  fit2 <- eBayes(fit)
  results <- topTable(fit2, coef = contrast)

  # Then plotting
  ggplot(results, aes(logFC, -log10(P.Value))) + ...
}

# GOOD: Visualization only reads pre-computed
create_volcano_plot <- function(de_results, contrast_name) {
  # de_results already computed and loaded
  ggplot(de_results, aes(logFC, -log10(P.Value))) + ...
}
```

### 8.3 Data Flow Between Layers

```
Processing Layer (1.x)
       |
   Checkpoints (RDS)
       |
Aggregation Layer (1.x)
       |
   Master Tables (CSV)
       |
   +-------+
   |       |
R Viz   Python Viz
(2.x)   (3.x)
```

---

## 9. DRY Principles

### 9.1 Parameterize, Don't Duplicate

```r
# BAD: Separate scripts per database
# run_gsea_hallmark.R
# run_gsea_kegg.R
# run_gsea_reactome.R

# GOOD: Single parameterized script
DATABASES <- list(
  H = c("H", ""),
  C2_KEGG = c("C2", "CP:KEGG"),
  C2_REACTOME = c("C2", "CP:REACTOME")
)

for (db_name in names(DATABASES)) {
  gsea_results[[db_name]] <- run_gsea(
    de_results,
    category = DATABASES[[db_name]][1],
    subcategory = DATABASES[[db_name]][2]
  )
}
```

### 9.2 Extract Common Patterns

```r
# BAD: Repeated plotting code
pdf("plot1.pdf"); plot1; dev.off()
pdf("plot2.pdf"); plot2; dev.off()

# GOOD: Helper function
save_plot <- function(plot, filename, width = 10, height = 8) {
  pdf(file.path(DIR_PLOTS, filename), width = width, height = height)
  print(plot)
  dev.off()
  message("Saved: ", filename)
}

save_plot(plot1, "plot1.pdf")
save_plot(plot2, "plot2.pdf")
```

### 9.3 Use Toolkit Functions

```r
# Source shared toolkit
source_toolkit()

# Use standardized functions
volcano <- create_standard_volcano(de_results, ...)
dotplot <- gsea_dotplot(gsea_results, ...)
```

---

## 10. Code Syntax Standards

### 10.1 R Style Guide

```r
# Function documentation (roxygen2)
#' Run GSEA analysis on DE results
#'
#' @param de_results Data frame with logFC, P.Value, etc.
#' @param rank_metric Column to use for ranking (default: "t")
#' @param species Species for msigdbr (default: "Mus musculus")
#' @return gseaResult object
#' @export
run_gsea <- function(de_results, rank_metric = "t",
                     species = "Mus musculus") {
  # Implementation
}

# Variable naming: snake_case
de_results <- topTable(fit, coef = contrast_name)
gsea_hallmark <- run_gsea(de_results, category = "H")

# Progress messages
message("[STEP] Loading data...")
message("[STEP] Running GSEA for ", db_name, "...")
message("[DONE] Saved ", n_files, " plots")

# Tidyverse pipes for data manipulation
results %>%
  filter(padj < 0.05) %>%
  arrange(desc(abs(NES))) %>%
  head(20)
```

### 10.2 Python Style Guide

```python
"""Module docstring explaining purpose."""

from pathlib import Path
from typing import Dict, List, Optional
import pandas as pd

def load_gsea_data(data_dir: Path) -> pd.DataFrame:
    """Load and validate GSEA master table.

    Args:
        data_dir: Path to tables directory

    Returns:
        DataFrame with validated GSEA results
    """
    df = pd.read_csv(data_dir / "master_gsea_table.csv")
    validate_schema(df, 'master_gsea_table')
    return df

# Path handling with pathlib
PROJECT_ROOT = Path(__file__).parent.parent.parent
DATA_DIR = PROJECT_ROOT / "03_results" / "tables"
```

### 10.3 Common Patterns

**Loading checkpoints:**
```r
# Always use load_or_compute for expensive operations
result <- load_or_compute(
  checkpoint_file = "checkpoint.rds",
  description = "Description",
  compute_fn = function() { ... }
)

# For cheap operations, direct load is OK
config <- readRDS("config.rds")
```

**Creating directories:**
```r
# Always ensure directory exists before writing
output_dir <- file.path(DIR_PLOTS, "GSEA", db_name)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}
```

**Error handling:**
```r
# Graceful fallback for missing data
result <- tryCatch(
  run_analysis(data),
  error = function(e) {
    warning("Analysis failed: ", e$message)
    return(NULL)
  }
)
```

---

## 11. Visualization Standards

### 11.1 Colorblind-Safe Palettes

Use **Okabe-Ito palette** for categorical data:

```yaml
# In pipeline.yaml
colors:
  databases:
    Hallmark: "#E69F00"      # Orange
    KEGG: "#56B4E9"          # Sky Blue
    Reactome: "#009E73"      # Bluish Green
    WikiPathways: "#F0E442"  # Yellow
    GO_BP: "#0072B2"         # Blue
    GO_MF: "#D55E00"         # Vermillion
    GO_CC: "#CC79A7"         # Reddish Purple
```

Use **Blue-White-Orange diverging** for NES/logFC:
```yaml
colors:
  diverging:
    negative: "#2166AC"  # Blue (down)
    neutral: "#F7F7F7"   # White
    positive: "#B35806"  # Orange (up)
```

### 11.2 Plot Quality Standards

```r
# Publication-ready settings
PLOT_WIDTH <- 10
PLOT_HEIGHT <- 8
PLOT_DPI <- 300

# Vector format for publication
ggsave("figure.pdf", plot, width = PLOT_WIDTH, height = PLOT_HEIGHT)

# Raster format for preview
ggsave("figure.png", plot, width = PLOT_WIDTH, height = PLOT_HEIGHT, dpi = PLOT_DPI)
```

### 11.3 Consistent Theme

```r
# Use custom minimal theme
custom_minimal_theme_with_grid <- function(base_size = 12) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.ticks = element_line(color = "black")
    )
}
```

---

## 12. Documentation Requirements

### 12.1 Per-Script Documentation

Every script must have a header:
```r
# {script_name} - {Brief description}
# Project: {PROJECT_ID}
# Phase: {1|2|3}
# Description: {What this script does}
# Dependencies: {Required checkpoints}
# Outputs: {Generated files}
# Updated: {YYYY-MM-DD}
```

### 12.2 Per-Folder READMEs

Every output folder should have a README:

```markdown
# {Folder Name}

Generated: {YYYY-MM-DD}

## Contents
- `file1.pdf`: Description
- `file2.csv`: Description

## Interpretation
- What these outputs show
- How to use them

## Generation
Rscript 02_analysis/2.1.visualizations.R
```

### 12.3 Project-Level Documentation

| File | Purpose |
|------|---------|
| `CLAUDE.md` | AI assistant context |
| `plan.md` | Analysis strategy |
| `tasks.md` | Execution tracker |
| `notes.md` | Research findings |

---

## Appendix A: Quick Reference

### A.1 File Locations

| What | Where |
|------|-------|
| Configuration | `02_analysis/config/pipeline.yaml` |
| Checkpoints | `03_results/checkpoints/*.rds` |
| Master tables | `03_results/tables/master_*.csv` |
| Plots | `03_results/plots/{theme}/` |
| Toolkit | `01_scripts/RNAseq-toolkit/` |

### A.2 Common Commands

```bash
# Run core pipeline
Rscript 02_analysis/1.1.core_pipeline.R

# Run all visualizations
for script in 02_analysis/2.*.R; do Rscript "$script"; done

# Generate interactive dashboard
python 02_analysis/3.1.pathway_explorer.py
```

### A.3 Troubleshooting

| Problem | Solution |
|---------|----------|
| GSEA slow | Check if checkpoint exists |
| Missing genes | Verify gene ID type (Symbol vs Entrez) |
| Plot errors | Ensure master tables are current |
| Import errors | Validate CSV schema |

---

**End of Guidelines**
