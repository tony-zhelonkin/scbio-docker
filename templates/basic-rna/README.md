# Basic RNA-seq Template

Standard single-cell RNA-seq project structure.

## Directory Structure

```
.
├── data/
│   ├── raw/           # Raw sequencing data (CellRanger outputs, etc.)
│   └── processed/     # Processed Seurat/AnnData objects
├── scripts/           # Analysis scripts
├── notebooks/         # Jupyter/Quarto notebooks
├── results/           # Figures, tables, reports
└── renv/              # R package environment (auto-generated)
```

## Getting Started

1. Copy this template to your project directory
2. Place raw data in `data/raw/`
3. Initialize R environment: `renv::init()`
4. Start analysis in `scripts/` or `notebooks/`

## Recommended Packages

- **R:** Seurat, ggplot2, dplyr, harmony
- **Python:** scanpy, anndata, scvi-tools (use `/opt/venvs/base`)
