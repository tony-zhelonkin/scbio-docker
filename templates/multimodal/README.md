# Multimodal Analysis Template

For projects combining RNA + ATAC, CITE-seq, or other multimodal data.

## Directory Structure

```
.
├── data/
│   ├── raw/
│   │   ├── rna/       # RNA-seq data
│   │   ├── atac/      # ATAC-seq data
│   │   └── adt/       # Antibody-derived tags (CITE-seq)
│   └── processed/
├── scripts/
│   ├── rna/           # RNA-specific scripts
│   ├── atac/          # ATAC-specific scripts
│   └── integration/   # Cross-modality integration
├── notebooks/
├── results/
└── renv/
```

## Getting Started

1. Copy this template to your project directory
2. Organize raw data by modality
3. Initialize R environment: `renv::init()`
4. Consider installing MOFA2, Signac, ArchR (if needed)

## Recommended Packages

- **R:** Seurat, Signac, MOFA2, harmony, Azimuth
- **Python:** muon, scvi-tools, scenic+ (use `/opt/venvs/atac` or `/opt/venvs/squid`)
