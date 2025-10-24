# Example DMATAC Template

Template for differential chromatin accessibility analysis.

## Directory Structure

```
.
├── data/
│   ├── raw/
│   │   ├── fragments/     # Fragment files per sample
│   │   └── metadata/      # Sample metadata
│   └── processed/
│       ├── peaks/         # Called peaks
│       └── matrices/      # Count matrices
├── scripts/
│   ├── 01_preprocessing/
│   ├── 02_peak_calling/
│   ├── 03_differential/
│   └── 04_motif_analysis/
├── notebooks/
├── results/
│   ├── figures/
│   └── tables/
└── renv/
```

## Getting Started

1. Copy this template to your project directory
2. Update `data/raw/metadata/` with sample information
3. Place fragment files in `data/raw/fragments/`
4. Follow numbered scripts in `scripts/`

## Recommended Tools

- **Peak calling:** MACS2 (Python), ArchR (R)
- **Differential analysis:** DESeq2, edgeR, chromVAR
- **Motif analysis:** JASPAR2024, TFBSTools, motifmatchr
- **Visualization:** ComplexHeatmap, dittoSeq, Gviz

## Notes

- For large datasets, consider using BPCells for on-disk storage
- Use ArchR for trajectory-based differential accessibility
