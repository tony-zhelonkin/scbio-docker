# ArchR-Focused Template

For scATAC-seq analysis with ArchR (R 4.4 environment).

## Important Notes

**This template assumes you're using the `dev-archr` service with the official ArchR image:**
```bash
docker pull greenleaflab/archr:1.0.3-base-r4.4
```

See `DEVOPS.md` for switching to the ArchR container.

## Directory Structure

```
.
├── data/
│   ├── raw/
│   │   └── fragments/  # Fragment files (.tsv.gz)
│   └── processed/
│       └── ArchRProject/  # ArchR output directory
├── scripts/
│   └── archr/          # ArchR scripts
├── notebooks/
├── results/
└── renv/               # May need project-specific packages
```

## Getting Started

1. Switch to ArchR container (see DEVOPS.md)
2. Copy this template to your project directory
3. Place fragment files in `data/raw/fragments/`
4. Initialize ArchR project:
   ```R
   library(ArchR)
   addArchRThreads(threads = 8)
   addArchRGenome("hg38")
   ```

## Recommended Workflow

- Use ArchR for peak calling, motif analysis, trajectory inference
- Export to Seurat/Signac for integration with RNA data (switch back to `dev-core`)
