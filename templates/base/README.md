# Base Project Template

Standard single-cell bioinformatics project structure for RNA-seq and multi-modal analyses.

## Directory Structure

After running `init-project.sh` + `setup-ai.sh`:

```
.
├── 00_data/
│   ├── raw/           # Raw sequencing data (CellRanger outputs, etc.)
│   ├── processed/     # Processed Seurat/AnnData objects
│   └── references/    # Reference files (GTF, genome builds)
├── 01_modules/        # Reusable analysis toolkits (submodules)
│   ├── RNAseq-toolkit/
│   └── SciAgent-toolkit/
├── 02_analysis/
│   ├── config/        # Project configuration files
│   │   ├── config.R
│   │   ├── pipeline.yaml
│   │   ├── color_config.R
│   │   └── analysis_config.yaml  # Created by setup-ai.sh
│   └── helpers/       # Project-specific scripts
├── 03_results/
│   ├── checkpoints/   # Cached intermediate objects
│   ├── plots/         # Generated figures
│   └── tables/        # Output tables
├── logs/              # Analysis logs
├── .claude/           # Claude Code configuration (created by setup-ai.sh)
├── CLAUDE.md          # AI context (created by setup-ai.sh)
├── context.md         # Scientific context (created by setup-ai.sh)
├── tasks.md           # Task tracker
└── notes.md           # Research notes
```

## Getting Started

1. Initialize project: `./init-project.sh ~/projects/my-analysis base`
2. Open in VS Code: `code ~/projects/my-analysis`
3. Reopen in container: Ctrl+Shift+P → "Dev Containers: Reopen in Container"
4. Run AI setup: `./01_modules/SciAgent-toolkit/scripts/setup-ai.sh`
5. Fill in `context.md` with your scientific question
6. Configure `02_analysis/config/pipeline.yaml`

## Recommended Packages

- **R:** Seurat, ggplot2, dplyr, harmony, edgeR, limma
- **Python:** scanpy, anndata, scvi-tools (use `/opt/venvs/base`)

## Workflow

1. Place raw data in `00_data/raw/` (or mount via docker-compose.yml)
2. Configure species/genome in `02_analysis/config/config.R`
3. Use methodology guidelines in `01_modules/SciAgent-toolkit/docs/guidelines/`
4. Cache intermediate results in `03_results/checkpoints/`
5. Generate figures in `03_results/plots/`
