# Quick Start: Enhanced Project Initialization

**Version:** 1.0 (2025-11-04)
**Branch:** dev (core) / dev-claude-integration (AI workflow)

---

## TL;DR

### Basic usage (minimal)
```bash
cd ~/pipeline/scbio-docker
./init-project.sh ~/projects/my-analysis basic-rna
```

### Recommended usage (interactive)
```bash
./init-project.sh ~/projects/my-analysis basic-rna --interactive
# Follow prompts for data mounts, git init, submodules
```

### Power user (all features)
```bash
./init-project.sh ~/projects/atac-study archr-focused \
  --data-mount atac:/scratch/data/DT-1234 \
  --data-mount rna:/scratch/data/DT-5678:ro \
  --git-init
```

---

## Available Templates

- **basic-rna** — Standard RNA-seq analysis (Seurat/scanpy workflow)
- **multimodal** — RNA + ATAC or CITE-seq (integration)
- **archr-focused** — ArchR scATAC-seq analysis (uses dev-archr service)
- **example-DMATAC** — Differential chromatin accessibility (from DC_Dictionary)

---

## Command-Line Options

| Option | Description | Example |
|--------|-------------|---------|
| `--interactive` | Prompt for all options | `--interactive` |
| `--data-mount KEY:PATH[:ro]` | Add data mount (repeatable) | `--data-mount atac:/scratch/data:ro` |
| `--git-init` | Initialize git repo | `--git-init` |
| `--submodules LIST` | Add submodules (comma-separated) | `--submodules R_GSEA_visualisations,scIBD` |

---

## What Gets Created

### Directory Structure
```
your-project/
├── 00_data/              # Data storage
│   ├── raw/              # Raw input data
│   ├── processed/        # Intermediate data
│   └── references/       # Genome references
├── 01_scripts/           # External tools, submodules
├── 02_analysis/          # Your analysis scripts
│   └── config/           # Pipeline configuration
├── 03_results/           # All outputs
│   ├── checkpoints/      # R/Python objects
│   ├── plots/            # Figures
│   └── tables/           # Result tables
├── logs/                 # Job logs (tmux outputs)
├── .devcontainer/        # Container config
├── .vscode/              # VS Code settings
├── README.md             # Project guide + workflow philosophy
├── plan.md               # Scientific strategy (~3000 words)
├── tasks.md              # Execution tracker (~2500 words)
├── .env                  # Environment variables (gitignored)
├── .env.example          # Template for team sharing
└── .gitignore            # Excludes data, checkpoints, logs
```

### Additional files (claude-integration branch only)
```
├── CLAUDE.md             # Minimal context file (~650 tokens)
├── WORKFLOW.md           # Claude Code usage philosophy
└── .claude/
    └── agents/
        ├── handoff-writer.md
        └── stage-reviewer.md
```

---

## Workflow After Initialization

### 1. Open in VS Code
```bash
cd ~/projects/your-project
code .
```

### 2. Reopen in Container
- Press `Ctrl+Shift+P` (or `Cmd+Shift+P` on Mac)
- Type: "Dev Containers: Reopen in Container"
- Wait for container to build/start

### 3. Verify Environment
The post-start script runs automatically and checks:
- R/Python availability
- Key packages (Seurat, scanpy)
- Library write permissions

### 4. Fill in Documentation

**plan.md:**
- Scientific question
- Data inventory
- Analysis strategy (Levels 1-4)
- Known issues

**tasks.md:**
- Break plan into stages (S1, S2, S3...)
- Define substeps for each stage
- Add testing criteria

### 5. Start Analyzing!
```r
# In R
source("02_analysis/config/checkpoint_helpers.R")
# Load data, analyze, save checkpoints
checkpoint_save(obj, "S1_preprocessed")
```

```python
# In Python
usepy base  # or squid, atac, comms
import scanpy as sc
# Analyze, save checkpoints
adata.write_h5ad("03_results/checkpoints/S1_preprocessed.h5ad")
```

---

## Switching Branches

### Currently on dev (no Claude integration)

**To get Claude integration files:**
```bash
cd ~/pipeline/scbio-docker
git checkout dev-claude-integration
# Now init-project.sh will include CLAUDE.md, WORKFLOW.md, .claude/agents/
```

### Currently on dev-claude-integration

**To go back to core only:**
```bash
git checkout dev
# Now init-project.sh creates projects without Claude files
```

---

## Claude Code Workflow (claude-integration branch)

### Recommended approach:

**Phase 1: Planning (in plan mode)**
1. Fill in rough draft of plan.md
2. Ask Claude to iterate: "Review and refine this analysis plan"
3. Ask Claude to create tasks.md: "Create detailed tasks from this plan"

**Phase 2: Execution (phase-by-phase)**
1. Clear chat (fresh context)
2. Execute Stage 1 substeps
3. Test after each substep
4. Mark ✅ in tasks.md
5. Clear chat, move to Stage 2
6. Repeat

**Key principle:** Context > Speed
- Spend 50 minutes on context, 10 on coding
- Never dump whole project at once
- Work one stage at a time

See `WORKFLOW.md` in generated projects for full details.

---

## Examples

### Example 1: Simple RNA-seq project
```bash
./init-project.sh ~/projects/pbmc-analysis basic-rna --interactive
# Prompts:
#   Mount label: rna
#   Host path: /scratch/data/pbmc-10k
#   Read-only: y
#   Git init: y
```

### Example 2: Multimodal project (RNA + ATAC)
```bash
./init-project.sh ~/projects/multiome-study multimodal \
  --data-mount rna:/scratch/data/DT-1579 \
  --data-mount atac:/scratch/data/DT-1634 \
  --git-init \
  --submodules R_GSEA_visualisations
```

### Example 3: ArchR-focused ATAC project
```bash
./init-project.sh ~/projects/atac-aging archr-focused \
  --data-mount atac:/scratch/data/aging-atlas \
  --data-mount refs:/data/genomes/mm10:ro \
  --git-init
```

---

## Troubleshooting

### "Template not found"
**Solution:** Check template name spelling. Available: `basic-rna`, `multimodal`, `archr-focused`, `example-DMATAC`

### Data mount not working in container
**Solution:**
1. Check host path exists: `ls /scratch/data/...`
2. Verify docker-compose.yml has correct mount
3. Rebuild container: Ctrl+Shift+P → "Dev Containers: Rebuild Container"

### Permission errors in container
**Solution:**
1. Check `.env` has correct UID/GID: `cat .env`
2. Should match host: `id -u` and `id -g`
3. Fix if needed, rebuild container

### Claude integration files missing
**Solution:**
```bash
cd ~/pipeline/scbio-docker
git branch  # Check you're on dev-claude-integration
git checkout dev-claude-integration
# Re-run init-project.sh
```

---

## Tips

### Tip 1: Use --interactive for first project
The interactive mode helps you learn what options are available.

### Tip 2: Create .env.example for team projects
Share `.env.example` with team, each person creates their own `.env` with actual UIDs.

### Tip 3: Test with /tmp first
```bash
./init-project.sh /tmp/test-project basic-rna --interactive
# Verify structure, then delete
rm -rf /tmp/test-project
```

### Tip 4: Commit early on main projects
```bash
cd ~/projects/my-project
git add .
git commit -m "Initial structure after Stage 1"
```

### Tip 5: Read the generated README.md
Each project gets a README.md with workflow philosophy and best practices.

---

## Time Investment vs Savings

**Time to learn:** ~15 minutes (read this guide)
**Time per project:** ~2 minutes (vs 15 minutes manual)
**Break-even:** After 2 projects
**Benefit:** Consistent structure, no forgotten steps, smooth Claude integration

---

## Next Steps

1. ✅ Try creating a test project: `./init-project.sh /tmp/test basic-rna --interactive`
2. ✅ Review generated files (README.md, plan.md, tasks.md)
3. ✅ Open in VS Code, reopen in container
4. ✅ Verify environment works
5. ✅ Create your first real project!

---

**Questions?** See:
- `INIT_PROJECT_ENHANCEMENT.md` — Full implementation details
- `templates/docs/README.md.template` — Project workflow guide
- `templates/claude/WORKFLOW.md` — Claude Code usage philosophy

---

**Happy analyzing!** 🚀
