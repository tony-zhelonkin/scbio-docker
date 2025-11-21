# Quick Start: init-scproject

Personal quick start. Use what’s useful.

---

## TL;DR

The `init-project.sh` script is available globally as `init-scproject` (symlinked to `~/.local/bin/`).

### Basic usage
```bash
# From anywhere (init-scproject is in PATH)
init-scproject ~/projects/my-analysis basic-rna
```

### Recommended usage (interactive)
```bash
init-scproject ~/projects/my-analysis basic-rna --interactive
# Follow prompts for data mounts, git init, submodules
```

### Power user (all features)
```bash
init-scproject ~/projects/atac-study archr-focused \
  --data-mount atac:/scratch/data/DT-1234 \
  --data-mount rna:/scratch/data/DT-5678:ro \
  --git-init \
  --submodules R_GSEA_visualisations
```

---

## Available Templates

| Template | Description | Service |
|----------|-------------|---------|
| `basic-rna` | Standard RNA-seq analysis (Seurat/scanpy) | dev-core |
| `multimodal` | RNA + ATAC or CITE-seq (integration) | dev-core |
| `archr-focused` | ArchR scATAC-seq analysis | dev-archr |
| `example-DMATAC` | Differential chromatin accessibility | dev-archr |

---

## Command-Line Options

| Option | Description | Example |
|--------|-------------|---------|
| `--interactive` | Prompt for all options (recommended) | `--interactive` |
| `--data-mount KEY:PATH[:ro]` | Add data mount (repeatable) | `--data-mount raw:/scratch/data` |
| `--git-init` | Initialize git repository | `--git-init` |
| `--submodules LIST` | Add submodules (comma-separated) | `--submodules scIBD,R_GSEA` |

**Data mount format:**
- `KEY` = Label for the mount (e.g., `atac`, `rna`, `refs`)
- `PATH` = Host path to data (e.g., `/scratch/data/DT-1234`)
- `:ro` = Optional read-only flag

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
├── .devcontainer/        # VS Code container config
│   ├── devcontainer.json
│   ├── docker-compose.yml
│   └── scripts/poststart_sanity.sh
├── .vscode/              # VS Code settings (Python + R)
│   └── settings.json
├── README.md             # Project guide + workflow philosophy
├── plan.md               # Scientific strategy (~3000 words)
├── tasks.md              # Execution tracker (~2500 words)
├── .env                  # Environment variables (gitignored)
├── .env.example          # Template for team sharing
└── .gitignore            # Excludes data, checkpoints, logs
```

### Additional Files (dev-claude-integration branch only)
```

---

## Workflow After Initialization

### 1. Initialize Project
```bash
# From anywhere (init-scproject is in PATH)
init-scproject ~/projects/my-analysis basic-rna --interactive
```

### 2. Open in VS Code
```bash
cd ~/projects/my-analysis
code .
```

### 3. Reopen in Container
- VS Code will prompt: **"Reopen in Container"**
- Or: `Ctrl+Shift+P` → "Dev Containers: Reopen in Container"

### 4. Verify Environment
The post-start script runs automatically and checks:
- R/Python availability
- Key packages (Seurat, scanpy)
- Library write permissions
- Virtual environment setup

### 5. Fill in Documentation

**plan.md:**
- Scientific question
- Data inventory
- Analysis strategy (Levels 1-4)
- Known issues

**tasks.md:**
- Break plan into stages (S1, S2, S3...)
- Define substeps for each stage
- Add testing criteria
- Track status: ⏸️ (pending), 🔄 (in progress), ✅ (complete), ⛔ (blocked)

### 6. Start Analyzing!

**In R:**
```r
radian  # Launch R terminal
source("02_analysis/config/checkpoint_helpers.R")
# Load data, analyze, save checkpoints
checkpoint_save(obj, "S1_preprocessed")
```

**In Python:**
```bash
usepy base  # or squid, atac, comms
```
```python
import scanpy as sc
# Analyze, save checkpoints
adata.write_h5ad("03_results/checkpoints/S1_preprocessed.h5ad")
```

---

## Branch Usage

### Using dev branch (core only)
```bash
cd ~/pipeline/scbio-docker
git checkout dev
init-scproject ~/projects/my-project basic-rna --interactive
# Creates project WITHOUT Claude integration files
```

### Using dev-claude-integration branch (AI workflow)
```bash
cd ~/pipeline/scbio-docker
git checkout dev-claude-integration
init-scproject ~/projects/my-project basic-rna --interactive
# Creates project WITH Claude integration files (CLAUDE.md, WORKFLOW.md, .claude/)
```

**The script auto-detects which branch you're on** and includes Claude files only if `templates/claude/` exists.

---

## Claude Code Workflow (dev-claude-integration branch only)

### Recommended Approach

**Phase 1: Planning (in plan mode)**
1. Fill in rough draft of `plan.md`
2. Ask Claude to iterate: *"Review and refine this analysis plan"*
3. Ask Claude to create `tasks.md`: *"Create detailed tasks from this plan"*

**Phase 2: Execution (phase-by-phase)**
1. Clear chat (fresh context)
2. Execute Stage 1 substeps
3. Test after each substep
4. Mark ✅ in `tasks.md`
5. Clear chat, move to Stage 2
6. Repeat

**Key Principle: Context > Speed**
- Spend 50 minutes on context, 10 on coding
- Never dump whole project at once
- Work one stage at a time
- Only ONE substep 🔄 at a time

See `WORKFLOW.md` in generated projects for full details.

---

## Examples

### Example 1: Simple RNA-seq Project
```bash
init-scproject ~/projects/pbmc-analysis basic-rna --interactive
# Prompts:
#   Mount label: rna
#   Host path: /scratch/data/pbmc-10k
#   Read-only: y
#   Git init: y
```

### Example 2: Multimodal Project (RNA + ATAC)
```bash
init-scproject ~/projects/multiome-study multimodal \
  --data-mount rna:/scratch/data/DT-1579 \
  --data-mount atac:/scratch/data/DT-1634 \
  --git-init \
  --submodules R_GSEA_visualisations
```

### Example 3: ArchR-Focused ATAC Project
```bash
init-scproject ~/projects/atac-aging archr-focused \
  --data-mount atac:/scratch/data/aging-atlas \
  --data-mount refs:/data/genomes/mm10:ro \
  --git-init
```

### Example 4: Test Project (Disposable)
```bash
init-scproject /tmp/test-project basic-rna --interactive
# Verify structure, then delete
rm -rf /tmp/test-project
```

---

## Troubleshooting

### Command not found: `init-scproject`
```bash
# Reload shell config
source ~/.bashrc

# Or verify PATH
echo $PATH | grep -o "$HOME/.local/bin"

# Check symlink exists
ls -l ~/.local/bin/init-scproject
```

### Template not found
**Solution:** Check template name spelling. Available: `basic-rna`, `multimodal`, `archr-focused`, `example-DMATAC`

### Data mount not working in container
**Solution:**
1. Check host path exists: `ls /scratch/data/...`
2. Verify `docker-compose.yml` has correct mount
3. Rebuild container: `Ctrl+Shift+P` → "Dev Containers: Rebuild Container"

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
# Re-run init-scproject
```

---

## Tips

### Tip 1: Use `--interactive` for your first project
Interactive mode helps you learn what options are available.

### Tip 2: Create `.env.example` for team projects
Share `.env.example` with team, each person creates their own `.env` with actual UIDs.

### Tip 3: Test with `/tmp` first
```bash
init-scproject /tmp/test-project basic-rna --interactive
# Verify structure looks good, then delete
rm -rf /tmp/test-project
```

### Tip 4: Commit early on real projects
```bash
cd ~/projects/my-project
git add .
git commit -m "Initial structure after Stage 1"
```

### Tip 5: Read the generated README.md
Each project gets a README.md with workflow philosophy and best practices specific to your template.

### Tip 6: Use tmux for long-running tasks
```bash
tmux new-session -s my-analysis radian
# Detach: Ctrl+B, then D
# Re-attach: tmux attach -t my-analysis
```

---

## Time Investment vs Savings

| Metric | Value |
|--------|-------|
| **Time to learn** | ~15 minutes (read this guide) |
| **Time per project** | ~2 minutes (vs 15 minutes manual) |
| **Break-even** | After 2 projects |
| **Benefits** | Consistent structure, no forgotten steps, smooth Claude integration |

**At 3 projects per week:**
- Weekly savings: ~39 minutes
- Monthly savings: ~2.6 hours
- Yearly savings: ~33.8 hours

---

## Next Steps

1. ✅ Try creating a test project: `init-scproject /tmp/test basic-rna --interactive`
2. ✅ Review generated files (README.md, plan.md, tasks.md)
3. ✅ Open in VS Code, reopen in container
4. ✅ Verify environment works (sanity checks run automatically)
5. ✅ Create your first real project!

---

## See Also

- **Full documentation:** `CLAUDE.md` — Architecture and technical details
- **DevOps guide:** `DEVOPS.md` — Container operations and workflows
- **Implementation details:** `INIT_PROJECT_ENHANCEMENT.md` — Enhancement summary
- **Template customization:** `templates/` — Modify templates for your needs
- **Branch management:** `BRANCH_MANAGEMENT.md` — How to sync dev branches

---

**Happy analyzing!** 🚀
