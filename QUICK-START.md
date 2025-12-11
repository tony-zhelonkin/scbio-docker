# Quick Start

Personal quick start.

---

## TL;DR

Use `./init-project.sh` to scaffold a new analysis project.
If you have the global wrapper `init-scproject` in PATH, you can use that instead.

### Basic usage
```bash
# From the repo root
./init-project.sh ~/projects/my-analysis basic-rna
```

Or, if `init-scproject` is available globally:
```bash
init-scproject ~/projects/my-analysis basic-rna
```

### Recommended usage (interactive)
```bash
./init-project.sh ~/projects/my-analysis basic-rna --interactive
# Follow prompts for data mounts, git init, resource limits
```

Or:
```bash
init-scproject ~/projects/my-analysis basic-rna --interactive
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

---

## Workflow After Initialization

### 1. Initialize Project
```bash
# From repo root
./init-project.sh ~/projects/my-analysis basic-rna --interactive
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

### Switching Services (core ↔ ArchR)
- Default: dev-core (R 4.5 + Bioc 3.21, Python base)
- ArchR: dev-archr (official ArchR image with R 4.4)
- Change `service` in `.devcontainer/devcontainer.json`, then Reopen in Container

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

---

## Examples

### Example 1: Simple RNA-seq Project
```bash
./init-project.sh ~/projects/pbmc-analysis basic-rna --interactive
# Prompts:
#   Mount label: rna
#   Host path: /scratch/data/pbmc-10k
#   Read-only: y
#   Git init: y
```

### Example 2: Multimodal Project (RNA + ATAC)
```bash
./init-project.sh ~/projects/multiome-study multimodal \
  --data-mount rna:/scratch/data/DT-1579 \
  --data-mount atac:/scratch/data/DT-1634 \
  --git-init
```

### Example 3: ArchR-Focused ATAC Project
```bash
./init-project.sh ~/projects/atac-aging archr-focused \
  --data-mount atac:/scratch/data/aging-atlas \
  --data-mount refs:/data/genomes/mm10:ro \
  --git-init
```

### Example 4: Test Project (Disposable)
```bash
./init-project.sh /tmp/test-project basic-rna --interactive
# Verify structure, then delete
rm -rf /tmp/test-project
```

---

## Troubleshooting

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

### Switch between core and ArchR
- Edit `.devcontainer/devcontainer.json` service
- Reopen in container

---

## Tips

### Tip 1: Use `--interactive` for your first project
Interactive mode helps you learn what options are available.

### Tip 2: Create `.env.example` for team projects
Share `.env.example` with team, each person creates their own `.env` with actual UIDs.

### Tip 3: Test with `/tmp` first
```bash
./init-project.sh /tmp/test-project basic-rna --interactive
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

---

## Next Steps

1. ✅ Create a test project: `./init-project.sh /tmp/test basic-rna --interactive`
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

---

## Adding Submodules (Manual Setup)

Submodules (like RNAseq-toolkit, SciAgent-toolkit) are added **manually** after project creation. This allows per-project branching and customization.

### After Creating Project
```bash
# 1. Initialize project
./init-project.sh ~/projects/my-analysis basic-rna --git-init

# 2. Test in VS Code container (verify it works)

# 3. Exit container, SSH to repo, add submodules
cd ~/projects/my-analysis
git submodule add git@github.com:tony-zhelonkin/RNAseq-toolkit.git 01_scripts/RNAseq-toolkit

# 4. Set per-project branch
git config -f .gitmodules submodule.01_scripts/RNAseq-toolkit.branch dev-myproject
git submodule update --remote

# 5. Commit and reopen container
git add .gitmodules 01_scripts/RNAseq-toolkit
git commit -m "Add RNAseq-toolkit submodule"
```

### For AI Tools (SciAgent-toolkit)
```bash
# Add SciAgent-toolkit submodule
git submodule add git@github.com:tony-zhelonkin/SciAgent-toolkit.git 01_scripts/SciAgent-toolkit

# Inside container, run setup
./01_scripts/SciAgent-toolkit/scripts/setup_mcp_infrastructure.sh
```

This separation keeps `init-project.sh` minimal and stable while allowing flexible per-project customization.
