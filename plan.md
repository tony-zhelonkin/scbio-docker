# scbio-docker: Project Vision and Development Plan

## Project Mission

Create a **production-ready, reproducible Docker-based development environment** for single-cell bioinformatics that balances:
- **Size efficiency** (~20GB vs typical 100GB+ bioinformatics containers)
- **Reproducibility** (pinned packages via renv.lock, CRAN snapshots)
- **Flexibility** (runtime package installation, layered Python venvs)
- **Shareability** (generic images, team-friendly UID handling)
- **Developer experience** (VS Code integration, tmux workflows, clear documentation)

## Target Users

1. **Bioinformaticians** analyzing single-cell RNA-seq, scATAC-seq, and multimodal data
2. **Computational biologists** requiring reproducible R + Python workflows
3. **Research teams** sharing containerized environments across HPC and local machines
4. **Students/trainees** learning single-cell analysis with consistent tooling

## Core Design Principles

### 1. Reproducibility Over Convenience

**Philosophy:** Stable baseline enables scientific reproducibility; experimental packages live in user space.

**Implementation:**
- Read-only system R library (`/usr/local/lib/R/library`) with ~80 core packages
- Packages pinned via `renv.lock` at build time
- CRAN snapshot dates (RSPM) for time-travel reproducibility
- User library (`~/R/...`) for runtime installs and project-specific packages

**Trade-off:** Users see "cannot update system packages" warnings, but gain:
- Identical baseline across all container instances
- No accidental version drift
- Safe testing of newer packages in user library

### 2. Size Efficiency Through Strategic Omissions

**Philosophy:** Pre-install 95% use-case packages; defer heavy/specialized tools to runtime.

**What's included:**
- Core single-cell: Seurat, Signac, scran, scater, edgeR, limma
- Pathway analysis: clusterProfiler, GSVA, fgsea, decoupleR
- Multi-modal: MOFA2, mixOmics, muscat, harmony
- Interop: anndataR, MuDataSeurat, reticulate
- Python base venv: scanpy, scvi-tools, muon, cellrank

**What's deferred to runtime:**
- Heavy annotation packages: BSgenome.*, EnsDb.*, org.*.eg.db (multi-GB each)
- Bulk RNA-seq aligners: STAR, BWA, Salmon, kallisto
- Pre-processing QC: FastQC, Trimmomatic, featureCounts
- Specialized Python tools: scenicplus, scGLUE (optional in layered venvs)

**Result:** 20GB base image vs 100GB+ "kitchen sink" containers

### 3. Layered Python Environments Over Full Duplication

**Philosophy:** Share core packages, specialize where needed.

**Architecture:**
```
/opt/venvs/base (~25GB)
  ├─ Core: scanpy, scvi-tools, muon, cellrank, scvelo, radian
  └─ Foundation: numpy, pandas, scipy, matplotlib, seaborn

/opt/venvs/squid (~3-5GB additional)
  └─ Inherits base via --system-site-packages
  └─ Adds: squidpy, spatialdata (spatial transcriptomics)

/opt/venvs/atac (~2-3GB additional)
  └─ Inherits base
  └─ Adds: snapatac2, episcanpy (scATAC-seq)

/opt/venvs/comms (~3-4GB additional)
  └─ Inherits base
  └─ Adds: liana, cellphonedb, pyscenic (communication & GRN)
```

**Trade-off:** Slight complexity in venv management vs 60GB+ saved (4 full venvs would be ~100GB)

### 4. Generic Images with Runtime UID Remapping

**Philosophy:** Build once, share with team; adapt to user at runtime.

**Generic build (default):**
- Image contains `devuser:1000` (arbitrary, consistent UID)
- Shareable via Docker registry (no baked-in user IDs)
- VS Code remaps 1000 → actual user UID via `updateRemoteUserUID: true`
- Files created in container are owned by actual user on host

**Personal build (optional):**
- Bakes your UID into image (`--build-arg USER_ID=$(id -u)`)
- No remapping needed, but image is user-specific
- Use case: Personal HPC workflows where sharing isn't needed

**Trade-off:** Slight VS Code config complexity vs massive shareability gain

### 5. Official Tools Over Custom Builds

**Philosophy:** Leverage upstream maintenance; wrap for compatibility only.

**Examples:**
- **ArchR:** Use official `greenleaflab/archr:1.0.3-base-r4.4` instead of custom build
  - Maintained by ArchR developers (fewer version conflicts)
  - Wrapper image adds `devuser` for UID consistency
  - Result: 15GB official image vs 30GB custom build

- **TinyTeX:** Minimal LaTeX distribution vs full texlive
  - ~200MB vs 5GB
  - Covers 99% of R Markdown/Quarto use cases

**Trade-off:** Occasionally need to wrap/patch vs maintaining full custom stacks

## Current State (v0.5.1)

### Achievements
- ✅ Multi-stage build: True 20GB final image (no layer bloat)
- ✅ Two-tier R library architecture documented
- ✅ ~80 core R packages pre-installed
- ✅ Python base venv + layered venvs (squid, atac, comms)
- ✅ Generic shareable images with UID remapping
- ✅ Project templates (basic-rna, multimodal, archr-focused)
- ✅ VS Code integration with httpgd graphics
- ✅ Official ArchR image integration

### Known Issues
1. **tidyverse meta-package missing** (components installed, but `library(tidyverse)` fails)
   - Root cause: `safe_install()` checks via `require()`, which succeeds if components present
   - Impact: Users must load ggplot2, dplyr, etc. individually
   - Severity: Low (all components work), but UX friction

2. **BiocManager update warnings confusing** (harmless but noisy)
   - Root cause: `BiocManager::install(update = TRUE)` checks system packages
   - Impact: Users think installation failed (it didn't)
   - Mitigation: Documented as EXPECTED, but still generates confusion

3. **Python layered venvs not auto-created** (created on first `usepy` call)
   - Root cause: On-demand creation to save image size
   - Impact: First use takes 2-5 minutes to install packages
   - Trade-off: 60GB saved vs slight first-use delay

## Strategic Direction (v0.5.x → v1.0)

### Phase 1: Robustness & UX Polish (v0.5.2 - v0.6.0)
**Goal:** Eliminate known issues, improve first-run experience

**Priorities:**
1. Fix R package installation robustness
   - Improve `safe_install()` to detect meta-package failures
   - Ensure tidyverse, BiocGenerics, other meta-packages install correctly
   - Add build-time verification script

2. Reduce BiocManager warning noise
   - Set `options(BiocManager.check_repositories = FALSE)` in `.Rprofile`
   - Or create wrapper functions (`bioc_install()`) with sane defaults
   - Document that `update = FALSE` is best practice for runtime installs

3. Improve Python venv first-run UX
   - Consider pre-building squid venv (most commonly used)
   - Add progress indicators to `usepy` switcher
   - Document expected creation times

4. Add smoke tests to CI
   - Validate core package imports (R and Python)
   - Check httpgd functionality
   - Verify UID remapping in VS Code
   - Test project template scaffolding

### Phase 2: Performance & Scale (v0.6.0 - v0.8.0)
**Goal:** Optimize for large datasets, HPC environments

**Priorities:**
1. BPCells integration testing
   - Validate on-disk matrix workflows
   - Document disk I/O best practices
   - Benchmark vs in-memory Seurat

2. HPC integration patterns
   - Nextflow/nf-core profile examples
   - Slurm/SGE job submission templates
   - Singularity conversion guide

3. GPU support (optional variant)
   - scVI-tools GPU acceleration
   - TensorFlow GPU for scGLUE
   - Document CUDA version requirements

4. Disk I/O optimization
   - Document VS Code file watcher exclusions
   - Recommend tmpfs for intermediate files
   - HDF5 compression strategies

### Phase 3: Advanced Workflows & Modularity (v0.8.0 - v1.0)
**Goal:** Enable complex multi-modal, multi-sample analyses

**Priorities:**
1. Enhanced R ↔ Python interop
   - Validate `.h5mu` round-trip fidelity
   - Document embedding preservation (scVI, PeakVI, SCGLUE)
   - Test graph/neighbor conversions

2. GRN inference workflows
   - SCENIC+ runtime installation guide
   - Memory/time benchmarks
   - Example notebooks (PBMC multiome)

3. Spatial transcriptomics support
   - squidpy venv validation
   - Visium/Xenium workflows
   - SpatialData integration

4. Package modularity
   - Consider splitting core utilities into R package
   - Shared plotting functions, QC metrics
   - Version and distribute via GitHub/r-universe

### Post-v1.0: Maintenance & Community

**Long-term goals:**
1. Quarterly releases tied to Bioconductor versions
2. Community-contributed templates
3. Registry hosting (Docker Hub, GitHub Container Registry)
4. Training materials and workshops
5. Slim/GPU/HPC variants as separate tags

## Non-Goals (Out of Scope)

**What this project will NOT do:**

1. **Include bulk RNA-seq alignment tools** (STAR, BWA, Salmon)
   - Rationale: Use dedicated pipeline containers (nf-core/rnaseq)
   - Defer to specialized tools for upstream processing

2. **Support every single-cell tool** ("kitchen sink" approach)
   - Rationale: Size explosion, dependency conflicts
   - Philosophy: 80/20 rule - cover 80% use cases with 20% packages

3. **Provide web interfaces** (RStudio Server, Jupyter Hub)
   - Rationale: VS Code Dev Containers is the primary UX
   - Alternative: Users can layer RStudio on top if needed

4. **Support Windows/Mac native Docker Desktop extensively**
   - Rationale: Optimized for Linux HPC environments
   - Windows/Mac work via VS Code Remote-Containers but not primary target

5. **Auto-update system packages at runtime**
   - Rationale: Breaks reproducibility
   - Philosophy: Updates go into next image version, not runtime patches

## Success Metrics

**How we measure progress:**

1. **Image size:** Maintain <25GB for base image
2. **Build time:** <30 minutes on GitHub Actions runners
3. **Documentation completeness:** No user questions about "known issues"
4. **Reproducibility:** Same renv.lock → identical package versions
5. **Community adoption:** GitHub stars, forks, issue engagement

## Contributing Philosophy

**What we welcome:**
- Bug reports with reproducible examples
- Documentation improvements
- Template contributions (new project types)
- Performance optimization PRs
- Compatibility testing (different OS, UID setups)

**Review standards:**
- Maintain image size targets
- Preserve reproducibility guarantees
- Add tests for new features
- Update documentation in same PR

## Related Projects & Inspiration

- **vscode-r** (Rami Krispin): VS Code + R integration patterns
- **rocker** (Rocker Project): R-focused Docker images
- **nf-core/modules**: Nextflow best practices
- **greenleaflab/ArchR**: Official ArchR Docker workflows
- **scverse**: Python single-cell ecosystem standards

## Project Initialization System (v0.5.1+)

**Problem:** Creating new analysis projects required manual copying of .devcontainer/, .vscode/, docker-compose.yml, and inconsistent directory structures across projects.

**Solution:** Enhanced `init-project.sh` with two-branch architecture:

### Branch Strategy

**dev branch (core functionality):**
- Enhanced directory structure: `00_data/`, `01_scripts/`, `02_analysis/`, `03_results/`
- Interactive data mount configuration (`--interactive`, `--data-mount`)
- Command-line options (`--git-init`, `--submodules`)
- Documentation templates (README.md, plan.md, tasks.md)
- .env generation and enhanced .gitignore
- **No AI dependencies** - Universal, shareable with any team

**dev-claude-integration branch (AI workflow extension):**
- CLAUDE.md template (minimal context file, ~650 tokens)
- WORKFLOW.md ("Context > Speed" philosophy)
- Agent stubs for future automation (handoff-writer, stage-reviewer)
- **Stacks on top of dev** - Always mergeable, graceful degradation

**Philosophy:** Project initialization is infrastructure, Claude Code integration is optional enhancement.

### Key Features

**Time savings:** ~13 min/project → ~2 min/project (via --interactive mode)
- At 3 projects/week: ~33 hours/year saved
- Consistent structure eliminates forgotten configuration steps

**Documentation templates:**
- README.md: Workflow philosophy ("50 minutes context, 10 minutes coding")
- plan.md: Scientific strategy (~3000 words, non-redundant with tasks.md)
- tasks.md: Execution tracker with status markers (⏸️ ⛔ 🔄 ✅)

**Data mount automation:**
```bash
./init-project.sh ~/projects/my-analysis basic-rna \
  --data-mount atac:/scratch/data/DT-1234 \
  --data-mount rna:/scratch/data/DT-5678:ro \
  --git-init
```

**Token economy (Claude integration):**
- Old CLAUDE.md: ~3000 tokens/session
- New CLAUDE.md: ~650 tokens/session
- **Savings:** 39 full conversation turns over 10 sessions

**Design rationale:**
- Minimal context files (CLAUDE.md) point to detailed docs (plan.md, tasks.md)
- Documentation non-redundancy: plan.md = strategy, tasks.md = execution, notes.md = research
- Agent stubs document future automation (handoff-writer, stage-reviewer)

### Usage

**For team sharing** (no AI):
```bash
git checkout dev
./init-project.sh ~/projects/my-project basic-rna --interactive
```

**For personal workflow** (with Claude Code):
```bash
git checkout dev-claude-integration
./init-project.sh ~/projects/my-project basic-rna --interactive
```

**Result:** Projects have consistent structure, ready-to-use documentation templates, and smooth container integration.

See `QUICK_START.md` and `INIT_PROJECT_ENHANCEMENT.md` for full details.

---

## Version History Highlights

- **v0.1-v0.3**: Initial builds, R 4.4 + Seurat + ArchR custom stacks
- **v0.4**: Two-service compose (dev-core + dev-archr), renv integration
- **v0.5.0**: Core packages approach, layered Python venvs, TinyTeX
- **v0.5.1**: Multi-stage build (20GB true size), documentation overhaul, **enhanced init-project.sh**
- **v0.6.x (planned)**: R package robustness, BiocManager UX fixes
- **v1.0 (goal)**: Production-ready, HPC-tested, community templates
