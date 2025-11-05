# scbio-docker: Development Tasks & Roadmap

## Current Status: v0.5.1 (Released + Enhanced)

**Date:** 2025-11-04
**Branch:** dev (core) / dev-claude-integration (AI workflow extension)
**Tag:** v0.5.1

### Completed in v0.5.1
- [x] Multi-stage build implementation (20GB true size)
- [x] Two-tier R library architecture documentation
- [x] Expected BiocManager warnings documented
- [x] CLAUDE.md created for Claude Code integration
- [x] Generic shareable images with UID remapping
- [x] **Enhanced init-project.sh with two-branch architecture** (2025-11-04)
  - [x] Core functionality on dev branch (templates/docs/, interactive config)
  - [x] Claude integration on dev-claude-integration branch (templates/claude/, agents)
  - [x] Documentation templates (README.md, plan.md, tasks.md)
  - [x] Interactive data mount configuration (--interactive, --data-mount)
  - [x] Command-line options (--git-init, --submodules)
  - [x] WORKFLOW.md documenting "Context > Speed" philosophy
  - [x] Agent stubs for future automation (handoff-writer, stage-reviewer)
  - [x] QUICK_START.md and INIT_PROJECT_ENHANCEMENT.md guides
  - **Time saved:** ~33 hours/year at 3 projects/week

### Known Issues Identified
1. **tidyverse meta-package missing** (components installed, `library(tidyverse)` fails)
2. **BiocManager update warnings confusing** (harmless but noisy)
3. **Python layered venvs not auto-created** (first-use delay)

---

## Phase 1: Robustness & UX Polish (v0.5.2 - v0.6.0)

### v0.5.2: R Package Installation Fixes

**Goal:** Ensure all declared core packages install correctly, eliminate silent failures

**Priority:** HIGH (User-facing installation issues)

**Estimated Effort:** 1-2 days

#### Task 1.1: Fix tidyverse Meta-Package Installation

**Problem Analysis:**
```r
# Current behavior
safe_install <- function(pkgs, installer, ...) {
  for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      installer(pkg, ...)
    }
  }
}

# When tidyverse components (ggplot2, dplyr, etc.) are already installed:
require("tidyverse")  # Returns TRUE (loads ggplot2 silently)
# -> installer never runs
# -> tidyverse meta-package never installs
```

**Root Cause:**
- `require("tidyverse")` succeeds if ANY tidyverse component loads
- Meta-package installation is skipped
- Result: `library(tidyverse)` fails with "package not found"

**Implementation Plan:**

**Option A: Check installed.packages() instead of require()**
```r
safe_install <- function(pkgs, installer, ...) {
  installed <- installed.packages()[, "Package"]
  for (pkg in pkgs) {
    if (!(pkg %in% installed)) {
      message(sprintf("Installing %s ...", pkg))
      tryCatch({
        installer(pkg, ...)
      }, error = function(e) {
        warning(sprintf("Failed to install %s: %s", pkg, e$message))
      })
    } else {
      message(sprintf("%s already installed", pkg))
    }
  }
}
```

**Option B: Force install meta-packages explicitly**
```r
# After component installations
meta_packages <- c("tidyverse", "BiocGenerics")
for (meta in meta_packages) {
  if (!(meta %in% installed.packages()[, "Package"])) {
    install.packages(meta, repos = "https://cloud.r-project.org")
  }
}
```

**Recommended:** Option A (more robust, catches all failures)

**Subtasks:**
- [ ] Modify `.devcontainer/install_R_core.R` to use installed.packages() check
- [ ] Add explicit logging for each package (installed vs skipped vs failed)
- [ ] Create post-install verification script
- [ ] Test build locally
- [ ] Extract updated renv.lock
- [ ] Update R-packages-manifest.csv

**Verification Script:**
```r
# .devcontainer/scripts/verify_core_packages.R
declared_packages <- c("tidyverse", "Seurat", "edgeR", ...)
installed <- installed.packages()[, "Package"]
missing <- setdiff(declared_packages, installed)

if (length(missing) > 0) {
  stop(sprintf("Missing packages: %s", paste(missing, collapse=", ")))
} else {
  message("All core packages verified OK")
}
```

#### Task 1.2: Improve Error Logging

**Goal:** Make build failures visible, not silent

**Subtasks:**
- [ ] Redirect install output to `/opt/settings/install_R_core.log`
- [ ] Log successful installs to `/opt/settings/R_core_success.txt`
- [ ] Log failures to `/opt/settings/R_core_failures.txt`
- [ ] Add summary at end of install_R_core.R:
  ```r
  cat(sprintf("\nCore R package installation complete\n"))
  cat(sprintf("  Successful: %d\n", length(successes)))
  cat(sprintf("  Failed: %d\n", length(failures)))
  if (length(failures) > 0) {
    cat(sprintf("  Failed packages: %s\n", paste(failures, collapse=", ")))
  }
  ```

#### Task 1.3: Add Build-Time Smoke Test

**Goal:** Catch package failures during build, not at runtime

**Subtasks:**
- [ ] Create `.devcontainer/scripts/smoke_test_R.R`:
  ```r
  # Test critical packages
  critical <- c("tidyverse", "Seurat", "Signac", "edgeR", "limma",
                "clusterProfiler", "GSVA", "anndataR")

  for (pkg in critical) {
    tryCatch({
      library(pkg, character.only = TRUE)
      message(sprintf("✓ %s loads OK", pkg))
    }, error = function(e) {
      stop(sprintf("✗ %s FAILED: %s", pkg, e$message))
    })
  }
  ```
- [ ] Add to Dockerfile.optimized before final stage:
  ```dockerfile
  RUN Rscript --vanilla /workspaces/scbio-docker/.devcontainer/scripts/smoke_test_R.R
  ```

**Files Modified:**
- `.devcontainer/install_R_core.R`
- `.devcontainer/scripts/verify_core_packages.R` (new)
- `.devcontainer/scripts/smoke_test_R.R` (new)
- `.devcontainer/Dockerfile.optimized`

**Testing Checklist:**
- [ ] Build completes successfully
- [ ] All 80+ core packages in manifest
- [ ] `library(tidyverse)` works in built image
- [ ] No unexpected failures in log
- [ ] renv.lock matches installed packages

---

### v0.5.3: BiocManager UX Improvements

**Goal:** Reduce confusion from "Installation paths not writeable" warnings

**Priority:** MEDIUM (Documented workaround exists, but UX friction remains)

**Estimated Effort:** 1 day

#### Task 2.1: Add BiocManager Options to .Rprofile

**Implementation:**

```r
# .devcontainer/.Rprofile (after user library setup)

# Reduce BiocManager warning noise
if (interactive()) {
  options(
    # Don't check for repository mismatches (we use CRAN + Bioc intentionally)
    BiocManager.check_repositories = FALSE
  )
}

# Helper message on first library load
if (interactive() && !file.exists("~/.bioc_welcome_shown")) {
  message("\n=== scbio-docker R Environment ===")
  message("Two-tier library design:")
  message("  System: /usr/local/lib/R/library (read-only, core packages)")
  message("  User:   ~/R/... (writable, runtime installs)")
  message("\nTo install packages:")
  message("  BiocManager::install('PACKAGE', update = FALSE)  # Recommended")
  message("  install.packages('PACKAGE')")
  message("\nWarnings about 'paths not writeable' are NORMAL and harmless.")
  message("See README.md for details.\n")

  # Create marker to show message only once per container
  file.create("~/.bioc_welcome_shown")
}
```

**Subtasks:**
- [ ] Update `.devcontainer/.Rprofile` with BiocManager options
- [ ] Test welcome message appears on first radian launch
- [ ] Verify message doesn't repeat on subsequent launches
- [ ] Update documentation to reference this feature

#### Task 2.2: Create Wrapper Functions (Optional)

**Goal:** Provide user-friendly alternatives to BiocManager::install()

**Implementation:**

```r
# .devcontainer/.Rprofile

# Wrapper for clean BiocManager installs
bioc_install <- function(pkgs, ..., update = FALSE) {
  BiocManager::install(pkgs, ..., update = update, ask = FALSE)
}

# Wrapper for GitHub installs
gh_install <- function(repo, ...) {
  remotes::install_github(repo, ..., upgrade = "never", quiet = FALSE)
}
```

**Decision Point:** This adds convenience but also non-standard functions.

**Recommendation:** Add to CLAUDE.md as "helpful patterns" but don't enforce in .Rprofile

**Files Modified:**
- `.devcontainer/.Rprofile`
- `CLAUDE.md` (document wrapper pattern)
- `README.md` (optional: add "Helper Functions" section)

**Testing Checklist:**
- [ ] Welcome message appears on first radian launch
- [ ] BiocManager.check_repositories = FALSE suppresses repo warnings
- [ ] Manual test: `BiocManager::install("AnnotationHub")` has minimal noise
- [ ] Wrapper functions work if implemented

---

### v0.6.0: Python Venv UX & CI Smoke Tests

**Goal:** Improve Python environment first-run experience, add automated testing

**Priority:** MEDIUM-HIGH

**Estimated Effort:** 2-3 days

#### Task 3.1: Pre-Build squidpy Venv (Most Common Use Case)

**Rationale:**
- squidpy (spatial transcriptomics) is commonly used
- Creating at build time adds ~3-5GB but saves 5-minute first-run delay

**Implementation:**

```dockerfile
# .devcontainer/Dockerfile.optimized (builder stage)

# Pre-create squidpy venv (most commonly used)
RUN /opt/scripts/create_layered_venv.sh squid /opt/environments/squid_requirements.txt && \
    /opt/venvs/squid/bin/python -c "import squidpy; print('squidpy OK')"
```

**Trade-off Analysis:**
- **Pros:** Better first-run UX, spatial analysis ready out-of-box
- **Cons:** +3-5GB image size, not everyone uses spatial
- **Decision:** Include in v0.6.0, document how to remove if not needed

**Subtasks:**
- [ ] Modify Dockerfile.optimized to pre-build squid venv
- [ ] Add smoke test for squidpy import
- [ ] Update documentation (squid is pre-built, atac/comms are on-demand)
- [ ] Measure final image size (should be ~23-25GB)

#### Task 3.2: Add Progress Indicators to usepy Switcher

**Goal:** Show user what's happening during venv creation

**Current Behavior:**
```bash
$ usepy atac
# Long pause (2-5 minutes), no feedback
# User wonders if it's frozen
```

**Improved Behavior:**
```bash
$ usepy atac
Creating ATAC venv (first use, ~2-3 minutes)...
  Installing snapatac2...
  Installing episcanpy...
  Installing dependencies...
Done! ATAC venv ready at /opt/venvs/atac
Activated: atac
```

**Implementation:**

```bash
# .devcontainer/scripts/usepy.sh (modified)

if [ ! -d "/opt/venvs/$venv_name" ]; then
  echo "Creating $venv_name venv (first use, ~2-5 minutes)..."
  create_layered_venv.sh "$venv_name" "$requirements_file" 2>&1 | \
    grep -E '(Installing|Collecting|Successfully)' || true
  echo "Done! $venv_name venv ready at /opt/venvs/$venv_name"
fi
```

**Subtasks:**
- [ ] Modify `usepy` to show creation progress
- [ ] Add estimated time to output
- [ ] Test with atac and comms venvs
- [ ] Update CLAUDE.md with expected creation times

#### Task 3.3: Add CI Smoke Tests (GitHub Actions)

**Goal:** Catch regressions before they reach users

**Workflow:** `.github/workflows/test.yml`

```yaml
name: Image Smoke Tests

on:
  push:
    branches: [main, dev]
  pull_request:

jobs:
  build-and-test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3

      - name: Build image
        run: ./build-optimized.sh --tag test:latest

      - name: Test R core packages
        run: |
          docker run --rm test:latest Rscript --vanilla \
            .devcontainer/scripts/smoke_test_R.R

      - name: Test Python base venv
        run: |
          docker run --rm test:latest bash -c \
            "source /opt/venvs/base/bin/activate && \
             python -c 'import scanpy, scvi, muon; print(\"Python OK\")'"

      - name: Test httpgd
        run: |
          docker run --rm test:latest R -e \
            "library(httpgd); print('httpgd OK')"

      - name: Test UID remapping simulation
        run: |
          docker run --rm -u 2000:2000 test:latest id
```

**Subtasks:**
- [ ] Create `.github/workflows/test.yml`
- [ ] Add R smoke test script
- [ ] Add Python smoke test script
- [ ] Configure GitHub Actions secrets (GITHUB_PAT if needed)
- [ ] Test workflow on dev branch
- [ ] Add status badge to README.md

**Files Created:**
- `.github/workflows/test.yml`
- `.devcontainer/scripts/smoke_test_R.R`
- `.devcontainer/scripts/smoke_test_python.sh`

---

## Phase 2: Performance & Scale (v0.6.0 - v0.8.0)

### v0.7.0: HPC Integration & BPCells Validation

**Goal:** Optimize for large datasets, HPC environments

**Priority:** MEDIUM (Power users, production workflows)

**Estimated Effort:** 1 week

#### Task 4.1: BPCells Integration Testing

**Subtasks:**
- [ ] Create example workflow: 1M+ cell dataset with BPCells
- [ ] Document on-disk matrix best practices
- [ ] Benchmark BPCells vs in-memory Seurat (time, RAM, disk I/O)
- [ ] Add to `templates/example-BPCells/` with README

#### Task 4.2: HPC Profile Examples

**Subtasks:**
- [ ] Create Nextflow profile for scbio-docker
- [ ] Slurm job submission template
- [ ] SGE/PBS examples
- [ ] Document UID passthrough for shared filesystems

#### Task 4.3: Singularity Conversion Guide

**Subtasks:**
- [ ] Test `docker2singularity` conversion
- [ ] Document bind mount patterns for HPC
- [ ] Create example Singularity definition file
- [ ] Add to DEVOPS.md

**Files Created:**
- `templates/example-BPCells/README.md`
- `hpc/nextflow.config`
- `hpc/slurm_template.sh`
- `hpc/singularity_conversion.md`

---

### v0.8.0: GPU Support & Advanced Workflows

**Goal:** Enable GPU-accelerated scVI, spatial transcriptomics

**Priority:** LOW-MEDIUM (Specialized use cases)

**Estimated Effort:** 1-2 weeks

#### Task 5.1: GPU Variant Image

**Subtasks:**
- [ ] Create `Dockerfile.gpu` based on nvidia/cuda base
- [ ] Install scvi-tools with GPU support
- [ ] Install TensorFlow GPU for scGLUE
- [ ] Document CUDA version requirements (11.8+)
- [ ] Add GPU smoke tests to CI

#### Task 5.2: Spatial Transcriptomics Workflows

**Subtasks:**
- [ ] Validate squidpy on Visium data
- [ ] Test SpatialData with Xenium/MERFISH
- [ ] Create example notebook: `templates/example-spatial/`
- [ ] Document image processing requirements

---

## Phase 3: Modularity & Community (v0.8.0 - v1.0)

### v0.9.0: R Package Utilities & Enhanced Interop

**Goal:** Package shared code, improve R ↔ Python round-trips

**Priority:** LOW (Nice-to-have, community contribution)

**Estimated Effort:** 2 weeks

#### Task 6.1: Create scbioUtils R Package

**Subtasks:**
- [ ] Extract common plotting functions
- [ ] Add QC metric calculators
- [ ] Package structure with DESCRIPTION, man/
- [ ] Host on GitHub + r-universe
- [ ] Version via semantic versioning

#### Task 6.2: Validate .h5mu Round-Trip Fidelity

**Subtasks:**
- [ ] Test Seurat → MuData → Seurat with embeddings
- [ ] Verify graph preservation (obsp → @graphs)
- [ ] Document lossy conversions (if any)
- [ ] Add to `templates/example-interop/`

---

### v1.0.0: Production Release

**Goal:** Stable, documented, community-ready

**Priority:** MILESTONE

**Estimated Effort:** 1 week (polish + release)

#### Release Checklist

**Documentation:**
- [ ] Complete API stability guarantees
- [ ] Migration guide from v0.5.x → v1.0
- [ ] Video walkthrough or workshop materials
- [ ] FAQ section in README.md

**Infrastructure:**
- [ ] Push to Docker Hub: `tonyzhelonkin/scbio-docker:1.0`
- [ ] Push to GitHub Container Registry: `ghcr.io/tony-zhelonkin/scbio-docker:1.0`
- [ ] Create release notes on GitHub
- [ ] Tag semantic versions (v1.0.0, v1.0.1, etc.)

**Community:**
- [ ] Announce on Bioconductor support forum
- [ ] Tweet/LinkedIn post
- [ ] Submit to awesome-single-cell lists
- [ ] Create contributing guide (CONTRIBUTING.md)

**Testing:**
- [ ] All CI tests passing
- [ ] Manual testing on 3+ HPC systems
- [ ] Validate on Mac/Windows Docker Desktop
- [ ] Community beta testing period (2 weeks)

---

## Backlog (Post-v1.0)

### Potential Future Features

**Image Variants:**
- [ ] Slim variant (no build tools, runtime-only, <15GB)
- [ ] Full variant (includes heavy annotation packages, ~40GB)
- [ ] ARM64 support (Apple Silicon, AWS Graviton)

**Integrations:**
- [ ] JupyterLab extension
- [ ] RStudio Server optional layer
- [ ] Dask for parallel Python workflows

**Advanced Workflows:**
- [ ] CITE-seq template
- [ ] TCR/BCR sequencing tools
- [ ] Long-read scRNA-seq (PacBio, ONT)

**Community Contributions:**
- [ ] Template gallery on GitHub Pages
- [ ] User-submitted workflows
- [ ] Plugin system for custom package sets

---

## Issue Tracking & Project Management

### GitHub Issues Organization

**Labels:**
- `bug` - Something isn't working
- `enhancement` - New feature or request
- `documentation` - Improvements or additions to docs
- `size-optimization` - Image size reduction
- `reproducibility` - renv, version pinning issues
- `hpc` - HPC-specific issues
- `good-first-issue` - Newcomer-friendly

**Milestones:**
- v0.5.2 - R Package Fixes
- v0.5.3 - BiocManager UX
- v0.6.0 - Python Venv & CI
- v0.7.0 - HPC Integration
- v1.0.0 - Production Release

### Development Workflow

**Branch Strategy:**
- `main` - Production releases only (v1.0+)
- `dev` - Active development (current: v0.5.2 work)
- `feature/*` - Feature branches (merge to dev)

**Release Process:**
1. Development on `dev` branch
2. Tag release candidate: `v0.5.2-rc1`
3. Testing period (1 week)
4. Merge to `main` (for v1.0+) or tag on `dev` (for v0.x)
5. Push to registry
6. Create GitHub release with notes

---

## Communication & Documentation

### User-Facing Docs

**README.md** - Quick start, features, troubleshooting
**DEVOPS.md** - Build/run operations, detailed commands
**CLAUDE.md** - AI assistant instructions (Claude Code)
**plan.md** - Project vision (this file)
**tasks.md** - Development roadmap (this file)

### Internal Docs

**CHANGELOG.md** (to create) - Version history
**CONTRIBUTING.md** (to create) - Contribution guidelines
**ARCHITECTURE.md** (to create) - Technical deep dive

---

## Success Metrics & Review Cadence

### Metrics to Track

**Image Quality:**
- Image size (target: <25GB)
- Build time (target: <30 min)
- Package count (R: ~100, Python: ~50)

**User Experience:**
- GitHub stars/forks
- Issue resolution time
- Documentation clarity (measured by repeat questions)

**Reproducibility:**
- renv.lock drift (should be zero between builds)
- CI test pass rate (target: 100%)

### Review Schedule

**Weekly (during active development):**
- Review open issues
- Update tasks.md status
- Merge ready PRs

**Monthly:**
- Security updates (base image, R/Python versions)
- Dependency updates (Bioconductor releases)
- User feedback review

**Quarterly:**
- Major version planning (v0.6, v0.7, etc.)
- Community contributions review
- Roadmap adjustments

---

## Notes & Context from Development History

### Issue: tidyverse Meta-Package Missing (2025-01-24)

**Discovered By:** User testing `library(tidyverse)` in radian

**Investigation:**
```r
# Components installed:
installed.packages()[, 'Package'] |> grep('ggplot2|dplyr|tidyr', value=T)
# [1] "ggplot2" "dplyr" "tidyr" "readr" "purrr" "tibble" "stringr" "forcats"

# Meta-package missing:
'tidyverse' %in% installed.packages()[, 'Package']
# [1] FALSE
```

**Root Cause:** `safe_install()` uses `require()` which succeeds if components present

**Fix Plan:** Check `installed.packages()` instead, force meta-package installs

**Status:** Planned for v0.5.2

### Issue: BiocManager Update Warnings (2025-01-24)

**Discovered By:** User confused by "Installation paths not writeable" message

**Investigation:** Warning is NORMAL - BiocManager checks system library for updates (cannot update because read-only by design)

**Mitigation:** Documented in README.md, CLAUDE.md, DEVOPS.md as EXPECTED behavior

**Future Enhancement:** Add `.Rprofile` options to suppress in v0.5.3

**Status:** Documented (v0.5.1), optional UX fix in v0.5.3

---

## Decision Log

### Decision: Two-Tier R Library Architecture (v0.4-v0.5)

**Date:** 2024-2025 development period

**Rationale:**
- Reproducibility: System packages pinned, user packages experimental
- Shareability: Generic images work for all users
- Disk efficiency: Shared core packages across containers

**Trade-offs:**
- Users see "cannot update" warnings (documented as expected)
- Slight complexity in explaining .libPaths()

**Outcome:** Kept in v0.5.1, extensively documented

### Decision: Official ArchR Image Over Custom Build (v0.5.0)

**Date:** v0.5.0 release

**Rationale:**
- Upstream maintenance (fewer version conflicts)
- Smaller size (15GB vs 30GB custom)
- R 4.4 stability (vs dev-core R 4.5)

**Trade-offs:**
- Need wrapper image for UID consistency
- Two-image strategy (dev-core + ArchR)

**Outcome:** Successful, documented in CLAUDE.md

### Decision: Layered Python Venvs Over Full Duplication (v0.5.0)

**Date:** v0.5.0 release

**Rationale:**
- 60GB+ savings (4 full venvs would be ~100GB)
- Most users only need 1-2 specialized venvs

**Trade-offs:**
- First-use creation delay (2-5 minutes)
- Slight complexity in explaining --system-site-packages

**Outcome:** Successful, considering pre-building squid in v0.6.0

---

## Appendix: Command Reference

### Quick Task Commands

**Start v0.5.2 development:**
```bash
git checkout dev
git pull origin dev
git checkout -b feature/fix-tidyverse-meta-package
```

**Build and test locally:**
```bash
./build-optimized.sh --tag test:v0.5.2
docker run --rm test:v0.5.2 R -e "library(tidyverse)"
```

**Run smoke tests:**
```bash
docker run --rm test:v0.5.2 Rscript --vanilla \
  .devcontainer/scripts/smoke_test_R.R
```

**Extract renv.lock after build:**
```bash
CID=$(docker create test:v0.5.2)
docker cp $CID:/opt/settings/renv.lock ./renv.lock
docker rm $CID
git add renv.lock
git commit -m "Update renv.lock (v0.5.2)"
```

**Tag and release:**
```bash
git tag -a v0.5.2 -m "Fix tidyverse meta-package installation"
git push origin dev
git push origin v0.5.2
```
