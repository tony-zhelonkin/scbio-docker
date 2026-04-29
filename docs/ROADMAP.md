# scbio-docker — Roadmap

Consolidated from `plan.md`, `direction.md`, and `tasks.md`. Originals preserved under `docs/archive/` for one cycle.

---

## Vision

A production-ready, reproducible Docker-based development environment for single-cell bioinformatics that balances:

- **Size efficiency** (~20 GB vs typical 100 GB+ bioinformatics containers)
- **Reproducibility** (pinned packages via `renv.lock`, CRAN snapshots, pinned Bioconductor)
- **Flexibility** (runtime package installation, layered Python venvs)
- **Shareability** (generic images, team-friendly UID handling)
- **Developer experience** (VS Code integration, tmux workflows, clear docs)

### Target users

1. Bioinformaticians analyzing scRNA-seq, scATAC-seq, and multimodal data
2. Computational biologists requiring reproducible R + Python workflows
3. Research teams sharing containerized environments across HPC and local machines
4. Students/trainees learning single-cell analysis with consistent tooling

### Core design principles

1. **Reproducibility over convenience.** Read-only system R library (~80 core packages, pinned via `renv.lock` + CRAN snapshots). User library (`~/R/...`) for runtime installs.
2. **Size efficiency through strategic omissions.** Pre-install the 95% use case; defer heavy/specialized tools (BSgenome.*, EnsDb.*, bulk aligners) to runtime.
3. **Layered Python environments over full duplication.** One `/opt/venvs/base` (~25 GB) shared across `squid`, `atac`, `comms` layered venvs (`--system-site-packages`).
4. **Generic images with runtime UID remapping.** Build once as `devuser:1000`; VS Code remaps to actual UID via `updateRemoteUserUID: true`. Personal builds remain available as opt-in.
5. **Official tools over custom builds.** Use `greenleaflab/archr:1.0.3-base-r4.4` directly (wrapper only for UID consistency); TinyTeX over full TeX.

### Strict separation of concerns

| Repository | Responsibility |
|------------|----------------|
| **scbio-docker** | Docker images, container spin-up, project directory structure. Image carries only AI **prerequisites** (`node`, `npm`, `uv`, Python `toml`). |
| **SciAgent-toolkit** (submodule) | All AI tooling: MCP servers, agents, skills, methodology guidelines, role activation. Installed at runtime via `setup-ai.sh`. |

Project flow:

1. `scbio-docker/scripts/init-project.sh` — creates project structure (no AI files).
2. User opens project in container (VS Code Dev Container).
3. `SciAgent-toolkit/scripts/setup-ai.sh` — installs AI tooling, creates `CLAUDE.md`, `.claude/`, `.mcp.json`, etc.

---

## Current Status

**In progress:** v0.5.4 — repo cleanup, AI strip from image, template consolidation, isolation hardening (see Roadmap below).

### v0.5.4 changes (Phase 4 landed)

- Stripped AI tooling from the image — image is now containerization-only. No Claude/Gemini CLI, no MCP servers, no ToolUniverse/Serena/PAL baked in. The `COPY toolkits/SciAgent-toolkit/{scripts,agents}` step is removed; the toolkit is attached per-project as a submodule and `setup-ai.sh` runs at runtime.
- Kept AI prerequisites (Node 20 LTS, `uv`/`uvx`, Python `toml`) so `SciAgent-toolkit/scripts/setup-ai.sh` can run cleanly inside the container.
- `docker/requirements/*.txt` audited and confirmed AI-clean (no `anthropic`, `openai`, `mcp-*`, `claude-*`, etc.). Only `toml` remains in `base.txt` as a deliberate AI prerequisite — it is small and broadly useful.
- `LABEL version="v0.5.4"` set on the runtime stage. `CLAUDE.md` and `README.md` updated to describe the runtime-install model.

**Previously released:** v0.5.3 (AI-prereqs added: Node 20, npm/npx, uv/uvx, Python `toml`).

### Achievements through v0.5.3

- Multi-stage build, true ~20 GB final image (no layer bloat)
- Two-tier R library architecture (system read-only, user writable)
- ~80 core R packages pre-installed (R 4.5 + Bioc 3.21)
- Python base venv + on-demand layered venvs (`squid`, `atac`, `comms`)
- Generic shareable images with UID remapping
- Project scaffolding via `init-project.sh`
- VS Code integration with httpgd graphics
- Official ArchR image integration via wrapper
- AI prerequisites pre-installed; AI tooling deferred to SciAgent-toolkit at runtime

### Known issues

1. **`tidyverse` meta-package missing.** Components install but `library(tidyverse)` fails because `safe_install()` checks via `require()` which returns TRUE on any component. Fix: check `installed.packages()` instead. Targeted in v0.5.2 work; verify in v0.5.4 build.
2. **BiocManager "paths not writeable" warnings.** Harmless (intentional: system library is read-only) but confusing. Documented; optional `BiocManager.check_repositories = FALSE` in `.Rprofile` planned.
3. **Layered venvs not auto-created.** Created on first `usepy` call (~2–5 min delay). Trade-off accepted; pre-building `squid` is under consideration.

---

## Roadmap

### v0.5.4 — Cleanup, AI strip, isolation (in progress)

Tracked in detail in the active refactor plan; phases:

1. **Root tidy.** Move scripts to `scripts/`, meta-docs to `docs/`, consolidate `plan.md` + `direction.md` + `tasks.md` into this `ROADMAP.md`. Add `VERSION` file.
2. **Template consolidation.** Collapse `templates/{base,config,devcontainer,docs,.vscode}/` into a single flat `templates/base/` matching DC_hum_verse layout. Add `docs/{plan, ai-generated, raw}` scaffolding and `.gitignore` with secrets/data patterns.
3. **`init-project.sh` refactor.** Read version from `VERSION`, drop hardcoded `v0.5.x`, generate from new template tree, drop unused multi-template UX (only `base` exists).
4. **Strip AI tooling from image.** Keep only Node 20, `uv`, Python `toml` as prerequisites. Remove any direct AI-tool installation (Claude/Gemini CLI, MCP server installs, ToolUniverse pre-install, Serena pre-build, PAL pre-install). Document split in `CLAUDE.md` and `README.md`.
5. **Filesystem isolation hardening.** Adopt portable patterns from `pi-coding-agent-container`: `tmpfs /tmp` (`noexec,nodev,nosuid`), `pids_limit: 4096`, commented opt-in `read_only` and `secrets:` blocks. Document in `docs/ISOLATION.md`. Skip LD_PRELOAD firewall and SetUID purge — they break interactive R/Python use.
6. **Deprecation cleanup.** Remove `.devcontainer/Dockerfile.archr` and `.devcontainer/install_R_archr.R`. Refresh `README.md` quick-start and version block.
7. **Submodule sync** (orchestrator-handled) — review `toolkits/SciAgent-toolkit` modifications with user.
8. **Build, smoke test, version bump** (orchestrator-handled). Smoke tests confirm AI strip (`! command -v claude`), Node/uv/toml availability, end-to-end `init-project.sh` produces correct tree.

### v0.6.0 — Robustness & UX polish

- **R install robustness.** Switch `safe_install()` to `installed.packages()` check; explicit per-package logging; build-time verification (`smoke_test_R.R`) covering tidyverse, Seurat, Signac, edgeR, limma, clusterProfiler, GSVA, anndataR.
- **BiocManager UX.** Set `BiocManager.check_repositories = FALSE` in `.Rprofile`; one-time welcome message explaining two-tier library design.
- **Python venv UX.** Pre-build `squid` venv in image (most common spatial use case, +3–5 GB; opt-out documented). Add progress indicators to `usepy` showing "Installing snapatac2…" etc.
- **CI smoke tests.** GitHub Actions workflow building image, validating R core packages, Python base imports, httpgd, UID remapping (`-u 2000:2000`).

### v0.7.0 — Performance & scale

- **BPCells integration.** Example workflow on 1M+ cell dataset; benchmark vs in-memory Seurat (time, RAM, disk I/O); `templates/example-BPCells/`.
- **HPC profiles.** Nextflow profile for scbio-docker; Slurm/SGE/PBS submission templates; UID passthrough on shared filesystems.
- **Singularity conversion guide.** Test `docker2singularity`; document bind mounts; example `.def` file.

### v0.8.0 — GPU & advanced workflows

- **GPU variant image.** `Dockerfile.gpu` on `nvidia/cuda` base; scvi-tools GPU; TensorFlow GPU for scGLUE; CUDA 11.8+; GPU smoke tests.
- **Spatial transcriptomics.** Validate `squidpy` on Visium; SpatialData on Xenium/MERFISH; `templates/example-spatial/`.

### v0.9.0 — Modularity & interop

- **`scbioUtils` R package.** Extract common plotting/QC functions; host on GitHub + r-universe; semantic versioning.
- **`.h5mu` round-trip validation.** Seurat → MuData → Seurat with embeddings; verify graph preservation (`obsp` ↔ `@graphs`); `templates/example-interop/`.

### v1.0.0 — Production release

- API stability guarantees; migration guide v0.5.x → v1.0; FAQ; video walkthrough.
- Push to Docker Hub (`tonyzhelonkin/scbio-docker:1.0`) and GHCR (`ghcr.io/tony-zhelonkin/scbio-docker:1.0`).
- Tagged GitHub release with notes.
- Community beta period (2 weeks); validate on 3+ HPC systems.

### Post-v1.0 backlog

- Slim variant (<15 GB, runtime-only)
- Full variant (~40 GB, includes heavy annotation packages)
- ARM64 support (Apple Silicon, AWS Graviton)
- JupyterLab extension layer
- Optional RStudio Server layer
- Dask integration for parallel Python workflows
- Templates: CITE-seq, TCR/BCR, long-read scRNA-seq (PacBio, ONT)
- Plugin system for custom package sets

---

## Non-goals

1. Bulk RNA-seq aligners (STAR, BWA, Salmon, kallisto) — defer to nf-core/rnaseq.
2. "Kitchen sink" coverage of every single-cell tool — 80/20 rule.
3. Web interfaces (RStudio Server, JupyterHub) baked in by default — VS Code Dev Containers is the primary UX.
4. Heavy Mac/Windows native Docker Desktop support — works via Remote-Containers, not a primary target.
5. Auto-updating system packages at runtime — breaks reproducibility.
6. Aggressive isolation (LD_PRELOAD firewall, SetUID purge, command wrappers) — incompatible with interactive R/Python.
7. Migrating existing analysis projects to the new template — they pin specific image versions and stay where they are.

---

## Open questions

- **Pre-bake `squid` venv?** Adds 3–5 GB to image but eliminates 5-min first-run delay. Lean toward yes in v0.6.0; opt-out documented.
- **Wrapper functions (`bioc_install()`, `gh_install()`)?** Add convenience but introduce non-standard API. Lean toward documenting patterns in `CLAUDE.md` rather than baking into `.Rprofile`.
- **Where does `scenicplus` live?** Currently a runtime install into `comms` venv. Move into SciAgent-toolkit's runtime requirements, or keep in scbio-docker `comms.txt`?
- **Registry hosting for v1.0** — Docker Hub vs GHCR vs both? Affects image-pull defaults in templates.
- **AI-only Python packages** — verify `docker/requirements/base.txt` is AI-clean as part of v0.5.4 Phase 4; if anything leaks, decide whether to defer to SciAgent-toolkit runtime or keep.

---

## Success metrics

- Image size: maintain <25 GB for base.
- Build time: <30 min on GitHub Actions runners.
- Reproducibility: identical `renv.lock` → identical package versions across rebuilds.
- Documentation: no recurring user questions about "known issues."
- Adoption: GitHub stars, forks, issue engagement.

## Review cadence

- **Weekly (active development):** review open issues, update roadmap status, merge ready PRs.
- **Monthly:** security updates (base image, R/Python patch versions), dependency updates (Bioconductor releases), user feedback.
- **Quarterly:** major version planning, community contributions, roadmap adjustments tied to Bioconductor cycle.
