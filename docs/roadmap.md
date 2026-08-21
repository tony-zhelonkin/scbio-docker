# scbio-docker — Roadmap

Forward-looking direction for the container substrate. Shipped work and per-version detail live in [changelog.md](changelog.md); this file tracks where the project is headed, not where it has been.

Consolidated from the historical `plan.md`, `direction.md`, and `tasks.md` (preserved under [archive/](archive/)).

---

## Vision

A production-ready, reproducible Docker development environment for single-cell bioinformatics that balances:

- **Size efficiency** — true ~20 GB image vs typical 100 GB+ bioinformatics containers.
- **Reproducibility** — pinned R packages via `renv.lock`, pinned CRAN/RSPM snapshot, pinned Bioconductor.
- **Flexibility** — runtime package installation, layered Python venvs.
- **Shareability** — generic images with team-friendly UID remapping.
- **Developer experience** — VS Code Dev Containers, tmux workflows, clear docs.

### Target users

1. Bioinformaticians analyzing scRNA-seq, scATAC-seq, and multimodal data.
2. Computational biologists needing reproducible R + Python workflows.
3. Research teams sharing containerized environments across HPC and local machines.
4. Students/trainees learning single-cell analysis with consistent tooling.

### Core design principles

1. **Reproducibility over convenience.** Read-only system R library (~80 core packages, pinned via `renv.lock` + CRAN/RSPM snapshot). Writable user library (`~/R/...`) for runtime installs.
2. **Size efficiency through strategic omissions.** Pre-install the 95% use case; defer heavy/specialized tools (`BSgenome.*`, `EnsDb.*`, bulk aligners) to runtime.
3. **Layered Python environments over full duplication.** One `/opt/venvs/base` shared across `squid`, `atac`, `comms` layered venvs (`python3.11 -m venv --system-site-packages`).
4. **Generic images with runtime UID remapping.** Build once as `devuser:1000`; VS Code remaps to the actual UID via `updateRemoteUserUID: true`. Personal builds remain an opt-in.
5. **Official tools over custom builds.** TinyTeX over full TeX; official upstream images preferred over bespoke rebuilds.
6. **Container substrate only.** The image ships AI *prerequisites* (`node`/`npm`, `uv`/`uvx`, `jq`, Python `toml`); all AI tooling is installed at runtime, per project — see the boundary with SciAgent-toolkit below.

### Strict separation of concerns

| Repository | Responsibility |
|------------|----------------|
| **scbio-docker** | Docker images, env specs, build scripts, devcontainer/compose templates, `init-container.sh`. Image carries only AI **prerequisites**. |
| **SciAgent-toolkit** (submodule, attached per-project) | Project tree and catalog, AI harness (agents/skills/commands), methodology guidelines. Never vendored into the image. |

Project flow: (1) `./init-project.sh <dir>` renders the container; (2) with SciAgent-toolkit vendored at `01_modules/SciAgent-toolkit/`, run `./01_modules/SciAgent-toolkit/bin/scio link` and then `./01_modules/SciAgent-toolkit/bin/scio craft` from `<dir>`; (3) open in VS Code and run `setup-ai.sh` once. See [ai-integration.md](ai-integration.md) and [devcontainer.md](devcontainer.md).

---

## Current status

Released versions and their changes are tracked in [changelog.md](changelog.md). The current baseline:

| Area | State |
|------|-------|
| R | 4.5.3 built from source (Cairo/BLAS/LAPACK/R-shlib), Bioconductor 3.22, CRAN/RSPM snapshot default 2026-04-15 |
| Python | 3.11 (deadsnakes), base venv `/opt/venvs/base` + on-demand `squid`/`atac`/`comms` layered venvs |
| R libraries | Two-tier: read-only system (~80 core, `renv`-pinned) + writable user library |
| Image | Multi-stage from `ubuntu:22.04`, true ~20 GB |
| AI model | Containerization-only image (AI prerequisites baked in); AI CLIs bootstrapped at runtime by `setup_ai_env.sh` on every container start |
| ArchR | Deprecated sidecar (`greenleaflab/archr:1.0.3-base-r4.4` + wrapper), on the path to removal |

### Known issues / rough edges

1. **BiocManager "paths not writeable" warnings.** Harmless (system library is intentionally read-only) but confusing to new users. Documented; an `.Rprofile` quieting option is under consideration.
2. **Layered venvs not auto-created.** Built on first `usepy` call (~2–5 min delay). Trade-off accepted; pre-baking `squid` is under consideration (see open questions).

### Deferred — concretely identified (2026-07)

Small, well-scoped items surfaced during the v0.5.10 pass; parked with enough context to pick up cold. Each has a matching `TODO(...)` marker at the code site.

- **Purge tmux build-only deps** (`bison`, `libevent-dev`) from the runtime image. The catch: `libevent-dev` drags in `libevent-2.1-7`, which tmux links against — a naive purge+autoremove breaks tmux. Fix: purge both but reinstall `libevent-2.1-7` in the same layer. Marker: `docker/base/Dockerfile` tmux block.
- **Consolidate the 4–5 runtime `apt-get update && install` layers** into fewer layers (single index refresh). Low risk but touches the whole runtime stage; batch with the next structural change. Marker: `docker/base/Dockerfile` top of runtime stage.
- **Bake `poststart_sanity.sh` into the image?** It is currently mounted at runtime, so bare `docker run <image> scripts/poststart_sanity.sh` fails (docs now show the `-v $PWD:/repo` mount form). Baking it to `/usr/local/bin` would make the cheat-sheet commands work as-is — decide vs. keeping it project-mounted.
- **CI smoke test still unbuilt.** Today's validation was manual (`docker run` core R/Python/kernel checks). The near-term GitHub Actions build+validate workflow (with `smoke_test_R.R`) would automate exactly that.

---

## Roadmap

Forward-looking only. Items already shipped have been removed — see [changelog.md](changelog.md).

### Near term — robustness & UX polish

- **R install robustness.** Explicit per-package logging; a build-time verification pass (`smoke_test_R.R`) covering tidyverse, Seurat, Signac, edgeR, limma, clusterProfiler, GSVA, anndataR.
- **BiocManager UX.** Quiet the read-only-library warnings; a one-time welcome message explaining the two-tier library design.
- **Python venv UX.** Consider pre-baking the `squid` venv (most common spatial use case, +3–5 GB; opt-out documented). Progress indicators on `usepy` showing "Installing snapatac2…".
- **CI smoke tests.** GitHub Actions workflow that builds the image and validates R core packages, Python base imports, httpgd, and UID remapping (`-u 2000:2000`).
- **Retire the ArchR sidecar.** ArchR upstream has been unmaintained for 2+ years and grows more obsolete month by month; the field has moved to per-language stacks (Seurat/Signac in R, scanpy/snapATAC in Python). **Maintenance is being dropped** — the `archr` compose profile and wrapper are frozen (no further updates) and `scdock-r-archr` will be removed once no project depends on it. Fold any remaining guidance into docs.

### Mid term — performance & scale

- **BPCells workflow example.** On-disk analysis of a 1M+ cell dataset; benchmark vs in-memory Seurat (time, RAM, disk I/O).
- **HPC profiles.** Nextflow profile; Slurm/SGE/PBS submission templates; UID passthrough on shared filesystems.
- **Singularity/Apptainer conversion guide.** Bind-mount recipes and an example `.def` file.

### Longer term — GPU & advanced workflows

- **GPU variant image.** `Dockerfile.gpu` on an `nvidia/cuda` base; GPU scvi-tools; GPU smoke tests. (The devcontainer template already exposes a `--gpu` flag for NVIDIA device passthrough on the CPU image.)
- **Spatial transcriptomics validation.** `squidpy` on Visium; SpatialData on Xenium/MERFISH; worked examples.

### Modularity & interop

- **`scbioUtils` R package.** Extract common plotting/QC helpers; host on GitHub + r-universe with semantic versioning.
- **`.h5mu` round-trip validation.** Seurat → MuData → Seurat with embeddings; verify graph preservation (`obsp` ↔ `@graphs`).

### v1.0.0 — production release

- API stability guarantees; a migration guide and FAQ.
- Publish to a public registry (Docker Hub and/or GHCR — see open questions).
- Tagged GitHub release with notes; short community beta validated on 3+ HPC systems.

### Post-v1.0 backlog

- Slim variant (<15 GB, runtime-only) and full variant (~40 GB, heavy annotation packages).
- ARM64 support (Apple Silicon, AWS Graviton).
- Optional JupyterLab / RStudio Server layers.
- Dask integration for parallel Python workflows.
- Templates: CITE-seq, TCR/BCR, long-read scRNA-seq (PacBio, ONT).

---

## Non-goals

1. Bulk RNA-seq aligners (STAR, BWA, Salmon, kallisto) — defer to nf-core/rnaseq.
2. "Kitchen sink" coverage of every single-cell tool — 80/20 rule.
3. Web interfaces (RStudio Server, JupyterHub) baked in by default — VS Code Dev Containers is the primary UX.
4. First-class Mac/Windows native Docker Desktop support — works via Remote-Containers, not a primary target.
5. Auto-updating system packages at runtime — breaks reproducibility.
6. Aggressive isolation (LD_PRELOAD firewall, SetUID purge, command wrappers) — incompatible with interactive R/Python. See [isolation.md](isolation.md).
7. Migrating existing analysis projects to new templates — they pin specific image versions and stay where they are.

---

## Open questions

- **Pre-bake `squid` venv?** Adds 3–5 GB but eliminates the first-run delay. Leaning yes, with a documented opt-out.
- **Where does `scenicplus` live?** Currently a runtime install into the `comms` venv. Keep in scbio-docker `comms.txt`, or move to SciAgent-toolkit runtime requirements?
- **Registry hosting for v1.0.** Docker Hub vs GHCR vs both — affects image-pull defaults in the templates.

---

## Success metrics

- Image size: maintain <25 GB for base.
- Build time: <30 min on CI runners.
- Reproducibility: identical `renv.lock` → identical package versions across rebuilds.
- Documentation: no recurring user questions about "known issues."
- Adoption: GitHub stars, forks, issue engagement.

## Review cadence

- **Active development:** review open issues, update roadmap status, merge ready PRs.
- **Monthly:** security updates (base image, R/Python patch versions), dependency updates.
- **Quarterly:** major-version planning tied to the Bioconductor release cycle.
