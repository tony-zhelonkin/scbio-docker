# AGENTS.md

Orientation for AI agents working in **scbio-docker**. Kept deliberately lean —
it carries the facts you can't cheaply re-derive each session and points to
`docs/` for depth. Do not duplicate `docs/` here; fix the doc and link it.

## What this repo is

A container substrate: Docker images + VS Code dev-container templates for
single-cell RNA-seq / epigenomics analysis in R and Python. Reproducible,
size-optimized, shareable. It is **containerization-only** — it ships no AI
tooling, only the prerequisites for it (see the boundary below).

Canonical version lives in [`VERSION`](VERSION). Version history lives in
[docs/changelog.md](docs/changelog.md) — the one place versions are written out.
Default branch: `main`.

## Ground truth (current image)

| Fact | Value |
|------|-------|
| Base image | `scdock-r-dev:$(cat VERSION)`, multi-stage from `ubuntu:22.04`, `docker/base/Dockerfile`, true ~20GB |
| R | 4.5.3 from source (Cairo/BLAS/LAPACK/R-shlib) + Bioconductor 3.22 |
| Python | 3.11 (deadsnakes); base venv `/opt/venvs/base` |
| R kernel / Py kernel | Jupyter `ir` and `python311-scagent` registered in-image |
| CLI bio tools | samtools/bcftools/htslib 1.21, bedtools 2.31.1, MACS3, scIBD |
| Baked tooling | tmux 3.7b (source), Quarto, jupyter-scatter, jq, pandoc, ripgrep, fd, TinyTeX |
| AI prereqs only | Node 20, `uv`/`uvx`, Python `toml`, `curl`, `jq` |
| ArchR | **legacy, maintenance dropped** — official `greenleaflab/archr:1.0.3-base-r4.4`, `dev-archr` behind the `archr` compose profile; frozen (no further updates) and on the path to removal |

## Where things live

- `docker/base/Dockerfile` — the build.
- `docker/base/R/{install_core.R,install_httpgd.R,install_renv_project.R}` — R installers (NOT under `.devcontainer/`).
- `docker/requirements/{base,squid,atac,comms}.txt` — Python stacks.
- `scripts/` — `build.sh`, `build-archr-wrapper.sh`, `init-container.sh`, `create_layered_venv.sh`, `poststart_sanity.sh`.
- `init-project.sh` → symlink to `scripts/init-container.sh`.
- `templates/devcontainer/` — the ONLY template tree (devcontainer.json/compose/.env/.vscode + `scripts/setup_ai_env.sh`).
- `.devcontainer/` — a rendered example, not the source of truth.
- `toolkits/SciAgent-toolkit/` — submodule (the AI harness; see boundary).
- `toolkits/refcache/` — submodule (own repo): shared reference-data cache tooling.
  The *mechanism* only; the bytes live on a host path passed via `--refcache`.
- `docs/` — all detailed docs; `docs/README.md` is the map.

Note: `scripts/provc` and `scripts/setup_pi.sh` are **workstation-private** tooling,
not scbio-docker features — don't document or surface them.

## The boundary (read before scaffolding anything)

scbio-docker owns the **container**; SciAgent-toolkit owns the **project + AI
harness**; refcache owns the **reference data**. They compose at seams and stay
independent.

| Repo | Owns |
|------|------|
| **scbio-docker** | Dockerfiles, env specs, build scripts, devcontainer/compose templates, `init-container.sh`, `.env` stub |
| **SciAgent-toolkit** | Project tree, analysis config, docs namespaces, AI harness (agents/skills/commands) |
| **refcache** | Reference-data fetchers, snapshot/verify/prune driver, fetcher image |

Ownership test: *does the content change when the **image** changes → scbio-docker;
when the **project / AI harness** changes → SciAgent-toolkit; when **upstream
reference data** changes → refcache.* SciAgent-toolkit is attached per-project at
`01_modules/SciAgent-toolkit/`; the image stays free of it. refcache is attached
at `toolkits/refcache/` and meets the container at one seam: the `:ro` mount.

Two-step workflow: `./init-project.sh <dir>` (render container) → `sciagent new
project --type analysis <dir>` (scaffold) → open in VS Code → run `setup-ai.sh` once.

## Command cheat-sheet

```bash
scripts/build.sh                 # generic build (devuser:1000, shareable); --personal for your UID
./init-project.sh <dir> [--data-mount k:PATH[:rw]] [--service dev-core|dev-archr] [--gpu]
usepy base|squid|atac|comms|scenic      # switch Python env (layered venvs created on first use)
r-base                           # radian on base R libs
docker run --rm -v "$PWD:/repo" -w /repo scdock-r-dev:$(cat VERSION) bash scripts/poststart_sanity.sh  # sanity script is mounted, not baked into the image
```

## Gotchas

- **Version-agnostic edits.** Don't hardcode `vX.Y.Z` in docs/scripts; read `VERSION`. Only `docs/changelog.md` spells out versions.
- **Two-tier R libs.** System `/usr/local/lib/R/library` is read-only (renv-pinned); runtime installs land in writable `~/R/...` with no sudo. The "installation paths not writeable" notice is expected — use `update = FALSE` to quiet it.
- **A green build does not mean a complete package set.** `safe_install()` tolerates per-package failures by design. It now verifies each package is loadable and writes `/opt/settings/install_failures.csv` — **check that file after every build**, not the exit code. Before this existed, the blind spot hid seven packages (`chromVAR` among them) for two releases.
- **C++ standard floor.** `/usr/local/lib/R/etc/Makevars.site` forces `CXX11STD`/`CXX14STD` to `-std=gnu++17`. Packages declaring `CXX_STD = CXX11` otherwise fail against current RcppArmadillo. Don't drop it.
- **UID remapping.** Generic image is `devuser:1000`; VS Code remaps via `updateRemoteUserUID`, Compose via `LOCAL_UID/LOCAL_GID`, `docker run` via `-u $(id -u):$(id -g)`.
- **`USE_ARCHR` is deprecated** — use the ArchR image/profile, not the toggle.
- Heavy annotation packages are NOT pre-installed, with documented exceptions: `org.Hs.eg.db`, `org.Mm.eg.db`, `EnsDb.Mmusculus.v79`, and **Azimuth's human reference chain** (`BSgenome.Hsapiens.UCSC.hg38`, `EnsDb.Hsapiens.v86`, `JASPAR2020`, ~700MB). The Azimuth chain was previously pulled in invisibly as a transitive dep; it is now declared explicitly in `install_core.R`. Everything else installs at runtime.

## Docs map

Start at [docs/README.md](docs/README.md). Most-used:
[architecture.md](docs/architecture.md) ·
[build.md](docs/build.md) ·
[environments.md](docs/environments.md) ·
[devcontainer.md](docs/devcontainer.md) ·
[ai-integration.md](docs/ai-integration.md) ·
[operations.md](docs/operations.md) ·
[repo-structure.md](docs/repo-structure.md) ·
[toolkits/refcache/README.md](toolkits/refcache/README.md) ·
[changelog.md](docs/changelog.md)
