# Changelog

All notable changes to the `scdock-r-dev` image and the surrounding
container substrate are documented here. This is the single source of
truth for version history — every other doc stays version-agnostic and
refers to the image as `scdock-r-dev:$(cat VERSION)`.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Newest release first. **[Unreleased]** collects changes queued for the next
image bump — fold them into a version heading when it ships.

## [Unreleased]

## [v0.5.10]

### Added
- **`bubblewrap` (`bwrap`) baked into the image** — the sandbox backend
  for the Codex CLI's `workspace-write` mode. Without it, sandboxed exec
  fails outright. `bwrap` still needs unprivileged user namespaces at run
  time; pair with `--security-opt seccomp=unconfined` (a commented
  `security_opt` line is now shipped on the compose dev service). See
  [ai-integration.md](ai-integration.md).
- **`pytest` added to `/opt/venvs/base`** — agents and analysis scripts
  assumed a test runner was importable/on `PATH`; it was missing.
- **`openpyxl==3.1.5` pinned in `docker/requirements/base.txt`** — the
  `.xlsx` reader backend for `pandas.read_excel`. It was present in
  `/opt/venvs/base` but never pinned, so a rebuild silently dropped it;
  the absence broke `score_eTreg` / the mt-hi script. Now reproducible.

### Changed
- **ArchR sidecar maintenance dropped.** ArchR upstream has been
  unmaintained for 2+ years; the `archr` compose profile + wrapper are
  frozen (no further updates) and `scdock-r-archr` is on the path to
  removal once no project depends on it.

### Fixed
- **Image `version` label no longer hardcoded** — it is fed from the
  `VERSION` file via a build-arg (`IMAGE_VERSION`), so it can't drift on
  a bump.
- **Node.js runtime layer** does an explicit `apt-get update` before
  installing `nodejs` rather than relying on the NodeSource script's
  refresh, hardening the layer against reordering.

## [v0.5.9]

### Added
- **tmux 3.7b built from source** (`/usr/local/bin` first on `PATH`) to
  get kitty graphics passthrough for inline image rendering in terminal
  R/Python sessions.
- Devcontainer template carries its own **dev-env layer wiring**: the
  `postStart` sequence optionally runs `/opt/dev-env/rollout.sh --copy
  --no-sync` **only if it is present and executable** (an optional
  personal layer — not a required feature).

### Changed
- Devcontainer template refresh so the dev-env hook is wired after the
  sanity and AI-env setup steps rather than inline.

## [v0.5.8]

### Added
- **Jupyter kernels registered in-image**: the Python kernel
  (`python311-scagent`, from `ipykernel` in the base venv) and the R
  kernel (`ir`, from IRkernel) are installed so notebooks resolve the
  right interpreters out of the box.

### Changed
- **renv refresh** — `renv.lock` / manifest regenerated against the
  current R 4.5.3 + Bioconductor 3.22 core set.
- Image dependency installs consolidated inline in the Dockerfile.

### Fixed
- **VS Code reconnect-freeze** on Remote-SSH: git scanning is capped
  (`git.detectSubmodules=false`, `autoRepositoryDetection=openEditors`,
  `repositoryScanMaxDepth=1`) and heavy I/O paths are excluded from the
  file watcher, preventing reconnect drops. See
  [vscode-remote-stability.md](vscode-remote-stability.md).

## [v0.5.7]

### Added
- Baked **`jupyter-scatter==1.0.1` + `ipykernel`** into the base venv
  for interactive scatter visualisation in notebooks.

## [v0.5.6]

### Added
- **`jq`, `pandoc`, and Quarto CLI (1.9.38)** baked into the image.
- `.devcontainer/scripts/setup_ai_env.sh` bootstraps AI CLIs at
  container start (idempotent, network-tolerant). See
  [ai-integration.md](ai-integration.md).

### Fixed
- **R side-pane terminal**: VS Code R extension settings adjusted
  (`r.alwaysUseActiveTerminal=false`, extension owns its radian
  side-pane) so R code routes to the correct terminal.

## [v0.5.5]

### Added
- **`ripgrep` + `fd`** pre-installed via apt. `fd-find` ships the binary
  as `fdfind`; a `/usr/local/bin/fd` symlink exposes the canonical `fd`
  name. `~/.local/bin` is on `PATH` in interactive shells so
  user-installed tooling resolves.
- **`init-container.sh` renders `.vscode/settings.json`** from the
  template (never clobbering an existing file). Fixes the Shift+Enter
  Python REPL opening on the system interpreter by setting
  `python.defaultInterpreterPath=/opt/venvs/base/bin/python`.

### Changed
- **`fzf` is no longer baked into the image** — it is owned by the
  sibling `dev-env` personal layer, keeping `nvim`/`fzf` out of the
  image to stay lean.

## [v0.5.4]

### Added
- **Expanded R core**: `IRkernel` (Jupyter R kernel), `Rsamtools`,
  `Rfast`, `data.table`, explicit `textshaping`, plus GitHub-installed
  `Zhen-Miao/PICsnATAC` + `Zhen-Miao/PACS`.
- **Filesystem isolation patterns** for the compose template (tmpfs
  `/tmp`, `pids_limit`, optional read-only root, secrets pattern). See
  [isolation.md](isolation.md).

### Changed
- **R 4.5.3 + Bioconductor 3.22** (from 4.5.0 + 3.21). `anndataR` now
  installs cleanly from Bioconductor.
- **Python moved to 3.11** via the deadsnakes PPA (from the Ubuntu
  default 3.10). Base venv and layered venvs
  (`squid`/`atac`/`comms`) build against `python3.11`.
- **Image is now containerization-only**: all AI tooling stripped from
  the build. The image keeps only the **prerequisites** downstream AI
  tooling needs — **Node.js 20 LTS**, **`uv`/`uvx`**, and the Python
  `toml` package — with the actual CLIs installed at runtime by
  `setup_ai_env.sh` / SciAgent-toolkit. See
  [ai-integration.md](ai-integration.md).

### Removed
- Claude/Gemini CLIs, MCP servers, ToolUniverse, Serena, PAL and any
  baked-in agents/skills — no longer part of the image.

## [v0.5.3]

### Added
- **AI prerequisites** kept in the image so downstream `setup-ai.sh`
  runs cleanly: **Node.js 20 LTS** (`node`/`npm`/`npx`), **`uv`/`uvx`**
  (Astral), and the Python `toml` package (Codex CLI config generation).

## [v0.5.2]

### Added
- **Ubuntu libraries**: `libhdf5-dev`, `libgsl-dev` (enables `hdf5r`,
  `DirichletMultinomial` compilation at runtime).
- **R packages**: chromVAR, motifmatchr, TFBSTools, JASPAR2022, SingleR,
  celldex, AnnotationHub, EnsDb.Mmusculus.v79, crescendo (GitHub).

### Fixed
- `safe_install()` now handles meta-packages correctly (`tidyverse`
  installs as expected).

## [v0.5.1]

### Changed
- **Multi-stage build** that discards build artifacts entirely — no
  layer bloat — while preserving build tools in the runtime stage for
  on-demand package compilation.
- **True ~20GB reported image** (previously ~200GB actual / 500GB
  reported).

## [v0.5.0]

### Added
- **Layered Python venvs**: a single base venv (`/opt/venvs/base`) plus
  on-demand `squid`/`atac`/`comms` venvs created with
  `--system-site-packages`, replacing four full venvs (~40GB total vs
  ~100GB). See [environments.md](environments.md).
- **Core + runtime R model**: ~80 core R packages pre-installed via
  `install_core.R`, with a read-only system library plus a writable user
  library for runtime installs.
- Devcontainer / compose templates and `init-container.sh` for quick
  container setup. See [devcontainer.md](devcontainer.md).

### Changed
- Aggressive cache cleanup (renv, pip, build artifacts).
- **TinyTeX** instead of a full TeX distribution.
- **Official ArchR image** (`greenleaflab/archr:1.0.3-base-r4.4`)
  instead of a custom ArchR build; now a deprecated legacy sidecar. See
  [architecture.md](architecture.md).

### Removed
- Bulk aligners (STAR, BWA, Bowtie2, Salmon, kallisto) and
  pre-processing tools (FastQC, Trimmomatic, Trim Galore, featureCounts,
  Picard). Retained: samtools, bcftools, bedtools, scIBD, MACS3.
- Heavy R annotation packages (`BSgenome.*`, `EnsDb.*`, `org.*.eg.db`)
  moved to optional / runtime install.
