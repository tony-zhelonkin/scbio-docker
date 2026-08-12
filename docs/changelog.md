# Changelog

All notable changes to the `scdock-r-dev` image and the surrounding
container substrate are documented here. This is the single source of
truth for version history — every other doc stays version-agnostic and
refers to the image as `scdock-r-dev:$(cat VERSION)`.

Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
Newest release first. **[Unreleased]** collects changes queued for the next
image bump — fold them into a version heading when it ships.

## [Unreleased]

## [v0.5.13]

### Added
- **`bulkiRNA` v0.4.0 pinned in the image**, replacing the `RNAseq-toolkit`
  submodule that projects used to `source()` file by file.
- **`gatom` and `mwcsr` for bulkiRNA's GATOM metabolic-network modules.**
  Each falls back to its `ctlab` GitHub source. A failure lands in
  `/opt/settings/install_failures.csv` and the build carries on.
- **An optional-dependency report**, written after the GitHub installs to
  `/opt/settings/bulkirna_optional_deps.csv`: every optional package with its
  presence, version and install command. The build proceeds whatever it finds,
  since these dependencies are optional by design.

  The report keeps its own file, and `install_failures.csv` keeps its meaning:
  requested and failed. `AGENTS.md` sends the reader there after every build, so
  its "no failures" line stays a real signal. A requested package or a
  `bulkiRNA` that fails to load counts as a genuine failure and lands there.

### Fixed
- **`safe_install()` recorded two failures for one package.** A hard install
  error was recorded by the `tryCatch`, then recorded again by the loadability
  check as "install returned without error" — a duplicate, and false. It now
  writes exactly one row per failed package, naming the stage that failed. That
  precision matters, because `install_failures.csv` is the report a reader is
  told to trust.
- **Pinned GitHub slugs now stay pinned through every fallback.**
  `install_gh_pkg()` separates `owner/repo@ref` into repository and ref parts,
  so package-name detection uses the repository while clone and `install_git`
  fallbacks use a valid repository URL and explicitly check out the ref.

## [v0.5.12]

> The reference-data cache tooling lives in its own repo, attached as the
> `toolkits/refcache` submodule, and keeps its own history — see
> [../toolkits/refcache/CHANGELOG.md](../toolkits/refcache/CHANGELOG.md).
> Its releases move independently of `VERSION`.

> This release skips `v0.5.11` in the changelog. A `v0.5.11` bump was prepared
> in-tree and then folded into this one before shipping, so it has no entry of
> its own. Note that a locally-built `scdock-r-dev:v0.5.11` image exists on the
> workstation (built 2026-08-07, from a tree whose `VERSION` still read
> `v0.5.10`); it corresponds to no changelog entry and should not be treated as
> a release.

### Added
- **MALLET 2.0.8 + a headless JRE (`openjdk-17-jre-headless`).** `pycisTopic`'s
  `run_cgs_models_mallet()` shells out to a MALLET install for collapsed-Gibbs
  LDA topic modelling; the pure-Python fallback is impractically slow on real
  cell × region matrices. Installed at `/opt/mallet` with `mallet` symlinked
  onto `PATH` and `MALLET_HOME` / `MALLET_MEMORY=16g` exported (the stock 1g
  JVM heap OOMs). This is the image's only Java consumer — hence the JRE, not
  a JDK.
- **Annotation databases baked in: `org.Hs.eg.db`, `org.Mm.eg.db`.** These
  were the two `org.*.eg.db` packages every project installed at runtime
  anyway (clusterProfiler / GSEA ID mapping), so the "install annotation
  packages at runtime" rule no longer applies to them. Other `org.*`,
  `BSgenome.*` and `EnsDb.*` packages remain runtime installs.
- **Ortholog conversion: `homologene`, `babelgene`, `orthogene`.**
  `babelgene` and `homologene` carry offline human↔mouse mappings (no
  network at call time); `orthogene` (Bioconductor) is the higher-level
  wrapper over both plus gprofiler, and also does ortholog-aware
  conversion of whole expression matrices. `biomaRt` stays the online
  fallback for non-mouse species.
- **TF motif databases for chromVAR: `JASPAR2024` and `chromVARmotifs`.**
  `JASPAR2022` was the only motif source in the image. `JASPAR2024` adds
  the current CORE release (note: it materialises its SQLite DB through
  `BiocFileCache` on first use, so the first call needs network).
  `chromVARmotifs` (`GreenleafLab/chromVARmotifs`, GitHub) provides the
  curated cisBP `human_pwms_v2` / `mouse_pwms_v2` PWM sets that the
  chromVAR and ArchR/Signac deviation workflows assume.
- **`FactoMineR` + `factoextra`** for PCA/MCA/MFA exploratory multivariate
  analysis and its ggplot2 visualisations.
- **`scenic` Python env (`usepy scenic`) — SCENIC+ stack, finally installable.**
  `pycistopic` / `pycistarget` sat commented out in `comms.txt` with a stale
  "not on PyPI for py3.10" note. They are not on PyPI *at all* (both 404 from
  the JSON API), so `docker/requirements/scenic.txt` installs them from git at
  pinned commits. Critically this env is **isolated, not layered**: pycisTopic
  pins `pandas == 1.5` against base's 2.2.3, and a `--system-site-packages`
  layer would shadow base's pandas for scanpy/anndata/muon — precisely the
  breakage `atac.txt` was split off to avoid. Like the other non-base envs it
  is built on first `usepy`, so it adds no image build risk.
  `scripts/create_layered_venv.sh` gained `--isolated` to match.
- **`init-container.sh --refcache PATH`** binds the shared reference-data cache
  `:ro` at `/refcache` and exports `REFCACHE_ROOT`, so analysis code resolves
  reference paths from an env var instead of a host layout or snapshot tag.

### Fixed
- **Seven R packages were silently absent, `chromVAR` and `motifmatchr` among
  them.** `safe_install()` turned install failures into deferred warnings that
  never reached the build log (6 `Installing ...` lines survived for 563
  packages, with zero recorded failures), so a green build hid an incomplete
  package set — for at least two releases. Each had its own cause:
  - `chromVAR` / `motifmatchr`: declare `CXX_STD = CXX11`, compiled at
    gnu++11, which current RcppArmadillo rejects. Fixed by a site-wide
    `Makevars.site` raising `CXX11STD`/`CXX14STD` to `-std=gnu++17`.
  - `WGCNA`: installed with CRAN-only `repos=`, which cannot resolve its
    Bioconductor deps `impute`/`preprocessCore`. Now via `BiocManager`.
  - `mbkmeans`: needs `ClusterR` → `gmp` → `libgmp-dev`, absent from the image.
  - `rliger`: the slug `welch-lab/liger` made the installer look for a package
    named `liger`; it is `rliger`. Its dep `RcppPlanc` needs cmake ≥ 3.24 while
    jammy apt ships 3.22.1, so cmake now comes from pip.
  - `Rfast` / `brms`: install once the above system deps are present.
- **`safe_install()` now reports.** It verifies each package is loadable and
  writes `/opt/settings/install_failures.csv`. **Check that file after a build
  rather than the exit code.** It immediately caught two further problems that
  would otherwise have shipped silently (below).
- **`igraph` failed to load in the runtime image**, taking `SeuratData`,
  `GeneNMF` and `leidenAlg` with it: adding `libglpk-dev` to the builder made
  igraph link GLPK, but the runtime stage had no `libglpk.so.40`. Both
  `libgmp-dev` and `libglpk-dev` are now installed in *both* stages.
- **GitHub installs are retried.** `api.github.com` intermittently drops HTTP/2
  streams from this host; a single attempt lost several packages per build.
  Now 3 attempts with backoff, then a `git`-transport fallback with
  `dependencies = FALSE` that bypasses the API's `Remotes:` resolution.
- **`seurat-disk` now installs before `azimuth`**, which depends on it.
- **The failure report over-reported.** `install_gh_pkg()` recorded a failure at
  the clone+`install_local` stage even when the subsequent `install_git`
  fallback then succeeded. Only the terminal loadability check records now.
- **`build.sh` no longer echoes the first 10 characters of `GITHUB_PAT`** into
  the build log.

### Changed
- **`snapatac2` 2.7.1 → 2.9.0** in `docker/requirements/atac.txt`. 2.9.0 is the
  newest release still publishing a cp311 wheel; 2.10.0 moved to
  `requires-python >=3.12` and cannot be used until the image leaves 3.11.
- **Azimuth's annotation chain is now explicit.** `BSgenome.Hsapiens.UCSC.hg38`,
  `EnsDb.Hsapiens.v86`, `JASPAR2020` and the shiny packages (~700MB) were
  already in the image as invisible transitive deps of Azimuth, contradicting
  the documented "no heavy annotation packages" rule. They are now declared in
  `install_core.R` and the exception is recorded in AGENTS.md and
  environments.md.

### Notes
- Running chromVAR still needs a `BSgenome.*` package for peak sequences;
  those are ~1GB each and stay runtime installs by design.
- The `install_failures.csv` baked into the shipped `v0.5.12` image lists
  `Azimuth`, `SeuratDisk`, `MuDataSeurat` and `PACS` as failed. That is the
  over-reporting bug above: all four are present and loadable. The report is
  correct from the next build onward.

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
