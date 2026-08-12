# Repository Structure

Where things live in scbio-docker, and what a rendered project looks like after
container init plus `sciagent new project`. scbio-docker is a container
substrate only — image definitions, environment specs, and the devcontainer
templates that wrap a project. Project scaffold and AI harness come from
SciAgent-toolkit (see the table at the bottom).

## scbio-docker layout

```
scbio-docker/
├── VERSION                         # Canonical image tag (source of truth)
├── README.md                       # Repo overview + doc index
├── QUICKSTART.md
├── AGENTS.md, CLAUDE.md            # AI context for this repo
├── renv.lock                       # Pinned R system library
├── R-packages-manifest.csv         # Core R package manifest
├── installed_R_core_packages.csv
├── docker/
│   ├── base/
│   │   ├── Dockerfile              # Multi-stage build (ubuntu:22.04 base)
│   │   └── R/
│   │       ├── install_core.R          # ~80 core R packages
│   │       ├── install_httpgd.R         # CRAN-first + GitHub fallback
│   │       └── install_renv_project.R   # renv restore/snapshot wrapper
│   └── requirements/
│       ├── base.txt               # Core Python stack (base venv)
│       ├── squid.txt              # Spatial transcriptomics (layered)
│       ├── atac.txt               # scATAC-seq tools (layered)
│       └── comms.txt              # Cell communication / GRN (layered)
├── scripts/
│   ├── build.sh                   # Base image build wrapper
│   ├── build-archr-wrapper.sh     # ArchR wrapper image (legacy sidecar)
│   ├── init-container.sh          # Renders devcontainer into a project dir
│   ├── create_layered_venv.sh     # Runtime layered venv helper
│   └── poststart_sanity.sh        # Container startup validation
├── init-project.sh -> scripts/init-container.sh   # Symlink
├── templates/
│   └── devcontainer/              # The ONLY template tree
│       ├── .devcontainer/
│       │   ├── devcontainer.json.template
│       │   ├── docker-compose.yml.template
│       │   ├── .env.example
│       │   └── scripts/
│       │       ├── setup_ai_env.sh
│       │       └── models.agentic.json
│       └── .vscode/
│           └── settings.json
├── .devcontainer/                 # A rendered example (not a template)
│   ├── devcontainer.json
│   ├── docker-compose.yml
│   ├── .Rprofile
│   └── Dockerfile.archr-wrapper
├── toolkits/
│   ├── SciAgent-toolkit/          # Git submodule (attached per-project)
│   └── refcache/                  # Git submodule: reference-data cache tooling
│       ├── refcache.sh            # Snapshot/verify/flip/prune driver
│       ├── sources/               # cistarget.sh, coresh.sh
│       ├── fetcher/Dockerfile     # One image that runs any source
│       └── CHANGELOG.md           # Versioned separately from the image
└── docs/                          # This documentation set
```

Notes:

- **R installers live under `docker/base/R/`**, not `.devcontainer/`. They run
  at image build time.
- **Python requirements live under `docker/requirements/`**. Only `base.txt` is
  baked into the image; `squid`/`atac`/`comms` are installed on demand into
  layered venvs.
- **`toolkits/refcache/` is a separate repo**, on the upstream-data clock: it
  keeps its own `CHANGELOG.md` and moves independently of `VERSION`. It holds
  the fetch *mechanism*; the bytes live on a host path supplied via
  `init-project.sh --refcache`. See
  [../toolkits/refcache/README.md](../toolkits/refcache/README.md).
- **`templates/` contains only `templates/devcontainer/`.** There is no
  `templates/base`, `templates/config`, or `templates/docs` — the project
  scaffold is owned by SciAgent-toolkit.
- **`.devcontainer/` at the repo root is a rendered example**, useful for
  reference. `init-container.sh` produces an equivalent tree in a target
  project directory.
- **`toolkits/SciAgent-toolkit/` is a git submodule**, tracked here but never
  vendored into the image; it is re-attached per project at
  `01_modules/SciAgent-toolkit/`.

## Project layout (after init + `sciagent new project`)

`init-container.sh <dir>` writes the container files; `sciagent new project`
adds the analysis tree and AI harness; `setup-ai.sh` (run once inside the
container) populates AI config.

```
my-project/
├── .devcontainer/                 # Rendered by init-container.sh
│   ├── devcontainer.json
│   ├── docker-compose.yml
│   ├── .env
│   └── scripts/
│       ├── poststart_sanity.sh
│       ├── setup_ai_env.sh
│       └── models.agentic.json
├── .vscode/
│   └── settings.json              # Never clobbers an existing file
├── 00_data/
│   └── <label>/                   # Data mounts land here (read-only default)
├── 01_modules/
│   └── SciAgent-toolkit/          # Git submodule (AI harness, guidelines)
├── 02_analysis/
│   └── config/
│       └── analysis_config.yaml   # Created by setup-ai.sh
├── 03_results/
│   └── checkpoints/
├── .claude/                       # Created by setup-ai.sh
│   ├── agents/
│   └── skills/
├── .mcp.json                      # Created by setup-ai.sh
├── CLAUDE.md, GEMINI.md, AGENTS.md  # Created by setup-ai.sh
├── context.md                     # Created by setup-ai.sh (you fill in)
└── README.md
```

The exact analysis tree is defined by SciAgent-toolkit, not this repo, so treat
the layout above as illustrative.

## Separation of concerns

Decide ownership by asking whether the content changes when the **image**
changes (scbio-docker) or when the **project / AI harness** changes
(SciAgent-toolkit).

| Repository | Owns |
|------------|------|
| **scbio-docker** | Dockerfiles, image definitions, R/Python env specs, build scripts, devcontainer/compose templates, `init-container.sh`, `.env` stub |
| **SciAgent-toolkit** | Project scaffold (directory tree), config templates, docs namespaces, AI harness (roles/skills/agents/commands), methodology guidelines |

## Key paths

| Purpose | Path |
|---------|------|
| Image tag | `cat VERSION` (e.g. `scdock-r-dev:$(cat VERSION)`) |
| Base image build | `scripts/build.sh` or `docker build -f docker/base/Dockerfile` |
| Container init | `./init-project.sh <dir>` (symlink to `scripts/init-container.sh`) |
| Project scaffold | `sciagent new project --type analysis <dir>` |
| AI setup (in container) | `./01_modules/SciAgent-toolkit/scripts/setup-ai.sh` |
| Guidelines | `01_modules/SciAgent-toolkit/docs/guidelines/` |

## Related docs

- [architecture.md](architecture.md) — image layering and design
- [build.md](build.md) — build modes and flags
- [environments.md](environments.md) — R and Python environment specs
- [devcontainer.md](devcontainer.md) — templates and `init-container.sh`
- [ai-integration.md](ai-integration.md) — runtime AI tooling setup
