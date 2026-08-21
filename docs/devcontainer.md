# Dev Container

`init-container.sh` renders a VS Code dev container (`.devcontainer/` +
`.vscode/settings.json`) into a project directory from the templates in
`templates/devcontainer/`. It is **container substrate only** — it knows nothing
about project structure, analysis trees, config schemas, or AI context. For the
project catalog and AI harness, vendor SciAgent-toolkit and use its CLI (see
[ai-integration.md](ai-integration.md)).

The command is exposed two ways, both identical (one is a symlink):

| Invocation | Path |
|------------|------|
| `scripts/init-container.sh` | canonical script |
| `init-project.sh` | repo-root symlink → `scripts/init-container.sh` |

## Two-step workflow

```bash
# 1. Render the dev container (this repo)
./init-project.sh ~/projects/atac-study \
    --data-mount atac:/scratch/data/DT-1234 \
    --data-mount scratch:/scratch/work/DT-5678:rw

# 2. Bind the toolkit catalog and render CRAFT (SciAgent-toolkit)
cd ~/projects/atac-study
./01_modules/SciAgent-toolkit/bin/scio link
./01_modules/SciAgent-toolkit/bin/scio craft

# 3. Open + reopen in container
code ~/projects/atac-study      # Ctrl+Shift+P -> "Dev Containers: Reopen in Container"
```

## What it renders

Into `<project-dir>/`:

```
.devcontainer/
  devcontainer.json          # from devcontainer.json.template
  docker-compose.yml         # from docker-compose.yml.template
  .env                       # generated (UID/GID, limits, MCP key stubs)
  scripts/
    poststart_sanity.sh      # copied from repo scripts/ (fallback: generated)
    setup_ai_env.sh          # copied from templates/.../scripts/
    models.agentic.json      # copied from templates/.../scripts/
.vscode/
  settings.json              # copied verbatim; NEVER clobbers an existing file
```

`.vscode/settings.json` is only written if the target does not already exist, so
your per-project edits survive re-runs.

## CLI

```
init-container.sh <project-dir> [OPTIONS]
```

| Option | Default | Meaning |
|--------|---------|---------|
| `--data-mount KEY:PATH[:rw]` | — | Add a host data mount (repeatable). **Read-only by default**; append `:rw` to opt out. |
| `--image-version vX.Y.Z` | contents of `VERSION` | Image tag applied to both compose services. |
| `--service dev-core\|dev-archr` | `dev-core` | Which compose service `devcontainer.json` targets. |
| `--max-cpus N` | `50` | CPU limit default written to `.env`. |
| `--max-memory NG` | `450g` | Memory limit default written to `.env`. |
| `--gpu` | off | Add an NVIDIA `devices` block (driver nvidia, count all, `[gpu]`). |

### Data mounts land under `00_data/<label>`

Each `--data-mount KEY:PATH[:rw]` becomes a compose volume mounting the host
`PATH` at `/workspaces/<PROJECT_NAME>/00_data/<KEY>`. Read-only unless `:rw` is
given (input data should stay read-only):

```
--data-mount atac:/scratch/data/DT-1234        -> /workspaces/<proj>/00_data/atac:ro
--data-mount scratch:/scratch/work/DT-5678:rw  -> /workspaces/<proj>/00_data/scratch
```

With no `--data-mount`, the compose file gets a commented placeholder to fill in
later.

### SSH agent auto-detection

`init-container.sh` inspects `SSH_AUTH_SOCK`. If it points at a **live socket**,
the socket is mounted read-only and the compose services set
`SSH_AUTH_SOCK=/ssh-agent` so agent-forwarded keys work inside the container
(e.g. `git push` over SSH). If no live socket is present, the mount token is
left empty and the compose file stays valid — no fallback hacks. See
[ssh-passthrough.md](ssh-passthrough.md).

### The generated `.env`

Written to `.devcontainer/.env` (alongside the compose file, where Docker
Compose auto-loads it). Existing `GEMINI_API_KEY` / `OPENAI_API_KEY` values are
**preserved** across re-runs; everything else is regenerated:

| Key | Value |
|-----|-------|
| `LOCAL_UID` / `LOCAL_GID` | your current `id -u` / `id -g` |
| `WORKSPACE_FOLDER` | `..` (project root relative to `.devcontainer/`) |
| `MAX_CPUS` / `MAX_MEMORY` | from `--max-cpus` / `--max-memory` |
| `OLLAMA_HOST` | `http://172.17.0.1:11434` (Docker bridge → host) |
| `GEMINI_API_KEY` / `OPENAI_API_KEY` | MCP key stubs, consumed later by SciAgent-toolkit |

## Template tree

Everything below lives in `templates/devcontainer/` and is version-agnostic;
`init-container.sh` substitutes the `{{TOKEN}}` placeholders at render time.

```
templates/devcontainer/
  .devcontainer/
    devcontainer.json.template
    docker-compose.yml.template
    .env.example
    scripts/
      setup_ai_env.sh
      models.agentic.json
  .vscode/
    settings.json
```

| Placeholder | Filled with |
|-------------|-------------|
| `{{PROJECT_NAME}}` | `basename` of the project dir |
| `{{SERVICE}}` | `--service` value |
| `{{IMAGE_VERSION}}` | `--image-version` / `VERSION` |
| `{{MAX_CPUS}}` / `{{MAX_MEMORY}}` | resource defaults |
| `{{SSH_AGENT_MOUNT}}` | agent socket volume line (or empty) |
| `{{DATA_MOUNTS}}` | rendered data-mount block |
| `{{GPU_DEVICES}}` | NVIDIA devices block (or empty) |

## `devcontainer.json.template`

Key settings:

- `remoteUser: devuser`, `updateRemoteUserUID: true` — the generic `devuser:1000`
  image is remapped to **your** UID at attach time, so files are owned correctly
  without a personal build. See [architecture.md](architecture.md) for the
  UID-remapping model.
- `dockerComposeFile: docker-compose.yml`, `service: {{SERVICE}}`,
  `workspaceFolder: /workspaces/{{PROJECT_NAME}}`, `shutdownAction: stopCompose`.
- VS Code customizations point the Python interpreter at `/opt/venvs/base/bin/python`
  and the R terminal at `/opt/venvs/base/bin/radian` (the fuller settings live in
  the copied `.vscode/settings.json`).
- **Extensions**: `rdebugger.r-debugger`, `reditorsupport.r`, `quarto.quarto`,
  `purocean.drawio-preview`, `redhat.vscode-yaml`, `yzhang.markdown-all-in-one`,
  `ms-azuretools.vscode-docker`, `ms-vscode-remote.remote-containers`,
  `ms-python.python`, `ms-toolsai.jupyter`.

### `postStartCommand` chain

Runs on every container start:

```
chmod +x .devcontainer/scripts/*.sh
  -> poststart_sanity.sh        # validate environment (see operations.md)
  -> setup_ai_env.sh            # install/refresh AI CLIs (see ai-integration.md)
  -> /opt/dev-env/rollout.sh    # OPTIONAL personal dotfiles layer, skipped if absent
```

The last step is an **optional personal dotfiles layer**: it runs only if
`/opt/dev-env/rollout.sh` is executable, and is silently skipped otherwise. It is
not a required feature of the dev container.

## `docker-compose.yml.template`

Two services share the same mount/user/limit shape:

| Service | Image | Started |
|---------|-------|---------|
| `dev-core` (default) | `scdock-r-dev:{{IMAGE_VERSION}}` | always |
| `dev-archr` | `scdock-r-archr:{{IMAGE_VERSION}}` | only under `profiles: ["archr"]` — `docker compose --profile archr up` |

`dev-archr` is a **deprecated legacy sidecar** (R 4.4 + ArchR 1.0.3); see
[architecture.md](architecture.md) and [changelog.md](changelog.md).

Shared per-service settings:

- `user: "${LOCAL_UID:-1000}:${LOCAL_GID:-1000}"` — runtime UID override from `.env`.
- `env_file: .env`; `SSH_AUTH_SOCK=/ssh-agent`; `OLLAMA_HOST` from `.env`.
- Workspace volume `${WORKSPACE_FOLDER:-.}:/workspaces/{{PROJECT_NAME}}`, plus the
  optional read-only `/opt/dev-env` personal layer, plus rendered SSH-agent and
  data-mount blocks.
- `tmpfs: /tmp` (`noexec,nodev,nosuid,mode=1777`).
- `deploy.resources.limits`: `cpus`/`memory` (from `.env`, defaulting to the
  rendered `{{MAX_CPUS}}`/`{{MAX_MEMORY}}`) and `pids: 4096`;
  `reservations` of 2 CPUs / 8G.

### Commented-out hardening + secrets

The template ships **commented** blocks for a hardened profile (`read_only`,
`cap_drop: [ALL]` + minimal `cap_add`, `security_opt: no-new-privileges`) and a
file-based secrets pattern (`~/.scbio-secrets/<name>` → `/run/secrets/<name>`,
`mode 0000`). Uncomment to enable. Rationale and full recipe live in
[isolation.md](isolation.md).

## `.vscode/settings.json`

Copied verbatim (no tokens), never clobbering an existing file. It encodes
image-path and stability decisions:

- **Python**: `python.defaultInterpreterPath=/opt/venvs/base/bin/python`,
  `python.terminal.activateEnvironment=false`,
  `python.REPL.enableREPLSmartSend=true`, `runInDedicatedTerminal=true`,
  `executeInFileDir=true`. Pointing at the base venv fixes the Shift+Enter REPL
  opening on the system Python (which cannot `import numpy`).
- **Jupyter**: `jupyter.kernels.filter=[]` so the baked `python311-scagent` and
  `ir` kernels stay visible.
- **R**: `r.rterm.linux=/opt/venvs/base/bin/radian`, `r.bracketedPaste=true`,
  `r.sessionWatcher=true`, `r.plot.useHttpgd=true`, `r.lsp.diagnostics=false`.
  `r.alwaysUseActiveTerminal=false` is deliberate — the R extension owns its own
  radian **side pane** so "Run Selection" / `.qmd` chunks never get piped into a
  focused bash terminal (which would throw `bash: syntax error`).
- **Terminal**: default profile `bash`.
- **Git scan caps**: `git.detectSubmodules=false`,
  `git.autoRepositoryDetection=openEditors`, `git.repositoryScanMaxDepth=1`.
- **Watcher/search excludes**: `files.watcherExclude` drops `00_data`,
  `01_modules`, `.venv`, `*-env`, `__pycache__`, `.ipynb_checkpoints`,
  `renv/library`, `03_results/checkpoints` — but keeps `03_results` itself
  watched so figure PNGs live-refresh.

The git and watcher caps exist to prevent Remote-SSH / Dev-Container reconnect
drops on large umbrella workspaces — details in
[vscode-remote-stability.md](vscode-remote-stability.md).

## See also

- [architecture.md](architecture.md) — image layering, UID remapping, ArchR status
- [ai-integration.md](ai-integration.md) — `setup_ai_env.sh`, AI CLIs, SciAgent-toolkit boundary
- [isolation.md](isolation.md) — hardening + secrets patterns
- [ssh-passthrough.md](ssh-passthrough.md) — SSH agent forwarding
- [vscode-remote-stability.md](vscode-remote-stability.md) — watcher/git tuning
- [operations.md](operations.md) — running containers, sanity checks
