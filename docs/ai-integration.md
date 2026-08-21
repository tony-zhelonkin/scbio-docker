# AI integration

The scbio-docker image is **containerization-only**: it bakes in no AI CLIs, no
MCP servers, and no agents. All AI tooling is installed at runtime, per-project,
by a script that ships with the devcontainer template. The image carries only
the **prerequisites** that downstream tooling needs.

| Prereq | Why |
|--------|-----|
| Node.js 20 LTS (`node`, `npm`, `npx`) | JS/TS CLIs + npm fallbacks |
| `uv` / `uvx` (Astral) | Python-based AI tools |
| Python `toml` (in `/opt/venvs/base`) | Codex CLI config generation |
| `jq` | settings deep-merge + Pi model config |
| `curl` | official installer scripts |
| `bubblewrap` (`bwrap`) | sandbox backend for Codex `workspace-write` mode |

### Codex sandbox (bubblewrap) requirement

Codex's `workspace-write` mode jails each exec with **bubblewrap** (`bwrap`),
which is baked into the image. `bwrap` sets up an unprivileged user namespace,
and Docker's default seccomp profile blocks the `clone`/`unshare` calls it
needs — so inside a stock container `bwrap` fails with *Operation not
permitted* even though the binary is present. Run the container with
`--security-opt seccomp=unconfined` (or an equivalent host userns-enabled
setup) for Codex sandboxing to work. The compose template ships this as a
commented `security_opt: ["seccomp=unconfined"]` line on the dev service —
uncomment it; add `--security-opt seccomp=unconfined` to ad-hoc `docker run`
invocations that use Codex.

## Boundary: scbio-docker vs SciAgent-toolkit

scbio-docker is a **container substrate** and nothing more. The AI *harness*
(agents, skills, methodology guidelines, project catalog) belongs to
**SciAgent-toolkit**, a sibling repo attached per-project — never vendored into
the image.

| Repository | Owns |
|------------|------|
| **scbio-docker** | Dockerfiles, image/env specs, build scripts, devcontainer + compose templates, `init-container.sh`, `.env` stub, `setup_ai_env.sh` |
| **SciAgent-toolkit** | Project tree and catalog, config templates, docs namespaces, AI harness (skills/agents/commands), methodology guidelines |

**Rule of ownership:** does the content change when the **image** changes
(scbio-docker) or when the **project / AI harness** changes (SciAgent-toolkit)?

SciAgent-toolkit is tracked here as a submodule at `toolkits/SciAgent-toolkit/`
and re-attached per-project at `01_modules/SciAgent-toolkit/`. It is **not**
copied into the image.

## Two-step workflow

```bash
# 1. Render the dev container into a target directory (scbio-docker)
./init-project.sh ~/projects/my-analysis \
    --data-mount atac:/scratch/data/DT-1234 \
    --service dev-core --max-cpus 50 --max-memory 450G

# 2. Bind the toolkit catalog and render CRAFT (SciAgent-toolkit)
cd ~/projects/my-analysis
./01_modules/SciAgent-toolkit/bin/scio link
./01_modules/SciAgent-toolkit/bin/scio craft

# 3. Open in VS Code, reopen in container
code ~/projects/my-analysis
```

Once inside the container, run the toolkit's setup once:

```bash
./01_modules/SciAgent-toolkit/scripts/setup-ai.sh
```

`init-project.sh` is a symlink to `scripts/init-container.sh`. See
[devcontainer.md](devcontainer.md) for the container rendering details.

## setup_ai_env.sh

`.devcontainer/scripts/setup_ai_env.sh` is rendered into every project by
`init-container.sh` and invoked from `devcontainer.json`'s `postStartCommand`.
It runs on **every container start** (the container home dir is not persisted
across rebuilds, so it re-establishes the toolchain each time).

Design contract:

- **Idempotent** — skips any CLI already on `PATH`; safe to run every start.
- **Network-tolerant** — a failed installer only WARNs; it never blocks startup
  (the script uses `set -uo pipefail`, not `-e`).

It prepends the installer target dirs to `PATH` for the run:

```
$HOME/.local/bin : $HOME/.opencode/bin : $HOME/.npm-global/bin : $PATH
```

### AI CLIs installed

Five CLIs, each skipped if already present. Curl installers hit GitHub's release
API (rate-limited to 60/hr per IP on a shared host), so most have an npm
fallback that is not rate-limited.

| CLI | Probe | Primary installer | Fallback |
|-----|-------|-------------------|----------|
| Claude Code | `claude` | `curl -fsSL https://claude.ai/install.sh \| bash` | — |
| OpenAI Codex | `codex` | `curl -fsSL https://chatgpt.com/codex/install.sh \| sh` | npm `@openai/codex` |
| opencode | `opencode` | `curl -fsSL https://opencode.ai/install \| bash` | npm `opencode-ai` |
| Pi agent | `pi` | npm `@mariozechner/pi-coding-agent` (prefix `$HOME/.local`) | — |
| Antigravity | `antigravity` / `agy` | `curl -fsSL https://antigravity.google/cli/install.sh \| bash` | — |

Notes:
- Pi is installed via npm directly: the official pi.dev installer is interactive
  and wants Node ≥ 22.19, while the npm package installs headless on Node 20.
- Antigravity ships as `antigravity` with the short alias `agy`; both are probed.

### Pi local-model config

`configure_pi_models()` reads `models.agentic.json` (shipped beside the script)
and writes `~/.pi/agent/models.json` via `jq`, rewriting the Ollama `baseUrl` to
`${OLLAMA_HOST:-http://172.17.0.1:11434}/v1` (the Docker bridge). If `jq` is
absent it writes the file verbatim (baseUrl not rewritten).

`models.agentic.json` is the single source of truth for the Pi roster. It lists
a **single Ollama provider** with **tool-capable (agentic) models only** —
reasoning-only models (FuseO1, DeepSeek-R1) are deliberately excluded because
they do not tool-call and 400 the agent. The `_comment` / `_policy` / `_bridge`
metadata keys are stripped on write.

### Claude power-user config (`~/.claude`, USER scope)

Seeded on every start; pre-existing keys are preserved via `jq` deep-merge with
the desired values winning on conflict.

- **`statusline.sh`** — renders `[model] [progress-bar] tokens (N% free)` from
  the real statusline payload.
- **`settings.json`** — `editorMode: vim`, `alwaysThinkingEnabled: true`,
  `effortLevel: xhigh`, `showThinkingSummaries: true`, `teammateMode: auto`,
  `autoMemoryEnabled: false`, `cleanupPeriodDays: 90`, and
  `env.CLAUDE_CODE_EXPERIMENTAL_AGENT_TEAMS: "1"`.
- **Attribution suppressed** — `attribution: { commit: "", pr: "", sessionUrl: false }`
  drops the `Co-Authored-By: Claude` commit trailer, the "Generated with Claude
  Code" PR footer, and the session-URL trailer, per SciAgent's no-AI-authorship
  policy.
- **Telemetry left ON** — disabling it trips the feature-flag layer that gates
  agent teams / 1M context, so `DISABLE_TELEMETRY` is intentionally not set.

### Shell-rc conveniences

Two marker-guarded blocks appended to `~/.bashrc` (each at most once per
container lifetime):

- **`si` alias** → `./01_modules/SciAgent-toolkit/bin/scio`. CWD-relative
  (not a PATH symlink) so it resolves to whichever project's vendored toolkit
  copy you `cd` into.
- **`_source_project_env`** — auto-sources the active project's
  `.devcontainer/.env` (and `.env`) into every interactive shell, loading MCP
  API-key stubs (`GEMINI_API_KEY`, `OPENAI_API_KEY`) and `OLLAMA_HOST`.

## See also

- [devcontainer.md](devcontainer.md) — template rendering and compose services
- [architecture.md](architecture.md) — image layering
- [repo-structure.md](repo-structure.md) — where these files live
