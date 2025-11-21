# Repository Guidelines

## Project Structure & Module Organization
- `.devcontainer/`: Dockerfiles, compose config, R/Python installers, and sanity scripts. Primary builds use `Dockerfile.optimized`; optional AI variant uses `Dockerfile.ai-enabled` (via SciAgent‑toolkit).
- `.environments/`: Python requirements for layered venvs (e.g., `base_requirements.txt`).
- `docs/`: Operations and AI tooling docs.
- `templates/`: Project scaffolds (`basic-rna`, `multimodal`, `archr-focused`, `example-DMATAC`).
- Root scripts: `build-optimized.sh`, `build-ai-enabled.sh`, `init-project.sh`.
- Version pins: `renv.lock`, `R-packages-manifest.csv`.

## Build, Test, and Development Commands
- Build optimized base image: `./build-optimized.sh [--github-pat <token>] [--tag scdock-r-dev:vX]`.
- Build AI-enabled image (optional): `./build-ai-enabled.sh [--with-tooluniverse] [--with-serena]`.
- Scaffold a project: `./init-project.sh ~/projects/my-analysis basic-rna`.
- Start via compose: `docker compose -f .devcontainer/docker-compose.yml up -d dev-core` (or `dev-archr`).
- Sanity checks (inside container): `.devcontainer/scripts/poststart_sanity.sh` → expect “OK” lines for Python, R, httpgd, venvs.

## Coding Style & Naming Conventions
- Shell: bash with `set -euo pipefail`; 2-space indent; functions `lower_snake_case`; scripts kebab-case (e.g., `build-archr-wrapper.sh`). Env vars `UPPER_SNAKE_CASE`.
- Docker/compose: keep args explicit (UID/GID, PAT); pin tags.
- R/Python in templates: follow tidyverse/PEP8 defaults; keep notebooks and scripts under `01_Scripts/` or `02_Analysis/` in projects.

## Testing Guidelines
- No unit tests; use smoke tests.
- After build/run, verify:
  - `radian` starts; `R -q -e "packageVersion('httpgd')"` prints a version.
  - `python -c "import scanpy"` exits 0; `which python` under `/opt/venvs/base` or selected venv.
  - Run `.devcontainer/scripts/poststart_sanity.sh` and confirm all OK.

## Commit & Pull Request Guidelines
- Commits: imperative, concise subjects; include scope or version when relevant (e.g., `build: multi-stage size cut (v0.5.1)`, `docs: README corrections`).
- PRs must include: summary, rationale, commands to reproduce (build/run), sanity-check output, and links to issues. Update `README.md`/`DEVOPS.md` when behavior changes.

## Security & Configuration Tips
- Do not commit secrets. Pass `GITHUB_PAT` via env or flags.
- Pin R via `renv.lock`; pin Python in `.environments/*.txt`. Avoid unpinned base image bumps without verification.
- Prefer the optimized multi-stage build for size and reproducibility.

## Branching Policy
- Active development occurs on `dev`; stable snapshots live on `main`.
- AI branches `dev-claude-integration` and `dev-gpt-codex-integration` are deprecated.
- AI integration is provided via the optional AI-enabled build and the external SciAgent‑toolkit. See `docs/AI_TOOLS.md`.
