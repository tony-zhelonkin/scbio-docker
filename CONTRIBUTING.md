Contributing

Thank you for your interest in contributing!

Workflow
- Branch from dev; open PRs into dev. Releases merge from dev → main.
- Keep changes focused and small where possible.
- Run a local build and sanity check before opening a PR.

Local build
- scripts/build.sh --tag scdock-r-dev:dev-local
- docker run --rm scdock-r-dev:dev-local bash -lc 'scripts/poststart_sanity.sh'

Structure
- docker/ — Dockerfiles and installers
- scripts/ — host-side scripts (build, init-project, sanity)
- .devcontainer/ — VS Code config only
- docs/ — long-form documentation

Coding conventions
- Shell: set -euo pipefail; POSIX-ish; readable functions.
- Dockerfile: multi-stage; copy only what’s needed; minimize layers; keep runtime toolchain for package installs.
- R installers: avoid silent failures; prefer installed.packages() checks for meta-packages.

Issues and PRs
- Use a clear title, short summary, and checklist of changed areas.
- Link to any related issues or design notes.

