Agent Guide

Scope
- This file applies to the entire repository.

Where to look first
- docker/base/Dockerfile — canonical build definition
- scripts/build.sh — build wrapper (supports PAT and UID/GID remap)
- .devcontainer/devcontainer.json — VS Code container config (uses prebuilt images)
- scripts/init-project.sh — scaffolds new analysis projects

Conventions
- Keep Docker build assets under docker/ and host-side scripts under scripts/.
- Do not add new build logic under .devcontainer/; keep it configuration-only.
- Prefer adding docs under docs/ and link from README.md.
- When moving files, update references across README.md and docs/.

Editing guidelines
- Minimal diffs; match existing style.
- Prefer small, reviewable PRs (scaffold first, then switch references, then delete legacy files).
- Validate with scripts/build.sh and run scripts/poststart_sanity.sh in a container.

Testing
- Smoke tests: ensure R (Seurat/edgeR/httpgd) and Python (scanpy/scvi) import.
- Sanity script location for built images: scripts/poststart_sanity.sh.