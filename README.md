# Single-Cell Docker Dev Environment

![Docker Image Version](https://img.shields.io/badge/Docker-v0.5.10-blue?style=flat-square)
![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)

Purpose-built Docker images and VS Code Dev Container config for single-cell analysis in R and Python. The goal is a clean, reproducible, and fast-to-start environment you can use locally or remotely without yak-shaving.

Who this helps
- Beginners who want a working R/Python stack without setup pain
- Researchers who switch across machines or share a consistent image with teammates
- Future me: predictable builds, minimal cognitive load, explicit docs

Key features
- ~20GB true image via multi-stage build (size-optimized)
- R 4.5.3 + Bioconductor 3.22 core stack (~80 packages); Python 3.11 base venv
- Layered Python venvs on demand: squid (spatial), atac, comms
- VS Code friendly: httpgd plotting, radian, language server, Jupyter R + Python kernels
- Containerization-only: AI tooling installs at runtime (see docs/ai-integration.md)
- Build tools retained to allow runtime installs when needed

Quick start
```bash
# 1) Build the base image
scripts/build.sh

# 2) (Optional) Pull official ArchR image for scATAC work
docker pull greenleaflab/archr:1.0.3-base-r4.4

# 3) Render a dev container into a project directory
./init-project.sh ~/projects/my-analysis

# 4) Scaffold the project structure (SciAgent-toolkit)
sciagent new project --type analysis ~/projects/my-analysis

# 5) Open + Reopen in Container
code ~/projects/my-analysis
```

Working in the container
- Default service: dev-core (R 4.5.3 + Python 3.11)
- ArchR (R 4.4) is a legacy sidecar, gated behind the `archr` compose profile — on the path to removal (see docs/architecture.md)
- R: run `radian` (or `r-base` wrapper)
- Python envs: `usepy base|squid|atac|comms` (creates layered venvs on demand)
- Sanity check: `.devcontainer/scripts/poststart_sanity.sh` (in a project) or run `scripts/poststart_sanity.sh` inside the image with a bind mount

Build modes (brief)
- Generic (default): shareable image with `devuser:1000`
- Personal: `scripts/build.sh --personal` bakes your UID/GID; not shareable but handy for local use

Documentation — see [docs/README.md](docs/README.md) for the full map.
- Quick start: QUICKSTART.md
- Architecture: docs/architecture.md
- Build guide: docs/build.md
- Environments (R/Python, runtime installs): docs/environments.md
- Dev container & templates: docs/devcontainer.md
- AI integration: docs/ai-integration.md
- Operations runbook: docs/operations.md
- Repo structure: docs/repo-structure.md
- Changelog: docs/changelog.md

AI integration (CLI agents)

The image is **containerization-only**. 
All context/LLM/agent management lives in a dedicated repo — [SciAgent-toolkit](https://github.com/tony-zhelonkin/SciAgent-toolkit).
The image carries only the **prerequisites** (Node.js 20, `uv`/`uvx`, Python `toml`) so
that SciAgent-toolkit can install AI tooling at runtime, per-project.

Thus my personal workflow is: 
- scbio-docker renders the dev container  
- SciAgent-toolkit scaffolds the project context management wrapping AI harness. 

They are intended to compose 

```bash
# 1. Render the dev container into a project directory (scbio-docker)
./init-project.sh ~/projects/my-analysis

# 2. Scaffold the project structure + AI harness (SciAgent-toolkit)
sciagent new project --type analysis ~/projects/my-analysis

# 3. Open in VS Code, Reopen in Container, then run AI setup (first time only)
./01_modules/SciAgent-toolkit/scripts/setup-ai.sh
```

License
- MIT. See LICENSE.

Acknowledgments
- Inspired by [Rami Krispin’s vscode-r](https://github.com/RamiKrispin/vscode-r) for VS Code + R workflow ideas.
