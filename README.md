# Single-Cell Docker Dev Environment

![Docker Image Version](https://img.shields.io/badge/Docker-v0.5.2-blue?style=flat-square)
![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)

Purpose-built Docker images and VS Code Dev Container config for single-cell analysis in R and Python. The goal is a clean, reproducible, and fast-to-start environment you can use locally or remotely without yak-shaving.

Who this helps
- Beginners who want a working R/Python stack without setup pain
- Researchers who switch across machines or share a consistent image with teammates
- Future me: predictable builds, minimal cognitive load, explicit docs

Key features
- ~20GB true image via multi-stage build (size-optimized)
- R 4.5 + Bioconductor 3.21 core stack; Python 3.10 base venv
- Layered Python venvs on demand: squid (spatial), atac, comms
- VS Code friendly: httpgd plotting, radian, language server
- Build tools retained to allow runtime installs when needed

Quick start
```bash
# 1) Build the base image
scripts/build.sh

# 2) (Optional) Pull official ArchR image for scATAC work
docker pull greenleaflab/archr:1.0.3-base-r4.4

# 3) Scaffold a new analysis project
./init-project.sh ~/projects/my-analysis basic-rna --interactive

# 4) Open + Reopen in Container
code ~/projects/my-analysis
```

Working in the container
- Default service: dev-core (R 4.5 + Python)
- Switch to ArchR (R 4.4) by setting service to dev-archr in `.devcontainer/devcontainer.json`
- R: run `radian` (or `r-base` wrapper)
- Python envs: `usepy base|squid|atac|comms` (creates layered venvs on demand)
- Sanity check: `.devcontainer/scripts/poststart_sanity.sh` (in a project) or run `scripts/poststart_sanity.sh` inside the image with a bind mount

Build modes (brief)
- Generic (default): shareable image with `devuser:1000`
- Personal: `scripts/build.sh --personal` bakes your UID/GID; not shareable but handy for local use

Documentation
- Quick start: QUICK-START.md
- Build guide: docs/build.md
- Architecture: docs/architecture.md
- DevOps and operations: docs/devops.md
- Runtime R/Python installs: docs/runtime-install.md
- Repo structure: docs/repo-structure.md
- Image size notes: docs/size-optimization.md
- Migration notes: docs/migration.md
- Branching model: docs/branching.md
- Changelog: docs/changelog.md

AI integration (CLI agents)
- This repository can be integrated with the `SciAgent-toolkit` to provide AI-powered assistance for your analysis.
- The `SciAgent-toolkit` is included as a Git submodule and can be installed into a dedicated Docker image.

To build the AI-enabled image, run the following command:
```bash
# Initialize the submodule (only needs to be done once)
git submodule update --init --recursive

# Build the AI-enabled image
./build-ai-enabled.sh
```
This will create a new Docker image tagged `scdock-ai-dev:v0.5.2` with the `SciAgent-toolkit` pre-installed.

To start a project using the AI-enabled image:
```bash
./init-project.sh ~/projects/my-analysis basic-rna --ai --interactive
```

License
- MIT. See LICENSE.

Acknowledgments
- Inspired by [Rami Krispin’s vscode-r](https://github.com/RamiKrispin/vscode-r) for VS Code + R workflow ideas.
