# Single-Cell Docker Dev Environment

![Docker Image Version](https://img.shields.io/badge/Docker-v0.5.1-blue?style=flat-square)
![License](https://img.shields.io/badge/License-MIT-green?style=flat-square)

A Docker-based single-cell bioinformatics environment for VS Code Dev Containers. Personal convenience toolkit shared publicly—use what’s useful and ignore the rest.

> Inspiration: strongly inspired by [Rami Krispin’s vscode-r repository](https://github.com/RamiKrispin/vscode-r). It was a great reference for R + VS Code setup.

## Why use it
- ~20GB true image via multi-stage build
- R 4.5 + Bioc 3.21 core stack, Python 3.10 base venv
- Layered venvs for spatial/ATAC/communication toolsets
- Build tools kept in runtime for on-demand installs
- Devcontainer-friendly, works locally/remotely

## Quick start
```bash
# Build base image
scripts/build.sh

# (optional) pull official ArchR image
docker pull greenleaflab/archr:1.0.3-base-r4.4

# Scaffold a project
./init-project.sh ~/projects/my-analysis basic-rna

# Open in VS Code and Reopen in Container
code ~/projects/my-analysis
```

## Features
- R: Seurat, edgeR, limma, clusterProfiler, GSVA, Signac, and friends
- Python: scanpy, scvi-tools, muon, squidpy (via layered venv), etc.
- CLI: samtools, bcftools, bedtools; scIBD
- TinyTeX + httpgd for plotting in VS Code

## Docs
- Quick start guide: QUICK-START.md
- Build: docs/build.md
- Architecture: docs/architecture.md
- DevOps details: DEVOPS.md
- Runtime package install: RUNTIME_INSTALL.md
- Repo layout: docs/repo-structure.md
- Image size notes: SIZE_OPTIMIZATION_SUMMARY.md
- Migration (restructure): docs/migration.md

## Contributing
Informal and lightweight. Open a PR or issue if you have a concrete improvement. This is optimized for personal use.

## License
MIT. See LICENSE.

## Acknowledgments
- Thanks to [Rami Krispin’s vscode-r repo](https://github.com/RamiKrispin/vscode-r) for inspiration.

