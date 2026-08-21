# Quick Start Guide

`scbio-docker` is a container substrate: it builds the Docker image and renders a
VS Code dev container into a project directory. The project catalog and AI harness
come from [SciAgent-toolkit](https://github.com/tony-zhelonkin/SciAgent-toolkit),
vendored at `01_modules/SciAgent-toolkit/`.

## 1. Build the image

```bash
scripts/build.sh                      # generic build (devuser:1000)
docker pull greenleaflab/archr:1.0.3-base-r4.4   # optional, for scATAC (ArchR)
```

## 2. Render a dev container into a project directory

```bash
./init-project.sh ~/projects/my-analysis \
    --data-mount atac:/scratch/data/DT-1234 \
    --data-mount rna:/scratch/data/DT-5678:ro
```

This writes only `.devcontainer/{devcontainer.json,docker-compose.yml,.env,scripts/}`.
See `init-project.sh --help` for `--service`, `--image-version`, `--max-cpus`, `--max-memory`.

## 3. Bind the toolkit catalog and render CRAFT (SciAgent-toolkit)

```bash
cd ~/projects/my-analysis
./01_modules/SciAgent-toolkit/bin/scio link
./01_modules/SciAgent-toolkit/bin/scio craft
```

## 4. Open in VS Code

```bash
code ~/projects/my-analysis
# Ctrl+Shift+P -> "Dev Containers: Reopen in Container"
```
