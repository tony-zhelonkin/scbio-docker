# Architecture (short overview)

Images
- Base image: built from docker/base/Dockerfile (multi-stage).
- ArchR: use official greenleaflab/archr:1.0.3-base-r4.4 when needed.

Dev Containers
- .devcontainer/devcontainer.json uses docker-compose to pick a service (dev-core or dev-archr).
- Both services share the same mounts and UID/GID passthrough.

Python envs
- Base venv at /opt/venvs/base (preinstalled).
- Layered venvs (squid, atac, comms) created on first use with system-site-packages.

R setup
- Core packages preinstalled during build via docker/base/R/install_core.R and renv snapshot.
- httpgd enabled for VS Code plotting via .devcontainer/.Rprofile.

Sanity
- scripts/poststart_sanity.sh checks R/Python/httpgd and venv paths.

That’s the gist—see DEVOPS.md for the longer version.

