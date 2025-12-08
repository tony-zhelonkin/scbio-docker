Build Guide (dev-restructure)

- Preferred: scripts/build.sh
  - Supports PAT via --github-pat
  - Supports UID/GID remap via --personal or --user-id/--group-id
- Direct docker build:
  - docker build . -f docker/base/Dockerfile -t scdock-r-dev:latest

Outputs
- Image tag: configurable via --tag
- Logs: build.log in repo root when using scripts/build.sh

Notes
- R installers are under docker/base/R/
- Python requirements are under docker/requirements/

