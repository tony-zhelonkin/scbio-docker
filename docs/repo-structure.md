Repository Structure (dev-restructure)

- docker/
  - base/
    - Dockerfile — primary multi-stage build
    - R/
      - install_core.R — core R stack installer
      - install_renv_project.R — renv init/restore wrapper
      - install_httpgd.R — CRAN-first + GitHub fallback
  - requirements/
    - base.txt, squid.txt, atac.txt, comms.txt — Python stacks
- scripts/
  - build.sh — wrapper to build docker/base/Dockerfile
  - create_layered_venv.sh — runtime layered venv helper
  - init-project.sh — wrapper (delegates to repo-root init-project.sh)
  - poststart_sanity.sh — sanity checks (copied into projects)
- .devcontainer/
  - devcontainer.json, docker-compose.yml — config only
- templates/ — project templates
- docs/ — consolidated documentation (core topics)
  - build.md, architecture.md, migration.md, repo-structure.md
  - additional docs at repo root: QUICK-START.md, DEVOPS.md, RUNTIME_INSTALL.md, SIZE_OPTIMIZATION_SUMMARY.md, HANDOFF.md, BRANCH_MANAGEMENT.md
- README.md — short overview + links into docs
- CONTRIBUTING.md, CODE_OF_CONDUCT.md, AGENTS.md — contribution guidance

Notes
- Old locations remain during the transition window. Prefer new paths in docs and scripts.
- Image builds should use scripts/build.sh or docker -f docker/base/Dockerfile.
