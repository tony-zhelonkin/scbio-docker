# Documentation

Detailed docs for **scbio-docker** — the Docker/dev-container substrate for
single-cell R + Python analysis. The canonical image version lives in the
[`VERSION`](../VERSION) file; version history lives in [changelog.md](changelog.md).

Start at the root [README.md](../README.md) and [QUICKSTART.md](../QUICKSTART.md).

## Map

| Doc | What it covers |
|-----|----------------|
| [architecture.md](architecture.md) | Image layering, venv/R-library model, design decisions, image size |
| [build.md](build.md) | Building the image: `scripts/build.sh`, build modes, UID/GID, renv lockfile |
| [environments.md](environments.md) | Python venvs (`usepy`), R two-tier libraries, runtime installs, R↔Python interop |
| [devcontainer.md](devcontainer.md) | `init-container.sh`, the `templates/devcontainer/` scaffold, compose services, VS Code settings |
| [ai-integration.md](ai-integration.md) | Containerization-only stance, SciAgent-toolkit boundary, `setup_ai_env.sh` AI-CLI bootstrap |
| [operations.md](operations.md) | Runbook: run/compose, ArchR (legacy), tmux R sessions, troubleshooting |
| [repo-structure.md](repo-structure.md) | Repository tree and where things live |
| [isolation.md](isolation.md) | Filesystem/process isolation and opt-in container hardening |
| [ssh-passthrough.md](ssh-passthrough.md) | SSH agent forwarding into the container |
| [vscode-remote-stability.md](vscode-remote-stability.md) | Avoiding Remote-SSH / Dev-Container reconnect freezes |
| [branching.md](branching.md) | Git branching and archive tags |
| [roadmap.md](roadmap.md) | Forward-looking direction |
| [changelog.md](changelog.md) | Cumulative version history (the one place versions live) |

Community health: [CONTRIBUTING.md](CONTRIBUTING.md) · [SECURITY.md](SECURITY.md) · [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md)

Historical planning snapshots are frozen under [`archive/`](archive/).

## Reading order

- **New users:** [../README.md](../README.md) → [../QUICKSTART.md](../QUICKSTART.md) → [architecture.md](architecture.md) → [operations.md](operations.md)
- **Building the image:** [build.md](build.md) → [architecture.md](architecture.md) → [environments.md](environments.md)
- **Setting up a project:** [devcontainer.md](devcontainer.md) → [ai-integration.md](ai-integration.md)
- **Contributors:** [CONTRIBUTING.md](CONTRIBUTING.md) → [repo-structure.md](repo-structure.md) → [branching.md](branching.md)

## Conventions

- **Version-agnostic** — commands read the tag from `VERSION`; explicit versions live only in [changelog.md](changelog.md).
- **One concern per doc** — docs link out rather than duplicate.
- **Scannable** — short intro, then tables and fenced commands.
