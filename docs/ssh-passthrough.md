# SSH passthrough — deferred

**Status (as of v0.5.4):** *Not implemented.* Tracked here as a known gap with a copy-pasteable recipe to enable when it lands.

## What's missing

| Layer | Current state |
|---|---|
| Image (`scdock-r-dev:v0.5.4`) | `git` is present. `ssh` and `gh` CLIs are **not** installed. |
| Compose template (`templates/base/.devcontainer/docker-compose.yml.template`) | No `SSH_AUTH_SOCK` env, no agent-socket mount, no `~/.ssh` mount, no `~/.config/gh` mount. |
| Existing project devcontainers (DC_hum_verse, AdaW_eWAT, GVDRP1, JBader, etc.) | None forward SSH or `gh` config. The pattern was never added. |

Practical impact inside any current container:
- `ssh git@github.com` → "ssh: command not found".
- `gh ...` → "gh: command not found".
- `git push` over an SSH remote → fails (no client, no key access).
- Workaround in current use: do everything SSH-/gh-related on the host, treat the container as a compute layer.

## Why it's deferred (not blocked)

- It widens the trust boundary of the container slightly (any in-container process can use the host's loaded SSH keys via the agent socket). For an interactive single-user container that's fine, but it's a deliberate decision to flip on.
- It requires an image rebuild (apt installs) plus template edits — should bundle with the next minor bump rather than amend v0.5.4 in place.
- No active project currently needs it; the host workflow has been sufficient.

## Recipe to enable (when it lands — target: v0.5.5 or v0.6.0)

### 1. Image: install SSH client and `gh` CLI

In `docker/base/Dockerfile`, add to the apt installs in the runtime stage:

```dockerfile
RUN apt-get update && apt-get install -y --no-install-recommends \
    openssh-client \
 && curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg \
        | gpg --dearmor -o /usr/share/keyrings/githubcli-archive-keyring.gpg \
 && chmod go+r /usr/share/keyrings/githubcli-archive-keyring.gpg \
 && echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" \
        | tee /etc/apt/sources.list.d/github-cli.list > /dev/null \
 && apt-get update && apt-get install -y --no-install-recommends gh \
 && rm -rf /var/lib/apt/lists/*
```

Adds ~30 MB. Bumps version to v0.5.5.

### 2. Compose template: forward SSH agent + (optionally) `gh` config

Edit `templates/base/.devcontainer/docker-compose.yml.template`. For BOTH `dev-core` and `dev-archr` services, add to the `environment:` and `volumes:` blocks:

```yaml
    environment:
      - SSH_AUTH_SOCK=/ssh-agent

    volumes:
      - ${WORKSPACE_FOLDER:-.}:/workspaces/{{PROJECT_NAME}}
      # SSH agent forwarding (falls back to /dev/null when no agent runs)
      - ${SSH_AUTH_SOCK:-/dev/null}:/ssh-agent
      # known_hosts so first-time GitHub connections don't prompt
      - ${HOME}/.ssh/known_hosts:/home/devuser/.ssh/known_hosts:ro
      # Optional: reuse host's gh CLI auth (mounted read-only)
      - ${HOME}/.config/gh:/home/devuser/.config/gh:ro
{{DATA_MOUNTS}}
```

The `${SSH_AUTH_SOCK:-/dev/null}` fallback is the standard pattern — when no agent is running on the host (CI, headless), the bind mount lands on `/dev/null` and SSH simply isn't available; nothing else breaks.

### 3. (No code change) Confirm the host has an SSH agent loaded

```bash
ssh-add -l                          # should list your GitHub key
gh auth status                      # should report logged-in
```

If `ssh-add -l` returns "Could not open a connection to your authentication agent", the user has not started one. Common patterns: `eval "$(ssh-agent)" && ssh-add ~/.ssh/id_ed25519_github`, or use `keychain` (the user's host already has `keychain 2.8.5` set up — see the `gh auth login` flow that ran earlier).

### 4. Rebuild and test

```bash
echo v0.5.5 > VERSION
scripts/build.sh --yes --github-pat "$(gh auth token)"
docker run --rm -v ${SSH_AUTH_SOCK}:/ssh-agent -e SSH_AUTH_SOCK=/ssh-agent \
    scdock-r-dev:v0.5.5 bash -lc 'ssh -T git@github.com || true; gh --version'
```

Expected: `Hi <username>! You've successfully authenticated...` from GitHub, and a `gh version 2.x.x` line.

## Decision log

- **Forward `~/.config/gh`?** *Optional.* Convenient (no re-auth in container) but mounts a token file. Keep it `:ro`. Skip if you prefer per-container `gh auth login`.
- **Mount `~/.ssh` directly?** *Don't.* It exposes private key files to the container filesystem. The agent-socket pattern keeps keys on the host while letting the container *use* them via the agent.
- **Why not VS Code's automatic agent forwarding?** VS Code does forward the agent automatically when launching a Dev Container — but only when the container is launched via VS Code, not via raw `docker compose up`. Explicit compose-level forwarding works in both cases.

## Cross-references

- Compose template that needs editing: `templates/base/.devcontainer/docker-compose.yml.template`
- Image Dockerfile: `docker/base/Dockerfile`
- Isolation policy (compatible — agent forwarding does not conflict with current `tmpfs /tmp` / `pids` defaults): [ISOLATION.md](ISOLATION.md)
- Roadmap entry: [ROADMAP.md](ROADMAP.md) under v0.6.0 → "SSH passthrough"
