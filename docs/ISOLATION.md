# Filesystem Isolation Patterns

scbio-docker borrows selective filesystem-isolation patterns from
[pi-coding-agent-container](https://github.com/) (a hardened agent runtime),
adapted for **interactive** scientific work. The defaults are conservative;
optional opt-in hardening is documented inline in the compose template.

Note on Docker Compose schemas: we set `pids` under
`deploy.resources.limits.pids` rather than as a top-level `pids_limit`. Recent
Compose versions reject a project that sets both, and treat an absent
`deploy.resources.limits.pids` as `0`, which conflicts with a top-level
`pids_limit`. Putting it inside `deploy.resources.limits` avoids the clash.

## Defaults (always on, transparent to users)

| Pattern | Why | Where |
|---------|-----|-------|
| `tmpfs /tmp` with `noexec,nodev,nosuid` | Stops downloaded binaries from running out of `/tmp`; standard tools don't need exec there. | Both `dev-core` and `dev-archr` services. |
| `deploy.resources.limits.pids: 4096` | Caps fork bombs; high enough that R `BiocParallel` and Python `multiprocessing` keep working. | Both services. |
| Host UID/GID mapping (`${LOCAL_UID:-1000}:${LOCAL_GID:-1000}`) | Files written from the container land owned by you, not root. | Both services (already present pre-isolation work). |

## Opt-in hardening (uncomment in docker-compose.yml)

- `read_only: true` root filesystem with writable tmpfs for `/tmp` and
  `/home/devuser`. Useful for shared/registry-pushed images. Breaks anything
  that writes outside bind-mounted workspaces.
- `cap_drop: [ALL]` + minimal `cap_add` set. Tighter than Docker's default. May
  break `apt-get install` or system-package compilation done inside a running
  container.
- `security_opt: ["no-new-privileges:true"]` — prevents `setuid` binaries from
  escalating. Safe to keep on if you don't need `sudo` inside the container.

## Secrets pattern

Avoid baking API keys into `.env`. Instead:

1. Store secrets outside the repo, e.g. `~/.scbio-secrets/anthropic_api_key`,
   `chmod 600`.
2. Wire them as Docker secrets in compose (commented example at the bottom of
   `templates/base/.devcontainer/docker-compose.yml.template`).
3. Apps read them from `/run/secrets/<name>`, mounted with mode `0000` (only
   accessible via container processes that have the right effective UID).

## Explicitly NOT adopted

| Pattern | Why skipped |
|---------|-------------|
| `LD_PRELOAD` filesystem firewall (pi-coding-agent's `fs-vault.so`) | Too aggressive for interactive use — breaks `htop`, `strace`, debuggers, and any tool that legitimately reads config. |
| SetUID binary purge (`chmod a-s` on entire FS) | Breaks `sudo`-based workflows for occasional package installs (which scientists do). |
| `gh-guard.sh` / command-level wrappers | Loss of CLI flexibility; mismatched threat model (you are the user, not an autonomous agent). |
| Application-level firewalls (`NODE_OPTIONS` monkeypatch) | R/Python image; would only cover Node tooling and break legitimate fs reads. |

## When to enable opt-in hardening

- **Pushing the image to a shared registry** other people will pull → enable
  `read_only: true`, `cap_drop`, `no-new-privileges`.
- **Running an autonomous agent (Claude Code, Aider, etc.) inside the container
  with elevated tool access** → all of the above + secrets pattern.
- **Personal interactive use on your own machine** → defaults are sufficient.
