# Filesystem Isolation

scbio-docker applies selective filesystem-isolation patterns adapted for
**interactive** scientific work. The defaults are conservative and always on;
optional opt-in hardening ships commented out in the compose template. All of it
lives in `templates/devcontainer/.devcontainer/docker-compose.yml.template`,
which [`init-container.sh`](devcontainer.md) renders into a project's
`.devcontainer/docker-compose.yml`.

## Compose schema note

The process cap is set as `deploy.resources.limits.pids` rather than a top-level
`pids_limit`. Recent Compose versions reject a project that sets both, and treat
an absent `deploy.resources.limits.pids` as `0`, which conflicts with a
top-level `pids_limit`. Keeping it inside `deploy.resources.limits` avoids the
clash.

## Defaults (always on)

| Pattern | Why | Where |
|---------|-----|-------|
| `tmpfs /tmp` → `noexec,nodev,nosuid,mode=1777` | Stops binaries downloaded to `/tmp` from executing; standard tooling doesn't need exec there. | Both `dev-core` and `dev-archr`. |
| `deploy.resources.limits.pids: 4096` | Caps fork bombs; high enough that R `BiocParallel` and Python `multiprocessing` keep working. | Both services. |
| Host UID/GID mapping (`${LOCAL_UID:-1000}:${LOCAL_GID:-1000}`) | Files written from the container land owned by you, not root. | Both services. |

## Opt-in hardening

These blocks are present but **commented out** in the compose template (between
the `--- Optional hardening ---` markers on each service). Enable by deleting the
leading `# ` on the lines you want, then rebuild/reopen the container.

- **`read_only: true`** — read-only root filesystem. Requires writable `tmpfs`
  for `/tmp` and (usually) `/home/devuser`; the template includes the matching
  tmpfs lines, also commented. Useful for shared/registry-pushed images. Breaks
  anything that writes outside bind-mounted workspaces.
- **`cap_drop: ["ALL"]` + minimal `cap_add`** — the template keeps
  `CHOWN, SETUID, SETGID, DAC_OVERRIDE`. Tighter than Docker's default. May
  break `apt-get install` or in-container package compilation.
- **`security_opt: ["no-new-privileges:true"]`** — prevents `setuid` binaries
  from escalating. Safe to keep on if you don't need `sudo` inside the
  container.

## Secrets pattern

Avoid baking API keys into `.env`. A commented Docker-secrets example sits at the
bottom of
`templates/devcontainer/.devcontainer/docker-compose.yml.template`:

1. Store each secret in a file outside the repo, e.g.
   `~/.scbio-secrets/anthropic_api_key`, `chmod 600`.
2. Declare it as a top-level `secrets:` entry pointing at that file.
3. Attach it to `dev-core`/`dev-archr` with `target` + `mode: 0000`; the app
   reads it from `/run/secrets/<name>`.

```yaml
secrets:
  anthropic_api_key:
    file: ${HOME}/.scbio-secrets/anthropic_api_key

# inside dev-core / dev-archr:
    secrets:
      - source: anthropic_api_key
        target: anthropic_api_key
        mode: 0000
```

## When to enable opt-in hardening

- **Pushing the image to a shared registry** others will pull → enable
  `read_only: true`, `cap_drop`, `no-new-privileges`.
- **Running an autonomous agent with elevated tool access inside the container**
  → all of the above plus the secrets pattern.
- **Personal interactive use on your own machine** → defaults are sufficient.

## Explicitly not adopted

| Pattern | Why skipped |
|---------|-------------|
| `LD_PRELOAD` filesystem firewall | Too aggressive for interactive use — breaks `htop`, `strace`, debuggers, and any tool that reads config. |
| SetUID binary purge (`chmod a-s` on the whole FS) | Breaks `sudo`-based workflows for occasional package installs. |
| Command-level CLI wrappers | Loss of flexibility; mismatched threat model (you are the user, not an autonomous agent). |
| Application-level firewalls (`NODE_OPTIONS` monkeypatch) | R/Python image; would only cover Node tooling and break legitimate fs reads. |
