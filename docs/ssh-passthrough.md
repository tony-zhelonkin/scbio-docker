# SSH agent passthrough

Git over SSH inside the container uses the **host's SSH agent** — no private keys are copied into the container. This is implemented in the devcontainer/compose templates and wired automatically by `scripts/init-container.sh`.

## How it works

| Piece | Behaviour |
|---|---|
| Compose services | Both `dev-core` and `dev-archr` set `SSH_AUTH_SOCK=/ssh-agent` in `environment:`. |
| `init-container.sh` | Auto-detects the host `SSH_AUTH_SOCK`. If it points at a live socket, the script emits a volume line mounting it **read-only** at `/ssh-agent`. |
| No agent on host | The `{{SSH_AGENT_MOUNT}}` token renders empty — no mount is added, the compose file stays valid, and nothing else breaks. |

The relevant detection logic in `scripts/init-container.sh` (`build_ssh_agent_mount`) mounts the socket only when it is a real socket:

```bash
local sock="${SSH_AUTH_SOCK:-}"
if [[ -n "$sock" && -S "$sock" ]]; then
    printf '      - %s:/ssh-agent:ro\n' "$sock"
fi
```

The resulting service block looks like:

```yaml
environment:
  - SSH_AUTH_SOCK=/ssh-agent
volumes:
  - /run/user/1000/keyring/ssh:/ssh-agent:ro   # host agent socket (path varies)
```

Because the container's `SSH_AUTH_SOCK` always points at `/ssh-agent`, `ssh` (and therefore `git` over an SSH remote) inside the container talks to whatever agent the host has loaded. Keys never leave the host — only the agent socket is shared, read-only.

## Prerequisite: an agent on the host

Passthrough forwards an agent; it does not create one. On the host, before rendering the container (or before opening it in VS Code):

```bash
ssh-add -l          # should list your key(s)
```

If it reports "Could not open a connection to your authentication agent" or "no identities", start one and load your key:

```bash
eval "$(ssh-agent)"
ssh-add ~/.ssh/id_ed25519
```

`init-container.sh` reads `SSH_AUTH_SOCK` at render time, so make sure the agent is running in the shell you run it from. VS Code Dev Containers also forward the agent automatically when it launches the container, so the compose-level mount and VS Code's forwarding both point the container at the same host agent.

## Verify inside the container

```bash
echo "$SSH_AUTH_SOCK"          # -> /ssh-agent
ssh-add -l                     # lists the same keys as on the host
ssh -T git@github.com          # -> "Hi <user>! You've successfully authenticated..."
git clone git@github.com:org/repo.git
```

If `ssh-add -l` lists your keys, passthrough is working. If it prints "Could not open a connection to your authentication agent", the mount was skipped (see below).

## When no agent is running

If `SSH_AUTH_SOCK` is unset or does not point at a live socket when `init-container.sh` runs:

- The `{{SSH_AGENT_MOUNT}}` placeholder is replaced with nothing — the rendered `docker-compose.yml` has **no** `/ssh-agent` mount.
- Inside the container `SSH_AUTH_SOCK` still equals `/ssh-agent`, but that path does not exist, so `ssh-add -l` / `ssh` cannot reach an agent.
- Nothing else is affected; the container starts normally and you can use HTTPS git remotes or a token instead.

To fix: start an agent and load your key on the host, then **re-render** the container (`init-container.sh <project-dir>`) so the mount line is emitted, and restart the container. There is no `/dev/null` fallback socket — the mount is simply present or absent.

## Notes

- The socket is mounted `:ro`. That prevents the container from replacing the socket file; it does not restrict use of the loaded keys — any process in the container can sign with them while the container runs. This matches the trust model of a single-user interactive dev container.
- Only the agent socket is shared. `~/.ssh` (private keys, config) is **not** mounted, by design.
- Known-hosts prompts on first GitHub connection are harmless; accept the fingerprint or pre-populate `~/.ssh/known_hosts` on the host.

## Cross-references

- Renderer: [`scripts/init-container.sh`](../scripts/init-container.sh) (`build_ssh_agent_mount`)
- Compose template: `templates/devcontainer/.devcontainer/docker-compose.yml.template`
- Devcontainer wiring: [devcontainer.md](devcontainer.md)
- Container isolation defaults (tmpfs `/tmp`, `pids`, optional hardening): [isolation.md](isolation.md)
