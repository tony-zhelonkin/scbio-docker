# VS Code Remote stability — why containers freeze, and the fixes

Diagnosed 2026-07-09 on the JCR workstation. Symptom: Remote-SSH / Dev-Container
sessions (esp. Meta-Aging, STING-JR) freeze 1-3 min after connecting; terminals
accumulate stray `source /opt/venvs/base/bin/activate` lines after each reconnect.

## Root cause (not host resources)

Host was healthy (30/502 GiB RAM, load ~3.5/72 cores, no OOM/SIGKILL). The freeze
is a **reconnect feedback loop**, not exhaustion:

1. Network/VPN blip drops the SSH channel. sshd is *not* the killer
   (`ClientAliveCountMax 720` = ~6 h tolerance); it just reports the dead socket.
2. VS Code reconnects but does **not** reap the old server stack. Each reconnect
   spawns a fresh `server-main → extensionHost → Pylance → R LSP → fileWatcher`
   plus a host-side `vscode-remote-containers-server-*.js` shim.
3. Stacks pile up (seen: 4 Pylance + 4 R LSP + 4 watchers in one STING-JR
   container; 40 orphaned host shims). Each new stack re-indexes + re-registers
   recursive watchers over large RNA-seq trees → CPU/inotify burst = the freeze.
4. **inotify is a host-wide budget per UID** (`max_user_instances=512`). Every
   container runs as the same host UID (`788715489`), so all watchers share one
   pool — orphan sprawl crowds it and new watchers stall.

The stray `activate` lines: image `/etc/bash.bashrc` already activates the venv
silently. `python.terminal.activateEnvironment: true` in per-repo settings made
the Python extension **re-inject** the command into every terminal on each
reconnect — stacking, and (dangerously) submitting as a prompt into Claude/codex
tabs. Machine-level `false` was overridden by workspace-level `true`.

## Fixes applied

- **`python.terminal.activateEnvironment: false`** everywhere (redundant with
  bashrc; stops terminal injection). Template + 25 repos.
- **Git + watcher caps** in every repo `.vscode/settings.json`:
  `git.repositoryScanMaxDepth: 1`, `git.autoRepositoryDetection: openEditors`,
  `git.detectSubmodules: false`, and `files.watcherExclude` over
  `00_data / 01_modules / .venv / *-env / renv/library / __pycache__` +
  `03_results/checkpoints`. Stops per-subrepo Git models + watcher thrash.
- SSH split kept (see `ssh-passthrough.md`): VS Code uses the **no-mux** host
  entry (`ControlPath none`); terminal uses the multiplexed entry. Sharing a mux
  socket with Remote-SSH is the classic "raw ssh fine, VS Code flaky" cause.

## Tooling (dev-env/scripts, on PATH)

- `vscode-sync-settings [--apply]` — re-applies the two policies above across all
  repos. Idempotent, dry-run by default. Run after scaffolding a new project.
- `vscode-cleanup` — survey/reap stale server sprawl:
  - `vscode-cleanup` — list dev containers: status, last VS Code activity, #stacks
  - `vscode-cleanup stale [DAYS]` — containers idle > N days (stop candidates)
  - `vscode-cleanup orphans -y` — reap host-side dead Remote-Containers shims
  - `vscode-cleanup trim <substr> -y` — keep newest server stack, kill duplicates
  - `vscode-cleanup restart|stop <substr> -y` — cleanest reset / free resources

## Habits

- **Close Remote Connection** (don't just shut the laptop) to avoid orphan stacks.
- Run long-lived Claude/codex sessions inside **tmux in the container** so
  reconnect junk hits the outer shell, not the agent, and drops don't SIGKILL it.
- Add laptop `~/.ssh/config`: `ServerAliveInterval 15`, `ServerAliveCountMax 8`.
- Keep only 1-2 heavy containers hot; `vscode-cleanup stale` then `stop` the rest.
