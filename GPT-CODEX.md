# GPT-CODEX.md

This handbook explains how to work on `scbio-docker` with the **Codex CLI (GPT-based)** agent. It mirrors `CLAUDE.md`, but is tuned for Codex’s execution model (shell-focused, deterministic actions, explicit approval flow).

## Environment Recap

- **Images:** `scdock-r-dev:v0.5.2` (core) and `scdock-r-archr:v0.5.2` (ArchR wrapper)
- **Languages:** R 4.5 + Bioconductor 3.21, Python 3.10 virtualenvs (`/opt/venvs/{base,squid,atac,comms}`)
- **Dev UX:** VS Code Dev Containers, tmux-friendly shell, `init-project.sh` scaffolding
- **Key dirs:** `.devcontainer/`, `.vscode/`, `templates/`, `.claude/` (Claude assets), `.gpt-codex/` (Codex assets on dev-gpt-codex-integration)

## Codex-Specific Guidance

1. **Plan-first** – Run `/plan` (or use the plan tool) whenever a task takes more than a couple edits. Keep plans short (≤4 steps) and update after each step.
2. **Shell usage** – Commands must run via `["bash","-lc", "..."]` with `workdir` set. Prefer `rg` for search, avoid destructive operations (`git reset --hard`, `rm -rf`) unless explicitly requested.
3. **Approvals** – When the sandbox blocks required writes/networking, rerun the same command with `with_escalated_permissions=true` plus a one-line justification, instead of messaging the user first.
4. **Context hygiene** – Don’t rely on previous turns. If you need design history, read `plan.md`, `tasks.md`, `INIT_PROJECT_ENHANCEMENT.md`, and the branch vignette (`BRANCH_MANAGEMENT.md`).
5. **Testing** – Default expectation is to run linters/tests relevant to the touched files (e.g., `./init-project.sh --help`, dry-run builds) unless the sandbox forbids it. Document any skipped tests with instructions for the user.
6. **Handoff discipline** – Update `HANDOFF_LATEST.md` (or a project-level `handoff.md`) with: what changed, remaining blockers, next commands to run.

## Working With `init-project.sh`

- Script lives at repository root. Use `./init-project.sh DEST TEMPLATE [flags]` or its symlink `init-scproject`.
- `dev` branch scaffolds core files only. Branches with AI integrations (`dev-claude-integration`, `dev-gpt-codex-integration`) add their respective templates automatically.
- Key flags: `--interactive`, `--data-mount LABEL:/host/path[:ro]`, `--git-init`, `--submodules list`.
- Validate generated projects by checking `.devcontainer/docker-compose.yml`, `.env`, `.gitignore`, and AI docs (e.g., `GPT-CODEX.md`, `.gpt-codex/agents/`).
- Pass `--ai codex` (or `--ai both`) to copy Codex docs, `.gpt-codex/agents`, the MCP helper script, and a project-specific `.mcp.json`. The script substitutes `{{PROJECT_NAME}}` and `/workspaces/<project>` automatically, so there’s no post-init editing. The shared template also enables `context7` by default—set `CONTEXT7_API_KEY` (API key workflow) or plan to authenticate via Claude Code if you need that MCP server.

## Bootstrap + MCP Configuration

- The base image now bundles **Node.js 20.x + npm/npx, `nvm`, and the Python `uvx` runner**, so Codex MCP servers (`context7`, `serena`, `sequential-thinking`) are available without re-installing dependencies per container.
- Every generated project includes `.devcontainer/scripts/install_ai_tooling.sh`, executed via `postCreateCommand`. It:
  1. Verifies `node`, `npm`, `npx`, and `uvx`.
  2. Installs the Claude CLI when Claude assets exist (harmless for Codex-only projects).
  3. Warns if `.mcp.json` is missing so you can regenerate it from `templates/ai-common/mcp.json.template`, and flags when `context7` is enabled but `CONTEXT7_API_KEY` is unset.
- Use the Codex-specific helper whenever you need to revalidate MCP setup:

```bash
bash -lc '.devcontainer/scripts/setup_codex_mcp.sh'
```

- Typical project scaffold:

```bash
./init-project.sh ~/projects/cardio-atac basic-rna --ai codex --git-init
```
- To disable `context7`, remove that block from `.mcp.json` or leave `CONTEXT7_API_KEY` blank (Codex will skip it); Claude users can still authenticate via `/mcp` when needed.

## AI Template Philosophy

- **Token budget:** Target 600–700 tokens per AI context file (Codex context window is smaller than Claude’s when running multi-shot coding sessions).
- **Docs separation:** `GPT-CODEX.md` carries only evergreen context (warnings, key commands). The analysis story lives in `plan.md`, execution steps in `tasks.md`, research scraps in `notes.md`, and session continuity in `handoff.md`.
- **Agent hooks:** `.gpt-codex/agents/` defines future automation (handoff writer, stage reviewer). They are stubs today but document the desired behavior so we can implement agents later.

## Workflow Expectations For Codex Sessions

1. **Kickoff**
   - Read `plan.md`, `tasks.md`, `HANDOFF_LATEST.md`.
   - Summarize the goal and constraints, propose a plan, wait for confirmation if scope is ambiguous.
2. **Execution**
   - Take one subtask at a time. After finishing, summarize changes, run validations, update docs.
   - When editing files, prefer `apply_patch` for single-file modifications; mention why if you fall back to other methods.
3. **Testing & Validation**
   - For scripts: run the relevant helper (`bash -lc "./init-project.sh --help"`). For docs: lint via markdown preview or `vale` (if configured). If testing is impossible, state exactly why and how to reproduce.
4. **Handoff**
   - Update `HANDOFF_LATEST.md` with: what changed, current branch, outstanding steps, commands/files to revisit.
   - Mention whether `templates/gpt-codex/` still needs edits or if init scaffolding was tested.

## Quick Reference Commands

```bash
# Branch hygiene
git status -sb
git log --oneline --graph --decorate -10

# Search
rg -n "claude" -g '*.md'

# Test init script (dry run)
./init-project.sh /tmp/test-gpt basic-rna --interactive <<'EOF'
# answers...
EOF

# Build image (multi-stage)
./build-optimized.sh --github-pat $GITHUB_PAT
```

## When To Escalate

- Need to run Docker builds, install packages, or access network resources blocked by sandbox.
- Have to inspect system paths outside the workspace.
- Need destructive git operations because the user explicitly asked (e.g., reset branch).

Always annotate the escalation with intent: “Need to run docker build to verify init script changes,” etc.

## Support Files

- `BRANCH_MANAGEMENT.md` – explains how `dev`, `dev-claude-integration`, and `dev-gpt-codex-integration` relate.
- `CLAUDE.md` – reference for Claude-specific workflows (read to mirror best practices).
- `GPT-CODEX.md` (this file) – keep updated as Codex capabilities evolve.
- `HANDOFF_LATEST.md` – session-by-session status; update before logging off.
