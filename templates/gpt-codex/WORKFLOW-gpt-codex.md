# Codex Workflow for Bioinformatics Projects

> “Codex writes fast, but only when the context is crisp. Spend the time to aim.”

---

## Core Principles

1. **Context before code** – load `plan.md`, `tasks.md`, and `HANDOFF.md` into the session before asking Codex to edit anything.
2. **Phase-based execution** – run one stage at a time; reset chat/context between stages to avoid drift.
3. **Explicit validation** – every code edit needs a verification step (command/run/log). If you skip a test, state why and how to reproduce it.
4. **Deterministic shell usage** – always specify `workdir`, prefer `rg` for search, and narrate any destructive action before executing.

---

## Documentation System

| File | Purpose | When to touch |
|------|---------|---------------|
| `plan.md` | Scientific narrative, known issues, hypotheses | After major discoveries or strategy changes |
| `tasks.md` | Numbered stages/substeps with ⏸️/🔄/✅/⛔ markers | In real time during execution |
| `notes.md` | Research links, troubleshooting, benchmarks | Whenever you look things up |
| `handoff.md` | Session summary (progress, blockers, next steps) | Every time you pause/finish |

Codex should never rely on conversation history; always re-read these files when returning.

---

## Recommended Session Flow

1. **Kickoff**
   - Read `GPT-CODEX.md`, `plan.md`, `tasks.md`, `HANDOFF.md`.
   - Restate the objective and propose a short plan (≤4 steps). Ask clarifying questions if scope is unclear.
2. **Plan Tool**
   - Use the Codex plan tool for any non-trivial work. Update the plan after each completed step so user logs stay fresh.
3. **Execution**
   - Use `apply_patch` for single-file edits; explain if you must fall back to other methods.
   - Run commands via `["bash","-lc", ...]` with `workdir` set. Mention important command output rather than dumping entire logs.
4. **Testing**
   - Run relevant scripts/tests (e.g., `./init-project.sh --help`, `bash -n script.sh`). If sandbox prevents it, document the exact command the user should run.
5. **Handoff**
   - Update `handoff.md` with accomplishments, blockers, next steps, file paths, and commands.
   - Mention pending approvals or tests that still need to run.

---

## Agent Stubs (Future Automation)

The `.gpt-codex/agents/` directory describes two future agents:

1. **handoff-writer** – scans `tasks.md`, outputs, and chat logs to auto-generate `handoff.md`.
2. **stage-reviewer** – validates a stage before it’s marked ✅ (checks files, logs, tasks status).

Until implemented, Codex must perform these duties manually and log work in `handoff.md`.

---

## Best Practices

- Keep diffs small; merge frequently.
- Prefer declarative descriptions (“Add GPT templates to init-project copy block”) before editing.
- If you detect unexpected repo changes, stop and ask how to proceed.
- For long-running work, use tmux (`tmux new -s stage1`) and log outputs to `logs/`.
- When editing docs, maintain ASCII text and short, focused sections.

---

## When Things Go Wrong

- **Sandbox denial** – rerun command with `with_escalated_permissions=true` and a one-line justification, then summarize results in the chat.
- **Merge conflicts** – describe the conflict, propose a resolution strategy, and wait for confirmation before editing.
- **Missing context** – request the relevant files or ask the user to clarify, rather than guessing.

---

By following this workflow, Codex instances stay aligned with the rest of the team, avoid context drift, and keep the single-cell environment reproducible.
