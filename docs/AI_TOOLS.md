# AI Integration via SciAgent-toolkit

This repository supports optional AI assistance (Claude Code or Codex CLI) through the external SciAgent-toolkit. The previous `dev-claude-integration` and `dev-gpt-codex-integration` branches are deprecated.

## Quick Start (concise)
- Build AI-enabled image (default installs Sequential Thinking MCP):
  - `./build-ai-enabled.sh`
  - Optional MCPs: `--with-tooluniverse`, `--with-serena`
- Scaffold a project with AI setup:
  - `./init-project.sh ~/projects/my-ai-analysis ai-enabled`
- Open the project in VS Code and reopen in container.
- Inside the container, install your assistant and verify MCPs:
  - Claude Code: `curl -fsSL https://claude.ai/install.sh | bash -s latest && claude login && claude mcp list`
  - Codex CLI: `npm install -g @openai/codex && codex login && codex mcp list`

## Notes and scope
- Single source of truth for MCP configuration, available tools, and advanced usage is the SciAgent-toolkit documentation:
  - https://github.com/tony-zhelonkin/SciAgent-toolkit
- This repo only provides the container baseline and a build wrapper for an AI-enabled image. It intentionally avoids duplicating MCP/tool specifics.
- If you plan to build the AI image, ensure the `.sciagent/` submodule is initialized: `git submodule update --init --recursive`.

## Minimal configuration references
- Sequential Thinking MCP is included by default (invoked via `npx`).
- Optional MCPs:
  - ToolUniverse (broad scientific tooling) → `--with-tooluniverse`
  - Serena (code intelligence) → `--with-serena`
- API keys (if needed for specific tools) should be provided via environment or a project `.env` referenced from compose. Do not commit secrets.

## Migration
See docs/MIGRATION_AI_INTEGRATION.md for a short deprecation notice and the recommended path from the old AI branches to the AI-enabled build.
