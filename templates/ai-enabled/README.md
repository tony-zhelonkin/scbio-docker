# AI-Enabled Project (Template)

This template uses the optional AI-enabled image for scbio-docker and integrates with SciAgent-toolkit.

Quick use:
- Open the project in VS Code → "Dev Containers: Reopen in Container".
- Inside the container, install your assistant of choice and verify MCPs:
  - `claude login && claude mcp list` or `codex login && codex mcp list`.

Docs (single sources of truth):
- AI integration (concise): ../../docs/AI_TOOLS.md
- MCP/agent configuration: https://github.com/tony-zhelonkin/SciAgent-toolkit

Notes:
- AI branches `dev-claude-integration` and `dev-gpt-codex-integration` are deprecated.
- Keep API keys in environment or a local `.env` (never commit secrets).
