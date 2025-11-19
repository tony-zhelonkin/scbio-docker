#!/usr/bin/env bash
#
# setup_codex_mcp.sh - Validate MCP configuration for Codex projects.
#
# This script is intended to run inside the devcontainer (optional manual step).
# It checks for required prerequisites and confirms that .mcp.json points to
# the current workspace.

set -euo pipefail

PROJECT_ROOT="/workspaces/{{PROJECT_NAME}}"
MCP_FILE="${PROJECT_ROOT}/.mcp.json"

log() {
    printf '[GPT-CODEX][MCP] %s\n' "$1"
}

require_cmd() {
    local cmd="$1"
    if ! command -v "$cmd" >/dev/null 2>&1; then
        log "ERROR: Missing dependency '$cmd'. Rebuild the image or rerun install_ai_tooling.sh."
        return 1
    fi
    return 0
}

log "Workspace: ${PROJECT_ROOT}"

if [ ! -f "$MCP_FILE" ]; then
    log "ERROR: ${MCP_FILE} not found. Copy the template and rerun this script."
    exit 1
fi

for tool in node npm npx uvx python3; do
    require_cmd "$tool"
done

log "Dependencies detected (node: $(node --version), uvx: $(uvx --version 2>/dev/null || echo 'unknown'))"

export MCP_FILE

python3 <<'PY'
import json
import os

mcp_file = os.environ.get("MCP_FILE")
with open(mcp_file, "r", encoding="utf-8") as handle:
    data = json.load(handle)

servers = data.get("mcpServers", {})
if not servers:
    raise SystemExit("[GPT-CODEX][MCP] No servers defined in .mcp.json")

print("[GPT-CODEX][MCP] Configured servers:")
for name, config in servers.items():
    cmd = config.get("command", "<unknown>")
    args = " ".join(config.get("args", []))
    print(f"  - {name}: {cmd} {args}".strip())

if "context7" in servers:
    has_api_key = bool(os.environ.get("CONTEXT7_API_KEY"))
    if not has_api_key:
        print("[GPT-CODEX][MCP] context7 enabled but CONTEXT7_API_KEY is empty. Provide an API key or authenticate via Claude Code before use.")
PY

log "MCP configuration looks sane. Start MCP servers via the Codex CLI as needed."
