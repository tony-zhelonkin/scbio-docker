#!/usr/bin/env bash
#
# install_ai_tooling.sh - One-time setup for AI assistants inside devcontainers.
# This script is invoked via VS Code's postCreateCommand so that we can:
#   - verify MCP prerequisites (node/npm/npx, uvx)
#   - bootstrap Claude CLI when Claude assets are present
#   - leave room for future Codex/Claude hybrid automation

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

log() {
    printf '[AI TOOLING] %s\n' "$1"
}

require_cmd() {
    local cmd="$1"
    if ! command -v "$cmd" >/dev/null 2>&1; then
        log "ERROR: '$cmd' is not available. Make sure the base image includes it."
        return 1
    fi
    return 0
}

log "Running AI tooling setup for project: ${PROJECT_ROOT}"

MISSING=0
for tool in node npm npx uvx; do
    if require_cmd "$tool"; then
        log "$tool version: $($tool --version 2>&1 | head -n 1)"
    else
        MISSING=1
    fi
done

if [ "$MISSING" -ne 0 ]; then
    log "Missing required tools. Please rebuild the base image and rerun this script."
    exit 1
fi

############################################################
# Claude CLI bootstrap (only if project ships Claude assets)
############################################################
if [ -d "${PROJECT_ROOT}/.claude" ] || [ -f "${PROJECT_ROOT}/CLAUDE.md" ]; then
    log "Claude integration detected."
    if command -v claude >/dev/null 2>&1; then
        log "Claude CLI already installed: $(claude --version 2>/dev/null || echo 'version check failed')"
    elif [ "${SKIP_CLAUDE_CLI_INSTALL:-0}" = "1" ]; then
        log "SKIP_CLAUDE_CLI_INSTALL=1, skipping Claude CLI installation."
    else
        log "Installing Claude CLI..."
        if curl -fsSL https://claude.ai/install.sh | bash -s -- latest; then
            log "Claude CLI installed successfully."
        else
            log "ERROR: Failed to install Claude CLI. Re-run this script after checking network access."
            exit 1
        fi
    fi
else
    log "No Claude assets found; skipping Claude CLI installation."
fi

############################################################
# Codex hook (placeholder for future automation requirements)
############################################################
if [ -d "${PROJECT_ROOT}/.gpt-codex" ] || [ -f "${PROJECT_ROOT}/GPT-CODEX.md" ]; then
    log "GPT-Codex assets detected."
    if [ -f "${PROJECT_ROOT}/.mcp.json" ]; then
        log "MCP config found at ${PROJECT_ROOT}/.mcp.json"
        if grep -q '"context7"' "${PROJECT_ROOT}/.mcp.json" 2>/dev/null; then
            if [ -z "${CONTEXT7_API_KEY:-}" ]; then
                log "NOTE: context7 is enabled but CONTEXT7_API_KEY is empty. Provide an API key or run OAuth via Claude Code when available."
            fi
        fi
    else
        log "WARNING: .mcp.json missing. Copy templates/ai-common/mcp.json.template and rerun this script."
    fi
else
    log "No GPT-Codex assets found; skipping Codex-specific setup."
fi

log "AI tooling setup complete."
