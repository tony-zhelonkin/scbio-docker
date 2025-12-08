#!/usr/bin/env bash
#
# Setup MCP configurations on container start
# This script automatically configures MCP servers for AI tools (Claude Code, Codex CLI)
#
set -euo pipefail

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo -e "${YELLOW}Setting up MCP configurations...${NC}"
echo "For full AI/MCP docs see SciAgent‑toolkit: https://github.com/tony-zhelonkin/SciAgent-toolkit"

# Detect which AI tools are available
HAS_CLAUDE=false
HAS_CODEX=false

if command -v claude &>/dev/null; then
    HAS_CLAUDE=true
fi

if command -v codex &>/dev/null; then
    HAS_CODEX=true
fi

# Install Claude Code config if needed
if [ "$HAS_CLAUDE" = true ]; then
    if [ ! -f "${HOME}/.config/claude/mcp-config.json" ]; then
        echo "Installing Claude Code MCP configuration..."
        if [ -f /opt/templates/mcp-configs/claude/install-claude-config.sh ]; then
            bash /opt/templates/mcp-configs/claude/install-claude-config.sh
        else
            # Fallback: create minimal config if install script not found
            mkdir -p "${HOME}/.config/claude"
            cat > "${HOME}/.config/claude/mcp-config.json" <<'EOF'
{
  "mcpServers": {
    "sequential-thinking": {
      "command": "npx",
      "args": [
        "-y",
        "@modelcontextprotocol/server-sequential-thinking"
      ]
    }
  }
}
EOF
            echo "Created minimal Claude MCP config"
        fi
    else
        echo -e "${GREEN}Claude Code MCP config already exists${NC}"
    fi
fi

# Install Codex CLI config if needed
if [ "$HAS_CODEX" = true ]; then
    if [ ! -f "${HOME}/.codex/config.toml" ]; then
        echo "Installing Codex CLI MCP configuration..."
        if [ -f /opt/templates/mcp-configs/codex/install-codex-config.sh ]; then
            bash /opt/templates/mcp-configs/codex/install-codex-config.sh
        else
            # Fallback: create minimal config if install script not found
            mkdir -p "${HOME}/.codex"
            cat > "${HOME}/.codex/config.toml" <<'EOF'
# Codex CLI MCP Server Configuration

[mcp_servers.sequential-thinking]
command = "npx"
args = [
  "-y",
  "@modelcontextprotocol/server-sequential-thinking"
]
startup_timeout_sec = 30
EOF
            echo "Created minimal Codex MCP config"
        fi
    else
        echo -e "${GREEN}Codex CLI MCP config already exists${NC}"
    fi
fi

# If neither tool is installed, provide guidance
if [ "$HAS_CLAUDE" = false ] && [ "$HAS_CODEX" = false ]; then
    echo -e "${YELLOW}No AI tools detected.${NC}"
    echo "To install AI tools inside this container:"
    echo ""
    echo "  Claude Code:"
    echo "    curl -fsSL https://claude.ai/install.sh | bash -s latest"
    echo "    source ~/.bashrc"
    echo "    claude login"
    echo ""
    echo "  Codex CLI:"
    echo "    npm install -g @openai/codex"
    echo "    codex login"
    echo ""
    echo "After installation, MCP configs will be set up automatically on next container start."
fi

echo -e "${GREEN}MCP setup complete${NC}"
