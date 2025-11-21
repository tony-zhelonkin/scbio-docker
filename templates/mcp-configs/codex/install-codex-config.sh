#!/usr/bin/env bash
#
# Install Codex CLI MCP configuration for container
#
set -euo pipefail

CONFIG_DIR="${HOME}/.codex"
TEMPLATE_FILE="/opt/templates/mcp-configs/codex/config.toml.template"
CONFIG_FILE="${CONFIG_DIR}/config.toml"

# Create config directory if it doesn't exist
mkdir -p "${CONFIG_DIR}"

# Check if config already exists
if [ -f "${CONFIG_FILE}" ]; then
    echo "Codex MCP config already exists at ${CONFIG_FILE}"
    echo "Skipping installation to preserve existing configuration."
    echo "To reinstall, remove the file and run this script again."
    exit 0
fi

# Copy template
if [ -f "${TEMPLATE_FILE}" ]; then
    cp "${TEMPLATE_FILE}" "${CONFIG_FILE}"
    echo "Codex CLI MCP configuration installed: ${CONFIG_FILE}"
else
    echo "Warning: Template file not found at ${TEMPLATE_FILE}"
    echo "Creating minimal configuration..."

    cat > "${CONFIG_FILE}" <<'EOF'
# Codex CLI MCP Server Configuration
# Generated for scbio-docker AI-enabled environment

[mcp_servers.sequential-thinking]
command = "npx"
args = [
  "-y",
  "@modelcontextprotocol/server-sequential-thinking"
]
startup_timeout_sec = 30
EOF
    echo "Minimal Codex MCP config created: ${CONFIG_FILE}"
fi

# Validate TOML (if tomli is available)
if python3 -c "import tomli" 2>/dev/null; then
    if python3 -c "import tomli; tomli.load(open('${CONFIG_FILE}', 'rb'))" 2>/dev/null; then
        echo "✓ Configuration is valid TOML"
    else
        echo "✗ Warning: Configuration may have syntax errors"
        exit 1
    fi
else
    echo "Note: Install tomli to validate TOML syntax (pip install tomli)"
fi

echo ""
echo "To use Codex CLI with MCP servers:"
echo "  1. Install Codex CLI (if not already installed):"
echo "     npm install -g @openai/codex"
echo "  2. Authenticate:"
echo "     codex login"
echo "  3. Verify MCP servers:"
echo "     codex mcp list"
echo ""
