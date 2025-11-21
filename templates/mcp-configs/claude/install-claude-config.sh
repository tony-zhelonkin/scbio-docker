#!/usr/bin/env bash
#
# Install Claude Code MCP configuration for container
#
set -euo pipefail

CONFIG_DIR="${HOME}/.config/claude"
TEMPLATE_FILE="/opt/templates/mcp-configs/claude/mcp-config.json.template"
CONFIG_FILE="${CONFIG_DIR}/mcp-config.json"

# Create config directory if it doesn't exist
mkdir -p "${CONFIG_DIR}"

# Check if config already exists
if [ -f "${CONFIG_FILE}" ]; then
    echo "Claude MCP config already exists at ${CONFIG_FILE}"
    echo "Skipping installation to preserve existing configuration."
    echo "To reinstall, remove the file and run this script again."
    exit 0
fi

# Copy template
if [ -f "${TEMPLATE_FILE}" ]; then
    cp "${TEMPLATE_FILE}" "${CONFIG_FILE}"
    echo "Claude Code MCP configuration installed: ${CONFIG_FILE}"
else
    echo "Warning: Template file not found at ${TEMPLATE_FILE}"
    echo "Creating minimal configuration..."

    cat > "${CONFIG_FILE}" <<'EOF'
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
    echo "Minimal Claude MCP config created: ${CONFIG_FILE}"
fi

# Validate JSON
if python3 -m json.tool "${CONFIG_FILE}" > /dev/null 2>&1; then
    echo "✓ Configuration is valid JSON"
else
    echo "✗ Warning: Configuration is not valid JSON"
    exit 1
fi

echo ""
echo "To use Claude Code with MCP servers:"
echo "  1. Install Claude Code (if not already installed):"
echo "     curl -fsSL https://claude.ai/install.sh | bash -s latest"
echo "  2. Authenticate:"
echo "     claude login"
echo "  3. Verify MCP servers:"
echo "     claude mcp list"
echo ""
