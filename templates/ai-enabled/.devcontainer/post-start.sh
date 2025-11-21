#!/bin/bash
#
# Post-start checks and MCP setup for AI-enabled environment
#

GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

echo "==================================================================="
echo "           AI-Enabled Container Environment Check"
echo "==================================================================="
echo ""
echo "R version:        $(R --version | head -n1)"
echo "Python version:   $(python3 --version)"
echo "Node.js version:  $(node --version)"
echo "npm version:      $(npm --version)"
echo "Working directory: $(pwd)"
echo "User:             $(whoami) (UID=$(id -u), GID=$(id -g))"
echo ""

# Check if MCP setup script exists
if [ -f /opt/scripts/setup-mcp.sh ]; then
    echo -e "${YELLOW}Setting up MCP configurations...${NC}"
    bash /opt/scripts/setup-mcp.sh
else
    echo -e "${YELLOW}Note: MCP setup script not found. MCP servers will be configured manually.${NC}"
    echo "To configure MCP servers:"
    echo "  - Claude Code: See ~/.config/claude/mcp-config.json"
    echo "  - Codex CLI:   See ~/.codex/config.toml"
fi

echo ""
echo "==================================================================="
echo "                     Environment Ready"
echo "==================================================================="
echo ""
echo "Quick commands:"
echo "  R / radian          - Start R session"
echo "  python3             - Start Python session"
echo "  usepy <env>         - Switch Python environment (base/squid/atac/comms)"
echo ""
echo "AI Tools (install if needed):"
echo "  claude --version    - Check if Claude Code is installed"
echo "  codex --version     - Check if Codex CLI is installed"
echo ""
echo "For help, see:"
echo "  README.md           - Project documentation"
echo "  ../../docs/AI_TOOLS.md - AI tools usage guide"
echo ""
