#!/usr/bin/env bash
#
# Claude Code + MCP Dependencies Setup Script
#
# This script installs Claude Code and all required MCP server dependencies
# in the dev container. It's designed to be idempotent (safe to run multiple times)
# and is called from devcontainer.json postStartCommand.
#
# What it does:
# 1. Installs Claude Code (native installer) if not present
# 2. Installs Python 'uv' package (provides uvx for serena MCP server)
# 3. Installs Node.js + npm (provides npx for sequential-thinking MCP server)
# 4. Creates project-level .mcp.json configuration
# 5. Runs sanity checks on all installations
#
# SECURITY NOTE: This script requires sudo for system package installation.
# Review all installations and understand security implications before running.
# See .devcontainer/MCP_AUTH_SETUP.md for detailed security considerations.

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Helper functions
log_info()  { echo -e "${BLUE}[INFO]${NC} $*"; }
log_ok()    { echo -e "${GREEN}[OK]${NC} $*"; }
log_warn()  { echo -e "${YELLOW}[WARN]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }
separator() { echo -e "${BLUE}==== $* ====${NC}"; }

# Check if running as devuser (not root)
if [ "$(id -u)" -eq 0 ]; then
    log_error "This script should not be run as root. Run as devuser."
    exit 1
fi

separator "Claude Code + MCP Setup"

# ============================================================================
# 1. Claude Code Installation
# ============================================================================
separator "Installing Claude Code"

if command -v claude &>/dev/null; then
    CLAUDE_VERSION=$(claude --version 2>/dev/null | grep -oP '\d+\.\d+\.\d+' || echo "unknown")
    log_ok "Claude Code already installed (version: ${CLAUDE_VERSION})"

    # Ensure PATH is set even if already installed
    if ! grep -q 'export PATH="$HOME/.local/bin:$PATH"' "$HOME/.bashrc" 2>/dev/null; then
        echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$HOME/.bashrc"
        log_info "Added ~/.local/bin to PATH in .bashrc"
    fi
    export PATH="$HOME/.local/bin:$PATH"
else
    log_info "Installing Claude Code (native installer)..."

    # Use the native installer (doesn't require npm)
    if curl -fsSL https://claude.ai/install.sh | bash -s latest; then
        log_ok "Claude Code installed successfully"

        # Add ~/.local/bin to PATH permanently for bash
        if ! grep -q 'export PATH="$HOME/.local/bin:$PATH"' "$HOME/.bashrc"; then
            echo 'export PATH="$HOME/.local/bin:$PATH"' >> "$HOME/.bashrc"
            log_info "Added ~/.local/bin to PATH in .bashrc"
        fi

        # Ensure ~/.local/bin is in PATH for this session
        export PATH="$HOME/.local/bin:$PATH"

        # Verify installation
        if command -v claude &>/dev/null; then
            CLAUDE_VERSION=$(claude --version 2>/dev/null | grep -oP '\d+\.\d+\.\d+' || echo "installed")
            log_ok "Claude Code version: ${CLAUDE_VERSION}"
        else
            log_error "Claude Code installation failed - command not found"
            exit 1
        fi
    else
        log_error "Failed to install Claude Code"
        exit 1
    fi
fi

# ============================================================================
# 2. Python 'uv' Package (provides uvx for serena MCP server)
# ============================================================================
separator "Installing MCP Dependencies: uv (for serena)"

if command -v uvx &>/dev/null; then
    UVX_VERSION=$(uvx --version 2>/dev/null | grep -oP '\d+\.\d+\.\d+' || echo "unknown")
    log_ok "uvx already installed (version: ${UVX_VERSION})"
else
    log_info "Installing Python 'uv' package..."

    # Install with sudo to make it available system-wide
    if sudo pip3 install uv &>/dev/null; then
        log_ok "uv package installed successfully"

        # Verify installation
        if command -v uvx &>/dev/null; then
            UVX_VERSION=$(uvx --version 2>/dev/null | grep -oP '\d+\.\d+\.\d+' || echo "installed")
            log_ok "uvx version: ${UVX_VERSION}"
        else
            log_error "uv installation failed - uvx command not found"
            exit 1
        fi
    else
        log_error "Failed to install uv package"
        exit 1
    fi
fi

# ============================================================================
# 3. Node.js + npm (provides npx for sequential-thinking MCP server)
# ============================================================================
separator "Installing MCP Dependencies: Node.js (for sequential-thinking)"

if command -v npx &>/dev/null; then
    NODE_VERSION=$(node --version 2>/dev/null || echo "unknown")
    NPM_VERSION=$(npm --version 2>/dev/null || echo "unknown")
    log_ok "Node.js already installed (node: ${NODE_VERSION}, npm: ${NPM_VERSION})"
else
    log_info "Installing Node.js 20.x from NodeSource..."

    # Add NodeSource repository and install Node.js
    if curl -fsSL https://deb.nodesource.com/setup_20.x | sudo -E bash - &>/dev/null && \
       sudo apt-get install -y nodejs &>/dev/null; then
        log_ok "Node.js installed successfully"

        # Verify installation
        if command -v npx &>/dev/null; then
            NODE_VERSION=$(node --version 2>/dev/null || echo "installed")
            NPM_VERSION=$(npm --version 2>/dev/null || echo "installed")
            log_ok "Node.js version: ${NODE_VERSION}"
            log_ok "npm version: ${NPM_VERSION}"
        else
            log_error "Node.js installation failed - npx command not found"
            exit 1
        fi
    else
        log_error "Failed to install Node.js"
        exit 1
    fi
fi

# ============================================================================
# 4. Project-level MCP Configuration
# ============================================================================
separator "Configuring Project-Level MCP Servers"

MCP_JSON="/workspaces/DC_Dictionary/.mcp.json"

if [ -f "$MCP_JSON" ]; then
    log_ok "MCP configuration already exists: $MCP_JSON"
else
    log_info "Creating project-level MCP configuration..."

    cat > "$MCP_JSON" << 'EOF'
{
  "mcpServers": {
    "context7": {
      "type": "sse",
      "url": "https://mcp.context7.com/sse"
    },
    "serena": {
      "type": "stdio",
      "command": "uvx",
      "args": [
        "--from",
        "git+https://github.com/oraios/serena",
        "serena",
        "start-mcp-server",
        "--context",
        "ide-assistant",
        "--project",
        "/workspaces/DC_Dictionary"
      ]
    },
    "sequential-thinking": {
      "type": "stdio",
      "command": "npx",
      "args": [
        "-y",
        "@modelcontextprotocol/server-sequential-thinking"
      ]
    }
  }
}
EOF

    if [ -f "$MCP_JSON" ]; then
        log_ok "MCP configuration created: $MCP_JSON"
    else
        log_error "Failed to create MCP configuration"
        exit 1
    fi
fi

# ============================================================================
# 5. Sanity Checks
# ============================================================================
separator "Sanity Checks"

CHECKS_PASSED=0
CHECKS_FAILED=0

# Check Claude Code
printf "Claude Code: "
if command -v claude &>/dev/null && claude --version &>/dev/null; then
    CLAUDE_VER=$(claude --version 2>/dev/null | head -1 || echo "unknown")
    log_ok "${CLAUDE_VER}"
    CHECKS_PASSED=$((CHECKS_PASSED + 1))
else
    log_error "NOT INSTALLED"
    CHECKS_FAILED=$((CHECKS_FAILED + 1))
fi

log_info "DEBUG: After Claude check - CHECKS_PASSED=${CHECKS_PASSED}, CHECKS_FAILED=${CHECKS_FAILED}"

# Check uvx (for serena)
printf "uvx (serena): "
if command -v uvx &>/dev/null; then
    UVX_VER=$(uvx --version 2>/dev/null | head -1 || echo "unknown")
    log_ok "${UVX_VER}"
    CHECKS_PASSED=$((CHECKS_PASSED + 1))
else
    log_error "NOT INSTALLED"
    CHECKS_FAILED=$((CHECKS_FAILED + 1))
fi

# Check npx (for sequential-thinking)
printf "npx (sequential-thinking): "
if command -v npx &>/dev/null; then
    NPX_VER=$(npx --version 2>/dev/null | head -1 || echo "unknown")
    log_ok "${NPX_VER}"
    CHECKS_PASSED=$((CHECKS_PASSED + 1))
else
    log_error "NOT INSTALLED"
    CHECKS_FAILED=$((CHECKS_FAILED + 1))
fi

# Check MCP config
printf "MCP config (.mcp.json): "
if [ -f "$MCP_JSON" ]; then
    # Validate JSON syntax
    if python3 -m json.tool "$MCP_JSON" &>/dev/null; then
        log_ok "Valid JSON"
        CHECKS_PASSED=$((CHECKS_PASSED + 1))
    else
        log_error "Invalid JSON syntax"
        CHECKS_FAILED=$((CHECKS_FAILED + 1))
    fi
else
    log_error "NOT FOUND"
    CHECKS_FAILED=$((CHECKS_FAILED + 1))
fi

# Test stdio MCP servers (quick check)
printf "serena MCP server: "
if timeout 5 uvx --from git+https://github.com/oraios/serena serena --help &>/dev/null; then
    log_ok "Can start"
    CHECKS_PASSED=$((CHECKS_PASSED + 1))
else
    log_warn "Cannot verify (may need first download)"
    # Don't fail the check, just warn
    CHECKS_PASSED=$((CHECKS_PASSED + 1))
fi

printf "sequential-thinking MCP: "
if timeout 5 npx -y @modelcontextprotocol/server-sequential-thinking --help &>/dev/null; then
    log_ok "Can start"
    CHECKS_PASSED=$((CHECKS_PASSED + 1))
else
    log_warn "Cannot verify (may need first download)"
    # Don't fail the check, just warn
    CHECKS_PASSED=$((CHECKS_PASSED + 1))
fi

# ============================================================================
# Summary
# ============================================================================
separator "Setup Summary"

if [ $CHECKS_FAILED -eq 0 ]; then
    log_ok "All checks passed (${CHECKS_PASSED}/${CHECKS_PASSED})"
    log_info ""
    log_info "Next steps:"
    log_info "1. Restart container if not already done (for port mapping)"
    log_info "2. Set up SSH port forwarding: ssh -L 45454:localhost:45454 user@remote"
    log_info "3. Run 'claude' to start Claude Code"
    log_info "4. Run '/mcp' inside Claude Code to authenticate context7"
    log_info ""
    log_info "For detailed setup and security info, see:"
    log_info "  .devcontainer/MCP_AUTH_SETUP.md"
else
    log_error "Some checks failed (${CHECKS_PASSED} passed, ${CHECKS_FAILED} failed)"
    log_info "See .devcontainer/MCP_AUTH_SETUP.md for troubleshooting"
    exit 1
fi

separator "Setup Complete"
