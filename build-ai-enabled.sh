#!/usr/bin/env bash
#
# Build AI-enabled variant of scbio-docker
# This script extends the base image with minimal MCP prerequisites.
# Full agent/MCP configuration lives in the external SciAgent-toolkit.
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

log_info()  { echo -e "${YELLOW}[INFO]${NC} $*"; }
log_ok()    { echo -e "${GREEN}[OK]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR]${NC} $*"; }

# Default values
BASE_IMAGE="scdock-r-dev:v0.5.1"
TAG="scdock-r-dev-ai:v0.5.1"
INSTALL_TOOLUNIVERSE=0
INSTALL_SERENA=0

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --base-image)
            BASE_IMAGE="$2"
            shift 2
            ;;
        --tag)
            TAG="$2"
            shift 2
            ;;
        --with-tooluniverse)
            INSTALL_TOOLUNIVERSE=1
            shift
            ;;
        --with-serena)
            INSTALL_SERENA=1
            shift
            ;;
        --help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Build AI-enabled variant of scbio-docker (MCP-ready)"
            echo ""
            echo "Options:"
            echo "  --base-image IMAGE        Base image to extend (default: $BASE_IMAGE)"
            echo "  --tag TAG                 Tag for the AI-enabled image (default: $TAG)"
            echo "  --with-tooluniverse       Install ToolUniverse MCP server (600+ scientific tools)"
            echo "  --with-serena             Install Serena MCP server (code intelligence)"
            echo "  --help                    Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0                                    # Build with Sequential Thinking only (default)"
            echo "  $0 --with-tooluniverse                # Include ToolUniverse"
            echo "  $0 --with-tooluniverse --with-serena  # Include all MCP servers"
            echo ""
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            echo "Run '$0 --help' for usage"
            exit 1
            ;;
    esac
done

echo "==================================================================="
echo "           Building AI-Enabled scbio-docker Image"
echo "==================================================================="
echo ""
echo "Base image:           $BASE_IMAGE"
echo "Output tag:           $TAG"
echo "Sequential Thinking:  YES (always installed)"
echo "ToolUniverse:         $([ $INSTALL_TOOLUNIVERSE -eq 1 ] && echo 'YES' || echo 'NO')"
echo "Serena:               $([ $INSTALL_SERENA -eq 1 ] && echo 'YES' || echo 'NO')"
echo ""

# Ensure submodule is initialized
if [ ! -f "${SCRIPT_DIR}/.sciagent/scripts/setup_mcp_infrastructure.sh" ]; then
    log_error "SciAgent-toolkit submodule not found!"
    echo ""
    echo "Please initialize the submodule first:"
    echo "  cd \"$SCRIPT_DIR\""
    echo "  git submodule update --init --recursive"
    echo ""
    echo "See .submodule-setup.md for detailed instructions."
    exit 1
fi

log_ok "SciAgent-toolkit submodule found"

# Check if base image exists
if ! docker image inspect "$BASE_IMAGE" &>/dev/null; then
    log_error "Base image '$BASE_IMAGE' not found!"
    echo ""
    echo "Please build the base image first:"
    echo "  cd \"$SCRIPT_DIR\""
    echo "  ./build-optimized.sh"
    echo ""
    exit 1
fi

log_ok "Base image '$BASE_IMAGE' exists"

# Build AI-enabled variant
log_info "Building AI-enabled image..."
echo ""

if docker build "${SCRIPT_DIR}" \
    -f "${SCRIPT_DIR}/.devcontainer/Dockerfile.ai-enabled" \
    --build-arg BASE_IMAGE="$BASE_IMAGE" \
    --build-arg INSTALL_TOOLUNIVERSE="$INSTALL_TOOLUNIVERSE" \
    --build-arg INSTALL_SERENA="$INSTALL_SERENA" \
    -t "$TAG"; then

    echo ""
    log_ok "Build complete: $TAG"
    echo ""
    echo "==================================================================="
    echo "                         Next Steps"
    echo "==================================================================="
    echo ""
    echo "1. Verify the image:"
    echo "   docker images | grep scdock-r-dev-ai"
    echo ""
    echo "2. Create a new AI-enabled project:"
    echo "   ./init-project.sh ~/projects/my-ai-analysis ai-enabled"
    echo ""
    echo "3. Open in VS Code and reopen in container:"
    echo "   code ~/projects/my-ai-analysis"
    echo "   # Then: Ctrl+Shift+P → 'Dev Containers: Reopen in Container'"
    echo ""
    echo "4. Inside the container, install your preferred assistant (Claude/Codex) and verify MCPs."
    echo "   See docs/AI_TOOLS.md and SciAgent-toolkit for commands and configuration."
    echo ""
    echo "Docs (single sources of truth):"
    echo "  - docs/AI_TOOLS.md (concise usage)"
    echo "  - docs/MIGRATION_AI_INTEGRATION.md (AI branches deprecated)"
    echo "  - SciAgent‑toolkit: https://github.com/tony-zhelonkin/SciAgent-toolkit"
    echo ""
else
    echo ""
    log_error "Build failed"
    exit 1
fi
