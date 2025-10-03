#!/usr/bin/env bash
# init-project.sh - Initialize a new project directory from templates
#
# Usage:
#   ./init-project.sh <project-dir> <template-name>
#
# Templates:
#   basic-rna       - Standard RNA-seq analysis
#   multimodal      - RNA + ATAC or CITE-seq
#   archr-focused   - ArchR scATAC-seq analysis
#   example-DMATAC  - Differential chromatin accessibility

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEMPLATES_DIR="${SCRIPT_DIR}/templates"

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

usage() {
    echo "Usage: $0 <project-dir> <template-name>"
    echo ""
    echo "Available templates:"
    echo "  basic-rna       - Standard RNA-seq analysis"
    echo "  multimodal      - RNA + ATAC or CITE-seq"
    echo "  archr-focused   - ArchR scATAC-seq analysis (uses dev-archr service)"
    echo "  example-DMATAC  - Differential chromatin accessibility"
    echo ""
    echo "Example:"
    echo "  $0 ~/projects/my-scrna-analysis basic-rna"
    exit 1
}

if [ $# -ne 2 ]; then
    usage
fi

PROJECT_DIR="$1"
TEMPLATE="$2"
TEMPLATE_PATH="${TEMPLATES_DIR}/${TEMPLATE}"

# Validate template
if [ ! -d "$TEMPLATE_PATH" ]; then
    echo -e "${RED}Error: Template '${TEMPLATE}' not found${NC}"
    echo ""
    usage
fi

# Check if project directory exists
if [ -d "$PROJECT_DIR" ]; then
    echo -e "${YELLOW}Warning: Directory '${PROJECT_DIR}' already exists${NC}"
    read -p "Continue and merge template? (y/N): " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Aborted."
        exit 0
    fi
else
    mkdir -p "$PROJECT_DIR"
fi

echo -e "${GREEN}Initializing project from '${TEMPLATE}' template...${NC}"

# Copy template structure
echo "Copying template files..."
cp -r "${TEMPLATE_PATH}"/* "$PROJECT_DIR/"

# Copy universal .vscode settings
echo "Configuring VS Code settings..."
mkdir -p "${PROJECT_DIR}/.vscode"
cp "${TEMPLATES_DIR}/.vscode/settings.json" "${PROJECT_DIR}/.vscode/settings.json"

# Create standard directories if they don't exist
for dir in data/raw data/processed scripts notebooks results; do
    mkdir -p "${PROJECT_DIR}/${dir}"
done

# Copy devcontainer configuration
echo "Setting up devcontainer configuration..."
mkdir -p "${PROJECT_DIR}/.devcontainer"

# Determine which service to use
if [ "$TEMPLATE" == "archr-focused" ]; then
    SERVICE="dev-archr"
    IMAGE="greenleaflab/archr:1.0.3-base-r4.4"
else
    SERVICE="dev-core"
    IMAGE="scdock-r-dev:v0.5.0"
fi

# Create devcontainer.json
cat > "${PROJECT_DIR}/.devcontainer/devcontainer.json" <<EOF
{
  "name": "scbio-${TEMPLATE}",
  "dockerComposeFile": "docker-compose.yml",
  "service": "${SERVICE}",
  "workspaceFolder": "/workspaces/project",
  "shutdownAction": "stopCompose",
  "customizations": {
    "vscode": {
      "extensions": [
        "REditorSupport.r",
        "ms-python.python",
        "ms-toolsai.jupyter",
        "quarto.quarto"
      ]
    }
  },
  "postStartCommand": "bash .devcontainer/post-start.sh"
}
EOF

# Create docker-compose.yml
cat > "${PROJECT_DIR}/.devcontainer/docker-compose.yml" <<EOF
services:
  dev-core:
    image: scdock-r-dev:v0.5.0
    user: "\${LOCAL_UID:-1000}:\${LOCAL_GID:-1000}"
    working_dir: /workspaces/project
    volumes:
      - \${WORKSPACE_FOLDER:-.}:/workspaces/project
      # Add your data mounts here:
      # - /path/to/data:/workspaces/project/data/raw:ro
    stdin_open: true
    tty: true
    command: /bin/bash

  dev-archr:
    image: greenleaflab/archr:1.0.3-base-r4.4
    user: "\${LOCAL_UID:-1000}:\${LOCAL_GID:-1000}"
    working_dir: /workspaces/project
    volumes:
      - \${WORKSPACE_FOLDER:-.}:/workspaces/project
      # Duplicate the same data mounts as dev-core for consistency
      # - /path/to/data:/workspaces/project/data/raw:ro
    stdin_open: true
    tty: true
    command: /bin/bash
EOF

# Create post-start.sh
cat > "${PROJECT_DIR}/.devcontainer/post-start.sh" <<'EOF'
#!/bin/bash
# Post-start sanity checks

echo "=== Container Environment Check ==="
echo "R version: $(R --version | head -n1)"
echo "Python version: $(python3 --version)"
echo "Working directory: $(pwd)"
echo "User: $(whoami) (UID=$(id -u), GID=$(id -g))"

# Check key packages
if command -v R &> /dev/null; then
    echo -n "Seurat installed: "
    R -q --vanilla -e 'if (requireNamespace("Seurat", quietly=TRUE)) cat("YES\n") else cat("NO\n")' 2>/dev/null | tail -n1
fi

if command -v python3 &> /dev/null; then
    echo -n "scanpy installed: "
    python3 -c "import scanpy; print('YES')" 2>/dev/null || echo "NO"
fi

echo "==================================="
EOF

chmod +x "${PROJECT_DIR}/.devcontainer/post-start.sh"

# Create .env file if it doesn't exist
if [ ! -f "${PROJECT_DIR}/.env" ]; then
    echo "Creating .env file..."
    cat > "${PROJECT_DIR}/.env" <<EOF
# Docker Compose environment variables
LOCAL_UID=$(id -u)
LOCAL_GID=$(id -g)
WORKSPACE_FOLDER=.
EOF
fi

# Create .gitignore if it doesn't exist
if [ ! -f "${PROJECT_DIR}/.gitignore" ]; then
    echo "Creating .gitignore..."
    cat > "${PROJECT_DIR}/.gitignore" <<'EOF'
# R
.Rproj.user
.Rhistory
.RData
.Ruserdata
renv/library/
renv/local/
renv/cellar/
renv/lock/
renv/python/
renv/staging/

# Python
__pycache__/
*.py[cod]
*$py.class
.venv/
venv/
*.egg-info/

# Jupyter
.ipynb_checkpoints/
*.ipynb_checkpoints

# Data (do not commit large files)
data/raw/*
data/processed/*
!data/raw/.gitkeep
!data/processed/.gitkeep

# Results (optional, adjust as needed)
results/figures/*
results/tables/*
!results/figures/.gitkeep
!results/tables/.gitkeep

# IDE
.vscode/*
!.vscode/settings.json
!.vscode/extensions.json
.idea/

# OS
.DS_Store
Thumbs.db

# Environment
.env.local
EOF
fi

# Create .gitkeep files
touch "${PROJECT_DIR}/data/raw/.gitkeep"
touch "${PROJECT_DIR}/data/processed/.gitkeep"
touch "${PROJECT_DIR}/results/.gitkeep"

echo ""
echo -e "${GREEN}✓ Project initialized successfully!${NC}"
echo ""
echo "Next steps:"
echo "  1. cd ${PROJECT_DIR}"
echo "  2. Review and update .devcontainer/docker-compose.yml (add data mounts)"
echo "  3. Open in VS Code: code ${PROJECT_DIR}"
echo "  4. Reopen in container: Ctrl+Shift+P → 'Dev Containers: Reopen in Container'"
echo ""
if [ "$TEMPLATE" == "archr-focused" ]; then
    echo -e "${YELLOW}Note: This project uses the ArchR image (R 4.4).${NC}"
    echo "      To switch between dev-core and dev-archr, see DEVOPS.md"
fi
echo ""
echo "Documentation:"
echo "  - Project structure: ${PROJECT_DIR}/README.md"
echo "  - Build/ops guide: ${SCRIPT_DIR}/DEVOPS.md"
echo "  - Architecture: ${SCRIPT_DIR}/CLAUDE.md"
