#!/usr/bin/env bash
# init-project.sh - Initialize a new project directory from templates
#
# Usage:
#   ./init-project.sh <project-dir> <template-name> [OPTIONS]
#
# Templates:
#   basic-rna       - Standard RNA-seq analysis
#   multimodal      - RNA + ATAC or CITE-seq
#   archr-focused   - ArchR scATAC-seq analysis
#   example-DMATAC  - Differential chromatin accessibility
#
# Options:
#   --data-mount KEY:PATH[:ro]    Add data mount (can be used multiple times)
#   --interactive                  Prompt for all configuration options
#   --git-init                     Initialize git repository
#   --submodules LIST              Add submodules to 01_scripts/ (comma-separated)
#
# Examples:
#   ./init-project.sh ~/projects/my-analysis basic-rna --interactive
#   ./init-project.sh ~/projects/atac-study archr-focused \
#       --data-mount atac:/scratch/data/DT-1234 \
#       --data-mount rna:/scratch/data/DT-5678:ro \
#       --git-init

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEMPLATES_DIR="${SCRIPT_DIR}/templates"

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Default options
INTERACTIVE=false
GIT_INIT=false
declare -a DATA_MOUNTS=()
declare -a SUBMODULES=()

usage() {
    echo "Usage: $0 <project-dir> <template-name> [OPTIONS]"
    echo ""
    echo "Available templates:"
    echo "  basic-rna       - Standard RNA-seq analysis"
    echo "  multimodal      - RNA + ATAC or CITE-seq"
    echo "  archr-focused   - ArchR scATAC-seq analysis (uses dev-archr service)"
    echo "  example-DMATAC  - Differential chromatin accessibility"
    echo ""
    echo "Options:"
    echo "  --data-mount KEY:PATH[:ro]    Add data mount (can be used multiple times)"
    echo "                                 KEY is a label, PATH is host path, :ro for read-only"
    echo "  --interactive                  Prompt for all configuration options"
    echo "  --git-init                     Initialize git repository"
    echo "  --submodules LIST              Add submodules to 01_scripts/ (comma-separated)"
    echo ""
    echo "Examples:"
    echo "  $0 ~/projects/my-analysis basic-rna --interactive"
    echo ""
    echo "  $0 ~/projects/atac-study archr-focused \\"
    echo "      --data-mount atac:/scratch/data/DT-1234 \\"
    echo "      --data-mount rna:/scratch/data/DT-5678:ro \\"
    echo "      --git-init"
    exit 1
}

# Parse command line arguments
if [ $# -lt 2 ]; then
    usage
fi

PROJECT_DIR="$1"
TEMPLATE="$2"
shift 2

# Parse options
while [[ $# -gt 0 ]]; do
    case $1 in
        --data-mount)
            DATA_MOUNTS+=("$2")
            shift 2
            ;;
        --interactive)
            INTERACTIVE=true
            shift
            ;;
        --git-init)
            GIT_INIT=true
            shift
            ;;
        --submodules)
            IFS=',' read -ra SUBMODULES <<< "$2"
            shift 2
            ;;
        *)
            echo -e "${RED}Error: Unknown option '$1'${NC}"
            usage
            ;;
    esac
done

TEMPLATE_PATH="${TEMPLATES_DIR}/${TEMPLATE}"

# Validate template
if [ ! -d "$TEMPLATE_PATH" ]; then
    echo -e "${RED}Error: Template '${TEMPLATE}' not found${NC}"
    echo ""
    usage
fi

# Extract project name from path
PROJECT_NAME=$(basename "$PROJECT_DIR")

# Interactive mode: prompt for options
if [ "$INTERACTIVE" = true ]; then
    echo -e "${BLUE}=== Interactive Configuration ===${NC}"
    echo ""

    # Data mounts
    echo "Configure data mounts (press Enter to skip each):"
    while true; do
        read -p "  Mount label (e.g., 'atac', 'rna'): " mount_label
        if [ -z "$mount_label" ]; then
            break
        fi
        read -p "  Host path (e.g., /scratch/data/DT-1234): " mount_path
        if [ -z "$mount_path" ]; then
            break
        fi
        read -p "  Read-only? (y/N): " -n 1 -r mount_ro
        echo
        if [[ $mount_ro =~ ^[Yy]$ ]]; then
            DATA_MOUNTS+=("${mount_label}:${mount_path}:ro")
        else
            DATA_MOUNTS+=("${mount_label}:${mount_path}")
        fi
    done

    # Git initialization
    if [ "$GIT_INIT" = false ]; then
        read -p "Initialize git repository? (y/N): " -n 1 -r git_reply
        echo
        if [[ $git_reply =~ ^[Yy]$ ]]; then
            GIT_INIT=true
        fi
    fi

    # Submodules
    read -p "Add submodules to 01_scripts/? (comma-separated, or Enter to skip): " submodules_input
    if [ -n "$submodules_input" ]; then
        IFS=',' read -ra SUBMODULES <<< "$submodules_input"
    fi

    echo ""
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

echo -e "${GREEN}Initializing project '${PROJECT_NAME}' from '${TEMPLATE}' template...${NC}"

# Copy template structure (if template has files)
if [ -n "$(ls -A "${TEMPLATE_PATH}" 2>/dev/null)" ]; then
    echo "Copying template files..."
    cp -r "${TEMPLATE_PATH}"/* "$PROJECT_DIR/" 2>/dev/null || true
fi

# Copy universal .vscode settings
echo "Configuring VS Code settings..."
mkdir -p "${PROJECT_DIR}/.vscode"
cp "${TEMPLATES_DIR}/.vscode/settings.json" "${PROJECT_DIR}/.vscode/settings.json"

# Create standard directories with new structure
echo "Creating project structure..."
for dir in 00_data/raw 00_data/processed 00_data/references \
           01_scripts 02_analysis/config 03_results/checkpoints \
           03_results/plots 03_results/tables logs; do
    mkdir -p "${PROJECT_DIR}/${dir}"
done

# Copy documentation templates
echo "Creating documentation..."
if [ -f "${TEMPLATES_DIR}/docs/README.md.template" ]; then
    cp "${TEMPLATES_DIR}/docs/README.md.template" "${PROJECT_DIR}/README.md"
    # Replace placeholders
    sed -i "s|{{PROJECT_NAME}}|${PROJECT_NAME}|g" "${PROJECT_DIR}/README.md"
    sed -i "s|{{PROJECT_PATH}}|${PROJECT_DIR}|g" "${PROJECT_DIR}/README.md"
    sed -i "s|{{DATE}}|$(date +%Y-%m-%d)|g" "${PROJECT_DIR}/README.md"
    sed -i "s|{{TEMPLATE_TYPE}}|${TEMPLATE}|g" "${PROJECT_DIR}/README.md"
    sed -i "s|{{IMAGE_VERSION}}|scdock-r-dev:v0.5.1|g" "${PROJECT_DIR}/README.md"
    sed -i "s|{{SCBIO_DOCKER_PATH}}|${SCRIPT_DIR}|g" "${PROJECT_DIR}/README.md"
fi

if [ -f "${TEMPLATES_DIR}/docs/plan.md.template" ]; then
    cp "${TEMPLATES_DIR}/docs/plan.md.template" "${PROJECT_DIR}/plan.md"
    sed -i "s|{{PROJECT_NAME}}|${PROJECT_NAME}|g" "${PROJECT_DIR}/plan.md"
    sed -i "s|{{DATE}}|$(date +%Y-%m-%d)|g" "${PROJECT_DIR}/plan.md"
fi

if [ -f "${TEMPLATES_DIR}/docs/tasks.md.template" ]; then
    cp "${TEMPLATES_DIR}/docs/tasks.md.template" "${PROJECT_DIR}/tasks.md"
    sed -i "s|{{PROJECT_NAME}}|${PROJECT_NAME}|g" "${PROJECT_DIR}/tasks.md"
    sed -i "s|{{DATE}}|$(date +%Y-%m-%d)|g" "${PROJECT_DIR}/tasks.md"
fi

# Copy .env.example
if [ -f "${TEMPLATES_DIR}/docs/.env.example" ]; then
    cp "${TEMPLATES_DIR}/docs/.env.example" "${PROJECT_DIR}/.env.example"
fi

# Copy Claude Code integration files (if available - dev-claude-integration branch)
if [ -d "${TEMPLATES_DIR}/claude" ]; then
    echo "Setting up Claude Code integration..."

    # Copy CLAUDE.md
    if [ -f "${TEMPLATES_DIR}/claude/CLAUDE.md.template" ]; then
        cp "${TEMPLATES_DIR}/claude/CLAUDE.md.template" "${PROJECT_DIR}/CLAUDE.md"
        sed -i "s|{{PROJECT_NAME}}|${PROJECT_NAME}|g" "${PROJECT_DIR}/CLAUDE.md"
        sed -i "s|{{DATE}}|$(date +%Y-%m-%d)|g" "${PROJECT_DIR}/CLAUDE.md"
        sed -i "s|{{TEMPLATE_TYPE}}|${TEMPLATE}|g" "${PROJECT_DIR}/CLAUDE.md"
        sed -i "s|{{STATUS}}|Planning|g" "${PROJECT_DIR}/CLAUDE.md"
        sed -i "s|{{IMAGE_VERSION}}|scdock-r-dev:v0.5.1|g" "${PROJECT_DIR}/CLAUDE.md"
    fi

    # Copy WORKFLOW.md
    if [ -f "${TEMPLATES_DIR}/claude/WORKFLOW.md" ]; then
        cp "${TEMPLATES_DIR}/claude/WORKFLOW.md" "${PROJECT_DIR}/WORKFLOW.md"
    fi

    # Copy .claude/agents/ directory
    if [ -d "${TEMPLATES_DIR}/claude/.claude" ]; then
        mkdir -p "${PROJECT_DIR}/.claude"
        cp -r "${TEMPLATES_DIR}/claude/.claude"/* "${PROJECT_DIR}/.claude/"
    fi
fi

# Copy devcontainer configuration
echo "Setting up devcontainer configuration..."
mkdir -p "${PROJECT_DIR}/.devcontainer"
mkdir -p "${PROJECT_DIR}/.devcontainer/scripts"

# Determine which service to use
if [ "$TEMPLATE" == "archr-focused" ]; then
    SERVICE="dev-archr"
    IMAGE="scdock-r-archr:v0.5.1"
else
    SERVICE="dev-core"
    IMAGE="scdock-r-dev:v0.5.1"
fi

# Create devcontainer.json
cat > "${PROJECT_DIR}/.devcontainer/devcontainer.json" <<EOF
{
  "name": "${PROJECT_NAME}",
  "dockerComposeFile": "docker-compose.yml",
  "service": "${SERVICE}",
  "workspaceFolder": "/workspaces/${PROJECT_NAME}",
  "remoteUser": "devuser",
  "updateRemoteUserUID": true,
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
  "postStartCommand": "bash .devcontainer/scripts/poststart_sanity.sh"
}
EOF

# Generate data mount volume lines
DATA_MOUNT_LINES=""
if [ ${#DATA_MOUNTS[@]} -gt 0 ]; then
    DATA_MOUNT_LINES+="      # Data mounts"$'\n'
    for mount in "${DATA_MOUNTS[@]}"; do
        IFS=':' read -ra MOUNT_PARTS <<< "$mount"
        mount_label="${MOUNT_PARTS[0]}"
        mount_path="${MOUNT_PARTS[1]}"
        mount_ro="${MOUNT_PARTS[2]:-}"

        if [ -n "$mount_ro" ]; then
            DATA_MOUNT_LINES+="      - ${mount_path}:/workspaces/${PROJECT_NAME}/00_data/${mount_label}:ro"$'\n'
        else
            DATA_MOUNT_LINES+="      - ${mount_path}:/workspaces/${PROJECT_NAME}/00_data/${mount_label}"$'\n'
        fi
    done
else
    DATA_MOUNT_LINES+="      # Add your data mounts here:"$'\n'
    DATA_MOUNT_LINES+="      # - /path/to/data:/workspaces/${PROJECT_NAME}/00_data/raw:ro"$'\n'
fi

# Create docker-compose.yml with dynamic mounts
cat > "${PROJECT_DIR}/.devcontainer/docker-compose.yml" <<EOF
services:
  dev-core:
    image: scdock-r-dev:v0.5.1
    user: "\${LOCAL_UID:-1000}:\${LOCAL_GID:-1000}"
    working_dir: /workspaces/${PROJECT_NAME}
    environment:
      - CONTEXT7_API_KEY=\${CONTEXT7_API_KEY}
    volumes:
      - \${WORKSPACE_FOLDER:-.}:/workspaces/${PROJECT_NAME}
${DATA_MOUNT_LINES}    stdin_open: true
    tty: true
    command: /bin/bash

  dev-archr:
    image: scdock-r-archr:v0.5.1
    user: "\${LOCAL_UID:-1000}:\${LOCAL_GID:-1000}"
    working_dir: /workspaces/${PROJECT_NAME}
    volumes:
      - \${WORKSPACE_FOLDER:-.}:/workspaces/${PROJECT_NAME}
${DATA_MOUNT_LINES}    stdin_open: true
    tty: true
    command: /bin/bash
EOF

# Copy poststart sanity script from scbio-docker repo
if [ -f "${SCRIPT_DIR}/.devcontainer/scripts/poststart_sanity.sh" ]; then
    cp "${SCRIPT_DIR}/.devcontainer/scripts/poststart_sanity.sh" \
       "${PROJECT_DIR}/.devcontainer/scripts/poststart_sanity.sh"
    chmod +x "${PROJECT_DIR}/.devcontainer/scripts/poststart_sanity.sh"
else
    # Create basic sanity check if original doesn't exist
    cat > "${PROJECT_DIR}/.devcontainer/scripts/poststart_sanity.sh" <<'SANITY_EOF'
#!/bin/bash
echo "=== Container Environment Check ==="
echo "R version: $(R --version | head -n1)"
echo "Python version: $(python3 --version)"
echo "Working directory: $(pwd)"
echo "User: $(whoami) (UID=$(id -u), GID=$(id -g))"
echo "==================================="
SANITY_EOF
    chmod +x "${PROJECT_DIR}/.devcontainer/scripts/poststart_sanity.sh"
fi

# Create .env file
if [ ! -f "${PROJECT_DIR}/.env" ]; then
    echo "Creating .env file..."
    cat > "${PROJECT_DIR}/.env" <<EOF
# Docker Compose environment variables
LOCAL_UID=$(id -u)
LOCAL_GID=$(id -g)
WORKSPACE_FOLDER=.

# MCP Server API Keys (optional)
CONTEXT7_API_KEY=
EOF
fi

# Create/update .gitignore
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
00_data/raw/*
00_data/processed/*
!00_data/raw/.gitkeep
!00_data/processed/.gitkeep

# Results
03_results/checkpoints/*
03_results/plots/*
03_results/tables/*
!03_results/checkpoints/.gitkeep
!03_results/plots/.gitkeep
!03_results/tables/.gitkeep

# Logs
logs/*
!logs/.gitkeep

# IDE
.vscode/*
!.vscode/settings.json
!.vscode/extensions.json
.idea/

# OS
.DS_Store
Thumbs.db

# Environment
.env
.env.local
EOF
fi

# Create .gitkeep files
touch "${PROJECT_DIR}/00_data/raw/.gitkeep"
touch "${PROJECT_DIR}/00_data/processed/.gitkeep"
touch "${PROJECT_DIR}/00_data/references/.gitkeep"
touch "${PROJECT_DIR}/03_results/checkpoints/.gitkeep"
touch "${PROJECT_DIR}/03_results/plots/.gitkeep"
touch "${PROJECT_DIR}/03_results/tables/.gitkeep"
touch "${PROJECT_DIR}/logs/.gitkeep"

# Add submodules if requested
if [ ${#SUBMODULES[@]} -gt 0 ]; then
    echo "Adding submodules to 01_scripts/..."
    cd "${PROJECT_DIR}/01_scripts"
    for submodule in "${SUBMODULES[@]}"; do
        # Trim whitespace
        submodule=$(echo "$submodule" | xargs)
        if [ -n "$submodule" ]; then
            echo "  - ${submodule}"
            # You may need to customize URLs based on your organization
            # This is a placeholder - adjust based on your needs
            git submodule add "https://github.com/YOUR_ORG/${submodule}.git" "${submodule}" 2>/dev/null || \
                echo "    (Note: Update submodule URL in git config manually)"
        fi
    done
    cd - > /dev/null
fi

# Git initialization
if [ "$GIT_INIT" = true ]; then
    echo "Initializing git repository..."
    cd "${PROJECT_DIR}"
    if [ ! -d ".git" ]; then
        git init
        git add .
        git commit -m "Initial project structure from ${TEMPLATE} template

Created with scbio-docker init-project.sh
Template: ${TEMPLATE}
Date: $(date +%Y-%m-%d)
"
        echo -e "${GREEN}✓ Git repository initialized${NC}"
    else
        echo -e "${YELLOW}Git repository already exists${NC}"
    fi
    cd - > /dev/null
fi

echo ""
echo -e "${GREEN}✓ Project initialized successfully!${NC}"
echo ""
echo -e "${BLUE}Project Summary:${NC}"
echo "  Name: ${PROJECT_NAME}"
echo "  Location: ${PROJECT_DIR}"
echo "  Template: ${TEMPLATE}"
echo "  Container service: ${SERVICE}"
if [ ${#DATA_MOUNTS[@]} -gt 0 ]; then
    echo "  Data mounts configured: ${#DATA_MOUNTS[@]}"
fi
if [ "$GIT_INIT" = true ]; then
    echo "  Git: initialized"
fi
echo ""
echo -e "${BLUE}Next steps:${NC}"
echo "  1. cd ${PROJECT_DIR}"
if [ ${#DATA_MOUNTS[@]} -eq 0 ]; then
    echo "  2. Edit .devcontainer/docker-compose.yml to add data mounts"
    echo "  3. Open in VS Code: code ${PROJECT_DIR}"
    echo "  4. Reopen in container: Ctrl+Shift+P → 'Dev Containers: Reopen in Container'"
else
    echo "  2. Open in VS Code: code ${PROJECT_DIR}"
    echo "  3. Reopen in container: Ctrl+Shift+P → 'Dev Containers: Reopen in Container'"
fi
echo "  $(( ${#DATA_MOUNTS[@]} > 0 ? 5 : 5 )). Review and fill in plan.md with your scientific question"
echo "  $(( ${#DATA_MOUNTS[@]} > 0 ? 6 : 6 )). Create detailed tasks in tasks.md"
echo ""
if [ "$TEMPLATE" == "archr-focused" ]; then
    echo -e "${YELLOW}Note: This project uses the ArchR image (R 4.4).${NC}"
    echo "      To switch between dev-core and dev-archr, edit .devcontainer/devcontainer.json"
fi
echo ""
echo -e "${BLUE}Documentation:${NC}"
echo "  - Project workflow: ${PROJECT_DIR}/README.md"
echo "  - Analysis plan: ${PROJECT_DIR}/plan.md"
echo "  - Task tracker: ${PROJECT_DIR}/tasks.md"
echo "  - Image build guide: ${SCRIPT_DIR}/DEVOPS.md"
echo "  - Architecture: ${SCRIPT_DIR}/CLAUDE.md"
echo ""
