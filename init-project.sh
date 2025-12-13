#!/usr/bin/env bash
# init-project.sh - Initialize a new scbio-dock project directory from templates
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
#
# Examples:
#   ./init-project.sh ~/projects/my-analysis basic-rna --interactive
#   ./init-project.sh ~/projects/atac-study archr-focused \
#       --data-mount atac:/scratch/data/DT-1234 \
#       --data-mount rna:/scratch/data/DT-5678:ro \
#       --git-init

set -euo pipefail

# Resolve symlinks to get actual script location
SCRIPT_PATH="${BASH_SOURCE[0]}"
while [ -L "$SCRIPT_PATH" ]; do
    SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"
    SCRIPT_PATH="$(readlink "$SCRIPT_PATH")"
    [[ $SCRIPT_PATH != /* ]] && SCRIPT_PATH="$SCRIPT_DIR/$SCRIPT_PATH"
done
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"
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

    # Resource limits
    read -p "Max CPUs (default: 50): " max_cpus
    MAX_CPUS="${max_cpus:-50}"
    read -p "Max Memory (default: 450G): " max_memory
    MAX_MEMORY="${max_memory:-450G}"

    # Species configuration
    echo ""
    echo "Species configuration:"
    read -p "Species (default: Mus musculus): " species_input
    SPECIES="${species_input:-Mus musculus}"

    # Auto-derive SPECIES_DB and GENOME_BUILD based on common species
    case "$SPECIES" in
        "Mus musculus"|"mouse"|"Mouse")
            SPECIES="Mus musculus"
            SPECIES_DB="MM"
            GENOME_BUILD="${GENOME_BUILD:-mm10}"
            ;;
        "Homo sapiens"|"human"|"Human")
            SPECIES="Homo sapiens"
            SPECIES_DB="HS"
            GENOME_BUILD="${GENOME_BUILD:-hg38}"
            ;;
        *)
            read -p "Species DB code (e.g., MM for mouse, HS for human): " SPECIES_DB
            read -p "Genome build (e.g., mm10, hg38): " GENOME_BUILD
            ;;
    esac

    echo ""
fi

# Non-interactive defaults (if not set by interactive mode)
: ${MAX_CPUS:=50}
: ${MAX_MEMORY:=450G}

# Species defaults based on template (can be overridden by interactive mode)
case "$TEMPLATE" in
    basic-rna|multimodal)
        : ${SPECIES:="Mus musculus"}
        : ${SPECIES_DB:="MM"}
        : ${GENOME_BUILD:="mm10"}
        ;;
    archr-focused|example-DMATAC)
        : ${SPECIES:="Mus musculus"}
        : ${SPECIES_DB:="MM"}
        : ${GENOME_BUILD:="mm10"}
        ;;
    *)
        : ${SPECIES:="Mus musculus"}
        : ${SPECIES_DB:="MM"}
        : ${GENOME_BUILD:="mm10"}
        ;;
esac

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

# Determine which service and image to use (needed for documentation templates)
if [ "$TEMPLATE" == "archr-focused" ]; then
    SERVICE="dev-archr"
    IMAGE="scdock-r-archr:v0.5.2"
else
    SERVICE="dev-core"
    IMAGE="scdock-r-dev:v0.5.2"
fi

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
    sed -i "s|{{IMAGE_VERSION}}|${IMAGE}|g" "${PROJECT_DIR}/README.md"
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

# Copy CLAUDE.md template (AI assistant context)
if [ -f "${TEMPLATES_DIR}/docs/CLAUDE.md.template" ]; then
    echo "Creating CLAUDE.md..."
    cp "${TEMPLATES_DIR}/docs/CLAUDE.md.template" "${PROJECT_DIR}/CLAUDE.md"
    sed -i "s|{{PROJECT_NAME}}|${PROJECT_NAME}|g" "${PROJECT_DIR}/CLAUDE.md"
    sed -i "s|{{DATE}}|$(date +%Y-%m-%d)|g" "${PROJECT_DIR}/CLAUDE.md"
    sed -i "s|{{TEMPLATE_TYPE}}|${TEMPLATE}|g" "${PROJECT_DIR}/CLAUDE.md"
    sed -i "s|{{SPECIES}}|${SPECIES}|g" "${PROJECT_DIR}/CLAUDE.md"
    sed -i "s|{{SPECIES_DB}}|${SPECIES_DB}|g" "${PROJECT_DIR}/CLAUDE.md"
    sed -i "s|{{GENOME_BUILD}}|${GENOME_BUILD}|g" "${PROJECT_DIR}/CLAUDE.md"
fi

# Copy notes.md template (research findings tracker)
if [ -f "${TEMPLATES_DIR}/docs/notes.md.template" ]; then
    echo "Creating notes.md..."
    cp "${TEMPLATES_DIR}/docs/notes.md.template" "${PROJECT_DIR}/notes.md"
    sed -i "s|{{PROJECT_NAME}}|${PROJECT_NAME}|g" "${PROJECT_DIR}/notes.md"
    sed -i "s|{{DATE}}|$(date +%Y-%m-%d)|g" "${PROJECT_DIR}/notes.md"
    sed -i "s|{{TEMPLATE_TYPE}}|${TEMPLATE}|g" "${PROJECT_DIR}/notes.md"
    sed -i "s|{{SPECIES}}|${SPECIES}|g" "${PROJECT_DIR}/notes.md"
    sed -i "s|{{SPECIES_DB}}|${SPECIES_DB}|g" "${PROJECT_DIR}/notes.md"
    sed -i "s|{{GENOME_BUILD}}|${GENOME_BUILD}|g" "${PROJECT_DIR}/notes.md"
fi

# Copy configuration templates
echo "Creating configuration files..."
if [ -f "${TEMPLATES_DIR}/config/config.R.template" ]; then
    cp "${TEMPLATES_DIR}/config/config.R.template" "${PROJECT_DIR}/02_analysis/config/config.R"
    sed -i "s|{{PROJECT_NAME}}|${PROJECT_NAME}|g" "${PROJECT_DIR}/02_analysis/config/config.R"
    sed -i "s|{{DATE}}|$(date +%Y-%m-%d)|g" "${PROJECT_DIR}/02_analysis/config/config.R"
    sed -i "s|{{TEMPLATE_TYPE}}|${TEMPLATE}|g" "${PROJECT_DIR}/02_analysis/config/config.R"
    sed -i "s|{{SPECIES}}|${SPECIES}|g" "${PROJECT_DIR}/02_analysis/config/config.R"
    sed -i "s|{{SPECIES_DB}}|${SPECIES_DB}|g" "${PROJECT_DIR}/02_analysis/config/config.R"
    sed -i "s|{{GENOME_BUILD}}|${GENOME_BUILD}|g" "${PROJECT_DIR}/02_analysis/config/config.R"
fi

if [ -f "${TEMPLATES_DIR}/config/pipeline.yaml.template" ]; then
    cp "${TEMPLATES_DIR}/config/pipeline.yaml.template" "${PROJECT_DIR}/02_analysis/config/pipeline.yaml"
    sed -i "s|{{PROJECT_NAME}}|${PROJECT_NAME}|g" "${PROJECT_DIR}/02_analysis/config/pipeline.yaml"
    sed -i "s|{{DATE}}|$(date +%Y-%m-%d)|g" "${PROJECT_DIR}/02_analysis/config/pipeline.yaml"
    sed -i "s|{{SPECIES}}|${SPECIES}|g" "${PROJECT_DIR}/02_analysis/config/pipeline.yaml"
    sed -i "s|{{SPECIES_DB}}|${SPECIES_DB}|g" "${PROJECT_DIR}/02_analysis/config/pipeline.yaml"
    sed -i "s|{{GENOME_BUILD}}|${GENOME_BUILD}|g" "${PROJECT_DIR}/02_analysis/config/pipeline.yaml"
fi

if [ -f "${TEMPLATES_DIR}/config/color_config.R.template" ]; then
    cp "${TEMPLATES_DIR}/config/color_config.R.template" "${PROJECT_DIR}/02_analysis/config/color_config.R"
    sed -i "s|{{DATE}}|$(date +%Y-%m-%d)|g" "${PROJECT_DIR}/02_analysis/config/color_config.R"
fi

# Copy guidelines to project (from canonical source in docs/)
echo "Copying analysis guidelines..."
mkdir -p "${PROJECT_DIR}/docs/guidelines"
if [ -f "${SCRIPT_DIR}/docs/guidelines/rnaseq-analysis-guidelines.md" ]; then
    cp "${SCRIPT_DIR}/docs/guidelines/rnaseq-analysis-guidelines.md" \
       "${PROJECT_DIR}/docs/guidelines/rnaseq-analysis-guidelines.md"
fi

# Copy devcontainer configuration
echo "Setting up devcontainer configuration..."
mkdir -p "${PROJECT_DIR}/.devcontainer"
mkdir -p "${PROJECT_DIR}/.devcontainer/scripts"

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
  "settings": {
    "files.associations": {
      "*.Rmd": "rmd"
    }
  },
  "customizations": {
    "vscode": {
      "extensions": [
        "rdebugger.r-debugger",
        "reditorsupport.r",
        "quarto.quarto",
        "purocean.drawio-preview",
        "redhat.vscode-yaml",
        "yzhang.markdown-all-in-one",
        "ms-azuretools.vscode-docker",
        "ms-vscode-remote.remote-containers",
        "ms-python.python",
        "ms-toolsai.jupyter"
      ]
    }
  },
  "postStartCommand": "bash -lc 'chmod +x .devcontainer/scripts/poststart_sanity.sh && .devcontainer/scripts/poststart_sanity.sh'"
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
    image: ${IMAGE}
    user: "\${LOCAL_UID:-1000}:\${LOCAL_GID:-1000}"
    working_dir: /workspaces/${PROJECT_NAME}
    env_file: .env
    environment:
      - CONTEXT7_API_KEY=\${CONTEXT7_API_KEY}
    ports:
      - "8787:8787"  # httpgd graphics server
    volumes:
      - \${WORKSPACE_FOLDER:-.}:/workspaces/${PROJECT_NAME}
${DATA_MOUNT_LINES}    stdin_open: true
    tty: true
    command: /bin/bash
    deploy:
      resources:
        limits:
          cpus: '\${MAX_CPUS:-${MAX_CPUS}}'
          memory: \${MAX_MEMORY:-${MAX_MEMORY}}
        reservations:
          cpus: '2'
          memory: 8G

  dev-archr:
    profiles: ["archr"]  # Only started with: docker compose --profile archr up
    image: scdock-r-archr:v0.5.2
    user: "\${LOCAL_UID:-1000}:\${LOCAL_GID:-1000}"
    working_dir: /workspaces/${PROJECT_NAME}
    env_file: .env
    environment:
      - CONTEXT7_API_KEY=\${CONTEXT7_API_KEY}
    ports:
      - "8787:8787"  # httpgd graphics server
    volumes:
      - \${WORKSPACE_FOLDER:-.}:/workspaces/${PROJECT_NAME}
${DATA_MOUNT_LINES}    stdin_open: true
    tty: true
    command: /bin/bash
    deploy:
      resources:
        limits:
          cpus: '\${MAX_CPUS:-${MAX_CPUS}}'
          memory: \${MAX_MEMORY:-${MAX_MEMORY}}
        reservations:
          cpus: '2'
          memory: 8G
EOF

# Copy poststart sanity script from scbio-docker repo
if [ -f "${SCRIPT_DIR}/scripts/poststart_sanity.sh" ]; then
    cp "${SCRIPT_DIR}/scripts/poststart_sanity.sh" \
       "${PROJECT_DIR}/.devcontainer/scripts/poststart_sanity.sh"
    chmod +x "${PROJECT_DIR}/.devcontainer/scripts/poststart_sanity.sh"
elif [ -f "${SCRIPT_DIR}/.devcontainer/scripts/poststart_sanity.sh" ]; then
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

# Create .env file in .devcontainer/
if [ ! -f "${PROJECT_DIR}/.devcontainer/.env" ]; then
    echo "Creating .env file..."
    cat > "${PROJECT_DIR}/.devcontainer/.env" <<EOF
# Docker Compose environment variables
LOCAL_UID=$(id -u)
LOCAL_GID=$(id -g)
WORKSPACE_FOLDER=..

# MCP Server API Keys (optional - only if using Claude Code)

# Context7 - Up-to-date library docs (optional - works without key)
# Get key for higher rate limits: https://context7.com/dashboard
CONTEXT7_API_KEY=

# PAL MCP Server - Multi-model AI collaboration
# Requires at least one key to function. Get keys from:
#   - Gemini: https://aistudio.google.com/apikey
#   - OpenAI: https://platform.openai.com/api-keys
GEMINI_API_KEY=
OPENAI_API_KEY=

# Resource Limits (adjust based on your system)
MAX_CPUS=${MAX_CPUS}
MAX_MEMORY=${MAX_MEMORY}
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
.devcontainer/.env
.devcontainer/.env.local
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
echo "  Species: ${SPECIES} (${SPECIES_DB}, ${GENOME_BUILD})"
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
    echo "  5. Review and fill in plan.md with your scientific question"
    echo "  6. Customize CLAUDE.md with project-specific context"
    echo "  7. Edit 02_analysis/config/pipeline.yaml for your experiment"
else
    echo "  2. Open in VS Code: code ${PROJECT_DIR}"
    echo "  3. Reopen in container: Ctrl+Shift+P → 'Dev Containers: Reopen in Container'"
    echo "  4. Review and fill in plan.md with your scientific question"
    echo "  5. Customize CLAUDE.md with project-specific context"
    echo "  6. Edit 02_analysis/config/pipeline.yaml for your experiment"
fi
echo ""
if [ "$TEMPLATE" == "archr-focused" ]; then
    echo -e "${YELLOW}Note: This project uses the ArchR image (R 4.4).${NC}"
    echo "      To switch between dev-core and dev-archr, edit .devcontainer/devcontainer.json"
fi
echo ""
echo -e "${BLUE}Documentation:${NC}"
echo "  - Project workflow: ${PROJECT_DIR}/README.md"
echo "  - AI context: ${PROJECT_DIR}/CLAUDE.md"
echo "  - Analysis plan: ${PROJECT_DIR}/plan.md"
echo "  - Task tracker: ${PROJECT_DIR}/tasks.md"
echo "  - Research notes: ${PROJECT_DIR}/notes.md"
echo "  - Analysis guidelines: ${PROJECT_DIR}/docs/guidelines/rnaseq-analysis-guidelines.md"
echo "  - Config: ${PROJECT_DIR}/02_analysis/config/"
echo ""
echo -e "${BLUE}scbio-docker references:${NC}"
echo "  - Image build guide: ${SCRIPT_DIR}/DEVOPS.md"
echo "  - Architecture: ${SCRIPT_DIR}/CLAUDE.md"
echo ""
