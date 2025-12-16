#!/usr/bin/env bash
# init-project.sh - Initialize a new scbio-dock project directory from templates
#
# Usage:
#   ./init-project.sh <project-dir> [template-name] [OPTIONS]
#
# Templates:
#   base            - Standard bioinformatics project (default)
#
# Options:
#   --data-mount KEY:PATH[:ro]    Add data mount (can be used multiple times)
#   --interactive                  Prompt for all configuration options
#   --git-init                     Initialize git repository
#   --with-submodules              Add RNAseq-toolkit and SciAgent-toolkit as git submodules
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
WITH_SUBMODULES=false
declare -a DATA_MOUNTS=()

# Submodule URLs (SSH - gh CLI handles authentication)
RNASEQ_TOOLKIT_URL="git@github.com:tony-zhelonkin/RNAseq-toolkit.git"
RNASEQ_TOOLKIT_BRANCH="dev"
SCIAGENT_TOOLKIT_URL="git@github.com:tony-zhelonkin/SciAgent-toolkit.git"
SCIAGENT_TOOLKIT_BRANCH="main"

usage() {
    echo "Usage: $0 <project-dir> [template-name] [OPTIONS]"
    echo ""
    echo "Available templates:"
    echo "  base            - Standard bioinformatics project (default)"
    echo ""
    echo "Options:"
    echo "  --data-mount KEY:PATH[:ro]    Add data mount (can be used multiple times)"
    echo "                                 KEY is a label, PATH is host path, :ro for read-only"
    echo "  --interactive                  Prompt for all configuration options"
    echo "  --git-init                     Initialize git repository"
    echo "  --with-submodules              Add RNAseq-toolkit and SciAgent-toolkit as git submodules"
    echo "                                 (implies --git-init, requires SSH key setup for GitHub)"
    echo ""
    echo "Examples:"
    echo "  $0 ~/projects/my-analysis"
    echo "  $0 ~/projects/my-analysis base --interactive"
    echo "  $0 ~/projects/atac-study base \\"
    echo "      --data-mount atac:/scratch/data/DT-1234 \\"
    echo "      --data-mount rna:/scratch/data/DT-5678:ro \\"
    echo "      --git-init --with-submodules"
    exit 1
}

# Parse command line arguments
if [ $# -lt 1 ]; then
    usage
fi

PROJECT_DIR="$1"
shift

# Template is optional, default to "base"
TEMPLATE="base"
if [ $# -gt 0 ] && [[ ! "$1" =~ ^-- ]]; then
    TEMPLATE="$1"
    shift
fi

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
        --with-submodules)
            WITH_SUBMODULES=true
            GIT_INIT=true  # Submodules require git
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

    # Submodules (only ask if git will be initialized)
    if [ "$GIT_INIT" = true ] && [ "$WITH_SUBMODULES" = false ]; then
        read -p "Add analysis toolkits as git submodules? (Y/n): " -n 1 -r submod_reply
        echo
        if [[ ! $submod_reply =~ ^[Nn]$ ]]; then
            WITH_SUBMODULES=true
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

# Species defaults (can be overridden by interactive mode)
: ${SPECIES:="Mus musculus"}
: ${SPECIES_DB:="MM"}
: ${GENOME_BUILD:="mm10"}

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

# Service and image defaults (ArchR is available via docker-compose profile)
SERVICE="dev-core"
IMAGE="scdock-r-dev:v0.5.2"

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
           01_modules 02_analysis/config 02_analysis/helpers \
           03_results/checkpoints 03_results/plots 03_results/tables logs; do
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

# Note: plan.md has been replaced by context.md, which is created by
# SciAgent-toolkit/scripts/setup-ai.sh along with other AI context files

if [ -f "${TEMPLATES_DIR}/docs/tasks.md.template" ]; then
    cp "${TEMPLATES_DIR}/docs/tasks.md.template" "${PROJECT_DIR}/tasks.md"
    sed -i "s|{{PROJECT_NAME}}|${PROJECT_NAME}|g" "${PROJECT_DIR}/tasks.md"
    sed -i "s|{{DATE}}|$(date +%Y-%m-%d)|g" "${PROJECT_DIR}/tasks.md"
fi

# Copy .env.example
if [ -f "${TEMPLATES_DIR}/docs/.env.example" ]; then
    cp "${TEMPLATES_DIR}/docs/.env.example" "${PROJECT_DIR}/.env.example"
fi

# Note: CLAUDE.md is now created by SciAgent-toolkit/scripts/setup-ai.sh
# This ensures AI context files are created when AI tools are set up, not at project init

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

# Note: Analysis guidelines are now in SciAgent-toolkit/docs/guidelines/
# They are not copied to each project - reference them directly from the submodule
# This ensures single source of truth and easier updates

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
    # Uncomment to expose httpgd graphics server (may conflict if port already in use)
    # ports:
    #   - "8787:8787"
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
    # Uncomment to expose httpgd graphics server (may conflict if port already in use)
    # ports:
    #   - "8787:8787"
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

# Add git submodules (requires git to be initialized)
if [ "$WITH_SUBMODULES" = true ] && [ -d "${PROJECT_DIR}/.git" ]; then
    echo "Adding analysis toolkits as git submodules..."
    cd "${PROJECT_DIR}"

    SUBMODULES_ADDED=false

    # Helper function to add submodule with gh fallback
    add_submodule_with_fallback() {
        local repo_owner="tony-zhelonkin"
        local repo_name="$1"
        local branch="$2"
        local target_path="$3"
        local ssh_url="git@github.com:${repo_owner}/${repo_name}.git"
        local https_url="https://github.com/${repo_owner}/${repo_name}.git"

        if [ -d "$target_path" ]; then
            echo -e "  ${YELLOW}⚠ ${repo_name} directory already exists, skipping${NC}"
            return 1
        fi

        # Try direct git submodule add first (uses SSH)
        if git submodule add -b "$branch" "$ssh_url" "$target_path" 2>/dev/null; then
            echo -e "  ${GREEN}✓ ${repo_name} added via git (SSH)${NC}"
            return 0
        fi

        # Fallback: use gh CLI with HTTPS (for agents without SSH access)
        echo "  SSH failed, trying gh CLI fallback (HTTPS)..."
        if command -v gh &> /dev/null && gh auth status &> /dev/null; then
            # Clone using gh with explicit HTTPS protocol
            if GIT_CONFIG_COUNT=1 \
               GIT_CONFIG_KEY_0="url.https://github.com/.insteadOf" \
               GIT_CONFIG_VALUE_0="git@github.com:" \
               gh repo clone "${repo_owner}/${repo_name}" "$target_path" -- -b "$branch" 2>/dev/null; then
                # Manually create .gitmodules entry (keep SSH URL for future clones by user)
                git config -f .gitmodules "submodule.${target_path}.path" "$target_path"
                git config -f .gitmodules "submodule.${target_path}.url" "$ssh_url"
                git config -f .gitmodules "submodule.${target_path}.branch" "$branch"
                # Register in .git/config
                git config "submodule.${target_path}.url" "$ssh_url"
                git config "submodule.${target_path}.active" "true"
                # Stage the submodule
                git add "$target_path"
                echo -e "  ${GREEN}✓ ${repo_name} added via gh CLI (HTTPS)${NC}"
                return 0
            fi
        fi

        echo -e "  ${RED}✗ Failed to add ${repo_name}${NC}"
        return 1
    }

    # Add RNAseq-toolkit submodule (dev branch)
    echo "  Adding RNAseq-toolkit (branch: ${RNASEQ_TOOLKIT_BRANCH})..."
    if add_submodule_with_fallback "RNAseq-toolkit" "${RNASEQ_TOOLKIT_BRANCH}" "01_modules/RNAseq-toolkit"; then
        SUBMODULES_ADDED=true
    fi

    # Add SciAgent-toolkit submodule
    echo "  Adding SciAgent-toolkit (branch: ${SCIAGENT_TOOLKIT_BRANCH})..."
    if add_submodule_with_fallback "SciAgent-toolkit" "${SCIAGENT_TOOLKIT_BRANCH}" "01_modules/SciAgent-toolkit"; then
        SUBMODULES_ADDED=true
    fi

    # Commit submodule additions
    if [ "$SUBMODULES_ADDED" = true ]; then
        if [ -f ".gitmodules" ]; then
            git add .gitmodules 01_modules/
            git commit -m "Add analysis toolkits as git submodules

- RNAseq-toolkit (${RNASEQ_TOOLKIT_BRANCH} branch): Reusable RNA-seq analysis functions
- SciAgent-toolkit: AI infrastructure and MCP server setup
"
            echo -e "${GREEN}✓ Submodules committed${NC}"
        fi
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
if [ "$WITH_SUBMODULES" = true ]; then
    echo "  Submodules: RNAseq-toolkit (${RNASEQ_TOOLKIT_BRANCH}), SciAgent-toolkit"
fi
echo ""
echo -e "${BLUE}Next steps:${NC}"
echo "  1. cd ${PROJECT_DIR}"
if [ ${#DATA_MOUNTS[@]} -eq 0 ]; then
    echo "  2. Edit .devcontainer/docker-compose.yml to add data mounts"
    echo "  3. Open in VS Code: code ${PROJECT_DIR}"
    echo "  4. Reopen in container: Ctrl+Shift+P → 'Dev Containers: Reopen in Container'"
    echo "  5. Run AI setup: ./01_modules/SciAgent-toolkit/scripts/setup-ai.sh"
    echo "  6. Fill in context.md with your scientific question"
    echo "  7. Edit 02_analysis/config/pipeline.yaml for your experiment"
else
    echo "  2. Open in VS Code: code ${PROJECT_DIR}"
    echo "  3. Reopen in container: Ctrl+Shift+P → 'Dev Containers: Reopen in Container'"
    echo "  4. Run AI setup: ./01_modules/SciAgent-toolkit/scripts/setup-ai.sh"
    echo "  5. Fill in context.md with your scientific question"
    echo "  6. Edit 02_analysis/config/pipeline.yaml for your experiment"
fi
echo ""
echo -e "${BLUE}Documentation:${NC}"
echo "  - Project workflow: ${PROJECT_DIR}/README.md"
echo "  - Task tracker: ${PROJECT_DIR}/tasks.md"
echo "  - Research notes: ${PROJECT_DIR}/notes.md"
echo "  - Config: ${PROJECT_DIR}/02_analysis/config/"
echo ""
echo -e "${BLUE}After running setup-ai.sh:${NC}"
echo "  - AI context: CLAUDE.md, GEMINI.md, AGENTS.md"
echo "  - Scientific context: context.md"
echo "  - Methodology: 01_modules/SciAgent-toolkit/docs/guidelines/"
echo ""
echo -e "${BLUE}scbio-docker references:${NC}"
echo "  - Image build guide: ${SCRIPT_DIR}/DEVOPS.md"
echo "  - Architecture: ${SCRIPT_DIR}/CLAUDE.md"
echo ""
