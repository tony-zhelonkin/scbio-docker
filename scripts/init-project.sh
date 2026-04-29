#!/usr/bin/env bash
# init-project.sh - Initialize a new scbio-dock project directory from templates
#
# Usage:
#   ./init-project.sh <project-dir> [template-name] [OPTIONS]
#
# Templates:
#   base            - Standard bioinformatics project (default, only one supported)
#
# Options:
#   --data-mount KEY:PATH[:ro]    Add data mount (can be used multiple times)
#   --interactive                  Prompt for all configuration options
#   --git-init                     Initialize git repository
#   --with-submodules              Add RNAseq-toolkit and SciAgent-toolkit as git submodules

set -euo pipefail

# Resolve symlinks to get actual script location
SCRIPT_PATH="${BASH_SOURCE[0]}"
while [ -L "$SCRIPT_PATH" ]; do
    SCRIPT_DIR_TMP="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"
    SCRIPT_PATH="$(readlink "$SCRIPT_PATH")"
    [[ $SCRIPT_PATH != /* ]] && SCRIPT_PATH="$SCRIPT_DIR_TMP/$SCRIPT_PATH"
done
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"

# Repo root is one level up from scripts/
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
TEMPLATES_DIR="${REPO_ROOT}/templates"

# Read version from VERSION file (single source of truth)
IMAGE_VERSION="$(tr -d '[:space:]' < "$REPO_ROOT/VERSION")"

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
    echo "  --image-version vX.Y.Z         Override image tag (default: read from VERSION file)"
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

# Only "base" is supported; reject anything else cleanly
if [ "$TEMPLATE" != "base" ]; then
    echo "Only 'base' template is currently supported." >&2
    exit 1
fi

# Parse options
while [[ $# -gt 0 ]]; do
    case $1 in
        --data-mount)
            DATA_MOUNTS+=("$2")
            shift 2
            ;;
        --image-version)
            IMAGE_VERSION="$2"
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

# Validate template exists on disk
if [ ! -d "$TEMPLATE_PATH" ]; then
    echo -e "${RED}Error: Template directory '${TEMPLATE_PATH}' not found${NC}"
    exit 1
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

echo -e "${GREEN}Initializing project '${PROJECT_NAME}' from '${TEMPLATE}' template (image ${IMAGE_VERSION})...${NC}"

# Service default (ArchR is available via docker-compose profile)
SERVICE="dev-core"

# Create standard directories matching the new universal template tree
echo "Creating project structure..."
for dir in 00_data/raw 00_data/processed 00_data/references \
           01_modules/.ref \
           02_analysis/config 02_analysis/helpers 02_analysis/scripts 02_analysis/notebooks \
           03_results/checkpoints 03_results/plots 03_results/tables \
           docs/raw docs/ai-generated/vignettes docs/ai-generated/research \
           docs/plan/phase-1 \
           logs; do
    mkdir -p "${PROJECT_DIR}/${dir}"
done

# Copy universal .vscode settings
echo "Configuring VS Code settings..."
mkdir -p "${PROJECT_DIR}/.vscode"
cp "${TEMPLATE_PATH}/.vscode/settings.json" "${PROJECT_DIR}/.vscode/settings.json"

# Copy project .gitignore from template
echo "Creating .gitignore..."
cp "${TEMPLATE_PATH}/.gitignore" "${PROJECT_DIR}/.gitignore"

# Generate project docs from templates
echo "Creating documentation..."
sed -e "s|{{PROJECT_NAME}}|${PROJECT_NAME}|g" \
    -e "s|{{DATE}}|$(date +%Y-%m-%d)|g" \
    "${TEMPLATE_PATH}/docs/README.md.template" > "${PROJECT_DIR}/docs/README.md"

sed -e "s|{{PROJECT_NAME}}|${PROJECT_NAME}|g" \
    -e "s|{{DATE}}|$(date +%Y-%m-%d)|g" \
    "${TEMPLATE_PATH}/docs/plan/README.md.template" > "${PROJECT_DIR}/docs/plan/README.md"

# Top-level project README (from templates/base/README.md if present)
if [ -f "${TEMPLATE_PATH}/README.md" ]; then
    sed -e "s|{{PROJECT_NAME}}|${PROJECT_NAME}|g" \
        -e "s|{{PROJECT_PATH}}|${PROJECT_DIR}|g" \
        -e "s|{{DATE}}|$(date +%Y-%m-%d)|g" \
        -e "s|{{TEMPLATE_TYPE}}|${TEMPLATE}|g" \
        -e "s|{{IMAGE_VERSION}}|scdock-r-dev:${IMAGE_VERSION}|g" \
        -e "s|{{SCBIO_DOCKER_PATH}}|${REPO_ROOT}|g" \
        "${TEMPLATE_PATH}/README.md" > "${PROJECT_DIR}/README.md"
fi

# Optional notes.md (research findings tracker)
if [ -f "${TEMPLATE_PATH}/docs/notes.md.template" ]; then
    echo "Creating notes.md..."
    sed -e "s|{{PROJECT_NAME}}|${PROJECT_NAME}|g" \
        -e "s|{{DATE}}|$(date +%Y-%m-%d)|g" \
        -e "s|{{TEMPLATE_TYPE}}|${TEMPLATE}|g" \
        -e "s|{{SPECIES}}|${SPECIES}|g" \
        -e "s|{{SPECIES_DB}}|${SPECIES_DB}|g" \
        -e "s|{{GENOME_BUILD}}|${GENOME_BUILD}|g" \
        "${TEMPLATE_PATH}/docs/notes.md.template" > "${PROJECT_DIR}/notes.md"
fi

# Copy configuration templates from new flat template tree
echo "Creating configuration files..."
CFG_SRC="${TEMPLATE_PATH}/02_analysis/config"
CFG_DST="${PROJECT_DIR}/02_analysis/config"

if [ -f "${CFG_SRC}/config.R.template" ]; then
    sed -e "s|{{PROJECT_NAME}}|${PROJECT_NAME}|g" \
        -e "s|{{DATE}}|$(date +%Y-%m-%d)|g" \
        -e "s|{{TEMPLATE_TYPE}}|${TEMPLATE}|g" \
        -e "s|{{SPECIES}}|${SPECIES}|g" \
        -e "s|{{SPECIES_DB}}|${SPECIES_DB}|g" \
        -e "s|{{GENOME_BUILD}}|${GENOME_BUILD}|g" \
        "${CFG_SRC}/config.R.template" > "${CFG_DST}/config.R"
fi

if [ -f "${CFG_SRC}/pipeline.yaml.template" ]; then
    sed -e "s|{{PROJECT_NAME}}|${PROJECT_NAME}|g" \
        -e "s|{{DATE}}|$(date +%Y-%m-%d)|g" \
        -e "s|{{SPECIES}}|${SPECIES}|g" \
        -e "s|{{SPECIES_DB}}|${SPECIES_DB}|g" \
        -e "s|{{GENOME_BUILD}}|${GENOME_BUILD}|g" \
        "${CFG_SRC}/pipeline.yaml.template" > "${CFG_DST}/pipeline.yaml"
fi

if [ -f "${CFG_SRC}/color_config.R.template" ]; then
    sed -e "s|{{DATE}}|$(date +%Y-%m-%d)|g" \
        "${CFG_SRC}/color_config.R.template" > "${CFG_DST}/color_config.R"
fi

# Set up devcontainer configuration
echo "Setting up devcontainer configuration..."
mkdir -p "${PROJECT_DIR}/.devcontainer"
mkdir -p "${PROJECT_DIR}/.devcontainer/scripts"

# Render devcontainer.json from template
sed -e "s|{{PROJECT_NAME}}|${PROJECT_NAME}|g" \
    -e "s|{{SERVICE}}|${SERVICE}|g" \
    "${TEMPLATE_PATH}/.devcontainer/devcontainer.json.template" \
    > "${PROJECT_DIR}/.devcontainer/devcontainer.json"

# Build the dynamic data-mount block for compose YAML
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

# Render docker-compose.yml from template
# Use python for the {{DATA_MOUNTS}} substitution (multi-line, sed-hostile).
# All other tokens are simple single-line replacements.
python3 - "$TEMPLATE_PATH" "$PROJECT_DIR" "$IMAGE_VERSION" "$PROJECT_NAME" \
    "$MAX_CPUS" "$MAX_MEMORY" "$DATA_MOUNT_LINES" <<'PYEOF'
import sys, pathlib
template_root, project_dir, image_version, project_name, max_cpus, max_memory, data_mounts = sys.argv[1:8]
src = pathlib.Path(template_root) / ".devcontainer" / "docker-compose.yml.template"
dst = pathlib.Path(project_dir) / ".devcontainer" / "docker-compose.yml"
content = src.read_text()
content = content.replace("{{IMAGE_VERSION}}", image_version)
content = content.replace("{{PROJECT_NAME}}", project_name)
content = content.replace("{{MAX_CPUS}}", max_cpus)
content = content.replace("{{MAX_MEMORY}}", max_memory)
content = content.replace("{{DATA_MOUNTS}}", data_mounts)
dst.write_text(content)
PYEOF

# Copy devcontainer scripts from template
if [ -d "${TEMPLATE_PATH}/.devcontainer/scripts" ]; then
    cp -r "${TEMPLATE_PATH}/.devcontainer/scripts"/. "${PROJECT_DIR}/.devcontainer/scripts/" 2>/dev/null || true
    chmod +x "${PROJECT_DIR}/.devcontainer/scripts"/*.sh 2>/dev/null || true
fi

# Copy poststart sanity script (fallback to repo-level scripts/ if not in template)
if [ ! -f "${PROJECT_DIR}/.devcontainer/scripts/poststart_sanity.sh" ]; then
    if [ -f "${REPO_ROOT}/scripts/poststart_sanity.sh" ]; then
        cp "${REPO_ROOT}/scripts/poststart_sanity.sh" \
           "${PROJECT_DIR}/.devcontainer/scripts/poststart_sanity.sh"
        chmod +x "${PROJECT_DIR}/.devcontainer/scripts/poststart_sanity.sh"
    elif [ -f "${REPO_ROOT}/.devcontainer/scripts/poststart_sanity.sh" ]; then
        cp "${REPO_ROOT}/.devcontainer/scripts/poststart_sanity.sh" \
           "${PROJECT_DIR}/.devcontainer/scripts/poststart_sanity.sh"
        chmod +x "${PROJECT_DIR}/.devcontainer/scripts/poststart_sanity.sh"
    else
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
fi

# Copy any documentation that lives next to the template's devcontainer
for doc in MCP_AUTH_SETUP.md PYTHON_VENV_GUIDE.md; do
    if [ -f "${TEMPLATE_PATH}/.devcontainer/${doc}" ]; then
        cp "${TEMPLATE_PATH}/.devcontainer/${doc}" "${PROJECT_DIR}/.devcontainer/"
    fi
done

# Create/update .env file in .devcontainer/ (preserve existing API keys if any)
ENV_FILE="${PROJECT_DIR}/.devcontainer/.env"
if [ -f "$ENV_FILE" ]; then
    echo "Updating existing .env file (preserving API keys)..."
    CONTEXT7_KEY=$(grep "^CONTEXT7_API_KEY=" "$ENV_FILE" 2>/dev/null | cut -d= -f2- || echo "")
    GEMINI_KEY=$(grep "^GEMINI_API_KEY=" "$ENV_FILE" 2>/dev/null | cut -d= -f2- || echo "")
    OPENAI_KEY=$(grep "^OPENAI_API_KEY=" "$ENV_FILE" 2>/dev/null | cut -d= -f2- || echo "")
else
    echo "Creating .env file..."
    CONTEXT7_KEY=""
    GEMINI_KEY=""
    OPENAI_KEY=""
fi

cat > "$ENV_FILE" <<EOF
# Docker Compose environment variables
LOCAL_UID=$(id -u)
LOCAL_GID=$(id -g)
WORKSPACE_FOLDER=..

# MCP Server API Keys (optional - only if using Claude Code)

# Context7 - Up-to-date library docs (optional - works without key)
# Get key for higher rate limits: https://context7.com/dashboard
CONTEXT7_API_KEY=${CONTEXT7_KEY}

# PAL MCP Server - Multi-model AI collaboration
# Requires at least one key to function. Get keys from:
#   - Gemini: https://aistudio.google.com/apikey
#   - OpenAI: https://platform.openai.com/api-keys
GEMINI_API_KEY=${GEMINI_KEY}
OPENAI_API_KEY=${OPENAI_KEY}

# Resource Limits (adjust based on your system)
MAX_CPUS=${MAX_CPUS}
MAX_MEMORY=${MAX_MEMORY}
EOF

# Create .gitkeep files (the template already provides these via copy patterns,
# but ensure presence for any directories we created above that the user didn't ship)
for keep in 00_data/raw 00_data/processed 00_data/references \
            01_modules/.ref \
            02_analysis/scripts 02_analysis/notebooks 02_analysis/helpers \
            03_results/checkpoints 03_results/plots 03_results/tables \
            docs/raw docs/ai-generated/vignettes docs/ai-generated/research \
            docs/plan/phase-1 \
            logs; do
    [ -e "${PROJECT_DIR}/${keep}/.gitkeep" ] || touch "${PROJECT_DIR}/${keep}/.gitkeep"
done

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
Image: scdock-r-dev:${IMAGE_VERSION}
Date: $(date +%Y-%m-%d)
"
        echo -e "${GREEN}Git repository initialized${NC}"
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

    add_submodule_with_fallback() {
        local repo_owner="tony-zhelonkin"
        local repo_name="$1"
        local branch="$2"
        local target_path="$3"
        local ssh_url="git@github.com:${repo_owner}/${repo_name}.git"

        if [ -d "$target_path" ]; then
            echo -e "  ${YELLOW}${repo_name} directory already exists, skipping${NC}"
            return 1
        fi

        if git submodule add -b "$branch" "$ssh_url" "$target_path" 2>/dev/null; then
            echo -e "  ${GREEN}${repo_name} added via git (SSH)${NC}"
            return 0
        fi

        echo "  SSH failed, trying gh CLI fallback (HTTPS)..."
        if command -v gh &> /dev/null && gh auth status &> /dev/null; then
            if GIT_CONFIG_COUNT=1 \
               GIT_CONFIG_KEY_0="url.https://github.com/.insteadOf" \
               GIT_CONFIG_VALUE_0="git@github.com:" \
               gh repo clone "${repo_owner}/${repo_name}" "$target_path" -- -b "$branch" 2>/dev/null; then
                git config -f .gitmodules "submodule.${target_path}.path" "$target_path"
                git config -f .gitmodules "submodule.${target_path}.url" "$ssh_url"
                git config -f .gitmodules "submodule.${target_path}.branch" "$branch"
                git config "submodule.${target_path}.url" "$ssh_url"
                git config "submodule.${target_path}.active" "true"
                git add "$target_path"
                echo -e "  ${GREEN}${repo_name} added via gh CLI (HTTPS)${NC}"
                return 0
            fi
        fi

        echo -e "  ${RED}Failed to add ${repo_name}${NC}"
        return 1
    }

    echo "  Adding RNAseq-toolkit (branch: ${RNASEQ_TOOLKIT_BRANCH})..."
    if add_submodule_with_fallback "RNAseq-toolkit" "${RNASEQ_TOOLKIT_BRANCH}" "01_modules/RNAseq-toolkit"; then
        SUBMODULES_ADDED=true
    fi

    echo "  Adding SciAgent-toolkit (branch: ${SCIAGENT_TOOLKIT_BRANCH})..."
    if add_submodule_with_fallback "SciAgent-toolkit" "${SCIAGENT_TOOLKIT_BRANCH}" "01_modules/SciAgent-toolkit"; then
        SUBMODULES_ADDED=true
    fi

    if [ "$SUBMODULES_ADDED" = true ]; then
        if [ -f ".gitmodules" ]; then
            git add .gitmodules 01_modules/
            git commit -m "Add analysis toolkits as git submodules

- RNAseq-toolkit (${RNASEQ_TOOLKIT_BRANCH} branch): Reusable RNA-seq analysis functions
- SciAgent-toolkit: AI infrastructure and MCP server setup
"
            echo -e "${GREEN}Submodules committed${NC}"
        fi
    fi

    cd - > /dev/null
fi

echo ""
echo -e "${GREEN}Project initialized successfully${NC}"
echo ""
echo -e "${BLUE}Project Summary:${NC}"
echo "  Name: ${PROJECT_NAME}"
echo "  Location: ${PROJECT_DIR}"
echo "  Template: ${TEMPLATE}"
echo "  Image: scdock-r-dev:${IMAGE_VERSION}"
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
    echo "  4. Reopen in container: Ctrl+Shift+P -> 'Dev Containers: Reopen in Container'"
    echo "  5. Run AI setup: ./01_modules/SciAgent-toolkit/scripts/setup-ai.sh"
    echo "  6. Fill in context.md with your scientific question"
    echo "  7. Edit 02_analysis/config/pipeline.yaml for your experiment"
else
    echo "  2. Open in VS Code: code ${PROJECT_DIR}"
    echo "  3. Reopen in container: Ctrl+Shift+P -> 'Dev Containers: Reopen in Container'"
    echo "  4. Run AI setup: ./01_modules/SciAgent-toolkit/scripts/setup-ai.sh"
    echo "  5. Fill in context.md with your scientific question"
    echo "  6. Edit 02_analysis/config/pipeline.yaml for your experiment"
fi
echo ""
echo -e "${BLUE}Documentation:${NC}"
echo "  - Project README: ${PROJECT_DIR}/README.md"
echo "  - Docs hub: ${PROJECT_DIR}/docs/README.md"
echo "  - Research plan: ${PROJECT_DIR}/docs/plan/README.md"
echo "  - Config: ${PROJECT_DIR}/02_analysis/config/"
echo ""
echo -e "${BLUE}scbio-docker references:${NC}"
echo "  - Repo root: ${REPO_ROOT}"
echo "  - Architecture: ${REPO_ROOT}/CLAUDE.md"
echo ""
