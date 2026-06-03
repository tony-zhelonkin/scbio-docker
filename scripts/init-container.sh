#!/usr/bin/env bash
# init-container.sh - Render a VS Code dev container (devcontainer.json +
# docker-compose.yml + .env) into a target project directory.
#
# This is container substrate only: it knows nothing about project structure,
# analysis trees, config schemas, docs, or AI context. For the project scaffold
# (directory tree, config, docs, AI harness) use SciAgent-toolkit:
#   sciagent new project --type analysis|software-tool <dir>
#
# Usage:
#   init-container.sh <project-dir> [OPTIONS]
#
# Options:
#   --data-mount KEY:PATH[:ro]    Add a data mount (repeatable)
#   --image-version vX.Y.Z        Image tag (default: read from VERSION)
#   --service dev-core|dev-archr  Compose service (default: dev-core)
#   --max-cpus N                  CPU limit default (default: 50)
#   --max-memory NG               Memory limit default (default: 450G)

set -euo pipefail

# --- Resolve script + repo paths (follow symlinks) ---------------------------
SCRIPT_PATH="${BASH_SOURCE[0]}"
while [ -L "$SCRIPT_PATH" ]; do
    SCRIPT_DIR_TMP="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"
    SCRIPT_PATH="$(readlink "$SCRIPT_PATH")"
    [[ $SCRIPT_PATH != /* ]] && SCRIPT_PATH="$SCRIPT_DIR_TMP/$SCRIPT_PATH"
done
SCRIPT_DIR="$(cd "$(dirname "$SCRIPT_PATH")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
TEMPLATES_DIR="${REPO_ROOT}/templates/devcontainer"

# --- Defaults ----------------------------------------------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
BLUE='\033[0;34m'
NC='\033[0m'

IMAGE_VERSION="$(tr -d '[:space:]' < "$REPO_ROOT/VERSION")"
SERVICE="dev-core"
MAX_CPUS="50"
MAX_MEMORY="450g"
declare -a DATA_MOUNTS=()

usage() {
    cat <<EOF
Usage: $0 <project-dir> [OPTIONS]

Renders .devcontainer/{devcontainer.json,docker-compose.yml,.env} into <project-dir>.

Options:
  --data-mount KEY:PATH[:ro]    Add a data mount (repeatable)
                                KEY is a label, PATH is a host path, :ro for read-only
  --image-version vX.Y.Z        Image tag (default: VERSION file -> ${IMAGE_VERSION})
  --service dev-core|dev-archr  Compose service (default: dev-core)
  --max-cpus N                  CPU limit default (default: 50)
  --max-memory NG               Memory limit default (default: 450G)

Example:
  $0 ~/projects/atac-study \\
      --data-mount atac:/scratch/data/DT-1234 \\
      --data-mount rna:/scratch/data/DT-5678:ro

For the project scaffold (tree, config, docs, AI harness):
  sciagent new project --type analysis <project-dir>
EOF
    exit 1
}

# --- Argument parsing --------------------------------------------------------
[ $# -lt 1 ] && usage
PROJECT_DIR="$1"
shift

while [[ $# -gt 0 ]]; do
    case $1 in
        --data-mount)    DATA_MOUNTS+=("$2"); shift 2 ;;
        --image-version) IMAGE_VERSION="$2"; shift 2 ;;
        --service)       SERVICE="$2"; shift 2 ;;
        --max-cpus)      MAX_CPUS="$2"; shift 2 ;;
        --max-memory)    MAX_MEMORY="$2"; shift 2 ;;
        *)
            echo -e "${RED}Error: Unknown option '$1'${NC}" >&2
            usage
            ;;
    esac
done

case "$SERVICE" in
    dev-core|dev-archr) ;;
    *)
        echo -e "${RED}Error: --service must be dev-core or dev-archr${NC}" >&2
        exit 1
        ;;
esac

if [ ! -d "$TEMPLATES_DIR" ]; then
    echo -e "${RED}Error: template dir '${TEMPLATES_DIR}' not found${NC}" >&2
    exit 1
fi

PROJECT_NAME="$(basename "$PROJECT_DIR")"
mkdir -p "$PROJECT_DIR/.devcontainer/scripts"

# --- Render devcontainer.json (simple single-line token substitution) --------
render_devcontainer_json() {
    sed -e "s|{{PROJECT_NAME}}|${PROJECT_NAME}|g" \
        -e "s|{{SERVICE}}|${SERVICE}|g" \
        "${TEMPLATES_DIR}/.devcontainer/devcontainer.json.template" \
        > "${PROJECT_DIR}/.devcontainer/devcontainer.json"
}

# --- Build the multi-line data-mount block for the compose YAML --------------
build_data_mount_block() {
    local lines=""
    if [ ${#DATA_MOUNTS[@]} -gt 0 ]; then
        lines+="      # Data mounts"$'\n'
        local mount label path ro
        for mount in "${DATA_MOUNTS[@]}"; do
            IFS=':' read -r label path ro <<< "$mount"
            if [ -n "${ro:-}" ]; then
                lines+="      - ${path}:/workspaces/${PROJECT_NAME}/00_data/${label}:ro"$'\n'
            else
                lines+="      - ${path}:/workspaces/${PROJECT_NAME}/00_data/${label}"$'\n'
            fi
        done
    else
        lines+="      # Add your data mounts here:"$'\n'
        lines+="      # - /path/to/data:/workspaces/${PROJECT_NAME}/00_data/raw:ro"$'\n'
    fi
    printf '%s' "$lines"
}

# --- Render docker-compose.yml (Python: {{DATA_MOUNTS}} is multi-line) --------
render_docker_compose() {
    local data_mount_block="$1"
    python3 - "$TEMPLATES_DIR" "$PROJECT_DIR" "$IMAGE_VERSION" "$PROJECT_NAME" \
        "$MAX_CPUS" "$MAX_MEMORY" "$data_mount_block" <<'PYEOF'
import sys, pathlib
tmpl, project_dir, image_version, project_name, max_cpus, max_memory, data_mounts = sys.argv[1:8]
src = pathlib.Path(tmpl) / ".devcontainer" / "docker-compose.yml.template"
dst = pathlib.Path(project_dir) / ".devcontainer" / "docker-compose.yml"
content = src.read_text()
content = content.replace("{{IMAGE_VERSION}}", image_version)
content = content.replace("{{PROJECT_NAME}}", project_name)
content = content.replace("{{MAX_CPUS}}", max_cpus)
content = content.replace("{{MAX_MEMORY}}", max_memory)
content = content.replace("{{DATA_MOUNTS}}", data_mounts)
dst.write_text(content)
PYEOF
}

# --- Copy devcontainer scripts + poststart sanity fallback -------------------
copy_devcontainer_scripts() {
    if [ -d "${TEMPLATES_DIR}/.devcontainer/scripts" ]; then
        cp -r "${TEMPLATES_DIR}/.devcontainer/scripts"/. \
              "${PROJECT_DIR}/.devcontainer/scripts/" 2>/dev/null || true
        chmod +x "${PROJECT_DIR}/.devcontainer/scripts"/*.sh 2>/dev/null || true
    fi

    local dst="${PROJECT_DIR}/.devcontainer/scripts/poststart_sanity.sh"
    [ -f "$dst" ] && return 0

    if [ -f "${REPO_ROOT}/scripts/poststart_sanity.sh" ]; then
        cp "${REPO_ROOT}/scripts/poststart_sanity.sh" "$dst"
    elif [ -f "${REPO_ROOT}/.devcontainer/scripts/poststart_sanity.sh" ]; then
        cp "${REPO_ROOT}/.devcontainer/scripts/poststart_sanity.sh" "$dst"
    else
        cat > "$dst" <<'SANITY_EOF'
#!/bin/bash
echo "=== Container Environment Check ==="
echo "R version: $(R --version | head -n1)"
echo "Python version: $(python3 --version)"
echo "Working directory: $(pwd)"
echo "User: $(whoami) (UID=$(id -u), GID=$(id -g))"
echo "==================================="
SANITY_EOF
    fi
    chmod +x "$dst"
}

# --- Create/update .env at project root (preserve existing MCP keys) ---------
# .env lives at the project root so Docker Compose v2 picks it up for variable
# interpolation (compose v2 reads .env from the working directory, not the
# compose file's directory).
write_env_file() {
    local env_file="${PROJECT_DIR}/.env"
    local context7_key="" gemini_key="" openai_key=""
    if [ -f "$env_file" ]; then
        context7_key=$(grep "^CONTEXT7_API_KEY=" "$env_file" 2>/dev/null | cut -d= -f2- || true)
        gemini_key=$(grep "^GEMINI_API_KEY=" "$env_file" 2>/dev/null | cut -d= -f2- || true)
        openai_key=$(grep "^OPENAI_API_KEY=" "$env_file" 2>/dev/null | cut -d= -f2- || true)
    fi

    cat > "$env_file" <<EOF
# Docker Compose environment variables
LOCAL_UID=$(id -u)
LOCAL_GID=$(id -g)
WORKSPACE_FOLDER=.

# Resource limits (override compose defaults)
MAX_CPUS=${MAX_CPUS}
MAX_MEMORY=${MAX_MEMORY}

# MCP Server API keys (consumed by SciAgent-toolkit's setup-ai.sh later)
# Context7 - library docs (works without a key; key raises rate limits)
CONTEXT7_API_KEY=${context7_key}
# PAL - multi-model AI collaboration (needs at least one of the below)
GEMINI_API_KEY=${gemini_key}
OPENAI_API_KEY=${openai_key}
EOF
}

# --- Main --------------------------------------------------------------------
echo -e "${GREEN}Rendering dev container for '${PROJECT_NAME}' (image scdock-r-dev:${IMAGE_VERSION}, service ${SERVICE})...${NC}"

render_devcontainer_json
render_docker_compose "$(build_data_mount_block)"
copy_devcontainer_scripts
write_env_file

echo ""
echo -e "${GREEN}Dev container rendered into ${PROJECT_DIR}/.devcontainer/${NC}"
echo "  - devcontainer.json"
echo "  - docker-compose.yml"
echo "  - scripts/poststart_sanity.sh"
echo -e "${GREEN}Project env file written to ${PROJECT_DIR}/.env${NC}"
echo ""
echo -e "${BLUE}Next steps:${NC}"
echo "  1. Scaffold the project structure (tree, config, docs, AI harness):"
echo "       sciagent new project --type analysis ${PROJECT_DIR}"
echo "  2. Open in VS Code: code ${PROJECT_DIR}"
echo "  3. Reopen in container: Ctrl+Shift+P -> 'Dev Containers: Reopen in Container'"
echo ""
