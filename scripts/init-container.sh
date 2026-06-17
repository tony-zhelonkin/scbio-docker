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
#   --data-mount KEY:PATH[:rw]    Add a data mount (repeatable; read-only by default)
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
GPU=false
declare -a DATA_MOUNTS=()

usage() {
    cat <<EOF
Usage: $0 <project-dir> [OPTIONS]

Renders .devcontainer/{devcontainer.json,docker-compose.yml,.env} into <project-dir>.

Options:
  --data-mount KEY:PATH[:rw]    Add a data mount (repeatable; read-only by default)
                                KEY is a label, PATH is a host path, :rw for read-write
  --image-version vX.Y.Z        Image tag (default: VERSION file -> ${IMAGE_VERSION})
  --service dev-core|dev-archr  Compose service (default: dev-core)
  --max-cpus N                  CPU limit default (default: 50)
  --max-memory NG               Memory limit default (default: 450G)
  --gpu                         Enable NVIDIA GPU passthrough (adds devices block)

Example:
  $0 ~/projects/atac-study \\
      --data-mount atac:/scratch/data/DT-1234 \\
      --data-mount scratch:/scratch/work/DT-5678:rw

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
        --gpu)           GPU=true; shift ;;
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

# --- Build the SSH agent mount line (empty when no agent socket present) -----
# Auto-detects SSH_AUTH_SOCK. If it's a live socket, returns a single volume
# line that forwards it into the container. If absent/not a socket, returns
# empty — the {{SSH_AGENT_MOUNT}} token is replaced with nothing, keeping the
# compose file valid without any fallback-socket hacks.
build_ssh_agent_mount() {
    local sock="${SSH_AUTH_SOCK:-}"
    if [[ -n "$sock" && -S "$sock" ]]; then
        printf '      - %s:/ssh-agent:ro\n' "$sock"
    fi
}

# --- Build the GPU devices block (empty when --gpu not passed) ---------------
build_gpu_devices() {
    if [[ "$GPU" == "true" ]]; then
        printf '          devices:\n            - driver: nvidia\n              count: all\n              capabilities: [gpu]\n'
    fi
}

# --- Build the multi-line data-mount block for the compose YAML --------------
build_data_mount_block() {
    local lines=""
    if [ ${#DATA_MOUNTS[@]} -gt 0 ]; then
        lines+="      # Data mounts"$'\n'
        local mount label path mode
        for mount in "${DATA_MOUNTS[@]}"; do
            IFS=':' read -r label path mode <<< "$mount"
            # Read-only by default (00_data is input data); pass :rw to opt out.
            if [ "${mode:-ro}" = "rw" ]; then
                lines+="      - ${path}:/workspaces/${PROJECT_NAME}/00_data/${label}"$'\n'
            else
                lines+="      - ${path}:/workspaces/${PROJECT_NAME}/00_data/${label}:ro"$'\n'
            fi
        done
    else
        lines+="      # Add your data mounts here:"$'\n'
        lines+="      # - /path/to/data:/workspaces/${PROJECT_NAME}/00_data/raw:ro"$'\n'
    fi
    printf '%s' "$lines"
}

# --- Render docker-compose.yml (Python: multi-line tokens passed as args) ----
render_docker_compose() {
    local data_mount_block="$1" ssh_agent_mount="$2" gpu_devices="$3"
    python3 - "$TEMPLATES_DIR" "$PROJECT_DIR" "$IMAGE_VERSION" "$PROJECT_NAME" \
        "$MAX_CPUS" "$MAX_MEMORY" "$data_mount_block" "$ssh_agent_mount" "$gpu_devices" <<'PYEOF'
import sys, pathlib
tmpl, project_dir, image_version, project_name, max_cpus, max_memory, data_mounts, ssh_agent, gpu_devices = sys.argv[1:10]
src = pathlib.Path(tmpl) / ".devcontainer" / "docker-compose.yml.template"
dst = pathlib.Path(project_dir) / ".devcontainer" / "docker-compose.yml"
content = src.read_text()
content = content.replace("{{IMAGE_VERSION}}", image_version)
content = content.replace("{{PROJECT_NAME}}", project_name)
content = content.replace("{{MAX_CPUS}}", max_cpus)
content = content.replace("{{MAX_MEMORY}}", max_memory)
content = content.replace("{{DATA_MOUNTS}}", data_mounts)
content = content.replace("{{SSH_AGENT_MOUNT}}", ssh_agent)
content = content.replace("{{GPU_DEVICES}}", gpu_devices)
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

# --- Create/update .devcontainer/.env (preserve existing MCP keys) -----------
# .env lives alongside the compose file in .devcontainer/ so Docker Compose
# picks it up for variable interpolation (compose --project-directory defaults
# to the compose file's directory).
write_env_file() {
    local env_file="${PROJECT_DIR}/.devcontainer/.env"
    local gemini_key="" openai_key=""
    if [ -f "$env_file" ]; then
        gemini_key=$(grep "^GEMINI_API_KEY=" "$env_file" 2>/dev/null | cut -d= -f2- || true)
        openai_key=$(grep "^OPENAI_API_KEY=" "$env_file" 2>/dev/null | cut -d= -f2- || true)
    fi

    cat > "$env_file" <<EOF
# Docker Compose environment variables
LOCAL_UID=$(id -u)
LOCAL_GID=$(id -g)
WORKSPACE_FOLDER=..

# Resource limits (override compose defaults)
MAX_CPUS=${MAX_CPUS}
MAX_MEMORY=${MAX_MEMORY}

# Local Ollama endpoint — Docker bridge IP (host as seen from inside container)
# Update if Ollama runs on a different host or port
OLLAMA_HOST=http://172.17.0.1:11434

# MCP Server API keys (consumed by SciAgent-toolkit's setup-ai.sh later)
# PAL - multi-model AI collaboration (needs at least one of the below)
GEMINI_API_KEY=${gemini_key}
OPENAI_API_KEY=${openai_key}
EOF
}

# --- Main --------------------------------------------------------------------
echo -e "${GREEN}Rendering dev container for '${PROJECT_NAME}' (image scdock-r-dev:${IMAGE_VERSION}, service ${SERVICE})...${NC}"

render_devcontainer_json
render_docker_compose "$(build_data_mount_block)" "$(build_ssh_agent_mount)" "$(build_gpu_devices)"
copy_devcontainer_scripts
write_env_file

echo ""
echo -e "${GREEN}Dev container rendered into ${PROJECT_DIR}/.devcontainer/${NC}"
echo "  - devcontainer.json"
echo "  - docker-compose.yml"
echo "  - .env"
echo "  - scripts/poststart_sanity.sh"
echo ""
echo -e "${BLUE}Next steps:${NC}"
echo "  1. Scaffold the project structure (tree, config, docs, AI harness):"
echo "       sciagent new project --type analysis ${PROJECT_DIR}"
echo "  2. Open in VS Code: code ${PROJECT_DIR}"
echo "  3. Reopen in container: Ctrl+Shift+P -> 'Dev Containers: Reopen in Container'"
echo ""
