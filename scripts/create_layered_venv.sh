#!/usr/bin/env bash
set -euo pipefail

# create_layered_venv.sh - Create Python venv with --system-site-packages
# Inherits from base venv, only installs additional packages

ISOLATED=0
if [ "${1:-}" = "--isolated" ]; then
    ISOLATED=1
    shift
fi

if [ $# -lt 2 ]; then
    echo "Usage: create_layered_venv.sh [--isolated] <venv_name> <requirements_file>"
    echo ""
    echo "Examples:"
    echo "  create_layered_venv.sh squid /opt/environments/squid.txt"
    echo "  create_layered_venv.sh atac /opt/environments/atac.txt"
    echo "  create_layered_venv.sh comms /opt/environments/comms.txt"
    echo "  create_layered_venv.sh --isolated scenic /opt/environments/scenic.txt"
    echo ""
    echo "Creates venv at /opt/venvs/<venv_name> inheriting from base venv."
    echo "--isolated drops --system-site-packages (needed when a stack pins a"
    echo "package differently from base, e.g. scenic's pandas 1.5)."
    exit 1
fi

ENV_NAME=$1
REQ_FILE=$2
VENV_DIR="/opt/venvs/${ENV_NAME}"

# Check if requirements file exists
if [ ! -f "$REQ_FILE" ]; then
    echo "Error: Requirements file not found: $REQ_FILE" >&2
    exit 1
fi

# Layered (inherits base) unless --isolated
if [ "$ISOLATED" -eq 1 ]; then
    echo "Creating ISOLATED venv: $ENV_NAME (no base inheritance)"
    python3.11 -m venv "$VENV_DIR"
else
    echo "Creating layered venv: $ENV_NAME (inherits from base)"
    python3.11 -m venv --system-site-packages "$VENV_DIR"
fi

# Upgrade pip/setuptools/wheel
"$VENV_DIR/bin/python" -m pip install --upgrade --no-cache-dir pip setuptools wheel

# Install only additional packages
echo "Installing packages from $REQ_FILE..."
"$VENV_DIR/bin/pip" install --no-cache-dir -r "$REQ_FILE"

# Freeze for reproducibility (includes inherited packages)
echo "Freezing environment..."
"$VENV_DIR/bin/pip" freeze > "/opt/environments/${ENV_NAME}_frozen.txt"

echo "✓ Layered venv created: $VENV_DIR"
echo "  Activate with: source $VENV_DIR/bin/activate"
echo "  Or use: usepy $ENV_NAME"
echo "  Frozen requirements: /opt/environments/${ENV_NAME}_frozen.txt"
