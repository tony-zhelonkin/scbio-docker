#!/usr/bin/env bash
set -euo pipefail

# create_layered_venv.sh - Create Python venv with --system-site-packages
# Inherits from base venv, only installs additional packages

if [ $# -lt 2 ]; then
    echo "Usage: create_layered_venv.sh <venv_name> <requirements_file>"
    echo ""
    echo "Examples:"
    echo "  create_layered_venv.sh squid /opt/environments/squid.txt"
    echo "  create_layered_venv.sh atac /opt/environments/atac.txt"
    echo "  create_layered_venv.sh comms /opt/environments/comms.txt"
    echo ""
    echo "Creates venv at /opt/venvs/<venv_name> inheriting from base venv"
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

# Create layered venv (inherits base packages)
echo "Creating layered venv: $ENV_NAME (inherits from base)"
python3 -m venv --system-site-packages "$VENV_DIR"

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
