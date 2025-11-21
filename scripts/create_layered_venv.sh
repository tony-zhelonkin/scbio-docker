#!/usr/bin/env bash
set -euo pipefail

# create_layered_venv.sh NAME REQUIREMENTS_TXT
# Creates a venv at /opt/venvs/NAME with --system-site-packages and installs given requirements.

if [ $# -lt 2 ]; then
  echo "usage: $0 <name> <requirements.txt>" >&2
  exit 1
fi

NAME="$1"
REQ="$2"

python3 -m venv --system-site-packages "/opt/venvs/${NAME}"
"/opt/venvs/${NAME}/bin/python" -m pip install --upgrade pip setuptools wheel
if [ -f "$REQ" ]; then
  "/opt/venvs/${NAME}/bin/pip" install --no-cache-dir -r "$REQ"
else
  echo "requirements file not found: $REQ" >&2
fi

