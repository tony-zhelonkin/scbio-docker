#!/usr/bin/env bash
set -euo pipefail

ok()  { echo "OK"; }
fail(){ echo "NOT OK"; }

sep(){ echo "==== $* ===="; }

sep "Devcontainer sanity checks"

printf "Python: "
if python -c "import sys; print(sys.version.split()[0])" >/dev/null 2>&1; then ok; else fail; fi

printf "R:      "
if R -q --vanilla -e "invisible(TRUE)" >/dev/null 2>&1; then ok; else fail; fi

printf "httpgd: "
if R -q --vanilla -e "suppressMessages(library(httpgd)); cat(as.character(packageVersion('httpgd')), '\n')" >/dev/null 2>&1; then ok; else fail; fi

printf "scanpy import: "
if python -c "import scanpy" >/dev/null 2>&1; then ok; else fail; fi

printf "default venv: "
if which python | grep -q "/opt/venvs/base/bin/python"; then echo "/opt/venvs/base" OK; else fail; fi

sep "Tips"
echo "usepy base|squid|atac   # switch Python env"
echo "r-base / r-archr       # regular vs ArchR R sessions"
