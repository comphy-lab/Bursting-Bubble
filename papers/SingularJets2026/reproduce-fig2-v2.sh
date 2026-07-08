#!/bin/bash
# Reproduce the manuscript Fig. 2 v2 panel.
set -euo pipefail
cd "$(dirname "$0")"

if [ -x ../../.venv/bin/python ]; then
  PY=(../../.venv/bin/python)
elif command -v uv >/dev/null 2>&1; then
  PY=(uv run --with numpy --with matplotlib --with scipy python)
else
  PY=(python3)
fi

"${PY[@]}" figure-scripts/make_fig2_v2.py --output fig2_v2.pdf

if command -v pdftoppm >/dev/null 2>&1; then
  pdftoppm -png -singlefile -r 300 fig2_v2.pdf fig2_v2
fi

echo "wrote fig2_v2.pdf"
