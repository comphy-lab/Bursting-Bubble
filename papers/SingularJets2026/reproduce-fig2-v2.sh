#!/bin/bash
# Reproduce the manuscript Fig. 2 v2 panel.
set -euo pipefail
cd "$(dirname "$0")"

if command -v uv >/dev/null 2>&1; then
  PY=(uv run --frozen --offline --project . python)
elif [ -x .venv/bin/python ]; then
  PY=(.venv/bin/python)
else
  PY=(python3)
fi

if [ "${SINGULARJETS_NO_TEX:-0}" = "1" ]; then
  "${PY[@]}" figure-scripts/make_fig2a_streamlines.py \
    --output fig2a_streamlines.pdf --no-tex
  "${PY[@]}" figure-scripts/make_fig2_v2.py \
    --output fig2_v2.pdf --no-tex
else
  "${PY[@]}" figure-scripts/make_fig2a_streamlines.py \
    --output fig2a_streamlines.pdf
  "${PY[@]}" figure-scripts/make_fig2_v2.py \
    --output fig2_v2.pdf
fi

for output in \
  fig2a_streamlines.pdf fig2a_streamlines.png \
  fig2_v2.pdf fig2_v2.png; do
  if [ ! -s "$output" ]; then
    echo "Missing or empty figure output: $output" >&2
    exit 1
  fi
done

echo "wrote fig2a_streamlines.{pdf,png} and fig2_v2.{pdf,png}"
