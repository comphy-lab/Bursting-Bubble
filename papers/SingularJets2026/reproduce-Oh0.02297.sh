#!/bin/bash
# Reproduce the Bo=0, Oh=0.02297 grid-convergence figure (L13 / L14).
# Run from the folder root:  bash reproduce-Oh0.02297.sh
# Requires python3 with numpy + matplotlib (+ LaTeX for the labels; if absent,
# set matplotlib.rcParams["text.usetex"]=False near the top of the plot script).
set -euo pipefail
cd "$(dirname "$0")"
D=data-Oh-0.02297

python3 plotJetMetricsTheory.py \
  --series 0.02297 13 $D/3010_L13_log.txt \
  --series 0.02297 14 $D/4012_L14_log.txt \
  --facet  0.02297 $D/facet_inception.txt \
  --bond 0 --rmin 0.003 --fit-window 0.008 0.025 \
  --grid-name 13 "L13 (3010)" \
  --grid-name 14 "L14 (4012)" \
  --out $D/gridconv_Oh0.02297_L13-L14
echo "wrote $D/gridconv_Oh0.02297_L13-L14.{png,pdf}"
