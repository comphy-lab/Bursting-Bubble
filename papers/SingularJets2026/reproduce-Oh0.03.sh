#!/bin/bash
# Reproduce the Bo=0, Oh=0.03 grid-convergence figure (L13 / L14 / two L15 runs).
# Run from the folder root:  bash reproduce-Oh0.03.sh
# Requires python3 with numpy + matplotlib (+ LaTeX for the labels; if absent,
# set matplotlib.rcParams["text.usetex"]=False near the top of the plot script).
set -euo pipefail
cd "$(dirname "$0")"
D=data-Oh-0.03

python3 plotJetMetricsTheory.py \
  --series 0.03 13 $D/3013_L13_log.txt \
  --series 0.03 14 $D/4015_L14_log.txt \
  --series 0.03 15 $D/5001_L15_focus13_log.txt \
  --series 0.03 16 $D/5003_L15_focus14_log.txt \
  --series 0.03 17 $D/5008_L15_focus15_log.txt \
  --facet  0.03 $D/facet_inception.txt \
  --bond 0 --rmin 0.003 --fit-window 0.008 0.025 \
  --grid-name 13 "L13 (3013)" \
  --grid-name 14 "L14 (4015)" \
  --grid-name 15 "L15 (5001, focus 13)" \
  --grid-name 16 "L15 (5003, focus 14)" \
  --grid-name 17 "L15 (5008, focus 15)" \
  --out $D/gridconv3_Oh0.03_L13-L14-L15
# NOTE: focus-15 series is now 5008 (survived inception with keStopMax=1e6; the
# earlier 5004 tripped the default 1e4 ke gate at inception, reaching only
# r_jet~0.015). 5008 reaches r_jet~0.05; re-run this once a to-completion
# focus-15 run extends further.
echo "wrote $D/gridconv3_Oh0.03_L13-L14-L15.{png,pdf}"
