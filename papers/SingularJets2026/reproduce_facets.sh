#!/bin/bash
# Reproduce the Bo=0, Oh=0.03 self-similar interface-COLLAPSE figure (Figure 3):
# pre/post-inception interface profiles, raw + rescaled by |t-t0|^alpha (our
# exponent, not 2/3). Run from the folder root:  bash reproduce_facets.sh
# Requires python3 with numpy + matplotlib.
set -euo pipefail
cd "$(dirname "$0")"
python3 facet-collapse-figure3/fig3_collapse.py
echo "wrote facet-collapse-figure3/fig3_Oh0.03_collapse.{png,pdf}"
