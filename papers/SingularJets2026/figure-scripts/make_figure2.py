#!/usr/bin/env python3
"""Write the exact vector Figure 2(a-c) render used by the paper.

The supplied PDF includes the four streamlines, velocity-magnitude, and
interface snapshots in panel (a), the jet-base data in panel (b), and the
Cattaneo comparison in panel (c). It is a vector render, not a new calculation.
The separate comparison notebook and source tables are attributed to Javier
Rodríguez-Rodríguez (Javi).
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "output" / "figure2.pdf")
    args = parser.parse_args()
    source = ROOT / "assets" / "main" / "figure2.pdf"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
