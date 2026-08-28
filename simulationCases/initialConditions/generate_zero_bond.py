#!/usr/bin/env python3
"""Bo = 0 bursting-bubble initial condition: sphere intersecting a plane.

Usage:
    ./generate_zero_bond.py [--delta 0.01] [--rmax 32] [--out ../DataFiles/Bo0.0000.dat]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_DIR = Path(__file__).resolve().parent
if str(_DIR) not in sys.path:
    sys.path.insert(0, str(_DIR))

from zero_bond import sphere_plane, write_basilisk_dat  # noqa: E402


def main(argv=None) -> int:
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--delta",
        type=float,
        default=0.01,
        help="fillet/neck regularisation (initial hole radius ~ 2 delta)",
    )
    ap.add_argument(
        "--rmax",
        type=float,
        default=32.0,
        help="radial extent of the flat free surface (>= domain size)",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=here.parent / "DataFiles" / "Bo0.0000.dat",
    )
    args = ap.parse_args(argv)

    shape = sphere_plane(delta=args.delta, rmax=args.rmax)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    write_basilisk_dat(shape, args.out)
    print(
        f"wrote {args.out}  (delta={args.delta:g}, Rf={shape.fillet_radius:.5f}, "
        f"{len(shape.axial)} points)"
    )
    print(
        f"  x_axial [{shape.axial.min():.4f}, {shape.axial.max():.4f}]  "
        f"(south pole {shape.south_pole_axial():.4f}, surface 0)"
    )
    print(f"  y_radial[{shape.radial.min():.4f}, {shape.radial.max():.4f}]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
