#!/usr/bin/env python3
"""Generate a Bond-number-dependent bursting-bubble initial polyline.

Writes ``Bo%5.4f.dat`` in the Basilisk convention used by Stage 1
(``distance.h`` / ``input_xy``): axial column, radial column, cavity in
−axial, far-field free surface at axial = 0.

Examples
--------
    ./generate_bond_shape.py --bond 0.001
    ./generate_bond_shape.py --bond 0.01,0.04 --out-dir ../DataFiles
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

_DIR = Path(__file__).resolve().parent
if str(_DIR) not in sys.path:
    sys.path.insert(0, str(_DIR))

from young_laplace import (  # noqa: E402
    EquilibriumShape,
    bond_filename,
    continuation_ladder,
    solve_equilibrium,
    write_basilisk_dat,
)


def _parse_bonds(text: str) -> list[float]:
    bonds = []
    for token in text.split(","):
        token = token.strip()
        if token:
            bonds.append(float(token))
    if not bonds:
        raise argparse.ArgumentTypeError("need at least one Bond number")
    return bonds


def _print_shape(shape: EquilibriumShape, path: Path) -> None:
    print(
        f"Bo={shape.bond:g}  Rb={shape.Rb:.6f}  Rc={shape.Rc:.6f}  "
        f"phic={shape.phic:.6f}  alpha_c={shape.alpha_c:.6f}  "
        f"hinf={shape.hinf:.6f}  vol_err={shape.volume_err:.3e}  "
        f"tail_err={shape.tail_err:.3e}"
    )
    print(
        f"  wrote {path}  ({len(shape.axial)} points, "
        f"south pole {shape.south_pole_axial():.5f}, "
        f"R in [{shape.radial.min():.3e}, {shape.radial.max():.3f}])"
    )


def main(argv=None) -> int:
    here = Path(__file__).resolve().parent
    default_out = here.parent / "DataFiles"

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--bond",
        type=_parse_bonds,
        required=True,
        help="Bond number, or a comma-separated list (R0-based)",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        default=default_out,
        help=f"directory for BoXXXX.dat (default: {default_out})",
    )
    ap.add_argument(
        "--rmax",
        type=float,
        default=32.0,
        help="radial extent of the written free surface (default 32)",
    )
    ap.add_argument(
        "--fillet-span",
        type=float,
        default=0.22,
        help="physical fillet length along the meniscus (default 0.22)",
    )
    ap.add_argument(
        "--no-continue",
        action="store_true",
        help="disable Bond-parameter continuation (cold start at each value)",
    )
    args = ap.parse_args(argv)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    previous = None
    for bond in args.bond:
        if bond == 0.0:
            print(
                "Bo = 0 is the sphere-plane generator: use generate_zero_bond.py",
                file=sys.stderr,
            )
            return 2
        if not args.no_continue:
            seed = 1e-3 if previous is None else max(1e-3, previous.bond)
            ladder = continuation_ladder(bond, seed=min(seed, bond))
            if len(ladder) > 1:
                print(
                    f"Bo={bond:g}: continuing through "
                    + ", ".join(f"{b:g}" for b in ladder)
                )
        shape = solve_equilibrium(
            bond,
            previous=previous,
            continue_in_bond=not args.no_continue,
            rmax_out=args.rmax,
            fillet_span=args.fillet_span,
        )
        path = args.out_dir / bond_filename(bond)
        write_basilisk_dat(shape, path)
        _print_shape(shape, path)
        previous = shape
    return 0


if __name__ == "__main__":
    sys.exit(main())
