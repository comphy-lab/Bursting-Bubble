#!/usr/bin/env python3
"""Reproduce the Lhuissier & Villermaux (2012) opening-angle comparison.

The abscissa is (Rc/R0) sqrt(Bo) = Rc/a. The ordinate is 2 α_c / π,
with α_c = π − φ_c. The small-Bo line is α_c = Rc / (2 √3 a).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def _configure_matplotlib(*, usetex: bool) -> None:
    matplotlib.rcParams["font.family"] = "serif"
    matplotlib.rcParams["font.serif"] = ["Computer Modern Roman"]
    matplotlib.rcParams["text.usetex"] = bool(usetex)
    if usetex:
        matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    else:
        matplotlib.rcParams["mathtext.fontset"] = "cm"


_DIR = Path(__file__).resolve().parent
if str(_DIR) not in sys.path:
    sys.path.insert(0, str(_DIR))

from young_laplace import EquilibriumShape, solve_equilibrium  # noqa: E402

VILLERMAUX_CSV = _DIR / "reference" / "Villermaux.csv"

# Continuation ladder covering the historical comparison set.
PLOT_BONDS = [
    1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 0.04, 0.10, 0.14, 0.25, 0.55,
    0.75, 1.00, 1.25, 2.00, 2.50, 3.50, 5.00, 8.88, 10.0, 11.25,
    15.0, 20.0, 25.0, 36.0, 50.0, 70.0, 100.0, 175.0, 222.0,
    325.0, 500.0, 900.0,
]
INSET_BONDS = [0.01, 1.0, 5.0, 222.0]


def _load_villermaux(path: Path) -> tuple[np.ndarray, np.ndarray]:
    data = np.loadtxt(path, delimiter=",")
    return data[:, 0], data[:, 1]


def _sweep(bonds, *, skip_failed: bool = False) -> list[EquilibriumShape]:
    shapes = []
    previous = None
    for bond in bonds:
        try:
            shape = solve_equilibrium(float(bond), previous=previous)
        except Exception as exc:
            if not skip_failed:
                raise
            print(f"Bo={bond:g} skipped: {exc}", file=sys.stderr)
            continue
        shapes.append(shape)
        previous = shape
        print(
            f"Bo={shape.bond:g}  2α_c/π={shape.opening_metric:.4f}  "
            f"Rc/a={shape.capillary_metric:.4f}  Rb={shape.Rb:.4f}"
        )
    if not shapes:
        raise RuntimeError("no Bond numbers produced a shape")
    return shapes


def _nearest(shapes: list[EquilibriumShape], bond: float) -> EquilibriumShape | None:
    if not shapes:
        return None
    return min(shapes, key=lambda s: abs(np.log(s.bond / bond)))


def _draw_shape(ax, shape: EquilibriumShape) -> None:
    """Physics coordinates with Z up from the south pole, mirrored."""
    segs = [
        (shape.R_bubble, shape.Z_bubble),
        (shape.R_cap, shape.Z_cap),
        (shape.R_tail, shape.Z_tail),
    ]
    for R, Z in segs:
        ax.plot(R, Z, "k-", lw=1.6)
        ax.plot(-R, Z, "k-", lw=1.6)
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(0.0, 2.15)
    ax.set_aspect("equal")
    ax.axis("off")


def plot_opening_angle(
    shapes: list[EquilibriumShape],
    out: Path,
    *,
    usetex: bool = True,
) -> None:
    xv, yv = _load_villermaux(VILLERMAUX_CSV)
    x = np.array([s.capillary_metric for s in shapes])
    y = np.array([s.opening_metric for s in shapes])

    xth = np.logspace(-2.2, 1.6, 400)
    yth = xth / (np.pi * np.sqrt(3.0))

    fig, ax = plt.subplots(figsize=(12, 12))
    ax.plot(
        xv, yv,
        linestyle="None",
        marker="o",
        markersize=14,
        markerfacecolor="k",
        markeredgecolor="k",
        label=(
            r"Lhuissier \& Villermaux (2012)"
            if usetex
            else "Lhuissier & Villermaux (2012)"
        ),
        zorder=3,
    )
    ax.plot(
        x, y,
        linestyle="None",
        marker="*",
        markersize=18,
        markerfacecolor="#1A64B3",
        markeredgecolor="k",
        markeredgewidth=0.6,
        label="Present results",
        zorder=4,
    )
    ax.plot(xth, yth, "--", color="#d62728", lw=3.0, zorder=2)
    ax.plot([1e-2, 30.0], [1.0, 1.0], "--", color="#d62728", lw=2.4, zorder=1)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(1e-2, 30.0)
    ax.set_ylim(1e-3, 2.0)
    ax.set_xlabel(
        r"$\sqrt{\mathcal{B}o\,(R_c/R_0)^2}$",
        fontsize=40,
        labelpad=15,
    )
    ax.set_ylabel(r"$2\alpha_c/\pi$", fontsize=40, labelpad=15)
    ax.tick_params(which="both", direction="out", width=3, labelsize=30, pad=10)
    ax.tick_params(which="major", length=12)
    ax.tick_params(which="minor", length=6)
    for spine in ax.spines.values():
        spine.set_linewidth(3)
    ax.minorticks_on()
    ax.set_box_aspect(1)
    ax.legend(loc="upper left", frameon=False, fontsize=26)
    ax.text(
        0.22, 0.42,
        r"$\alpha_c=\frac{1}{2\sqrt{3}}\,\frac{R_c}{a}$",
        transform=ax.transAxes,
        color="#d62728",
        fontsize=28,
    )

    inset_pos = [
        [0.20, 0.16, 0.22, 0.13],
        [0.34, 0.36, 0.22, 0.13],
        [0.46, 0.52, 0.22, 0.13],
        [0.58, 0.66, 0.22, 0.13],
    ]
    for bond, pos in zip(INSET_BONDS, inset_pos):
        shape = _nearest(shapes, bond)
        if shape is None or abs(np.log(shape.bond / bond)) > 0.8:
            continue
        ax_in = fig.add_axes(pos)
        _draw_shape(ax_in, shape)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", pad_inches=0.1, dpi=300)
    png = out.with_suffix(".png")
    fig.savefig(png, bbox_inches="tight", pad_inches=0.1, dpi=200)
    plt.close(fig)
    print(f"wrote {out}")
    print(f"wrote {png}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--out",
        type=Path,
        default=_DIR / "opening_angle.pdf",
    )
    ap.add_argument(
        "--no-usetex",
        action="store_true",
        help="render math with matplotlib mathtext (no external LaTeX)",
    )
    ap.add_argument(
        "--skip-failed",
        action="store_true",
        help="omit Bond numbers that fail to converge instead of aborting",
    )
    args = ap.parse_args(argv)
    usetex = not args.no_usetex
    _configure_matplotlib(usetex=usetex)
    shapes = _sweep(PLOT_BONDS, skip_failed=args.skip_failed)
    plot_opening_angle(shapes, args.out, usetex=usetex)
    return 0


if __name__ == "__main__":
    sys.exit(main())
