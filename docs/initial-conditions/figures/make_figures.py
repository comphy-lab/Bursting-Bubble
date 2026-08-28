#!/usr/bin/env python3
"""Publication figures for the Young-Laplace initial-condition report."""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

matplotlib.rcParams["font.family"] = "serif"
matplotlib.rcParams["font.serif"] = ["Computer Modern Roman"]
matplotlib.rcParams["text.usetex"] = True
matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"

HERE = Path(__file__).resolve().parent
IC = HERE.parents[2] / "simulationCases" / "initialConditions"
if str(IC) not in sys.path:
    sys.path.insert(0, str(IC))

from plot_opening_angle import (  # noqa: E402
    _configure_matplotlib,
    _sweep,
    plot_opening_angle,
    PLOT_BONDS,
)
from young_laplace import solve_equilibrium  # noqa: E402
from zero_bond import sphere_plane  # noqa: E402

BLUE = "#1A64B3"
ORANGE = "#fc8d59"
COLORS = ["#1A64B3", "#2ca02c", "#ff7f0e", "#9467bd", "#d62728"]


def _style_axes(ax, *, equal: bool = True) -> None:
    ax.tick_params(which="both", direction="out", width=3, labelsize=30, pad=10)
    ax.tick_params(which="major", length=12)
    ax.tick_params(which="minor", length=6)
    for spine in ax.spines.values():
        spine.set_linewidth(3)
    ax.minorticks_on()
    if equal:
        ax.set_aspect("equal")


def _save(fig, name: str) -> None:
    pdf = HERE / name
    png = pdf.with_suffix(".png")
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.1, dpi=300)
    fig.savefig(png, bbox_inches="tight", pad_inches=0.1, dpi=200)
    plt.close(fig)
    print(f"wrote {pdf}")
    print(f"wrote {png}")


def make_geometry() -> None:
    shape = solve_equilibrium(1.0)
    bonds = [1e-3, 0.1, 1.0, 10.0]
    family = []
    previous = None
    for bond in bonds:
        previous = solve_equilibrium(bond, previous=previous)
        family.append(previous)
    zero = sphere_plane(delta=0.01, rmax=4.0)

    fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(24, 11))

    for R, Z, lw in (
        (shape.R_bubble, shape.Z_bubble, 3.2),
        (shape.R_cap, shape.Z_cap, 2.4),
        (shape.R_tail, shape.Z_tail, 2.4),
    ):
        ax0.plot(R, Z, color="k", lw=lw)
        ax0.plot(-R, Z, color="k", lw=lw)
    xc = shape.R_bubble[-1]
    yc = shape.Z_bubble[-1]
    ax0.plot([xc, -xc], [yc, yc], linestyle="none", marker="o",
             markersize=10, color=ORANGE, zorder=4)
    ax0.axhline(shape.hinf, color=BLUE, lw=1.6, ls="--")
    ax0.annotate(
        r"$\varphi_c$",
        xy=(xc, yc),
        xytext=(xc + 0.55, yc + 0.35),
        fontsize=28,
        arrowprops=dict(arrowstyle="->", lw=1.6, color=ORANGE),
        color=ORANGE,
    )
    ax0.text(1.55, shape.hinf + 0.08, r"$h_\infty$", fontsize=28, color=BLUE)
    ax0.text(0.12, 0.18, r"south pole", fontsize=24)
    ax0.set_xlim(-2.4, 2.4)
    ax0.set_ylim(-0.15, 2.35)
    ax0.set_xlabel(r"$\mathcal{R}/R_0$", fontsize=40, labelpad=15)
    ax0.set_ylabel(r"$\mathcal{Z}/R_0$", fontsize=40, labelpad=15)
    _style_axes(ax0)
    ax0.set_title(r"(a) $\mathcal{B}o=1$", fontsize=30, pad=12)

    curves = [(zero.radial, zero.axial, r"$\mathcal{B}o=0$")]
    for shape_i, bond in zip(family, bonds):
        mask = shape_i.radial <= 4.0
        label = rf"$\mathcal{{B}}o={bond:g}$"
        curves.append((shape_i.radial[mask], shape_i.axial[mask], label))
    for (R, Z, label), color in zip(curves, COLORS):
        ax1.plot(R, Z, color=color, lw=2.8, label=label)
        ax1.plot(-R, Z, color=color, lw=2.8)
    ax1.axhline(0.0, color="0.4", lw=1.2, ls=":")
    ax1.set_xlim(-3.2, 3.2)
    ax1.set_ylim(-2.3, 0.35)
    ax1.set_xlabel(r"radial $/R_0$", fontsize=40, labelpad=15)
    ax1.set_ylabel(r"axial $/R_0$", fontsize=40, labelpad=15)
    _style_axes(ax1)
    ax1.legend(loc="lower right", frameon=False, fontsize=22)
    ax1.set_title(r"(b) Basilisk convention", fontsize=30, pad=12)

    fig.tight_layout()
    _save(fig, "geometry.pdf")


def make_zero_bond() -> None:
    delta = 0.08
    shape = sphere_plane(delta=delta, rmax=2.4, n=800)
    fig, ax = plt.subplots(figsize=(12, 12))
    ax.plot(shape.radial, shape.axial, "k-", lw=3.2)
    ax.plot(-shape.radial, shape.axial, "k-", lw=3.2)

    # After the write convention, south pole is near -2 and surface at 0.
    south = shape.south_pole_axial()
    circle = plt.Circle(
        (0.0, south + 1.0),
        1.0,
        fill=False,
        lw=1.8,
        ls="--",
        color=BLUE,
        zorder=1,
    )
    ax.add_patch(circle)
    ax.plot(0.0, south, "o", color=ORANGE, markersize=10, zorder=4)
    ax.annotate(
        r"south pole",
        xy=(0.0, south),
        xytext=(0.55, south - 0.15),
        fontsize=26,
        arrowprops=dict(arrowstyle="->", lw=1.6),
    )
    ax.annotate(
        r"fillet $\sim 2\delta$",
        xy=(2.0 * delta, 0.0),
        xytext=(0.85, 0.28),
        fontsize=26,
        color=ORANGE,
        arrowprops=dict(arrowstyle="->", lw=1.6, color=ORANGE),
    )
    ax.annotate(
        r"free surface",
        xy=(1.6, 0.0),
        xytext=(1.15, 0.42),
        fontsize=26,
        arrowprops=dict(arrowstyle="->", lw=1.6),
    )
    ax.text(-1.85, south + 1.0, r"unit sphere", fontsize=26, color=BLUE)
    ax.set_xlim(-2.3, 2.3)
    ax.set_ylim(-2.35, 0.65)
    ax.set_xlabel(r"radial $/R_0$", fontsize=40, labelpad=15)
    ax.set_ylabel(r"axial $/R_0$", fontsize=40, labelpad=15)
    _style_axes(ax)
    handles = [
        Line2D([0], [0], color="k", lw=3.2, label=rf"$\delta={delta:g}$"),
        Line2D([0], [0], color=BLUE, lw=1.8, ls="--", label="unit circle"),
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=26)
    _save(fig, "zero_bond.pdf")


def make_opening_angle() -> None:
    out = HERE / "opening_angle.pdf"
    _configure_matplotlib(usetex=True)
    shapes = _sweep(PLOT_BONDS, skip_failed=False)
    plot_opening_angle(shapes, out, usetex=True)


def main() -> int:
    HERE.mkdir(parents=True, exist_ok=True)
    make_geometry()
    make_zero_bond()
    make_opening_angle()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
