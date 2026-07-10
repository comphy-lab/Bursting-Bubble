#!/usr/bin/env python3
"""Generate the schematic for the Self-similar Worthington jets paper."""

from __future__ import annotations

import argparse
import math
import shutil
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
from matplotlib.patches import Arc, Circle, FancyArrowPatch, PathPatch
from matplotlib.path import Path as MplPath


REPO_ROOT = Path(__file__).resolve().parents[1]

INTERFACE = "#000000"
AXIS = "#5f6770"
BETA = "#c23b22"
THETA = "#7b3294"
R_COLOR = "#0072b2"
RP_COLOR = "#d55e00"
Z_COLOR = "#009e73"
ELL_COLOR = "#cc79a7"
VJ_COLOR = "#1b1b1b"


def configure_matplotlib(use_tex: bool) -> None:
    rc = {
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "mathtext.fontset": "cm",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.linewidth": 0.8,
    }
    if use_tex and shutil.which("latex"):
        rc["text.usetex"] = True
        rc["text.latex.preamble"] = r"\usepackage{amsmath}"
    else:
        rc["text.usetex"] = False
    mpl.rcParams.update(rc)


def cubic_path(points: list[tuple[float, float]]) -> MplPath:
    """Build a path from one start point plus 3-point cubic Bezier segments."""
    if (len(points) - 1) % 3 != 0:
        raise ValueError("Bezier point list must be start + 3*n points")
    codes = [MplPath.MOVETO] + [MplPath.CURVE4] * (len(points) - 1)
    return MplPath(points, codes)


def draw_curve(
    ax: plt.Axes,
    points: list[tuple[float, float]],
    *,
    lw: float = 1.12,
    ls: str = "-",
    alpha: float = 1.0,
    zorder: int = 3,
) -> None:
    patch = PathPatch(
        cubic_path(points),
        facecolor="none",
        edgecolor=INTERFACE,
        lw=lw,
        linestyle=ls,
        capstyle="round",
        joinstyle="round",
        alpha=alpha,
        zorder=zorder,
    )
    ax.add_patch(patch)


def arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    color: str = INTERFACE,
    lw: float = 0.92,
    mutation_scale: float = 7.5,
    style: str = "-|>",
    alpha: float = 1.0,
    zorder: int = 5,
) -> None:
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle=style,
            color=color,
            lw=lw,
            mutation_scale=mutation_scale,
            shrinkA=0,
            shrinkB=0,
            capstyle="round",
            joinstyle="round",
            alpha=alpha,
            zorder=zorder,
        )
    )


def pol2cart(radius: float, angle_deg_from_x: float) -> tuple[float, float]:
    angle = math.radians(angle_deg_from_x)
    return radius * math.cos(angle), radius * math.sin(angle)


def label(
    ax: plt.Axes,
    x: float,
    y: float,
    text: str,
    *,
    color: str,
    fontsize: float = 11.0,
    ha: str = "center",
    va: str = "center",
    zorder: int = 10,
) -> None:
    ax.text(
        x,
        y,
        text,
        color=color,
        fontsize=fontsize,
        ha=ha,
        va=va,
        zorder=zorder,
        path_effects=[
            pe.SimpleLineShadow(offset=(0.55, -0.55), alpha=0.22, rho=0.95),
            pe.Stroke(linewidth=2.4, foreground=(0.97, 0.96, 0.92, 0.62)),
            pe.Normal(),
        ],
    )


def build_figure(output_dir: Path, use_tex: bool = True) -> None:
    configure_matplotlib(use_tex)

    output_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(3.0, 3.5))
    ax.set_aspect("equal")
    ax.set_xlim(-3.35, 3.35)
    ax.set_ylim(-0.35, 4.75)
    ax.axis("off")

    # Symmetry axis and conical far-field guide.
    ax.plot([0.0, 0.0], [-0.30, 4.55], color=AXIS, lw=0.75, dashes=(5, 5), zorder=1)
    beta_deg = 38.0
    beta_ray_angle = 90.0 - beta_deg
    cone_x = lambda y: y * math.tan(math.radians(beta_deg))
    for sign in (-1, 1):
        y_end = 2.18
        x_end = sign * cone_x(y_end)
        ax.plot(
            [0.0, x_end],
            [0.0, y_end],
            color=BETA,
            lw=0.65,
            dashes=(4, 4),
            alpha=0.78,
            zorder=1,
        )

    # Free surface: outer conical cavity plus the rising Worthington jet.
    cone_y0, cone_y1 = 1.72, 4.20
    for sign in (-1, 1):
        ax.plot(
            [sign * cone_x(cone_y0), sign * cone_x(cone_y1)],
            [cone_y0, cone_y1],
            color=INTERFACE,
            lw=1.12,
            solid_capstyle="round",
            zorder=3,
        )
    draw_curve(
        ax,
        [
            (-cone_x(cone_y0), cone_y0),
            (-1.14, 1.48),
            (-0.89, 1.22),
            (-0.60, 1.19),
            (-0.41, 1.18),
            (-0.33, 1.36),
            (-0.33, 1.60),
            (-0.33, 1.95),
            (-0.20, 2.58),
            (-0.15, 3.08),
            (-0.13, 3.55),
            (-0.10, 3.92),
            (-0.05, 4.22),
        ],
    )
    draw_curve(
        ax,
        [
            (cone_x(cone_y0), cone_y0),
            (1.14, 1.48),
            (0.89, 1.22),
            (0.60, 1.19),
            (0.41, 1.18),
            (0.33, 1.36),
            (0.33, 1.60),
            (0.33, 1.95),
            (0.20, 2.58),
            (0.15, 3.08),
            (0.13, 3.55),
            (0.10, 3.92),
            (0.05, 4.22),
        ],
    )

    # Local jet radius and axial velocity at the jet base.
    arrow(ax, (-0.88, 1.10), (-0.10, 1.10), color=ELL_COLOR, lw=0.82, mutation_scale=6.8, style="<->")
    label(ax, -1.00, 1.00, r"$r_j$", color=ELL_COLOR, fontsize=11, va="center")

    arrow(ax, (0.0, 1.10), (0.0, 1.76), color=VJ_COLOR, lw=1.02, mutation_scale=8.0)
    label(ax, 0.07, 1.24, r"$v_j$", color=VJ_COLOR, fontsize=12.0, ha="left")

    # Cylindrical and spherical coordinate directions.
    arrow(ax, (0.0, 0.0), (0.0, 0.82), color=Z_COLOR, lw=0.86, mutation_scale=7.0)
    label(ax, -0.30, 0.62, r"$z$", color=Z_COLOR, fontsize=10.5, ha="right")

    arrow(ax, (0.0, 0.0), (1.05, 0.0), color=RP_COLOR, lw=0.86, mutation_scale=7.0)
    label(ax, 0.58, -0.25, r"$r_p$", color=RP_COLOR, fontsize=11, va="top")

    theta_deg = 58.0
    theta_ray_angle = 90.0 - theta_deg
    r_end = pol2cart(2.45, theta_ray_angle)
    arrow(ax, (0.0, 0.0), r_end, color=R_COLOR, lw=0.88, mutation_scale=7.2)
    label(ax, r_end[0] + 0.12, r_end[1] + 0.03, r"$r$", color=R_COLOR, fontsize=12, ha="left")

    # Angle annotations, measured from the symmetry axis.
    ax.add_patch(Arc((0.0, 0.0), 0.78, 0.78, theta1=beta_ray_angle, theta2=90, color=BETA, lw=0.72))
    label(ax, 0.13, 0.50, r"$\beta$", color=BETA, fontsize=10.8, ha="center")

    ax.add_patch(Arc((0.0, 0.0), 1.86, 1.86, theta1=theta_ray_angle, theta2=90, color=THETA, lw=0.72))
    label(ax, 0.92, 0.80, r"$\theta$", color=THETA, fontsize=10.8, ha="left")

    ax.plot(0.0, 0.0, marker="o", ms=2.25, color=INTERFACE, zorder=6)

    pdf_path = output_dir / "sketch_paper_worthington.pdf"
    png_path = output_dir / "sketch_paper_worthington.png"
    fig.savefig(pdf_path, bbox_inches="tight", pad_inches=0.03, dpi=300)
    fig.savefig(png_path, bbox_inches="tight", pad_inches=0.03, dpi=450)
    plt.close(fig)

    print(f"Wrote {pdf_path}")
    print(f"Wrote {png_path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT,
        help="Directory for sketch_paper_worthington.pdf/png",
    )
    parser.add_argument("--no-tex", action="store_true", help="Use Matplotlib mathtext instead of LaTeX")
    args = parser.parse_args()
    build_figure(args.output_dir, use_tex=not args.no_tex)


if __name__ == "__main__":
    main()
