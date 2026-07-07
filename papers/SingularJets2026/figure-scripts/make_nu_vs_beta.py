#!/usr/bin/env python3
"""Build the End Matter nu(beta) and alpha(Oh) figure.

This is a publication-plot script: it uses Matplotlib with LaTeX typography and
SciPy's Legendre function to solve P_nu(-cos(beta)) = 0.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy import optimize, special
except ImportError as exc:  # pragma: no cover - depends on local environment
    raise SystemExit(
        "make_nu_vs_beta.py needs numpy, scipy, matplotlib, and a working LaTeX installation."
    ) from exc


ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().with_name("beta_alpha_Oh.csv")

APS_DOUBLE_COL = 6.75

APS = {
    "LabelFont": 10,
    "AxesFont": 9,
    "LegendFont": 8,
    "PanelFont": 10,
}

LINE = {
    "linewidth": 1.5,
    "theory_linewidth": 1.2,
    "spine_width": 1.0,
    "tick_width": 0.8,
    "tick_length_major": 5,
    "tick_length_minor": 2.5,
    "markersize": 4.8,
    "markeredgewidth": 0.7,
}

BLUE = "#1f77b4"
ORANGE = "#d95f02"
RED = "#d7191c"
GREY = "#6f6f6f"
LIGHT_GREY = "#b7b7b7"


def configure_matplotlib(use_tex: bool = True) -> None:
    matplotlib.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Computer Modern Roman"],
            "mathtext.fontset": "cm",
            "text.usetex": use_tex,
            "text.latex.preamble": r"\usepackage{amsmath}",
            "axes.linewidth": LINE["spine_width"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.unicode_minus": False,
        }
    )


def legendre_condition(nu: float, beta_deg: float) -> float:
    """Return P_nu(-cos(beta))."""

    return float(special.lpmv(0, nu, -np.cos(np.deg2rad(beta_deg))))


def nu_for_beta(beta_deg: float, previous: float | None = None) -> float:
    """Solve P_nu(-cos(beta)) = 0 for 0 < nu < 1."""

    if beta_deg <= 0.0:
        return 0.0

    center = 0.0 if previous is None else previous
    lo = max(1.0e-10, center - 0.04)
    hi = min(0.999999, center + 0.10)
    flo = legendre_condition(lo, beta_deg)
    fhi = legendre_condition(hi, beta_deg)

    while not (np.isfinite(flo) and np.isfinite(fhi) and flo * fhi <= 0.0):
        lo = max(1.0e-10, lo - 0.05)
        hi = min(0.999999, hi + 0.05)
        flo = legendre_condition(lo, beta_deg)
        fhi = legendre_condition(hi, beta_deg)
        if lo <= 1.0e-10 and hi >= 0.999999 and flo * fhi > 0.0:
            raise RuntimeError(f"No root bracketed for beta={beta_deg:g}")

    return float(
        optimize.brentq(
            lambda nu: legendre_condition(nu, beta_deg),
            lo,
            hi,
            xtol=1.0e-12,
            rtol=1.0e-12,
            maxiter=80,
        )
    )


def nu_curve() -> tuple[np.ndarray, np.ndarray]:
    """Compute a dense curve with extra resolution near beta=0."""

    beta = np.r_[0.0, np.geomspace(1.0e-3, 1.0, 90), np.linspace(1.05, 60.0, 260)]
    nu = np.empty_like(beta)
    previous: float | None = None
    for i, b in enumerate(beta):
        nu[i] = nu_for_beta(float(b), previous)
        previous = float(nu[i])
    return beta, nu


def read_alpha_data(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open(newline="") as f:
        for row in csv.DictReader(f):
            rows.append(
                {
                    "bo": float(row["Bo"]),
                    "level": float(row["MAXlevel"]),
                    "oh": float(row["Oh"]),
                    "beta": float(row["beta_deg"]),
                    "nu": float(row["nu"]),
                    "alpha": float(row["alpha"]),
                }
            )
    return rows


def style_axes(ax: plt.Axes) -> None:
    ax.tick_params(
        axis="both",
        which="major",
        labelsize=APS["AxesFont"],
        width=LINE["tick_width"],
        length=LINE["tick_length_major"],
        direction="out",
        pad=3,
    )
    ax.tick_params(
        which="minor",
        width=LINE["tick_width"] * 0.7,
        length=LINE["tick_length_minor"],
        direction="out",
    )
    for spine in ax.spines.values():
        spine.set_linewidth(LINE["spine_width"])
    ax.minorticks_on()
    ax.set_box_aspect(0.78)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.16,
        1.08,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=APS["PanelFont"],
        fontweight="bold",
    )


def build_figure(data_path: Path, output: Path, use_tex: bool = True) -> None:
    configure_matplotlib(use_tex=use_tex)
    rows = read_alpha_data(data_path)
    beta, nu = nu_curve()

    fig, (ax_nu, ax_alpha) = plt.subplots(1, 2, figsize=(APS_DOUBLE_COL, 2.7))
    fig.set_facecolor("white")

    ax_nu.plot(beta, nu, color="black", lw=LINE["linewidth"], zorder=2)
    ax_nu.plot(
        [r["beta"] for r in rows],
        [r["nu"] for r in rows],
        linestyle="None",
        marker="o",
        ms=3.2,
        mfc="white",
        mec=GREY,
        mew=0.55,
        zorder=3,
    )
    ax_nu.plot(49.3, 0.5, "o", ms=5.0, color=RED, mec=RED, zorder=4)
    ax_nu.annotate(
        r"Taylor cone",
        xy=(49.3, 0.5),
        xytext=(33.5, 0.590),
        textcoords="data",
        color=RED,
        fontsize=APS["LegendFont"],
        arrowprops=dict(arrowstyle="-", color=RED, lw=0.8, shrinkA=2, shrinkB=4),
    )
    ax_nu.set_xlim(0.0, 60.0)
    ax_nu.set_ylim(0.0, 0.64)
    ax_nu.set_xticks(np.arange(0, 61, 10))
    ax_nu.set_yticks(np.arange(0.0, 0.61, 0.1))
    ax_nu.set_xlabel(r"$\beta\,(^\circ)$", fontsize=APS["LabelFont"], labelpad=4)
    ax_nu.set_ylabel(r"$\nu$", fontsize=APS["LabelFont"], labelpad=4)
    style_axes(ax_nu)
    add_panel_label(ax_nu, r"(a)")

    series = [
        (0.0, 13.0, "o", BLUE, r"$Bo=0$, L13"),
        (0.0, 14.0, "^", BLUE, r"$Bo=0$, L14"),
        (0.0, 15.0, "D", BLUE, r"$Bo=0$, L15"),
        (1.0e-3, 13.0, "s", ORANGE, r"$Bo=0.001$, L13"),
    ]
    for bo, level, marker, colour, label in series:
        group = [
            row
            for row in rows
            if abs(row["bo"] - bo) < 1.0e-12 and abs(row["level"] - level) < 1.0e-12
            and 0.02 < row["oh"] < 0.04
        ]
        if not group:
            continue
        group.sort(key=lambda row: row["oh"])
        ax_alpha.plot(
            [row["oh"] for row in group],
            [row["alpha"] for row in group],
            linestyle="None",
            marker=marker,
            ms=LINE["markersize"],
            mfc="white",
            mec=colour,
            mew=LINE["markeredgewidth"],
            label=label,
            zorder=3,
        )

    ax_alpha.axhline(2.0 / 3.0, color=GREY, lw=LINE["theory_linewidth"], ls=(0, (4, 3)), zorder=1)
    ax_alpha.text(
        0.0367,
        2.0 / 3.0 + 0.0010,
        r"$2/3$",
        ha="left",
        va="center",
        fontsize=APS["LegendFont"],
        color=GREY,
        bbox=dict(facecolor="white", edgecolor="none", pad=0.5),
    )
    ax_alpha.set_xlim(0.020, 0.040)
    ax_alpha.set_ylim(0.620, 0.670)
    ax_alpha.set_xticks([0.02, 0.03, 0.04])
    ax_alpha.set_yticks([0.62, 0.64, 0.66])
    ax_alpha.set_xlabel(r"$Oh$", fontsize=APS["LabelFont"], labelpad=4)
    ax_alpha.set_ylabel(r"$\alpha$", fontsize=APS["LabelFont"], labelpad=4)
    style_axes(ax_alpha)
    add_panel_label(ax_alpha, r"(b)")
    ax_alpha.legend(
        loc="upper left",
        bbox_to_anchor=(0.03, 0.86),
        ncol=2,
        frameon=False,
        fontsize=APS["LegendFont"],
        handlelength=1.0,
        handletextpad=0.35,
        columnspacing=0.8,
        borderpad=0.1,
        labelspacing=0.25,
    )

    fig.subplots_adjust(left=0.075, right=0.995, bottom=0.20, top=0.90, wspace=0.38)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", pad_inches=0.04, dpi=300)
    plt.close(fig)
    print(f"Wrote {output}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DATA, help="CSV file with beta, nu, and alpha measurements")
    parser.add_argument("--output", type=Path, default=ROOT / "nu_vs_beta.pdf", help="Output PDF path")
    parser.add_argument("--no-tex", action="store_true", help="Use Matplotlib mathtext instead of LaTeX")
    args = parser.parse_args()
    build_figure(args.data, args.output, use_tex=not args.no_tex)


if __name__ == "__main__":
    main()
