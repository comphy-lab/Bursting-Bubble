#!/usr/bin/env python3
"""Plot the cone-angle map by solving P_nu(-cos(beta)) = 0.

For each semiangle beta, the degree nu is the real root in (0, 1) of

    P_nu(-cos(beta)) = 0.

The root is bracketed on that interval and found with Brent's method. The
theoretical curve needs no fitted parameters. Open markers show the DNS
semiangles from the supplied CSV. Taylor's cone is marked at beta=49.3 degrees,
nu=1/2.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq
from scipy.special import lpmv

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "figure4b.csv"


def nu_at_beta(beta_deg: float) -> float:
    """Return the first real degree nu in (0, 1) satisfying the cone condition."""
    if not 0.0 < beta_deg < 90.0:
        raise ValueError("beta must be strictly between 0 and 90 degrees")
    x = -np.cos(np.deg2rad(beta_deg))
    return float(brentq(lambda nu: lpmv(0, nu, x), 1e-10, 1.0 - 1e-10,
                        xtol=1e-12, rtol=1e-12, maxiter=100))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DATA)
    parser.add_argument("--output", type=Path, default=ROOT / "output" / "figure4a.pdf")
    parser.add_argument("--no-tex", action="store_true")
    args = parser.parse_args()
    matplotlib.rcParams.update({
        "font.family": "serif", "font.serif": ["Computer Modern Roman"],
        "mathtext.fontset": "cm", "text.usetex": not args.no_tex,
        "text.latex.preamble": r"\usepackage{amsmath}",
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    beta = np.r_[0.0, np.geomspace(1.0e-3, 1.0, 90), np.linspace(1.05, 60.0, 260)]
    nu = np.r_[0.0, [nu_at_beta(float(value)) for value in beta[1:]]]
    with args.data.open(newline="") as stream:
        measured = list(csv.DictReader(stream))
    beta_dns = np.array([float(row["beta_deg"]) for row in measured])
    nu_dns = np.array([float(row["nu"]) for row in measured])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(3.3, 2.6))
    ax.plot(beta, nu, color="black", linewidth=1.35)
    ax.plot(beta_dns, nu_dns, linestyle="none", marker="o", markersize=2.8,
            markerfacecolor="white", markeredgecolor="0.45", markeredgewidth=0.65)
    ax.plot(49.3, 0.5, "o", color="#b33b2e", markersize=4)
    ax.annotate("Taylor cone", xy=(49.3, 0.5), xytext=(33.5, 0.59), color="#b33b2e",
                fontsize=7, arrowprops={"arrowstyle": "-", "color": "#b33b2e", "lw": 0.7})
    ax.set(xlim=(0, 60), ylim=(0, 0.65), xlabel=r"$\beta\,(^\circ)$", ylabel=r"$\nu$")
    ax.tick_params(direction="out", labelsize=8, width=0.8, length=3)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    fig.tight_layout(pad=0.5)
    fig.savefig(args.output, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
