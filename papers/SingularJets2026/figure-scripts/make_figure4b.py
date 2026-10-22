#!/usr/bin/env python3
"""Render the measured cone-angle exponents in Figure 4(b) from CSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data" / "figure4b.csv"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DATA)
    parser.add_argument("--output", type=Path, default=ROOT / "output" / "figure4b.pdf")
    parser.add_argument("--no-tex", action="store_true")
    args = parser.parse_args()
    matplotlib.rcParams.update({
        "font.family": "serif", "font.serif": ["Computer Modern Roman"],
        "mathtext.fontset": "cm", "text.usetex": not args.no_tex,
        "text.latex.preamble": r"\usepackage{amsmath}",
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })
    groups: dict[tuple[float, float], list[tuple[float, float]]] = {}
    with args.data.open(newline="") as stream:
        for row in csv.DictReader(stream):
            if not 0.02 < float(row["Oh"]) < 0.04:
                continue
            key = (float(row["Bo"]), float(row["MAXlevel"]))
            groups.setdefault(key, []).append((float(row["Oh"]), float(row["alpha"])))
    styles = [("o", "#0072B2", r"$Bo=0$, L13"), ("^", "#0072B2", r"$Bo=0$, L14"),
              ("D", "#0072B2", r"$Bo=0$, L15"), ("s", "#D55E00", r"$Bo=0.001$, L13")]
    labels = [(0.0, 13.0), (0.0, 14.0), (0.0, 15.0), (1e-3, 13.0)]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(3.3, 2.6))
    for key, (marker, colour, label) in zip(labels, styles):
        points = sorted(groups.get(key, []))
        if points:
            ax.plot([p[0] for p in points], [p[1] for p in points], linestyle="none",
                    marker=marker, markersize=3.5, markerfacecolor="white",
                    markeredgecolor=colour, markeredgewidth=0.8, label=label)
    ax.axhline(2 / 3, color="0.45", linewidth=1, linestyle="--")
    ax.text(0.0367, 2 / 3 + 0.001, r"$2/3$", color="0.45", fontsize=6.5,
            ha="left", va="center", bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.5})
    ax.set(xlim=(0.02, 0.04), ylim=(0.62, 0.67), xlabel=r"$Oh$", ylabel=r"$\alpha$")
    ax.legend(frameon=False, fontsize=5.8, loc="upper left", bbox_to_anchor=(0.02, 0.84), ncol=2)
    ax.tick_params(direction="out", labelsize=8, width=0.8, length=3)
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
    fig.tight_layout(pad=0.5)
    fig.savefig(args.output, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
