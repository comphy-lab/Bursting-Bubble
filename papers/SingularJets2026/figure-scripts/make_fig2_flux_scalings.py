#!/usr/bin/env python3
"""Build the Fig. 2 scaffold from the Oh=0.03 grid-convergence logs.

Paper notation differs from the solver-column names:
  q_j  = q_l
  Q_j  = 2*pi*q_jet = pi*r_j*q_j for the plug-like jet profile
  We_j = q_l**2/r_j

The cone and inertio-capillary prefactors are fit in the near-inception
inertial window. The PRF 2023 prefactor is fit in the finite-radius plateau
window where q_j is approximately constant and We_j ~ r_j^{-1}.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
except ImportError as exc:  # pragma: no cover - depends on local environment
    raise SystemExit(
        "make_fig2_flux_scalings.py needs numpy, matplotlib, and a working LaTeX installation."
    ) from exc


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DATA_DIR = ROOT / "data-Oh-0.03"

APS_DOUBLE_COL = 6.75
ALPHA = 0.629
MAX_RJ = 1.0
CONE_FIT_WINDOW = (0.005, 0.023952)
PRF_FIT_WINDOW = (0.11, 0.19)
CONE_DRAW_WINDOW = (0.005, 0.10)
PRF_DRAW_WINDOW = (0.052, 0.60)

APS = {
    "LabelFont": 10,
    "AxesFont": 8.5,
    "LegendFont": 7.5,
    "PanelFont": 10,
}

LINE = {
    "linewidth": 1.35,
    "theory_linewidth": 1.45,
    "spine_width": 0.9,
    "tick_width": 0.75,
    "tick_length_major": 4.2,
    "tick_length_minor": 2.1,
    "markersize": 2.65,
    "markeredgewidth": 0.28,
}

LEVEL_COLOURS = {
    13: "#0072B2",
    14: "#D55E00",
    15: "#009E73",
}
GREY = "#666666"
LIGHT_GREY = "#d9d9d9"
BLACK = "#111111"


@dataclass(frozen=True)
class RunSpec:
    filename: str
    label: str
    level: int
    focus: int | None
    marker: str

    @property
    def colour(self) -> str:
        return LEVEL_COLOURS[self.level]


RUNS = (
    RunSpec("3013_L13_log.txt", r"L13", 13, None, "o"),
    RunSpec("4015_L14_log.txt", r"L14", 14, None, "s"),
    RunSpec("5001_L15_focus13_log.txt", r"L15, focus 13", 15, 13, "D"),
    RunSpec("5003_L15_focus14_log.txt", r"L15, focus 14", 15, 14, "P"),
    RunSpec("5008_L15_focus15_log.txt", r"L15, focus 15", 15, 15, "X"),
)

COLUMNS = (
    "i",
    "dt",
    "t",
    "ke",
    "maxlevel",
    "r_b",
    "z_b",
    "r_base",
    "z_base",
    "q_jet",
    "q_l",
)


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


def read_log(path: Path) -> dict[str, np.ndarray]:
    rows_by_time: dict[float, list[float]] = {}
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("MAXlevel") or line.startswith("i "):
                continue
            parts = line.split()
            if len(parts) != len(COLUMNS):
                continue
            row = [float(value) for value in parts]
            rows_by_time[round(row[2], 8)] = row

    if not rows_by_time:
        raise ValueError(f"No numeric rows found in {path}")

    rows = [rows_by_time[key] for key in sorted(rows_by_time)]
    array = np.asarray(rows, dtype=float)
    return {name: array[:, idx] for idx, name in enumerate(COLUMNS)}


def reconnection_time(data: dict[str, np.ndarray], pin_r: float = 0.005) -> float | None:
    mask = (data["t"] > 0.40) & (data["r_base"] < pin_r)
    if not np.any(mask):
        return None
    return float(np.max(data["t"][mask]))


def processed_series(data: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    r_j = data["r_base"]
    q_j = data["q_l"]
    Q_j = 2.0 * np.pi * data["q_jet"]
    We_j = q_j**2 / r_j
    incept_t = reconnection_time(data)

    # Post-inception branch used in the grid-convergence diagnostic. The
    # reconnection-time filter follows the reference gridconv3 workflow and
    # removes isolated pre-jet spikes from lower-resolution runs.
    mask = (
        np.isfinite(r_j)
        & np.isfinite(q_j)
        & np.isfinite(Q_j)
        & np.isfinite(We_j)
        & (r_j > 0.0)
        & (r_j <= MAX_RJ)
        & (q_j > 0.0)
        & (Q_j > 0.0)
        & (We_j > 0.0)
    )
    if incept_t is not None:
        mask &= data["t"] > incept_t

    order = np.argsort(r_j[mask])
    return {
        "r_j": r_j[mask][order],
        "q_j": q_j[mask][order],
        "Q_j": Q_j[mask][order],
        "We_j": We_j[mask][order],
    }


def log_bin_indices(x: np.ndarray, y: np.ndarray, target: int = 75) -> np.ndarray:
    if len(x) <= target:
        return np.arange(len(x), dtype=int)
    valid = np.isfinite(x) & np.isfinite(y) & (x > 0.0) & (y > 0.0)
    if not np.any(valid):
        return np.array([], dtype=int)

    valid_idx = np.flatnonzero(valid)
    edges = np.geomspace(np.min(x[valid]), np.max(x[valid]), target + 1)
    chosen: list[int] = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        in_bin = valid_idx[(x[valid_idx] >= lo) & (x[valid_idx] < hi)]
        if len(in_bin) == 0:
            continue
        centre = np.sqrt(lo * hi)
        chosen.append(int(in_bin[np.argmin(np.abs(np.log(x[in_bin] / centre)))]))
    return np.asarray(sorted(set(chosen)), dtype=int)


def style_axes(ax: plt.Axes) -> None:
    ax.tick_params(
        axis="both",
        which="major",
        labelsize=APS["AxesFont"],
        width=LINE["tick_width"],
        length=LINE["tick_length_major"],
        direction="out",
        pad=2.5,
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
    ax.set_box_aspect(1.0)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.19,
        1.06,
        label,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=APS["PanelFont"],
        fontweight="bold",
    )


def pool_series(series_by_run: list[tuple[RunSpec, dict[str, np.ndarray]]]) -> dict[str, np.ndarray]:
    keys = ("r_j", "q_j", "Q_j", "We_j")
    return {key: np.concatenate([series[key] for _, series in series_by_run]) for key in keys}


def normalisation(
    series: dict[str, np.ndarray],
    quantity: str,
    slope: float,
    window: tuple[float, float],
) -> float:
    fit = (series["r_j"] >= window[0]) & (series["r_j"] <= window[1])
    if not np.any(fit):
        raise ValueError(f"No points in the fit window for {quantity}")
    return float(
        np.exp(np.mean(np.log(series[quantity][fit]) - slope * np.log(series["r_j"][fit])))
    )


def run_weighted_normalisation(
    series_by_run: list[tuple[RunSpec, dict[str, np.ndarray]]],
    quantity: str,
    slope: float,
    window: tuple[float, float],
) -> float:
    prefactors = []
    for _, series in series_by_run:
        fit = (series["r_j"] >= window[0]) & (series["r_j"] <= window[1])
        if not np.any(fit):
            continue
        prefactors.append(
            np.mean(np.log(series[quantity][fit]) - slope * np.log(series["r_j"][fit]))
        )
    if not prefactors:
        raise ValueError(f"No runs have points in the fit window for {quantity}")
    return float(np.exp(np.median(prefactors)))


def reference_run_normalisation(
    series_by_run: list[tuple[RunSpec, dict[str, np.ndarray]]],
    quantity: str,
    slope: float,
    window: tuple[float, float],
    level: int = 15,
    focus: int = 15,
) -> float:
    for run, series in series_by_run:
        if run.level != level or run.focus != focus:
            continue
        fit = (series["r_j"] >= window[0]) & (series["r_j"] <= window[1])
        if not np.any(fit):
            raise ValueError(f"No {run.label} points in the fit window for {quantity}")
        return float(
            np.exp(np.mean(np.log(series[quantity][fit]) - slope * np.log(series["r_j"][fit])))
        )
    raise ValueError(f"No L{level}, focus {focus} run found for {quantity}")


def marker_alpha(run: RunSpec) -> float:
    if run.level == 15 and run.focus == 15:
        return 0.94
    if run.level == 15:
        return 0.74
    return 0.48


def draw_theory(
    ax: plt.Axes,
    series_by_run: list[tuple[RunSpec, dict[str, np.ndarray]]],
    quantity: str,
    slopes: tuple[float, float, float],
    show_labels: bool,
    cone_fit_window: tuple[float, float],
) -> None:
    cone_slope, ic_slope, prf_slope = slopes
    r_cone = np.geomspace(*CONE_DRAW_WINDOW, 120)
    r_prf = np.geomspace(*PRF_DRAW_WINDOW, 120)

    cone_prefactor = reference_run_normalisation(
        series_by_run, quantity, cone_slope, cone_fit_window
    )
    ax.plot(
        r_cone,
        cone_prefactor * r_cone**cone_slope,
        color=BLACK,
        ls="-",
        lw=LINE["theory_linewidth"],
        zorder=8,
        label=rf"cone ($\alpha={ALPHA:.3f}$)" if show_labels else None,
    )

    if quantity == "We_j":
        ax.axhline(
            1.0,
            color=GREY,
            ls="--",
            lw=LINE["theory_linewidth"],
            zorder=8,
            label=r"inertio-capillary" if show_labels else None,
        )
    else:
        ic_prefactor = reference_run_normalisation(
            series_by_run, quantity, ic_slope, cone_fit_window
        )
        ax.plot(
            r_cone,
            ic_prefactor * r_cone**ic_slope,
            color=GREY,
            ls="--",
            lw=LINE["theory_linewidth"],
            zorder=8,
            label=r"inertio-capillary" if show_labels else None,
        )

    prf_prefactor = run_weighted_normalisation(
        series_by_run, quantity, prf_slope, PRF_FIT_WINDOW
    )
    ax.plot(
        r_prf,
        prf_prefactor * r_prf**prf_slope,
        color=GREY,
        ls=":",
        lw=LINE["theory_linewidth"] + 0.25,
        zorder=8,
        label=r"Gordillo \& Blanco-Rodr\'iguez 2023 [25]" if show_labels else None,
    )


def build_figure(
    data_dir: Path,
    output: Path,
    use_tex: bool = True,
    cone_fit_window: tuple[float, float] = CONE_FIT_WINDOW,
) -> None:
    configure_matplotlib(use_tex=use_tex)

    series_by_run: list[tuple[RunSpec, dict[str, np.ndarray]]] = []
    for run in RUNS:
        series_by_run.append((run, processed_series(read_log(data_dir / run.filename))))

    fig, axes = plt.subplots(1, 3, figsize=(APS_DOUBLE_COL, 2.46))
    fig.set_facecolor("white")

    panels = (
        ("Q_j", r"$Q_j$", (3.0 * ALPHA - 1.0) / ALPHA, 1.5, 1.0),
        ("q_j", r"$q_j$", (2.0 * ALPHA - 1.0) / ALPHA, 0.5, 0.0),
        ("We_j", r"$We_j$", (3.0 * ALPHA - 2.0) / ALPHA, 0.0, -1.0),
    )

    for ax, (quantity, ylabel, cone_slope, ic_slope, prf_slope), panel_label in zip(
        axes, panels, (r"(a)", r"(b)", r"(c)")
    ):
        ax.axvspan(*cone_fit_window, color=LIGHT_GREY, alpha=0.22, lw=0, zorder=0)
        draw_theory(
            ax,
            series_by_run=series_by_run,
            quantity=quantity,
            slopes=(cone_slope, ic_slope, prf_slope),
            show_labels=(quantity == "Q_j"),
            cone_fit_window=cone_fit_window,
        )

        for run, series in series_by_run:
            idx = log_bin_indices(series["r_j"], series[quantity], target=68)
            ax.plot(
                series["r_j"][idx],
                series[quantity][idx],
                linestyle="None",
                marker=run.marker,
                ms=LINE["markersize"],
                mfc=run.colour,
                mec="black",
                mew=LINE["markeredgewidth"],
                alpha=marker_alpha(run),
                label=run.label if quantity == "Q_j" else None,
                zorder=3,
            )

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(0.005, MAX_RJ)
        ax.set_xlabel(r"$r_j$", fontsize=APS["LabelFont"], labelpad=3)
        ax.set_ylabel(ylabel, fontsize=APS["LabelFont"], labelpad=2)
        style_axes(ax)
        add_panel_label(ax, panel_label)

    axes[0].set_ylim(0.01, 4.0)
    axes[1].set_ylim(0.35, 3.2)
    axes[2].set_ylim(0.06, 260.0)
    axes[1].yaxis.set_label_coords(-0.18, 0.5)

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="upper center",
        bbox_to_anchor=(0.53, 1.05),
        ncol=4,
        frameon=False,
        fontsize=APS["LegendFont"],
        handlelength=1.5,
        columnspacing=0.75,
        handletextpad=0.35,
    )

    fig.subplots_adjust(left=0.07, right=0.99, bottom=0.18, top=0.76, wspace=0.50)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight", pad_inches=0.035, dpi=300)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--output", type=Path, default=ROOT / "fig2_flux_scalings.pdf")
    parser.add_argument("--no-tex", action="store_true", help="Use mathtext instead of LaTeX.")
    parser.add_argument(
        "--cone-fit-window",
        type=float,
        nargs=2,
        metavar=("RMIN", "RMAX"),
        default=CONE_FIT_WINDOW,
        help="Near-inception fit window used for the cone and inertio-capillary prefactors.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_figure(
        args.data_dir,
        args.output,
        use_tex=not args.no_tex,
        cone_fit_window=tuple(args.cone_fit_window),
    )
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
