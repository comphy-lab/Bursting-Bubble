#!/usr/bin/env python3
"""Build the two-column Fig. 2 v2 layout.

Panel (a) shows the four-frame velocity/streamline diagnostic. Panel (b)
uses the approved Q_j processing from make_fig2_flux_scalings.py. Panel (c)
first constructs an intermediate Q_j branch and then evaluates
We_j = Q_j^2/(pi^2 r_j^3).
"""

from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from scipy.interpolate import PchipInterpolator

from capsule_utils import atomic_savefig


sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
APS_DOUBLE_COL = 6.75
FIG_HEIGHT = 2.62


def load_script_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


fig2a = load_script_module("fig2a_streamlines", SCRIPT_DIR / "make_fig2a_streamlines.py")
flux = load_script_module("fig2_flux_scalings", SCRIPT_DIR / "make_fig2_flux_scalings.py")

SHORT_LEGEND_LABELS = {
    rf"cone ($\alpha={flux.ALPHA:.3f}$)": rf"present theory, $\alpha={flux.ALPHA:.3f}$",
    r"inertio-capillary": r"inertio-capillary, $\alpha=2/3$",
    r"L13": r"Level 13, focus 13",
    r"L14": r"Level 14, focus 13",
    r"L15, focus 13": r"Level 15, focus 13",
    r"L15, focus 14": r"Level 15, focus 14",
    r"L15, focus 15": r"Level 15, focus 15",
}


def tune_flux_style() -> None:
    flux.APS.update(
        {
            "LabelFont": 8.8,
            "AxesFont": 7.5,
            "LegendFont": 6.2,
            "PanelFont": 9.4,
        }
    )
    flux.LINE.update(
        {
            "theory_linewidth": 1.15,
            "spine_width": 0.75,
            "tick_width": 0.6,
            "tick_length_major": 3.2,
            "tick_length_minor": 1.7,
            "markersize": 2.15,
            "markeredgewidth": 0.22,
        }
    )


def load_streamline_fields(args: argparse.Namespace):
    return fig2a.load_archived_inputs(args.fig2a_data_dir, tuple(args.snapshots))


def draw_panel_a(fig: plt.Figure, bbox: tuple[float, float, float, float], args: argparse.Namespace) -> None:
    fields, all_segments = load_streamline_fields(args)
    speed_norm = Normalize(vmin=0.0, vmax=args.vmax)
    cmap_speed = plt.get_cmap("Blues").copy()
    cmap_speed.set_bad((1, 1, 1, 0))

    left, bottom, width, height = bbox
    gap_x = 0.006
    gap_y = 0.004
    cbar_height = 0.018
    cbar_gap = 0.006
    frame_width = (width - gap_x) / 2.0
    frame_height = (height - cbar_height - cbar_gap - gap_y) / 2.0
    y_top = bottom + cbar_height + cbar_gap + frame_height + gap_y

    axes = [
        fig.add_axes([left, y_top, frame_width, frame_height]),
        fig.add_axes([left + frame_width + gap_x, y_top, frame_width, frame_height]),
        fig.add_axes([left, bottom + cbar_height + cbar_gap, frame_width, frame_height]),
        fig.add_axes(
            [
                left + frame_width + gap_x,
                bottom + cbar_height + cbar_gap,
                frame_width,
                frame_height,
            ]
        ),
    ]
    frame_labels = [r"(i)", r"(ii)", r"(iii)", r"(iv)"]
    for idx, (ax, field, segs, snap, label) in enumerate(
        zip(axes, fields, all_segments, args.snapshots, frame_labels)
    ):
        fig2a.draw_frame(
            ax,
            field,
            segs,
            snap,
            args.zmin,
            args.zmax,
            args.rmax,
            speed_norm,
            cmap_speed,
            label,
        )
        ax.set_anchor("S" if idx < 2 else "N")

    cbar_width = width * 0.86
    cbar_left = left + 0.5 * (width - cbar_width)
    cax = fig.add_axes([cbar_left, bottom, cbar_width, cbar_height])
    cb = fig.colorbar(
        ScalarMappable(norm=speed_norm, cmap=cmap_speed),
        cax=cax,
        orientation="horizontal",
    )
    cb.set_label(r"$|\mathbf{u}|$", fontsize=6.7, labelpad=0.7)
    cb.set_ticks([0, 25, 50])
    cb.ax.tick_params(labelsize=6.6, length=2.0, width=0.45, pad=0.7)
    cb.outline.set_linewidth(0.45)


def draw_flux_panel(
    ax: plt.Axes,
    series_by_run,
    quantity: str,
    ylabel: str,
    slopes: tuple[float, float, float],
    panel_label: str,
    cone_fit_window: tuple[float, float],
    show_labels: bool,
    resample_markers: bool = True,
) -> None:
    ax.axvspan(*cone_fit_window, color=flux.LIGHT_GREY, alpha=0.22, lw=0, zorder=0)
    draw_theory_v2(
        ax,
        series_by_run=series_by_run,
        quantity=quantity,
        slopes=slopes,
        show_labels=show_labels,
        cone_fit_window=cone_fit_window,
    )

    for run, series in series_by_run:
        if resample_markers:
            idx = flux.log_bin_indices(series["r_j"], series[quantity], target=52)
        else:
            idx = np.arange(len(series["r_j"]), dtype=int)
        ax.plot(
            series["r_j"][idx],
            series[quantity][idx],
            linestyle="None",
            marker=run.marker,
            ms=flux.LINE["markersize"],
            mfc=run.colour,
            mec="black",
            mew=flux.LINE["markeredgewidth"],
            alpha=flux.marker_alpha(run),
            label=run.label if show_labels else None,
            zorder=3,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.005, flux.MAX_RJ)
    ax.set_xlabel(r"$r_j$", fontsize=flux.APS["LabelFont"], labelpad=2)
    ax.set_ylabel(ylabel, fontsize=flux.APS["LabelFont"], labelpad=1.5)
    flux.style_axes(ax)


def draw_theory_v2(
    ax: plt.Axes,
    series_by_run,
    quantity: str,
    slopes: tuple[float, float, float],
    show_labels: bool,
    cone_fit_window: tuple[float, float],
) -> None:
    cone_slope, ic_slope, prf_slope = slopes
    r_cone = np.geomspace(*flux.CONE_DRAW_WINDOW, 120)
    r_prf = np.geomspace(*flux.PRF_DRAW_WINDOW, 120)

    cone_prefactor = flux.reference_run_normalisation(
        series_by_run, quantity, cone_slope, cone_fit_window
    )
    ax.plot(
        r_cone,
        cone_prefactor * r_cone**cone_slope,
        color=flux.BLACK,
        ls="-",
        lw=flux.LINE["theory_linewidth"],
        zorder=8,
        label=rf"cone ($\alpha={flux.ALPHA:.3f}$)" if show_labels else None,
    )

    if quantity == "We_j":
        ax.axhline(
            1.0,
            color=flux.GREY,
            ls="--",
            lw=flux.LINE["theory_linewidth"],
            zorder=8,
            label=r"inertio-capillary" if show_labels else None,
        )
    else:
        ic_prefactor = flux.reference_run_normalisation(
            series_by_run, quantity, ic_slope, cone_fit_window
        )
        ax.plot(
            r_cone,
            ic_prefactor * r_cone**ic_slope,
            color=flux.GREY,
            ls="--",
            lw=flux.LINE["theory_linewidth"],
            zorder=8,
            label=r"inertio-capillary" if show_labels else None,
        )

    prf_prefactor = flux.run_weighted_normalisation(
        series_by_run, quantity, prf_slope, flux.PRF_FIT_WINDOW
    )
    ax.plot(
        r_prf,
        prf_prefactor * r_prf**prf_slope,
        color=flux.BLACK,
        ls=":",
        lw=flux.LINE["theory_linewidth"] + 0.2,
        zorder=9,
        label=flux.literature_label() if show_labels else None,
    )


def _unique_log_points(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return strictly increasing log-space points for PCHIP."""
    order = np.argsort(x)
    x_sorted = x[order]
    y_sorted = y[order]
    unique_x: list[float] = []
    unique_y: list[float] = []
    start = 0
    while start < len(x_sorted):
        stop = start + 1
        while stop < len(x_sorted) and np.isclose(
            x_sorted[stop], x_sorted[start], rtol=0.0, atol=1e-13
        ):
            stop += 1
        unique_x.append(float(x_sorted[start]))
        unique_y.append(float(np.median(y_sorted[start:stop])))
        start = stop
    return np.asarray(unique_x), np.asarray(unique_y)


def constrained_q_interpolation(
    series: dict[str, np.ndarray],
    r_eval: np.ndarray,
    slope: float,
    fit_window: tuple[float, float],
    anchor_r: float,
    blend_start_r: float,
) -> tuple[np.ndarray, float]:
    """Build the intermediate Q_j branch used before evaluating We_j."""
    r = series["r_j"]
    q = series["Q_j"]
    fit = (
        np.isfinite(r)
        & np.isfinite(q)
        & (r >= fit_window[0])
        & (r <= fit_window[1])
        & (r > 0.0)
        & (q > 0.0)
    )
    if not np.any(fit):
        raise ValueError("No Q_j data in the requested fit window")

    prefactor = float(np.exp(np.mean(np.log(q[fit]) - slope * np.log(r[fit]))))
    r_eval = np.asarray(r_eval)

    valid = np.isfinite(r) & np.isfinite(q) & (r > 0.0) & (q > 0.0)
    if np.count_nonzero(valid) < 2:
        return prefactor * r_eval**slope, prefactor

    x = np.log(r[valid])
    y = np.log(q[valid])
    x, y = _unique_log_points(x, y)
    if len(x) < 2:
        return prefactor * r_eval**slope, prefactor

    if not (0.0 < blend_start_r < anchor_r):
        raise ValueError("--interp-blend-start-r must be positive and smaller than --interp-anchor-r")

    interpolant = PchipInterpolator(x, y, extrapolate=True)
    log_r_eval = np.log(r_eval)
    log_q_asymptote = np.log(prefactor) + slope * log_r_eval
    log_q_data = interpolant(log_r_eval)

    t = (log_r_eval - np.log(blend_start_r)) / (np.log(anchor_r) - np.log(blend_start_r))
    t = np.clip(t, 0.0, 1.0)
    data_weight = t * t * (3.0 - 2.0 * t)
    q_eval = np.exp((1.0 - data_weight) * log_q_asymptote + data_weight * log_q_data)
    return q_eval, prefactor


def build_interpolated_we_series(
    series_by_run,
    q_slope: float,
    fit_window: tuple[float, float],
    anchor_r: float,
    blend_start_r: float,
    marker_target: int,
):
    interpolated = []
    prefactors = []
    for run, series in series_by_run:
        idx = flux.log_bin_indices(series["r_j"], series["Q_j"], target=marker_target)
        r_eval = series["r_j"][idx]
        q_eval, prefactor = constrained_q_interpolation(
            series,
            r_eval,
            q_slope,
            fit_window,
            anchor_r,
            blend_start_r,
        )
        we_eval = q_eval**2 / (np.pi**2 * r_eval**3)
        interpolated.append(
            (
                run,
                {
                    "r_j": r_eval,
                    "Q_j": q_eval,
                    "q_j": q_eval / (np.pi * r_eval),
                    "We_j": we_eval,
                },
            )
        )
        prefactors.append((run, prefactor))
    return interpolated, prefactors


def build_figure(args: argparse.Namespace) -> None:
    flux.configure_matplotlib(use_tex=not args.no_tex)
    tune_flux_style()

    series_by_run = [
        (run, flux.processed_series(flux.read_log(args.data_dir / run.filename)))
        for run in flux.RUNS
    ]

    fig = plt.figure(figsize=(APS_DOUBLE_COL, FIG_HEIGHT))
    fig.set_facecolor("white")

    draw_panel_a(fig, (0.014, 0.105, 0.365, 0.790), args)

    legend_top = 0.160
    theory_legend_ax = fig.add_axes([0.405, 0.020, 0.342, legend_top - 0.020])
    theory_legend_ax.axis("off")
    symbol_left_legend_ax = fig.add_axes([0.658, 0.098, 0.172, legend_top - 0.098])
    symbol_left_legend_ax.axis("off")
    symbol_right_legend_ax = fig.add_axes([0.828, 0.062, 0.167, legend_top - 0.062])
    symbol_right_legend_ax.axis("off")
    ax_b = fig.add_axes([0.430, 0.325, 0.258, 0.600])
    ax_c = fig.add_axes([0.735, 0.325, 0.255, 0.600])

    q_slope = (3.0 * flux.ALPHA - 1.0) / flux.ALPHA
    we_slope = (3.0 * flux.ALPHA - 2.0) / flux.ALPHA
    we_series_by_run, q_prefactors = build_interpolated_we_series(
        series_by_run,
        q_slope,
        tuple(args.cone_fit_window),
        args.interp_anchor_r,
        args.interp_blend_start_r,
        args.interp_marker_target,
    )
    resample_we_markers = False
    print(
        "Q_j interpolation prefactors for "
        f"Q_j ~ A*r_j^{q_slope:.6g} as r_j -> 0, "
        f"blended to PCHIP data over "
        f"{args.interp_blend_start_r:g} <= r_j <= {args.interp_anchor_r:g}:"
    )
    for run, prefactor in q_prefactors:
        print(f"  {SHORT_LEGEND_LABELS.get(run.label, run.label)}: A = {prefactor:.6g}")

    draw_flux_panel(
        ax_b,
        series_by_run,
        "Q_j",
        r"$Q_j$",
        (q_slope, 1.5, 1.0),
        r"(b)",
        tuple(args.cone_fit_window),
        show_labels=True,
    )
    draw_flux_panel(
        ax_c,
        we_series_by_run,
        "We_j",
        r"$We_j$",
        (we_slope, 0.0, -1.0),
        r"(c)",
        tuple(args.cone_fit_window),
        show_labels=False,
        resample_markers=resample_we_markers,
    )

    ax_b.set_ylim(0.01, 4.0)
    ax_c.set_ylim(0.75, 260.0)
    ax_b.yaxis.set_label_coords(-0.155, 0.5)
    ax_c.yaxis.set_label_coords(-0.190, 0.5)

    fig.text(0.006, 0.875, r"(a)", ha="left", va="bottom",
             fontsize=flux.APS["PanelFont"], fontweight="bold")
    fig.text(0.395, 0.875, r"(b)", ha="left", va="bottom",
             fontsize=flux.APS["PanelFont"], fontweight="bold")
    fig.text(0.690, 0.875, r"(c)", ha="left", va="bottom",
             fontsize=flux.APS["PanelFont"], fontweight="bold")

    handles, labels = ax_b.get_legend_handles_labels()
    labels = [SHORT_LEGEND_LABELS.get(label, label) for label in labels]
    theory_legend_ax.legend(
        handles[:3],
        labels[:3],
        loc="upper left",
        bbox_to_anchor=(0.0, 1.0),
        ncol=1,
        frameon=False,
        fontsize=flux.APS["LegendFont"],
        handlelength=1.10,
        handletextpad=0.22,
        labelspacing=0.32,
        borderaxespad=0.0,
    )
    symbol_left_legend_ax.legend(
        handles[3:5],
        labels[3:5],
        loc="upper left",
        bbox_to_anchor=(0.0, 1.0),
        ncol=1,
        frameon=False,
        fontsize=flux.APS["LegendFont"],
        handlelength=0.55,
        handletextpad=0.12,
        labelspacing=0.30,
        borderaxespad=0.0,
    )
    symbol_right_legend_ax.legend(
        handles[5:],
        labels[5:],
        loc="upper left",
        bbox_to_anchor=(0.0, 1.0),
        ncol=1,
        frameon=False,
        fontsize=flux.APS["LegendFont"],
        handlelength=0.55,
        handletextpad=0.12,
        labelspacing=0.30,
        borderaxespad=0.0,
    )

    atomic_savefig(fig, args.output, dpi=300)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=flux.DEFAULT_DATA_DIR)
    parser.add_argument("--output", type=Path, default=ROOT / "fig2_v2.pdf")
    parser.add_argument("--fig2a-data-dir", type=Path, default=fig2a.DEFAULT_DATA_DIR)
    parser.add_argument("--snapshots", nargs="+", default=list(fig2a.DEFAULT_SNAPSHOTS))
    parser.add_argument("--zmin", type=float, default=-1.72)
    parser.add_argument("--zmax", type=float, default=-0.82)
    parser.add_argument("--rmax", type=float, default=0.58)
    parser.add_argument("--nr", type=int, default=190)
    parser.add_argument("--vmax", type=float, default=fig2a.DEFAULT_VMAX)
    parser.add_argument(
        "--we-from-q-interp",
        action="store_true",
        help=argparse.SUPPRESS,
    )
    parser.add_argument(
        "--interp-anchor-r",
        type=float,
        default=0.1,
        help="Radius by which the Q_j interpolation has fully returned to the PCHIP data branch.",
    )
    parser.add_argument(
        "--interp-blend-start-r",
        type=float,
        default=0.005,
        help="Radius below which the Q_j interpolation uses the asymptotic A*r_j^1.41 branch.",
    )
    parser.add_argument(
        "--interp-marker-target",
        type=int,
        default=52,
        help="Log-binned Q_j marker count used before sampling interpolated We_j.",
    )
    parser.add_argument(
        "--cone-fit-window",
        type=float,
        nargs=2,
        metavar=("RMIN", "RMAX"),
        default=flux.CONE_FIT_WINDOW,
    )
    parser.add_argument("--no-tex", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_figure(args)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
