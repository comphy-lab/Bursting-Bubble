#!/usr/bin/env python3
"""Build the two-column Fig. 2 v2 layout.

Panel (a) shows the four-frame velocity/streamline diagnostic. Panel (b)
uses the approved Q_j processing from make_fig2_flux_scalings.py. Panel (c)
first constructs an intermediate Q_j branch and then evaluates
We_j = Q_j^2/(pi^2 r_j^3), with the approved empirical base-grid ceiling
applied to its extrapolated r_j-to-zero limit.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from scipy.interpolate import PchipInterpolator

from capsule_utils import atomic_savefig


sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
APS_DOUBLE_COL = 6.75
FIG_HEIGHT = 2.78


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
GRID_WE_METADATA = flux.FIG2_METADATA["grid_informed_weber"]
GRID_WE_CASES = {
    str(item["case"]): item for item in GRID_WE_METADATA["cases"]
}

SHORT_LEGEND_LABELS = {
    rf"cone ($\alpha={flux.ALPHA:.3f}$)": rf"present theory, $\alpha={flux.ALPHA:.3f}$",
    r"inertio-capillary": r"inertio-capillary, $\alpha=2/3$",
    r"(12,13)": r"$(12,13)$",
    r"(12,14)": r"$(12,14)$",
    r"(13,14)": r"$(13,14)$",
    r"(13,15)": r"$(13,15)$",
    r"(14,15)": r"$(14,15)$",
    r"(15,15)": r"$(15,15)$",
    r"(14,16)": r"$(14,16)$",
    r"(15,16)": r"$(15,16)$",
    r"(16\to14,16)": r"$(16\!\to\!14,16)$",
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


def case_id(run) -> str:
    """Return the four-digit case identifier encoded in a run filename."""
    return run.filename.split("_", 1)[0]


def grid_informed_runs():
    """Return the eight unbridged runs used by the empirical grid correlation."""
    excluded = {str(case) for case in GRID_WE_METADATA["excluded_cases"]}
    runs = tuple(run for run in flux.RUNS if case_id(run) not in excluded)
    if {case_id(run) for run in runs} != set(GRID_WE_CASES):
        raise RuntimeError("Figure 2 grid-informed case set differs from metadata")
    return runs


def grid_weber_ceiling(level: int) -> float:
    """Evaluate the fixed eight-case correlation without refitting it."""
    correlation = GRID_WE_METADATA["correlation"]
    delta = float(correlation["domain_size_over_R"]) / 2**level
    return float(correlation["prefactor"] * delta ** correlation["exponent"])


def cap_flux(radius, volume_flux, weber_ceiling):
    """Apply the empirical Weber ceiling through the equivalent flux bound."""
    radius = np.asarray(radius, dtype=float)
    volume_flux = np.asarray(volume_flux, dtype=float)
    if np.any(~np.isfinite(radius)) or np.any(radius <= 0.0):
        raise ValueError("Grid-informed Weber processing requires positive finite radii")
    if np.any(~np.isfinite(volume_flux)) or np.any(volume_flux < 0.0):
        raise ValueError("Grid-informed Weber processing requires finite outward flux")
    if not np.isfinite(weber_ceiling) or weber_ceiling <= 0.0:
        raise ValueError("Grid-informed Weber ceiling must be positive and finite")
    return np.minimum(
        volume_flux, np.pi * radius**1.5 * np.sqrt(weber_ceiling)
    )


def apply_grid_informed_weber_caps(
    series_by_run,
    prefactors,
    q_slope: float,
    blend_start_r: float,
    data_dir: Path = flux.DEFAULT_DATA_DIR,
):
    """Cap the existing constrained-Q series and verify its hidden r->0 limit."""
    adjusted = []
    diagnostics = []
    we_slope = 2.0 * q_slope - 3.0
    display_minimum = float(GRID_WE_METADATA["display"]["x_limits"][0])
    for (run, series), (same_run, q_prefactor) in zip(
        series_by_run, prefactors, strict=True
    ):
        if run != same_run:
            raise RuntimeError("Q-interpolation prefactor order differs from run order")
        case = case_id(run)
        evidence = GRID_WE_CASES[case]
        source = data_dir / run.filename
        if hashlib.sha256(source.read_bytes()).hexdigest() != evidence["sha256"]:
            raise RuntimeError(f"Figure 2 source hash differs for case {case}")
        level = int(evidence["base_level_at_native_peak"])
        maximum_weber = grid_weber_ceiling(level)
        capped_q = cap_flux(series["r_j"], series["Q_j"], maximum_weber)
        capped_weber = np.minimum(series["We_j"], maximum_weber)
        capped_line_flux = capped_q / (np.pi * series["r_j"])
        if not np.allclose(
            capped_weber,
            capped_q**2 / (np.pi**2 * series["r_j"] ** 3),
            rtol=1.0e-12,
        ):
            raise RuntimeError(f"Flux and Weber caps disagree for case {case}")
        changed = capped_weber < series["We_j"]
        if np.any(changed & (series["r_j"] >= display_minimum)):
            raise RuntimeError("Empirical grid cap became active in the displayed range")

        we_prefactor = q_prefactor**2 / np.pi**2
        saturation_radius = float((maximum_weber / we_prefactor) ** (1.0 / we_slope))
        if not 1.0e-6 < saturation_radius < blend_start_r:
            raise RuntimeError(f"Unexpected saturation radius for case {case}")
        for test_radius in (1.0e-8, 1.0e-12):
            original_q = q_prefactor * test_radius**q_slope
            limited_q = cap_flux([test_radius], [original_q], maximum_weber)[0]
            if not np.isclose(
                limited_q**2 / (np.pi**2 * test_radius**3),
                maximum_weber,
                rtol=1.0e-12,
            ):
                raise RuntimeError(f"r_j-to-zero limit check failed for case {case}")
        adjusted.append(
            (
                run,
                {
                    "r_j": series["r_j"].copy(),
                    "Q_j": capped_q,
                    "q_j": capped_line_flux,
                    "We_j": capped_weber,
                },
            )
        )
        diagnostics.append(
            {
                "case": case,
                "level": level,
                "weber_ceiling": maximum_weber,
                "saturation_radius": saturation_radius,
                "displayed_values_changed": int(
                    np.count_nonzero(changed & (series["r_j"] >= display_minimum))
                ),
            }
        )
    return adjusted, diagnostics


def grid_legend_runs(runs):
    """Order the legend by measured peak-base level, then numerical pair."""
    return sorted(
        runs,
        key=lambda run: (
            int(GRID_WE_CASES[case_id(run)]["base_level_at_native_peak"]),
            run.focus,
            run.level,
            int(case_id(run)),
        ),
    )


def grid_legend_label(run) -> str:
    """Append the independently established native-peak base level."""
    level = int(GRID_WE_CASES[case_id(run)]["base_level_at_native_peak"])
    return rf"${run.label}_{{{level}}}$"


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
    marker_target: int = 36,
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
            idx = flux.log_bin_indices(
                series["r_j"], series[quantity], target=marker_target
            )
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

    runs = grid_informed_runs()
    series_by_run = [
        (run, flux.processed_series(flux.read_log(args.data_dir / run.filename)))
        for run in runs
    ]

    fig = plt.figure(figsize=(APS_DOUBLE_COL, FIG_HEIGHT))
    fig.set_facecolor("white")

    draw_panel_a(fig, (0.014, 0.105, 0.365, 0.790), args)

    legend_top = 0.195
    theory_legend_ax = fig.add_axes([0.405, 0.025, 0.310, legend_top - 0.025])
    theory_legend_ax.axis("off")
    symbol_legend_ax = fig.add_axes([0.735, 0.012, 0.263, 0.185])
    symbol_legend_ax.axis("off")
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
    we_series_by_run, grid_diagnostics = apply_grid_informed_weber_caps(
        we_series_by_run,
        q_prefactors,
        q_slope,
        args.interp_blend_start_r,
        data_dir=args.data_dir,
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
    print("Grid-informed Weber limits (inactive in the displayed range):")
    for item in grid_diagnostics:
        print(
            f"  case {item['case']}: L_B={item['level']} "
            f"We_max={item['weber_ceiling']:.6g} "
            f"r_sat={item['saturation_radius']:.6g}"
        )

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
    ax_c.set_xlim(*GRID_WE_METADATA["display"]["x_limits"])
    ax_c.set_ylim(*GRID_WE_METADATA["display"]["y_limits"])
    ax_c.set_xticks(
        [0.005, 0.01, 0.1, 1.0],
        [r"$0.005$", r"$10^{-2}$", r"$10^{-1}$", r"$10^0$"],
    )
    ax_c.get_xticklabels()[0].set_fontsize(6.0)
    ax_c.get_xticklabels()[0].set_ha("right")
    ax_b.yaxis.set_label_coords(-0.155, 0.5)
    ax_c.yaxis.set_label_coords(-0.153, 0.5)

    grid_information = [r"$L_B:\ \Delta_b/R\ \mapsto\ We_{\max}$"]
    for level in (13, 14, 15):
        delta = GRID_WE_METADATA["correlation"]["domain_size_over_R"] / 2**level
        grid_information.append(
            rf"${level}:\ {delta/1e-3:.3f}\!\times\!10^{{-3}}\ \mapsto\ "
            rf"{grid_weber_ceiling(level):.0f}$"
        )
    grid_information.extend(
        [r"Grid caps; inactive in this range", r"$L_B$: base level at native peak"]
    )
    ax_c.text(
        0.045,
        0.40,
        "\n".join(grid_information),
        transform=ax_c.transAxes,
        fontsize=4.7,
        ha="left",
        va="top",
        linespacing=1.35,
        bbox={"facecolor": "none", "edgecolor": "none", "pad": 0},
    )

    fig.text(0.006, 0.875, r"(a)", ha="left", va="bottom",
             fontsize=flux.APS["PanelFont"], fontweight="bold")
    fig.text(0.395, 0.875, r"(b)", ha="left", va="bottom",
             fontsize=flux.APS["PanelFont"], fontweight="bold")
    fig.text(0.690, 0.875, r"(c)", ha="left", va="bottom",
             fontsize=flux.APS["PanelFont"], fontweight="bold")

    handles, raw_labels = ax_b.get_legend_handles_labels()
    labels = [SHORT_LEGEND_LABELS.get(label, label) for label in raw_labels]
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
    handle_by_label = dict(zip(raw_labels[3:], handles[3:]))
    ordered_runs = grid_legend_runs(runs)
    groups = [
        [
            run
            for run in ordered_runs
            if GRID_WE_CASES[case_id(run)]["base_level_at_native_peak"] == level
        ]
        for level in (13, 14, 15)
    ]
    if [len(group) for group in groups] != [3, 2, 3]:
        raise RuntimeError("Grid-informed legend groups must have sizes 3, 2, 3")
    legend_slots = [group + [None] * (3 - len(group)) for group in groups]
    column_slots = [
        legend_slots[row][column] for column in range(3) for row in range(3)
    ]
    blank = Line2D([], [], linestyle="none", marker="None", color="none")
    symbol_legend_ax.legend(
        [handle_by_label[run.label] if run is not None else blank for run in column_slots],
        [grid_legend_label(run) if run is not None else "" for run in column_slots],
        loc="upper left",
        bbox_to_anchor=(0.0, 1.0),
        ncol=3,
        frameon=False,
        fontsize=5.7,
        handlelength=0.45,
        handletextpad=0.08,
        columnspacing=0.35,
        labelspacing=0.20,
        borderaxespad=0.0,
        title=r"Levels $(L_{\rm pre},L_{\rm post})_{L_B}$",
        title_fontsize=6.2,
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
        default=36,
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
