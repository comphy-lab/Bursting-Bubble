#!/usr/bin/env python3
"""Build Reynolds/Weber self-consistency and cutoff-resolution diagnostics.

The self-similar figure derives ``Re_j = q_j/Oh`` from the same processed
jet-base series used by Fig. 2 and compares it with ``We_j``.  The cutoff
figure deliberately remains separate: it evaluates ``R_kappa``, ``Re_m`` and
``We_m`` at the earliest available post-inception sample and labels those
values as grid-censored whenever ``R_kappa`` is of order the local cell size.
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import asdict, dataclass
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, NullFormatter
import numpy as np

import make_fig2_flux_scalings as flux


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OH = 0.03
TAU_MIN = 1.0e-6
PLOT_R_MAX = 0.06


@dataclass(frozen=True)
class CutoffEstimate:
    """One earliest-post-inception, resolution-censored cutoff estimate."""

    label: str
    level: int
    source: str
    tau: float
    delta: float
    curvature_radius: float
    curvature_radius_cells: float
    u_z_tip: float
    re_m: float
    we_m: float
    tip_cell_offset_cells: float


def similarity_exponents(alpha: float) -> tuple[float, float]:
    """Return the fixed-shape exponents for ``Re_j(r_j)`` and ``We_j(r_j)``."""
    if not math.isfinite(alpha) or alpha <= 0.5:
        raise ValueError("alpha must be finite and greater than 1/2")
    return (2.0 * alpha - 1.0) / alpha, (3.0 * alpha - 2.0) / alpha


def load_jet_series(data_dir: Path, oh: float) -> list[tuple[flux.RunSpec, dict[str, np.ndarray]]]:
    """Load the Fig. 2 jet-base series and derive the local Reynolds number."""
    if not math.isfinite(oh) or oh <= 0.0:
        raise ValueError("Oh must be finite and positive")
    output: list[tuple[flux.RunSpec, dict[str, np.ndarray]]] = []
    for run in flux.RUNS:
        series = flux.processed_series(flux.read_log(data_dir / run.filename))
        series["Re_j"] = series["q_j"] / oh
        output.append((run, series))
    return output


def cutoff_estimate(
    label: str,
    level: int,
    t0: float,
    path: Path,
    oh: float,
    tau_min: float = TAU_MIN,
) -> CutoffEstimate:
    """Evaluate the first finite tip sample at or after ``tau_min``."""
    if level <= 0 or not math.isfinite(t0):
        raise ValueError("level must be positive and t0 finite")
    if not math.isfinite(oh) or oh <= 0.0 or tau_min < 0.0:
        raise ValueError("Oh must be positive and tau_min nonnegative")
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    required = {
        "time", "inverse_mean_curvature", "u_z_tip", "delta_tip",
        "tip_cell_offset_cells",
    }
    if not rows or not required.issubset(rows[0]):
        raise ValueError(f"Missing cutoff columns in {path}")

    candidates: list[tuple[float, float, float, float, float]] = []
    for row in rows:
        values = tuple(float(row[name]) for name in required)
        if not all(math.isfinite(value) for value in values):
            continue
        tau = float(row["time"]) - t0
        radius = float(row["inverse_mean_curvature"])
        delta = float(row["delta_tip"])
        if tau >= tau_min and radius > 0.0 and delta > 0.0:
            candidates.append(
                (
                    tau,
                    radius,
                    abs(float(row["u_z_tip"])),
                    delta,
                    float(row["tip_cell_offset_cells"]),
                )
            )
    if not candidates:
        raise ValueError(f"No finite post-inception cutoff sample in {path}")
    tau, radius, velocity, delta, offset = min(candidates, key=lambda item: item[0])
    return CutoffEstimate(
        label=label,
        level=level,
        source=str(path.resolve()),
        tau=tau,
        delta=delta,
        curvature_radius=radius,
        curvature_radius_cells=radius / delta,
        u_z_tip=velocity,
        re_m=radius * velocity / oh,
        we_m=radius * velocity**2,
        tip_cell_offset_cells=offset,
    )


def atomic_json(path: Path, value: object) -> None:
    """Write a JSON receipt atomically."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def build_selfsimilar_figure(
    series_by_run: list[tuple[flux.RunSpec, dict[str, np.ndarray]]],
    oh: float,
    output: Path,
) -> dict[str, object]:
    """Plot resolved jet-base Reynolds and Weber numbers without connected DNS lines."""
    flux.configure_matplotlib(use_tex=True)
    re_slope, we_slope = similarity_exponents(flux.ALPHA)
    re_prefactor = flux.reference_run_normalisation(
        series_by_run, "Re_j", re_slope, flux.CONE_FIT_WINDOW
    )
    we_prefactor = flux.reference_run_normalisation(
        series_by_run, "We_j", we_slope, flux.CONE_FIT_WINDOW
    )
    radius = np.geomspace(flux.CONE_FIT_WINDOW[0], flux.CONE_FIT_WINDOW[1], 240)

    fig, axes = plt.subplots(1, 2, figsize=(6.75, 3.05))
    data_handles: list[Line2D] = []
    ranges: dict[str, dict[str, float]] = {}
    for run, series in series_by_run:
        window = (series["r_j"] >= flux.CONE_FIT_WINDOW[0]) & (series["r_j"] <= PLOT_R_MAX)
        indices = flux.log_bin_indices(series["r_j"][window], series["Re_j"][window], target=55)
        r_plot = series["r_j"][window][indices]
        for axis, quantity in zip(axes, ("Re_j", "We_j"), strict=True):
            axis.plot(
                r_plot,
                series[quantity][window][indices],
                linestyle="None",
                marker=run.marker,
                ms=flux.LINE["markersize"] + 0.25,
                mfc=run.colour,
                mec="black",
                mew=flux.LINE["markeredgewidth"],
                alpha=flux.marker_alpha(run),
                zorder=3,
            )
        data_handles.append(
            Line2D([], [], linestyle="None", marker=run.marker, ms=4.2,
                   mfc=run.colour, mec="black", mew=0.35, label=run.label)
        )
        fit = (series["r_j"] >= flux.CONE_FIT_WINDOW[0]) & (
            series["r_j"] <= flux.CONE_FIT_WINDOW[1]
        )
        ranges[run.label] = {
            "re_min": float(np.min(series["Re_j"][fit])),
            "re_max": float(np.max(series["Re_j"][fit])),
            "we_min": float(np.min(series["We_j"][fit])),
            "we_max": float(np.max(series["We_j"][fit])),
        }

    for axis in axes:
        axis.axvspan(*flux.CONE_FIT_WINDOW, color=flux.LIGHT_GREY, alpha=0.30, lw=0)
        axis.axhline(1.0, color="0.45", ls="--", lw=1.0, zorder=1)
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_xlim(0.0045, PLOT_R_MAX)
        axis.set_xlabel(r"$r_j$")
        flux.style_axes(axis)
        axis.set_xticks((0.005, 0.01, 0.02, 0.04), ("0.005", "0.01", "0.02", "0.04"))
        axis.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10)*0.1))
        axis.xaxis.set_minor_formatter(NullFormatter())
    axes[0].plot(radius, re_prefactor * radius**re_slope, color="black", lw=1.5)
    axes[1].plot(radius, we_prefactor * radius**we_slope, color="black", lw=1.5)
    axes[0].set_ylabel(r"$Re_j=q_j/Oh$")
    axes[1].set_ylabel(r"$We_j=q_j^2/r_j$")
    axes[0].set_ylim(0.7, 80.0)
    axes[1].set_ylim(0.7, 300.0)
    flux.add_panel_label(axes[0], r"(a)")
    flux.add_panel_label(axes[1], r"(b)")
    for axis, reference_label in zip(axes, (r"$Re_j=1$", r"$We_j=1$"), strict=True):
        axis.legend(
            handles=(
                Line2D([], [], color="black", lw=1.5, label="present similarity"),
                Line2D([], [], color="0.45", ls="--", lw=1.0, label=reference_label),
            ),
            frameon=False,
            fontsize=7.5,
            loc="lower right",
        )
    fig.legend(handles=data_handles, loc="upper center", bbox_to_anchor=(0.52, 1.02),
               ncol=6, frameon=False, fontsize=7.5, handletextpad=0.25,
               columnspacing=0.65)
    fig.text(0.505, 0.01, rf"$Oh={oh:g}$; shaded: main asymptotic window",
             ha="center", va="bottom", fontsize=8)
    fig.subplots_adjust(left=0.10, right=0.99, bottom=0.20, top=0.82, wspace=0.31)
    flux.atomic_savefig(fig, output, dpi=300, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return {
        "oh": oh,
        "alpha": flux.ALPHA,
        "fit_window": list(flux.CONE_FIT_WINDOW),
        "re_slope": re_slope,
        "we_slope": we_slope,
        "re_prefactor": re_prefactor,
        "we_prefactor": we_prefactor,
        "run_ranges_in_fit_window": ranges,
    }


def build_cutoff_figure(estimates: Sequence[CutoffEstimate], output: Path) -> None:
    """Plot the available cutoff estimates as unconnected, grid-censored markers."""
    flux.configure_matplotlib(use_tex=True)
    ordered = sorted(estimates, key=lambda item: item.level)
    levels = np.asarray([item.level for item in ordered], dtype=float)
    colours = plt.get_cmap("viridis")(np.linspace(0.12, 0.88, len(ordered)))
    values = (
        ("curvature_radius_cells", r"$R_\kappa/\Delta_{\mathrm{tip}}$", (0.0, 1.5)),
        ("re_m", r"$Re_m=R_\kappa u_{z,\mathrm{tip}}/Oh$", (0.0, 6.0)),
        ("we_m", r"$We_m=u_{z,\mathrm{tip}}^2R_\kappa$", (0.0, 48.0)),
    )
    fig, axes = plt.subplots(1, 3, figsize=(6.75, 2.75))
    for axis, (field, ylabel, ylim), panel in zip(axes, values, ("(a)", "(b)", "(c)"), strict=True):
        axis.axhline(1.0, color="0.45", ls="--", lw=1.0)
        for level, colour, estimate in zip(levels, colours, ordered, strict=True):
            value = float(getattr(estimate, field))
            axis.plot(level, value, linestyle="None", marker="o", ms=6.0,
                      mfc="white", mec=colour, mew=1.5)
            axis.annotate(f"{value:.2f}", (level, value), xytext=(0, 6),
                          textcoords="offset points", ha="center", va="bottom", fontsize=7.5)
        axis.set_xlim(levels.min() - 0.6, levels.max() + 0.6)
        axis.set_ylim(*ylim)
        axis.set_xticks(levels, [f"L{int(level)}" for level in levels])
        axis.set_xlabel("actual tip level")
        axis.set_ylabel(ylabel)
        flux.style_axes(axis)
        flux.add_panel_label(axis, panel)
    fig.text(0.51, 0.01,
             "Earliest available post-inception sample; open symbols are grid-censored",
             ha="center", va="bottom", fontsize=8)
    fig.subplots_adjust(left=0.075, right=0.995, bottom=0.24, top=0.88, wspace=0.52)
    flux.atomic_savefig(fig, output, dpi=300, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=flux.DEFAULT_DATA_DIR)
    parser.add_argument("--oh", type=float, default=DEFAULT_OH)
    parser.add_argument(
        "--cutoff-series", action="append", nargs=4,
        metavar=("LABEL", "LEVEL", "T0", "CSV"), required=True,
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Build both figures and their machine-readable evidence receipt."""
    args = build_parser().parse_args(argv)
    output_dir = args.output_dir.resolve()
    series = load_jet_series(args.data_dir.resolve(), args.oh)
    estimates = [
        cutoff_estimate(label, int(level), float(t0), Path(path).resolve(), args.oh)
        for label, level, t0, path in args.cutoff_series
    ]
    selfsimilar = build_selfsimilar_figure(
        series, args.oh, output_dir / "re-we-self-similar.pdf"
    )
    build_cutoff_figure(estimates, output_dir / "cutoff-resolution-diagnostic.pdf")
    atomic_json(
        output_dir / "re-we-consistency.json",
        {
            "schema_version": 1,
            "claim": {
                "self_similar": "existing jet-base DNS supports Re_j >> 1 in the main fit window",
                "cutoff": "grid-censored diagnostic; Re_m approximately one is not established",
            },
            "self_similar": selfsimilar,
            "cutoff_method": "earliest finite tip sample with tau >= 1e-6",
            "cutoff_estimates": [asdict(estimate) for estimate in estimates],
        },
    )
    print(f"RE_WE_COMPLETE output={output_dir} cutoff_series={len(estimates)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
