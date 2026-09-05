#!/usr/bin/env python3
r"""
# Diagnose tip-curvature and curvature-Weber scaling

Compare extracted tip metrics across $Oh$ and refinement. The script first
checks whether the inverse-curvature radius spans a configurable number of local
cells. It then reports the minimum radius inside a common post-inception time
window, rejects it from model comparison when grid- or window-boundary-limited,
and compares two fixed-shape hypotheses with fitted prefactors:

$$
R_m/\ell_\mu \propto Oh^{(2-3\alpha)/(2\alpha-1)},\qquad
We_\kappa \propto Oh^{(3\alpha-2)/(2\alpha-1)},
$$

and the inertio-capillary alternatives, which are constant in $Oh$. If
$\alpha=\alpha[\beta(Oh)]$ varies between series, the conical prediction is
evaluated pointwise rather than represented as one global power law.

This is a diagnostic comparison. It does not establish that the measured
$R_\kappa=1/|\kappa|$ is the theoretical cutoff radius $R_m$.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Sequence

import numpy as np


ONLINE_COLUMNS = (
    "i", "dt", "time", "jet_formed", "tip_pinched", "liquid_components",
    "tip_status", "z_tip", "r_tip", "z_cell", "r_cell", "kappa_mean",
    "u_z_tip", "u_r_tip", "speed_tip", "delta_tip", "level_tip", "f_tip",
    "r_base", "z_base", "q_jet", "q_l",
)


def cutoff_exponents(alpha: float) -> tuple[float, float, float]:
    """Return exponents for `r_m`, `R_m/ell_mu`, and `We_kappa`."""
    if not math.isfinite(alpha) or alpha <= 0.5:
        raise ValueError("alpha must be finite and greater than 1/2")
    radius = alpha / (2.0 * alpha - 1.0)
    normalised_radius = radius - 2.0
    weber = -normalised_radius
    return radius, normalised_radius, weber


def read_online_metrics(path: Path) -> dict[str, np.ndarray]:
    """Read valid, connected, pre-pinch rows from `tip_metrics.log`."""
    raw: list[dict[str, float]] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        if len(fields) != len(ONLINE_COLUMNS):
            raise ValueError(
                f"Expected {len(ONLINE_COLUMNS)} tip-log columns in {path}, found {len(fields)}"
            )
        values = np.asarray([float(field) for field in fields])
        if not np.all(np.isfinite(values)):
            raise ValueError(f"Non-finite tip-log row in {path}")
        row = dict(zip(ONLINE_COLUMNS, values, strict=True))
        if row["jet_formed"] == 1 and row["tip_pinched"] == 0 and row["tip_status"] == 3:
            raw.append(row)
    if not raw:
        raise ValueError(f"No valid connected pre-pinch tip rows in {path}")
    output = {name: np.asarray([row[name] for row in raw]) for name in ONLINE_COLUMNS}
    kappa = np.abs(output["kappa_mean"])
    output["curvature_radius"] = 1.0 / kappa
    output["curvature_radius_cells"] = output["curvature_radius"] / output["delta_tip"]
    output["we_curvature_uz"] = output["u_z_tip"]**2 * output["curvature_radius"]
    output["we_curvature_speed"] = output["speed_tip"]**2 * output["curvature_radius"]
    output["tip_cell_offset_cells"] = np.hypot(
        output["z_cell"] - output["z_tip"], output["r_cell"] - output["r_tip"]
    ) / output["delta_tip"]
    return output


def read_metrics(path: Path) -> dict[str, np.ndarray]:
    """Read finite metrics from an extractor CSV or online sidecar log."""
    with path.open(encoding="utf-8") as probe:
        first = probe.readline()
    if first.strip() == "# tip-metrics-v1":
        output = read_online_metrics(path)
        required = {
            "time", "z_tip", "curvature_radius", "curvature_radius_cells",
            "we_curvature_uz", "we_curvature_speed", "u_z_tip", "delta_tip",
            "tip_cell_offset_cells",
        }
        return {name: output[name] for name in required}

    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No rows in {path}")
    required = {
        "time",
        "z_tip",
        "inverse_mean_curvature",
        "speed_tip",
        "u_z_tip",
        "delta_tip",
        "tip_cell_offset_cells",
    }
    missing = required - set(rows[0])
    if missing:
        raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
    output: dict[str, np.ndarray] = {}
    for name in required:
        values = np.asarray([float(row[name]) for row in rows], dtype=float)
        if not np.all(np.isfinite(values)):
            raise ValueError(f"Non-finite {name} in {path}")
        output[name] = values
    output["curvature_radius"] = output["inverse_mean_curvature"]
    output["curvature_radius_cells"] = output["curvature_radius"] / output["delta_tip"]
    output["we_curvature_uz"] = output["u_z_tip"]**2 * output["curvature_radius"]
    output["we_curvature_speed"] = output["speed_tip"]**2 * output["curvature_radius"]
    order = np.argsort(output["time"])
    if np.any(np.diff(output["time"][order]) <= 0.0):
        raise ValueError(f"Times are not unique and increasing in {path}")
    return {name: values[order] for name, values in output.items()}


def summarise_series(
    label: str,
    oh: float,
    level: int,
    t0: float,
    alpha: float,
    path: Path,
    tau_window: tuple[float, float],
    min_cells: float,
    max_tip_offset_cells: float,
) -> tuple[dict[str, object], dict[str, np.ndarray]]:
    """Summarise one series and retain arrays for plotting."""
    if not math.isfinite(oh) or oh <= 0.0:
        raise ValueError("Oh must be finite and positive")
    if level <= 0:
        raise ValueError("level must be positive")
    if not math.isfinite(t0):
        raise ValueError("t0 must be finite")
    _, radius_exponent, weber_exponent = cutoff_exponents(alpha)
    data = read_metrics(path)
    tau = data["time"] - t0
    selected = (tau >= tau_window[0]) & (tau <= tau_window[1])
    if np.count_nonzero(selected) < 3:
        raise ValueError(f"Fewer than three rows in tau window for {label}")
    resolved = (
        (data["curvature_radius_cells"] >= min_cells)
        & (data["tip_cell_offset_cells"] <= max_tip_offset_cells)
    )
    admissible = selected & resolved
    raw_indices = np.flatnonzero(selected)
    raw_minimum = raw_indices[np.argmin(data["curvature_radius"][selected])]
    measurement = raw_minimum

    z_speed = np.gradient(data["z_tip"], data["time"], edge_order=1)
    crosscheck = np.abs(z_speed - data["u_z_tip"]) / np.maximum(
        np.abs(data["u_z_tip"]), 1.0e-30
    )
    selected_positions = np.flatnonzero(selected)
    minimum_position = int(np.flatnonzero(selected_positions == measurement)[0])
    boundary_minimum = minimum_position in {0, len(selected_positions) - 1}
    radius = float(data["curvature_radius"][measurement])
    summary: dict[str, object] = {
        "label": label,
        "oh": oh,
        "level": level,
        "t0": t0,
        "alpha": alpha,
        "csv": str(path),
        "tau_window": list(tau_window),
        "row_count": int(np.count_nonzero(selected)),
        "resolved_row_count": int(np.count_nonzero(admissible)),
        "resolved_fraction": float(np.count_nonzero(admissible) / np.count_nonzero(selected)),
        "min_cells": min_cells,
        "max_tip_offset_cells": max_tip_offset_cells,
        "minimum_is_resolved": bool(resolved[raw_minimum]),
        "minimum_is_window_boundary": boundary_minimum,
        "measurement_time": float(data["time"][measurement]),
        "measurement_tau": float(tau[measurement]),
        "curvature_radius": radius,
        "curvature_radius_cells": float(data["curvature_radius_cells"][measurement]),
        "normalised_curvature_radius": radius / oh**2,
        "we_curvature_uz": float(data["we_curvature_uz"][measurement]),
        "we_curvature_speed": float(data["we_curvature_speed"][measurement]),
        "median_tip_kinematic_relative_mismatch": float(np.median(crosscheck[selected])),
        "radius_exponent": radius_exponent,
        "weber_exponent": weber_exponent,
        "radius_theory_factor": oh**radius_exponent,
        "weber_theory_factor": oh**weber_exponent,
    }
    data.update(
        {
            "tau": tau,
            "selected": selected,
            "resolved": resolved,
            "z_speed": z_speed,
            "kinematic_relative_mismatch": crosscheck,
        }
    )
    return summary, data


def fitted_prefactor(values: np.ndarray, factors: np.ndarray) -> tuple[float, float]:
    """Fit one positive prefactor in log space and return RMS log residual."""
    if values.size < 2 or np.any(values <= 0.0) or np.any(factors <= 0.0):
        return math.nan, math.nan
    log_prefactor = float(np.mean(np.log(values) - np.log(factors)))
    residual = np.log(values) - (log_prefactor + np.log(factors))
    return math.exp(log_prefactor), float(np.sqrt(np.mean(residual**2)))


def model_comparison(summaries: Sequence[dict[str, object]]) -> dict[str, object]:
    """Compare conical and constant-in-Oh shapes using resolved minima only."""
    valid = [
        row
        for row in summaries
        if bool(row["minimum_is_resolved"])
        and not bool(row["minimum_is_window_boundary"])
    ]
    result: dict[str, object] = {"resolved_series_count": len(valid)}
    for observable, factor_name in (
        ("normalised_curvature_radius", "radius_theory_factor"),
        ("we_curvature_uz", "weber_theory_factor"),
    ):
        values = np.asarray([float(row[observable]) for row in valid])
        conical = np.asarray([float(row[factor_name]) for row in valid])
        constant = np.ones_like(conical)
        cpref, crms = fitted_prefactor(values, conical)
        ipref, irms = fitted_prefactor(values, constant)
        result[observable] = {
            "conical_prefactor": cpref,
            "conical_log_rms": crms,
            "constant_prefactor": ipref,
            "constant_log_rms": irms,
        }
    return result


def atomic_json(path: Path, value: object) -> None:
    """Atomically write one JSON document."""
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


def configure_matplotlib(use_tex: bool) -> None:
    """Apply compact APS double-column typography."""
    import matplotlib

    matplotlib.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Computer Modern Roman", "DejaVu Serif"],
            "mathtext.fontset": "cm",
            "text.usetex": use_tex,
            "text.latex.preamble": r"\usepackage{amsmath}",
            "axes.linewidth": 1.0,
            "axes.labelsize": 10,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def format_log_decade(value: float, _position: float | None = None) -> str:
    """Render a logarithmic decade as an unambiguous ``1e+/-N`` label."""
    if not math.isfinite(value) or value <= 0.0:
        return ""
    exponent = int(round(math.log10(value)))
    if not math.isclose(value, 10.0**exponent, rel_tol=1.0e-8):
        return ""
    return f"1e{exponent:+d}"


def format_plain_tick(value: float, _position: float | None = None) -> str:
    """Render compact decimal labels without an implicit exponent offset."""
    if not math.isfinite(value):
        return ""
    return f"{value:.3g}"


def plot_report(
    output_stem: Path,
    summaries: Sequence[dict[str, object]],
    arrays: Sequence[dict[str, np.ndarray]],
    comparison: dict[str, object],
    use_tex: bool,
) -> None:
    """Render the time and cross-Oh resolution diagnostic."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

    configure_matplotlib(use_tex)
    fig, axes = plt.subplots(2, 2, figsize=(7.2, 5.7))
    colours = plt.get_cmap("viridis")(np.linspace(0.12, 0.88, len(summaries)))
    for colour, summary, data in zip(colours, summaries, arrays, strict=True):
        chosen = data["selected"]
        resolved = chosen & data["resolved"]
        unresolved = chosen & ~data["resolved"]
        label = str(summary["label"])
        axes[0, 0].plot(data["tau"][resolved], data["curvature_radius"][resolved] / float(summary["oh"])**2,
                        color=colour, lw=1.2, label=label)
        axes[0, 0].plot(data["tau"][unresolved], data["curvature_radius"][unresolved] / float(summary["oh"])**2,
                        "x", color=colour, ms=3.5, mew=0.8)
        axes[0, 1].plot(data["tau"][resolved], data["we_curvature_uz"][resolved],
                        color=colour, lw=1.2, label=label)
        axes[0, 1].plot(data["tau"][unresolved], data["we_curvature_uz"][unresolved],
                        "x", color=colour, ms=3.5, mew=0.8)

    valid = sorted(
        (
            row
            for row in summaries
            if bool(row["minimum_is_resolved"])
            and not bool(row["minimum_is_window_boundary"])
        ),
        key=lambda row: float(row["oh"]),
    )
    if len(valid) >= 2:
        oh = np.asarray([float(row["oh"]) for row in valid])
        radius = np.asarray([float(row["normalised_curvature_radius"]) for row in valid])
        weber = np.asarray([float(row["we_curvature_uz"]) for row in valid])
        axes[1, 0].plot(oh, radius, "o", color="#0072B2", ms=5, label="resolved DNS")
        axes[1, 1].plot(oh, weber, "o", color="#D55E00", ms=5, label="resolved DNS")
        rcomp = comparison["normalised_curvature_radius"]
        wcomp = comparison["we_curvature_uz"]
        rfactor = np.asarray([float(row["radius_theory_factor"]) for row in valid])
        wfactor = np.asarray([float(row["weber_theory_factor"]) for row in valid])
        axes[1, 0].plot(oh, float(rcomp["conical_prefactor"])*rfactor, "-", color="black",
                        lw=1.2, label="conical shape")
        axes[1, 0].axhline(float(rcomp["constant_prefactor"]), color="0.45", ls="--",
                           lw=1.0, label="inertio-capillary")
        axes[1, 1].plot(oh, float(wcomp["conical_prefactor"])*wfactor, "-", color="black",
                        lw=1.2, label="conical shape")
        axes[1, 1].axhline(float(wcomp["constant_prefactor"]), color="0.45", ls="--",
                           lw=1.0, label="inertio-capillary")
        axes[1, 0].set(xlabel=r"$Oh$", ylabel=r"$R_{\kappa,\min}/\ell_\mu$")
        axes[1, 1].set(xlabel=r"$Oh$", ylabel=r"$We_{\kappa,\min}$")
        time_axes = (axes[0, 0], axes[0, 1])
        log_y_axes = tuple(axes.flat)
    else:
        velocity_handles: list[Line2D] = []
        for colour, summary, data in zip(colours, summaries, arrays, strict=True):
            chosen = data["selected"]
            resolved = chosen & data["resolved"]
            unresolved = chosen & ~data["resolved"]
            label = str(summary["label"])
            axes[1, 0].plot(data["tau"][resolved], data["curvature_radius_cells"][resolved],
                            color=colour, lw=1.2, label=label)
            axes[1, 0].plot(data["tau"][unresolved], data["curvature_radius_cells"][unresolved],
                            "x", color=colour, ms=3.5, mew=0.8)
            kinematic = chosen & (data["kinematic_relative_mismatch"] <= 0.02)
            resolved_velocity = kinematic & data["resolved"]
            unresolved_velocity = kinematic & ~data["resolved"]
            axes[1, 1].plot(
                data["tau"][resolved_velocity],
                np.abs(data["u_z_tip"][resolved_velocity]),
                color=colour,
                lw=1.2,
            )
            axes[1, 1].plot(
                data["tau"][unresolved_velocity],
                np.abs(data["u_z_tip"][unresolved_velocity]),
                "x",
                color=colour,
                ms=3.5,
                mew=0.8,
            )
            axes[1, 1].plot(
                data["tau"][resolved_velocity],
                np.abs(data["z_speed"][resolved_velocity]),
                "--",
                color=colour,
                alpha=0.65,
                lw=1.0,
            )
            velocity_handles.append(Line2D([], [], color=colour, lw=1.2, label=label))
        axes[1, 0].axhline(float(summaries[0]["min_cells"]), color="0.45", ls="--",
                           lw=1.0, label="analysis gate")
        axes[1, 0].axhline(1.0, color="0.65", ls=":", lw=1.0,
                           label=r"empirical $\Delta$ floor")
        velocity_handles.extend(
            (
                Line2D([], [], color="0.35", ls="--", lw=1.0,
                       label=r"$|\mathrm{d}z_{\mathrm{tip}}/\mathrm{d}t|$"),
                Line2D([], [], color="0.35", marker="x", ls="", ms=4,
                       label="outside analysis gate"),
            )
        )
        axes[1, 1].legend(handles=velocity_handles, frameon=False, fontsize=7.5,
                          loc="upper right")
        axes[1, 0].set(ylabel=r"$R_\kappa/\Delta_{\mathrm{tip}}$")
        axes[1, 1].set(ylabel="tip speed")
        time_axes = tuple(axes.flat)
        log_y_axes = (axes[0, 0], axes[1, 0])

    axes[0, 0].set(ylabel=r"$R_\kappa/\ell_\mu$")
    axes[0, 1].set(ylabel=r"$We_\kappa=u_{z,\mathrm{tip}}^2R_\kappa$")
    for axis in time_axes:
        axis.set_xscale("log")
        axis.set_xlabel(r"time after inception, $\tau$")
        axis.xaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0,)))
        axis.xaxis.set_major_formatter(FuncFormatter(format_log_decade))
        axis.xaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10)*0.1))
        axis.xaxis.set_minor_formatter(NullFormatter())

    time_values = [
        data["tau"][data["selected"] & (data["tau"] > 0.0)]
        for data in arrays
    ]
    time_min = min(float(np.min(values)) for values in time_values if values.size)
    time_max = max(float(np.max(values)) for values in time_values if values.size)
    for axis in time_axes:
        axis.set_xlim(time_min, time_max)

    for axis in log_y_axes:
        axis.set_yscale("log")
        axis.yaxis.set_major_locator(LogLocator(base=10.0, subs=(1.0, 2.0, 3.0, 4.0, 6.0)))
        axis.yaxis.set_major_formatter(FuncFormatter(format_plain_tick))
        axis.yaxis.set_minor_formatter(NullFormatter())

    if len(valid) < 2:
        axes[0, 1].set_yscale("linear")
        axes[0, 1].set_ylim(bottom=0.0)
        axes[0, 1].yaxis.set_major_formatter(FuncFormatter(format_plain_tick))
        axes[1, 1].set_yscale("linear")
        axes[1, 1].yaxis.set_major_formatter(FuncFormatter(format_plain_tick))

    for index, axis in enumerate(axes.flat):
        axis.tick_params(which="major", direction="out", labelsize=8.5, width=0.8, length=4)
        axis.tick_params(which="minor", direction="out", width=0.6, length=2)
        for spine in axis.spines.values():
            spine.set_linewidth(1.0)
        axis.text(0.02, 0.96, f"({chr(97 + index)})", transform=axis.transAxes,
                  ha="left", va="top", fontsize=10.5)
        if axis.lines and axis is not axes[1, 1]:
            if axis is axes[1, 0] and len(valid) < 2:
                axis.legend(frameon=False, fontsize=7.5, loc="upper left",
                            bbox_to_anchor=(0.02, 0.88))
            else:
                axis.legend(frameon=False, fontsize=7.5)
    fig.tight_layout()
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    for suffix in (".pdf", ".png"):
        fig.savefig(output_stem.with_suffix(suffix), dpi=300, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--series",
        action="append",
        nargs=6,
        metavar=("LABEL", "OH", "LEVEL", "T0", "ALPHA", "CSV"),
        required=True,
        help="series label, Oh, level, inception time, alpha and extractor CSV",
    )
    parser.add_argument("--tau-window", nargs=2, type=float, default=(0.0, 0.01))
    parser.add_argument(
        "--min-cells", type=float, default=2.0,
        help="analysis-quality threshold for R_kappa=1/abs(kappa) (default: 2)",
    )
    parser.add_argument("--max-tip-offset-cells", type=float, default=1.0)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-stem", type=Path, required=True)
    parser.add_argument("--no-tex", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry point."""
    args = build_parser().parse_args(argv)
    tau_window = tuple(float(value) for value in args.tau_window)
    if not all(math.isfinite(value) for value in tau_window) or tau_window[0] < 0.0:
        raise ValueError("tau-window must be finite and nonnegative")
    if tau_window[0] >= tau_window[1]:
        raise ValueError("tau-window must be increasing")
    if not math.isfinite(args.min_cells) or args.min_cells <= 0.0:
        raise ValueError("min-cells must be finite and positive")
    if not math.isfinite(args.max_tip_offset_cells) or args.max_tip_offset_cells <= 0.0:
        raise ValueError("max-tip-offset-cells must be finite and positive")
    summaries: list[dict[str, object]] = []
    arrays: list[dict[str, np.ndarray]] = []
    for label, oh, level, t0, alpha, csv_path in args.series:
        summary, data = summarise_series(
            label,
            float(oh),
            int(level),
            float(t0),
            float(alpha),
            Path(csv_path).resolve(),
            tau_window,
            args.min_cells,
            args.max_tip_offset_cells,
        )
        summaries.append(summary)
        arrays.append(data)
    comparison = model_comparison(summaries)
    report = {
        "schema_version": 1,
        "claim": "resolution diagnostic; R_kappa-to-R_m identification remains unverified",
        "series": summaries,
        "model_comparison": comparison,
    }
    atomic_json(args.output_json.resolve(), report)
    plot_report(args.output_stem.resolve(), summaries, arrays, comparison, not args.no_tex)
    print(
        f"TIP_SCALING_COMPLETE series={len(summaries)} "
        f"resolved={comparison['resolved_series_count']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
