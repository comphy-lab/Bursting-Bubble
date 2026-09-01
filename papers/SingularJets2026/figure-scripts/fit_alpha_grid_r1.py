#!/usr/bin/env python3
"""Fit an R1 flux exponent from explicitly selected resolution series.

Only ``Q_j = 2 pi q_jet`` is fitted. The q and We slopes are algebraic
consequences. ``alpha_flux`` is not the beta-derived geometry alpha and must
never overwrite or be conflated with it.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import tempfile
from dataclasses import dataclass
from io import StringIO
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

COLUMNS = ("i", "dt", "t", "ke", "maxlevel", "r_b", "z_b", "r_base",
           "z_base", "q_jet", "q_l")
DEFAULT_ALPHA = 0.629
INERTIO_CAPILLARY_ALPHA = 2.0 / 3.0


@dataclass(frozen=True)
class LogSeries:
    path: Path
    level: int
    values: Mapping[str, np.ndarray]


@dataclass(frozen=True)
class RunSeries:
    """One named run in physical time order."""
    label: str
    level: int
    path: Path
    r_j: np.ndarray
    Q_j: np.ndarray


@dataclass(frozen=True)
class BinnedPoint:
    label: str
    level: int
    bin_index: int
    r_j: float
    Q_j: float
    raw_count: int


@dataclass(frozen=True)
class Regression:
    slope: float
    intercepts: Mapping[str, float]
    rms: float
    point_count: int


def read_log(path: Path | str, *, time_digits: int = 8) -> LogSeries:
    """Parse/deduplicate an 11-column log; the last restarted row wins."""
    path = Path(path)
    rows: dict[float, list[float]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            parts = line.split()
            if len(parts) != len(COLUMNS):
                continue
            try:
                row = [float(value) for value in parts]
            except ValueError:
                continue
            if np.all(np.isfinite(row)):
                rows[round(row[2], time_digits)] = row
    if not rows:
        raise ValueError(f"No finite 11-column rows found in {path}")
    array = np.asarray([rows[t] for t in sorted(rows)], float)
    match = re.search(r"(?:^|[_-])L(\d+)(?:[_-]|$)", path.stem, re.I)
    level = int(match.group(1)) if match else int(np.max(array[:, 4]))
    return LogSeries(path, level,
                     {name: array[:, index] for index, name in enumerate(COLUMNS)})


def reconnection_time(log: LogSeries, *, pin_radius: float = 0.005,
                      search_after: float = 0.40) -> float | None:
    data = log.values
    mask = (data["t"] > search_after) & (data["r_base"] < pin_radius)
    return float(np.max(data["t"][mask])) if np.any(mask) else None


def post_inception_points(log: LogSeries, *, pin_radius: float = 0.005,
                          search_after: float = 0.40) -> tuple[np.ndarray, np.ndarray]:
    """Derive r_j and Q_j while retaining time order for block inference."""
    data = log.values
    r_j = np.asarray(data["r_base"], float)
    Q_j = 2.0 * np.pi * np.asarray(data["q_jet"], float)
    mask = np.isfinite(r_j) & np.isfinite(Q_j) & (r_j > 0) & (Q_j > 0)
    inception = reconnection_time(log, pin_radius=pin_radius, search_after=search_after)
    if inception is not None:
        mask &= data["t"] > inception
    return r_j[mask], Q_j[mask]


def load_selected_series(specs: Sequence[Sequence[str]], *, pin_radius: float,
                         search_after: float) -> list[RunSeries]:
    """Load explicit LABEL LEVEL LOG triples without pooling cases by level."""
    labels: set[str] = set()
    paths: set[Path] = set()
    result = []
    for label, level_text, path_text in specs:
        if not label or label in labels:
            raise ValueError(f"Series labels must be non-empty and unique: {label!r}")
        path = Path(path_text).resolve()
        if path in paths:
            raise ValueError(f"A log may be selected only once: {path}")
        try:
            level = int(level_text)
        except ValueError as error:
            raise ValueError(f"Invalid level for {label}: {level_text}") from error
        log = read_log(path)
        if level < 1 or log.level != level:
            raise ValueError(f"Explicit L{level} disagrees with inferred L{log.level}: {path}")
        r_j, Q_j = post_inception_points(log, pin_radius=pin_radius,
                                         search_after=search_after)
        if not r_j.size:
            raise ValueError(f"No positive post-inception points in {path}")
        result.append(RunSeries(label, level, path, r_j, Q_j))
        labels.add(label)
        paths.add(path)
    if not result:
        raise ValueError("At least one explicit --series is required")
    return result


def common_radius_support(series: Sequence[RunSeries]) -> tuple[float, float]:
    if not series:
        raise ValueError("No series supplied")
    lower = max(float(run.r_j.min()) for run in series)
    upper = min(float(run.r_j.max()) for run in series)
    if not lower < upper:
        raise ValueError("Selected series have no common radius support")
    return lower, upper


def log_radius_bins(r_j: np.ndarray, Q_j: np.ndarray, *, lower: float,
                    upper: float, bins: int, label: str,
                    level: int) -> list[BinnedPoint]:
    """Use deterministic geometric means; both window boundaries are included."""
    if not 0 < lower < upper or bins < 1:
        raise ValueError("Require 0 < lower < upper and positive bins")
    r_j, Q_j = np.asarray(r_j, float), np.asarray(Q_j, float)
    valid = (np.isfinite(r_j) & np.isfinite(Q_j) & (r_j >= lower) &
             (r_j <= upper) & (r_j > 0) & (Q_j > 0))
    r_j, Q_j = r_j[valid], Q_j[valid]
    if not r_j.size:
        return []
    edges = np.geomspace(lower, upper, bins + 1)
    ids = np.clip(np.searchsorted(edges, r_j, side="right") - 1, 0, bins - 1)
    output = []
    for index in range(bins):
        chosen = ids == index
        if np.any(chosen):
            output.append(BinnedPoint(
                label, level, index,
                float(np.exp(np.mean(np.log(r_j[chosen])))),
                float(np.exp(np.mean(np.log(Q_j[chosen])))),
                int(np.count_nonzero(chosen))))
    return output


def bin_runs(series: Sequence[RunSeries], *, lower: float, upper: float,
             bins: int, min_occupied_bins: int) -> dict[str, list[BinnedPoint]]:
    """Bin each run independently, including distinct cases at the same level."""
    result = {run.label: log_radius_bins(
        run.r_j, run.Q_j, lower=lower, upper=upper, bins=bins,
        label=run.label, level=run.level) for run in series}
    deficient = {label: len(points) for label, points in result.items()
                 if len(points) < min_occupied_bins}
    if deficient:
        detail = ", ".join(f"{label}: {count}" for label, count in deficient.items())
        raise ValueError(f"Insufficient occupied bins (minimum {min_occupied_bins}; {detail})")
    return result


def _xy(points: Sequence[BinnedPoint]) -> tuple[np.ndarray, np.ndarray]:
    return np.log([point.r_j for point in points]), np.log([point.Q_j for point in points])


def ols(points: Sequence[BinnedPoint], *, fixed_slope: float | None = None) -> Regression:
    if len(points) < (2 if fixed_slope is None else 1):
        raise ValueError("Too few binned points")
    x, y = _xy(points)
    if fixed_slope is None:
        denominator = float(np.sum((x - x.mean()) ** 2))
        if denominator <= 0:
            raise ValueError("Radii have zero variance")
        slope = float(np.sum((x-x.mean())*(y-y.mean())) / denominator)
    else:
        slope = float(fixed_slope)
    intercept = float(np.mean(y-slope*x))
    residual = y-intercept-slope*x
    return Regression(slope, {points[0].label: intercept},
                      float(np.sqrt(np.mean(residual**2))), len(points))


def pooled_run_intercepts(points_by_run: Mapping[str, Sequence[BinnedPoint]], *,
                          fixed_slope: float | None = None,
                          aggregation: str = "point-weighted") -> Regression:
    """Fit a common slope with a distinct intercept for every explicit run.

    ``point-weighted`` is the historical pooled least-squares estimator.  In
    ``equal-run`` mode each run contributes its mean squared residual, so a
    finely sampled trajectory cannot dominate a coarser trajectory merely by
    occupying more radius bins.
    """
    arrays = {label: _xy(points) for label, points in points_by_run.items()}
    if not arrays:
        raise ValueError("No binned points")
    if aggregation not in {"point-weighted", "equal-run"}:
        raise ValueError(f"Unknown aggregation: {aggregation}")
    weights = {label: (1.0 if aggregation == "point-weighted" else 1.0/len(x))
               for label, (x, _) in arrays.items()}
    if fixed_slope is None:
        numerator = sum(weights[label]*float(np.sum((x-x.mean())*(y-y.mean())))
                        for label, (x, y) in arrays.items())
        denominator = sum(weights[label]*float(np.sum((x-x.mean())**2))
                          for label, (x, _) in arrays.items())
        if denominator <= 0:
            raise ValueError("Pooled radii have zero within-run variance")
        slope = numerator / denominator
    else:
        slope = float(fixed_slope)
    intercepts = {label: float(np.mean(y-slope*x)) for label, (x, y) in arrays.items()}
    residuals_by_run = {label: y-intercepts[label]-slope*x
                        for label, (x, y) in arrays.items()}
    residuals = np.concatenate(list(residuals_by_run.values()))
    if aggregation == "equal-run":
        rms = float(np.sqrt(np.mean([
            np.mean(residual**2) for residual in residuals_by_run.values()])))
    else:
        rms = float(np.sqrt(np.mean(residuals**2)))
    return Regression(float(slope), intercepts,
                      rms, int(residuals.size))


def alpha_from_flux_slope(slope: float) -> float:
    if not np.isfinite(slope) or slope >= 3:
        raise ValueError(f"Nonphysical flux slope for alpha_flux: {slope}")
    alpha = 1.0 / (3.0-slope)
    if not np.isfinite(alpha) or alpha <= 0:
        raise ValueError(f"Nonpositive alpha_flux from slope: {slope}")
    return alpha


def slope_from_alpha(alpha: float) -> float:
    if not np.isfinite(alpha) or alpha <= 0:
        raise ValueError("alpha must be positive")
    return 3.0 - 1.0/alpha


def derived_slopes(s_Q: float) -> dict[str, float]:
    return {"s_Q": float(s_Q), "s_q": float(s_Q-1), "s_We": float(2*s_Q-3)}


def _block_indices(count: int, width: int, rng: np.random.Generator) -> np.ndarray:
    width = min(width, count)
    result: list[int] = []
    while len(result) < count:
        start = int(rng.integers(0, count))
        result.extend((start+offset) % count for offset in range(width))
    return np.asarray(result[:count])


def temporal_block_inference(
    series: Sequence[RunSeries], *, lower: float, upper: float, bins: int,
    min_occupied_bins: int, replicates: int, block_size: int, seed: int,
    null_alpha: float = INERTIO_CAPILLARY_ALPHA,
    aggregation: str = "point-weighted",
) -> dict[str, object]:
    """Resample raw time blocks, re-bin, and test the alpha_flux=2/3 null."""
    if replicates < 1 or block_size < 1:
        raise ValueError("replicates and block_size must be positive")
    observed = pooled_run_intercepts(bin_runs(
        series, lower=lower, upper=upper, bins=bins,
        min_occupied_bins=min_occupied_bins), aggregation=aggregation)
    alpha_from_flux_slope(observed.slope)
    null_slope = slope_from_alpha(null_alpha)
    pair_rng, null_rng = np.random.default_rng(seed), np.random.default_rng(seed+1)
    pair_slopes: list[float] = []
    null_slopes: list[float] = []
    attempts, max_attempts = 0, max(100, 25*replicates)
    while (len(pair_slopes) < replicates or len(null_slopes) < replicates) and attempts < max_attempts:
        attempts += 1
        if len(pair_slopes) < replicates:
            sampled = []
            for run in series:
                idx = _block_indices(run.r_j.size, block_size, pair_rng)
                sampled.append(RunSeries(run.label, run.level, run.path,
                                         run.r_j[idx], run.Q_j[idx]))
            try:
                slope = pooled_run_intercepts(bin_runs(
                    sampled, lower=lower, upper=upper, bins=bins,
                    min_occupied_bins=min_occupied_bins),
                    aggregation=aggregation).slope
                alpha_from_flux_slope(slope)
            except ValueError:
                pass
            else:
                pair_slopes.append(slope)
        if len(null_slopes) < replicates:
            synthetic = []
            for run in series:
                eligible = (run.r_j >= lower) & (run.r_j <= upper)
                x, y = np.log(run.r_j[eligible]), np.log(run.Q_j[eligible])
                if x.size < 2:
                    raise ValueError(f"Too few raw points in {run.label}")
                intercept = float(np.mean(y-null_slope*x))
                free_intercept = float(np.mean(y-observed.slope*x))
                # Centre the unrestricted-model residuals at the null slope;
                # retaining the null-fit trend would erase power under H1.
                residual = y-free_intercept-observed.slope*x
                residual -= residual.mean()
                idx = _block_indices(x.size, block_size, null_rng)
                synthetic.append(RunSeries(run.label, run.level, run.path,
                    np.exp(x), np.exp(intercept+null_slope*x+residual[idx])))
            try:
                slope = pooled_run_intercepts(bin_runs(
                    synthetic, lower=lower, upper=upper, bins=bins,
                    min_occupied_bins=min_occupied_bins),
                    aggregation=aggregation).slope
                alpha_from_flux_slope(slope)
            except ValueError:
                pass
            else:
                null_slopes.append(slope)
    if len(pair_slopes) < replicates or len(null_slopes) < replicates:
        raise ValueError("Too few valid temporal-block bootstrap replicates")
    pair, null = np.asarray(pair_slopes), np.asarray(null_slopes)
    alphas = np.asarray([alpha_from_flux_slope(value) for value in pair])
    distance = abs(observed.slope-null_slope)
    p_value = (1 + np.count_nonzero(np.abs(null-null_slope) >= distance-1e-14))/(replicates+1)
    return {
        "replicates": replicates, "block_size_raw_observations": block_size,
        "seed": seed, "s_Q_ci95": np.percentile(pair, [2.5, 97.5]).tolist(),
        "alpha_flux_ci95": np.percentile(alphas, [2.5, 97.5]).tolist(),
        "null_test_alpha_flux": null_alpha, "null_test_s_Q": null_slope,
        "null_test_p_two_sided": float(p_value),
        "null_test_method": "null-centred circular temporal-block residual bootstrap",
        "aggregation": aggregation,
    }


def fit_window(series: Sequence[RunSeries], *, lower: float, upper: float,
               bins: int, min_occupied_bins: int,
               aggregation: str = "point-weighted"):
    binned = bin_runs(series, lower=lower, upper=upper, bins=bins,
                      min_occupied_bins=min_occupied_bins)
    pooled = pooled_run_intercepts(binned, aggregation=aggregation)
    alpha = alpha_from_flux_slope(pooled.slope)
    runs = {run.label: run for run in series}
    fixed = {}
    for name, candidate in (("alpha_flux_0p629", DEFAULT_ALPHA),
                            ("alpha_flux_2_over_3", INERTIO_CAPILLARY_ALPHA)):
        result = pooled_run_intercepts(binned, fixed_slope=slope_from_alpha(candidate),
                                       aggregation=aggregation)
        fixed[name] = {"alpha_flux_fixed": candidate, "s_Q": result.slope,
                       "rms_log_Q": result.rms, "point_count": result.point_count}
    summary = {
        "lower": lower, "upper": upper, "aggregation": aggregation,
        "occupied_bins": {label: len(points) for label, points in binned.items()},
        "pooled": {"alpha_flux": alpha, **derived_slopes(pooled.slope),
                   "rms_log_Q": pooled.rms, "point_count": pooled.point_count,
                   "log_intercepts": dict(pooled.intercepts)},
        "per_run": {}, "fixed_alpha_flux_comparisons": fixed,
        "aggregation_comparison": {},
    }
    for mode in ("point-weighted", "equal-run"):
        try:
            result = pooled_run_intercepts(binned, aggregation=mode)
            candidate_alpha = alpha_from_flux_slope(result.slope)
        except ValueError as error:
            # A diagnostic alternative must not erase a valid result from the
            # explicitly selected aggregation mode.
            summary["aggregation_comparison"][mode] = {
                "valid": False, "reason": str(error)}
        else:
            summary["aggregation_comparison"][mode] = {
                "valid": True, "reason": "", "alpha_flux": candidate_alpha,
                **derived_slopes(result.slope), "rms_log_Q": result.rms,
                "point_count": result.point_count,
            }
    for label, points in binned.items():
        result = ols(points)
        summary["per_run"][label] = {
            "level": runs[label].level, "alpha_flux": alpha_from_flux_slope(result.slope),
            **derived_slopes(result.slope), "rms_log_Q": result.rms,
            "point_count": result.point_count, "log_intercept": result.intercepts[label]}
    return summary, binned


def window_sensitivity(series: Sequence[RunSeries], *, lower_grid: Sequence[float],
                       upper_grid: Sequence[float], bins: int,
                       min_occupied_bins: int, replicates: int,
                       block_size: int, seed: int,
                       aggregation: str = "point-weighted"):
    support = common_radius_support(series)
    rows, index = [], 0
    for requested_lower in sorted(set(lower_grid)):
        for requested_upper in sorted(set(upper_grid)):
            lower, upper = max(requested_lower, support[0]), min(requested_upper, support[1])
            row: dict[str, object] = {
                "requested_lower": requested_lower, "requested_upper": requested_upper,
                "effective_lower": lower, "effective_upper": upper,
                "aggregation": aggregation}
            if not lower < upper:
                row.update(valid=False, reason="empty common-window intersection")
            else:
                try:
                    summary, _ = fit_window(series, lower=lower, upper=upper, bins=bins,
                                            min_occupied_bins=min_occupied_bins,
                                            aggregation=aggregation)
                    inference = temporal_block_inference(
                        series, lower=lower, upper=upper, bins=bins,
                        min_occupied_bins=min_occupied_bins, replicates=replicates,
                        block_size=block_size, seed=seed+1009*index,
                        aggregation=aggregation)
                except ValueError as error:
                    row.update(valid=False, reason=str(error))
                else:
                    pooled, fixed = summary["pooled"], summary["fixed_alpha_flux_comparisons"]
                    row.update(valid=True, reason="", occupied_bins=summary["occupied_bins"],
                        N_binned=pooled["point_count"], alpha_flux=pooled["alpha_flux"],
                        alpha_flux_ci95=inference["alpha_flux_ci95"],
                        p_alpha_flux_2_over_3=inference["null_test_p_two_sided"],
                        s_Q=pooled["s_Q"], s_q=pooled["s_q"], s_We=pooled["s_We"],
                        rms_free=pooled["rms_log_Q"],
                        rms_alpha_flux_0p629=fixed["alpha_flux_0p629"]["rms_log_Q"],
                        rms_alpha_flux_2_over_3=fixed["alpha_flux_2_over_3"]["rms_log_Q"])
            rows.append(row)
            index += 1
    return rows, support


def _atomic_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
        handle.write(text)
    os.replace(temporary, path)


def write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    columns = ("requested_lower", "requested_upper", "effective_lower", "effective_upper",
        "aggregation",
        "valid", "reason", "occupied_bins", "N_binned", "alpha_flux",
        "alpha_flux_ci95", "p_alpha_flux_2_over_3", "s_Q", "s_q", "s_We",
        "rms_free", "rms_alpha_flux_0p629", "rms_alpha_flux_2_over_3")
    buffer = StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=columns, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        flat = dict(row)
        for key in ("occupied_bins", "alpha_flux_ci95"):
            if key in flat:
                flat[key] = json.dumps(flat[key], sort_keys=True, separators=(",", ":"))
        writer.writerow(flat)
    _atomic_text(path, buffer.getvalue())


def make_diagnostic_pdf(path: Path, binned, fit, sensitivity, *, use_tex: bool) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, MaxNLocator
    matplotlib.rcParams.update({"font.family": "serif",
        "font.serif": ["Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm", "text.usetex": use_tex,
        "text.latex.preamble": r"\usepackage{amsmath}", "axes.linewidth": 1.2,
        "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, (ax, sensitivity_ax) = plt.subplots(1, 2, figsize=(7.1, 3.15))
    colours = plt.get_cmap("viridis")(np.linspace(.15, .85, len(binned)))
    pooled = fit["pooled"]
    for colour, (label, points) in zip(colours, sorted(binned.items())):
        x, y = np.asarray([p.r_j for p in points]), np.asarray([p.Q_j for p in points])
        ax.loglog(x, y, "o", ms=4, color=colour, label=label)
        line_x = np.geomspace(x.min(), x.max(), 100)
        ax.loglog(line_x, np.exp(pooled["log_intercepts"][label])*line_x**pooled["s_Q"],
                  color=colour, lw=1.5)
    valid = [row for row in sensitivity if row.get("valid")]
    if valid:
        centres = [np.sqrt(row["effective_lower"]*row["effective_upper"]) for row in valid]
        sensitivity_ax.scatter(centres, [row["alpha_flux"] for row in valid], s=28)
    sensitivity_ax.axhline(DEFAULT_ALPHA, color="#444", label=r"$\alpha_{\rm flux}=0.629$")
    sensitivity_ax.axhline(INERTIO_CAPILLARY_ALPHA, color="#777", ls="--",
                           label=r"$\alpha_{\rm flux}=2/3$")
    sensitivity_ax.xaxis.set_major_locator(MaxNLocator(nbins=3))
    sensitivity_ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:.3f}"))
    ax.set(xlabel=r"$r_j$", ylabel=r"$Q_j$")
    sensitivity_ax.set(xlabel="geometric window centre", ylabel=r"$\alpha_{\rm flux}$")
    for axis in (ax, sensitivity_ax):
        axis.tick_params(which="both", direction="out", width=1.1)
        axis.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.stem}.", suffix=".pdf", dir=path.parent)
    os.close(descriptor)
    try:
        fig.savefig(temporary, format="pdf", dpi=300, bbox_inches="tight", pad_inches=.1)
        os.replace(temporary, path)
    finally:
        plt.close(fig)
        Path(temporary).unlink(missing_ok=True)


def analyse(args: argparse.Namespace) -> dict[str, object]:
    series = load_selected_series(args.series, pin_radius=args.pin_radius,
                                  search_after=args.inception_search_after)
    support = common_radius_support(series)
    lower, upper = max(args.window[0], support[0]), min(args.window[1], support[1])
    if not lower < upper:
        raise ValueError("Requested window has no common-run intersection")
    fit, binned = fit_window(series, lower=lower, upper=upper, bins=args.bins,
                             min_occupied_bins=args.min_occupied_bins,
                             aggregation=args.aggregation)
    fit["pooled"]["temporal_block_inference"] = temporal_block_inference(
        series, lower=lower, upper=upper, bins=args.bins,
        min_occupied_bins=args.min_occupied_bins, replicates=args.bootstrap,
        block_size=args.block_size, seed=args.seed,
        aggregation=args.aggregation)
    for offset, run in enumerate(series, start=1):
        fit["per_run"][run.label]["temporal_block_inference"] = (
            temporal_block_inference(
                [run],
                lower=lower,
                upper=upper,
                bins=args.bins,
                min_occupied_bins=args.min_occupied_bins,
                replicates=args.bootstrap,
                block_size=args.block_size,
                seed=args.seed + 100_003 * offset,
                aggregation=args.aggregation,
            )
        )
    sensitivity, _ = window_sensitivity(series,
        lower_grid=args.lower_grid or [args.window[0]],
        upper_grid=args.upper_grid or [args.window[1]], bins=args.bins,
        min_occupied_bins=args.min_occupied_bins, replicates=args.bootstrap,
        block_size=args.block_size, seed=args.seed,
        aggregation=args.aggregation)
    report = {"schema_version": 2, "physics_group": args.physics_group,
        "method": {"independent_observable": "Q_j = 2*pi*q_jet",
            "regression": "OLS(log(Q_j),log(r_j)); run-specific intercepts",
            "derived_relations": {"alpha_flux": "1/(3-s_Q)", "s_q": "s_Q-1",
                                  "s_We": "2*s_Q-3"},
            "case_isolation": "each --series is a separate trajectory",
            "aggregation": args.aggregation,
            "aggregation_note": ("each run has equal objective weight" if
                args.aggregation == "equal-run" else
                "each occupied binned point has equal objective weight"),
            "serial_correlation": "raw temporal blocks are resampled and re-binned",
            "inference_scope": ("confidence intervals are conditional on the selected "
                                "trajectories and do not include between-run or "
                                "resolution heterogeneity"),
            "interpretation": "alpha_flux is not the beta-derived geometry alpha"},
        "inputs": [{"label": run.label, "level": run.level, "path": str(run.path),
                    "post_inception_rows": int(run.r_j.size)} for run in series],
        "common_radius_support": list(support), "requested_window": list(args.window),
        "fit": fit, "binned_points": {
            label: [{"level": p.level, "bin_index": p.bin_index, "r_j": p.r_j,
                     "Q_j": p.Q_j, "raw_count": p.raw_count} for p in points]
            for label, points in binned.items()}, "window_sensitivity": sensitivity}
    _atomic_text(args.output_json, json.dumps(report, indent=2, sort_keys=True)+"\n")
    write_csv(args.output_csv, sensitivity)
    if args.output_pdf:
        no_tex = os.environ.get("SINGULARJETS_NO_TEX", "").lower() in {"1", "true", "yes", "on"}
        make_diagnostic_pdf(args.output_pdf, binned, fit, sensitivity,
                            use_tex=not (args.no_tex or no_tex))
    return report


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--physics-group", required=True)
    parser.add_argument("--series", action="append", nargs=3,
                        metavar=("LABEL", "LEVEL", "LOG"), required=True)
    parser.add_argument("--window", type=float, nargs=2, default=(.005, .023952))
    parser.add_argument("--lower-grid", type=float, nargs="+")
    parser.add_argument("--upper-grid", type=float, nargs="+")
    parser.add_argument("--bins", type=int, default=24)
    parser.add_argument("--min-occupied-bins", type=int, default=5)
    parser.add_argument("--aggregation", choices=("point-weighted", "equal-run"),
                        default="point-weighted")
    parser.add_argument("--bootstrap", type=int, default=2000)
    parser.add_argument("--block-size", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260901)
    parser.add_argument("--pin-radius", type=float, default=.005)
    parser.add_argument("--inception-search-after", type=float, default=.40)
    parser.add_argument("--output-json", type=Path, default=Path("fit-alpha-grid-r1.json"))
    parser.add_argument("--output-csv", type=Path, default=Path("fit-alpha-grid-r1.csv"))
    parser.add_argument("--output-pdf", type=Path)
    parser.add_argument("--no-tex", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.bins < 1 or args.min_occupied_bins < 2:
        parser.error("--bins positive; --min-occupied-bins at least 2")
    if args.bootstrap < 1 or args.block_size < 1:
        parser.error("--bootstrap and --block-size must be positive")
    if not 0 < args.window[0] < args.window[1]:
        parser.error("--window must satisfy 0 < LOWER < UPPER")
    analyse(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
