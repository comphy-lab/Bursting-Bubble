#!/usr/bin/env python3
"""Quantitatively compare candidate exponents for the Figure 3 collapse.

Each snapshot contributes equally.  Facet segments are chained, uniformly
resampled in arclength, rescaled by ``|t-t0|**alpha``, and compared using a
scale-normalised symmetric Chamfer RMS.  Candidate exponents are evaluated on
identical snapshot sets for every t0-fit-window/offset sensitivity scenario.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import tempfile
from collections import defaultdict
from functools import lru_cache
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
from scipy.spatial import cKDTree


def read_index(path: Path) -> tuple[float | None, dict[float, tuple[float, float]]]:
    """Return the recorded t0 and ``time -> (z_base, r_j)`` mapping."""
    t0, rows = None, {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "t0" and len(fields) >= 2:
                t0 = float(fields[1])
                continue
            if len(fields) < 3:
                continue
            try:
                time, radius, z_base = map(float, fields[:3])
            except ValueError:
                continue
            if np.all(np.isfinite((time, radius, z_base))):
                rows[time] = (z_base, radius)
    if not rows:
        raise ValueError(f"No finite index rows in {path}")
    return t0, rows


def fit_t0(rows: Mapping[float, tuple[float, float]], alpha: float,
           radius_window: tuple[float, float]) -> float:
    """Fit the virtual origin from ``r_j**(1/alpha) = A (t-t0)``."""
    times = np.asarray(sorted(rows))
    radii = np.asarray([rows[time][1] for time in times])
    chosen = ((radii >= radius_window[0]) & (radii <= radius_window[1]) &
              np.isfinite(radii) & (radii > 0))
    if np.count_nonzero(chosen) < 3:
        raise ValueError(f"Fewer than three t0-fit points in {radius_window}")
    slope, intercept = np.polyfit(times[chosen], radii[chosen]**(1.0/alpha), 1)
    if not np.isfinite(slope) or slope <= 0:
        raise ValueError("Nonpositive t0-fit slope")
    return float(-intercept/slope)


def read_segments(path: Path) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    """Read getFacet ``z r`` segment pairs."""
    points = []
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        if len(fields) >= 2:
            try:
                point = (float(fields[0]), float(fields[1]))
            except ValueError:
                continue
            if np.all(np.isfinite(point)):
                points.append(point)
    if len(points) % 2:
        raise ValueError(f"Odd number of facet endpoints in {path}")
    return list(zip(points[::2], points[1::2]))


def chain_segments(segments, *, gap: float = 1.2e-2):
    """Chain unordered facet segments with the Figure 3 spatial-hash rule."""
    used = [False] * len(segments)
    cell = max(gap, 1e-3)
    buckets = defaultdict(list)

    def key(point):
        return tuple(int(math.floor(value/cell)) for value in point)

    for index, (first, second) in enumerate(segments):
        buckets[key(first)].append((index, 0))
        buckets[key(second)].append((index, 1))

    def nearest(point, excluded):
        base = key(point)
        best, best_distance = None, gap
        for dz in (-1, 0, 1):
            for dr in (-1, 0, 1):
                for index, end in buckets.get((base[0]+dz, base[1]+dr), ()):
                    if used[index] or index == excluded:
                        continue
                    candidate = segments[index][end]
                    distance = math.hypot(candidate[0]-point[0], candidate[1]-point[1])
                    if distance < best_distance:
                        best, best_distance = (index, end), distance
        return best

    paths = []
    for start in range(len(segments)):
        if used[start]:
            continue
        used[start] = True
        first, second = segments[start]
        path = [first, second]
        for at_head, current in ((False, second), (True, first)):
            excluded = start
            while (match := nearest(current, excluded)) is not None:
                index, end = match
                used[index] = True
                current = segments[index][1-end]
                path.insert(0, current) if at_head else path.append(current)
                excluded = index
        paths.append(np.asarray(path, float))
    return paths


@lru_cache(maxsize=None)
def chained_file(path: Path) -> tuple[np.ndarray, ...]:
    """Cache topology work; sensitivity scenarios reuse the same raw facets."""
    return tuple(chain_segments(read_segments(path)))


def resample_path(path: np.ndarray, count: int) -> np.ndarray:
    distances = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(path, axis=0), axis=1))]
    keep = np.r_[True, np.diff(distances) > 0]
    path, distances = path[keep], distances[keep]
    if len(path) < 2 or distances[-1] <= 0:
        return np.empty((0, 2))
    targets = np.linspace(0, distances[-1], count)
    return np.column_stack([np.interp(targets, distances, path[:, axis])
                            for axis in range(2)])


def profile(path: Path, *, z_base: float, scale: float, z_window: tuple[float, float],
            samples: int, component_fraction: float = 0.3) -> np.ndarray:
    paths = chained_file(path)
    if not paths:
        raise ValueError(f"No facet paths in {path}")
    largest = max(len(item) for item in paths)
    paths = [item for item in paths if len(item) >= component_fraction*largest]
    scaled_paths = []
    for item in paths:
        scaled = np.column_stack((np.abs(item[:, 1])/scale,
                                  (item[:, 0]-z_base)/scale))
        # The comparison window belongs to the collapsed coordinates.  A raw
        # axial crop would select a different physical shape at every tau.
        scaled = scaled[(scaled[:, 1] >= z_window[0]) &
                        (scaled[:, 1] <= z_window[1])]
        if len(scaled) >= 2:
            scaled_paths.append(scaled)
    if not scaled_paths:
        raise ValueError(f"No facets remain inside z window for {path}")
    if len(scaled_paths) > samples:
        raise ValueError("Fewer samples than retained interface components")

    # Equal component representation to within one point, with a hard total
    # budget.  Truncating a concatenated over-allocation would systematically
    # discard the components appearing last.
    base, remainder = divmod(samples, len(scaled_paths))
    output = []
    for index, scaled in enumerate(scaled_paths):
        sampled = resample_path(scaled, base+(index < remainder))
        if sampled.size:
            output.append(sampled)
    return np.concatenate(output)


def normalised_chamfer(first: np.ndarray, second: np.ndarray) -> float:
    first_to_second = cKDTree(second).query(first, workers=1)[0]
    second_to_first = cKDTree(first).query(second, workers=1)[0]
    directed = np.r_[first_to_second**2, second_to_first**2]
    size = 0.5*(np.mean(np.sum(first**2, axis=1))+
                np.mean(np.sum(second**2, axis=1)))
    if size <= 0:
        raise ValueError("Degenerate zero-size profiles")
    return float(np.sqrt(np.mean(directed)/size))


def collapse_rms(profiles: Sequence[np.ndarray]) -> float:
    if len(profiles) < 2:
        raise ValueError("At least two profiles are required")
    distances = [normalised_chamfer(profiles[i], profiles[j])
                 for i in range(len(profiles)) for j in range(i+1, len(profiles))]
    return float(np.sqrt(np.mean(np.square(distances))))


def analyse(args: argparse.Namespace) -> dict[str, object]:
    facets = args.facets.resolve()
    recorded_t0, post = read_index(facets/"index.txt")
    _, pre = read_index(facets/"index_pre.txt")
    candidates = [(label, float(alpha)) for label, alpha in args.candidate]
    scenarios, scenario_id = [], 0
    for window in args.t0_fit_window:
        fitted = {label: fit_t0(post, alpha, tuple(window))
                  for label, alpha in candidates}
        for offset in args.t0_offset:
            t0s = {label: value+offset for label, value in fitted.items()}
            common = {}
            for phase, rows, prefix, sign in (
                ("pre", pre, "facetpremain", -1),
                ("post", post, "facetmain", +1),
            ):
                times = []
                for time in sorted(rows):
                    path = facets/f"{prefix}_{time:.6f}.txt"
                    if not path.is_file():
                        continue
                    taus = [sign*(time-t0s[label]) for label, _ in candidates]
                    if all(args.tau_window[0] <= tau <= args.tau_window[1]
                           for tau in taus):
                        times.append(time)
                if len(times) < args.min_profiles:
                    raise ValueError(f"Only {len(times)} common {phase} profiles in scenario")
                common[phase] = times

            for label, alpha in candidates:
                row = {"scenario": scenario_id, "candidate": label, "alpha": alpha,
                       "t0_fit_radius_window": list(window), "t0_offset": offset,
                       "t0": t0s[label], "recorded_t0": recorded_t0}
                for phase, rows, prefix, sign in (
                    ("pre", pre, "facetpremain", -1),
                    ("post", post, "facetmain", +1),
                ):
                    profiles = []
                    for time in common[phase]:
                        tau = sign*(time-t0s[label])
                        profiles.append(profile(
                            facets/f"{prefix}_{time:.6f}.txt",
                            z_base=rows[time][0], scale=tau**alpha,
                            z_window=tuple(args.z_window), samples=args.samples))
                    row[f"{phase}_profile_count"] = len(profiles)
                    row[f"{phase}_normalised_chamfer_rms"] = collapse_rms(profiles)
                scenarios.append(row)
            scenario_id += 1

    report = {
        "schema_version": 1,
        "method": {
            "metric": "equal-snapshot pairwise symmetric Chamfer RMS",
            "normalisation": "pairwise pooled RMS distance from shifted origin",
            "t0": "alpha-specific r_j**(1/alpha) linear fit plus explicit offsets",
            "comparison": "identical snapshot set within each sensitivity scenario",
        },
        "facets": str(facets), "candidates": dict(candidates),
        "tau_window": list(args.tau_window), "z_window": list(args.z_window),
        "scenarios": scenarios,
    }
    atomic_text(args.output_json, json.dumps(report, indent=2, sort_keys=True)+"\n")
    if args.output_pdf:
        diagnostic_pdf(args.output_pdf, scenarios, use_tex=not args.no_tex)
    return report


def atomic_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
        handle.write(content)
    os.replace(temporary, path)


def diagnostic_pdf(path: Path, rows: Sequence[Mapping[str, object]], *, use_tex: bool) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams.update({"font.family": "serif",
        "font.serif": ["Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm", "text.usetex": use_tex,
        "text.latex.preamble": r"\usepackage{amsmath}", "axes.linewidth": 1.2,
        "pdf.fonttype": 42, "ps.fonttype": 42})
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.0), sharex=True)
    labels = sorted({str(row["candidate"]) for row in rows})
    colours = plt.get_cmap("viridis")(np.linspace(.15, .85, len(labels)))
    for axis, phase in zip(axes, ("pre", "post")):
        for colour, label in zip(colours, labels):
            selected = [row for row in rows if row["candidate"] == label]
            axis.plot([row["scenario"] for row in selected],
                      [row[f"{phase}_normalised_chamfer_rms"] for row in selected],
                      "o-", color=colour, lw=1.5, ms=4, label=label)
        axis.set(xlabel="sensitivity scenario", ylabel="normalised collapse RMS",
                 title=phase)
        axis.tick_params(which="both", direction="out", width=1.1)
        for spine in axis.spines.values():
            spine.set_linewidth(1.2)
        axis.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, format="pdf", dpi=300, bbox_inches="tight", pad_inches=.1)
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--facets", type=Path, required=True)
    parser.add_argument("--candidate", action="append", nargs=2,
                        metavar=("LABEL", "ALPHA"), default=[])
    parser.add_argument("--t0-fit-window", action="append", nargs=2, type=float,
                        default=[])
    parser.add_argument("--t0-offset", nargs="+", type=float,
                        default=(-5e-5, 0.0, 5e-5))
    parser.add_argument("--tau-window", nargs=2, type=float, default=(1.5e-4, 1e-2))
    parser.add_argument("--z-window", nargs=2, type=float, default=(-.5, 8.0))
    parser.add_argument("--samples", type=int, default=160)
    parser.add_argument("--min-profiles", type=int, default=5)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-pdf", type=Path)
    parser.add_argument("--no-tex", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.candidate:
        args.candidate = [("alpha_beta", "0.629"), ("two_thirds", str(2/3))]
    if not args.t0_fit_window:
        args.t0_fit_window = [(0.01, 0.05), (0.01, 0.06), (0.012, 0.06)]
    if len({label for label, _ in args.candidate}) != len(args.candidate):
        parser.error("candidate labels must be unique")
    if any(float(alpha) <= 0 for _, alpha in args.candidate):
        parser.error("candidate alpha values must be positive")
    if not 0 < args.tau_window[0] < args.tau_window[1]:
        parser.error("--tau-window must satisfy 0 < LOWER < UPPER")
    if not args.z_window[0] < args.z_window[1] or args.samples < 12:
        parser.error("invalid z window or too few samples")
    analyse(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
