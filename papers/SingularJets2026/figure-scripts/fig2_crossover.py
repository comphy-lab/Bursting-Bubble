"""Quality cutoff and phenomenological crossover used by Fig. 2(b,c)."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class RefinementCut:
    """Boundary before a certified one-level descending mesh staircase."""

    index: int
    first_drop_time: float
    last_retained_time: float


def first_refinement_loss(data, pre_cap: int, reference_time: float | None):
    """Find the first post-reference staircase that falls below ``pre_cap``.

    The returned index is the first excluded row.  The backward walk accepts
    only exact one-level decreases and retains the last row on the preceding
    high-level plateau.  A multi-level jump is ambiguous and fails closed.
    """
    if reference_time is None:
        raise ValueError("A branch reference time is required")
    level = np.asarray(data["maxlevel"], dtype=float)
    time = np.asarray(data["t"], dtype=float)
    if level.ndim != 1 or level.shape != time.shape or not len(level):
        raise ValueError("Level and time columns must be non-empty aligned vectors")
    if np.any(~np.isfinite(level)) or np.any(~np.isfinite(time)):
        raise ValueError("Level and time columns must be finite")
    if np.any(np.diff(time) <= 0.0):
        raise ValueError("Time rows must be strictly increasing")
    if not isinstance(pre_cap, (int, np.integer)) or pre_cap < 1:
        raise ValueError("The pre-inception ceiling must be a positive integer")
    if not np.isfinite(reference_time):
        raise ValueError("The branch reference time must be finite")
    if np.any(level != np.rint(level)):
        raise ValueError("Requested mesh ceilings must be integer valued")
    level = level.astype(int)

    below = np.flatnonzero((time > reference_time) & (level < pre_cap))
    if not len(below):
        return None

    staircase = int(below[0])
    while (
        staircase > 0
        and time[staircase - 1] > reference_time
        and level[staircase - 1] > level[staircase]
    ):
        decrease = level[staircase - 1] - level[staircase]
        if decrease != 1:
            raise ValueError("Refinement-loss staircase contains a multi-level jump")
        staircase -= 1
    if level[staircase] < pre_cap or time[staircase] <= reference_time:
        raise ValueError("No post-reference high-level row precedes refinement loss")
    first_excluded = staircase + 1
    if first_excluded >= len(time) or level[first_excluded - 1] <= level[first_excluded]:
        raise ValueError("Could not certify the start of the refinement-loss staircase")
    return RefinementCut(
        index=first_excluded,
        first_drop_time=float(time[first_excluded]),
        last_retained_time=float(time[first_excluded - 1]),
    )


def truncate_at_refinement_loss(data, cut: RefinementCut | None):
    """Return a column-aligned view ending before the certified mesh drop."""
    if cut is None:
        return data
    return {name: values[: cut.index] for name, values in data.items()}


def paired_log_samples(series, log_bin_indices, target: int = 36):
    """Sample native Q once and derive Weber from those exact same rows."""
    indices = log_bin_indices(series["r_j"], series["Q_j"], target=target)
    radius = np.asarray(series["r_j"])[indices].copy()
    volume_flux = np.asarray(series["Q_j"])[indices].copy()
    weber = volume_flux**2 / (np.pi**2 * radius**3)
    return {
        "r_j": radius,
        "Q_j": volume_flux,
        "q_j": volume_flux / (np.pi * radius),
        "We_j": weber,
        "source_index": indices.copy(),
    }


def crossover_flux(radius, *, amplitude_cone: float, amplitude_gb: float,
                   cone_slope: float, sharpness: float):
    """Evaluate the stable generalized-mean crossover in volume flux."""
    radius = np.asarray(radius, dtype=float)
    parameters = (amplitude_cone, amplitude_gb, cone_slope, sharpness)
    if np.any(~np.isfinite(radius)) or np.any(radius <= 0.0):
        raise ValueError("Crossover radii must be positive and finite")
    if any(not np.isfinite(value) or value <= 0.0 for value in parameters):
        raise ValueError("Crossover parameters must be positive and finite")
    log_radius = np.log(radius)
    log_cone = np.log(amplitude_cone) + cone_slope * log_radius
    log_gb = np.log(amplitude_gb) + log_radius
    return np.exp(
        -np.logaddexp(-sharpness * log_cone, -sharpness * log_gb) / sharpness
    )


def crossover_log_slope(radius, *, amplitude_cone: float, amplitude_gb: float,
                        cone_slope: float, sharpness: float):
    """Return d(ln Q)/d(ln r) for the analytic crossover."""
    radius = np.asarray(radius, dtype=float)
    parameters = (amplitude_cone, amplitude_gb, cone_slope, sharpness)
    if np.any(~np.isfinite(radius)) or np.any(radius <= 0.0):
        raise ValueError("Crossover radii must be positive and finite")
    if any(not np.isfinite(value) or value <= 0.0 for value in parameters):
        raise ValueError("Crossover parameters must be positive and finite")
    log_radius = np.log(radius)
    log_cone_weight = -sharpness * (
        np.log(amplitude_cone) + cone_slope * log_radius
    )
    log_gb_weight = -sharpness * (np.log(amplitude_gb) + log_radius)
    gb_fraction = np.exp(
        log_gb_weight - np.logaddexp(log_cone_weight, log_gb_weight)
    )
    return cone_slope + (1.0 - cone_slope) * gb_fraction
