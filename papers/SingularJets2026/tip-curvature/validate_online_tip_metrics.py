#!/usr/bin/env python3
r"""
# Acceptance gate for online tip metrics

Validate one completed short-window run before dependent sweep jobs are
released. The gate requires a persistent pinch termination, enough valid
connected pre-pinch rows at the requested tip-cell level, an interior minimum,
a resolved minimum radius, and a sustained rebound after that minimum.

This is a software/data-quality gate. Passing it does not identify the measured
$R_\kappa$ with the theoretical cutoff radius $R_m$.
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path
import tempfile
from typing import Sequence

import numpy as np

from analyse_tip_scaling import ONLINE_COLUMNS


def read_rows(path: Path) -> list[dict[str, float]]:
    """Read strict numeric rows from a `tip-metrics-v1` sidecar."""
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].strip() != "# tip-metrics-v1":
        raise ValueError(f"Missing tip-metrics-v1 schema header: {path}")
    rows: list[dict[str, float]] = []
    for line in lines:
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        if len(fields) != len(ONLINE_COLUMNS):
            raise ValueError(
                f"Expected {len(ONLINE_COLUMNS)} columns in {path}, found {len(fields)}"
            )
        values = np.asarray([float(field) for field in fields])
        if not np.all(np.isfinite(values)):
            raise ValueError(f"Non-finite row in {path}")
        rows.append(dict(zip(ONLINE_COLUMNS, values, strict=True)))
    if not rows:
        raise ValueError(f"No numeric rows in {path}")
    if any(second["time"] <= first["time"] for first, second in zip(rows, rows[1:])):
        raise ValueError(f"Times are not strictly increasing in {path}")
    return rows


def validate(
    rows: Sequence[dict[str, float]],
    target_level: int,
    min_valid_rows: int,
    min_cells: float,
    rebound_factor: float,
    rebound_rows: int,
) -> dict[str, object]:
    """Evaluate the deterministic short-window acceptance predicate."""
    if target_level <= 0 or min_valid_rows < 3 or min_cells <= 0.0:
        raise ValueError("target level, row count and cell threshold must be positive")
    if rebound_factor <= 1.0 or rebound_rows <= 0:
        raise ValueError("rebound factor must exceed one and rebound rows be positive")
    valid = [
        row
        for row in rows
        if row["jet_formed"] == 1
        and row["tip_pinched"] == 0
        and row["tip_status"] == 3
        and int(row["level_tip"]) == target_level
    ]
    reasons: list[str] = []
    if rows[-1]["tip_pinched"] != 1:
        reasons.append("run did not reach the persistent tip-pinch latch")
    if len(valid) < min_valid_rows:
        reasons.append(f"only {len(valid)} valid target-level rows; need {min_valid_rows}")

    minimum_index = -1
    minimum_radius = math.nan
    minimum_cells_value = math.nan
    rebounded = False
    if valid:
        radii = np.asarray([2.0 / abs(row["kappa_mean"]) for row in valid])
        cells = np.asarray(
            [radius / row["delta_tip"] for radius, row in zip(radii, valid, strict=True)]
        )
        minimum_index = int(np.argmin(radii))
        minimum_radius = float(radii[minimum_index])
        minimum_cells_value = float(cells[minimum_index])
        if minimum_index in {0, len(valid) - 1}:
            reasons.append("minimum lies on the valid time-window boundary")
        if minimum_cells_value < min_cells:
            reasons.append(
                f"minimum spans {minimum_cells_value:.3f} cells; need {min_cells:.3f}"
            )
        tail = radii[minimum_index + 1 :]
        if tail.size >= rebound_rows:
            rebounded = bool(np.median(tail[-rebound_rows:]) >= rebound_factor * minimum_radius)
        if not rebounded:
            reasons.append("no sustained post-minimum curvature-radius rebound")

    return {
        "accepted": not reasons,
        "reasons": reasons,
        "row_count": len(rows),
        "valid_target_level_rows": len(valid),
        "target_level": target_level,
        "minimum_index": minimum_index,
        "minimum_radius": minimum_radius,
        "minimum_radius_cells": minimum_cells_value,
        "minimum_interior": minimum_index > 0 and minimum_index < len(valid) - 1,
        "rebounded": rebounded,
        "tip_pinched": bool(rows[-1]["tip_pinched"] == 1),
    }


def atomic_json(path: Path, value: object) -> None:
    """Atomically write one validation receipt."""
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


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--target-level", type=int, required=True)
    parser.add_argument("--min-valid-rows", type=int, default=50)
    parser.add_argument("--min-cells", type=float, default=4.0)
    parser.add_argument("--rebound-factor", type=float, default=1.1)
    parser.add_argument("--rebound-rows", type=int, default=20)
    parser.add_argument("--output-json", type=Path, required=True)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry point."""
    args = build_parser().parse_args(argv)
    report = validate(
        read_rows(args.log.resolve()),
        args.target_level,
        args.min_valid_rows,
        args.min_cells,
        args.rebound_factor,
        args.rebound_rows,
    )
    atomic_json(args.output_json.resolve(), report)
    print(
        f"TIP_METRICS_GATE accepted={int(bool(report['accepted']))} "
        f"valid={report['valid_target_level_rows']} "
        f"minimum_cells={report['minimum_radius_cells']}"
    )
    return 0 if report["accepted"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
