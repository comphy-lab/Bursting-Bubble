"""Tests for tip-curvature scaling diagnostics."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np


PAPER_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PAPER_ROOT / "tip-curvature/analyse_tip_scaling.py"
sys.path.insert(0, str(SCRIPT.parent))

from analyse_tip_scaling import cutoff_exponents, summarise_series  # noqa: E402


def write_metrics(path: Path, radius_cells: tuple[float, ...] = (5.0, 4.0, 6.0)) -> None:
    """Write the minimum extractor columns for one synthetic series."""
    fieldnames = [
        "time",
        "z_tip",
        "apex_radius",
        "apex_radius_cells",
        "we_apex_uz",
        "we_apex_speed",
        "u_z_tip",
        "delta_tip",
        "tip_cell_offset_cells",
    ]
    rows = []
    for index, cells in enumerate(radius_cells):
        time = 0.50 + 0.01 * index
        radius = 0.001 * cells
        rows.append(
            {
                "time": time,
                "z_tip": time * 10.0,
                "apex_radius": radius,
                "apex_radius_cells": cells,
                "we_apex_uz": 100.0 * radius,
                "we_apex_speed": 100.0 * radius,
                "u_z_tip": 10.0,
                "delta_tip": 0.001,
                "tip_cell_offset_cells": 0.5,
            }
        )
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def test_cutoff_exponents_pin_0629_and_two_thirds() -> None:
    radius, normalised, weber = cutoff_exponents(0.629)
    ic_radius, ic_normalised, ic_weber = cutoff_exponents(2.0 / 3.0)

    assert math.isclose(radius, 2.4379844961240315)
    assert math.isclose(normalised, 0.43798449612403145)
    assert math.isclose(weber, -0.43798449612403145)
    assert math.isclose(ic_radius, 2.0)
    assert math.isclose(ic_normalised, 0.0, abs_tol=1e-14)
    assert math.isclose(ic_weber, 0.0, abs_tol=1e-14)


def test_summary_selects_resolved_interior_minimum(tmp_path: Path) -> None:
    path = tmp_path / "metrics.csv"
    write_metrics(path)

    summary, arrays = summarise_series(
        "synthetic", 0.03, 15, 0.49, 0.629, path, (0.0, 0.04), 4.0, 1.0
    )

    assert summary["minimum_is_resolved"] is True
    assert summary["minimum_is_window_boundary"] is False
    assert math.isclose(float(summary["measurement_time"]), 0.51)
    assert math.isclose(float(summary["apex_radius"]), 0.004)
    assert np.all(arrays["resolved"])
    assert float(summary["median_tip_kinematic_relative_mismatch"]) < 1e-12


def test_summary_flags_grid_limited_series(tmp_path: Path) -> None:
    path = tmp_path / "unresolved.csv"
    write_metrics(path, (1.0, 1.5, 2.0))

    summary, _ = summarise_series(
        "unresolved", 0.03, 16, 0.49, 0.629, path, (0.0, 0.04), 4.0, 1.0
    )

    assert summary["minimum_is_resolved"] is False
    assert summary["resolved_row_count"] == 0


def test_cli_writes_json_and_two_figure_formats(tmp_path: Path) -> None:
    path = tmp_path / "metrics.csv"
    write_metrics(path)
    output_json = tmp_path / "report.json"
    output_stem = tmp_path / "tip-scaling"

    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--series",
            "synthetic",
            "0.03",
            "15",
            "0.49",
            "0.629",
            str(path),
            "--tau-window",
            "0",
            "0.04",
            "--output-json",
            str(output_json),
            "--output-stem",
            str(output_stem),
            "--no-tex",
        ],
        capture_output=True,
        text=True,
        timeout=60,
    )

    assert result.returncode == 0, result.stderr
    assert json.loads(output_json.read_text())["series"][0]["minimum_is_resolved"]
    assert output_stem.with_suffix(".pdf").stat().st_size > 0
    assert output_stem.with_suffix(".png").stat().st_size > 0
