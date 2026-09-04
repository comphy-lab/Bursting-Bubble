"""Tests for the online tip-metrics acceptance gate."""

from __future__ import annotations

from pathlib import Path
import sys

import numpy as np


PAPER_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = PAPER_ROOT / "tip-curvature"
sys.path.insert(0, str(SCRIPT_DIR))

from validate_online_tip_metrics import read_rows, validate  # noqa: E402


def metric_row(
    index: int,
    radius: float,
    *,
    pinched: int = 0,
    level: int = 18,
    status: int = 3,
) -> str:
    """Return one syntactically valid online-log row."""
    kappa = 2.0 / radius
    values = (
        index,
        1.0e-6,
        0.5 + index * 1.0e-6,
        1,
        pinched,
        1 + pinched,
        status,
        1.0,
        0.0,
        0.9999,
        0.0001,
        kappa if status == 3 else -1000.0,
        10.0 if status == 3 else -1000.0,
        0.0 if status == 3 else -1000.0,
        10.0 if status == 3 else -1000.0,
        0.001 if status == 3 else -1000.0,
        level if status == 3 else -1,
        0.5 if status == 3 else -1000.0,
        0.02,
        -1.0,
        0.01,
        1.0,
    )
    return " ".join(str(value) for value in values)


def write_log(path: Path, radii: np.ndarray, *, final_pinch: bool = True) -> None:
    """Write one synthetic versioned sidecar."""
    lines = ["# tip-metrics-v1", "# segment case=6401"]
    lines.extend(metric_row(index, float(radius)) for index, radius in enumerate(radii))
    if final_pinch:
        lines.append(metric_row(len(radii), float(radii[-1]), pinched=1, status=0))
    path.write_text("\n".join(lines) + "\n")


def test_accepts_resolved_interior_minimum_and_rebound(tmp_path: Path) -> None:
    path = tmp_path / "tip_metrics.log"
    radii = np.r_[np.linspace(0.01, 0.005, 31), np.linspace(0.0052, 0.008, 30)]
    write_log(path, radii)

    report = validate(read_rows(path), 18, 50, 4.0, 1.1, 20)

    assert report["accepted"] is True
    assert report["minimum_interior"] is True
    assert report["minimum_radius_cells"] == 5.0
    assert report["tip_pinched"] is True


def test_rejects_grid_limited_or_incomplete_run(tmp_path: Path) -> None:
    path = tmp_path / "tip_metrics.log"
    radii = np.r_[np.linspace(0.01, 0.002, 31), np.linspace(0.0022, 0.004, 30)]
    write_log(path, radii, final_pinch=False)

    report = validate(read_rows(path), 18, 50, 4.0, 1.1, 20)

    assert report["accepted"] is False
    assert any("minimum spans" in reason for reason in report["reasons"])
    assert any("tip-pinch" in reason for reason in report["reasons"])
