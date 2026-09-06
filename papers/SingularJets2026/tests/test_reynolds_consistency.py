"""Tests for Reynolds/Weber self-consistency diagnostics."""

from __future__ import annotations

import csv
import math
from pathlib import Path
import sys


SCRIPT_DIR = Path(__file__).resolve().parents[1] / "figure-scripts"
sys.path.insert(0, str(SCRIPT_DIR))

from make_reynolds_consistency import cutoff_estimate, similarity_exponents  # noqa: E402


def test_similarity_exponents_for_alpha_0629() -> None:
    re_slope, we_slope = similarity_exponents(0.629)

    assert math.isclose(re_slope, 0.4101748807631161)
    assert math.isclose(we_slope, -0.17965023847376795)


def test_cutoff_estimate_uses_earliest_admissible_time(tmp_path: Path) -> None:
    path = tmp_path / "tip.csv"
    fields = (
        "time", "inverse_mean_curvature", "u_z_tip", "delta_tip",
        "tip_cell_offset_cells",
    )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(
            (
                dict(zip(fields, (0.5000005, 0.002, 8.0, 0.001, 0.5), strict=True)),
                dict(zip(fields, (0.5000020, 0.003, -10.0, 0.001, 0.4), strict=True)),
                dict(zip(fields, (0.5000030, 0.004, 12.0, 0.001, 0.3), strict=True)),
            )
        )

    estimate = cutoff_estimate("L15", 15, 0.5, path, 0.03)

    assert math.isclose(estimate.tau, 2.0e-6, rel_tol=1.0e-9)
    assert math.isclose(estimate.curvature_radius_cells, 3.0)
    assert math.isclose(estimate.u_z_tip, 10.0)
    assert math.isclose(estimate.re_m, 1.0)
    assert math.isclose(estimate.we_m, 0.3)
