"""Static contracts for the logged logarithmic-Oh sweep."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np


PAPER_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = PAPER_ROOT.parents[1]


def test_sweep_has_seven_log_points_and_one_resolution_anchor() -> None:
    with (PAPER_ROOT / "tip-curvature/logspace-sweep.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    assert [int(row["index"]) for row in rows] == list(range(8))
    assert len({row["case_no"] for row in rows}) == 8
    primary = rows[:7]
    oh = np.asarray([float(row["oh"]) for row in primary])
    ratios = oh[1:] / oh[:-1]
    assert np.max(np.abs(ratios / np.mean(ratios) - 1.0)) < 5.0e-5
    assert all(row["pre_level"] == "12" and row["post_level"] == "18" for row in primary)
    assert rows[7]["oh"] == rows[0]["oh"]
    assert rows[7]["pre_level"] == "12" and rows[7]["post_level"] == "19"


def test_online_sidecar_does_not_change_established_log_schema() -> None:
    source = (REPO_ROOT / "simulationCases/burstingBubble-drillResolution.c").read_text()

    assert "i dt t ke maxlevel r_b z_b r_base z_base q_jet q_l\\n" in source
    assert "# tip-metrics-v1\\n" in source
    assert "drillHoldMaxUntilTipPinch" in source
    assert 'fopen("tip_metrics.log", "a+")' in source


def test_fresh_case_runner_refuses_implicit_restart() -> None:
    runner = (REPO_ROOT / "runSnelliusDrillOhSweep.sbatch").read_text()

    assert "refusing an implicit restart" in runner
    assert 'find "${CASE_DIR}" -mindepth 1 -maxdepth 1 -print -quit' in runner
    assert 'mkdir -p "${CASE_DIR}"' not in runner
