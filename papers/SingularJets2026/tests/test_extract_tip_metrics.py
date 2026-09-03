"""Regression tests for the parallel tip-metrics extractor."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import subprocess
import sys


PAPER_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PAPER_ROOT / "tip-curvature/extract_tip_metrics.py"
HELPER_SOURCE = PAPER_ROOT.parents[1] / "postProcess/getTipMetrics.c"


def make_case(tmp_path: Path) -> Path:
    """Create a synthetic dump directory with path metacharacters."""
    case = tmp_path / "case with spaces;literal-$value"
    intermediate = case / "intermediate"
    intermediate.mkdir(parents=True)
    for time in (0.50, 0.51, 0.52, 0.53):
        (intermediate / f"snapshot-{time:.2f}").write_text("synthetic\n")
    return case


def make_helper(tmp_path: Path, fail_time: str | None = None) -> Path:
    """Create a deterministic stand-in for the compiled Basilisk helper."""
    helper = tmp_path / "getTipMetrics fake"
    helper.write_text(
        "#!/usr/bin/env python3\n"
        "import pathlib, sys\n"
        "time = float(pathlib.Path(sys.argv[1]).name.rsplit('-', 1)[1])\n"
        f"fail_time = {fail_time!r}\n"
        "if fail_time is not None and fail_time in pathlib.Path(sys.argv[1]).name:\n"
        "    raise SystemExit(7)\n"
        "kappa = -100.0 - 10.0*(time - 0.5)\n"
        "print('TIP_METRICS', time, 1.0 + time, 0.0, 1.0 + time, 0.0005, "
        "kappa, 10.0, 0.0, 10.0, 0.001, 14, 0.5, 1, 10, 2, 1, 0, file=sys.stderr)\n"
    )
    helper.chmod(0o755)
    return helper


def run_extractor(
    case: Path, helper: Path, output: Path, cpus: int, *extra: str
) -> subprocess.CompletedProcess[str]:
    """Run the extractor with fixed scientific provenance metadata."""
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "192"
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--case",
            str(case),
            "--out",
            str(output),
            "--gettip",
            str(helper),
            "--cpus",
            str(cpus),
            "--oh",
            "0.03",
            "--bond",
            "0",
            "--pre-level",
            "14",
            "--post-level",
            "15",
            "--t0",
            "0.49",
            "--t0-protocol",
            "synthetic bracket",
            "--lineage",
            "synthetic-case",
            "--skip-video",
            *extra,
        ],
        capture_output=True,
        text=True,
        env=environment,
        timeout=60,
    )


def files(path: Path) -> dict[str, bytes]:
    """Return installed output bytes keyed by relative path."""
    return {
        item.relative_to(path).as_posix(): item.read_bytes()
        for item in sorted(path.rglob("*"))
        if item.is_file()
    }


def test_serial_and_parallel_outputs_are_byte_identical(tmp_path: Path) -> None:
    case = make_case(tmp_path)
    helper = make_helper(tmp_path)
    serial = tmp_path / "serial"
    parallel = tmp_path / "parallel"

    serial_result = run_extractor(case, helper, serial, 1)
    parallel_result = run_extractor(case, helper, parallel, 4)

    assert serial_result.returncode == 0, serial_result.stderr
    assert parallel_result.returncode == 0, parallel_result.stderr
    assert files(serial) == files(parallel)
    with (serial / "tip_metrics.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 4
    assert float(rows[0]["apex_radius"]) == 0.02
    assert float(rows[0]["apex_radius_cells"]) == 20.0
    assert float(rows[0]["we_apex_uz"]) == 2.0
    assert float(rows[0]["tip_cell_offset_cells"]) == 0.5


def test_worker_failure_publishes_no_partial_result(tmp_path: Path) -> None:
    case = make_case(tmp_path)
    helper = make_helper(tmp_path, "0.52")
    output = tmp_path / "failed"

    result = run_extractor(case, helper, output, 4)

    assert result.returncode != 0
    assert files(output) == {}


def test_manifest_records_scientific_provenance(tmp_path: Path) -> None:
    case = make_case(tmp_path)
    helper = make_helper(tmp_path)
    output = tmp_path / "output"

    result = run_extractor(case, helper, output, 1, "--bridged")

    assert result.returncode == 0, result.stderr
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["row_count"] == 4
    assert manifest["provenance"] == {
        "bond": 0.0,
        "bridged": True,
        "lineage": "synthetic-case",
        "oh": 0.03,
        "post_level": 15,
        "pre_level": 14,
        "t0": 0.49,
        "t0_protocol": "synthetic bracket",
    }
    assert all({"name", "size", "mtime_ns"} <= set(row) for row in manifest["snapshots"])


def test_duplicate_snapshot_times_are_rejected(tmp_path: Path) -> None:
    case = make_case(tmp_path)
    (case / "intermediate/snapshot-0.500").write_text("duplicate\n")
    helper = make_helper(tmp_path)

    result = run_extractor(case, helper, tmp_path / "duplicate", 1)

    assert result.returncode != 0
    assert "Duplicate snapshot time" in result.stderr


def test_help_and_cpu_validation(tmp_path: Path) -> None:
    help_result = subprocess.run(
        [sys.executable, str(SCRIPT), "--help"], capture_output=True, text=True
    )
    case = make_case(tmp_path)
    helper = make_helper(tmp_path)
    cpu_result = run_extractor(case, helper, tmp_path / "cpu", 0)

    assert help_result.returncode == 0
    assert "--cpus" in help_result.stdout and "--CPUs" in help_result.stdout
    assert "--max-frames" in help_result.stdout and "--skip-video" in help_result.stdout
    assert cpu_result.returncode != 0
    assert "cpus must be positive" in cpu_result.stderr


def test_c_helper_uses_axis_intercept_not_broad_tip_band() -> None:
    source = HELPER_SOURCE.read_text()

    assert "fabs(r) <= tolerance" in source
    assert "#define R_TIP" not in source
    assert "kappa_tip = kappa[]" in source
    assert "closest valid curvature" not in source
