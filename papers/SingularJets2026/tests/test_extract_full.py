"""Focused regression tests for the full-facet extractor."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys


PAPER_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = PAPER_ROOT / "facet-collapse-figure3/extract_full.py"


def make_case(tmp_path: Path) -> Path:
    case = tmp_path / "case with spaces;literal-$value"
    intermediate = case / "intermediate"
    intermediate.mkdir(parents=True)
    for time in (0.10, 0.15, 0.30, 0.35):
        (intermediate / f"snapshot-{time:.2f}").write_text("snapshot\n")
    (case / "log").write_text(
        "0 0 0.10 0 0 0 0 1.0 2.0\n"
        "0 0 0.15 0 0 0 0 1.5 2.5\n"
        "0 0 0.30 0 0 0 0 3.0 4.0\n"
        "0 0 0.35 0 0 0 0 3.5 4.5\n"
    )
    return case


def make_helper(tmp_path: Path, mode: str = "ok") -> Path:
    helper = tmp_path / f"getFacet {mode}"
    action = {
        "ok": "",
        "fail": "if '0.30' in sys.argv[1]: raise SystemExit(7)",
        "malformed": "if '0.30' in sys.argv[1]: sys.stderr.write('0 0\\n'); raise SystemExit(0)",
        "nan": "if '0.30' in sys.argv[1]: sys.stderr.write('0 nan\\n1 0.2\\n\\n'); raise SystemExit(0)",
    }[mode]
    helper.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        f"{action}\n"
        "sys.stderr.write('0.0 0.1\\n1.0 0.2\\n\\n2.0 0.7\\n3.0 0.8\\n\\n')\n"
    )
    helper.chmod(0o755)
    return helper


def run_extractor(
    case: Path, helper: Path, out: Path, cpus: int, *extra: str
) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = "192"
    return subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "--case",
            str(case),
            "--out",
            str(out),
            "--getfacet",
            str(helper),
            "--t0",
            "0.20",
            "--pre-lo",
            "0.05",
            "--cpus",
            str(cpus),
            *extra,
        ],
        capture_output=True,
        text=True,
        env=environment,
        timeout=60,
    )


def directory_bytes(path: Path) -> dict[str, bytes]:
    return {
        item.relative_to(path).as_posix(): item.read_bytes()
        for item in sorted(path.rglob("*"))
        if item.is_file()
    }


def test_serial_and_parallel_outputs_are_byte_identical(tmp_path: Path) -> None:
    case = make_case(tmp_path)
    helper = make_helper(tmp_path)
    serial = tmp_path / "serial output"
    parallel = tmp_path / "parallel output"

    serial_result = run_extractor(case, helper, serial, 1)
    parallel_result = run_extractor(case, helper, parallel, 4)

    assert serial_result.returncode == 0, serial_result.stderr
    assert parallel_result.returncode == 0, parallel_result.stderr
    assert directory_bytes(serial) == directory_bytes(parallel)
    assert len(directory_bytes(serial)) == 7


def test_worker_failures_do_not_publish_partial_results(tmp_path: Path) -> None:
    case = make_case(tmp_path)
    for mode in ("fail", "malformed", "nan"):
        helper = make_helper(tmp_path, mode)
        out = tmp_path / f"facets {mode}"

        result = run_extractor(case, helper, out, 4)

        assert result.returncode != 0
        assert directory_bytes(out) == {}


def test_refuses_to_overwrite_without_force(tmp_path: Path) -> None:
    case = make_case(tmp_path)
    helper = make_helper(tmp_path)
    out = tmp_path / "facets"
    first = run_extractor(case, helper, out, 1)
    before = directory_bytes(out)

    second = run_extractor(case, helper, out, 4)

    assert first.returncode == 0, first.stderr
    assert second.returncode != 0
    assert "use --force" in second.stderr
    assert directory_bytes(out) == before


def test_force_removes_stale_managed_facets_after_success(tmp_path: Path) -> None:
    case = make_case(tmp_path)
    helper = make_helper(tmp_path)
    out = tmp_path / "facets"
    first = run_extractor(case, helper, out, 1)
    stale = out / "facetmain_9.999999.txt"
    stale.write_text("stale\n", encoding="utf-8")

    forced = run_extractor(case, helper, out, 4, "--force")

    assert first.returncode == 0, first.stderr
    assert forced.returncode == 0, forced.stderr
    assert not stale.exists()
    manifest = (out / "manifest.json").read_text(encoding="utf-8")
    assert "9.999999" not in manifest


def test_rejects_zero_duplicate_and_t0_snapshots(tmp_path: Path) -> None:
    case = make_case(tmp_path)
    helper = make_helper(tmp_path)
    for snapshot in (case / "intermediate").iterdir():
        snapshot.unlink()
    empty_result = run_extractor(case, helper, tmp_path / "empty", 1)
    assert empty_result.returncode != 0
    assert "No snapshots found" in empty_result.stderr

    (case / "intermediate" / "snapshot-0.1").write_text("snapshot\n")
    (case / "intermediate" / "snapshot-0.10").write_text("snapshot\n")
    duplicate_result = run_extractor(case, helper, tmp_path / "duplicate", 1)
    assert duplicate_result.returncode != 0
    assert "Duplicate snapshot time" in duplicate_result.stderr

    (case / "intermediate" / "snapshot-0.10").unlink()
    (case / "intermediate" / "snapshot-0.20").write_text("snapshot\n")
    t0_result = run_extractor(case, helper, tmp_path / "at t0", 1)
    assert t0_result.returncode != 0
    assert "equals t0" in t0_result.stderr


def test_validates_cpu_and_log_gap_arguments(tmp_path: Path) -> None:
    case = make_case(tmp_path)
    helper = make_helper(tmp_path)

    cpu_result = run_extractor(case, helper, tmp_path / "cpu", 0)
    (case / "log").write_text(
        (case / "log").read_text().replace("0 0 0.10 ", "0 0 0.11 ")
    )
    gap_result = run_extractor(
        case, helper, tmp_path / "gap", 1, "--max-log-gap", "1e-9"
    )

    assert cpu_result.returncode != 0
    assert "cpus must be positive" in cpu_result.stderr
    assert gap_result.returncode != 0
    assert "Nearest log row" in gap_result.stderr


def test_help_exposes_both_cpu_spellings() -> None:
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--help"], capture_output=True, text=True
    )

    assert result.returncode == 0
    assert "--cpus" in result.stdout
    assert "--CPUs" in result.stdout
