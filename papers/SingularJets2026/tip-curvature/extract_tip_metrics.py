#!/usr/bin/env python3
r"""
# Parallel extraction of incipient-jet tip metrics

Run `getTipMetrics` once per Basilisk snapshot and install one deterministic
CSV plus a provenance manifest. Worker processes share one precompiled helper;
no compilation or plotting occurs inside the snapshot loop.

The raw helper curvature is the axisymmetric mean curvature. Derived columns
use the locally spherical/paraboloidal apex relation
$R_\kappa=2/|\kappa|$ and report $R_\kappa/\Delta$ explicitly. The output is a
resolution diagnostic, not proof that $R_\kappa$ equals the theoretical cutoff
radius $R_m$.
"""

from __future__ import annotations

import argparse
import csv
import glob
import hashlib
import json
import math
import multiprocessing
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import Sequence


DEFAULT_CPUS = 4
RAW_COLUMNS = (
    "time",
    "z_tip",
    "r_tip",
    "z_cell",
    "r_cell",
    "kappa_mean",
    "u_z_tip",
    "u_r_tip",
    "speed_tip",
    "delta_tip",
    "level_tip",
    "f_tip",
    "liquid_components",
    "curvature_height",
    "curvature_fit",
    "curvature_average",
    "curvature_centroid",
)
INTEGER_COLUMNS = {
    "level_tip",
    "liquid_components",
    "curvature_height",
    "curvature_fit",
    "curvature_average",
    "curvature_centroid",
}
OUTPUT_COLUMNS = (
    "snapshot",
    *RAW_COLUMNS,
    "inverse_mean_curvature",
    "apex_radius",
    "kappa_delta",
    "apex_radius_cells",
    "we_apex_uz",
    "we_apex_speed",
    "tip_cell_offset",
    "tip_cell_offset_cells",
)


def snapshot_time(path: Path) -> float:
    """Return the finite numeric suffix of a `snapshot-<time>` path."""
    try:
        value = float(path.name.rsplit("-", 1)[1])
    except (IndexError, ValueError) as exc:
        raise ValueError(f"Invalid snapshot name: {path}") from exc
    if not math.isfinite(value):
        raise ValueError(f"Non-finite snapshot time: {path}")
    return value


def file_sha256(path: Path) -> str:
    """Return a streaming SHA-256 digest."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def discover_snapshots(
    case: Path,
    time_min: float | None,
    time_max: float | None,
    max_frames: int | None,
) -> list[Path]:
    """Discover, filter and deterministically order case snapshots."""
    snapshots = sorted(
        (Path(item) for item in glob.glob(str(case / "intermediate" / "snapshot-*"))),
        key=snapshot_time,
    )
    selected: list[Path] = []
    seen: set[float] = set()
    for snapshot in snapshots:
        time = snapshot_time(snapshot)
        if time in seen:
            raise ValueError(f"Duplicate snapshot time {time}: {snapshot}")
        seen.add(time)
        if time_min is not None and time < time_min:
            continue
        if time_max is not None and time > time_max:
            continue
        selected.append(snapshot)
    if max_frames is not None:
        selected = selected[:max_frames]
    if not selected:
        raise RuntimeError(f"No snapshots selected from {case / 'intermediate'}")
    return selected


def parse_helper_output(stderr: str, snapshot: Path) -> dict[str, float | int | str]:
    """Parse the unique `TIP_METRICS` line emitted by the helper."""
    records = [line.split() for line in stderr.splitlines() if line.startswith("TIP_METRICS ")]
    if len(records) != 1:
        raise RuntimeError(
            f"Expected one TIP_METRICS row for {snapshot}, found {len(records)}"
        )
    fields = records[0][1:]
    if len(fields) != len(RAW_COLUMNS):
        raise RuntimeError(
            f"Malformed TIP_METRICS row for {snapshot}: expected "
            f"{len(RAW_COLUMNS)} values, found {len(fields)}"
        )
    row: dict[str, float | int | str] = {"snapshot": snapshot.name}
    for name, field in zip(RAW_COLUMNS, fields, strict=True):
        try:
            value = int(field) if name in INTEGER_COLUMNS else float(field)
        except ValueError as exc:
            raise RuntimeError(f"Non-numeric {name} for {snapshot}: {field}") from exc
        if isinstance(value, float) and not math.isfinite(value):
            raise RuntimeError(f"Non-finite {name} for {snapshot}: {field}")
        row[name] = value

    kappa = abs(float(row["kappa_mean"]))
    delta = float(row["delta_tip"])
    if kappa <= 0.0 or delta <= 0.0:
        raise RuntimeError(f"Nonpositive curvature magnitude or grid spacing for {snapshot}")
    inverse = 1.0 / kappa
    apex_radius = 2.0 * inverse
    uz_tip = float(row["u_z_tip"])
    speed_tip = float(row["speed_tip"])
    row.update(
        {
            "inverse_mean_curvature": inverse,
            "apex_radius": apex_radius,
            "kappa_delta": kappa * delta,
            "apex_radius_cells": apex_radius / delta,
            "we_apex_uz": uz_tip**2 * apex_radius,
            "we_apex_speed": speed_tip**2 * apex_radius,
            "tip_cell_offset": math.hypot(
                float(row["z_cell"]) - float(row["z_tip"]),
                float(row["r_cell"]) - float(row["r_tip"]),
            ),
        }
    )
    row["tip_cell_offset_cells"] = float(row["tip_cell_offset"]) / delta
    return row


def extract_task(task: tuple[str, str]) -> dict[str, float | int | str]:
    """Run the serial helper for one snapshot in an isolated worker."""
    snapshot_text, helper_text = task
    snapshot = Path(snapshot_text)
    environment = os.environ.copy()
    for variable in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        environment[variable] = "1"
    result = subprocess.run(
        [helper_text, snapshot_text],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
        env=environment,
    )
    if result.returncode:
        raise RuntimeError(
            f"getTipMetrics failed for {snapshot} with rc={result.returncode}: "
            f"{result.stderr[-1000:]}"
        )
    row = parse_helper_output(result.stderr, snapshot)
    named_time = snapshot_time(snapshot)
    if not math.isclose(float(row["time"]), named_time, rel_tol=0.0, abs_tol=5.1e-7):
        raise RuntimeError(
            f"Snapshot/helper time mismatch for {snapshot}: {row['time']}"
        )
    return row


def atomic_write_text(path: Path, payload: str) -> None:
    """Atomically replace `path` with a non-empty text payload."""
    if not payload:
        raise ValueError(f"Refusing to write empty output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent, text=True
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def csv_payload(rows: Sequence[dict[str, float | int | str]]) -> str:
    """Serialize rows with stable column and row ordering."""
    import io

    handle = io.StringIO(newline="")
    writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, lineterminator="\n")
    writer.writeheader()
    for row in sorted(rows, key=lambda item: float(item["time"])):
        writer.writerow(row)
    return handle.getvalue()


def extract_case(
    case: Path,
    output: Path,
    helper: Path,
    cpus: int,
    time_min: float | None,
    time_max: float | None,
    max_frames: int | None,
    force: bool,
    provenance: dict[str, object],
) -> int:
    """Extract a case and return the installed row count."""
    if cpus <= 0:
        raise ValueError("cpus must be positive")
    if max_frames is not None and max_frames <= 0:
        raise ValueError("max-frames must be positive")
    for name, value in (("time-min", time_min), ("time-max", time_max)):
        if value is not None and not math.isfinite(value):
            raise ValueError(f"{name} must be finite")
    if time_min is not None and time_max is not None and time_min > time_max:
        raise ValueError("time-min must not exceed time-max")
    if not case.is_dir():
        raise FileNotFoundError(f"Case directory does not exist: {case}")
    if not helper.is_file() or not os.access(helper, os.X_OK):
        raise FileNotFoundError(f"getTipMetrics helper is not executable: {helper}")

    csv_path = output / "tip_metrics.csv"
    manifest_path = output / "manifest.json"
    if (csv_path.exists() or manifest_path.exists()) and not force:
        raise FileExistsError(f"Refusing to overwrite {output}; use --force")

    snapshots = discover_snapshots(case, time_min, time_max, max_frames)
    tasks = [(str(snapshot), str(helper)) for snapshot in snapshots]
    context = multiprocessing.get_context("spawn")
    with context.Pool(processes=min(cpus, len(tasks))) as pool:
        rows = pool.map(extract_task, tasks, chunksize=1)

    output.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".tip-metrics-", dir=output))
    try:
        staged_csv = staging / csv_path.name
        staged_manifest = staging / manifest_path.name
        atomic_write_text(staged_csv, csv_payload(rows))
        manifest = {
            "schema_version": 1,
            "case": str(case),
            "configuration": {
                "max_frames": max_frames,
                "time_max": time_max,
                "time_min": time_min,
            },
            "definition": {
                "apex_radius": "2/abs(axisymmetric mean curvature)",
                "tip": "highest near-axis facet endpoint on largest liquid component",
                "we_apex_speed": "speed_tip**2 * apex_radius",
                "we_apex_uz": "u_z_tip**2 * apex_radius",
            },
            "provenance": provenance,
            "helper": str(helper),
            "helper_sha256": file_sha256(helper),
            "row_count": len(rows),
            "snapshots": [
                {
                    "name": snapshot.name,
                    "size": snapshot.stat().st_size,
                    "mtime_ns": snapshot.stat().st_mtime_ns,
                }
                for snapshot in snapshots
            ],
        }
        atomic_write_text(
            staged_manifest, json.dumps(manifest, indent=2, sort_keys=True) + "\n"
        )
        os.replace(staged_csv, csv_path)
        # The manifest is the completion receipt and is installed last.
        os.replace(staged_manifest, manifest_path)
    finally:
        shutil.rmtree(staging, ignore_errors=True)
    return len(rows)


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True, help="simulation case directory")
    parser.add_argument("--out", type=Path, required=True, help="output directory")
    parser.add_argument("--gettip", type=Path, required=True, help="compiled getTipMetrics helper")
    parser.add_argument(
        "--cpus", "--CPUs", type=int, default=DEFAULT_CPUS, help="worker processes (default: 4)"
    )
    parser.add_argument("--time-min", type=float, help="inclusive minimum snapshot time")
    parser.add_argument("--time-max", type=float, help="inclusive maximum snapshot time")
    parser.add_argument("--oh", type=float, required=True)
    parser.add_argument("--bond", type=float, required=True)
    parser.add_argument("--pre-level", type=int, required=True)
    parser.add_argument("--post-level", type=int, required=True)
    parser.add_argument("--t0", type=float, required=True, help="recorded jet-inception time")
    parser.add_argument("--t0-protocol", required=True, help="brief inception-time definition")
    parser.add_argument("--lineage", required=True, help="case/restart lineage label")
    parser.add_argument("--bridged", action="store_true", help="mark a lower-level topology bridge")
    parser.add_argument("--max-frames", type=int, help="maximum selected snapshots")
    parser.add_argument("--force", action="store_true", help="replace managed outputs")
    parser.add_argument(
        "--skip-video",
        action="store_true",
        help="accepted for the common snapshot-smoke contract; this extractor never renders video",
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry point."""
    args = build_parser().parse_args(argv)
    count = extract_case(
        args.case.resolve(),
        args.out.resolve(),
        args.gettip.resolve(),
        args.cpus,
        args.time_min,
        args.time_max,
        args.max_frames,
        args.force,
        {
            "bond": args.bond,
            "bridged": args.bridged,
            "lineage": args.lineage,
            "oh": args.oh,
            "post_level": args.post_level,
            "pre_level": args.pre_level,
            "t0": args.t0,
            "t0_protocol": args.t0_protocol,
        },
    )
    print(f"TIP_METRICS_COMPLETE rows={count} out={args.out.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
