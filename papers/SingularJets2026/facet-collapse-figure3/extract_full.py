#!/usr/bin/env python3
"""Extract full interface facets from a simulation case.

The helper is run once per snapshot. Facets are staged by worker processes,
then installed and indexed deterministically by the parent process only after
every extraction succeeds.
"""

from __future__ import annotations

import argparse
import bisect
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


DEFAULT_RMAX = 0.6
DEFAULT_CPUS = 4
DEFAULT_MAX_LOG_GAP = 1.0e-3


def snapshot_time(path: Path) -> float:
    """Return the numeric suffix of a ``snapshot-<time>`` path."""
    try:
        time = float(path.name.rsplit("-", 1)[1])
    except (IndexError, ValueError) as exc:
        raise ValueError(f"Invalid snapshot name: {path}") from exc
    if not math.isfinite(time):
        raise ValueError(f"Non-finite snapshot time: {path}")
    return time


def atomic_write(path: Path, payload: str) -> None:
    """Atomically replace *path* with a non-empty text payload."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.name}.", dir=path.parent, text=True
    )
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        if temporary.stat().st_size == 0:
            raise RuntimeError(f"Refusing to install empty output: {path}")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def read_log(case: Path) -> list[tuple[float, float, float]]:
    """Read ``(time, jet radius, base height)`` rows from the case log."""
    rows: list[tuple[float, float, float]] = []
    with (case / "log").open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 9:
                continue
            try:
                row = (float(fields[2]), float(fields[7]), float(fields[8]))
            except ValueError:
                continue
            if not all(math.isfinite(value) for value in row):
                raise RuntimeError(f"Non-finite numeric row in {case / 'log'}")
            rows.append(row)
    rows.sort()
    if not rows:
        raise RuntimeError(f"No numeric rows found in {case / 'log'}")
    return rows


def nearest_log(
    time: float,
    rows: Sequence[tuple[float, float, float]],
    max_gap: float,
) -> tuple[float, float, float]:
    """Return the log row nearest to *time*."""
    times = [row[0] for row in rows]
    index = bisect.bisect_left(times, time)
    candidates = [
        candidate
        for candidate in (index - 1, index, index + 1)
        if 0 <= candidate < len(rows)
    ]
    nearest = min(candidates, key=lambda candidate: abs(times[candidate] - time))
    gap = abs(times[nearest] - time)
    if gap > max_gap:
        raise RuntimeError(
            f"Nearest log row is {gap:.8e} from snapshot t={time:.8e}; "
            f"maximum is {max_gap:.8e}"
        )
    return rows[nearest]


def file_sha256(path: Path) -> str:
    """Return a streaming SHA-256 digest for *path*."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_facet_record(lines: Sequence[str], snapshot: Path) -> tuple[float, ...]:
    """Parse one blank-line-delimited two-point facet record."""
    if len(lines) != 2:
        raise RuntimeError(f"Malformed facet record in {snapshot}: expected two lines")
    fields = [line.split() for line in lines]
    if any(len(point) != 2 for point in fields):
        raise RuntimeError(f"Malformed facet record in {snapshot}: expected two columns")
    try:
        z1, r1 = map(float, fields[0])
        z2, r2 = map(float, fields[1])
    except ValueError as exc:
        raise RuntimeError(f"Malformed numeric facet record in {snapshot}") from exc
    values = (z1, r1, z2, r2)
    if not all(math.isfinite(value) for value in values):
        raise RuntimeError(f"Non-finite facet coordinate in {snapshot}")
    if r1 < 0 or r2 < 0:
        raise RuntimeError(f"Negative facet radius in {snapshot}")
    return values


def filter_facets(raw: Path, destination: Path, snapshot: Path, rmax: float) -> None:
    """Stream, validate and atomically install a filtered facet payload."""
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.name}.", dir=destination.parent, text=True
    )
    temporary = Path(temporary_name)
    kept = 0
    record: list[str] = []
    try:
        with raw.open() as source, os.fdopen(descriptor, "w") as target:
            for line in source:
                stripped = line.strip()
                if stripped:
                    record.append(stripped)
                    continue
                if not record:
                    continue
                z1, r1, z2, r2 = parse_facet_record(record, snapshot)
                if r1 < rmax and r2 < rmax:
                    target.write(f"{z1} {r1}\n{z2} {r2}\n\n")
                    kept += 1
                record = []
            if record:
                z1, r1, z2, r2 = parse_facet_record(record, snapshot)
                if r1 < rmax and r2 < rmax:
                    target.write(f"{z1} {r1}\n{z2} {r2}\n\n")
                    kept += 1
            target.flush()
            os.fsync(target.fileno())
        if kept == 0:
            raise RuntimeError(f"No valid facets retained for {snapshot}")
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def extract_task(task: tuple[str, str, str, float]) -> tuple[str, float]:
    """Extract one snapshot into its unique staging path."""
    snapshot_text, helper_text, staged_text, rmax = task
    snapshot = Path(snapshot_text)
    staged = Path(staged_text)
    raw = staged.with_suffix(".raw")
    environment = os.environ.copy()
    for variable in (
        "OMP_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "MKL_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS",
        "NUMEXPR_NUM_THREADS",
    ):
        environment[variable] = "1"
    try:
        with raw.open("w") as stderr:
            subprocess.run(
                [helper_text, snapshot_text],
                stdout=subprocess.DEVNULL,
                stderr=stderr,
                text=True,
                check=True,
                env=environment,
            )
        if raw.stat().st_size == 0:
            raise RuntimeError(f"getFacet returned no facet payload for {snapshot}")
        filter_facets(raw, staged, snapshot, rmax)
    finally:
        raw.unlink(missing_ok=True)
    return staged_text, snapshot_time(snapshot)


def discover_snapshots(case: Path) -> list[Path]:
    """Find snapshots and reject duplicate times or output-name collisions."""
    snapshots = sorted(
        (Path(path) for path in glob.glob(str(case / "intermediate" / "snapshot-*"))),
        key=snapshot_time,
    )
    if not snapshots:
        raise RuntimeError(f"No snapshots found in {case / 'intermediate'}")

    seen_times: dict[float, Path] = {}
    seen_labels: dict[str, Path] = {}
    for snapshot in snapshots:
        time = snapshot_time(snapshot)
        label = f"{time:.6f}"
        if time in seen_times:
            raise RuntimeError(
                f"Duplicate snapshot time {time}: {seen_times[time]} and {snapshot}"
            )
        if label in seen_labels:
            raise RuntimeError(
                "Snapshot times collide at six-decimal output precision: "
                f"{seen_labels[label]} and {snapshot}"
            )
        seen_times[time] = snapshot
        seen_labels[label] = snapshot
    return snapshots


def select_tasks(
    snapshots: Sequence[Path], t0: float, pre_lo: float, max_frames: int | None
) -> list[tuple[Path, str, str]]:
    """Select pre/post snapshots and return one combined deterministic task list."""
    at_t0 = [path for path in snapshots if snapshot_time(path) == t0]
    if at_t0:
        raise RuntimeError(
            f"Snapshot time equals t0 ({t0}); inception snapshots are neither pre nor post"
        )
    pre = [path for path in snapshots if pre_lo <= snapshot_time(path) < t0]
    post = [path for path in snapshots if snapshot_time(path) > t0]
    if max_frames is not None:
        pre = pre[:max_frames]
        post = post[:max_frames]
    if not pre:
        raise RuntimeError("No snapshots selected for facetpremain")
    if not post:
        raise RuntimeError("No snapshots selected for facetmain")
    return [
        *((path, "facetpremain", "index_pre.txt") for path in pre),
        *((path, "facetmain", "index.txt") for path in post),
    ]


def run_extraction(
    case: Path,
    out: Path,
    getfacet: Path,
    t0: float,
    pre_lo: float,
    rmax: float,
    cpus: int,
    max_frames: int | None,
    max_log_gap: float,
    force: bool,
) -> tuple[int, int]:
    """Extract a case and return the pre/post frame counts."""
    if cpus <= 0:
        raise ValueError("cpus must be positive")
    if max_frames is not None and max_frames <= 0:
        raise ValueError("max-frames must be positive")
    if pre_lo >= t0:
        raise ValueError("pre-lo must be less than t0")
    if rmax <= 0:
        raise ValueError("rmax must be positive")
    if not all(math.isfinite(value) for value in (t0, pre_lo, rmax, max_log_gap)):
        raise ValueError("time, radius and log-gap arguments must be finite")
    if max_log_gap <= 0:
        raise ValueError("max-log-gap must be positive")
    if not case.is_dir():
        raise FileNotFoundError(f"Case directory does not exist: {case}")
    if not getfacet.is_file() or not os.access(getfacet, os.X_OK):
        raise FileNotFoundError(f"getFacet helper is not executable: {getfacet}")

    snapshots = discover_snapshots(case)
    selected = select_tasks(snapshots, t0, pre_lo, max_frames)
    log_rows = read_log(case)
    out.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix=".extract-full-", dir=out))
    try:
        worker_tasks: list[tuple[str, str, str, float]] = []
        destinations: list[Path] = []
        for snapshot, prefix, _ in selected:
            destination = out / f"{prefix}_{snapshot_time(snapshot):.6f}.txt"
            staged = staging / destination.name
            destinations.append(destination)
            worker_tasks.append((str(snapshot), str(getfacet), str(staged), rmax))

        final_paths = [*destinations, out / "index_pre.txt", out / "index.txt", out / "manifest.json"]
        managed_existing = {
            *out.glob("facetpremain_*.txt"),
            *out.glob("facetmain_*.txt"),
            *(path for path in final_paths[-3:] if path.exists()),
        }
        existing = sorted(managed_existing)
        if existing and not force:
            names = ", ".join(path.name for path in existing[:5])
            raise FileExistsError(f"Refusing to overwrite existing outputs: {names}; use --force")
        stale_outputs = sorted(managed_existing - set(final_paths))

        context = multiprocessing.get_context("spawn")
        with context.Pool(processes=min(cpus, len(worker_tasks))) as pool:
            results = pool.map(extract_task, worker_tasks, chunksize=1)

        index_rows: dict[str, list[str]] = {
            "index_pre.txt": [f"t0 {t0:.6f}\n"],
            "index.txt": [f"t0 {t0:.6f}\n"],
        }
        for (snapshot, _, index_name), (_, time) in zip(selected, results, strict=True):
            if time != snapshot_time(snapshot):
                raise RuntimeError(f"Worker returned the wrong time for {snapshot}")
            _, radius, base = nearest_log(time, log_rows, max_log_gap)
            index_rows[index_name].append(f"{time:.6f} {radius:.8e} {base:.8e}\n")
        for index_name in ("index_pre.txt", "index.txt"):
            atomic_write(staging / index_name, "".join(index_rows[index_name]))

        manifest = {
            "schema_version": 1,
            "configuration": {
                "max_log_gap": max_log_gap,
                "pre_lo": pre_lo,
                "rmax": rmax,
                "t0": t0,
            },
            "helper_sha256": file_sha256(getfacet),
            "snapshots": [
                {
                    "file": snapshot.name,
                    "output": destination.name,
                    "phase": prefix,
                    "time": snapshot_time(snapshot),
                    "log_time": nearest_log(snapshot_time(snapshot), log_rows, max_log_gap)[0],
                }
                for (snapshot, prefix, _), destination in zip(
                    selected, destinations, strict=True
                )
            ],
        }
        atomic_write(
            staging / "manifest.json", json.dumps(manifest, indent=2, sort_keys=True) + "\n"
        )

        staged_paths = [
            *(Path(staged_text) for staged_text, _ in results),
            staging / "index_pre.txt",
            staging / "index.txt",
            staging / "manifest.json",
        ]
        if force and final_paths[-1].exists():
            final_paths[-1].unlink()
        for staged_path, final_path in zip(staged_paths[:-1], final_paths[:-1], strict=True):
            os.replace(staged_path, final_path)
        for stale_path in stale_outputs:
            stale_path.unlink()
        # The manifest is the completion receipt and is installed only after
        # every selected output is current and every stale managed output is gone.
        os.replace(staged_paths[-1], final_paths[-1])
    finally:
        shutil.rmtree(staging, ignore_errors=True)

    pre_count = sum(prefix == "facetpremain" for _, prefix, _ in selected)
    post_count = len(selected) - pre_count
    return pre_count, post_count


def build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True, help="simulation case directory")
    parser.add_argument("--out", type=Path, help="output directory (default: CASE/facets_full)")
    parser.add_argument("--getfacet", type=Path, help="helper executable (default: CASE/getFacet)")
    parser.add_argument("--t0", type=float, required=True, help="jet-inception time")
    parser.add_argument("--pre-lo", type=float, help="lower pre-inception time (default: T0-0.05)")
    parser.add_argument("--rmax", type=float, default=DEFAULT_RMAX, help="maximum retained radius")
    parser.add_argument(
        "--cpus", "--CPUs", type=int, default=DEFAULT_CPUS, help="worker processes (default: 4)"
    )
    parser.add_argument(
        "--max-frames", type=int, help="maximum snapshots to process in each pre/post phase"
    )
    parser.add_argument(
        "--max-log-gap",
        type=float,
        default=DEFAULT_MAX_LOG_GAP,
        help="maximum snapshot-to-log time gap (default: 1e-3)",
    )
    parser.add_argument("--force", action="store_true", help="replace matching existing outputs")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Command-line entry point."""
    args = build_parser().parse_args(argv)
    case = args.case.resolve()
    out = (args.out if args.out is not None else case / "facets_full").resolve()
    getfacet = (
        args.getfacet if args.getfacet is not None else case / "getFacet"
    ).resolve()
    pre_lo = args.pre_lo if args.pre_lo is not None else args.t0 - 0.05
    pre_count, post_count = run_extraction(
        case,
        out,
        getfacet,
        args.t0,
        pre_lo,
        args.rmax,
        args.cpus,
        args.max_frames,
        args.max_log_gap,
        args.force,
    )
    print(f"DONE pre={pre_count} post={post_count} -> {out}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
