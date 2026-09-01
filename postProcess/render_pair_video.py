#!/usr/bin/env python3
"""Render interface/jet-base and interface/AMR video frames in parallel.

Each snapshot is independent. Workers restore one Basilisk dump, extract the
interface and jet-base position, render the adaptive mesh into process-local
scratch, and atomically install one deterministic PNG frame. Encoding is
optional so an expensive full-node frame job can be followed by a small
dependency job for ``ffmpeg``.

The positional interface remains compatible with the historical script::

    render_pair_video.py CASE_DIR LDOMAIN OUT [ZBOT ZTOP RMAX FPS]

Use ``--cpus``, ``--frames-dir`` and ``--skip-video-encode`` for production
batch work.
"""

from __future__ import annotations

import argparse
import bisect
import glob
import math
import multiprocessing as mp
import os
import subprocess as sp
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

# Spawned workers import this module afresh. Assign process-specific caches
# before importing matplotlib and pin nested numerical libraries to one thread.
os.environ["MPLCONFIGDIR"] = os.path.join(tempfile.gettempdir(), f"mpl_{os.getpid()}")
os.environ["TEXMFVAR"] = os.path.join(tempfile.gettempdir(), f"texmf-var_{os.getpid()}")
os.environ["TEXMFCONFIG"] = os.path.join(tempfile.gettempdir(), f"texmf-config_{os.getpid()}")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import matplotlib

matplotlib.use("Agg")
matplotlib.rcParams["font.family"] = "serif"
matplotlib.rcParams["mathtext.fontset"] = "cm"
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from PIL import Image, ImageDraw, ImageFont

GREEN = (0.0, 0.5, 0.0)


@dataclass(frozen=True)
class RenderConfig:
    """Immutable settings shared by spawned rendering workers."""

    case_dir: str
    helper_dir: str
    frames_dir: str
    ldomain: float
    zbot: float
    ztop: float
    rmax: float


@dataclass(frozen=True)
class FrameTask:
    """One deterministic snapshot-to-frame rendering task."""

    index: int
    snapshot: str
    time: float
    kinetic_energy: float
    maxlevel: int
    target: str


_CONFIG: RenderConfig | None = None


def parse_arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse the backward-compatible positional and production options."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_dir")
    parser.add_argument("ldomain", type=float)
    parser.add_argument("out")
    parser.add_argument("zbot", nargs="?", type=float, default=-2.0)
    parser.add_argument("ztop", nargs="?", type=float, default=1.5)
    parser.add_argument("rmax", nargs="?", type=float, default=1.1)
    parser.add_argument("fps", nargs="?", type=int, default=18)
    parser.add_argument("--cpus", "--CPUs", type=int, default=1, dest="cpus")
    parser.add_argument("--frames-dir", default=None)
    parser.add_argument("--helper-dir", default=None)
    parser.add_argument("--max-frames", type=int, default=0)
    parser.add_argument("--skip-video-encode", action="store_true")
    args = parser.parse_args(argv)
    if args.cpus < 1:
        parser.error("--cpus/--CPUs must be positive")
    if args.max_frames < 0:
        parser.error("--max-frames cannot be negative")
    if args.ldomain <= 0 or args.ztop <= args.zbot or args.rmax <= 0:
        parser.error("invalid physical rendering bounds")
    return args


def run_helper(command: Sequence[str], case_dir: str) -> list[str]:
    """Run one compiled helper and return its stderr payload as lines."""

    process = sp.run(command, cwd=case_dir, capture_output=True, text=True)
    if process.returncode != 0:
        detail = process.stderr.strip() or process.stdout.strip()
        raise RuntimeError(f"{' '.join(command)} failed ({process.returncode}): {detail}")
    return process.stderr.splitlines()


def facets(helper: str, relative_snapshot: str, case_dir: str):
    """Extract and mirror the interface line segments."""

    lines = run_helper([helper, relative_snapshot], case_dir)
    segments = []
    skip = False
    for index, line in enumerate(lines):
        fields = line.split()
        if not fields:
            skip = False
            continue
        if skip or index + 1 >= len(lines):
            continue
        following = lines[index + 1].split()
        try:
            z1, r1 = float(fields[0]), float(fields[1])
            z2, r2 = float(following[0]), float(following[1])
        except (IndexError, ValueError):
            continue
        segments.extend([((r1, z1), (r2, z2)), ((-r1, z1), (-r2, z2))])
        skip = True
    return segments


def jet_base(helper: str, relative_snapshot: str, case_dir: str) -> tuple[float, float]:
    """Return the axial and radial jet-base coordinates, when present."""

    for line in run_helper([helper, relative_snapshot], case_dir):
        fields = line.split()
        if len(fields) < 6:
            continue
        try:
            return float(fields[1]), float(fields[2])
        except ValueError:
            continue
    return -1000.0, -1000.0


def read_log(case_dir: str) -> tuple[list[float], list[tuple[float, int]]]:
    """Read time, kinetic energy and maxlevel from the solver log."""

    rows = []
    with open(os.path.join(case_dir, "log"), encoding="utf-8") as handle:
        for line in handle:
            columns = line.split()
            if not columns or not columns[0].lstrip("-").isdigit():
                continue
            rows.append((float(columns[2]), float(columns[3]), int(columns[4])))
    rows.sort()
    return [row[0] for row in rows], [(row[1], row[2]) for row in rows]


def nearest_log_value(
    time: float, log_times: Sequence[float], log_values: Sequence[tuple[float, int]]
) -> tuple[float, int]:
    """Return the nearest log annotation in logarithmic search time."""

    if not log_times:
        return float("nan"), -1
    insertion = bisect.bisect_left(log_times, time)
    candidates = [max(0, insertion - 1), min(len(log_times) - 1, insertion)]
    selected = min(candidates, key=lambda idx: abs(log_times[idx] - time))
    return log_values[selected]


def discover_tasks(case_dir: str, frames_dir: str, max_frames: int) -> list[FrameTask]:
    """Validate snapshots and create stable output indices."""

    snapshots = sorted(
        glob.glob(os.path.join(case_dir, "intermediate", "snapshot-*")),
        key=lambda path: float(path.rsplit("snapshot-", 1)[-1]),
    )
    if max_frames:
        snapshots = snapshots[:max_frames]
    if not snapshots:
        raise RuntimeError(f"No snapshots found below {case_dir}/intermediate")
    times = [float(path.rsplit("snapshot-", 1)[-1]) for path in snapshots]
    if len(times) != len(set(times)):
        raise RuntimeError("Duplicate snapshot times would overwrite deterministic frames")
    empty = [path for path in snapshots if os.path.getsize(path) == 0]
    if empty:
        raise RuntimeError(f"Zero-byte snapshot found: {empty[0]}")
    log_times, log_values = read_log(case_dir)
    tasks = []
    for index, (path, time) in enumerate(zip(snapshots, times)):
        kinetic_energy, maxlevel = nearest_log_value(time, log_times, log_values)
        tasks.append(
            FrameTask(
                index=index,
                snapshot=path,
                time=time,
                kinetic_energy=kinetic_energy,
                maxlevel=maxlevel,
                target=os.path.join(frames_dir, f"{index:06d}.png"),
            )
        )
    return tasks


def initialise_worker(config: RenderConfig) -> None:
    """Install immutable worker state after process spawn."""

    global _CONFIG
    _CONFIG = config


def render_frame(task: FrameTask) -> str:
    """Render one frame atomically with process-local mesh scratch."""

    if _CONFIG is None:
        raise RuntimeError("worker configuration is unavailable")
    config = _CONFIG
    if os.path.isfile(task.target) and os.path.getsize(task.target) > 0:
        return task.target

    relative = os.path.join("intermediate", os.path.basename(task.snapshot))
    get_facet = os.path.join(config.helper_dir, "getFacet")
    get_base = os.path.join(config.helper_dir, "getBase")
    get_view = os.path.join(config.helper_dir, "getView2D")
    segments = facets(get_facet, relative, config.case_dir)
    zbase, rbase = jet_base(get_base, relative, config.case_dir)

    zcentre = 0.5 * (config.zbot + config.ztop)
    height = config.ztop - config.zbot
    fov = math.degrees(2.0 * math.atan(config.rmax / (3.0 * config.ldomain)))
    width = 1100
    pixels_high = max(200, int(round(width * (2.0 * config.rmax) / height)))
    raw = os.path.join(tempfile.gettempdir(), f"sbbv_{os.getpid()}_{task.index}.png")
    grid = None
    try:
        for _ in range(3):
            if os.path.exists(raw):
                os.remove(raw)
            process = sp.run(
                [
                    get_view,
                    relative,
                    raw,
                    f"{fov:.4f}",
                    f"{-zcentre / config.ldomain:.6f}",
                    "0",
                    str(width),
                    str(pixels_high),
                    "0",
                ],
                cwd=config.case_dir,
                stdout=sp.DEVNULL,
                stderr=sp.DEVNULL,
            )
            if process.returncode == 0 and os.path.isfile(raw):
                with Image.open(raw) as image:
                    grid = image.rotate(90, expand=True).copy()
                break
        if grid is None:
            raise RuntimeError(f"getView2D failed for {relative}")

        figure, (left, right) = plt.subplots(1, 2, figsize=(11, 8.5))
        if segments:
            left.add_collection(LineCollection(segments, linewidths=2.0, colors=[GREEN]))
        left.plot([0, 0], [config.zbot, config.ztop], "-.", color="grey", lw=0.7)
        if rbase > -900:
            left.plot([rbase, -rbase], [zbase, zbase], "o", color="red", ms=11, zorder=10)
        left.set(xlim=(-config.rmax, config.rmax), ylim=(config.zbot, config.ztop))
        left.set_aspect("equal")
        left.axis("off")
        right.imshow(
            grid,
            extent=[-config.rmax, config.rmax, config.zbot, config.ztop],
            aspect="equal",
            origin="upper",
        )
        right.set(xlim=(-config.rmax, config.rmax), ylim=(config.zbot, config.ztop))
        right.axis("off")
        # Reserve a stable text band. Labels are added with Pillow after the
        # Matplotlib render because concurrent Matplotlib text can lose glyphs.
        figure.subplots_adjust(left=0.035, right=0.965, bottom=0.045, top=0.82, wspace=0.14)
        temporary = f"{task.target}.tmp.{os.getpid()}"
        try:
            figure.savefig(temporary, format="png", dpi=120)
            font_path = os.environ.get("VIDEO_FONT_PATH", "DejaVuSerif.ttf")
            title_font = ImageFont.truetype(font_path, 28)
            panel_font = ImageFont.truetype(font_path, 23)
            with Image.open(temporary) as rendered:
                annotated = rendered.convert("RGB")
            draw = ImageDraw.Draw(annotated)
            width, _ = annotated.size
            draw.text(
                (width // 2, 22),
                f"t/τγ = {task.time:.4f}    ke = {task.kinetic_energy:.3f}"
                f"    maxlevel = {task.maxlevel}",
                fill="black",
                font=title_font,
                anchor="ma",
            )
            draw.text(
                (int(0.27 * width), 80),
                "interface + jet-base marker",
                fill="black",
                font=panel_font,
                anchor="ma",
            )
            draw.text(
                (int(0.73 * width), 80),
                "interface + adaptive mesh",
                fill="black",
                font=panel_font,
                anchor="ma",
            )
            annotated.save(temporary, format="PNG")
            os.replace(temporary, task.target)
        finally:
            plt.close(figure)
            if os.path.exists(temporary):
                os.remove(temporary)
    finally:
        if os.path.exists(raw):
            os.remove(raw)
    return task.target


def encode_video(frames_dir: str, fps: int, output: str) -> None:
    """Encode deterministic frames into one H.264 video atomically."""

    destination = Path(output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(f".{destination.name}.tmp.{os.getpid()}.mp4")
    command = [
        "ffmpeg",
        "-y",
        "-framerate",
        str(fps),
        "-i",
        os.path.join(frames_dir, "%06d.png"),
        "-vf",
        "pad=ceil(iw/2)*2:ceil(ih/2)*2",
        "-c:v",
        "libx264",
        "-pix_fmt",
        "yuv420p",
        "-crf",
        "20",
        str(temporary),
    ]
    sp.run(command, check=True)
    os.replace(temporary, destination)


def main(argv: Sequence[str] | None = None) -> int:
    """Render all requested frames and optionally encode the MP4."""

    args = parse_arguments(argv)
    case_dir = os.path.abspath(args.case_dir)
    helper_dir = os.path.abspath(args.helper_dir or os.path.dirname(__file__))
    frames_dir = os.path.abspath(args.frames_dir or os.path.join(case_dir, "Video_pair"))
    os.makedirs(frames_dir, exist_ok=True)
    tasks = discover_tasks(case_dir, frames_dir, args.max_frames)
    config = RenderConfig(
        case_dir=case_dir,
        helper_dir=helper_dir,
        frames_dir=frames_dir,
        ldomain=args.ldomain,
        zbot=args.zbot,
        ztop=args.ztop,
        rmax=args.rmax,
    )
    workers = min(args.cpus, len(tasks))
    print(
        f"rendering {len(tasks)} frames with {workers} workers, "
        f"window z[{args.zbot:.2f},{args.ztop:.2f}] r+-{args.rmax:.2f}",
        flush=True,
    )
    with mp.get_context("spawn").Pool(
        processes=workers, initializer=initialise_worker, initargs=(config,)
    ) as pool:
        pool.map(render_frame, tasks, chunksize=1)
    missing = [
        task.target
        for task in tasks
        if not os.path.isfile(task.target) or os.path.getsize(task.target) == 0
    ]
    if missing:
        raise RuntimeError(f"Frame generation incomplete; first missing frame: {missing[0]}")
    if not args.skip_video_encode:
        encode_video(frames_dir, args.fps, args.out)
        print(f"wrote {args.out} ({len(tasks)} frames)", flush=True)
    else:
        print(f"wrote {len(tasks)} frames below {frames_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
