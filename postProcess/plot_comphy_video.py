#!/usr/bin/env python3
"""CoMPhy field video from bursting-bubble snapshots.

Mathtext only (publication-plots): usetex deadlocks under multiprocessing.
Left panel is log10 dissipation (Newtonian) or log10 tr(A) in PuOr_r (VE);
right panel is |u| in Blues.
"""

from __future__ import annotations

import argparse
import glob
import multiprocessing as mp
import os
import shutil
import subprocess as sp
import sys
import tempfile
from functools import partial
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), f"mpl_{os.getpid()}"))

from plot_comphy_fields import (  # noqa: E402
    configure_matplotlib,
    load_facets,
    load_fields,
    render_split_frame,
)


def discover_all(case_dir: Path) -> list[tuple[float, Path]]:
    paths = [Path(p) for p in glob.glob(str(case_dir / "intermediate" / "snapshot-*"))]
    timed = [(float(p.name.rsplit("snapshot-", 1)[-1]), p) for p in paths]
    timed.sort()
    return timed


def _init_worker() -> None:
    os.environ["MPLCONFIGDIR"] = os.path.join(tempfile.gettempdir(), f"mpl_{os.getpid()}")
    configure_matplotlib(usetex=False)


def _render_one(item: tuple[int, float, str], args: argparse.Namespace) -> str:
    index, time, snapshot = item
    case_dir = Path(args.case).resolve()
    snap = Path(snapshot)
    fields = load_fields(snap, case_dir, args.zmin, args.zmax, args.rmax, args.nr)
    facets = load_facets(snap, case_dir)
    frame = Path(args.frames) / f"{index:05d}.png"
    render_split_frame(args, fields, facets, time, frame, dpi=args.dpi)
    return str(frame)


def ffmpeg_bin() -> str:
    nearby = str(Path(sys.executable).with_name("ffmpeg"))
    for candidate in (os.environ.get("FFMPEG"), shutil.which("ffmpeg"), nearby):
        if candidate and Path(candidate).is_file():
            return candidate
    raise RuntimeError("ffmpeg not found; set FFMPEG or install ffmpeg on PATH")


def encode_video(frames_dir: Path, out: Path, framerate: int, output_fps: int) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        ffmpeg_bin(),
        "-y",
        "-framerate",
        str(framerate),
        "-i",
        str(frames_dir / "%05d.png"),
        "-vf",
        "pad=ceil(iw/2)*2:ceil(ih/2)*2",
        "-c:v",
        "libx264",
        "-r",
        str(output_fps),
        "-pix_fmt",
        "yuv420p",
        str(out),
    ]
    proc = sp.run(cmd, stdout=sp.PIPE, stderr=sp.PIPE, text=True, check=False)
    if proc.returncode != 0:
        raise RuntimeError(f"ffmpeg failed ({proc.returncode}): {proc.stderr[-2000:]}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True)
    parser.add_argument("--mode", choices=("auto", "newtonian", "ve"), default="auto")
    parser.add_argument("--zmin", type=float, default=-2.5)
    parser.add_argument("--zmax", type=float, default=3.0)
    parser.add_argument("--rmax", type=float, default=2.0)
    parser.add_argument("--nr", type=int, default=256)
    parser.add_argument("--vel-vmin", type=float, default=0.0)
    parser.add_argument("--vel-vmax", type=float, default=5.0)
    parser.add_argument("--diss-vmin", type=float, default=-3.0)
    parser.add_argument("--diss-vmax", type=float, default=2.0)
    parser.add_argument("--tra-vmin", type=float, default=0.45)
    parser.add_argument("--tra-vmax", type=float, default=2.0)
    parser.add_argument("--label", type=str, default=None)
    parser.add_argument("--frames", type=str, required=True)
    parser.add_argument("--out", type=str, required=True)
    parser.add_argument("--cpus", type=int, default=4)
    parser.add_argument("--framerate", type=int, default=12)
    parser.add_argument("--output-fps", type=int, default=24)
    parser.add_argument("--dpi", type=int, default=120)
    parser.add_argument("--tmax", type=float, default=1.0)
    args = parser.parse_args()

    case_dir = Path(args.case).resolve()
    snaps = [(t, p) for t, p in discover_all(case_dir) if t <= args.tmax + 1e-9]
    if not snaps:
        raise SystemExit(f"no snapshots <= {args.tmax} in {case_dir / 'intermediate'}")
    frames_dir = Path(args.frames)
    frames_dir.mkdir(parents=True, exist_ok=True)
    for stale in frames_dir.glob("[0-9][0-9][0-9][0-9][0-9].png"):
        stale.unlink()
    items = [(i, t, str(p)) for i, (t, p) in enumerate(snaps)]
    print(f"rendering {len(items)} frames from {case_dir} -> {frames_dir}", flush=True)
    worker = partial(_render_one, args=args)
    if args.cpus <= 1:
        _init_worker()
        for item in items:
            print(worker(item), flush=True)
    else:
        with mp.Pool(processes=args.cpus, initializer=_init_worker) as pool:
            for frame in pool.imap(worker, items):
                print(frame, flush=True)
    out = Path(args.out)
    encode_video(frames_dir, out, args.framerate, args.output_fps)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
