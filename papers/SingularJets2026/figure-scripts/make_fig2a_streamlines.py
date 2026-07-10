#!/usr/bin/env python3
"""Build the four-frame streamline diagnostic for Fig. 2(a).

The default path reads the committed regular-grid fields and interface segments
and is fully offline.  Maintainers may explicitly refresh those extracts from
the checksum-verified raw snapshots archived on Dropbox.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess as sp
import tempfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

from capsule_utils import atomic_savefig, load_metadata, require_nonempty


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_DIR = Path(__file__).resolve().parent
METADATA = load_metadata()
FIG2A_METADATA = METADATA["figure_2a"]
DEFAULT_DATA_DIR = ROOT / "data-fig2a"
DEFAULT_SNAPSHOT_MANIFEST = DEFAULT_DATA_DIR / "raw-snapshots.json"
DEFAULT_SNAPSHOT_CACHE = Path(tempfile.gettempdir()) / "sbb_fig2a_5003"
DEFAULT_SNAPSHOTS = tuple(FIG2A_METADATA["snapshots"])

T0 = 0.49443
DEFAULT_VMAX = 50.0
INTERFACE_COLOR = "#D55E00"

APS = {
    "TitleFont": 7.5,
    "ColorbarFont": 7.0,
}


@dataclass(frozen=True)
class Field:
    z: np.ndarray
    r: np.ndarray
    f: np.ndarray
    uz: np.ndarray
    ur: np.ndarray
    speed: np.ndarray


def configure_matplotlib(use_tex: bool = True) -> None:
    matplotlib.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Computer Modern Roman", "DejaVu Serif"],
            "mathtext.fontset": "cm",
            "text.usetex": use_tex,
            "text.latex.preamble": r"\usepackage{amsmath}",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "axes.unicode_minus": False,
        }
    )


def run(cmd: list[str], *, cwd: Path | None = None, capture_stdout: bool = False) -> str:
    result = sp.run(
        cmd,
        cwd=cwd,
        text=True,
        stdout=sp.PIPE if capture_stdout else None,
        stderr=sp.PIPE,
    )
    if result.returncode != 0:
        raise SystemExit(
            "command failed:\n"
            f"  {' '.join(cmd)}\n"
            f"stderr:\n{result.stderr}"
        )
    return result.stdout if capture_stdout else result.stderr


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_download(path: Path, entry: dict) -> None:
    require_nonempty(path)
    if path.stat().st_size != int(entry["bytes"]):
        raise ValueError(f"Size mismatch for {path}")
    if sha256(path) != entry["sha256"]:
        raise ValueError(f"SHA-256 mismatch for {path}")


def ensure_snapshot_cache(
    case_dir: Path,
    manifest_path: Path,
    snapshots: tuple[str, ...],
) -> None:
    """Download the optional raw inputs only when a maintainer requests refresh."""
    with manifest_path.open(encoding="utf-8") as handle:
        manifest = json.load(handle)
    entries = {entry["snapshot"]: entry for entry in manifest["snapshots"]}
    case_dir.mkdir(parents=True, exist_ok=True)
    intermediate = case_dir / "intermediate"
    intermediate.mkdir(exist_ok=True)

    for snapshot in snapshots:
        if snapshot not in entries:
            raise KeyError(f"Snapshot {snapshot} is absent from {manifest_path}")
        entry = entries[snapshot]
        target = intermediate / f"snapshot-{snapshot}"
        if target.exists():
            validate_download(target, entry)
            continue
        descriptor, temporary_name = tempfile.mkstemp(
            prefix=f".{target.name}.", dir=target.parent
        )
        os.close(descriptor)
        temporary = Path(temporary_name)
        try:
            url = entry.get("url")
            if not url:
                raise RuntimeError(
                    f"No public Dropbox URL is recorded for snapshot {snapshot}; "
                    "the committed extracts remain fully offline, but maintainer refresh "
                    "requires an authenticated Dropbox link update"
                )
            with urllib.request.urlopen(url) as response, temporary.open("wb") as output:
                shutil.copyfileobj(response, output)
            validate_download(temporary, entry)
            os.replace(temporary, target)
        finally:
            temporary.unlink(missing_ok=True)


def compile_field_helper() -> Path:
    qcc = shutil.which("qcc")
    if qcc is None:
        raise SystemExit("qcc was not found on PATH; cannot compile the Basilisk field extractor.")

    source = SCRIPT_DIR / "extract_fig2a_fields.c"
    helper = Path(tempfile.gettempdir()) / "extract_fig2a_fields"
    if helper.exists() and helper.stat().st_mtime >= source.stat().st_mtime:
        return helper

    run(
        [qcc, "-O2", "-disable-dimensions", source.name, "-o", str(helper), "-lm"],
        cwd=SCRIPT_DIR,
    )
    return helper


def compile_facet_helper() -> Path:
    qcc = shutil.which("qcc")
    if qcc is None:
        raise SystemExit("qcc was not found on PATH; cannot compile the Basilisk facet extractor.")

    source = ROOT.parents[1] / "postProcess" / "getFacet.c"
    if not source.is_file():
        raise SystemExit(f"Missing postProcess source: {source}")
    helper = Path(tempfile.gettempdir()) / "sbb_getFacet"
    if helper.exists() and helper.stat().st_mtime >= source.stat().st_mtime:
        return helper

    run(
        [qcc, "-O2", "-disable-dimensions", source.name, "-o", str(helper), "-lm"],
        cwd=source.parent,
    )
    return helper


def parse_field(payload: str, nr: int) -> Field:
    data = np.loadtxt(payload.splitlines())
    if data.ndim != 2 or data.shape[1] != 6:
        raise ValueError("field extractor returned malformed data")
    nz = data.shape[0] // nr
    if nz * nr != data.shape[0]:
        raise ValueError("field extractor row count is not divisible by nr")

    z = data[:, 0].reshape(nz, nr)
    r = data[:, 1].reshape(nz, nr)
    f = data[:, 2].reshape(nz, nr)
    uz = data[:, 3].reshape(nz, nr)
    ur = data[:, 4].reshape(nz, nr)
    speed = data[:, 5].reshape(nz, nr)
    field = Field(z=z, r=r, f=f, uz=uz, ur=ur, speed=speed)
    validate_field(field)
    return field


def validate_field(field: Field) -> None:
    shape = field.z.shape
    if not shape or any(array.shape != shape for array in field.__dict__.values()):
        raise ValueError("Archived field arrays do not share one non-empty shape")
    if not np.all(np.isfinite(field.z)) or not np.all(np.isfinite(field.r)):
        raise ValueError("Archived grid coordinates contain non-finite values")
    if not np.all(np.diff(field.z[:, 0]) > 0.0) or not np.all(np.diff(field.r[0, :]) > 0.0):
        raise ValueError("Archived field coordinates are not strictly increasing")
    if not np.any(np.isfinite(field.speed)):
        raise ValueError("Archived velocity field contains no finite values")


def extract_field(
    helper: Path,
    case_dir: Path,
    snapshot: str,
    zmin: float,
    zmax: float,
    rmax: float,
    nr: int,
) -> Field:
    payload = run(
        [
            str(helper),
            f"intermediate/snapshot-{snapshot}",
            f"{zmin:.8g}",
            "0",
            f"{zmax:.8g}",
            f"{rmax:.8g}",
            str(nr),
        ],
        cwd=case_dir,
        capture_stdout=True,
    )
    return parse_field(payload, nr)


def facets(
    helper: Path,
    case_dir: Path,
    snapshot: str,
) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    raw = sp.run(
        [str(helper), f"intermediate/snapshot-{snapshot}"],
        cwd=case_dir,
        check=True,
        text=True,
        stdout=sp.PIPE,
        stderr=sp.PIPE,
    ).stderr.splitlines()

    segments: list[tuple[tuple[float, float], tuple[float, float]]] = []
    skip = False
    for idx, line in enumerate(raw):
        first = line.split()
        if not first:
            skip = False
            continue
        if skip or idx + 1 >= len(raw):
            continue
        second = raw[idx + 1].split()
        if len(first) != 2 or len(second) != 2:
            continue
        try:
            z1, r1 = float(first[0]), float(first[1])
            z2, r2 = float(second[0]), float(second[1])
        except ValueError:
            continue
        segments.append(((r1, z1), (r2, z2)))
        segments.append(((-r1, z1), (-r2, z2)))
        skip = True
    if not segments:
        raise ValueError(f"Facet extractor returned no interface segments for {snapshot}")
    return segments


def _atomic_npz(path: Path, **arrays: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{path.stem}.", suffix=".npz", dir=path.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        np.savez_compressed(temporary, **arrays)
        require_nonempty(temporary)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def field_path(data_dir: Path, snapshot: str) -> Path:
    return data_dir / "fields" / f"snapshot-{snapshot}.npz"


def facet_path(data_dir: Path, snapshot: str) -> Path:
    return data_dir / "interfaces" / f"snapshot-{snapshot}.npz"


def refresh_archived_inputs(args: argparse.Namespace) -> None:
    ensure_snapshot_cache(args.snapshot_cache, args.snapshot_manifest, tuple(args.snapshots))
    field_helper = compile_field_helper()
    facet_helper = compile_facet_helper()
    for snapshot in args.snapshots:
        field = extract_field(
            field_helper,
            args.snapshot_cache,
            snapshot,
            args.zmin,
            args.zmax,
            args.rmax,
            args.nr,
        )
        segments = facets(facet_helper, args.snapshot_cache, snapshot)
        _atomic_npz(field_path(args.data_dir, snapshot), **field.__dict__)
        _atomic_npz(
            facet_path(args.data_dir, snapshot),
            segments=np.asarray(segments, dtype=float),
        )


def load_archived_inputs(
    data_dir: Path,
    snapshots: tuple[str, ...],
) -> tuple[list[Field], list[list[tuple[tuple[float, float], tuple[float, float]]]]]:
    fields: list[Field] = []
    all_segments = []
    for snapshot in snapshots:
        with np.load(field_path(data_dir, snapshot)) as payload:
            field = Field(**{name: payload[name] for name in Field.__dataclass_fields__})
        validate_field(field)
        with np.load(facet_path(data_dir, snapshot)) as payload:
            segment_array = np.asarray(payload["segments"], dtype=float)
        if segment_array.ndim != 3 or segment_array.shape[1:] != (2, 2) or len(segment_array) == 0:
            raise ValueError(f"Malformed or empty archived interface segments for {snapshot}")
        if not np.all(np.isfinite(segment_array)):
            raise ValueError(f"Non-finite archived interface segments for {snapshot}")
        fields.append(field)
        all_segments.append(
            [(tuple(first), tuple(second)) for first, second in segment_array]
        )
    return fields, all_segments


def stream_arrays(field: Field, mirror: bool = False):
    r = field.r[0, :]
    z = field.z[:, 0]
    ur = np.ma.masked_invalid(field.ur)
    uz = np.ma.masked_invalid(field.uz)
    if mirror:
        return -r[::-1], z, -ur[:, ::-1], uz[:, ::-1]
    return r, z, ur, uz


def filtered_segments(
    segments: list[tuple[tuple[float, float], tuple[float, float]]],
    rmax: float,
    zmin: float,
    zmax: float,
) -> list[tuple[tuple[float, float], tuple[float, float]]]:
    kept = []
    for (r1, z1), (r2, z2) in segments:
        if max(abs(r1), abs(r2)) > 1.02 * rmax:
            continue
        if max(z1, z2) < zmin or min(z1, z2) > zmax:
            continue
        kept.append(((r1, z1), (r2, z2)))
    return kept


def draw_frame(
    ax: plt.Axes,
    field: Field,
    segments,
    snapshot: str,
    zmin: float,
    zmax: float,
    rmax: float,
    speed_norm: Normalize,
    cmap_speed,
    frame_label: str | None = None,
) -> None:
    speed = np.ma.masked_invalid(field.speed)
    mirrored_speed = np.ma.concatenate([speed[:, ::-1], speed], axis=1)

    ax.imshow(
        mirrored_speed,
        origin="lower",
        extent=[-rmax, rmax, zmin, zmax],
        cmap=cmap_speed,
        norm=speed_norm,
        interpolation="bilinear",
        rasterized=True,
        zorder=0,
    )

    for mirror in (True, False):
        xs, ys, uu, vv = stream_arrays(field, mirror=mirror)
        ax.streamplot(
            xs,
            ys,
            uu,
            vv,
            color="#6f6f6f",
            density=0.62,
            linewidth=0.25,
            arrowsize=0.32,
            minlength=0.08,
            maxlength=2.0,
            zorder=2,
        )

    interface = LineCollection(
        filtered_segments(segments, rmax, zmin, zmax),
        colors=INTERFACE_COLOR,
        linewidths=0.68,
        zorder=4,
    )
    ax.add_collection(interface)
    ax.axvline(0, color="0.72", lw=0.35, zorder=1)
    ax.set_xlim(-rmax, rmax)
    ax.set_ylim(zmin, zmax)
    ax.set_aspect("equal")
    ax.axis("off")

    if frame_label is not None:
        ax.text(
            0.04,
            0.94,
            frame_label,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=APS["TitleFont"],
            color="black",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.2},
            zorder=5,
        )


def build_figure(args: argparse.Namespace) -> None:
    configure_matplotlib(use_tex=not args.no_tex)
    if args.refresh_data:
        refresh_archived_inputs(args)
    fields, all_segments = load_archived_inputs(args.data_dir, tuple(args.snapshots))

    if args.vmax is None:
        all_speed = np.concatenate(
            [field.speed[np.isfinite(field.speed)].ravel() for field in fields]
        )
        vmax_speed = float(np.nanpercentile(all_speed, args.speed_percentile))
    else:
        vmax_speed = args.vmax
    speed_norm = Normalize(vmin=0.0, vmax=max(vmax_speed, 1e-6))

    cmap_speed = plt.get_cmap("Blues").copy()
    cmap_speed.set_bad((1, 1, 1, 0))

    fig, axes = plt.subplots(
        2,
        2,
        figsize=(3.42, 3.0),
        constrained_layout=False,
    )
    axes = np.ravel(axes)

    frame_labels = [r"(i)", r"(ii)", r"(iii)", r"(iv)"]
    for ax, field, segs, snap, label in zip(axes, fields, all_segments, args.snapshots, frame_labels):
        draw_frame(
            ax,
            field,
            segs,
            snap,
            args.zmin,
            args.zmax,
            args.rmax,
            speed_norm,
            cmap_speed,
            label,
        )

    for ax in axes[len(args.snapshots) :]:
        ax.axis("off")

    fig.subplots_adjust(left=0.01, right=0.99, top=0.995, bottom=0.16, wspace=0.015, hspace=0.025)

    cax1 = fig.add_axes([0.22, 0.065, 0.56, 0.035])
    sm_speed = ScalarMappable(norm=speed_norm, cmap=cmap_speed)
    cb1 = fig.colorbar(sm_speed, cax=cax1, orientation="horizontal")
    cb1.set_label(r"$|\mathbf{u}|$", fontsize=APS["ColorbarFont"], labelpad=1.5)
    cb1.ax.tick_params(labelsize=APS["ColorbarFont"], length=2.5, width=0.5, pad=1.0)
    cb1.outline.set_linewidth(0.5)

    atomic_savefig(fig, args.output, bbox_inches="tight", pad_inches=0.02, dpi=300)
    plt.close(fig)

    if args.frames_dir is not None:
        save_individual_frames(
            args,
            fields,
            all_segments,
            speed_norm,
            cmap_speed,
        )


def save_individual_frames(
    args: argparse.Namespace,
    fields: list[Field],
    all_segments,
    speed_norm: Normalize,
    cmap_speed,
) -> None:
    args.frames_dir.mkdir(parents=True, exist_ok=True)
    frame_labels = [r"(i)", r"(ii)", r"(iii)", r"(iv)"]
    for field, segs, snap, label in zip(fields, all_segments, args.snapshots, frame_labels):
        fig, ax = plt.subplots(figsize=(1.65, 1.65), constrained_layout=False)
        draw_frame(
            ax,
            field,
            segs,
            snap,
            args.zmin,
            args.zmax,
            args.rmax,
            speed_norm,
            cmap_speed,
            label,
        )
        fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
        atomic_savefig(
            fig,
            args.frames_dir / f"fig2a_streamlines_snapshot_{snap}.pdf",
            bbox_inches="tight",
            pad_inches=0.01,
            dpi=300,
        )
        plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--snapshot-cache", type=Path, default=DEFAULT_SNAPSHOT_CACHE)
    parser.add_argument("--snapshot-manifest", type=Path, default=DEFAULT_SNAPSHOT_MANIFEST)
    parser.add_argument("--snapshots", nargs="+", default=list(DEFAULT_SNAPSHOTS))
    parser.add_argument("--output", type=Path, default=ROOT / "fig2a_streamlines.pdf")
    parser.add_argument(
        "--frames-dir",
        type=Path,
        default=ROOT / "fig2a_streamlines_frames",
        help="Directory for per-snapshot PDF frame renders; set to '' to skip.",
    )
    parser.add_argument("--zmin", type=float, default=-1.72)
    parser.add_argument("--zmax", type=float, default=-0.82)
    parser.add_argument("--rmax", type=float, default=0.58)
    parser.add_argument("--nr", type=int, default=190)
    parser.add_argument("--vmax", type=float, default=DEFAULT_VMAX)
    parser.add_argument("--speed-percentile", type=float, default=99.2)
    parser.add_argument(
        "--refresh-data",
        action="store_true",
        help="Maintainer-only: download raw snapshots and atomically rebuild committed extracts.",
    )
    parser.add_argument("--no-tex", action="store_true")
    args = parser.parse_args()
    if args.frames_dir == Path(""):
        args.frames_dir = None
    return args


def main() -> None:
    build_figure(parse_args())


if __name__ == "__main__":
    main()
