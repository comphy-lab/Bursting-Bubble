#!/usr/bin/env python3
"""CoMPhy-style field snapshot for bursting-bubble cases.

Split axisymmetric frame: left is log10 dissipation (Newtonian) or
log10 tr(A) (VE); right is velocity magnitude in Blues. Interface is
the CoMPhy cyan overlay used in Bursting-Bubble/Video.py.

Single-process PDF path uses LaTeX typography (publication-plots).
"""

from __future__ import annotations

import argparse
import glob
import os
import subprocess as sp
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.ticker import StrMethodFormatter


def configure_matplotlib(*, usetex: bool) -> None:
    """LaTeX for single stills; mathtext for video/multiprocessing."""
    matplotlib.rcParams["font.family"] = "serif"
    if usetex:
        matplotlib.rcParams["font.serif"] = ["Computer Modern Roman"]
        matplotlib.rcParams["text.usetex"] = True
        matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"
    else:
        matplotlib.rcParams["mathtext.fontset"] = "cm"
        matplotlib.rcParams["text.usetex"] = False

SCRIPT_DIR = Path(__file__).resolve().parent
HELPER_GETFACET = SCRIPT_DIR / "getFacet"
HELPER_GETDATA = SCRIPT_DIR / "getData"
INTERFACE_COLOR = "#00B2FF"


def run_helper(cmd: list[str], cwd: Path) -> list[str]:
    proc = sp.run(cmd, cwd=cwd, stdout=sp.PIPE, stderr=sp.PIPE, check=False)
    if proc.returncode != 0:
        raise RuntimeError(
            f"{' '.join(cmd)} failed ({proc.returncode}):\n"
            f"{proc.stderr.decode('utf-8', errors='replace')}"
        )
    return proc.stderr.decode("utf-8").split("\n")


def discover_snapshot(case_dir: Path, target_time: float) -> tuple[float, Path]:
    paths = sorted(Path(p) for p in glob.glob(str(case_dir / "intermediate" / "snapshot-*")))
    if not paths:
        raise FileNotFoundError(f"no snapshots in {case_dir / 'intermediate'}")
    timed = [(float(p.name.rsplit("snapshot-", 1)[-1]), p) for p in paths]
    timed.sort()
    return min(timed, key=lambda item: abs(item[0] - target_time))


def load_facets(snapshot: Path, case_dir: Path) -> list[tuple]:
    lines = run_helper([str(HELPER_GETFACET), str(snapshot.relative_to(case_dir))], cwd=case_dir)
    segs = []
    skip = False
    for i, line in enumerate(lines):
        parts = line.split()
        if not parts:
            skip = False
            continue
        if skip or i + 1 >= len(lines):
            continue
        nxt = lines[i + 1].split()
        if len(parts) < 2 or len(nxt) < 2:
            continue
        r1, z1 = float(parts[1]), float(parts[0])
        r2, z2 = float(nxt[1]), float(nxt[0])
        segs.append(((r1, z1), (r2, z2)))
        segs.append(((-r1, z1), (-r2, z2)))
        skip = True
    return segs


def load_fields(snapshot: Path, case_dir: Path, zmin: float, zmax: float, rmax: float, nr: int):
    lines = run_helper(
        [
            str(HELPER_GETDATA),
            str(snapshot.relative_to(case_dir)),
            str(zmin),
            "0",
            str(zmax),
            str(rmax),
            str(nr),
        ],
        cwd=case_dir,
    )
    Z, R, diss, vel, tra = [], [], [], [], []
    for line in lines:
        parts = line.split()
        if len(parts) < 4:
            continue
        Z.append(float(parts[0]))
        R.append(float(parts[1]))
        diss.append(float(parts[2]))
        vel.append(float(parts[3]))
        tra.append(float(parts[4]) if len(parts) > 4 else -1.0)
    R = np.asarray(R)
    Z = np.asarray(Z)
    nz = int(len(Z) / nr)
    shape = (nz, nr)
    return {
        "R": R.reshape(shape),
        "Z": Z.reshape(shape),
        "diss": np.asarray(diss).reshape(shape),
        "vel": np.asarray(vel).reshape(shape),
        "tra": np.asarray(tra).reshape(shape),
    }


def nice_vmax(values: np.ndarray, floor: float, *, percentile: float = 99.0) -> float:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return floor
    active = finite[finite > 0.05 * np.nanmax(finite)]
    sample = active if active.size > 32 else finite
    raw = float(np.percentile(sample, percentile))
    if raw <= floor:
        return floor
    mag = 10 ** np.floor(np.log10(max(raw, 1e-6)))
    return float(np.ceil(raw / mag) * mag)


def add_colorbar(fig, ax, mappable, *, align: str, label: str) -> None:
    left, bottom, width, height = ax.get_position().bounds
    if align == "left":
        cax = fig.add_axes([left - 0.070, bottom, 0.032, height])
        cbar = fig.colorbar(mappable, cax=cax)
        cbar.ax.yaxis.set_ticks_position("left")
        cbar.ax.yaxis.set_label_position("left")
    else:
        cax = fig.add_axes([left + width + 0.018, bottom, 0.032, height])
        cbar = fig.colorbar(mappable, cax=cax)
    cbar.set_label(label, fontsize=28, labelpad=10)
    cbar.ax.tick_params(labelsize=22, width=2.4, length=10, pad=6)
    cbar.outline.set_linewidth(2.4)
    cbar.ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.2f}"))


def resolve_mode(mode: str, fields: dict) -> str:
    if mode != "auto":
        return mode
    return "ve" if float(np.nanmax(fields["tra"])) > 0.2 else "newtonian"


def left_field_spec(mode: str, fields: dict, args: argparse.Namespace):
    if mode == "ve":
        return (
            fields["tra"],
            r"$\log_{10}\mathrm{tr}(\mathsf{A})$",
            args.tra_vmin,
            args.tra_vmax if args.tra_vmax is not None else 2.0,
            "PuOr_r",
        )
    return (
        fields["diss"],
        r"$\log_{10}\left(2\mu_r(\boldsymbol{\mathcal{D}:\mathcal{D}})\right)$",
        args.diss_vmin,
        args.diss_vmax if args.diss_vmax is not None else 2.0,
        "hot_r",
    )


def render_split_frame(args: argparse.Namespace, fields: dict, facets, time: float, out: Path, *, dpi: int) -> None:
    mode = resolve_mode(args.mode, fields)
    left, left_label, left_vmin, left_vmax, left_cmap = left_field_spec(mode, fields, args)
    vel_vmax = args.vel_vmax if args.vel_vmax is not None else 5.0
    fig, ax = plt.subplots(figsize=(14.0, 12.5))
    rmin, rmax = float(fields["R"].min()), float(fields["R"].max())
    zmin, zmax = float(fields["Z"].min()), float(fields["Z"].max())
    im_left = ax.imshow(
        left,
        cmap=left_cmap,
        interpolation="bilinear",
        origin="lower",
        extent=[-rmin, -rmax, zmin, zmax],
        vmin=left_vmin,
        vmax=left_vmax,
        rasterized=True,
    )
    im_right = ax.imshow(
        fields["vel"],
        cmap="Blues",
        interpolation="bilinear",
        origin="lower",
        extent=[rmin, rmax, zmin, zmax],
        vmin=args.vel_vmin,
        vmax=vel_vmax,
        rasterized=True,
    )
    if facets:
        ax.add_collection(
            LineCollection(
                facets, linewidths=2.2, colors=INTERFACE_COLOR, linestyle="solid", zorder=5
            )
        )
    ax.axvline(0.0, color="0.45", linestyle="-.", linewidth=1.6, zorder=4)
    ax.set_xlim(-args.rmax, args.rmax)
    ax.set_ylim(args.zmin, args.zmax)
    ax.set_aspect("equal")
    ax.set_xlabel(r"$r/R_0$", fontsize=32, labelpad=12)
    ax.set_ylabel(r"$z/R_0$", fontsize=32)
    ax.yaxis.set_label_coords(-0.38, 0.5)
    ax.tick_params(which="both", direction="out", width=3, labelsize=22, pad=8)
    ax.tick_params(which="major", length=12)
    ax.tick_params(which="minor", length=6)
    for spine in ax.spines.values():
        spine.set_linewidth(3)
    ax.minorticks_on()
    label = args.label if args.label else Path(args.case).name
    ax.set_title(rf"{label}: $t/\tau_0 = {time:.3f}$", fontsize=28, pad=14)
    fig.subplots_adjust(left=0.30, right=0.80, bottom=0.12, top=0.90)
    add_colorbar(fig, ax, im_left, align="left", label=left_label)
    add_colorbar(fig, ax, im_right, align="right", label=r"$\|\boldsymbol{u}\|$")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight", pad_inches=0.12, dpi=dpi)
    plt.close(fig)


def plot_case(args: argparse.Namespace) -> Path:
    configure_matplotlib(usetex=True)
    case_dir = Path(args.case).resolve()
    time, snapshot = discover_snapshot(case_dir, args.time)
    fields = load_fields(snapshot, case_dir, args.zmin, args.zmax, args.rmax, args.nr)
    facets = load_facets(snapshot, case_dir)
    mode = resolve_mode(args.mode, fields)
    out = Path(args.out) if args.out else case_dir / f"field_t{time:.4f}.pdf"
    render_split_frame(args, fields, facets, time, out, dpi=300)
    png = out.with_suffix(".png")
    render_split_frame(args, fields, facets, time, png, dpi=200)
    _, left_label, left_vmin, left_vmax, _ = left_field_spec(mode, fields, args)
    vel_vmax = args.vel_vmax if args.vel_vmax is not None else 5.0
    print(f"wrote {out} and {png} (mode={mode}, |u|_max={vel_vmax:g}, left=[{left_vmin:g},{left_vmax:g}])")
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", required=True, help="case directory containing intermediate/")
    parser.add_argument("--time", type=float, default=1.0)
    parser.add_argument("--mode", choices=("auto", "newtonian", "ve"), default="auto")
    parser.add_argument("--zmin", type=float, default=-2.5)
    parser.add_argument("--zmax", type=float, default=3.0)
    parser.add_argument("--rmax", type=float, default=2.0)
    parser.add_argument("--nr", type=int, default=512)
    parser.add_argument("--vel-vmin", type=float, default=0.0)
    parser.add_argument("--vel-vmax", type=float, default=5.0)
    parser.add_argument("--diss-vmin", type=float, default=-3.0)
    parser.add_argument("--diss-vmax", type=float, default=2.0)
    parser.add_argument("--tra-vmin", type=float, default=0.45)
    parser.add_argument("--tra-vmax", type=float, default=2.0)
    parser.add_argument("--label", type=str, default=None)
    parser.add_argument("--out", type=str, default=None)
    args = parser.parse_args()
    if not HELPER_GETDATA.is_file() or not HELPER_GETFACET.is_file():
        raise SystemExit("compile postProcess/getData and getFacet before plotting")
    plot_case(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
