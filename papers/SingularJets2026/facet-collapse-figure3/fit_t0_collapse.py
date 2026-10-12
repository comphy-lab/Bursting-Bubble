#!/usr/bin/env python3
"""Virtual origins tau_0^- (before inception) and tau_0^+ (after) that best
collapse the Figure 3 profiles for a given exponent alpha.

Each profile (largest connected component of the interface, shifted by z_base)
is rescaled by |t - tau_0|^alpha, cropped to the rescaled panel and resampled
in arclength; the score is the pairwise symmetric Chamfer RMS between profiles,
measured in panel units of a 420 x 464 panel. tau_0 is scanned on a grid that
keeps every frame on its own side of the origin, and the minimum is reported.

Frame sets are fixed by the Letter's tau_0^L (the r_j^(1/0.629) fit over
0.01 < r_j < 0.06), so every candidate tau_0 is scored on the same snapshots:
after inception 1.5e-4 <= t - tau_0^L <= 1e-2; before inception
1e-3 <= tau_0^L - t <= 2e-2 (closer frames show the bubble-entrapment stall,
earlier cavities are not yet self-similar).

Run:  python3 fit_t0_collapse.py [--alpha 0.629 0.6667]
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

import compare_collapse_r1 as cc

HERE = Path(__file__).resolve().parent
FAC = HERE / "facets"
POST_TAU = (1.5e-4, 1.0e-2)
PRE_TAU = (1.0e-3, 2.0e-2)
BOX = {"post": ((-4.0, 4.0), (-0.5, 8.0)), "pre": ((-4.0, 4.0), (-0.3, 2.5))}
PANEL = (420.0, 464.0)
SAMPLES = 240


def frames(phase: str, tau_letter: float):
    """[(t, [array(r, z - z_base)])] for the fixed frame set of one phase."""
    _, post = cc.read_index(FAC / "index.txt")
    _, pre = cc.read_index(FAC / "index_pre.txt")
    if phase == "post":
        rows, prefix = post, "facetmain"
        times = [t for t in sorted(post) if POST_TAU[0] <= t - tau_letter <= POST_TAU[1]]
    else:
        rows, prefix = pre, "facetpremain"
        times = [t for t in sorted(pre) if PRE_TAU[0] <= tau_letter - t <= PRE_TAU[1]]
    out = []
    for t in times:
        paths = cc.chained_file(FAC / f"{prefix}_{t:.6f}.txt")
        keep = max(len(p) for p in paths)
        out.append((t, [np.column_stack([np.abs(p[:, 1]), p[:, 0] - rows[t][0]])
                        for p in paths if len(p) >= 0.3 * keep]))
    return out


def panel_points(comps, scale: float, phase: str) -> np.ndarray:
    """Rescale, map to panel units, crop to the panel and resample in arclength."""
    (x0, x1), (y0, y1) = BOX[phase]
    sx, sy = PANEL[0] / (x1 - x0), PANEL[1] / (y1 - y0)
    runs = []
    for c in comps:
        x, y = c[:, 0] / scale, c[:, 1] / scale
        inside = (x <= x1) & (y >= y0) & (y <= y1)
        edges = np.flatnonzero(np.diff(inside.astype(int))) + 1
        for sx_, sy_, k in zip(np.split(x, edges), np.split(y, edges), np.split(inside, edges)):
            if k[0] and len(sx_) >= 2:
                runs.append(np.column_stack([sx_ * sx, (sy_ - y0) * sy]))
    if not runs:
        return np.empty((0, 2))
    lengths = [np.sum(np.linalg.norm(np.diff(r, axis=0), axis=1)) for r in runs]
    out = []
    for r, length in zip(runs, lengths):
        s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(r, axis=0), axis=1))]
        keep = np.r_[True, np.diff(s) > 0]
        r, s = r[keep], s[keep]
        if len(r) < 2:
            continue
        target = np.linspace(0.0, s[-1], max(2, int(round(SAMPLES * length / sum(lengths)))))
        out.append(np.column_stack([np.interp(target, s, r[:, 0]), np.interp(target, s, r[:, 1])]))
    return np.concatenate(out) if out else np.empty((0, 2))


def score(fr, phase: str, alpha: float, tau0: float) -> float:
    sign = 1.0 if phase == "post" else -1.0
    pts = []
    for t, comps in fr:
        tau = sign * (t - tau0)
        if tau <= 0:
            return float("nan")
        p = panel_points(comps, tau ** alpha, phase)
        if len(p) < 12:
            return float("nan")
        pts.append(p)
    trees = [cKDTree(p) for p in pts]
    d2 = [np.mean(np.r_[trees[j].query(pts[i])[0], trees[i].query(pts[j])[0]] ** 2)
          for i in range(len(pts)) for j in range(i + 1, len(pts))]
    return float(np.sqrt(np.mean(d2)))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--alpha", type=float, nargs="+", default=[0.629])
    ap.add_argument("--step", type=float, default=2.5e-6)
    a = ap.parse_args()
    _, post = cc.read_index(FAC / "index.txt")
    tau_letter = round(cc.fit_t0(post, 0.629, (0.01, 0.06)), 6)
    fr = {ph: frames(ph, tau_letter) for ph in ("post", "pre")}
    grids = {"post": np.arange(0.49410, fr["post"][0][0] - 2.0e-5 + 1e-12, a.step),
             "pre": np.arange(fr["pre"][-1][0] + 4.0e-5, 0.49500 + 1e-12, a.step)}
    result = {"tau_letter": tau_letter, "frames": {ph: len(v) for ph, v in fr.items()}, "tau0": {}}
    for alpha in a.alpha:
        for ph in ("pre", "post"):
            s = np.array([score(fr[ph], ph, alpha, t0) for t0 in grids[ph]])
            i = int(np.nanargmin(s))
            result["tau0"][f"{ph} {alpha:.4f}"] = {"tau0": round(float(grids[ph][i]), 7), "rms": round(float(s[i]), 3)}
            print(f"alpha={alpha:.4f} {ph:4s} tau0={grids[ph][i]:.7f}  rms={s[i]:.2f}", flush=True)
    (HERE / "t0_collapse.json").write_text(json.dumps(result, indent=1) + "\n")


if __name__ == "__main__":
    main()
