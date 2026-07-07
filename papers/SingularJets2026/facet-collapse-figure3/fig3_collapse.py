#!/usr/bin/env python3
"""Figure 3 — self-similar collapse of the Bo=0, Oh=0.03 Worthington-jet
interface profiles (Duchemin/Cattaneo-style, 4 panels).

Rows:  (top) PRE-inception cavity collapse ; (bottom) POST-inception jet.
Cols:  (left) raw (r/R0, z/R0) ; (right) rescaled by the self-similar time
       scale tau=|t-t0| with OUR exponent alpha (not 2/3):
         x = (r/R0) / (|t-t0|/t_ic)^alpha
         y = (z - z_base)/R0 / (|t-t0|/t_ic)^alpha
Lengths are in R0 and time in t_ic (the raw simulation units), so the ONLY
change from an inertio-capillary collapse is the exponent alpha.

Data (in facets/): the FULL raw interface (getFacet) per snapshot. The facet
segments are chained into continuous solid polylines and only the largest
connected component is drawn, so the entrapped bubble + shed droplets (their own
small components) drop out while the whole L15 interface — including the neck —
is kept. index files give, per snapshot time t: r_j (getBase jet radius) and
z_base (getBase main-body base, the shift reference — NOT a raw min(z), which
the bubble would pollute).

Run:  python3 fig3_collapse.py     (needs numpy + matplotlib)
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.collections import LineCollection
from matplotlib.ticker import FixedLocator, FuncFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
FAC = os.path.join(HERE, "facets")
ALPHA = 0.629          # cone exponent alpha(beta) for Oh=0.03, beta=38.4 deg
ZWIN = 0.30            # axial window above the base to keep (drop far free surface)

APS_SINGLE_COL = 3.375
FIG_HEIGHT = 3.45
LABEL_FONT = 8.5
TICK_FONT = 8
PANEL_FONT = 9
LINEWIDTH = 0.45

matplotlib.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": TICK_FONT,
    "axes.labelsize": LABEL_FONT,
    "xtick.labelsize": TICK_FONT,
    "ytick.labelsize": TICK_FONT,
    "text.usetex": True,
    "text.latex.preamble": r"\usepackage{amsmath}",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.unicode_minus": False,
})

def read_index(fn):
    """{t: (z_base, r_j)} from an index file (cols: t  r_j  z_base [ ...])."""
    out = {}
    for ln in open(os.path.join(FAC, fn)):
        if ln.startswith(("t0", "incept", "DONE")) or not ln.strip():
            continue
        p = ln.split()
        out[float(p[0])] = (float(p[2]), float(p[1]))
    return out

def load_segments(prefix, t):
    """Read getFacetMain output ('z r' pairs, one blank line between pairs) as
    VOF facet SEGMENTS [((z1,r1),(z2,r2)), ...] — the exact method the render
    scripts (render_pair.py) use: each facet is one short line drawn via a
    LineCollection, so NO curve reconstruction and NO spurious connections."""
    lines = open(os.path.join(FAC, f"{prefix}_{t:.6f}.txt")).read().splitlines()
    segs = []; i = 0
    while i < len(lines) - 1:
        p = lines[i].split()
        if not p:
            i += 1; continue
        q = lines[i + 1].split()
        try:
            z1, r1 = float(p[0]), float(p[1]); z2, r2 = float(q[0]), float(q[1])
        except (ValueError, IndexError):
            i += 1; continue
        segs.append(((z1, r1), (z2, r2))); i += 3   # skip the blank line after the pair
    return segs

def chain(segs, gap=1.2e-2):
    """Stitch the individual VOF facet segments into ordered, CONTINUOUS
    polylines so each interface draws as a SOLID line (not dashed facets).

    Greedy nearest-endpoint walk: from a path end, hop to the nearest unused
    facet endpoint and append that facet; break the path only when the nearest
    endpoint is farther than `gap` (a genuine hole in the surface). Facets are
    ~3e-4 long and share endpoints to ~1e-7 where the surface is well resolved,
    but the wide crater rim at early times is only sparsely sampled; `gap`=1.2e-2
    stitches those sparse stretches into solid lines while staying well below the
    inter-branch spacing. Axisymmetric output is single-valued in r (no r<0
    twin), so a larger gap cannot bridge across the axis; the entrapped bubble
    (already removed by getFacetMain) and detached droplets never join by a long
    spurious jump. Spatial-hash neighbour lookup keeps it O(n).
    Returns a list of [(z,r), ...] paths."""
    import math
    from collections import defaultdict
    n = len(segs); used = [False] * n
    H = max(gap, 1e-3)                                    # bucket size for neighbour search
    def bk(p):
        return (int(math.floor(p[0] / H)), int(math.floor(p[1] / H)))
    buckets = defaultdict(list)
    for i, (a, b) in enumerate(segs):
        buckets[bk(a)].append((i, 0)); buckets[bk(b)].append((i, 1))
    def nearest(p, excl):
        bx, by = bk(p); best = None; bestd = gap
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for j, e in buckets.get((bx + dx, by + dy), ()):
                    if used[j] or j == excl:
                        continue
                    q = segs[j][e]
                    d = math.hypot(q[0] - p[0], q[1] - p[1])
                    if d < bestd:
                        bestd = d; best = (j, e)
        return best
    lines = []
    for start in range(n):
        if used[start]:
            continue
        used[start] = True
        a, b = segs[start]; path = [a, b]
        for head, cur, excl in ((False, b, start), (True, a, start)):  # grow both ends
            c = cur; ex = excl
            while True:
                nb = nearest(c, ex)
                if nb is None:
                    break
                j, e = nb; used[j] = True; other = segs[j][1 - e]
                path.insert(0, other) if head else path.append(other)
                c = other; ex = j
        lines.append(path)
    return lines

def fit_t0(post):
    """t0 from r_j = C (t-t0)^alpha  =>  r_j^(1/alpha) linear in t (window in r_j)."""
    ts = np.array(sorted(post)); rj = np.array([post[t][1] for t in ts])
    m = (rj > 0.01) & (rj < 0.06)
    A, b = np.polyfit(ts[m], rj[m] ** (1.0 / ALPHA), 1)
    return -b / A

pre = read_index("index_pre.txt")
post = read_index("index.txt")
T0 = fit_t0(post)

fig = plt.figure(figsize=(APS_SINGLE_COL, FIG_HEIGHT))
panel_w = 0.365
panel_h = panel_w * APS_SINGLE_COL / FIG_HEIGHT
left = 0.145
gap_x = 0.100
bottom = 0.205
gap_y = 0.055
aRawPost = fig.add_axes([left, bottom, panel_w, panel_h])
aScPost = fig.add_axes([left + panel_w + gap_x, bottom, panel_w, panel_h])
aRawPre = fig.add_axes([left, bottom + panel_h + gap_y, panel_w, panel_h])
aScPre = fig.add_axes([left + panel_w + gap_x, bottom + panel_h + gap_y, panel_w, panel_h])
cax = fig.add_axes([0.265, 0.075, 0.545, 0.026])
CBAR_LO, CBAR_HI = 1e-4, 0.05            # fixed colorbar range for |t-t0|/t_ic
norm = colors.LogNorm(CBAR_LO, CBAR_HI); cmap = cm.viridis

def brk(x, y, cut):
    """Insert NaN breaks where a single step exceeds `cut` (in the panel's own
    units), so matplotlib lifts the pen instead of drawing a long straight
    chord. Real interface facets rescale to tiny steps; a chain mis-join near a
    sub-resolution nub, magnified by the 1/tau^alpha rescaling, makes one huge
    step — this removes only that spurious segment, not the smooth profile."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    d = np.hypot(np.diff(x), np.diff(y))
    xb = [x[0]]; yb = [y[0]]
    for k in range(len(d)):
        if d[k] > cut:
            xb.append(np.nan); yb.append(np.nan)
        xb.append(x[k + 1]); yb.append(y[k + 1])
    return np.array(xb), np.array(yb)

SC_CUT = 1.0        # max allowed step in rescaled coords (real facets rescale to <0.15)
TAU_FLOOR = 1.5e-4  # drop frames with |t-t0| below this from BOTH panels of a row, so the raw
                    # (uncompensated) and rescaled (compensated) panels show the SAME instants.
                    # These near-origin frames sit between physical nucleation and the fitted
                    # self-similar virtual origin t0: the jet is a grid-scale nub whose overhang
                    # the 1/tau^alpha blow-up turns into spurious teardrops/crowns. (The data
                    # still resolves closer to t0; we just don't plot those instants.)

def draw(prefix, times, zbmap, sign, axRaw, axSc, tau_floor=0.0):
    for t in sorted(times):
        tau = (T0 - t) if sign > 0 else (t - T0)
        if tau <= tau_floor:                       # same instant set in raw and rescaled
            continue
        zref = zbmap[t][0]                         # tag-based main-body base
        s = tau ** ALPHA; col = cmap(norm(tau))
        comps = chain(load_segments(prefix, t))    # connected components of the raw interface
        if not comps:
            continue
        keep = max(len(c) for c in comps)          # main body = largest component;
        comps = [c for c in comps if len(c) >= 0.3 * keep]   # drop entrapped bubble + droplets
        for path in comps:
            path = [(z, r) for (z, r) in path if (z - zref) <= ZWIN]  # keep jet/cavity window
            if len(path) < 2:
                continue
            zz = np.array([z for z, _ in path]); rr = np.array([r for _, r in path])
            axRaw.plot(rr, zz, color=col, lw=LINEWIDTH, solid_capstyle="round")   # solid, right half
            axRaw.plot(-rr, zz, color=col, lw=LINEWIDTH, solid_capstyle="round")  # mirror, left half
            xs, ys = brk(rr / s, (zz - zref) / s, SC_CUT)   # break rescaling-magnified mis-joins
            axSc.plot(xs, ys, color=col, lw=LINEWIDTH, solid_capstyle="round")
            axSc.plot(-xs, ys, color=col, lw=LINEWIDTH, solid_capstyle="round")

draw("facetpremain", pre, pre, +1, aRawPre, aScPre, tau_floor=TAU_FLOOR)
draw("facetmain",   post, post, -1, aRawPost, aScPost, tau_floor=TAU_FLOOR)

aRawPre.set_ylabel(r"$z/R$", labelpad=1)
aScPre.set_ylabel(r"$\eta$", labelpad=1)
aRawPost.set_xlabel(r"$r/R$", labelpad=1); aRawPost.set_ylabel(r"$z/R$", labelpad=1)
aScPost.set_xlabel(r"$\xi$", labelpad=1); aScPost.set_ylabel(r"$\eta$", labelpad=1)
for a in (aScPre, aScPost): a.set_xlim(-8, 8); a.set_ylim(-0.5, 8)
for a in (aRawPre, aRawPost): a.set_xlim(-0.35, 0.35)
aRawPre.set_ylim(-1.80, -1.33)     # pre-inception cavity z-range (LineCollection: set explicitly)
aRawPost.set_ylim(-1.52, -1.00)    # post-inception jet z-range

for label, ax in zip((r"(a)", r"(b)", r"(c)", r"(d)"),
                     (aRawPre, aScPre, aRawPost, aScPost)):
    ax.set_box_aspect(1)
    ax.text(0.03, 0.95, label, transform=ax.transAxes, ha="left", va="top",
            fontsize=PANEL_FONT, fontweight="bold")
    ax.tick_params(axis="both", which="major", direction="out", width=0.6,
                   length=2.6, pad=1.5)
    ax.tick_params(axis="both", which="minor", direction="out", width=0.4,
                   length=1.5)
    ax.minorticks_on()
    for spine in ax.spines.values():
        spine.set_linewidth(0.65)

for ax in (aRawPre, aRawPost):
    ax.xaxis.set_major_locator(FixedLocator([-0.3, 0.0, 0.3]))
for ax in (aScPre, aScPost):
    ax.xaxis.set_major_locator(FixedLocator([-6, 0, 6]))
    ax.yaxis.set_major_locator(FixedLocator([0, 4, 8]))
aRawPre.yaxis.set_major_locator(FixedLocator([-1.8, -1.6, -1.4]))
aRawPost.yaxis.set_major_locator(FixedLocator([-1.5, -1.25, -1.0]))

def signed_label(value, decimals=1):
    sign = r"\mbox{-}" if value < 0 else ""
    return rf"${sign}{abs(value):.{decimals}f}$"

def signed_int_label(value):
    sign = r"\mbox{-}" if value < 0 else ""
    return rf"${sign}{abs(int(value))}$"

for ax in (aRawPre, aRawPost):
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: signed_label(v, 1)))
for ax in (aScPre, aScPost):
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: signed_int_label(v)))
aRawPre.yaxis.set_major_formatter(FuncFormatter(lambda v, _: signed_label(v, 1)))
aRawPost.yaxis.set_major_formatter(FuncFormatter(lambda v, _: signed_label(v, 2)))
for ax in (aRawPre, aScPre):
    ax.tick_params(labelbottom=False)

sm = cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
cb = fig.colorbar(sm, cax=cax, orientation="horizontal", extend="both")
cb.set_ticks([1e-4, 1e-3, 1e-2, 0.05])
cb.ax.xaxis.set_major_formatter(FuncFormatter(
    lambda v, _: (r"$10^{\mbox{-}4}$" if abs(v-1e-4) < 1e-6 else
                  r"$10^{\mbox{-}3}$" if abs(v-1e-3) < 1e-5 else
                  r"$10^{\mbox{-}2}$" if abs(v-1e-2) < 1e-4 else
                  "0.05" if abs(v-0.05) < 1e-4 else "")))
cb.ax.tick_params(axis="x", which="major", direction="out", width=0.6, length=2.6,
                  labelsize=TICK_FONT, pad=1.5)
cb.ax.text(1.075, 0.5, r"$\Delta\tau$", transform=cb.ax.transAxes,
           ha="left", va="center", fontsize=LABEL_FONT)
cb.outline.set_linewidth(0.65)

out = os.path.join(HERE, "fig3_Oh0.03_collapse")
plt.savefig(out + ".png", dpi=600)
plt.savefig(out + ".pdf")
print("t0(fitted)=%.6f  alpha=%.3f  pre=%d post=%d" % (T0, ALPHA, len(pre), len(post)))
print("wrote fig3_Oh0.03_collapse.{png,pdf}")
