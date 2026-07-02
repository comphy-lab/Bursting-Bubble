#!/usr/bin/env python3
"""Render the two check-in images for a running drill case:
   1) <prefix>_wide.png   -- interface (mirrored) + probe marker, whole box
   2) <prefix>_grid.png   -- adaptive mesh (cells) zoomed on the current probe

Both rotated 90deg CCW to the project convention (axial vertical/up,
radial horizontal). Marker POSITION (r_b, z_b) comes straight from the
solver's own `log` (authoritative -- it's what the solver's AMR trigger
actually used). Marker COLOR (regime) is re-derived with one getJetFoot call
on the matching snapshot, using the solver's own instantaneous criterion --
NOT a threshold on r_b, because post-inception the jet-base radius keeps
*shrinking* as the jet lengthens (the q_jet ~ r^0.95 slender-jet finding), so
r_b alone cannot distinguish "still focusing" from "jet formed, now slender".

Usage: monitor_snapshot.py <case_dir> <Ldomain> <out_prefix>
"""
import os, sys, glob, subprocess as sp
from PIL import Image
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

CASE = os.path.abspath(sys.argv[1])
LDOMAIN = float(sys.argv[2])
PREFIX = sys.argv[3]
HERE = os.path.dirname(os.path.abspath(__file__))
GETFACET = os.path.join(HERE, "getFacet")
GETFOOT = os.path.join(HERE, "getJetFoot")
GETVIEW2D = os.path.join(HERE, "getView2D")
AXIS_BAND = 0.04; R_AXIS_K = 0.05

def stderr_lines(cmd, cwd):
    p = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, cwd=cwd)
    _, e = p.communicate()
    return e.decode("utf-8").split("\n")

def get_facets(rel):
    out = stderr_lines([GETFACET, rel], CASE); segs = []; skip = False
    for i in range(len(out)):
        p = out[i].split()
        if not p: skip = False; continue
        if not skip and i + 1 < len(out):
            q = out[i + 1].split()
            try:
                r1, z1 = float(p[1]), float(p[0]); r2, z2 = float(q[1]), float(q[0])
            except Exception:
                continue
            segs.append(((r1, z1), (r2, z2))); segs.append(((-r1, z1), (-r2, z2))); skip = True
    return segs

def get_regime(rel):
    """Solver's own instantaneous inception test, applied post-hoc: jet-like
    if the max-curvature point sits near axis AND the lowest point is lifted
    off it. Returns 1 (jet base, red) or 2 (cavity focus, blue)."""
    for line in stderr_lines([GETFOOT, rel], CASE):
        s = line.split()
        if len(s) >= 10:
            try:
                _, zlow, rlow, zmaxk, rmaxk = (float(v) for v in s[:5])
            except Exception:
                continue
            if 0.0 <= rmaxk < R_AXIS_K and rlow > AXIS_BAND:
                return 1
            return 2
    return 2

# --- read the solver's own log for the latest state ---
with open(os.path.join(CASE, "log")) as fh:
    lines = [ln for ln in fh if ln.strip() and ln.split()[0].lstrip("-").isdigit()]
last = lines[-1].split()
i, dt, t, ke, maxlevel, rb, zb = (int(last[0]), float(last[1]), float(last[2]),
                                   float(last[3]), int(last[4]), float(last[5]), float(last[6]))
print("latest: i=%d t=%.5f ke=%.4f maxlevel=%d r_b=%.5f z_b=%.5f" % (i, t, ke, maxlevel, rb, zb))

# --- find the matching snapshot dump (nearest by t) ---
snaps = sorted(glob.glob(os.path.join(CASE, "intermediate", "snapshot-*")),
               key=lambda f: float(f.rsplit("snapshot-", 1)[-1]))
if not snaps:
    print("NO SNAPSHOTS YET"); sys.exit(0)
latest_snap = min(snaps, key=lambda f: abs(float(f.rsplit("snapshot-", 1)[-1]) - t))
snap_t = float(latest_snap.rsplit("snapshot-", 1)[-1])
rel = os.path.join("intermediate", os.path.basename(latest_snap))
print("using snapshot t=%.4f (%s)" % (snap_t, os.path.basename(latest_snap)))

regime = get_regime(rel)
print("regime=%d (%s)" % (regime, "jet base" if regime == 1 else "cavity focus"))

# --- image 1: wide interface + marker (matplotlib, r horizontal / z vertical) ---
segs = get_facets(rel)
fig, ax = plt.subplots(figsize=(6.0, 8.5))
ax.plot([0, 0], [-LDOMAIN, LDOMAIN], "-.", color="grey", lw=1.0)
if segs:
    ax.add_collection(LineCollection(segs, linewidths=2.0, colors="green"))
if rb > -900:
    col = "red" if regime == 1 else "blue"
    ax.plot(rb, zb, "o", color=col, ms=13, zorder=10)
    ax.plot(-rb, zb, "o", color=col, ms=13, zorder=10)
allz = [p[1] for s in segs for p in s] or [zb]
zmin, zmax = min(allz) - 0.3, max(allz) + 0.3
rmax = max(1.0, max(abs(rb) * 1.5, 0.8))
ax.set_xlim(-rmax, rmax); ax.set_ylim(zmin, zmax)
ax.set_aspect("equal"); ax.axis("off")
ax.set_title(r"$t/\tau_\gamma = %.4f$,  ke$=%.3f$,  maxlevel$=%d$" % (t, ke, maxlevel), fontsize=13)
wide_png = PREFIX + "_wide.png"
fig.savefig(wide_png, bbox_inches="tight", dpi=150)
plt.close(fig)
print("wrote %s" % wide_png)

# --- image 2: zoomed adaptive-mesh view via getView2D, centered on the probe ---
tx = -zb / LDOMAIN
ty = -rb / LDOMAIN
raw_png = PREFIX + "_grid_raw.png"
sp.run([GETVIEW2D, rel, raw_png, "5", str(tx), str(ty), "1000", "1000"], cwd=CASE, check=True,
       stdout=sp.DEVNULL, stderr=sp.DEVNULL)
grid_png = PREFIX + "_grid.png"
Image.open(os.path.join(CASE, raw_png)).rotate(90, expand=True).save(grid_png)
os.remove(os.path.join(CASE, raw_png))
print("wrote %s" % grid_png)
