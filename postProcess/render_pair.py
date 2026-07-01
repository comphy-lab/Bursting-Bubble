#!/usr/bin/env python3
"""Two-panel check-in for a drill run, ONE row, TWO columns:
   col 1: interface (green, mirrored) + red marker at the robust jet base
   col 2: same interface (green) + the adaptive mesh (cells)
Both panels share the same framing (a zoom on the jet/base region) and the same
interface colour. Suptitle carries t, ke, and maxlevel.

Marker uses getBase (tag.h main-body isolation -> outer free surface -> lowest
point), so it tracks the jet base and ignores satellite bubbles.

Usage: render_pair.py <case_dir> <Ldomain> <out.png>
"""
import os, sys, glob, subprocess as sp
from PIL import Image
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

CASE = os.path.abspath(sys.argv[1]); LDOMAIN = float(sys.argv[2]); OUT = sys.argv[3]
HERE = os.path.dirname(os.path.abspath(__file__))
GETFACET = os.path.join(HERE, "getFacet"); GETBASE = os.path.join(HERE, "getBase")
GETVIEW2D = os.path.join(HERE, "getView2D")
GREEN = (0.0, 0.5, 0.0)

def stderr_lines(cmd):
    p = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, cwd=CASE); _, e = p.communicate()
    return e.decode().split("\n")

def get_facets(rel):
    out = stderr_lines([GETFACET, rel]); segs = []; skip = False
    for i in range(len(out)):
        p = out[i].split()
        if not p: skip = False; continue
        if not skip and i + 1 < len(out):
            q = out[i+1].split()
            try: z1, r1 = float(p[0]), float(p[1]); z2, r2 = float(q[0]), float(q[1])
            except Exception: continue
            segs.append(((r1, z1), (r2, z2))); segs.append(((-r1, z1), (-r2, z2))); skip = True
    return segs

def get_base(rel):
    for ln in stderr_lines([GETBASE, rel]):
        s = ln.split()
        if len(s) >= 6:
            try:
                t, zb, rb, zt, rt, nout = float(s[0]), float(s[1]), float(s[2]), float(s[3]), float(s[4]), int(s[5])
                return t, zb, rb, zt, rt
            except Exception: pass
    return None

# --- latest solver-log state (t, ke, maxlevel) ---
with open(os.path.join(CASE, "log")) as fh:
    rows = [ln.split() for ln in fh if ln.strip() and ln.split()[0].lstrip("-").isdigit()]
last = rows[-1]
t_log, ke, maxlevel = float(last[2]), float(last[3]), int(last[4])

# --- matching snapshot ---
snaps = sorted(glob.glob(os.path.join(CASE, "intermediate", "snapshot-*")),
               key=lambda f: float(f.rsplit("snapshot-", 1)[-1]))
if not snaps:
    print("NO SNAPSHOTS YET"); sys.exit(0)
snap = min(snaps, key=lambda f: abs(float(f.rsplit("snapshot-", 1)[-1]) - t_log))
rel = os.path.join("intermediate", os.path.basename(snap))
snap_t = float(snap.rsplit("snapshot-", 1)[-1])

segs = get_facets(rel)
gb = get_base(rel)
t, zb, rb, zt, rt = gb if gb else (snap_t, -1000., -1000., -1000., -1000.)
print("t=%.4f ke=%.3f maxlevel=%d base=(%.4f,%.4f) tip_z=%.4f" % (t, ke, maxlevel, rb, zb, zt))

# --- shared framing window (zoom on jet + base + shoulders) ---
zbot = (zb if zb > -900 else min(z for _, z in [p for s in segs for p in s])) - 0.45
ztop = (zt if zt > -900 else 0.0) + 0.35
rmax = max(0.55, (abs(rb) if rb > -900 else 0.3) * 1.8 + 0.15)
zc = 0.5 * (zbot + ztop)

# --- col 2: adaptive-mesh render via getView2D, framed on the SAME window ---
# getView2D native orientation: z horizontal, r vertical; fov is the vertical
# (r) field of view. Empirical: visible half-r ~ D*L0*tan(fov/2) with D~1.78
# (calibrated vs known-fov renders). So to show r in [-rmax, rmax]:
#   fov = 2*atan(rmax / (D*L0)); and to show z-extent = winH horizontally the
#   render aspect W/H = winH/(2*rmax). After a 90deg CCW rotation z becomes
#   vertical, matching col1 (r horizontal, z vertical).
import math
winH = ztop - zbot
D = 2.55                                        # calibrated so col2 r-window ~ col1
fov = math.degrees(2.0 * math.atan(rmax / (D * LDOMAIN)))
W = 1100
H = max(200, int(round(W * (2.0 * rmax) / winH)))
raw = OUT + ".gv2d_raw.png"
sp.run([GETVIEW2D, rel, raw, "%.4f" % fov, "%.6f" % (-zc / LDOMAIN), "0", str(W), str(H), "0"],
       cwd=CASE, check=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
grid_img = Image.open(raw).rotate(90, expand=True)
os.remove(raw)

# --- assemble the 1x2 figure ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 8.5))

# col1: interface + marker
if segs:
    ax1.add_collection(LineCollection(segs, linewidths=2.2, colors=[GREEN]))
ax1.plot([0, 0], [zbot, ztop], "-.", color="grey", lw=0.8)
if rb > -900:
    ax1.plot(rb, zb, "o", color="red", ms=13, zorder=10)
    ax1.plot(-rb, zb, "o", color="red", ms=13, zorder=10)
ax1.set_xlim(-rmax, rmax); ax1.set_ylim(zbot, ztop)
ax1.set_aspect("equal"); ax1.axis("off")
ax1.set_title("interface + jet-base marker", fontsize=13)

# col2: interface + grid (raster from getView2D, already green)
ax2.imshow(grid_img)
ax2.axis("off")
ax2.set_title("interface + adaptive mesh", fontsize=13)

fig.suptitle(r"$t/\tau_\gamma = %.4f$    ke $= %.3f$    maxlevel $= %d$" % (t, ke, maxlevel),
             fontsize=15, y=0.98)
fig.savefig(OUT, bbox_inches="tight", dpi=140)
plt.close(fig)
print("wrote %s" % OUT)
