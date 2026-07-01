#!/usr/bin/env python3
"""Assemble a 2-column drill video: for every snapshot, one row / two columns --
   col 1: interface (green, mirrored) + red jet-base marker (getBase)
   col 2: same green interface + adaptive mesh (getView2D)
Same interface colour both panels; suptitle carries t, ke, maxlevel. A FIXED
framing window is used for temporal stability across frames.

Usage: render_pair_video.py <case_dir> <Ldomain> <out.mp4> [ZBOT ZTOP RMAX FPS]
Defaults frame the cavity+jet: ZBOT=-2.0 ZTOP=1.5 RMAX=1.1 FPS=18.
"""
import os, sys, glob, math, subprocess as sp
from PIL import Image
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

CASE = os.path.abspath(sys.argv[1]); LDOMAIN = float(sys.argv[2]); OUT = sys.argv[3]
ZBOT = float(sys.argv[4]) if len(sys.argv) > 4 else -2.0
ZTOP = float(sys.argv[5]) if len(sys.argv) > 5 else 1.5
RMAX = float(sys.argv[6]) if len(sys.argv) > 6 else 1.1
FPS  = int(sys.argv[7]) if len(sys.argv) > 7 else 18
HERE = os.path.dirname(os.path.abspath(__file__))
GETFACET = os.path.join(HERE, "getFacet"); GETBASE = os.path.join(HERE, "getBase")
GETVIEW2D = os.path.join(HERE, "getView2D")
GREEN = (0.0, 0.5, 0.0)
FRAMES = os.path.join(CASE, "Video_pair"); os.makedirs(FRAMES, exist_ok=True)

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
            try: return float(s[1]), float(s[2])   # zb, rb
            except Exception: pass
    return -1000., -1000.

# log lookup: snapshot t -> (ke, maxlevel) via nearest row
rows = []
with open(os.path.join(CASE, "log")) as fh:
    for ln in fh:
        c = ln.split()
        if c and c[0].lstrip("-").isdigit():
            rows.append((float(c[2]), float(c[3]), int(c[4])))
rows.sort()
def logat(t):
    if not rows: return (float("nan"), -1)
    best = min(rows, key=lambda r: abs(r[0] - t))
    return best[1], best[2]

snaps = sorted(glob.glob(os.path.join(CASE, "intermediate", "snapshot-*")),
               key=lambda f: float(f.rsplit("snapshot-", 1)[-1]))
zc = 0.5 * (ZBOT + ZTOP); winH = ZTOP - ZBOT
D = 3.0
fov = math.degrees(2.0 * math.atan(RMAX / (D * LDOMAIN)))
W = 1100; H = max(200, int(round(W * (2.0 * RMAX) / winH)))

print("rendering %d frames, window z[%.2f,%.2f] r+-%.2f" % (len(snaps), ZBOT, ZTOP, RMAX))
for k, snap in enumerate(snaps):
    rel = os.path.join("intermediate", os.path.basename(snap))
    t = float(snap.rsplit("snapshot-", 1)[-1])
    ke, maxlevel = logat(t)
    segs = get_facets(rel); zb, rb = get_base(rel)
    raw = os.path.join(FRAMES, "raw.png")
    sp.run([GETVIEW2D, rel, raw, "%.4f" % fov, "%.6f" % (-zc / LDOMAIN), "0", str(W), str(H), "0"],
           cwd=CASE, check=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
    grid_img = Image.open(raw).rotate(90, expand=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 8.5))
    if segs: ax1.add_collection(LineCollection(segs, linewidths=2.0, colors=[GREEN]))
    ax1.plot([0, 0], [ZBOT, ZTOP], "-.", color="grey", lw=0.7)
    if rb > -900:
        ax1.plot(rb, zb, "o", color="red", ms=11, zorder=10); ax1.plot(-rb, zb, "o", color="red", ms=11, zorder=10)
    ax1.set_xlim(-RMAX, RMAX); ax1.set_ylim(ZBOT, ZTOP); ax1.set_aspect("equal"); ax1.axis("off")
    ax1.set_title("interface + jet-base marker", fontsize=12)
    ax2.imshow(grid_img, extent=[-RMAX, RMAX, ZBOT, ZTOP], aspect="equal", origin="upper")
    ax2.set_xlim(-RMAX, RMAX); ax2.set_ylim(ZBOT, ZTOP); ax2.axis("off")
    ax2.set_title("interface + adaptive mesh", fontsize=12)
    fig.suptitle(r"$t/\tau_\gamma = %.4f$    ke $= %.3f$    maxlevel $= %d$" % (t, ke, maxlevel), fontsize=14, y=0.97)
    fig.savefig(os.path.join(FRAMES, "%05d.png" % k), bbox_inches="tight", dpi=120); plt.close(fig)
os.path.exists(os.path.join(FRAMES, "raw.png")) and os.remove(os.path.join(FRAMES, "raw.png"))

sp.run(["ffmpeg", "-y", "-framerate", str(FPS), "-i", os.path.join(FRAMES, "%05d.png"),
        "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-c:v", "libx264", "-pix_fmt", "yuv420p",
        "-crf", "20", OUT], check=True, stdout=sp.DEVNULL, stderr=sp.DEVNULL)
print("wrote %s (%d frames)" % (OUT, len(snaps)))
