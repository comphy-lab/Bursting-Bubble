#!/usr/bin/env python3
"""One-off render of the most recent snapshot, reusing Video.py's pipeline.

Run from the project root (where postProcess/ and simulationCases/ live):
    ~/miniconda3/envs/default/bin/python render_one.py [case] [out.png]
"""
import sys, os, glob

CASE = sys.argv[1] if len(sys.argv) > 1 else "simulationCases/1000"
OUT  = os.path.abspath(sys.argv[2] if len(sys.argv) > 2 else "latest_frame.png")

sys.path.insert(0, os.path.abspath("postProcess"))
import Video as V
from Video import RuntimeConfig, SnapshotInfo

case_dir = os.path.abspath(CASE)
snaps = sorted(glob.glob(os.path.join(case_dir, "intermediate", "snapshot-*")))
if not snaps:
    sys.exit("no snapshots found in %s/intermediate" % case_dir)
latest = max(snaps, key=lambda f: float(os.path.basename(f).replace("snapshot-", "")))
t = float(os.path.basename(latest).replace("snapshot-", ""))
# use the actual filename — snapshot names may carry 4 or 6 decimals
rel = os.path.relpath(latest, case_dir)

# Focused view on the collapsing cavity / jet region.
zmin, zmax, rmax = -2.5, 1.5, 1.5
grids_per_r = 256
nr = int(grids_per_r * rmax)

config = RuntimeConfig(
    cpus=1, n_snapshots=1, grids_per_r=grids_per_r, tsnap=0.01,
    zmin=zmin, zmax=zmax, rmax=rmax, case_dir=CASE,
    output_dir=os.path.join(CASE, "Video"), skip_video_encode=True,
    framerate=90, output_fps=30,
    d2_vmin=-2.0, d2_vmax=2.0, vel_vmin=0.0, vel_vmax=1.0,
)

facets = V.get_facets(rel, case_dir)
field = V.get_field(rel, case_dir, zmin, zmax, rmax, nr)
snap = SnapshotInfo(index=0, time=t, source=latest, target=OUT)
V.plot_snapshot(field, facets, config.bounds, snap, config, V.PLOT_STYLE)
print("WROTE %s  t=%.4f  facets=%d  nz=%d" % (OUT, t, len(facets), field.nz))
