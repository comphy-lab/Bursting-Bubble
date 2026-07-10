#!/usr/bin/env python3
"""Re-extract Fig-3 interface facets for case 5003 with the FULL raw interface
(getFacet), not getFacetMain. getFacetMain's main-body tagging was dropping
~75% of the neck facets right after inception, giving a segmented neck. Here we
keep every facet (full L15 resolution); the entrapped bubble + droplets are
removed later in the plot by keeping only the largest connected component.

Writes facetpremain_<t>.txt (pre-inception) and facetmain_<t>.txt (post) plus
index_pre.txt / index.txt (t r_j z_base, from the run log) into facets_full/.
Facets are filtered to r<0.6 (well outside the r<0.35 view) to keep files small.
"""
import glob
import os
import subprocess
import tempfile

CASE = "/gpfs/work2/0/nctt0620/vatsal/2026-07-03-Singular-Bursting-Bubbles-Bo0-OhSweep-L14/simulationCases/5003"
GETFACET = os.path.join(CASE, "getFacet")
T0 = 0.494416          # inception time (fitted); pre < T0 < post
PRE_LO = T0 - 0.05     # cavity-collapse window to keep
RMAX = 0.6
OUT = os.path.join(CASE, "facets_full")
os.makedirs(OUT, exist_ok=True)

def tof(p):
    return float(os.path.basename(p).rsplit("-", 1)[1])

snaps = sorted(glob.glob(os.path.join(CASE, "intermediate", "snapshot-*")), key=tof)
pre = [p for p in snaps if PRE_LO <= tof(p) < T0]
post = [p for p in snaps if tof(p) > T0]

# ---- index from the run log (cols: 1=i 2=dt 3=t ... 8=r_base(r_j) 9=z_base) ----
logrows = []
for ln in open(os.path.join(CASE, "log")):
    p = ln.split()
    if len(p) < 9:
        continue
    try:
        logrows.append((float(p[2]), float(p[7]), float(p[8])))  # t, r_j, z_base
    except ValueError:
        continue
logrows.sort()
logts = [r[0] for r in logrows]
if not logrows:
    raise RuntimeError("No numeric rows found in the case log")


def atomic_write(path, payload):
    descriptor, temporary = tempfile.mkstemp(
        prefix=".%s." % os.path.basename(path), dir=os.path.dirname(path), text=True
    )
    try:
        with os.fdopen(descriptor, "w") as handle:
            handle.write(payload)
            handle.flush()
            os.fsync(handle.fileno())
        if os.path.getsize(temporary) == 0:
            raise RuntimeError("Refusing to install empty output: %s" % path)
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)

def nearest_log(t):
    import bisect
    i = bisect.bisect_left(logts, t)
    cand = [j for j in (i - 1, i, i + 1) if 0 <= j < len(logrows)]
    j = min(cand, key=lambda k: abs(logts[k] - t))
    return logrows[j]

def extract(path, prefix):
    t = tof(path)
    r = subprocess.run([GETFACET, path], capture_output=True, text=True, check=True)
    if not r.stderr.strip():
        raise RuntimeError("getFacet returned no facet payload for %s" % path)
    lines = r.stderr.splitlines()
    out = []; i = 0
    while i < len(lines) - 1:
        a = lines[i].split()
        if not a:
            i += 1; continue
        b = lines[i + 1].split()
        try:
            z1, r1 = float(a[0]), float(a[1]); z2, r2 = float(b[0]), float(b[1])
        except (ValueError, IndexError):
            i += 1; continue
        if r1 < RMAX and r2 < RMAX:
            out.append("%s %s\n%s %s\n\n" % (z1, r1, z2, r2))
        i += 3
    if not out:
        raise RuntimeError("No valid facets parsed for %s" % path)
    atomic_write(os.path.join(OUT, "%s_%.6f.txt" % (prefix, t)), "".join(out))
    return t

def run(frames, prefix, idxname):
    if not frames:
        raise RuntimeError("No snapshots selected for %s" % prefix)
    idx = ["t0 %.6f\n" % T0]
    for k, p in enumerate(frames):
        t = extract(p, prefix)
        _, rj, zb = nearest_log(t)
        idx.append("%.6f %.8e %.8e\n" % (t, rj, zb))
        if k % 10 == 0:
            print("%s %d/%d t=%.6f" % (prefix, k, len(frames), t), flush=True)
    atomic_write(os.path.join(OUT, idxname), "".join(idx))

run(pre, "facetpremain", "index_pre.txt")
run(post, "facetmain", "index.txt")
print("DONE pre=%d post=%d -> %s" % (len(pre), len(post), OUT), flush=True)
