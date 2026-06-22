"""
VideoFoot.py — jet-base / cavity-focus probe tracking + base flux for bursting bubbles.

Pipeline (case 1000):
  1. GATHER: run postProcess/getJetFoot on every snapshot ->
       t, (z_low,r_low), (z_maxk,r_maxk), and base fluxes (q,q_l) for each.
  2. LATCH (time-ordered hysteresis): rule 2 (max-|k|, cavity focus) until
     inception, rule 1 (jet base = lowest off-axis point) after, latched.
     Inception = first frame where the max-|k| point reaches the axis
     (0 <= r_maxk < R_AXIS_K) AND the lowest point is off-axis (r_low > AXIS_BAND).
     The probe (z_b,r_b) and its flux (q,q_l) follow the active regime.
  3. r_PROBE(t), z_PROBE(t), q(t), q_l(t) plot + foot.dat (the observables).
  4. RENDER: facets + marker (blue=rule2 focus, red=rule1 base), mirrored;
     ffmpeg -> MP4.

Usage:
  python VideoFoot.py --caseToProcess simulationCases/1000 [--CPUs 6] [--no-video]
"""
import os, tempfile, argparse, subprocess as sp
import multiprocessing as mp
from functools import partial
os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "mpl_foot"))
os.environ.setdefault("OMP_NUM_THREADS", "1")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GETFACET = os.path.join(SCRIPT_DIR, "getFacet")
GETFOOT  = os.path.join(SCRIPT_DIR, "getJetFoot")

# --- regime / inception tuning (geometry-tuned for case 1000) ---
AXIS_BAND = 0.04   # r_low above this => lowest point genuinely off-axis
R_AXIS_K  = 0.05   # r_maxk below this => focusing point reached the axis
SENT      = -999.  # sentinel guard


def run_helper(cmd, cwd):
    p = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, cwd=cwd)
    _, e = p.communicate()
    return e.decode("utf-8").split("\n")


def get_facets(rel, case_dir):
    out = run_helper([GETFACET, rel], case_dir); segs = []; skip = False
    if len(out) > 1:
        for i in range(len(out)):
            p = out[i].split()
            if not p: skip = False; continue
            if not skip and i + 1 < len(out):
                q = out[i + 1].split()
                    r1, z1 = float(p[1]), float(p[0]); r2, z2 = float(q[1]), float(q[0])
                except Exception:
                    continue
                segs.append(((r1, z1), (r2, z2))); segs.append(((-r1, z1), (-r2, z2))); skip = True
    return segs


def get_candidates(rel, case_dir):
    # t z_low r_low z_maxk r_maxk q_low ql_low q_maxk ql_maxk z_jet
    for line in run_helper([GETFOOT, rel], case_dir):
        s = line.split()
        if len(s) >= 10:
            try:
                return tuple(float(v) for v in s[:10])
            except Exception:
                pass
    return None


def gather_one(idx, case_dir, tsnap):
    t = tsnap * idx
    rel = os.path.join("intermediate", "snapshot-%.4f" % t)
    if not os.path.exists(os.path.join(case_dir, rel)):
        return None
    return (idx, get_candidates(rel, case_dir))


def latch_regimes(frames):
    """frames sorted: (idx,t,zlow,rlow,zmaxk,rmaxk,qlow,qllow,qmaxk,qlmaxk).
    Returns dict idx->(zb,rb,q,ql,regime) and incept_t."""
    jet = False; incept_t = None; out = {}
    for fr in frames:
        idx, t, zlow, rlow, zmaxk, rmaxk, qlow, qllow, qmaxk, qlmaxk, zjet = fr
        if not jet and (0.0 <= rmaxk < R_AXIS_K) and (rlow > AXIS_BAND):
            jet = True; incept_t = t
        if jet:
            out[idx] = (zlow, rlow, qlow, qllow, 1)
        else:
            out[idx] = (zmaxk, rmaxk, qmaxk, qlmaxk, 2)
    return out, incept_t


def render_one(item, case_dir, out_dir, bounds, tsnap):
    idx, zb, rb, regime = item
    t = tsnap * idx
    rel = os.path.join("intermediate", "snapshot-%.4f" % t)
    target = os.path.join(out_dir, "%08d.png" % int(round(t * 1000)))
    if not os.path.exists(os.path.join(case_dir, rel)):
        return
    zmin, zmax, rmax = bounds
    segs = get_facets(rel, case_dir)
    fig, ax = plt.subplots(figsize=(6.4, 9.0))
    ax.plot([0, 0], [zmin, zmax], "-.", color="grey", lw=1.2)
    ax.add_collection(LineCollection(segs, linewidths=2.2, colors="green"))
    if zb is not None and zb > SENT and rb > SENT:
        c = "red" if regime == 1 else "blue"
        ax.plot(rb, zb, "o", color=c, ms=15, zorder=10)
        ax.plot(-rb, zb, "o", color=c, ms=15, zorder=10)
    ax.set_xlim(-rmax, rmax); ax.set_ylim(zmin, zmax); ax.set_aspect("equal"); ax.axis("off")
    ax.set_title(r"$t/\tau_\gamma = %.3f$" % t, fontsize=16)
    plt.savefig(target, bbox_inches="tight", dpi=160); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--caseToProcess", default="simulationCases/1000")
    ap.add_argument("--CPUs", "--cpus", type=int, default=6, dest="cpus")
    ap.add_argument("--nGFS", type=int, default=151)
    ap.add_argument("--tsnap", type=float, default=0.01)
    ap.add_argument("--ZMIN", type=float, default=-2.2)
    ap.add_argument("--ZMAX", type=float, default=2.0)
    ap.add_argument("--RMAX", type=float, default=1.5)
    ap.add_argument("--no-video", action="store_true")
    a = ap.parse_args()

    case_dir = os.path.abspath(a.caseToProcess)
    out_dir = os.path.join(a.caseToProcess, "Video_foot")
    os.makedirs(out_dir, exist_ok=True)
    bounds = (a.ZMIN, a.ZMAX, a.RMAX)

    # 1. gather
    with mp.Pool(a.cpus) as pool:
        got = pool.map(partial(gather_one, case_dir=case_dir, tsnap=a.tsnap), range(a.nGFS))
    frames = []
    for g in got:
        if g is None or g[1] is None:
            continue
        idx, c = g
        frames.append((idx, c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], c[8], c[9]))
    frames.sort(key=lambda r: r[0])
    print("[gather] %d frames" % len(frames))

    # 2. latch
    chosen, incept_t = latch_regimes(frames)
    print("[latch] inception t ~ %s" % (("%.3f" % incept_t) if incept_t is not None else "none"))

    # 3. foot.dat + observables plot
    datp = os.path.join(a.caseToProcess, "foot.dat")
    with open(datp, "w") as fh:
        fh.write("# t z_b r_b q_jet q_l regime(1=base,2=focus) z_low r_low z_maxk r_maxk z_jet\n")
        for fr in frames:
            idx, ti = fr[0], fr[1]
            zb, rb, q, ql, reg = chosen[idx]
            fh.write("%.4f %.6e %.6e %.6e %.6e %d %.6e %.6e %.6e %.6e %.6e\n"
                     % (ti, zb, rb, q, ql, reg, fr[2], fr[3], fr[4], fr[5], fr[10]))
    print("[dat] wrote %s" % datp)

    ts = [fr[1] for fr in frames]
    reg = [chosen[fr[0]][4] for fr in frames]
    series = {
        r"$r_{\rm probe}$": [chosen[fr[0]][1] for fr in frames],
        r"$z_{\rm probe}$": [chosen[fr[0]][0] for fr in frames],
        r"$q_{\rm jet}=\int_0^{r_b}\!u_z\,r\,dr$": [chosen[fr[0]][2] for fr in frames],
        r"$q_\ell=\int_0^{r_b}\!u_z\,dr$": [chosen[fr[0]][3] for fr in frames],
    }
    def split(xs, ys, want):
        return ([x for x, g in zip(xs, reg) if g == want],
                [y for y, g in zip(ys, reg) if g == want])
    fig, axes = plt.subplots(len(series), 1, figsize=(8, 11), sharex=True)
    for ax, (lab, ys) in zip(axes, series.items()):
        x2, y2 = split(ts, ys, 2); x1, y1 = split(ts, ys, 1)
        ax.plot(x2, y2, ".", color="blue", ms=5, label="rule 2: max|k| (cavity focus)")
        ax.plot(x1, y1, ".", color="red", ms=5, label="rule 1: jet base")
        if incept_t is not None:
            ax.axvline(incept_t, color="k", ls="--", lw=1, label="inception ~%.2f" % incept_t)
        ax.set_ylabel(lab, fontsize=12); ax.grid(alpha=0.3)
    axes[0].legend(fontsize=9, loc="best")
    axes[0].set_title("Jet-base probe + base flux (case 1000)", fontsize=13)
    axes[-1].set_xlabel(r"$t/\tau_\gamma$", fontsize=13)
    probe_png = os.path.join(a.caseToProcess, "footprobe.png")
    plt.savefig(probe_png, bbox_inches="tight", dpi=150); plt.close(fig)
    print("[plot] wrote %s" % probe_png)

    # 3b. flux vs r_jet (= probe radius); blue (rule 2) and red (rule 1) in SEPARATE panels
    rb_all = [chosen[fr[0]][1] for fr in frames]
    qj_all = [chosen[fr[0]][2] for fr in frames]
    ql_all = [chosen[fr[0]][3] for fr in frames]

    def xy(want, ys):
        out = [(x, y) for x, y, g in zip(rb_all, ys, reg)
               if g == want and x > SENT and y > SENT]
        return [p[0] for p in out], [p[1] for p in out]

    fig2, ax2 = plt.subplots(2, 2, figsize=(11, 8.5))
    panels = [
        (ax2[0][0], 2, qj_all, "blue",  r"$q_{\rm jet}$", "rule 2: cavity focus"),
        (ax2[0][1], 1, qj_all, "red",   r"$q_{\rm jet}$", "rule 1: jet base"),
        (ax2[1][0], 2, ql_all, "blue",  r"$q_\ell$",      "rule 2: cavity focus"),
        (ax2[1][1], 1, ql_all, "red",   r"$q_\ell$",      "rule 1: jet base"),
    ]
    for ax, want, ys, col, ylab, ttl in panels:
        xs, yy = xy(want, ys)
        ax.plot(xs, yy, ".", color=col, ms=6)
        ax.set_xlabel(r"$r_{\rm jet}$", fontsize=12)
        ax.set_ylabel(ylab, fontsize=13)
        ax.set_title(ttl, fontsize=11)
        ax.grid(alpha=0.3)
    fig2.suptitle(r"Base flux vs $r_{\rm jet}$ (case 1000)", fontsize=13)
    scaling_png = os.path.join(a.caseToProcess, "footscaling.png")
    plt.savefig(scaling_png, bbox_inches="tight", dpi=150); plt.close(fig2)
    print("[plot] wrote %s" % scaling_png)

    if a.no_video:
        return

    # 4. render + encode
    items = [(fr[0], chosen[fr[0]][0], chosen[fr[0]][1], chosen[fr[0]][4]) for fr in frames]
    with mp.Pool(a.cpus) as pool:
        pool.map(partial(render_one, case_dir=case_dir, out_dir=out_dir, bounds=bounds, tsnap=a.tsnap), items)
    print("[render] done")
    case_no = os.path.basename(os.path.normpath(a.caseToProcess))
    mp4 = os.path.join(a.caseToProcess, case_no + "_foot.mp4")
    cmd = ["ffmpeg", "-y", "-framerate", "30", "-pattern_type", "glob",
           "-i", os.path.join(out_dir, "*.png"),
           "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", "-c:v", "libx264",
           "-r", "30", "-pix_fmt", "yuv420p", mp4]
    r = sp.run(cmd, capture_output=True, text=True)
    print("[ffmpeg] rc=%d -> %s" % (r.returncode, mp4))
    if r.returncode != 0:
        print(r.stderr[-800:])


if __name__ == "__main__":
    main()
