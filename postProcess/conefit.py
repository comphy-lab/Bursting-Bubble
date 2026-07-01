"""
conefit.py — measure the cavity cone semi-angle beta at jet inception, then
nu(beta), alpha(beta) and the predicted q_jet / q_l exponents; compare with the
measured log-log slopes. Pure Python (no scipy/numpy).

beta : fit a line r = m z + c to the conical cavity wall at inception;
       semi-angle from the AXIS is beta = atan(m)  (r = (z - z_apex) tan beta).
nu   : solve P_nu(cos beta) = 0  via P_nu(x) = 2F1(-nu, nu+1; 1; (1-x)/2).
alpha= 1/(2 - nu).  q_jet ~ r^((3a-1)/a),  q_l(paper)=r_jet v_jet ~ r^((2a-1)/a).
"""
import os, glob, math, subprocess as sp, argparse

def run(cmd, cwd):
    p = sp.Popen(cmd, stdout=sp.PIPE, stderr=sp.PIPE, cwd=cwd)
    _, e = p.communicate(); return e.decode().split("\n")

def facet_points(rel, case_dir):
    pts = []
    for ln in run([os.path.join(case_dir, "../../postProcess/getFacet"), rel], case_dir):
        s = ln.split()
        if len(s) == 2:
            try:
                z, r = float(s[0]), float(s[1])
                if r > 0: pts.append((z, r))
            except: pass
    return pts

def linfit(xs, ys):
    n = len(xs); sx = sum(xs); sy = sum(ys)
    sxx = sum(x*x for x in xs); sxy = sum(x*y for x, y in zip(xs, ys))
    m = (n*sxy - sx*sy) / (n*sxx - sx*sx)
    c = (sy - m*sx) / n
    # R^2
    yb = sy/n; ss_tot = sum((y-yb)**2 for y in ys)
    ss_res = sum((y-(m*x+c))**2 for x, y in zip(xs, ys))
    r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0.0
    return m, c, r2

def P_nu(nu, x):                      # Legendre P_nu(x) via Gauss 2F1
    z = (1.0 - x)/2.0
    a, b, c = -nu, nu+1.0, 1.0
    term = 1.0; s = 1.0
    for n in range(1, 2000):
        term *= (a+n-1)*(b+n-1)/((c+n-1)*n) * z
        s += term
        if abs(term) < 1e-15*max(1.0, abs(s)): break
    return s

def solve_nu(beta_deg):
    # cone surface at polar angle theta = 180-beta from the upward axis:
    # potential vanishes there => P_nu(cos(180-beta)) = P_nu(-cos beta) = 0.
    # (Taylor cone beta=49.3deg <-> nu=1/2; beta=90deg <-> nu=1.)
    x = -math.cos(math.radians(beta_deg))
    # scan upward from nu~0 for the FIRST sign change (smallest positive root)
    step = 0.005; nu = step; fprev = P_nu(nu, x); lo = None
    while nu < 2.0:
        nu2 = nu + step; f2 = P_nu(nu2, x)
        if fprev * f2 <= 0:
            lo, hi, flo = nu, nu2, fprev; break
        nu, fprev = nu2, f2
    if lo is None: return None
    for _ in range(200):                # bisect the bracket
        mid = 0.5*(lo+hi); fm = P_nu(mid, x)
        if flo*fm <= 0: hi = mid
        else: lo = mid; flo = fm
    return 0.5*(lo+hi)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default="simulationCases/1000")
    ap.add_argument("--dat", default="simulationCases/1000/foot.dat")
    ap.add_argument("--zlo", type=float, default=-1.30)   # cone-fit window (axial)
    ap.add_argument("--zhi", type=float, default=-0.25)
    ap.add_argument("--rlo", type=float, default=0.20)    # exclude near-axis jet column
    ap.add_argument("--rhi", type=float, default=1.10)    # exclude flat outer surface
    a = ap.parse_args()
    case_dir = os.path.abspath(a.case)

    # inception time + clean-window q fits from foot.dat
    rows = []
    with open(a.dat) as fh:
        for ln in fh:
            if ln.startswith('#') or not ln.strip(): continue
            rows.append([float(v) for v in ln.split()])
    reg1 = [r for r in rows if int(r[5]) == 1]
    if not reg1:
        raise SystemExit("conefit: no jet-base rows (regime==1) in %s — inception never "
                         "latched, nothing to fit. Check the run / latch thresholds." % a.dat)
    incept_t = reg1[0][0]
    # nearest actual snapshot: names may carry 4 or 6 decimals
    snaps = glob.glob(os.path.join(case_dir, "intermediate", "snapshot-*"))
    if not snaps:
        raise SystemExit("conefit: no snapshots in %s/intermediate" % case_dir)
    rel = os.path.relpath(
        min(snaps, key=lambda f: abs(float(f.rsplit("snapshot-", 1)[-1]) - incept_t)),
        case_dir)

    # --- cone fit on the cavity wall at inception ---
    pts = facet_points(rel, case_dir)
    wall = [(z, r) for (z, r) in pts if a.zlo <= z <= a.zhi and a.rlo <= r <= a.rhi]
    if len(wall) < 2:
        raise SystemExit("conefit: cone-fit window {z in [%.2f,%.2f], r in [%.2f,%.2f]} has %d "
                         "facet points (<2) at t=%.4f — widen --zlo/--zhi/--rlo/--rhi."
                         % (a.zlo, a.zhi, a.rlo, a.rhi, len(wall), incept_t))
    zs = [z for z, r in wall]; rs = [r for z, r in wall]
    m, c, r2 = linfit(zs, rs)
    beta = math.degrees(math.atan(m))
    z_apex = -c/m
    nu = solve_nu(beta)
    alpha = 1.0/(2.0 - nu)
    exp_qjet = (3*alpha - 1)/alpha
    exp_ql   = (2*alpha - 1)/alpha

    # --- measured log-log slopes over the clean jet window r_jet in [0.09,0.40] ---
    win = [r for r in reg1 if 0.09 <= r[2] <= 0.40 and r[3] > 0 and r[4] > 0]
    lr = [math.log(r[2]) for r in win]
    mqj, _, r2qj = linfit(lr, [math.log(r[3]) for r in win])
    mql, _, r2ql = linfit(lr, [math.log(r[4]) for r in win])

    print("=== cone fit at inception t=%.3f (window z[%.2f,%.2f] r[%.2f,%.2f], %d pts) ==="
          % (incept_t, a.zlo, a.zhi, a.rlo, a.rhi, len(wall)))
    print("  slope dr/dz = %.4f   beta = %.2f deg   z_apex = %.3f   R^2 = %.4f"
          % (m, beta, z_apex, r2))
    print("  nu(beta) = %.4f   alpha = 1/(2-nu) = %.4f" % (nu, alpha))
    print("  PREDICTED  q_jet ~ r^%.3f   q_l(paper) ~ r^%.3f" % (exp_qjet, exp_ql))
    print("=== measured log-log slopes (r_jet in [0.09,0.40], %d pts) ===" % len(win))
    print("  q_jet ~ r^%.3f (R^2=%.3f)   q_l ~ r^%.3f (R^2=%.3f)" % (mqj, r2qj, mql, r2ql))
    # sanity check of nu solver against paper table
    print("=== nu solver check (paper: 35->0.383, 40->0.423, 45->0.463, 49.3->0.5) ===")
    for b in (35, 40, 45, 49.3):
        print("  beta=%.1f -> nu=%.4f" % (b, solve_nu(b)))

if __name__ == "__main__":
    main()
