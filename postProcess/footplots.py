"""
footplots.py — publication figures for the jet-base probe (PRL, 160 mm full width).

Reads foot.dat (written by VideoFoot.py):
  t z_b r_b q_jet q_l regime(1=base,2=focus) z_low r_low z_maxk r_maxk z_jet

Figure 1 (3 stacked panels, FULL process): r_probe(t), z_probe(t), z_jet(t).
Figure 2 (2 panels side by side, JET BASE ONLY, windowed inception -> max jet
          height): q_jet(r_jet), q_l(r_jet).

Run on machine-ts (has LaTeX + matplotlib):
  ~/miniconda3/envs/default/bin/python postProcess/footplots.py --dat simulationCases/1000/foot.dat --out simulationCases/1000
"""
import os, glob, argparse, math, subprocess as sp
import matplotlib
matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['font.serif'] = ['Computer Modern Roman']
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, FixedLocator, FuncFormatter

MM = 1.0 / 25.4
W = 160 * MM            # PRL full width = 160 mm
SENT = -999.
LAB, TICK, LEG = 10, 9, 8
BLUE, RED = '#1A64B3', '#C81E1E'
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GETFACET = os.path.join(SCRIPT_DIR, "getFacet")


def style(ax):
    ax.tick_params(which='major', labelsize=TICK, width=0.8, length=4, direction='out', pad=3)
    ax.tick_params(which='minor', width=0.6, length=2, direction='out')
    for s in ax.spines.values():
        s.set_linewidth(0.9)
    ax.minorticks_on()


def load(datp):
    R = []
    with open(datp) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            R.append([float(x) for x in line.split()])
    return R


# ---- cone fit + Legendre nu(beta), alpha(beta) (pure Python, no scipy) ----
def facet_points(rel, case_dir):
    p = sp.Popen([GETFACET, rel], stdout=sp.PIPE, stderr=sp.PIPE, cwd=case_dir)
    _, e = p.communicate()
    pts = []
    for ln in e.decode().split("\n"):
        s = ln.split()
        if len(s) == 2:
            try:
                z, r = float(s[0]), float(s[1])
                if r > 0:
                    pts.append((z, r))
            except Exception:
                pass
    return pts


def linfit(xs, ys):
    n = len(xs); sx = sum(xs); sy = sum(ys)
    sxx = sum(x*x for x in xs); sxy = sum(x*y for x, y in zip(xs, ys))
    m = (n*sxy - sx*sy) / (n*sxx - sx*sx); c = (sy - m*sx) / n
    return m, c


def P_nu(nu, x):
    z = (1.0 - x)/2.0; a, b, c = -nu, nu+1.0, 1.0; term = 1.0; s = 1.0
    for n in range(1, 2000):
        term *= (a+n-1)*(b+n-1)/((c+n-1)*n) * z; s += term
        if abs(term) < 1e-15*max(1.0, abs(s)):
            break
    return s


def solve_nu(beta_deg):           # P_nu(-cos beta) = 0, first positive root
    x = -math.cos(math.radians(beta_deg)); step = 0.005
    nu = step; fprev = P_nu(nu, x); lo = None
    while nu < 2.0:
        nu2 = nu + step; f2 = P_nu(nu2, x)
        if fprev*f2 <= 0:
            lo, hi, flo = nu, nu2, fprev; break
        nu, fprev = nu2, f2
    if lo is None:
        return None
    for _ in range(200):
        mid = 0.5*(lo+hi); fm = P_nu(mid, x)
        if flo*fm <= 0: hi = mid
        else: lo = mid; flo = fm
    return 0.5*(lo+hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dat', default='simulationCases/1000/foot.dat')
    ap.add_argument('--out', default='simulationCases/1000')
    ap.add_argument('--incept-t', type=float, default=None, dest='incept_t',
                    help='override the inception time (jet phase = t >= incept_t). The '
                         'foot.dat regime latch is case-1000-tuned and may never fire at '
                         'other Oh; take the true inception from the drill solver log '
                         '(first positive q_jet).')
    a = ap.parse_args()
    R = load(a.dat)
    t   = [r[0] for r in R]; zb = [r[1] for r in R]; rb = [r[2] for r in R]
    qj  = [r[3] for r in R]; ql = [r[4] for r in R]; reg = [int(r[5]) for r in R]
    zjet = [r[10] for r in R]
    N = len(t)

    # ---------- Figure 1: time series, full process ----------
    fig, ax = plt.subplots(3, 1, figsize=(W, 5.8), sharex=True)

    def regime_scatter(axi, ys):
        i2 = [i for i in range(N) if reg[i] == 2 and ys[i] > SENT]
        i1 = [i for i in range(N) if reg[i] == 1 and ys[i] > SENT]
        axi.plot([t[i] for i in i2], [ys[i] for i in i2], 'o', ms=2.4, color=BLUE,
                 label=r'rule 2: max$\,|\kappa|$ (focus)')
        axi.plot([t[i] for i in i1], [ys[i] for i in i1], 'o', ms=2.4, color=RED,
                 label=r'rule 1: jet base')

    regime_scatter(ax[0], rb); ax[0].set_ylabel(r'$r_{\mathrm{probe}}$', fontsize=LAB, labelpad=4)
    regime_scatter(ax[1], zb); ax[1].set_ylabel(r'$z_{\mathrm{probe}}$', fontsize=LAB, labelpad=4)
    ij = [i for i in range(N) if zjet[i] > SENT]
    ax[2].plot([t[i] for i in ij], [zjet[i] for i in ij], '-', color='k', lw=1.3)
    ax[2].set_ylabel(r'$z_{\mathrm{jet}}$', fontsize=LAB, labelpad=4)
    ax[2].set_xlabel(r'$t/\tau_\gamma$', fontsize=LAB, labelpad=4)
    for axi in ax:
        style(axi)
    ax[0].legend(fontsize=LEG, loc='upper center', frameon=False, ncol=2,
                 handletextpad=0.3, columnspacing=1.0)
    for axi, lab in zip(ax, ['(a)', '(b)', '(c)']):
        axi.text(-0.085, 1.04, lab, transform=axi.transAxes, fontsize=LAB)
    fig.set_facecolor('white'); fig.subplots_adjust(hspace=0.13)
    f1pdf = os.path.join(a.out, 'fig_foot_timeseries.pdf')
    f1png = os.path.join(a.out, 'fig_foot_timeseries.png')
    fig.savefig(f1pdf, bbox_inches='tight', dpi=300, pad_inches=0.05)
    fig.savefig(f1png, bbox_inches='tight', dpi=200, pad_inches=0.05)
    plt.close(fig)

    # ---------- window: inception (first rule 1) -> max jet height ----------
    if a.incept_t is not None:
        incept_t = a.incept_t          # user override (case-1000-tuned latch may never fire)
    else:
        jet_i = [i for i in range(N) if reg[i] == 1]
        incept_t = t[jet_i[0]] if jet_i else None
    valid = [(zjet[i], t[i]) for i in range(N) if zjet[i] > SENT]
    tmax = max(valid)[1] if valid else None
    if incept_t is None or tmax is None:
        raise SystemExit("footplots: no jet-base phase (regime==1) or no valid z_jet in %s — "
                         "cannot build the jet-base scaling figure (the time-series figure was "
                         "still written). If the latch never fired, rerun with --incept-t <t> "
                         "from the solver log's q_jet onset." % a.dat)

    def inwin(i):
        injet = (t[i] >= incept_t) if a.incept_t is not None else (reg[i] == 1)
        return (injet and rb[i] > SENT
                and t[i] >= incept_t and tmax is not None and t[i] <= tmax)
    W_i = [i for i in range(N) if inwin(i)]
    xr  = [rb[i] for i in W_i]
    yqj = [qj[i] for i in W_i]
    yql = [ql[i] for i in W_i]

    # ---------- Figure 2: q_jet, q_l vs r_jet (jet base only, windowed) — log-log ----------
    def pos(xs, ys):  # keep strictly positive pairs (log drops <=0 silently)
        px, py = [], []
        for x, y in zip(xs, ys):
            if x > 0. and y > 0.:
                px.append(x); py.append(y)
        return px, py
    xq, yq = pos(xr, yqj)
    xl, yl = pos(xr, yql)
    ndrop = (len(xr) - len(xq), len(xr) - len(xl))
    if any(ndrop):
        print("log-log dropped non-positive pts: q_jet=%d q_l=%d" % ndrop)

    # ---- cone fit at inception: beta, nu(beta), alpha; measured q_jet slope ----
    case_dir = os.path.abspath(a.out)
    # nearest actual snapshot: names may carry 4 or 6 decimals
    _snaps = glob.glob(os.path.join(case_dir, "intermediate", "snapshot-*"))
    if not _snaps:
        raise SystemExit("footplots: no snapshots in %s/intermediate" % case_dir)
    rel = os.path.relpath(
        min(_snaps, key=lambda f: abs(float(f.rsplit("snapshot-", 1)[-1]) - incept_t)),
        case_dir)
    pts = facet_points(rel, case_dir)
    wall = [(z, r) for (z, r) in pts if -1.30 <= z <= -0.25 and 0.20 <= r <= 1.10]
    if len(wall) < 2:
        raise SystemExit("footplots: cone-fit window has %d facet points (<2) at t=%.4f — "
                         "adjust the fit window or check the inception snapshot." % (len(wall), incept_t))
    cm, cc = linfit([z for z, r in wall], [r for z, r in wall])
    beta = math.degrees(math.atan(cm)); z_apex = -cc/cm
    nu = solve_nu(beta); alpha = 1.0/(2.0 - nu); exp_pred = (3*alpha - 1)/alpha
    fitpts = [(x, y) for x, y in zip(xq, yq) if 0.09 <= x <= 0.40]
    mqj, cqj = linfit([math.log(x) for x, y in fitpts], [math.log(y) for x, y in fitpts])
    print("beta=%.1f nu=%.3f alpha=%.3f  q_jet pred r^%.2f  measured r^%.2f"
          % (beta, nu, alpha, exp_pred, mqj))

    fig2, ax2 = plt.subplots(1, 2, figsize=(W, 3.2))

    # (a) q_jet vs r_jet: data + theory prediction (solid) + best fit (dashed)
    rmin = min(xq); rmax = max(xq)
    rr = [rmin * (rmax/rmin)**(k/60.0) for k in range(61)]
    # theory line q = K r^exp_pred: fit the prefactor K (slope fixed at exp_pred)
    # by least squares over the SAME window as the free best-fit (blue) line.
    logK = sum(math.log(y) - exp_pred*math.log(x) for x, y in fitpts) / len(fitpts)
    y_theory = [math.exp(logK) * r**exp_pred for r in rr]    # theory: K r^((3a-1)/a)
    y_fit    = [math.exp(cqj + mqj*math.log(r)) for r in rr] # free least-squares fit
    ax2[0].plot(rr, y_theory, '-',  color='k',  lw=1.3,
                label=r'theory $\propto r^{%.2f}$' % exp_pred)
    ax2[0].plot(rr, y_fit,    '--', color=BLUE, lw=1.3,
                label=r'fit $\propto r^{%.2f}$' % mqj)
    ax2[0].plot(xq, yq, 'o', ms=3.4, color=RED, mec='k', mew=0.3, label='data')
    ax2[0].set_xlabel(r'$r_{\mathrm{jet}}$', fontsize=LAB, labelpad=2)
    ax2[0].set_ylabel(r'$q_{\mathrm{jet}}=\int_0^{r_b}\! u_z\, r\,\mathrm{d}r$', fontsize=LAB, labelpad=2)
    ax2[0].legend(fontsize=7, loc='lower right', frameon=False,
                  handlelength=1.5, borderpad=0.2, labelspacing=0.25)
    ax2[0].text(0.04, 0.96,
                r'$\beta=%.1f^\circ,\ \nu=%.2f,\ \alpha=%.2f$' % (beta, nu, alpha),
                transform=ax2[0].transAxes, fontsize=7.5, va='top')

    # (b) q_l vs r_jet
    ax2[1].plot(xl, yl, 'o', ms=3.4, color=RED, mec='k', mew=0.3)
    ax2[1].set_xlabel(r'$r_{\mathrm{jet}}$', fontsize=LAB, labelpad=2)
    ax2[1].set_ylabel(r'$q_\ell=\int_0^{r_b}\! u_z\,\mathrm{d}r$', fontsize=LAB, labelpad=2)
    ax2[1].yaxis.set_label_coords(-0.16, 0.5)              # clear of panel (a)

    gfmt = FuncFormatter(lambda v, pos: '%g' % v)
    CAND = [0.02, 0.03, 0.05, 0.07, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0]
    def ticks_in(lo, hi):
        return [t for t in CAND if lo <= t <= hi]
    for axi in ax2:
        axi.set_xscale('log'); axi.set_yscale('log'); style(axi)
    # explicit, data-driven limits + ticks (avoids ScalarFormatter sub-decade rounding)
    axlo, axhi = min(xq)*0.85, max(xq)*1.15
    aylo, ayhi = min(yq)*0.90, max(yq)*2.6
    ax2[0].set_xlim(axlo, axhi); ax2[0].set_ylim(aylo, ayhi)
    ax2[0].xaxis.set_major_locator(FixedLocator(ticks_in(axlo, axhi)))
    ax2[0].yaxis.set_major_locator(FixedLocator(ticks_in(aylo, ayhi)))
    bylo, byhi = min(yl)*0.90, max(yl)*1.15
    ax2[1].set_xlim(axlo, axhi); ax2[1].set_ylim(bylo, byhi)
    ax2[1].xaxis.set_major_locator(FixedLocator(ticks_in(axlo, axhi)))
    ax2[1].yaxis.set_major_locator(FixedLocator(ticks_in(bylo, byhi)))
    for axi in ax2:
        for axis in (axi.xaxis, axi.yaxis):
            axis.set_major_formatter(gfmt)
            axis.set_minor_formatter(NullFormatter())

    # ---- inset on (a), upper-left below the beta/nu/alpha text: interface + fitted cone ----
    axin = ax2[0].inset_axes([0.07, 0.55, 0.38, 0.36])
    zz = [z for z, r in pts if -1.75 <= z <= 0.25 and r <= 1.25]
    rr2 = [r for z, r in pts if -1.75 <= z <= 0.25 and r <= 1.25]
    axin.plot(rr2, zz, '.', ms=0.7, color='green')
    axin.plot([-r for r in rr2], zz, '.', ms=0.7, color='green')
    zl = [z_apex, -0.20]; rl = [cm*z + cc for z in zl]
    axin.plot(rl, zl, '--', color=RED, lw=1.0)
    axin.plot([-x for x in rl], zl, '--', color=RED, lw=1.0)
    axin.plot([0, 0], [-1.75, 0.45], '-.', color='grey', lw=0.6)
    axin.set_aspect('equal'); axin.set_xlim(-1.25, 1.25); axin.set_ylim(-1.75, 0.55)
    axin.set_xticks([]); axin.set_yticks([])

    for axi, lab in zip(ax2, ['(a)', '(b)']):
        axi.text(-0.22, 1.03, lab, transform=axi.transAxes, fontsize=LAB)
    fig2.set_facecolor('white'); fig2.subplots_adjust(wspace=0.40)
    f2pdf = os.path.join(a.out, 'fig_foot_scaling.pdf')
    f2png = os.path.join(a.out, 'fig_foot_scaling.png')
    fig2.savefig(f2pdf, bbox_inches='tight', dpi=300, pad_inches=0.05)
    fig2.savefig(f2png, bbox_inches='tight', dpi=200, pad_inches=0.05)
    plt.close(fig2)

    print("inception=%.3f  t_maxheight=%.3f  window_points=%d"
          % (incept_t if incept_t else -1, tmax if tmax else -1, len(xr)))
    print("wrote", f1pdf, f2pdf)


if __name__ == '__main__':
    main()
