#!/usr/bin/env python3
"""
plotJetMetricsTheory.py — jet-observable scaling figure for the self-similar
Worthington-jet study (CoMPhy Lab). Self-contained: all quantities are defined
here in terms of what the drill solver logs, independent of any particular
paper's notation (which is still in flux).

## Definitions (our own)

The drill solver's `log` (simulationCases/<case>/log), after two header lines,
carries whitespace-separated rows:
    i dt t ke maxlevel r_b z_b r_base z_base q_jet q_l
with, at the jet-base plane,
    q_jet = INT v_z r dr          (flow rate feeding the jet)
    q_l   = INT v_z dr            (flow rate per unit length ~ v_j r_j)
    r_base = r_j                  (robust getBase jet-base radius; NOT r_b, the
                                   AMR-tracking probe which can latch onto
                                   satellites late in the jet phase)

From these we plot two observables vs the jet radius r_j:

  (a)  q_j  = INT v_z r dr = q_jet
  (b)  We_j = v_j^2 r_j = q_l^2 / r_j        (with v_j := q_l / r_j)

## Theory

A conical cavity of semiangle beta collapses purely inertially; nu(beta) is the
smallest positive root of the Legendre condition P_nu(-cos beta) = 0 (this sign
convention reproduces the Taylor-cone anchor beta=49.3deg -> nu=0.5), and
alpha(beta) = 1/(2 - nu(beta)). The self-similar scalings are

    q_j (r_j) ~ r_j ^ ((3 alpha - 1) / alpha)
    We_j(r_j) ~ r_j ^ ((3 alpha - 2) / alpha)

**Our theory line** uses alpha from a cone fit at inception (per Oh).

**The inertio-capillary line** is *the same formula with alpha forced to 2/3*
(the classical inertio-capillary balance r_j ~ tau^{2/3}). That gives
    q_j  ~ r_j ^ 1.5      and      We_j ~ r_j ^ 0  (i.e. We_j = constant).
The purely-inertial cone theory has alpha < 2/3, so its q_j is shallower and
its We_j *decreases* with r_j rather than staying constant — the data
distinguishes the two.

**The PRF 2023 line** (Gordillo & Rodriguez-Rodriguez, Phys. Rev. Fluids 2023)
is *the same formula with alpha forced to 1/2* — the purely-inertial cone with
a CONSTANT far-field volume flux (their high-Laplace, La>2500 branch). That
gives, in the notation JR uses (Q_j the full flow rate ~ volume/time; q_l the
flow rate per unit outer perimeter, [q_l] = L^2/T):
    q_l  ~ r_j ^ 0   (i.e. q_l = const)   and   We_j ~ r_j ^ -1.
It is the constant that JR flagged: q_l is what is constant at short times, not
Q_j (= our q_j = INT v_z r dr).

## Fit windows (each scaling holds over a different r_j band)

The running logarithmic slope d ln q / d ln r_j of the pooled data is a smooth,
MONOTONIC crossover, not a set of clean plateaus: it starts near the cone
exponent at the deepest r_j (approaching inception) and decreases through the
PRF exponent at intermediate r_j. So each theory line is fit — and drawn — only
over the r_j band where its own scaling actually holds:

  * cone + inertio-capillary  -> near-inception window --fit-window (default
    [0.040, 0.054]), the r_j -> 0 self-similar asymptote where the data slope
    matches alpha(beta);
  * PRF 2023 (q_l = const / We_j ~ r_j^-1) -> finite window --prf-window
    (default [0.11, 0.19]), the q_l plateau / We_j ~ r_j^-1 crossing away from
    inception.

Theory-line SLOPES are fixed by the three alphas; PREFACTORS are least-squares
fit in log space over the matching window (never the full range). For the We_j
panels the inertio-capillary prediction is r_j-INDEPENDENT and order unity, so
it is drawn as a horizontal line at We_j = O(1) = 1 (NOT fit to the data) — the
physical claim the data refutes.

## Note on the two We_j definitions

We_j is shown two ways (rows of the 2x2 figure): We_j = q_l^2/r_j (from
q_l = INT v_z dr) and We_j = q_j^2/r_j^3 (from q_j = INT v_z r dr). They come
out with identical SHAPE, offset only by a constant ~4x (q_l^2/r_j ~=
4 q_j^2/r_j^3). That factor is exactly what a self-similar velocity profile
gives (INT v dr = v_j r_j, INT v r dr = v_j r_j^2/2, so the ratio is fixed):
the jet profile IS self-similar and the choice of definition is cosmetic — it
only shifts where the data crosses the We_j = O(1) line. q_l^2/r_j keeps We_j
further above unity (a more dramatic refutation of the inertio-capillary
prediction) and is the paper-consistent form (their q_j = INT v_z dr).

## Grouping

Data are grouped by (Oh, grid=MAXlevel): colour encodes Oh, marker encodes the
grid level. Theory lines are drawn per Oh in the matching colour (cone theory
solid, inertio-capillary alpha=2/3 dashed, PRF 2023 alpha=1/2 dotted).

## Usage

    python3 plotJetMetricsTheory.py \\
        --series 0.029 12 simulationCases/1007/log \\
        --series 0.029 13 simulationCases/1006/log \\
        --series 0.03  12 simulationCases/1012/log \\
        --facet  0.029 facets_oh029_inception.txt \\
        --facet  0.03  facets_oh030_inception.txt \\
        --out figures/jet_metrics

Each --series is (Oh, grid, log-path). A cone fit needs the inception facet
cloud per Oh: pass --facet OH FILE (a "z r" dump from postProcess/getFacet), or
let the script run getFacet itself against a case directory's intermediate/
snapshots (requires a case dir, not a bare log, and a local getFacet binary).

@author Vatsal Sanjay (vatsal.sanjay@comphy-lab.org) / CoMPhy Lab
"""
import os
import re
import math
import argparse
import subprocess as sp
from collections import defaultdict

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FuncFormatter
from matplotlib.lines import Line2D

matplotlib.rcParams["font.family"] = "serif"
matplotlib.rcParams["font.serif"] = ["Computer Modern Roman"]
matplotlib.rcParams["text.usetex"] = True
matplotlib.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
GETFACET = os.path.join(SCRIPT_DIR, "getFacet")

ALPHA_IC = 2.0 / 3.0                     # inertio-capillary balance
ALPHA_PRF = 0.5                          # Gordillo & Rodriguez-Rodriguez, PRF 2023
                                         # (La>2500 constant-far-field-flux branch)

# Okabe-Ito colourblind-safe, keyed by Oh in sorted order.
OH_COLOURS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9"]
# marker per grid level
GRID_MARKERS = {12: "o", 13: "s", 14: "^", 15: "D", 11: "v", 16: "P"}


# ============================== Legendre / cone-fit math ====================
# Identical to postProcess/conefit.py (validated against the Taylor-cone anchor
# beta=49.3 -> nu=0.5). Duplicated so this script has no import-path dependency.

def _linfit(xs, ys):
    n = len(xs)
    sx, sy = sum(xs), sum(ys)
    sxx = sum(x * x for x in xs)
    sxy = sum(x * y for x, y in zip(xs, ys))
    m = (n * sxy - sx * sy) / (n * sxx - sx * sx)
    c = (sy - m * sx) / n
    yb = sy / n
    ss_tot = sum((y - yb) ** 2 for y in ys)
    ss_res = sum((y - (m * x + c)) ** 2 for x, y in zip(xs, ys))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return m, c, r2


def _P_nu(nu, x):
    z = (1.0 - x) / 2.0
    a, b, c = -nu, nu + 1.0, 1.0
    term = 1.0
    s = 1.0
    for n in range(1, 2000):
        term *= (a + n - 1) * (b + n - 1) / ((c + n - 1) * n) * z
        s += term
        if abs(term) < 1e-15 * max(1.0, abs(s)):
            break
    return s


def _solve_nu(beta_deg):
    x = -math.cos(math.radians(beta_deg))
    step = 0.005
    nu = step
    fprev = _P_nu(nu, x)
    lo = None
    while nu < 2.0:
        nu2 = nu + step
        f2 = _P_nu(nu2, x)
        if fprev * f2 <= 0:
            lo, hi, flo = nu, nu2, fprev
            break
        nu, fprev = nu2, f2
    if lo is None:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        fm = _P_nu(mid, x)
        if flo * fm <= 0:
            hi = mid
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


# ============================== log / facet parsing ==========================

def parse_log(path):
    """Parse a drill-solver `log` into (rows, oh, maxlevel_hdr).

    rows: (N,6) float array [t, r_base, q_jet, q_l, maxlevel, z_base], time-
    ordered, sentinels dropped, and restart-overlap duplicates removed (a
    mid-run restart re-writes rows for an already-present t range; keep the
    last row seen per t).
    """
    oh = None
    ml_hdr = None
    keyed = {}
    with open(path) as fh:
        for ln in fh:
            if oh is None:
                m = re.search(r"Oh\s+([0-9.eE+-]+)", ln)
                if m:
                    oh = float(m.group(1))
                m2 = re.search(r"MAXlevel\s+(\d+)", ln)
                if m2:
                    ml_hdr = int(m2.group(1))
            p = ln.split()
            if not p:
                continue
            try:
                int(p[0])
            except ValueError:
                continue
            if len(p) < 11:
                continue
            try:
                t, ml = float(p[2]), int(p[4])
                r_base, z_base = float(p[7]), float(p[8])
                q_jet, q_l = float(p[9]), float(p[10])
            except ValueError:
                continue
            if r_base <= -900 or z_base <= -900:
                continue
            keyed[round(t, 8)] = (t, r_base, q_jet, q_l, ml, z_base)
    if not keyed:
        raise SystemExit("plotJetMetricsTheory: no numeric rows parsed from %s" % path)
    rows = np.array(sorted(keyed.values(), key=lambda r: r[0]))
    return rows, oh, ml_hdr


def reconnection_time(rows, pin_r=0.005, t_min=0.40):
    """Last t > t_min with the base still pinned on-axis — end of the singular
    focusing approach; the jet exists only beyond this point."""
    t, r_base = rows[:, 0], rows[:, 1]
    mask = (t > t_min) & (r_base < pin_r)
    if not mask.any():
        return None
    return float(t[mask].max())


def jet_phase(rows, incept_t, rmin=0.04):
    t, r_base, q_jet, q_l = rows[:, 0], rows[:, 1], rows[:, 2], rows[:, 3]
    m = (t > incept_t) & (r_base > rmin) & (q_jet > 0) & (q_l > 0)
    return r_base[m], q_jet[m], q_l[m]


def facet_points_from_file(path):
    pts = []
    with open(path) as fh:
        for ln in fh:
            s = ln.split()
            if len(s) == 2:
                try:
                    z, r = float(s[0]), float(s[1])
                    if r > 0:
                        pts.append((z, r))
                except ValueError:
                    pass
    return pts


def facet_points_from_case(case_dir, snap_t):
    import glob
    snaps = glob.glob(os.path.join(case_dir, "intermediate", "snapshot-*"))
    if not snaps or not os.path.exists(GETFACET):
        return None
    rel = os.path.relpath(
        min(snaps, key=lambda f: abs(float(f.rsplit("snapshot-", 1)[-1]) - snap_t)),
        case_dir)
    out, _ = sp.Popen([GETFACET, rel], cwd=case_dir,
                       stdout=sp.PIPE, stderr=sp.PIPE).communicate()
    pts = []
    for ln in out.decode(errors="ignore").split("\n"):
        s = ln.split()
        if len(s) == 2:
            try:
                z, r = float(s[0]), float(s[1])
                if r > 0:
                    pts.append((z, r))
            except ValueError:
                pass
    return pts


def cone_fit(pts, zlo, zhi, rlo, rhi):
    wall = [(z, r) for (z, r) in pts if zlo <= z <= zhi and rlo <= r <= rhi]
    if len(wall) < 2:
        return None
    m, c, r2 = _linfit([z for z, r in wall], [r for z, r in wall])
    beta = math.degrees(math.atan(m))
    nu = _solve_nu(beta)
    if nu is None:
        return None
    return dict(beta=beta, nu=nu, alpha=1.0 / (2.0 - nu), r2=r2, n=len(wall))


# ============================== fitting / plotting helpers ===================

def powerfit_prefactor(r, q, slope, rlo, rhi):
    m = (r >= rlo) & (r <= rhi) & (q > 0)
    if m.sum() < 2:
        return None, 0
    return math.exp(np.mean(np.log(q[m]) - slope * np.log(r[m]))), int(m.sum())


def thin_logspace(r, q, n=250):
    idx = np.argsort(r)
    r, q = r[idx], q[idx]
    if len(r) <= n:
        return r, q
    targets = np.logspace(np.log10(r.min()), np.log10(r.max()), n)
    keep = np.unique(np.searchsorted(r, targets).clip(0, len(r) - 1))
    return r[keep], q[keep]


def nice_ticks(vmin, vmax):
    cand = [1e-3, 2e-3, 3e-3, 5e-3, 1e-2, 2e-2, 3e-2, 5e-2, 0.1, 0.2, 0.3, 0.5,
            1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 300, 500, 1000]
    ticks = [c for c in cand if vmin * 0.9 <= c <= vmax * 1.1]
    if len(ticks) < 3:
        ticks = list(np.geomspace(vmin, vmax, 4))
    return ticks


def style_log_axis(ax, xlim, ylim):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    fmt = FuncFormatter(lambda v, _: f"{v:g}")
    ax.xaxis.set_major_locator(FixedLocator(nice_ticks(*xlim)))
    ax.yaxis.set_major_locator(FixedLocator(nice_ticks(*ylim)))
    ax.xaxis.set_minor_locator(FixedLocator([]))
    ax.yaxis.set_minor_locator(FixedLocator([]))
    ax.xaxis.set_major_formatter(fmt)
    ax.yaxis.set_major_formatter(fmt)
    ax.tick_params(which="major", direction="out", width=1.7, length=8, labelsize=15, pad=5)
    ax.tick_params(which="minor", length=0)
    for s in ax.spines.values():
        s.set_linewidth(1.7)


# ==================================== main ====================================

def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--series", nargs=3, action="append", required=True,
                     metavar=("OH", "GRID", "LOG"),
                     help="One data series: Oh value, grid (MAXlevel), and log path "
                          "(or a case directory containing `log`). Repeatable.")
    ap.add_argument("--facet", nargs=2, action="append", default=[],
                     metavar=("OH", "FILE"),
                     help="Inception facet 'z r' dump for the cone fit of a given Oh. "
                          "Repeatable (one per Oh). If omitted for an Oh, getFacet is "
                          "run against that Oh's first series case dir.")
    ap.add_argument("--incept-t", nargs=2, action="append", default=[],
                     metavar=("OH", "T"), dest="incept_t",
                     help="Override reconnection/inception time for an Oh (else "
                          "auto-detected per series: last t>0.40 with r_base<0.005).")
    ap.add_argument("--fit-window", nargs=2, type=float, default=[0.040, 0.054],
                     dest="fit_window",
                     help="r_j window for the cone / inertio-capillary prefactor "
                          "fit — the near-inception (r_j->0) self-similar asymptote "
                          "where the data slope matches alpha(beta).")
    ap.add_argument("--prf-window", nargs=2, type=float, default=[0.11, 0.19],
                     dest="prf_window",
                     help="r_j window for the PRF 2023 (alpha=1/2) prefactor fit — "
                          "the FINITE intermediate band where q_l=const / We_j~r_j^-1 "
                          "hold (the q_l plateau), away from inception.")
    ap.add_argument("--cone-window", nargs=4, type=float,
                     default=[-1.60, -0.40, 0.10, 1.00], dest="cone_window",
                     metavar=("ZLO", "ZHI", "RLO", "RHI"))
    ap.add_argument("--out", required=True, help="Output stem (.png and .pdf).")
    args = ap.parse_args()

    facet_by_oh = {float(o): f for o, f in args.facet}
    incept_by_oh = {float(o): float(t) for o, t in args.incept_t}
    fit_rlo, fit_rhi = args.fit_window
    prf_rlo, prf_rhi = args.prf_window
    zlo, zhi, rlo_c, rhi_c = args.cone_window

    # ---- load series ----
    series = []
    for oh_s, grid_s, log_s in args.series:
        oh = float(oh_s)
        grid = int(grid_s)
        case_dir = log_s if os.path.isdir(log_s) else os.path.dirname(os.path.abspath(log_s))
        log_path = log_s if os.path.isfile(log_s) else os.path.join(log_s, "log")
        if not os.path.exists(log_path):
            raise SystemExit("plotJetMetricsTheory: no log for --series entry %r" % log_s)
        rows, oh_hdr, _ = parse_log(log_path)
        incept_t = incept_by_oh.get(oh, reconnection_time(rows))
        if incept_t is None:
            raise SystemExit("plotJetMetricsTheory: no reconnection time for %s — pass "
                              "--incept-t %g <t>." % (log_s, oh))
        r_j, q_jet, q_l = jet_phase(rows, incept_t)
        if r_j.size < 5:
            raise SystemExit("plotJetMetricsTheory: only %d jet rows for %s (incept_t=%.4f)."
                              % (r_j.size, log_s, incept_t))
        # Two estimators of the local Weber number We_j = v_j^2 r_j, using the
        # two velocity-profile moments the solver logs:
        #   We_ql = q_l^2 / r_j       with v_j := q_l / r_j    (0th moment INT v dr)
        #   We_qj = q_jet^2 / r_j^3   with v_j := q_jet / r_j^2 (1st moment INT v r dr)
        # They estimate the SAME physical We_j and differ only by a
        # velocity-profile-shape constant, so comparing them tests whether the
        # profile is self-similar.
        We_ql = q_l ** 2 / r_j
        We_qj = q_jet ** 2 / r_j ** 3
        series.append(dict(oh=oh, grid=grid, case_dir=case_dir, incept_t=incept_t,
                           r_j=r_j, q_j=q_jet, q_l=q_l, We_ql=We_ql, We_qj=We_qj))

    ohs = sorted({s["oh"] for s in series})
    oh_colour = {oh: OH_COLOURS[i % len(OH_COLOURS)] for i, oh in enumerate(ohs)}

    # ---- cone fit + alpha per Oh ----
    alpha_by_oh = {}
    for oh in ohs:
        oh_series = [s for s in series if s["oh"] == oh]
        pts = None
        if oh in facet_by_oh:
            pts = facet_points_from_file(facet_by_oh[oh])
        else:
            pts = facet_points_from_case(oh_series[0]["case_dir"], oh_series[0]["incept_t"])
        cone = cone_fit(pts, zlo, zhi, rlo_c, rhi_c) if pts else None
        alpha_by_oh[oh] = cone
        if cone:
            print("Oh=%.4g: cone fit beta=%.2fdeg nu=%.4f alpha=%.4f (R^2=%.4f, n=%d)"
                  % (oh, cone["beta"], cone["nu"], cone["alpha"], cone["r2"], cone["n"]))
        else:
            print("Oh=%.4g: CONE FIT FAILED (no facet data) — theory line skipped" % oh)

    # ---- exponents. Both We_j estimators share the same We_j exponent. ------
    def exps(alpha):
        e_we = (3 * alpha - 2) / alpha
        return {"q_j": (3 * alpha - 1) / alpha,
                "q_l": (2 * alpha - 1) / alpha,
                "We_qj": e_we, "We_ql": e_we}
    ic_exp = exps(ALPHA_IC)   # q_j->1.5, q_l->0.5, We->0 (constant)
    print("inertio-capillary (alpha=2/3): q_j~r^%.3f  q_l~r^%.3f  We_j~r^%.3f (constant)"
          % (ic_exp["q_j"], ic_exp["q_l"], ic_exp["We_ql"]))
    prf_exp = exps(ALPHA_PRF)  # q_j->1, q_l->0 (const), We->-1
    print("PRF 2023 (alpha=1/2):          q_j~r^%.3f  q_l~r^%.3f (const)  We_j~r^%.3f"
          % (prf_exp["q_j"], prf_exp["q_l"], prf_exp["We_ql"]))

    # ============================ figure =====================================
    # 2x2. Row 1: q_j (INT v_z r dr) and its We_j = q_j^2/r_j^3.
    #      Row 2: q_l (INT v_z dr)   and its We_j = q_l^2/r_j.
    # The two We_j columns estimate the SAME physical We_j = v_j^2 r_j from two
    # velocity-profile moments; agreement up to a constant confirms a
    # self-similar jet profile.
    fig, axes = plt.subplots(2, 2, figsize=(13.8, 11.4))
    panels = [
        (axes[0][0], "q_j",  r"$q_j = \int v_z\, r\, \mathrm{d}r$", "(a)", False),
        (axes[0][1], "We_qj", r"$We_j = q_j^2 / r_j^3$", "(b)", True),
        (axes[1][0], "q_l",  r"$q_\ell = \int v_z\, \mathrm{d}r$", "(c)", False),
        (axes[1][1], "We_ql", r"$We_j = q_\ell^2 / r_j$", "(d)", True),
    ]

    for ax, key, ylabel, tag, is_we in panels:
        for s in series:
            rt, qt = thin_logspace(s["r_j"], s[key], n=220)
            ax.plot(rt, qt, GRID_MARKERS.get(s["grid"], "o"), ms=5.2,
                    mfc=oh_colour[s["oh"]], mec="k", mew=0.3, lw=0, zorder=3)

        # Each theory line is drawn ONLY over the r_j band where its scaling
        # actually holds (established from the running d ln q/d ln r_j slope of
        # the pooled data), and its prefactor is least-squares fit over that
        # same band. The cone self-similar solution is the r_j->0 asymptote near
        # inception; the PRF 2023 scalings are a finite intermediate band.
        # Drawn spans are extended well past each fit window so the SLOPE is
        # legible (a scaling line needs length to read); the prefactor is still
        # fit only within the window. The PRF line in particular must reach far
        # enough on both sides of its crossing to separate from the data cloud
        # and show the r_j^-1 (or const) slope.
        cone_draw = (fit_rlo / 1.4, fit_rhi * 3.5)    # near-inception, extended
        prf_draw = (prf_rlo / 2.1, prf_rhi * 3.1)     # finite PRF band, extended

        # cone-theory lines, per Oh (slope fixed by the fitted alpha; prefactor
        # least-squares over the near-inception window [fit_rlo, fit_rhi])
        for oh in ohs:
            cone = alpha_by_oh[oh]
            if cone is None:
                continue
            col = oh_colour[oh]
            r_all = np.concatenate([s["r_j"] for s in series if s["oh"] == oh])
            q_all = np.concatenate([s[key] for s in series if s["oh"] == oh])
            rline = np.geomspace(*cone_draw, 40)
            s_our = exps(cone["alpha"])[key]
            K_our, _ = powerfit_prefactor(r_all, q_all, s_our, fit_rlo, fit_rhi)
            if K_our is not None:
                ax.plot(rline, K_our * rline ** s_our, "-", color=col, lw=2.0,
                        alpha=0.9, zorder=4)

        # PRF 2023 (Gordillo & Rodriguez-Rodriguez, La>2500 constant-far-field-
        # flux branch) = the same formula with alpha=1/2: q_l = const (r_j^0),
        # We_j ~ r_j^-1, q_j ~ r_j. Dotted, per Oh, prefactor fit over the FINITE
        # PRF window [prf_rlo, prf_rhi] where those scalings hold (the q_l
        # plateau / We_j ~ r_j^-1 crossing), NOT the r_j->0 cone window. The two
        # Jose asked for are (c) [q_l = const] and (b),(d) [We_j ~ r_j^-1].
        for oh in ohs:
            if alpha_by_oh[oh] is None:
                continue
            col = oh_colour[oh]
            r_all_oh = np.concatenate([s["r_j"] for s in series if s["oh"] == oh])
            q_all_oh = np.concatenate([s[key] for s in series if s["oh"] == oh])
            rline = np.geomspace(*prf_draw, 40)
            K_prf, _ = powerfit_prefactor(r_all_oh, q_all_oh, prf_exp[key],
                                          prf_rlo, prf_rhi)
            if K_prf is not None:
                ax.plot(rline, K_prf * rline ** prf_exp[key], ls=(0, (1, 1.1)),
                        color=col, lw=2.8, alpha=0.95, zorder=5)

        r_all = np.concatenate([s["r_j"] for s in series])
        q_all = np.concatenate([s[key] for s in series])

        if is_we:
            # inertio-capillary: We_j = O(1), i.e. r_j-independent AND order
            # unity (NOT fit to the data). Drawn as a single horizontal line at
            # We_j = 1 — the physical claim the data is meant to refute.
            ax.axhline(1.0, ls="--", color="0.35", lw=2.0, zorder=1)
            ylo = min(q_all.min() * 0.7, 0.6)
            yhi = q_all.max() * 1.5
        else:
            # q_j inertio-capillary: still a power law (alpha=2/3 -> r^1.5),
            # a near-inception competitor to the cone, so fit + drawn over the
            # same near-inception window, per Oh.
            for oh in ohs:
                if alpha_by_oh[oh] is None:
                    continue
                col = oh_colour[oh]
                r_all_oh = np.concatenate([s["r_j"] for s in series if s["oh"] == oh])
                q_all_oh = np.concatenate([s[key] for s in series if s["oh"] == oh])
                rline = np.geomspace(*cone_draw, 40)
                K_ic, _ = powerfit_prefactor(r_all_oh, q_all_oh, ic_exp[key], fit_rlo, fit_rhi)
                if K_ic is not None:
                    ax.plot(rline, K_ic * rline ** ic_exp[key], "--", color=col,
                            lw=2.0, alpha=0.9, zorder=4)
            ylo = q_all.min() * 0.7
            yhi = q_all.max() * 1.4

        style_log_axis(ax, (r_all.min() * 0.85, r_all.max() * 1.2), (ylo, yhi))
        ax.set_xlabel(r"$r_j$", fontsize=19, labelpad=5)
        ax.set_ylabel(ylabel, fontsize=18, labelpad=7)
        ax.set_title(tag, loc="left", fontsize=17, pad=8)

    ax_a, ax_b = axes[0][0], axes[0][1]
    ax_c, ax_d = axes[1][0], axes[1][1]
    a_ref = next((alpha_by_oh[o]["alpha"] for o in ohs if alpha_by_oh[o]), None)

    # ---- legends -----------------------------------------------------------
    # (a): Oh colours + grid markers + Bo (Bo is a legend entry so the
    # forthcoming Bo=0 runs can be distinguished from Bo=1e-3).
    def oh_label(oh):
        cone = alpha_by_oh[oh]
        if cone is not None:
            return r"$Oh = %.4g,\ \beta = %.1f^\circ$" % (oh, cone["beta"])
        return r"$Oh = %.4g$" % oh
    oh_handles = [Line2D([0], [0], marker="o", ls="", mfc=oh_colour[oh], mec="k",
                          mew=0.3, ms=8, label=oh_label(oh)) for oh in ohs]
    grids = sorted({s["grid"] for s in series})
    grid_handles = [Line2D([0], [0], marker=GRID_MARKERS.get(g, "o"), ls="", mfc="0.6",
                           mec="k", mew=0.3, ms=8, label=r"L%d" % g) for g in grids]
    bo_handle = [Line2D([0], [0], ls="", label=r"$Bo = 10^{-3}$")]
    leg1 = ax_a.legend(handles=bo_handle + oh_handles + grid_handles, fontsize=11.5,
                       loc="lower right", frameon=False, handletextpad=0.4,
                       labelspacing=0.32)
    ax_a.add_artist(leg1)

    cone_lbl = (r"cone ($\alpha\!\approx\!%.2f$)" % a_ref) if a_ref else "cone"
    # q-panel line legend (power-law inertio-capillary, alpha=2/3).
    # (a) q_j is small at small r_j -> upper-left is clear.
    # (c) q_l is O(1) even at small r_j (non-monotonic) -> lower-left is clear.
    def q_handles(prf_lbl):
        return [
            Line2D([0], [0], color="0.3", ls="-", lw=2.0, label=cone_lbl),
            Line2D([0], [0], color="0.3", ls="--", lw=2.0,
                   label=r"inertio-capillary ($\alpha=2/3$)"),
            Line2D([0], [0], color="0.3", ls=":", lw=2.3, label=prf_lbl),
        ]
    # PRF 2023 = alpha=1/2. On (a) [q_j] that is q_j ~ r_j; on (c) [q_l] it is
    # the q_l = const claim Jose asked for (the flat dotted line).
    ax_a.legend(handles=q_handles(r"PRF 2023 ($\alpha=1/2$)"), fontsize=11,
                loc="upper left", frameon=False, handletextpad=0.6, labelspacing=0.3)
    ax_c.legend(handles=q_handles(r"PRF 2023 ($q_\ell=$ const)"), fontsize=11,
                loc="lower left", frameon=False, handletextpad=0.6, labelspacing=0.3)
    # We-panel line legend (inertio-capillary We_j = O(1); PRF 2023 We_j ~
    # r_j^-1): on (b) and (d)
    for ax in (ax_b, ax_d):
        ax.legend(handles=[
            Line2D([0], [0], color="0.3", ls="-", lw=2.0, label=cone_lbl),
            Line2D([0], [0], color="0.35", ls="--", lw=2.0,
                   label=r"inertio-capillary ($We_j=O(1)$)"),
            Line2D([0], [0], color="0.3", ls=":", lw=2.3,
                   label=r"PRF 2023 ($We_j \propto r_j^{-1}$)"),
        ], fontsize=11, loc="lower left", frameon=False, handletextpad=0.6, labelspacing=0.3)

    plt.tight_layout(w_pad=2.4, h_pad=2.2)
    plt.savefig(args.out + ".png", dpi=300, bbox_inches="tight", pad_inches=0.08)
    plt.savefig(args.out + ".pdf", dpi=300, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    print("\nWROTE: %s.png\nWROTE: %s.pdf" % (args.out, args.out))


if __name__ == "__main__":
    main()
