"""Young-Laplace equilibrium cavity for a bubble at a free surface.

Length is scaled by the equivalent-sphere radius R0 (gas volume 4π/3).
Bond is Bo = ρ g R0² / γ. That is the Bursting-Bubble / MATLAB
``InitialCondition.m`` convention. It is not the Aberny sandbox
convention, which uses the bottom curvature radius as the length.

The submerged profile (R(φ), Z(φ)) and the outer meniscus are shot
together so that:

* the spherical cap matches slope and curvature at the contact
* the meniscus approaches Z → h∞ as R → ∞
* the enclosed gas volume equals the unit equivalent sphere

φ is the tangent angle from the horizontal, measured from the south
pole (φ = 0) toward the top (φ = π). The opening angle in Lhuissier
& Villermaux (2012) is α_c = π − φ_c.

Basilisk ``Bo%5.4f.dat`` files use (axial, radial) with the cavity in
−axial, the far-field free surface at axial = 0, and a circular fillet
replacing the curvature singularity at the contact.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

__all__ = [
    "EquilibriumShape",
    "SolveError",
    "bond_filename",
    "continuation_ladder",
    "rb_bracket",
    "solve_equilibrium",
    "sweep_bonds",
    "write_basilisk_dat",
]


class SolveError(RuntimeError):
    """Raised when the nested Bond-number shooting fails."""


def _trapz(y, x):
    if hasattr(np, "trapezoid"):
        return float(np.trapezoid(y, x))
    return float(np.trapz(y, x))


def bond_filename(bond: float) -> str:
    """Match MATLAB ``sprintf('Bo%5.4f.dat', Bond)`` and the solvers."""
    return f"Bo{bond:5.4f}.dat"


def continuation_ladder(
    target: float,
    seed: float = 1e-3,
    max_ratio: float = 2.0,
) -> list[float]:
    """Increasing Bond numbers from a cheap seed up to ``target``.

    Small-Bo shapes are nearly spherical (Rb ≈ 1). Walking Bo upward,
    using each hit as the next (Rb, φ_c) guess, is much more stable than
    a cold start at large Bond with a guessed RbMax.
    """
    target = float(target)
    seed = float(seed)
    if target <= seed * 1.01:
        return [target]
    n = max(1, int(np.ceil(np.log(target / seed) / np.log(max_ratio))))
    raw = np.geomspace(seed, target, n + 1)
    raw[-1] = target
    return [float(x) for x in raw]


def rb_bracket(bond: float, previous_rb: Optional[float] = None) -> tuple[float, float]:
    """Bond-dependent bottom-curvature bracket.

    The MATLAB driver never shipped a table of initial Rc; each working
    copy used a different ``RbGuess`` (``sqrt(2)`` at Bo = 10^{-3},
    ``sqrt(25)`` at Bo = 6, ``sqrt(1000)`` at Bo = 10). Continuation
    from a previous solution is preferred; otherwise the upper bound
    grows with Bond so the volume residual can change sign.
    """
    if previous_rb is not None and np.isfinite(previous_rb) and previous_rb > 0.0:
        lo = max(0.45, 0.55 * previous_rb)
        hi = max(previous_rb * 1.85, previous_rb + 0.35)
        return lo, hi
    lo = 0.70
    hi = 1.15 + 2.6 * (bond ** 0.38) + 0.35 * bond ** 0.15
    return lo, max(hi, 1.45)


def _bubble_rhs(phi, y, bond, rb):
    R, Z = y
    R = float(np.clip(R, 1e-14, 80.0))
    Z = float(np.clip(Z, -8.0, 80.0))
    rb = max(float(rb), 1e-8)
    den = R * (2.0 / rb + bond * Z) - np.sin(phi)
    if abs(den) < 1e-18:
        den = 1e-18
    dR = R * np.cos(phi) / den
    # sinφ form, not tanφ·dR/dφ, so the equator is regular
    dZ = R * np.sin(phi) / den
    return (dR, dZ)


def _tail_rhs(R, y, bond, hinf):
    Z, slope = y
    R = max(float(R), 1e-14)
    f2 = 1.0 + slope * slope
    dZ = slope
    dslope = bond * (Z - hinf) * f2 ** 1.5 - (slope / R) * f2
    return (dZ, dslope)


def _interp_state(phi, y, phic):
    """Linear interpolate (R, Z) at φ_c. MATLAB used a broken denominator."""
    phi = np.asarray(phi)
    y = np.asarray(y)
    if phic <= phi[0]:
        return float(y[0, 0]), float(y[0, 1]), 0
    if phic >= phi[-1]:
        return float(y[-1, 0]), float(y[-1, 1]), int(len(phi) - 1)
    ind = int(np.searchsorted(phi, phic) - 1)
    den = phi[ind + 1] - phi[ind]
    ratio = 0.0 if den == 0.0 else (phic - phi[ind]) / den
    R = y[ind, 0] + ratio * (y[ind + 1, 0] - y[ind, 0])
    Z = y[ind, 1] + ratio * (y[ind + 1, 1] - y[ind, 1])
    return float(R), float(Z), ind


def _event_bubble_leave(_phi, y, _bond, _rb):
    return 8.0 - max(abs(y[0]), abs(y[1]))


_event_bubble_leave.terminal = True


def _event_tail_leave(_R, y, _bond, hinf):
    zlim = 5.0 * max(abs(hinf), 4.0) + 10.0
    slim = 1.0e5
    return min(zlim - abs(y[0]), slim - abs(y[1]))


_event_tail_leave.terminal = True


def integrate_bubble(bond: float, rb: float, n_eval: int = 1200):
    """Integrate the submerged cavity from a regularised south-pole start."""
    phi0 = 1e-5
    y0 = (rb * np.sin(phi0), rb * (1.0 - np.cos(phi0)))
    sol = solve_ivp(
        _bubble_rhs,
        (phi0, np.pi - 1e-8),
        y0,
        method="RK45",
        args=(bond, rb),
        rtol=1e-7,
        atol=1e-9,
        max_step=0.05,
        events=_event_bubble_leave,
        dense_output=True,
    )
    if sol.y.size == 0 or sol.sol is None:
        raise SolveError(f"bubble ODE produced no samples (Bo={bond}, Rb={rb})")
    phi_end = float(sol.t[-1])
    if phi_end < 0.5 * np.pi + 0.05:
        raise SolveError(f"bubble ODE stalled at φ={phi_end:.4f} (Bo={bond}, Rb={rb})")
    phi = np.linspace(phi0, phi_end, n_eval)
    y = sol.sol(phi).T
    finite = np.isfinite(y).all(axis=1) & (y[:, 0] > 0.0)
    if np.count_nonzero(finite) < 20:
        raise SolveError(f"bubble ODE left the physical region (Bo={bond}, Rb={rb})")
    return phi[finite], y[finite], sol


def integrate_tail(bond, hinf, r0, z0, slope0, rmax, n_eval=800):
    if rmax <= r0 + 1e-6 or not np.isfinite(slope0):
        return None
    sol = solve_ivp(
        _tail_rhs,
        (r0, rmax),
        (z0, slope0),
        method="RK45",
        args=(bond, hinf),
        rtol=1e-7,
        atol=1e-9,
        max_step=max((rmax - r0) / 80.0, 5e-3),
        events=_event_tail_leave,
        dense_output=True,
    )
    if sol.y.size == 0 or sol.sol is None:
        return None
    if sol.t[-1] < r0 + 1e-3:
        return None
    R = np.linspace(r0, float(sol.t[-1]), n_eval)
    Z = sol.sol(R)[0]
    if not np.all(np.isfinite(Z)):
        return None
    return R, Z


def _spherical_cap(xc, yc, phic, n=800):
    Rc = xc / np.sin(phic)
    R = np.linspace(xc, 0.0, n)
    root = Rc * Rc - R * R
    root[root < 0.0] = 0.0
    Z = yc + xc / np.tan(phic) + np.sqrt(root)
    return R, Z, Rc


def _volume(R_sub, Z_sub, R_cap, Z_cap) -> float:
    vol1 = 0.75 * _trapz(R_sub ** 2, Z_sub)
    vol2 = 0.75 * _trapz(R_cap ** 2, Z_cap)
    return vol1 + vol2


def _tail_rmax(bond: float, requested: Optional[float]) -> float:
    """Finite radius at which the meniscus is asked to meet h∞.

    MATLAB ``InitialCondition.m`` shoots to ``TailxMax = 8`` (sometimes 4)
    and only afterwards pads the free surface out to 32. Checking at
    infinity is neither necessary nor stable.
    """
    if requested is not None:
        return float(requested)
    return float(np.clip(6.0 + 2.0 / max(np.sqrt(bond), 0.25), 6.0, 12.0))


def _phi_bracket(
    bond: float,
    previous_phic: Optional[float] = None,
) -> tuple[float, float]:
    # Lhuissier small-Bo estimate α_c ≈ (1/(2√3)) sqrt(Bo) when Rc ~ R0.
    # φ_c = π − α_c sits very close to π at small Bond, so the bracket
    # must open near the pole rather than at the equator.
    hard_lo = 0.5 * np.pi + 1e-4
    hard_hi = np.pi - 1e-5
    if previous_phic is not None and np.isfinite(previous_phic):
        # φ_c falls as Bo rises; bias the window downward from the last hit.
        lo = max(hard_lo, float(previous_phic) - 0.40)
        hi = min(hard_hi, float(previous_phic) + 0.12)
        if hi - lo < 0.05:
            lo = max(hard_lo, hi - 0.20)
        return lo, hi
    alpha_est = 0.5 * np.sqrt(max(bond, 0.0) / 3.0)
    window = max(8.0 * max(alpha_est, 0.02), 0.20)
    lo = max(hard_lo, np.pi - window)
    lo = min(lo, np.pi - 5e-3)
    return lo, hard_hi


def _match_radius(bond, xc, rmax):
    """Station used for the meniscus residual.

    The decaying meniscus mode is swamped by the growing mode beyond a
    few capillary lengths, so Z(R_max) − h∞ is discontinuous in φ_c.
    Matching at xc + O(1/sqrt(Bo)) keeps the residual smooth.
    """
    cap = 2.5 / max(np.sqrt(bond), 1e-3)
    return float(min(rmax, xc + np.clip(cap, 1.5, 6.0)))


def _shoot_tail(bond, phi, y, rb, rmax, tail_tol, previous_phic=None):
    phi_lo, phi_hi = _phi_bracket(bond, previous_phic=previous_phic)
    phi_lo = max(phi_lo, float(phi[0]) + 1e-4)
    phi_hi = min(phi_hi, float(phi[-1]) - 1e-4)
    if phi_lo >= phi_hi:
        raise SolveError(f"empty φ_c bracket (Bo={bond})")

    def residual(phic):
        xc, yc, _ = _interp_state(phi, y, phic)
        if xc <= 0.0 or not np.isfinite(xc):
            return 1e3
        Rc = xc / np.sin(phic)
        hinf = (2.0 / bond) * (2.0 / Rc - 1.0 / rb)
        slope = np.tan(phic)
        if not np.isfinite(slope) or abs(slope) > 1e6:
            return 1e3 if phic < 0.75 * np.pi else -1e3
        r_match = _match_radius(bond, xc, rmax)
        tail = integrate_tail(bond, hinf, xc, yc, slope, r_match)
        if tail is None:
            return 1e3
        return float(tail[1][-1] - hinf)

    r_lo = residual(phi_lo)
    r_hi = residual(phi_hi)
    if not (np.isfinite(r_lo) and np.isfinite(r_hi)):
        raise SolveError(f"tail residual is not finite (Bo={bond})")
    if r_lo * r_hi > 0.0 or abs(r_lo) + abs(r_hi) > 0.5:
        # Prefer the sign-change whose residuals are actually small.
        # A growing-mode jump is also a sign change but |f| stays O(1).
        grid = np.linspace(phi_lo, phi_hi, 28)
        vals = [residual(p) for p in grid]
        best = None
        for a, b, fa, fb in zip(grid[:-1], grid[1:], vals[:-1], vals[1:]):
            if not (np.isfinite(fa) and np.isfinite(fb)):
                continue
            if fa * fb > 0.0:
                continue
            score = abs(fa) + abs(fb)
            if best is None or score < best[0]:
                best = (score, a, b)
        if best is None:
            raise SolveError(
                f"tail residual does not change sign (Bo={bond}, "
                f"r({phi_lo:.4f})={r_lo:.3e}, r({phi_hi:.4f})={r_hi:.3e})"
            )
        phi_lo, phi_hi = best[1], best[2]

    phic = brentq(residual, phi_lo, phi_hi, xtol=max(tail_tol * 0.1, 1e-12))
    xc, yc, ind = _interp_state(phi, y, phic)
    Rc = xc / np.sin(phic)
    hinf = (2.0 / bond) * (2.0 / Rc - 1.0 / rb)
    tail_err = abs(residual(phic))
    tail = integrate_tail(bond, hinf, xc, yc, np.tan(phic), rmax)
    if tail is None:
        raise SolveError(f"tail ODE failed at the converged φ_c (Bo={bond})")
    R_t, Z_t = tail
    # Keep the near-field meniscus (Z is not yet h∞), then drop from the
    # first growing-mode spike after the profile has approached h∞.
    far = np.isfinite(Z_t) & (np.abs(Z_t - hinf) < 0.35 * max(abs(hinf), 1.0) + 0.25)
    entered = np.flatnonzero(far)
    if entered.size >= 8:
        start_far = int(entered[0])
        rest = far[start_far:]
        cut = len(Z_t) if bool(rest.all()) else start_far + int(np.argmin(rest))
        if cut >= 8:
            R_t = R_t[:cut]
            Z_t = Z_t[:cut]
    tail = (R_t, Z_t)
    return {
        "phic": float(phic),
        "xc": float(xc),
        "yc": float(yc),
        "ind": int(ind),
        "Rc": float(Rc),
        "hinf": float(hinf),
        "R_tail": tail[0],
        "Z_tail": tail[1],
        "tail_err": tail_err,
    }


def _cut_index_by_arclength(R, Z, span, from_start=True):
    ds = np.sqrt(np.diff(R) ** 2 + np.diff(Z) ** 2)
    s = np.concatenate([[0.0], np.cumsum(ds)])
    if s[-1] <= span:
        return (1 if from_start else len(R) - 2)
    if from_start:
        return int(np.clip(np.searchsorted(s, span), 1, len(R) - 2))
    return int(np.clip(len(R) - 1 - np.searchsorted(s, span), 1, len(R) - 2))


def _blend_fillet(R_sub, Z_sub, phi_sub, R_tail, Z_tail, fillet_span):
    """C¹ Hermite blend between cavity and meniscus.

    MATLAB's circular-fillet radius can change sign with the sample
    index. A cubic matching position and tangent at both cuts is
    unambiguous and keeps the neck a single forward march.
    """
    if fillet_span <= 0.0:
        return None
    j = _cut_index_by_arclength(R_tail, Z_tail, fillet_span, from_start=True)
    i = _cut_index_by_arclength(R_sub[::-1], Z_sub[::-1], fillet_span, from_start=True)
    i = len(R_sub) - 1 - i
    i = int(np.clip(i, 1, len(R_sub) - 2))
    j = int(np.clip(j, 1, len(R_tail) - 2))
    p0 = np.array([float(R_sub[i]), float(Z_sub[i])])
    p1 = np.array([float(R_tail[j]), float(Z_tail[j])])
    phif = float(phi_sub[i])
    t0 = np.array([np.cos(phif), np.sin(phif)])
    t1 = np.array([float(R_tail[j + 1] - R_tail[j - 1]),
                   float(Z_tail[j + 1] - Z_tail[j - 1])])
    n1 = np.linalg.norm(t1)
    if n1 < 1e-14:
        t1 = np.array([1.0, 0.0])
    else:
        t1 = t1 / n1
    chord = np.linalg.norm(p1 - p0)
    if chord < 1e-12:
        return None
    u = np.linspace(0.0, 1.0, 400)[:, None]
    h00 = 2 * u ** 3 - 3 * u ** 2 + 1
    h10 = u ** 3 - 2 * u ** 2 + u
    h01 = -2 * u ** 3 + 3 * u ** 2
    h11 = u ** 3 - u ** 2
    pts = h00 * p0 + h10 * (chord * t0) + h01 * p1 + h11 * (chord * t1)
    return {
        "i_cut": i,
        "j_cut": j,
        "R": pts[:, 0],
        "Z": pts[:, 1],
        "r": 0.25 * chord,
    }


@dataclass
class EquilibriumShape:
    bond: float
    Rb: float
    Rc: float
    phic: float
    hinf: float
    volume: float
    volume_err: float
    tail_err: float
    R_bubble: np.ndarray
    Z_bubble: np.ndarray
    R_cap: np.ndarray
    Z_cap: np.ndarray
    R_tail: np.ndarray
    Z_tail: np.ndarray
    axial: np.ndarray
    radial: np.ndarray
    fillet_radius: float = 0.0
    notes: dict = field(default_factory=dict)

    @property
    def alpha_c(self) -> float:
        return float(np.pi - self.phic)

    @property
    def opening_metric(self) -> float:
        """2 α_c / π, the Lhuissier–Villermaux ordinate."""
        return 2.0 * self.alpha_c / np.pi

    @property
    def capillary_metric(self) -> float:
        """(Rc/R0) sqrt(Bo) = Rc/a, the Lhuissier–Villermaux abscissa."""
        return self.Rc * np.sqrt(self.bond)

    def south_pole_axial(self) -> float:
        return float(self.axial[0])


def _solve_at_bond(
    bond: float,
    previous: Optional[EquilibriumShape] = None,
    *,
    tail_rmax: Optional[float] = None,
    rmax_out: float = 32.0,
    fillet_span: float = 0.22,
    vol_tol: float = 1e-8,
    tail_tol: float = 1e-6,
    max_vol_iter: int = 60,
) -> EquilibriumShape:
    """One Bond-number hit, optionally seeded by the previous shape."""
    if bond <= 0.0:
        raise ValueError("Bo > 0 is required; use generate_zero_bond.py for Bo = 0")
    if bond < 1e-6:
        raise ValueError(
            f"Bo = {bond} is too small for the 2/Bo meniscus match; "
            "use the zero-Bond sphere-plane generator"
        )

    rmax = _tail_rmax(bond, tail_rmax)
    prev_rb = None if previous is None else previous.Rb
    prev_phic = None if previous is None else previous.phic
    lo, hi = rb_bracket(bond, prev_rb)

    def packed(rb):
        phi, y, _ = integrate_bubble(bond, rb)
        tail = _shoot_tail(
            bond, phi, y, rb, rmax, tail_tol, previous_phic=prev_phic
        )
        R_sub = y[: tail["ind"] + 1, 0]
        Z_sub = y[: tail["ind"] + 1, 1]
        R_cap, Z_cap, _ = _spherical_cap(tail["xc"], tail["yc"], tail["phic"])
        vol = _volume(R_sub, Z_sub, R_cap, Z_cap)
        return phi, y, tail, R_sub, Z_sub, R_cap, Z_cap, vol

    def volume_mismatch(rb):
        try:
            *_, vol = packed(rb)
        except (SolveError, ValueError, FloatingPointError):
            return None
        return vol - 1.0

    def first_value(rb):
        val = volume_mismatch(rb)
        return val

    vlo = first_value(lo)
    vhi = first_value(hi)
    expand = 0
    while expand < 10:
        if vlo is not None and vhi is not None and vlo * vhi <= 0.0:
            break
        if vlo is not None and vhi is not None and vlo < 0.0 and vhi < 0.0:
            lo, vlo = hi, vhi
            hi = hi * 1.8
            vhi = first_value(hi)
        elif vlo is not None and vhi is not None and vlo > 0.0 and vhi > 0.0:
            hi, vhi = lo, vlo
            lo = max(0.40, lo * 0.7)
            vlo = first_value(lo)
        else:
            # one or both ends failed: probe a wider log grid
            probes = np.geomspace(max(lo, 0.5), max(hi, 1.2) * (1.6 ** (expand + 1)), 6)
            found = False
            last_ok = None
            for rb_p in probes:
                val = first_value(float(rb_p))
                if val is None:
                    continue
                if last_ok is not None and last_ok[1] * val <= 0.0:
                    lo, vlo = last_ok
                    hi, vhi = float(rb_p), val
                    found = True
                    break
                last_ok = (float(rb_p), val)
            if not found:
                if last_ok is not None:
                    lo, vlo = last_ok
                    hi = lo * 1.8
                    vhi = first_value(hi)
                else:
                    hi *= 1.8
                    vhi = first_value(hi)
        expand += 1
    if vlo is None or vhi is None or vlo * vhi > 0.0:
        raise SolveError(
            f"volume residual does not change sign (Bo={bond}, "
            f"Rb in [{lo:.4g}, {hi:.4g}], ΔV=({vlo}, {vhi}))"
        )

    def volume_for_root(rb):
        val = volume_mismatch(rb)
        if val is None:
            raise SolveError(
                f"volume residual evaluation failed (Bo={bond}, Rb={rb:.6g})"
            )
        return val

    rb = brentq(volume_for_root, lo, hi, xtol=max(vol_tol, 1e-12), maxiter=max_vol_iter)
    phi, y, tail, R_sub, Z_sub, R_cap, Z_cap, vol = packed(rb)
    tail_ok = max(50.0 * tail_tol, 5e-2 * max(abs(tail["hinf"]), 1.0), 0.12)
    if tail["tail_err"] > tail_ok:
        raise SolveError(f"tail residual {tail['tail_err']:.3e} after volume match (Bo={bond})")
    if abs(vol - 1.0) > max(100.0 * vol_tol, 1e-5):
        raise SolveError(f"volume residual {vol-1:.3e} after match (Bo={bond})")

    R_tail, Z_tail = tail["R_tail"], tail["Z_tail"]
    # Hold the far field exactly at h∞ once the shooting has converged.
    if R_tail[-1] < rmax_out:
        R_tail = np.concatenate([R_tail, np.linspace(R_tail[-1], rmax_out, 400)[1:]])
        Z_tail = np.concatenate([Z_tail, np.full(399, tail["hinf"])])
    else:
        Z_tail = Z_tail.copy()
        Z_tail[-1] = tail["hinf"]

    phi_sub = phi[: tail["ind"] + 1]
    fillet = _blend_fillet(R_sub, Z_sub, phi_sub, R_tail, Z_tail, fillet_span)
    z_shift = Z_tail[-1]
    if fillet is None:
        R_poly = np.concatenate([R_sub, R_tail])
        Z_poly = np.concatenate([Z_sub, Z_tail])
        fillet_r = 0.0
    else:
        R_poly = np.concatenate(
            [R_sub[: fillet["i_cut"] + 1], fillet["R"], R_tail[fillet["j_cut"] :]]
        )
        Z_poly = np.concatenate(
            [Z_sub[: fillet["i_cut"] + 1], fillet["Z"], Z_tail[fillet["j_cut"] :]]
        )
        fillet_r = fillet["r"]

    axial = Z_poly - z_shift
    radial = np.maximum(R_poly, 1e-8)
    radial[0] = 1e-8

    return EquilibriumShape(
        bond=float(bond),
        Rb=float(rb),
        Rc=float(tail["Rc"]),
        phic=float(tail["phic"]),
        hinf=float(tail["hinf"]),
        volume=float(vol),
        volume_err=float(vol - 1.0),
        tail_err=float(tail["tail_err"]),
        R_bubble=R_sub,
        Z_bubble=Z_sub,
        R_cap=R_cap,
        Z_cap=Z_cap,
        R_tail=R_tail,
        Z_tail=Z_tail,
        axial=axial,
        radial=radial,
        fillet_radius=float(fillet_r),
        notes={
            "tail_rmax": rmax,
            "fillet_span": fillet_span,
            "rb_bracket": (lo, hi),
            "bracket_expansions": expand,
            "continued_from": None if previous is None else previous.bond,
        },
    )


def solve_equilibrium(
    bond: float,
    previous: Optional[EquilibriumShape] = None,
    *,
    continue_in_bond: bool = True,
    continuation_seed: float = 1e-3,
    continuation_max_ratio: float = 2.0,
    **kwargs,
) -> EquilibriumShape:
    """Solve the floating-bubble Young-Laplace system at one Bond number.

    When ``continue_in_bond`` is true, Bond is walked from a cheap seed
    (or from ``previous``) up to the target. Each accepted shape seeds
    the next (Rb, φ_c) brackets, so a cold start at large Bo does not
    need a guessed RbMax.
    """
    if previous is not None and previous.bond > float(bond) * 1.01:
        previous = None
    if not continue_in_bond:
        return _solve_at_bond(bond, previous=previous, **kwargs)

    seed = continuation_seed
    if previous is not None:
        seed = max(seed, float(previous.bond))
    ladder = continuation_ladder(
        float(bond), seed=min(seed, float(bond)), max_ratio=continuation_max_ratio
    )

    current = previous
    walked = []
    for step in ladder:
        if current is not None and abs(
            np.log(max(step, 1e-30) / max(current.bond, 1e-30))
        ) < np.log(1.02):
            continue
        current = _solve_at_bond(step, previous=current, **kwargs)
        walked.append(step)
    if current is None:
        raise SolveError(f"continuation produced no shape (Bo={bond})")
    current.notes["continuation"] = walked
    return current


def sweep_bonds(bonds, previous: Optional[EquilibriumShape] = None, **kwargs):
    """Solve a Bond list in increasing order, continuing from the last hit."""
    shapes = []
    current = previous
    for bond in np.asarray(bonds, dtype=float):
        shape = solve_equilibrium(float(bond), previous=current, **kwargs)
        shapes.append(shape)
        current = shape
    return shapes


def write_basilisk_dat(shape: EquilibriumShape, path) -> None:
    """Write the Stage-1 polyline: column 1 axial, column 2 radial."""
    with open(path, "w") as fh:
        for x, y in zip(shape.axial, shape.radial):
            fh.write("% .7e   % .7e\n" % (x, y))
