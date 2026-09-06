/**
Vendored from comphy-lab/MultiRheoFlow `two-phaseVE.h`, reconciled with
upstream `7d9c3df` (2026-08-30, PR #7 "adapt-basilisk-v2026-08-30").

DELIBERATE DIVERGENCE FROM UPSTREAM. Upstream now calls `set_prolongation()`
unconditionally, which requires Basilisk newer than the `sf.dirty` removal.
This copy keeps the `VE_USE_SET_PROLONGATION` shim so the same source builds
against old and new trees; `append_solver_qcc_flags()` in
`src-local/parse_params.sh` defines the macro by probing the host tree for
`void set_prolongation`. Do not "simplify" this back to upstream without
first confirming every campaign host is on a new enough Basilisk.

Keep the G-lambda properties in lockstep with MultiRheoFlow.
*/

/**
# Two-Phase Viscoelastic Solver

Modified from Basilisk `two-phase.h` and `two-phase-generic.h` to
support viscoelastic flows with log-conformation rheology.

## Change Log

- 2024-10-17: Add support for VE simulations.

## Two-Phase Interfacial Flows

The interface is tracked with VOF. The volume fraction is `f = 1` in
fluid 1 and `f = 0` in fluid 2. Densities and viscosities are `rho1`,
`mu1`, `rho2`, `mu2`.
*/

#include "vof.h"

scalar f[], * interfaces = {f};

double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0.;
double G1 = 0., G2 = 0.; // elastic moduli
double lambda1 = 0., lambda2 = 0.; // relaxation times
double TOLelastic = 1e-2; // tolerance for elastic modulus #TOFIX: this must always be a very small number.

/**
Auxiliary fields define the specific volume $\alpha = 1/\rho$ and the
cell-centered density.
*/

face vector alphav[];
scalar rhov[];
scalar Gpd[];
scalar lambdapd[];

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;
  Gp = Gpd;
  lambda = lambdapd;

  /**
  If the viscosity is non-zero, we need to allocate the face-centered
  viscosity field. */

  mu = new face vector;
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). */

#ifndef rho
# define rho(f) (clamp(f,0.,1.)*(rho1 - rho2) + rho2)
#endif
#ifndef mu
// for Arithmetic mean, use this
# define mu(f)  (clamp(f,0.,1.)*(mu1 - mu2) + mu2)
#endif

/**
We have the option of using some "smearing" of the density/viscosity
jump. */

#ifdef FILTERED
scalar sf[];
#else
# define sf f
#endif

event tracer_advection (i++) {

  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. */

#ifndef sf
#if dimension <= 2
  foreach()
    sf[] = (4.*f[] +
	    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
	    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
#else // dimension == 3
  foreach()
    sf[] = (8.*f[] +
	    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
	    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] +
		f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
		f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
	    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
	    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
#endif
#endif

#if TREE
#ifdef VE_USE_SET_PROLONGATION
  set_prolongation (sf, refine_bilinear);
#else
  sf.prolongation = refine_bilinear;
  sf.dirty = true; // boundary conditions need to be updated
#endif
#endif
}

event properties (i++) {
  
  foreach_face() {
    double ff = (sf[] + sf[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff);
    face vector muv = mu;
    muv.x[] = fm.x[]*mu(ff);
  }

  foreach(){
    rhov[] = cm[]*rho(sf[]);

    Gpd[] = 0.;
    lambdapd[] = 0.;

    if (clamp(sf[], 0., 1.) > TOLelastic){
      Gpd[] += G1*clamp(sf[], 0., 1.);
      lambdapd[] += lambda1*clamp(sf[], 0., 1.);
    }
    if (clamp((1-sf[]), 0., 1.) > TOLelastic){
      Gpd[] += G2*clamp((1-sf[]), 0., 1.);
      lambdapd[] += lambda2*clamp((1-sf[]), 0., 1.);
    }
  }

#if TREE
#ifdef VE_USE_SET_PROLONGATION
  set_prolongation (sf, fraction_refine);
#else
  sf.prolongation = fraction_refine;
  sf.dirty = true; // boundary conditions need to be updated
#endif
#endif
}

/**
## Elastic-wave stability condition

`tension.h` limits the timestep to the capillary-wave period. Nothing limited
it to the *elastic* wave, and that omission is what destroyed the first
viscoelastic campaign: seven runs at `Oh >= 0.024` blew up within one timestep,
kinetic energy jumping from ~4 to between 1e9 and 1e100, always at the
cavity-focus instant.

The mechanism, measured rather than guessed. In the extensional flow at the jet
base the axial conformation reaches `A11 ~ 3e5` — an extension ratio of ~550,
which Oldroyd-B permits because it has no finite extensibility. The polymeric
stress `Gp*A11` then supports a shear wave whose speed is
$$
c_e = \sqrt{\frac{G_p \,\mathrm{tr}\,\mathbf{A}}{\rho}}
$$
because a stretched dumbbell stiffens along its stretch direction: the modulus
governing perturbations is `Gp*A`, not `Gp`. At the moment of failure `c_e`
was ~52, giving `Delta/c_e = 4.70e-5`, while the solver was stepping at
4.76e-5 — a ratio of 1.014. It was sitting exactly on the stability boundary
and stepped over it.

Restarting that same case with a step five times smaller walked straight
through: the two runs agree to 1 part in 1e4 up to the failing step, after
which one explodes by nine orders of magnitude and the other decays smoothly.

So this is an explicit-scheme CFL condition on a wave family the code did not
account for, not a limitation of Oldroyd-B. Imposing it here makes every
viscoelastic run safe by construction instead of depending on a hand-tuned
`DT` that must be re-guessed whenever `Ec` or the stretch changes.

`CFL_elastic` defaults to 0.25 because that is the margin demonstrated to work
(`dt/dt_limit = 0.213` in the successful restart), not a round number chosen
for looks. `tr(A)` is used rather than the largest eigenvalue: it bounds the
eigenvalue from above, so the criterion errs safe, and it costs no
decomposition. The condition is inert for an unstretched polymer — at
`A = I` it gives a limit far above the capillary one — so it only binds where
the stretch is genuinely large.
*/

double CFL_elastic = 0.25;

event stability (i++)
{
  if (CFL_elastic <= 0.)
    return 0;
  double dtelastic = HUGE;
  foreach (reduction(min:dtelastic)) {
    if (Gp[] > 0.) {
      double trA = A11[] + A22[];
#if AXI
      trA += AThTh[];
#endif
      if (trA > 0. && isfinite(trA)) {
        double rhom = rho(f[]);
        if (rhom > 0.) {
          double ce = sqrt (Gp[]*trA/rhom);
          if (ce > 0.) {
            double dte = CFL_elastic*Delta/ce;
            if (dte < dtelastic) dtelastic = dte;
          }
        }
      }
    }
  }
  if (dtelastic < dtmax)
    dtmax = dtelastic;
}

/**
## Conformation-source stability condition

The elastic-wave condition above bounds the speed at which a *stretched*
polymer transmits information. It says nothing about how fast the
log-conformation update itself is allowed to change `Psi = log A`, and that is
a separate explicit source. `log-conform-viscoelastic-scalar-2D.h` advances
`Psi` with forward Euler,
$$
\Psi^{n+1} = \Psi^n + \Delta t\,\bigl[2\mathbf{B} + (\Omega\Psi - \Psi\Omega)\bigr],
\qquad
\Psi_{\theta\theta}^{n+1} = \Psi_{\theta\theta}^n + \Delta t\,\frac{2u_r}{r},
$$
and then exponentiates. An increment of order unity therefore multiplies the
conformation by `e` in a single step. The relaxation is not the hazard: it is
integrated analytically (`intFactor = exp(-dt/lambda)`), so it is
unconditionally stable and small `lambda` is harmless in itself.

The axisymmetric hoop term is the one nothing bounded. Its rate is `2 u_r / r`,
evaluated in the first cell off the axis at `r = Delta/2`, and the advective
CFL is blind to it: `timestep()` sees the metric-weighted face velocity, which
carries a factor `r` and so vanishes exactly where this rate diverges.

Measured on case 2330 (`De = 0.02`, `Ec = 0.009`, `Oh = 0.024`, level 12), which
blew up at `t = 0.54267` with `ke` going from 4.61 to 4.15e11 in one step. The
peak per-step increment sits in one on-axis interfacial cell at
`(z, r) = (-1.638, Delta/2)` and grows

| t | dt | max abs dPsi_qq |
|---|----|-----------------|
| 0.500 | 2.00e-5 | 0.029 |
| 0.520 | 2.94e-5 | 0.119 |
| 0.540 | 4.00e-5 | 0.123 |
| 0.541 | 4.17e-5 | 0.262 |
| 0.542 | 4.17e-5 | 0.525 |

reaching order unity within the ~16 steps that remained. At the same instants
the elastic-wave condition was SATISFIED and saturated (`dt/dt_elastic` = 0.97
to 0.99), so it neither caused nor caught this.

The `De` dependence is indirect but decisive, and explains why only the lowest
`De` on the line died. At fixed `Ec` a weaker polymer stretches less, so the
elastic-wave limit is looser and `dt` is larger: at `t = 0.542`,
`max tr(A)` = 2.2e4 for `De = 0.02` against 1.7e5 for `De = 0.055` and 1.8e5
for `De = 0.1`, giving `dt` = 4.17e-5 against 1.41e-5. Combined with a ~20x
larger local hoop rate, the resulting increment is
`0.53` (`De = 0.02`) against `0.0078` (`De = 0.055`) and `0.0011`
(`De = 0.1`) — a factor of 480 across the line, monotone in `De`. The elastic
condition was incidentally protecting the higher-`De` cases from a constraint
it does not represent.

`CFL_conform` is that missing bound: the largest log-conformation increment
permitted in one step. `0` disables it, which is the pre-existing behaviour
exactly, so no completed run is affected. `0.1` is the value the measurements
support — it is inert for the cases that already ran (it would have allowed
`dt <= 1.8e-4` for case 2314 at `t = 0.542`, an order of magnitude above the
1.41e-5 that case was taking) and cuts case 2330's step by 5.3x at the instant
it failed. It is a bound on the increment, not a proof of stability; the
velocity-gradient terms are included alongside the hoop term because they are
the same class of explicit source, though only the hoop term is measured here
to have reached order unity.

The condition is restricted to cells carrying polymer (`Gp > 0`), matching the
elastic-wave condition: that is where a corrupted `A` feeds the momentum
equation through `T = Gp (A - I)`.
*/

double CFL_conform = 0.;

event stability (i++)
{
  if (CFL_conform <= 0.)
    return 0;
  double dtconf = HUGE;
  foreach (reduction(min:dtconf)) {
    if (Gp[] <= 0.)
      continue;
    double g = 0.;
#if AXI
    g = fabs (2.*u.y[])/max (y, 1e-20);          // hoop stretch, dPsi_qq
#endif
    double exx = fabs (u.x[1] - u.x[-1])/Delta;  // = |2 du_x/dx|
    double eyy = fabs (u.y[0,1] - u.y[0,-1])/Delta;
    double exy = (fabs (u.x[0,1] - u.x[0,-1]) +
                  fabs (u.y[1] - u.y[-1]))/(2.*Delta);
    if (exx > g) g = exx;
    if (eyy > g) g = eyy;
    if (exy > g) g = exy;
    if (g > 0. && isfinite (g)) {
      double dtg = CFL_conform/g;
      if (dtg < dtconf) dtconf = dtg;
    }
  }
  if (dtconf < dtmax)
    dtmax = dtconf;
}
