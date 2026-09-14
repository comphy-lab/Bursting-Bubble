/**
Vendored from comphy-lab/MultiRheoFlow `log-conform-viscoelastic-scalar-2D.h`, reconciled with upstream
`7d9c3df` (2026-08-30). Keep in lockstep with MultiRheoFlow when the
Oldroyd-B / log-conformation solver changes.

DELIBERATE DIVERGENCE FROM UPSTREAM. This copy aborts on a non-positive or
non-finite axisymmetric hoop conformation `Aqq`, and uses `MPI_Abort` rather
than `exit(1)` so a single bad rank does not leave the job hanging. Upstream
has neither. The guard exists because the high-Weissenberg thinning thread is
exactly where the log-conformation update loses positive-definiteness, and a
silent `log` of a non-positive `Aqq` is worse than a stop. Worth upstreaming;
do not drop it when re-syncing.
*/

/**
# Log-Conformation (Scalar 2D/Axi)

Scalar log-conformation implementation for 2D and axisymmetric
viscoelastic flows.

## Key Features

- Conformation tensor components stored as scalars.
- Supports 2D and axisymmetric configurations.
- Mirrors `log-conform-viscoelastic-scalar-3D.h`.

## Change Log

- 2024-10-18: Initial 2D/axi implementation.
- 2024-11-03: Axisymmetric mirror of 3D scalar version.
- 2024-11-14: Infinite Deborah number support.
- 2024-11-23: Documentation updates.

## Future Work

- Convert to tensor formulation for consistency and maintainability.
- Track: comphy-lab/Viscoelastic3D#11, #5.

## Author

Vatsal Sanjay (vatsal.sanjay@comphy-lab.org)
CoMPhy Lab
Last updated: 2024-11-23
*/

/**
# The log-conformation method for viscoelastic constitutive models

## Introduction

Viscoelastic fluids exhibit both viscous and elastic behaviour when
subjected to deformation. Therefore these materials are governed by
the Navier--Stokes equations enriched with an extra *elastic* stress
$Tij$
$$
\rho\left[\partial_t\mathbf{u}+\nabla\cdot(\mathbf{u}\otimes\mathbf{u})\right] =
- \nabla p + \nabla\cdot(2\mu_s\mathbf{D}) + \nabla\cdot\mathbf{T}
+ \rho\mathbf{a}
$$
where $\mathbf{D}=[\nabla\mathbf{u} + (\nabla\mathbf{u})^T]/2$ is the
deformation tensor and $\mu_s$ is the solvent viscosity of the
viscoelastic fluid.

The *polymeric* stress $\mathbf{T}$ represents memory effects due to
the polymers. Several constitutive rheological models are available in
the literature where the polymeric stress $\mathbf{T}$ is typically a
function $\mathbf{f_s}(\cdot)$ of the conformation tensor $\mathbf{A}$ such as
$$
\mathbf{T} = G_p \mathbf{f_s}(\mathbf{A})
$$
where $G_p$ is the elastic modulus and $\mathbf{f_s}(\cdot)$ is the relaxation function.

The conformation tensor $\mathbf{A}$ is related to the deformation of
the polymer chains. $\mathbf{A}$ is governed by the equation
$$
D_t \mathbf{A} - \mathbf{A} \cdot \nabla \mathbf{u} - \nabla
\mathbf{u}^{T} \cdot \mathbf{A} =
-\frac{\mathbf{f_r}(\mathbf{A})}{\lambda}
$$
where $D_t$ denotes the material derivative and
$\mathbf{f_r}(\cdot)$ is the relaxation function. Here, $\lambda$ is the relaxation time.

In the case of an Oldroyd-B viscoelastic fluid, $\mathbf{f}_s
 (\mathbf{A}) = \mathbf{f}_r (\mathbf{A}) = \mathbf{A} -\mathbf{I}$,
and the above equations can be combined to avoid the use of
$\mathbf{A}$
$$
\mathbf{T} + \lambda (D_t \mathbf{T} -
\mathbf{T} \cdot \nabla \mathbf{u} -
\nabla \mathbf{u}^{T} \cdot \mathbf{T})  = 2 G_p\lambda \mathbf{D}
$$

[Comminal et al. (2015)](#comminal2015) gathered the functions
$\mathbf{f}_s (\mathbf{A})$ and $\mathbf{f}_r (\mathbf{A})$ for
different constitutive models.

## Parameters

The primary parameters are the relaxation time
$\lambda$ and the elastic modulus $G_p$. The solvent viscosity
$\mu_s$ is defined in the [Navier-Stokes
solver](navier-stokes/centered.h).

Gp and lambda are defined in [two-phaseVE.h](two-phaseVE.h).
*/

/**
## The log conformation approach

The numerical resolution of viscoelastic fluid problems often faces the
[High-Weissenberg Number
Problem](http://www.ma.huji.ac.il/~razk/iWeb/My_Site/Research_files/Visco1.pdf).
This is a numerical instability appearing when strongly elastic flows
create regions of high stress and fine features. This instability
poses practical limits to the values of the relaxation time of the
viscoelastic fluid, $\lambda$.  [Fattal \& Kupferman (2004,
2005)](#fattal2004) identified the exponential nature of the solution
as the origin of the instability. They proposed to use the logarithm
of the conformation tensor $\Psi = \log \, \mathbf{A}$ rather than the
viscoelastic stress tensor to circumvent the instability.

The constitutive equation for the log of the conformation tensor is
$$
D_t \Psi = (\Omega \cdot \Psi -\Psi \cdot \Omega) + 2 \mathbf{B} +
\frac{e^{-\Psi} \mathbf{f}_r (e^{\Psi})}{\lambda}
$$
where $\Omega$ and $\mathbf{B}$ are tensors that result from the
decomposition of the transpose of the tensor gradient of the
velocity
$$
(\nabla \mathbf{u})^T = \Omega + \mathbf{B} + N
\mathbf{A}^{-1}
$$

The antisymmetric tensor $\Omega$ requires only the memory of a scalar
in 2D since,
$$
\Omega = \left(
\begin{array}{cc}
0 & \Omega_{12} \\
-\Omega_{12} & 0
\end{array}
\right)
$$

For 3D, $\Omega$ is a skew-symmetric tensor given by

$$
\Omega = \left(
\begin{array}{ccc}
0 & \Omega_{12} & \Omega_{13} \\
-\Omega_{12} & 0 & \Omega_{23} \\
-\Omega_{13} & -\Omega_{23} & 0
\end{array}
\right)
$$

The log-conformation tensor, $\Psi$, is related to the
polymeric stress tensor $\mathbf{T}$, by the strain function
$\mathbf{f}_s (\mathbf{A})$
$$
\Psi = \log \, \mathbf{A} \quad \mathrm{and} \quad \mathbf{T} =
\frac{G_p}{\lambda} \mathbf{f}_s (\mathbf{A})
$$
where $Tr$ denotes the trace of the tensor and $L$ is an additional
property of the viscoelastic fluid.

We will use the Bell--Collela--Glaz scheme to advect the log-conformation
tensor $\Psi$. */

/*
TODO:
- Perhaps, instead of the Bell--Collela--Glaz scheme, we can use the conservative form of the advection equation and transport the log-conformation tensor with the VoF color function, similar to [http://basilisk.fr/src/navier-stokes/conserving.h](http://basilisk.fr/src/navier-stokes/conserving.h)
*/

#include "bcg.h"

(const) scalar Gp = unity; // elastic modulus
(const) scalar lambda = unity; // relaxation time

/**
## Finite extensibility: FENE-P

Oldroyd-B is a Hookean dumbbell: a coil stretches without bound and the
extensional viscosity grows without limit. FENE-P replaces the spring with a
finitely extensible one through the Peterlin closure. We use the
equilibrium-normalised form standard in the filament-thinning literature,

$$
f(\mathrm{tr}\,\mathbf{A}) = \frac{L^2 - d}{L^2 - \mathrm{tr}\,\mathbf{A}},
\qquad \mathbf{T} = G_p\,(f\mathbf{A} - \mathbf{I}),
\qquad \left.\partial_t\mathbf{A}\right|_{relax}
      = -\frac{f\mathbf{A} - \mathbf{I}}{\lambda}
$$

`L2` is the squared ratio of the fully extended dumbbell length to its
equilibrium r.m.s. extension, so `tr(A)` runs from `d` at equilibrium to `L2`
at full extension.

**Why this normalisation.** `f(d) = 1`, so the equilibrium is `A = I` and
`T = 0` exactly as for Oldroyd-B. Every initial condition, boundary condition,
gas reset, adaptation tolerance and post-processing script written for
Oldroyd-B keeps its meaning, and the zero-shear polymeric viscosity is
`eta_p = Gp*lambda` exactly, so `Ec` and `De` are unchanged. `L2 = HUGE` (the
default) gives `f == 1` identically and recovers Oldroyd-B bit for bit.

The other common convention (Bird et al. 1987, and Basilisk's stock
`fene-p.h`) writes `f = 1/(1 - tr(A)/b)` with equilibrium `A = I*b/(b+d)`.
It is the same model: `A_here = A_Bird*(b+d)/b`, `L2 = b + d`,
`lambda_here = lambda_Bird*b/(b+d)`, same `G`. Convert before comparing
against anything built on the stock header; above `L2 ~ 1e3` the two agree to
better than a part in 300.

**This is FENE-P, not FENE-CR.** FENE-CR is `T = Gp*f*(A - I)` and has a
constant shear viscosity. The difference is the `I` inside versus outside the
`f`, and it is a different fluid.

`d` is the number of degrees of freedom of the dumbbell. Axisymmetric runs
carry the real hoop component `AThTh`, so the trace is the true
three-dimensional one and `d = 3`. A planar run carries only two components;
`d = 2` there describes a two-dimensional dumbbell, which is what Basilisk's
stock header also does in planar flow. Production is axisymmetric; the planar
case exists for the Poiseuille cross-check.
*/

#if AXI
# define FENEP_NDOF 3.
#else
# define FENEP_NDOF 2.
#endif

double L2 = HUGE;   // finite extensibility; HUGE = Oldroyd-B

static inline int fenep_active (void) {
  /**
  `HUGE` is `1e30f`, which is finite, so `isfinite()` alone does not exclude
  the Oldroyd-B sentinel. It happens that `f` then evaluates to exactly 1 in
  double precision and the answer is unchanged, but that is a coincidence of
  rounding, not a design: the test below keeps the Hookean path genuinely
  free of the FENE-P solve. */
  return isfinite (L2) && L2 < HUGE;
}

/**
`f` evaluated defensively: the argument is clamped just below `L2` so that a
diagnostic reading `A` between steps -- after adaptation, prolongation or a
restart -- can never produce a non-positive or infinite `f`. The relaxation
below does not rely on this clamp; it guarantees `tr(A) < L2` structurally. */

static inline double fenep_f (double s) {
  if (!fenep_active())
    return 1.;
  double smax = L2*(1. - 1e-9);
  if (!(s < smax)) s = smax;
  return (L2 - FENEP_NDOF)/(L2 - s);
}

/**
### The relaxation trace, solved implicitly in `f`

Oldroyd-B relaxes by an exact integral because `dA/dt = -(A - I)/lambda` is
linear. FENE-P is not: `f` depends on `tr(A)`, which is what is relaxing. We
freeze `f` over the step but choose it *self-consistently* with the end-of-step
trace, which keeps the update exact in `A` for that `f` and, crucially, keeps
`tr(A) < L2` for any input -- including an input that the closure-blind stretch
substep has already pushed past `L2`.

Solve for `s1`:
$$
s_1 = \frac{d}{f(s_1)}
    + \left(s_0 - \frac{d}{f(s_1)}\right) e^{-f(s_1) h}, \qquad h = \Delta t/\lambda .
$$
The residual is positive at `s = 0` and tends to `-L2` as `s -> L2`, so a root
is always bracketed by `[0, L2)`. Fixed-point iteration converges in two or
three steps at our `h ~ 1e-4`; the bracket is maintained so that a failure
falls back to bisection rather than to a wrong answer.

Freezing `f` at the *old* trace instead -- what the stock header does -- is
cheaper and stable while `tr(A) < L2`, but if the stretch substep overshoots
then `f(s0) <= 0`, `A` loses positive-definiteness and the stress becomes a
large compressive force in the momentum equation. That is the failure mode
case 2330 showed with an unbounded hoop source. Build with
`-DFENEP_EXPLICIT_F` to reproduce it deliberately in a test; never in
production. */

static double fenep_relax_trace (double s0, double h)
{
  const double d = FENEP_NDOF;
  double lo = 0., hi = L2*(1. - 1e-12);
  double s = (s0 < hi ? s0 : hi);
  if (!(s > 0.)) s = 0.;
  for (int it = 0; it < 60; it++) {
    double f = (L2 - d)/(L2 - s);
    double Aeq = d/f;
    double snew = Aeq + (s0 - Aeq)*exp (-f*h);
    double g = snew - s;                 // g(lo) > 0, g(hi) < 0
    if (g > 0.) lo = s; else hi = s;
    if (fabs (g) <= 1e-12*(1. + fabs (s)))
      return snew < hi ? snew : hi;
    s = (snew > lo && snew < hi) ? snew : 0.5*(lo + hi);
  }
  return 0.5*(lo + hi);
}

scalar A11[], A12[], A22[]; // conformation tensor
scalar T11[], T12[], T22[]; // stress tensor
#if AXI
scalar AThTh[], T_ThTh[];
#endif

event defaults (i = 0) {
  if (is_constant (a.x))
    a = new face vector;

  /*
  initialize A and T
  */
  for (scalar s in {A11, A22}) {
    foreach () {
      s[] = 1.;
    }
  }
  for (scalar s in {T11, T12, T22, A12}) {
    foreach(){
      s[] = 0.;
    }
  }
#if AXI
  foreach(){
    T_ThTh[] = 0;
    AThTh[] = 1.;
  }
#endif

  for (scalar s in {T11, T12, T22}) {
    if (s.boundary[left] != periodic_bc) {
        s[left] = neumann(0);
	      s[right] = neumann(0);
      }
  }

  for (scalar s in {A11, A12, A22}) {
    if (s.boundary[left] != periodic_bc) {
        s[left] = neumann(0);
	      s[right] = neumann(0);
    }
  }

#if AXI
  T12[bottom] = dirichlet (0.);
  A12[bottom] = dirichlet (0.);
#endif
}

/**
## Useful functions in 2D

The first step is to implement a routine to calculate the eigenvalues
and eigenvectors of the conformation tensor $\mathbf{A}$.

These structs ressemble Basilisk vectors and tensors but are just
arrays not related to the grid. */

typedef struct { double x, y;}   pseudo_v;
typedef struct { pseudo_v x, y;} pseudo_t;

// Function to initialize pseudo_v
static inline void init_pseudo_v(pseudo_v *v, double value) {
    v->x = value;
    v->y = value;
}

// Function to initialize pseudo_t
static inline void init_pseudo_t(pseudo_t *t, double value) {
    init_pseudo_v(&t->x, value);
    init_pseudo_v(&t->y, value);
}

static void diagonalization_2D (pseudo_v * Lambda, pseudo_t * R, pseudo_t * A)
{
  /**
  The eigenvalues are saved in vector $\Lambda$ computed from the
  trace and the determinant of the symmetric conformation tensor
  $\mathbf{A}$. */

  if (sq(A->x.y) < 1e-15) {
    R->x.x = R->y.y = 1.;
    R->y.x = R->x.y = 0.;
    Lambda->x = A->x.x; Lambda->y = A->y.y;
    return;
  }

  double T = A->x.x + A->y.y; // Trace of the tensor
  double D = A->x.x*A->y.y - sq(A->x.y); // Determinant

  /**
  The eigenvectors, $\mathbf{v}_i$ are saved by columns in tensor
  $\mathbf{R} = (\mathbf{v}_1|\mathbf{v}_2)$. */

  R->x.x = R->x.y = A->x.y;
  R->y.x = R->y.y = -A->x.x;
  double s = 1.;
  for (int i = 0; i < dimension; i++) {
    double * ev = (double *) Lambda;
    ev[i] = T/2 + s*sqrt(sq(T)/4. - D);
    s *= -1;
    double * Rx = (double *) &R->x;
    double * Ry = (double *) &R->y;
    Ry[i] += ev[i];
    double mod = sqrt(sq(Rx[i]) + sq(Ry[i]));
    Rx[i] /= mod;
    Ry[i] /= mod;
  }
}

/**
The stress tensor depends on previous instants and has to be
integrated in time. In the log-conformation scheme the advection of
the stress tensor is circumvented, instead the conformation tensor,
$\mathbf{A}$ (or more precisely the related variable $\Psi$) is
advanced in time.

In what follows we will adopt a scheme similar to that of [Hao \& Pan
(2007)](#hao2007). We use a split scheme, solving successively

a) the upper convective term:
$$
\partial_t \Psi = 2 \mathbf{B} + (\Omega \cdot \Psi -\Psi \cdot \Omega)
$$
b) the advection term:
$$
\partial_t \Psi + \nabla \cdot (\Psi \mathbf{u}) = 0
$$
c) the model term (but set in terms of the conformation
tensor $\mathbf{A}$). In an Oldroyd-B viscoelastic fluid, the model is
$$
\partial_t \mathbf{A} = -\frac{\mathbf{f}_r (\mathbf{A})}{\lambda}
$$
*/

event tracer_advection(i++)
{
  scalar Psi11 = A11;
  scalar Psi12 = A12;
  scalar Psi22 = A22;
#if AXI
  scalar Psiqq = AThTh;
#endif

  /**
  ### Computation of $\Psi = \log \mathbf{A}$ and upper convective term */

  foreach() {
    /**
      We assume that the stress tensor $\mathbf{\tau}_p$ depends on the
      conformation tensor $\mathbf{A}$ as follows
      $$
      \mathbf{\tau}_p = G_p (\mathbf{A}) =
      G_p (\mathbf{A} - I)
      $$
    */

    pseudo_t A;

    A.x.x = A11[]; A.y.y = A22[];
    A.x.y = A12[];

#if AXI
    double Aqq = AThTh[];
    if (Aqq <= 0. || !isfinite(Aqq)) {
      fprintf(ferr, "Invalid axisymmetric conformation: Aqq = %g\n", Aqq);
      fprintf(ferr, "x = %g, y = %g\n", x, y);
#if _MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
#endif
      exit(1);
    }
    Psiqq[] = log (Aqq);
#endif

    /**
    The conformation tensor is diagonalized through the
    eigenvector tensor $\mathbf{R}$ and the eigenvalues diagonal
    tensor, $\Lambda$. */

    pseudo_v Lambda;
    init_pseudo_v(&Lambda, 0.0);
    pseudo_t R;
    init_pseudo_t(&R, 0.0);
    diagonalization_2D (&Lambda, &R, &A);

    /*
    Check for negative eigenvalues -- this should never happen. If it does, print the location and value of the offending eigenvalue.
    Please report this bug by opening an issue on the GitHub repository.
    */
    if (Lambda.x <= 0. || Lambda.y <= 0.) {
      fprintf(ferr, "Negative eigenvalue detected: Lambda.x = %g, Lambda.y = %g\n", Lambda.x, Lambda.y);
      fprintf(ferr, "x = %g, y = %g\n", x, y);
#if _MPI
      MPI_Abort(MPI_COMM_WORLD, 1);
#endif
      exit(1);
    }

    /**
    $\Psi = \log \mathbf{A}$ is easily obtained after diagonalization,
    $\Psi = R \cdot \log(\Lambda) \cdot R^T$. */

    Psi12[] = R.x.x*R.y.x*log(Lambda.x) + R.y.y*R.x.y*log(Lambda.y);
    Psi11[] = sq(R.x.x)*log(Lambda.x) + sq(R.x.y)*log(Lambda.y);
    Psi22[] = sq(R.y.y)*log(Lambda.y) + sq(R.y.x)*log(Lambda.x);

    /**
    We now compute the upper convective term $2 \mathbf{B} +
    (\Omega \cdot \Psi -\Psi \cdot \Omega)$.

    The diagonalization will be applied to the velocity gradient
    $(\nabla u)^T$ to obtain the antisymmetric tensor $\Omega$ and
    the traceless, symmetric tensor, $\mathbf{B}$. If the conformation
    tensor is $\mathbf{I}$, $\Omega = 0$ and $\mathbf{B}= \mathbf{D}$.

    Otherwise, compute M = R * (nablaU)^T * R^T, where nablaU is the velocity gradient tensor. Then,

    1. Calculate omega using the off-diagonal elements of M and eigenvalues:
       omega = (Lambda.y*M.x.y + Lambda.x*M.y.x)/(Lambda.y - Lambda.x)
       This represents the rotation rate in the eigenvector basis.

    2. Transform omega back to physical space to get OM:
       OM = (R.x.x*R.y.y - R.x.y*R.y.x)*omega
       This gives us the rotation tensor Omega in the original coordinate system.

    3. Compute B tensor components using M and R: B is related to M and R through:

       In 2D:
       $$
       B_{xx} = R_{xx}^2 M_{xx} + R_{xy}^2 M_{yy} \\
       B_{xy} = R_{xx}R_{yx} M_{xx} + R_{xy}R_{yy} M_{yy} \\
       B_{yx} = B_{xy} \\
       B_{yy} = -B_{xx}
       $$

       Where:
       - R is the eigenvector matrix of the conformation tensor
       - M is the velocity gradient tensor in the eigenvector basis
       - The construction ensures B is symmetric and traceless
    */

    pseudo_t B;
    init_pseudo_t(&B, 0.0);
    double OM = 0.;
    if (fabs(Lambda.x - Lambda.y) <= 1e-20) {
      B.x.y = (u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(4.*Delta);
      foreach_dimension()
        B.x.x = (u.x[1,0] - u.x[-1,0])/(2.*Delta);
    } else {
      pseudo_t M;
      init_pseudo_t(&M, 0.0);
      foreach_dimension() {
        M.x.x = (sq(R.x.x)*(u.x[1] - u.x[-1]) +
        sq(R.y.x)*(u.y[0,1] - u.y[0,-1]) +
        R.x.x*R.y.x*(u.x[0,1] - u.x[0,-1] +
        u.y[1] - u.y[-1]))/(2.*Delta);

        M.x.y = (R.x.x*R.x.y*(u.x[1] - u.x[-1]) +
        R.x.y*R.y.x*(u.y[1] - u.y[-1]) +
        R.x.x*R.y.y*(u.x[0,1] - u.x[0,-1]) +
        R.y.x*R.y.y*(u.y[0,1] - u.y[0,-1]))/(2.*Delta);
      }
      double omega = (Lambda.y*M.x.y + Lambda.x*M.y.x)/(Lambda.y - Lambda.x);
      OM = (R.x.x*R.y.y - R.x.y*R.y.x)*omega;

      B.x.y = M.x.x*R.x.x*R.y.x + M.y.y*R.y.y*R.x.y;
      foreach_dimension()
        B.x.x = M.x.x*sq(R.x.x)+M.y.y*sq(R.x.y);
    }

    /**
    We now advance $\Psi$ in time, adding the upper convective
    contribution. */

    double s = -Psi12[];
    Psi12[] += dt * (2. * B.x.y + OM * (Psi22[] - Psi11[]));
    s *= -1;
    Psi11[] += dt * 2. * (B.x.x + s * OM);
    s *= -1;
    Psi22[] += dt * 2. * (B.y.y + s * OM);

    /**
    In the axisymmetric case, the governing equation for $\Psi_{\theta
    \theta}$ only involves that component,
    $$
    \Psi_{\theta \theta}|_t - 2 L_{\theta \theta} =
    \frac{\mathbf{f}_r(e^{-\Psi_{\theta \theta}})}{\lambda}
    $$
    with $L_{\theta \theta} = u_y/y$. Therefore step (a) for
    $\Psi_{\theta \theta}$ is */

#if AXI
    Psiqq[] += dt*2.*u.y[]/max(y, 1e-20);
#endif

}

  /**
  ### Advection of $\Psi$

  We proceed with step (b), the advection of the log of the
  conformation tensor $\Psi$. */

#if AXI
  advection ({Psi11, Psi12, Psi22, Psiqq}, uf, dt);
#else
  advection ({Psi11, Psi12, Psi22}, uf, dt);
#endif

  /**
  ### Convert back to Aij */

  foreach() {
    /**
    It is time to undo the log-conformation, again by
    diagonalization, to recover the conformation tensor $\mathbf{A}$
    and to perform step (c).*/

    pseudo_t A = {{Psi11[], Psi12[]}, {Psi12[], Psi22[]}}, R;
    init_pseudo_t(&R, 0.0);
    pseudo_v Lambda;
    init_pseudo_v(&Lambda, 0.0);
    diagonalization_2D (&Lambda, &R, &A);
    Lambda.x = exp(Lambda.x), Lambda.y = exp(Lambda.y);

    A.x.y = R.x.x*R.y.x*Lambda.x + R.y.y*R.x.y*Lambda.y;
    foreach_dimension()
      A.x.x = sq(R.x.x)*Lambda.x + sq(R.x.y)*Lambda.y;
#if AXI
      double Aqq = exp(Psiqq[]);
#endif

    /**
    We perform now step (c) by integrating
    $\mathbf{A}_t = -\mathbf{f}_r (\mathbf{A})/\lambda$ to obtain
    $\mathbf{A}^{n+1}$. This step is analytic,
    $$
    \int_{t^n}^{t^{n+1}}\frac{d \mathbf{A}}{\mathbf{I}- \mathbf{A}} =
    \frac{\Delta t}{\lambda}
    $$
    */

    double intFactor, Aeq = 1.;
    if (lambda[] == 0.)
      intFactor = 0.;              // no polymer here: reset to equilibrium
    else if (lambda[] == 1e30)
      intFactor = 1.;              // infinite Deborah: no relaxation at all
    else if (!fenep_active())
      intFactor = exp(-dt/lambda[]);            // Oldroyd-B, exact
    else {
      double s0 = A.x.x + A.y.y;
#if AXI
      s0 += Aqq;
#endif
      double h = dt/lambda[];
#ifdef FENEP_EXPLICIT_F
      double fstar = fenep_f (s0);              // test-only: stock's freeze
#else
      double fstar = (L2 - FENEP_NDOF)/(L2 - fenep_relax_trace (s0, h));
#endif
      intFactor = exp(-fstar*h);
      Aeq = 1./fstar;
    }

#if AXI
      Aqq = Aeq*(1. - intFactor) + intFactor*Aqq;
#endif

    A.x.y *= intFactor;
    foreach_dimension()
      A.x.x = Aeq*(1. - intFactor) + A.x.x*intFactor;

    /**
      Then the Conformation tensor $\mathcal{A}_p^{n+1}$ is restored from
      $\mathbf{A}^{n+1}$.  */

    /**
    The stress uses `f` at the end-of-step trace. For Oldroyd-B `fs == 1` and
    this is the original `T = Gp*(A - I)`. */

    double sEnd = A.x.x + A.y.y;
#if AXI
      sEnd += Aqq;
#endif
    double fs = fenep_f (sEnd);

    A12[] = A.x.y;
    T12[] = Gp[]*fs*A.x.y;
#if AXI
      AThTh[] = Aqq;
      T_ThTh[] = Gp[]*(fs*Aqq - 1.);
#endif

    A11[] = A.x.x;
    T11[] = Gp[]*(fs*A.x.x - 1.);
    A22[] = A.y.y;
    T22[] = Gp[]*(fs*A.y.y - 1.);
  }
}

/**
### Divergence of the viscoelastic stress tensor

The viscoelastic stress tensor $\mathbf{\tau}_p$ is defined at cell centers
while the corresponding force (acceleration) will be defined at cell
faces. Two terms contribute to each component of the momentum
equation. For example the $x$-component in Cartesian coordinates has
the following terms: $\partial_x \mathbf{\tau}_{p_{xx}} + \partial_y
\mathbf{\tau}_{p_{xy}}$. The first term is easy to compute since it can be
calculated directly from center values of cells sharing the face. The
other one is harder. It will be computed from vertex values. The
vertex values are obtained by averaging centered values.  Note that as
a result of the vertex averaging cells `[]` and `[-1,0]` are not
involved in the computation of shear. */

event acceleration (i++)
{
  face vector av = a;

  foreach_face(x){
    if (fm.x[] > 1e-20) {

      double shearX = (T12[0,1]*cm[0,1] + T12[-1,1]*cm[-1,1] -
      T12[0,-1]*cm[0,-1] - T12[-1,-1]*cm[-1,-1])/4.;

      av.x[] += (shearX + cm[]*T11[] - cm[-1]*T11[-1])*
      alpha.x[]/(sq(fm.x[])*Delta);

    }
  }

  foreach_face(y){
    if (fm.y[] > 1e-20) {

      double shearY = (T12[1,0]*cm[1,0] + T12[1,-1]*cm[1,-1] -
      T12[-1,0]*cm[-1,0] - T12[-1,-1]*cm[-1,-1])/4.;

      av.y[] += (shearY + cm[]*T22[] - cm[0,-1]*T22[0,-1])*
      alpha.y[]/(sq(fm.y[])*Delta);

    }
  }

#if AXI
  foreach_face(y)
    if (y > 1e-20)
      av.y[] -= (T_ThTh[] + T_ThTh[0,-1])*alpha.y[]/sq(y)/2.;
#endif
}
