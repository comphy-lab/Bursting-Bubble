/**
# Bursting Bubbles in Newtonian Fluids — DRILL adaptive resolution

This is `burstingBubble-adaptiveResolution.c` (the jet-base probe) with the
probe promoted from a passive diagnostic into the *driver* of a
feature-tracking adaptive-mesh + adaptive-time trigger. It "drills" resolution
into the collapsing region as the singularity approaches and lets it relax
elsewhere. The refinement philosophy is taken from
`comphy-lab/ElasticPinchOff` (`LiquidOutThinning.c`), generalised from a
neck-radius scalar to the bursting-bubble probe.

## What changed relative to the logging-only reference

1. The mesh ceiling is no longer the fixed `MAXlevel`. A *local* ceiling
   `maxlevelLocal` tracks the resolved length of the active feature and is
   passed to `adapt_wavelet`. `MAXlevel` becomes the hard cap; the far field
   still coarsens to `MINlevel` through the wavelet criterion.
2. Event layout mirrors the reference's ordering `adapt → probe`: the
   `adapt(i++)` event is the reference adapt with `maxlevelLocal`; a separate
   `drillProbe(i++)` event (declared after) runs the probe and sets
   `maxlevelLocal` for the NEXT step's adapt (a one-step lag, negligible).
   Keeping the probe out of the adapt event (not between `curvature()` and
   `adapt_wavelet()`) matches the proven reference structure.
3. Snapshot cadence (`tsnap`) is tightened as the mesh refines, so the fast
   jet/pinch window is sampled densely — adaptive time-resolution of output.
4. `logWriting` reuses the probe/level state via the `g_*` globals.

## The drill (ElasticPinchOff-style, two regimes)

Tracked length `s`:
- pre-inception (cavity-focus collapse): `s = 1/|kappa|_max` — the local radius
  of curvature at the focusing cavity base. Diverging curvature => shrinking `s`.
- post-inception (jet growth): `s = r_b` — the jet-base radius (`rlow`). On the
  slender/constant-flux branch this keeps shrinking as the jet lengthens, so
  resolution demand keeps climbing through the whole jet phase.

Target level = smallest `L in [drillMaxlevelStart, MAXlevel]` such that `s`
spans at least `nCells` cells: `s >= nCells * Ldomain / 2^L`. `nCells` differs
by regime (`drillNcellsK`, `drillNcellsJet`) because a curving front and a
slender column are geometrically different features. Hysteresis: refine
immediately to the demanded level (never under-resolve a growing feature),
coarsen at most one level per step (kills refine/coarsen churn). Because
`adapt_wavelet` itself moves the grid by at most one level per call, the actual
grid ramps smoothly regardless.

Terminal relaxation (OFF by default): once the first tip droplet sheds
(`jetFormed && n_components > 1`), if `drillRelaxLevel > 0` the ceiling relaxes
to that level to run the long post-pinch tail cheaply. Disabled by default
because the target observable R_j x Q_L is measured at the receding base, which
keeps demanding resolution — see the profile discussion.

Why not just port ElasticPinchOff verbatim: pinch-off tracks ONE scalar whose
global min IS the neck (`statsf(Y).min` is safe by construction), has ONE
singular event, and never resets `broken`. Bursting-bubble has multiple
interface features (flat free surface, wall meniscus, cavity wall, shed drop),
so it needs the MainPhase-tagged, near-axis-restricted probe; it has (at least)
two singular events in sequence governed by DIFFERENT scalars; and the jet base
keeps demanding resolution after inception rather than stopping.

## Adaptive time

The solver timestep is already capillary-wave-limited by `tension.h`
(`dt ~ sqrt(rho_m * Delta_min^3 / (pi*sigma))`), so as `maxlevelLocal` ramps and
`Delta_min` shrinks the timestep tightens automatically — the mesh drill drives
the time drill for free. The only manual time control is the staged snapshot
cadence above.

## Probe (unchanged from the reference)

Over the MAIN connected liquid body (detached drops excluded), interfacial
cells f in (1e-6, 1-1e-6) with y < RCAV:
  - (z_low, r_low) : globally lowest interfacial point (min axial x).
  - (z_maxk, r_maxk): point of max |curvature| with x < ZSURF_CURV.
Computed MPI-SAFELY (tag() + MPI_Allreduce for MainPhase; cross-rank reductions
for the argmin/argmax value then coordinates). Inception latch (never resets):
jetFormed = 1 when rmaxk in [0, R_AXIS_K) and rlow > AXIS_BAND. Sentinel -1000
when no probe exists.

Set `drillAMR=0` to pin `maxlevelLocal = MAXlevel` and reproduce the
fixed-level reference run bit-for-bit (for A/B validation).

## Status (Jul 2026)

Validated serial, OpenMP, and MPI: the full regime sweep fires end-to-end and
pre-focus kinetic energy matches the fixed-level-12 reference to ~1 %.

MPI note: build the MPI binary WITHOUT -D_GNU_SOURCE. That flag enables
Basilisk's FP trap on Linux, which fires SPURIOUSLY on its own `undefined`
NaN-sentinel in a transient ghost cell left by the drill's coarsen/refine/
rebalance and aborts an otherwise-correct run with SIGFPE. With the trap off,
np=2 and np=4 give bit-identical, physically-correct results (the ke
blow-up/decay checks remain the real-divergence guard). See DRILL-RESOLUTION.md.

Coords: x = axial (= z), y = radial (= r >= 0). Dump carries only f and u.

@file burstingBubble-drillResolution.c
@author Vatsal Sanjay (vatsal.sanjay@comphy-lab.org) / CoMPhy Lab
@version 2.0
@date Jan 04, 2025 (jet-base logging Jun 2026; drill AMR/dt trigger Jul 2026)
*/

#include "axi.h"
#include "navier-stokes/centered.h"

/**
## Solver Configuration

- `FILTERED`: Enable density and viscosity jump smoothing
*/
#define FILTERED 1 // Smear density and viscosity jumps
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"

#if !_MPI
#include "distance.h"
#endif

/**
## Extra includes for jet-base diagnostics

`tag.h` provides connected-component tagging (main liquid body vs detached
drops); `curvature.h` provides `curvature()`; `fractions.h` provides the VOF
fraction helpers. These are logging-only and do not touch the solver state.
(`two-phase.h`/`tension.h` already pull `curvature.h`/`fractions.h` in via
include guards; the explicit includes here just make the dependency obvious.)
*/
#include "fractions.h"
#include "curvature.h"
#include "tag.h"

/**
## Runtime Parameters

All configuration is read at runtime from `case.params` via the C-side
parameter layer. Adaptive space (levels, wavelet tolerances) and adaptive
time (CFL, dtmax ceiling, solver tolerance) are now tunable knobs rather
than compile-time constants — see `src-local/params.h`.
*/
#include "params.h"

struct SimulationParams params;

// Boundary conditions - outflow on the right boundary
u.n[right] = neumann(0.);
p[right] = dirichlet(0.);

// Boundary conditions - solid wall at left boundary (bottom)
f[left] = dirichlet(1.0);      // Liquid at wall
u.n[left] = dirichlet(0.0);    // No-slip normal
u.t[left] = dirichlet(0.0);    // No-slip tangential

// Mesh control (set from params in main)
int MAXlevel, MINlevel;

// Drill trigger state:
//   maxlevelLocal -> the LOCAL, time-varying refinement ceiling handed to
//     adapt_wavelet each step (<= MAXlevel). Seeded in init(), ramped in
//     drillProbe(). This is what makes the run "drill" resolution into the
//     collapsing region and relax it elsewhere.
//   drill* mirrors of the params knobs, populated in main() for terse use.
int maxlevelLocal;
int drillAMR, drillMaxlevelStart, drillRelaxLevel, drillTsnapStages;
double drillNcellsK, drillNcellsJet, drillTsnapMinFactor;

// Probe state exported by drillProbe(i++) and consumed by logWriting(i++):
//   g_rb, g_zb  -> selected jet-base / cavity-focus probe coordinates
//   g_jetFormed -> inception latch (0 before jet, 1 after)
double g_rb = -1000., g_zb = -1000.;
int    g_jetFormed = 0;

// Science observable exported by drillProbe(i++) (inlined getBase.c logic):
//   g_rbase, g_zbase -> lowest OUTER-free-surface cell of the main liquid
//     body (latch-free, satellite-proof; this is the r_b of R_j x Q_L)
//   g_qjet = INT_0^{r_base} u_z r dr, g_ql = INT_0^{r_base} u_z dr at z_base
double g_rbase = -1000., g_zbase = -1000.;
double g_qjet = -1000., g_ql = -1000.;

// Physical parameters (set from params in main):
//   Oh  -> Ohnesorge number (liquid)
//   Oha -> Ohnesorge number (gas) = OhRatio * Oh
double Oh, Oha, Bond, tmax;

// Domain parameters:
//   zWall   -> distance from bubble south pole to bottom wall
//   Ldomain -> computed domain size: min(zWall + 6.0, 16.0)
double zWall, Ldomain;

// Adaptive-resolution controls (set from params in main):
//   fErr/VelErr/KErr -> wavelet error tolerances (VOF, velocity, curvature)
double fErr, VelErr, KErr;

// tsnap -> snapshot/restart dump interval. Needs a non-zero static initial
// value: Basilisk classifies event expressions (e.g. `t += tsnap`) before
// main() runs, and a zero increment would be misread as a second condition.
// main() overrides this with params.tsnap for the actual firing interval.
double tsnap = 1e-2;

char nameOut[80], dumpFile[80];

/**
## Jet-base probe constants (match getJetFoot.c exactly)

Geometry-tuned for case 1000: origin(-6,0), L0=10, free surface near z=0.
*/
#define RCAV       1.20    // exclude the flat outer free surface (r spans to L0)
#define ZSURF_CURV 0.0     // max|k| search restricted below this axial level
#define R_AXIS_K   0.05    // inception latch: max|k| point must be near the axis
#define AXIS_BAND  0.04    // inception latch: lowest point must be off the axis

/**
## Main Function

Reads parameters, configures the domain and fluid properties, and starts
the run.

- Parses `case.params` (or legacy positional CLI) into `params`
- Validates the configuration before allocating the grid
- Sets up the physical domain with appropriate dimensions
- Configures fluid properties for both phases
- Maps the adaptive space/time knobs onto Basilisk's solver globals
*/
int main(int argc, char *argv[]) {
  // Parse and validate runtime configuration
  if (params_init(argc, argv, &params) != 0)
    return 1;
  if (!validate_params(&params)) {
    fprintf(ferr, "ERROR: Invalid parameters. Aborting.\n");
    return 1;
  }

  // Map physical parameters onto module globals
  MAXlevel = params.MAXlevel;
  MINlevel = params.MINlevel;
  Oh = params.Oh;
  Oha = params.OhRatio * params.Oh;
  Bond = params.Bond;
  tmax = params.tmax;
  zWall = params.zWall;
  fErr = params.fErr;
  VelErr = params.VelErr;
  KErr = params.KErr;
  tsnap = params.tsnap;

  // Drill trigger knobs
  drillAMR            = params.drillAMR;
  drillMaxlevelStart  = params.drillMaxlevelStart;
  drillNcellsK        = params.drillNcellsK;
  drillNcellsJet      = params.drillNcellsJet;
  drillRelaxLevel     = params.drillRelaxLevel;
  drillTsnapStages    = params.drillTsnapStages;
  drillTsnapMinFactor = params.drillTsnapMinFactor;

  // Calculate domain size: Ldomain = min(zWall + 6.0, 16.0)
  // zWall = distance from bubble south pole to bottom wall
  // +2.0 buffer below bubble, +4.0 space above for jet
  Ldomain = fmin(zWall + 6.0, 16.0);

  L0 = Ldomain;
  origin(-2.0 - zWall, 0.);

  init_grid(1 << params.init_grid_level);

  /**
  ## Adaptive Time Control

  Set the advective CFL and the timestep ceiling. Surface tension is
  time-explicit, so `tension.h` reduces the timestep each step to the
  capillary-wave period `T = sqrt(rho_m * Delta_min^3 / (pi * sigma))`;
  `dtmax` is therefore a safety ceiling and the effective step is adaptive
  (it scales with the finest cell size and the resolved physics).
  */
  CFL = params.CFL;
  dtmax = params.dtmax;
  TOLERANCE = params.TOLERANCE;

  // Create a folder named intermediate where all the simulation snapshots are stored.
  char comm[80];
  sprintf(comm, "mkdir -p intermediate");
  system(comm);

  // Name of the restart file. See writingFiles event.
  sprintf(dumpFile, "restart");

  /**
  ## Physical Properties Configuration

  Sets up the material properties for both phases:
  - `rho1`, `rho2`: Density of liquid and gas phases
  - `mu1`, `mu2`: Dynamic viscosity of liquid and gas phases
  */
  rho1 = 1., rho2 = 1e-3;
  mu1 = Oh, mu2 = Oha;

  f.sigma = 1.0;

  if (pid() == 0)
    print_params(&params, ferr);

  run();
}

/**
## Initialization Event

Sets up the initial conditions for the simulation.

The function attempts to restore from a dump file first. If that fails:
- For MPI runs: Ends with an error
- For non-MPI runs: Tries to load an initial shape from a data file,
  creates a distance field, and initializes the volume fraction
*/
event init(t = 0) {
#if _MPI // This is for supercomputers without OpenMP support
  if (!restore(file = dumpFile)) {
    fprintf(ferr, "Cannot restore from dump file '%s': MPI runs must start from a restart dump (distance.h init is MPI-incompatible). Aborting.\n", dumpFile);
    exit(1);
  }
#else  // Note that distance.h is incompatible with OpenMPI. So, the below code should not be used with MPI
  if (!restore(file = dumpFile)) {
    char filename[60];
    sprintf(filename, "DataFiles/Bo%5.4f.dat", Bond);
    FILE *fp = fopen(filename, "rb");
    if (fp == NULL) {
      fprintf(ferr, "There is no file named %s\n", filename);
      // Try in folder one level up
      sprintf(filename, "../DataFiles/Bo%5.4f.dat", Bond);
      fp = fopen(filename, "rb");
      if (fp == NULL) {
        fprintf(ferr, "There is no file named %s\n", filename);
        return 1;
      }
    }
    coord *InitialShape;
    InitialShape = input_xy(fp);
    fclose(fp);
    scalar d[];
    distance(d, InitialShape);

    while (adapt_wavelet((scalar *){f, d}, (double[]){1e-8, 1e-8}, MAXlevel).nf);

    // The distance function is defined at the center of each cell, we have
    // to calculate the value of this function at each vertex.
    vertex scalar phi[];
    foreach_vertex() {
      phi[] = -(d[] + d[-1] + d[0,-1] + d[-1,-1])/4.;
    }

    // We can now initialize the volume fraction of the domain.
    fractions(phi, f);
  }
#endif

  /**
  Seed the local refinement ceiling. We adopt the grid's current depth so a
  restart from a genuinely-refined dump does not get coarsened on step 1;
  drillProbe then coarsens (at most one level per step) wherever the feature
  no longer demands it. On a fresh run the init loop above refined the initial
  shape to MAXlevel, so the seed is MAXlevel and the drill coarsens down toward
  drillMaxlevelStart over the first few steps. With drillAMR off, pin MAXlevel.
  */
  maxlevelLocal = drillAMR
    ? max(drillMaxlevelStart, min(MAXlevel, depth()))
    : MAXlevel;
}

/**
## Adaptive Mesh Refinement (drilled ceiling)

Byte-for-byte the reference `burstingBubble.c` adapt event EXCEPT that the
ceiling is `maxlevelLocal` (the drilled, feature-tracking level) instead of the
fixed `MAXlevel`. `maxlevelLocal` was set by the `drillProbe` event on the
PREVIOUS step (a one-step lag, negligible since the feature evolves slowly per
step); on the first step it is the seed from `init`.

Ordering: this event runs FIRST (declared before `drillProbe`), so
`adapt_wavelet` sees the clean post-solver state — exactly the position the
reference solver adapts in. The probe's `tag()` / reduction passes run AFTER
adapt, in `drillProbe`, mirroring the reference (whose probe lives in
`logWriting`, after adapt). Keeping the probe out of the adapt event matches
the proven reference structure. (A separate MPI restart FPE on the dynamic
level path is still open — see the file header Status note.)
*/
event adapt(i++) {
  scalar KAPPA[];
  curvature(f, KAPPA);

  adapt_wavelet((scalar *){f, u.x, u.y, KAPPA},
    (double[]){fErr, VelErr, VelErr, KErr},
    maxlevelLocal, MINlevel);
}

/**
## Drill probe + level decision

Runs once per step, AFTER `adapt` (declared after it). Computes curvature and
the jet-base / cavity-focus probe once, derives the demanded refinement level
from the tracked feature length, updates `maxlevelLocal` (refine-fast /
coarsen-slow hysteresis) for the NEXT step's `adapt`, and tightens the snapshot
cadence. Probe state is exported to `logWriting` via the `g_*` globals.
*/
event drillProbe(i++) {
  /**
  ### Curvature

  KAPPA is LOCAL to this event (computed, used, discarded each step), matching
  the reference solver's memory pattern rather than persisting a field full of
  the `nodata` sentinel across MPI tree restructurings. Good practice; not, on
  its own, sufficient for MPI stability — see the file header's MPI note.
  */
  scalar KAPPA[];
  curvature(f, KAPPA);

  /**
  ### Main connected liquid region (MPI-safe)

  tag() merges components across ranks; the per-region size sum is reduced so
  MainPhase is the GLOBALLY largest region. `n` is the component count — used
  below as the tip-pinch signal (n > 1 after inception => a droplet has shed).
  */
  scalar dtag[];
  foreach() dtag[] = (f[] > 1e-4);
  int n = tag(dtag);
  int MainPhase = 0;
  if (n > 0) {
    double *sz = calloc(n, sizeof(double));
    foreach(serial) if (dtag[] > 0) sz[((int) dtag[]) - 1] += 1.;  // serial: avoid OpenMP race
#if _MPI
    MPI_Allreduce(MPI_IN_PLACE, sz, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    double sm = -1.;
    for (int j = 0; j < n; j++) if (sz[j] > sm) { sm = sz[j]; MainPhase = j + 1; }
    free(sz);
  }

  /**
  ### Candidates (MPI-safe argmin/argmax, identical to the reference probe)
  */
  double zlow = HUGE, kmax = -1.;
  foreach(reduction(min:zlow) reduction(max:kmax)) {
    if (f[] <= 1e-6 || f[] >= 1. - 1e-6 || dtag[] != MainPhase || y > RCAV) continue;
    if (x < zlow) zlow = x;
    if (x < ZSURF_CURV && KAPPA[] != nodata) {
      double ak = fabs(KAPPA[]);
      if (ak > kmax) kmax = ak;
    }
  }
  double rlow = HUGE, zk = HUGE;
  foreach(reduction(min:rlow) reduction(min:zk)) {
    if (f[] <= 1e-6 || f[] >= 1. - 1e-6 || dtag[] != MainPhase || y > RCAV) continue;
    if (zlow != HUGE && x == zlow && y < rlow) rlow = y;
    if (kmax >= 0. && x < ZSURF_CURV && KAPPA[] != nodata
        && fabs(KAPPA[]) == kmax && x < zk) zk = x;
  }
  double rk = HUGE;
  foreach(reduction(min:rk)) {
    if (f[] <= 1e-6 || f[] >= 1. - 1e-6 || dtag[] != MainPhase || y > RCAV) continue;
    if (kmax >= 0. && zk != HUGE && x < ZSURF_CURV && KAPPA[] != nodata
        && fabs(KAPPA[]) == kmax && x == zk && y < rk) rk = y;
  }
  if (zlow == HUGE) { zlow = -1000.; rlow = -1000.; }
  if (rlow == HUGE) rlow = -1000.;
  if (kmax < 0. || zk == HUGE) { zk = -1000.; rk = -1000.; }
  if (rk == HUGE) rk = -1000.;

  /**
  ### Inception latch + probe selection (identical to the reference)
  */
  static int jetFormed = 0;
  if (!jetFormed && rk >= 0 && rk < R_AXIS_K && rlow > AXIS_BAND)
    jetFormed = 1;

  double zb, rb;
  if (jetFormed) { zb = zlow; rb = rlow; }   // rule 1: jet base
  else           { zb = zk;   rb = rk;   }   // rule 2: cavity focus
  if (rb <= 0.) { zb = -1000.; rb = -1000.; }

  /**
  ### The drill: demanded level from the tracked feature length

  Tracked length `s`: pre-inception it is the cavity-focus radius of curvature
  `1/kmax`; post-inception it is the jet-base radius `rb`. Target level is the
  smallest L in [drillMaxlevelStart, MAXlevel] such that `s` spans >= nCells
  cells at level L, i.e. `s >= nCells * Ldomain / 2^L`. Integer while-loop
  (no log2/pow) so it is exact and branch-clean.
  */
  int Ltarget = maxlevelLocal;
  if (drillAMR) {
    // Terminal relaxation: latch on the first tip-droplet shed after inception.
    static int tipPinched = 0;
    if (!tipPinched && jetFormed && n > 1) tipPinched = 1;

    if (tipPinched && drillRelaxLevel > 0) {
      Ltarget = drillRelaxLevel;              // run the post-pinch tail cheaply
    } else {
      double s     = jetFormed ? rb : (kmax > 0. ? 1.0 / kmax : -1.0);
      double nCell = jetFormed ? drillNcellsJet : drillNcellsK;
      int L = drillMaxlevelStart;
      if (s > 0.) {
        double need = nCell * Ldomain / s;    // 2^L must reach this
        while (L < MAXlevel && (double)(1 << L) < need) L++;
      }
      Ltarget = L;
    }

    // Hysteresis: refine immediately (never under-resolve a growing feature),
    // coarsen at most one level per step (kills refine/coarsen churn).
    if      (Ltarget > maxlevelLocal) maxlevelLocal = Ltarget;
    else if (Ltarget < maxlevelLocal) maxlevelLocal--;
  } else {
    maxlevelLocal = MAXlevel;
  }

#if DRILL_DEBUG
  if (pid() == 0) {
    double sdbg = jetFormed ? rb : (kmax > 0. ? 1.0/kmax : -1.0);
    fprintf(ferr, "DRILL i=%d jetFormed=%d n=%d kmax=%g s=%g Ltarget=%d maxlevelLocal=%d\n",
            i, jetFormed, n, kmax, sdbg, Ltarget, maxlevelLocal);
  }
#endif

  /**
  ### Adaptive time-resolution of output

  Tighten tsnap as the mesh refines so the fast jet/pinch window is densely
  sampled: tsnap = params.tsnap * 2^(start - maxlevelLocal), floored at
  drillTsnapMinFactor * params.tsnap. (The solver timestep itself is already
  adaptive via the capillary-wave limit in tension.h and tightens as Delta_min
  shrinks with maxlevelLocal — nothing to set for it here.)
  */
  if (drillAMR && drillTsnapStages) {
    double fac = ldexp(1.0, drillMaxlevelStart - maxlevelLocal); // 2^(start-level)
    if (fac < drillTsnapMinFactor) fac = drillTsnapMinFactor;
    if (fac > 1.0) fac = 1.0;                                    // never coarser than base
    tsnap = params.tsnap * fac;
  }

  /**
  ### Robust base observable + base fluxes (inlined getBase.c)

  The science observable (R_j x Q_L) wants the OUTER-free-surface base, not
  the AMR probe above: post-inception, satellite bubbles / shed droplets
  string along the axis below the base and the "globally lowest interfacial
  point" latches onto them. getBase.c's protocol is inlined here verbatim:
  MainLiq = largest PURE-liquid component (drops detached droplets), MainGas
  = largest PURE-gas component (drops entrained bubbles); an outer-surface
  cell face-touches both; the base is the lowest such cell. All MPI-safe
  (Allreduce'd tallies + cross-rank reductions), matching the validated
  serial post-processing tool. This is LOGGING ONLY — the AMR ceiling above
  keeps tracking the sharpest feature (thin jet / diverging curvature),
  which demands far more resolution than the base radius would.
  */
  scalar dl[], dg[];
  foreach() {
    dl[] = (f[] > 1. - 1e-4);
    dg[] = (f[] < 1e-4);
  }
  int nliq = tag(dl), ngas = tag(dg);
  int MainLiq = 0, MainGas = 0;
  if (nliq > 0) {
    double *sz = calloc(nliq, sizeof(double));
    foreach(serial) if (dl[] > 0) sz[(int)dl[] - 1] += 1.;   // serial: no OpenMP race
#if _MPI
    MPI_Allreduce(MPI_IN_PLACE, sz, nliq, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    double sm = -1.;
    for (int j = 0; j < nliq; j++) if (sz[j] > sm) { sm = sz[j]; MainLiq = j + 1; }
    free(sz);
  }
  if (ngas > 0) {
    double *sz = calloc(ngas, sizeof(double));
    foreach(serial) if (dg[] > 0) sz[(int)dg[] - 1] += 1.;
#if _MPI
    MPI_Allreduce(MPI_IN_PLACE, sz, ngas, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    double sm = -1.;
    for (int j = 0; j < ngas; j++) if (sz[j] > sm) { sm = sz[j]; MainGas = j + 1; }
    free(sz);
  }
  boundary((scalar *){dl, dg});   // sync tag labels into ghosts (incl. across ranks)

  // Pass 1: lowest outer-surface cell (min axial x)
  double zbase = HUGE;
  foreach(reduction(min:zbase)) {
    if (f[] <= 1e-6 || f[] >= 1. - 1e-6 || y > RCAV) continue;
    bool touchGas = ((int)dg[1,0] == MainGas) || ((int)dg[-1,0] == MainGas) ||
                    ((int)dg[0,1] == MainGas) || ((int)dg[0,-1] == MainGas);
    bool touchLiq = ((int)dl[1,0] == MainLiq) || ((int)dl[-1,0] == MainLiq) ||
                    ((int)dl[0,1] == MainLiq) || ((int)dl[0,-1] == MainLiq);
    if (touchGas && touchLiq && x < zbase) zbase = x;
  }
  // Pass 2: radius there (min-y tiebreak)
  double rbase = HUGE;
  foreach(reduction(min:rbase)) {
    if (f[] <= 1e-6 || f[] >= 1. - 1e-6 || y > RCAV) continue;
    bool touchGas = ((int)dg[1,0] == MainGas) || ((int)dg[-1,0] == MainGas) ||
                    ((int)dg[0,1] == MainGas) || ((int)dg[0,-1] == MainGas);
    bool touchLiq = ((int)dl[1,0] == MainLiq) || ((int)dl[-1,0] == MainLiq) ||
                    ((int)dl[0,1] == MainLiq) || ((int)dl[0,-1] == MainLiq);
    if (touchGas && touchLiq && zbase != HUGE && x == zbase && y < rbase) rbase = y;
  }

  // Base fluxes: single cell-row at the z_base plane, dr = Delta, y < r_base
  // (getJetFoot.c definition: no 2*pi factor, no f-weighting)
  //   q_jet = sum u_z * y * Delta  [L^3/T] ;  q_l = sum u_z * Delta  [L^2/T]
  double qjet = 0., ql = 0.;
  int haveBase = (zbase != HUGE && rbase != HUGE && rbase > 0.);
  if (haveBase) {
    foreach(reduction(+:qjet) reduction(+:ql)) {
      if (fabs(x - zbase) < 0.5*Delta && y > 0. && y < rbase) {
        qjet += u.x[] * y * Delta;
        ql   += u.x[] * Delta;
      }
    }
  }
  if (!haveBase) { zbase = -1000.; rbase = -1000.; qjet = -1000.; ql = -1000.; }

  /**
  ### Export probe/level state for logWriting

  The updated `maxlevelLocal` is consumed by the `adapt` event on the NEXT
  step. No mesh operation happens here — that is deliberate (see the `adapt`
  event header: the probe must not sit between curvature and adapt_wavelet).
  */
  g_zb = zb; g_rb = rb; g_jetFormed = jetFormed;
  g_zbase = zbase; g_rbase = rbase; g_qjet = qjet; g_ql = ql;
}

/**
## Output Management

Creates periodic snapshots of the simulation state.
- Dumps restart files for simulation recovery
- Saves intermediate snapshots at regular intervals defined by `tsnap`
*/
event writingFiles(t = 0; t += tsnap; t <= tmax) {
  dump(file = dumpFile);
  // 6 decimal places: the staged tsnap can drop below 1e-4 near inception,
  // where 4-decimal names would collide and silently overwrite snapshots.
  sprintf(nameOut, "intermediate/snapshot-%8.6f", t);
  dump(file = nameOut);
}

/**
## Simulation Termination

Writes a final summary of the simulation parameters when the simulation ends.
*/
event end(t = end) {
  if (pid() == 0)
    fprintf(ferr, "Level %d, Oh %2.1e, Oha %2.1e, Bo %4.3f, zWall %g, Ldomain %g\n",
            MAXlevel, Oh, Oha, Bond, zWall, Ldomain);
}

/**
## Simulation Logging (drill trigger + jet-base tracking)

Records key simulation data at each timestep:
- Iteration number, timestep size, current simulation time, kinetic energy
- The active local refinement ceiling `maxlevelLocal` (the drill state)
- AMR probe location r_b, z_b (the sharp-feature tracker driving the drill)
- The science observable: robust outer-surface base r_base, z_base (inlined
  getBase.c) and the base fluxes q_jet = INT u_z r dr, q_l = INT u_z dr
  through the z_base plane — the on-the-fly R_j x Q_L data.

Nothing is recomputed here — all state was computed in `drillProbe(i++)`
(which runs after `adapt`, earlier in the same step) and exported through
the `g_*` globals. Log columns:
`i dt t ke maxlevel r_b z_b r_base z_base q_jet q_l`.
Time is printed to 8 decimals and observables to 6 significant decimals so
post-hoc q_jet(r_jet) / q_l(r_jet) fits are not precision-starved.

Also performs the original safety checks (kinetic-energy blow-up / too-small),
which are left fully intact.
*/
event logWriting(i++) {
  // Calculate kinetic energy
  double ke = 0.;
  foreach(reduction(+:ke)) {
    ke += (2*pi*y)*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }

  // Probe/level/observable state exported by drillProbe (this step, post-adapt grid)
  double rb = g_rb, zb = g_zb;
  double rbase = g_rbase, zbase = g_zbase, qjet = g_qjet, ql = g_ql;

  if (pid() == 0) {
    static FILE *fp;
    if (i == 0) {
      fprintf(ferr, "MAXlevel %d, Oh %2.1e, Oha %2.1e, Bo %4.3f, zWall %g, Ldomain %g, drillAMR %d\n",
              MAXlevel, Oh, Oha, Bond, zWall, Ldomain, drillAMR);
      fprintf(ferr, "i dt t ke maxlevel r_b z_b r_base z_base q_jet q_l\n");
      fp = fopen("log", "w");
      fprintf(fp, "MAXlevel %d, Oh %2.1e, Oha %2.1e, Bo %4.3f, zWall %g, Ldomain %g, drillAMR %d\n",
              MAXlevel, Oh, Oha, Bond, zWall, Ldomain, drillAMR);
      fprintf(fp, "i dt t ke maxlevel r_b z_b r_base z_base q_jet q_l\n");
      fprintf(fp, "%d %.6e %.8f %.6e %d %.6e %.6e %.6e %.6e %.6e %.6e\n",
              i, dt, t, ke, maxlevelLocal, rb, zb, rbase, zbase, qjet, ql);
      fclose(fp);
    } else {
      fp = fopen("log", "a");
      fprintf(fp, "%d %.6e %.8f %.6e %d %.6e %.6e %.6e %.6e %.6e %.6e\n",
              i, dt, t, ke, maxlevelLocal, rb, zb, rbase, zbase, qjet, ql);
      fclose(fp);
    }
    fprintf(ferr, "%d %.6e %.8f %.6e %d %.6e %.6e %.6e %.6e %.6e %.6e\n",
            i, dt, t, ke, maxlevelLocal, rb, zb, rbase, zbase, qjet, ql);

    assert(ke > -1e-10);

    // Check for energy blowup (numerical instability)
    if (ke > 1e2 && i > 1e1) {
      fprintf(ferr, "The kinetic energy blew up. Stopping simulation\n");
      fp = fopen("log", "a");
      fprintf(fp, "The kinetic energy blew up. Stopping simulation\n");
      fclose(fp);
      dump(file = dumpFile);
      return 1;
    }
    assert(ke < 1e2);

    // Check for energy dissipation below threshold
    if (ke < 1e-6 && i > 1e1) {
      fprintf(ferr, "kinetic energy too small now! Stopping!\n");
      dump(file = dumpFile);
      fp = fopen("log", "a");
      fprintf(fp, "kinetic energy too small now! Stopping!\n");
      fclose(fp);
      return 1;
    }
  }
}
