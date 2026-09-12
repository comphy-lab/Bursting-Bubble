/**
# Bursting Bubbles in Newtonian Fluids — DRILL adaptive resolution

This is the jet-base / cavity-focus probe promoted from a passive diagnostic
into the *driver* of a
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

HIGH-MAXLEVEL WARNING (case 1005, Oh=0.029, MAXlevel=14): the cavity-focus
collapse is a genuine singularity, and chasing it with the full MAXlevel is
self-defeating — 1/|kappa|max -> 0 demands MAXlevel exactly at the singular
instant, the CFL condition then chases the diverging focusing velocity in
ever-smaller cells (dt 2e-6 -> 7e-10 while t froze at 0.46887), and the
topology change finally blew ke up 5.3 -> 1100 in two steps (sub-grid gas
wisps at reconnection). The level-12 run (case 1004) stepped straight over
the same instant. Regularisation knobs (case 1006 onward): drillMaxlevelFocus
caps the PRE-inception ramp (12 validated; full MAXlevel is released at the
inception latch, where the jet is fast but smooth), and drillRemoveGasSize
absorbs sub-grid gas wisps each step (bubbles only — liquid droplets are
physics, never touched).

POST-LATCH PINCH EVENTS (case 1009, MAXlevel=13): the released full-depth
ceiling is GLOBAL, so the entrained satellite / retracting floor remnant
BELOW the base gets resolved at full depth too — and its own pinch
singularity reproduces the 1005 pathology (t frozen at 0.47421, dt -> 7.6e-10,
ke through the relaxed gate) while the jet itself is fine. Post-inception the
ceiling is therefore REGIONAL (adapt_wavelet_limited, vendored): full
maxlevelLocal at/above z_base - DRILL_BASE_BUFFER, drillMaxlevelFocus below
it. The latch state is also persisted across restarts via the `drillstate`
file (see init); `drillAssumeJet=1` covers legacy post-inception dumps.

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

`FILTERED` is a compile-time switch from `case.params`. The runner
adds `-DFILTERED` when `FILTERED=1`. `two-phase.h` smears density
and viscosity jumps only when that macro is defined.
*/
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
`adapt_wavelet_limited.h` (vendored from the Pairetti Basilisk sandbox, see
its provenance block) provides adapt_wavelet with a POSITION-DEPENDENT max
level — the regional ceiling used post-inception (case-1009 lesson).
*/
#include "adapt_wavelet_limited.h"

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
int drillMaxlevelFocus, drillRemoveGasSize, drillMinlevelJet;
double drillNcellsK, drillNcellsJet, drillTsnapMinFactor;

// Probe state exported by drillProbe(i++) and consumed by logWriting(i++):
//   g_rb, g_zb  -> selected jet-base / cavity-focus probe coordinates
//   g_jetFormed -> inception latch (0 before jet, 1 after)
double g_rb = -1000., g_zb = -1000.;
int    g_jetFormed = 0;

// Drill latch state. File-scope (not event-local statics) so init() can
// RESTORE it across restarts — statics used to reset on restore, and a
// mid-jet restart with r_base in (AXIS_BAND, RBASE_JET) could never
// re-latch. writingFiles mirrors these three to the "drillstate" file.
int jetFormed = 0;    // inception latch (never clears during a run)
int drillArmed = 0;   // arm/fire latch: armed in the final singular approach
int baseOffAxis = 0;  // consecutive off-axis-base steps since arming
int tipPinched = 0;   // first tip-droplet-shed latch (terminal relaxation)
int tipPinchSteps = 0; // consecutive n>1 steps toward the tipPinched latch
                       // (transient, not persisted — re-counts after restart)

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

// Robust arm/fire reconnection latch (case-1006 lesson; see drillProbe):
#define ARM_BAND    0.005  // base pinned on-axis (final singular approach)
#define LATCH_STEPS 25     // consecutive off-axis steps before firing
#define RBASE_JET   0.15   // unambiguous developed-jet base radius (self-heal
                           // for restarts that jump straight into the jet phase;
                           // the pre-reconnection dimple ring peaked at ~0.11)

// Regional ceiling (case-1009 lesson; see drillMLFun/adapt):
#define DRILL_BASE_BUFFER 0.05  // full-resolution zone starts this far BELOW
                                // z_base, keeping the base + flux plane inside

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
  drillMaxlevelFocus  = params.drillMaxlevelFocus;
  drillMinlevelJet    = params.drillMinlevelJet;
  drillRemoveGasSize  = params.drillRemoveGasSize;

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
  the ceiling is therefore a safety limit and the effective step is adaptive
  (it scales with the finest cell size and the resolved physics).

  The ceiling must be written to `DT`, not to `dtmax`. `centered.h` carries
  `event set_dtmax (i++,last) dtmax = DT;`, so anything assigned to `dtmax`
  here is discarded on the first step and the knob silently does nothing.
  */
  CFL = params.CFL;
  DT = params.dtmax;
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
  bool restored = false;
#if _MPI // This is for supercomputers without OpenMP support
  if (!restore(file = dumpFile)) {
    fprintf(ferr, "Cannot restore from dump file '%s': MPI runs must start from a restart dump (distance.h init is MPI-incompatible). Aborting.\n", dumpFile);
    exit(1);
  }
  restored = true;
#else  // Note that distance.h is incompatible with OpenMPI. So, the below code should not be used with MPI
  if (restore(file = dumpFile))
    restored = true;
  else {
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
  Restore the drill latch state across restarts (case-1009 lesson): the
  latches were event-local statics and reset on every restore, so a mid-jet
  restart could never re-latch (r_base already off-axis => arming impossible)
  and the focus cap silently re-bound the whole jet. `drillstate` is written
  by writingFiles alongside each dump. Only read after a successful restore —
  a fresh run in a dirty directory must not inherit stale latches. All ranks
  read the same tiny file; no communication needed. `drillAssumeJet=1` covers
  LEGACY post-inception snapshots that predate the drillstate file (it only
  ever sets the latch, never clears it).
  */
  if (restored) {
    FILE *fs = fopen("drillstate", "r");
    if (fs) {
      int jf, da, tp;
      if (fscanf(fs, "%d %d %d", &jf, &da, &tp) == 3) {
        jetFormed = jf; drillArmed = da; tipPinched = tp;
        fprintf(ferr, "drillstate restored: jetFormed=%d drillArmed=%d tipPinched=%d\n",
                jetFormed, drillArmed, tipPinched);
      }
      fclose(fs);
    }
    if (params.drillAssumeJet && !jetFormed) {
      jetFormed = 1;
      fprintf(ferr, "drillAssumeJet=1: forcing jetFormed=1 at restore\n");
    }
  }

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
/**
### Regional ceiling function (case-1009 lesson)

Below the jet base, satellites and the retracting floor remnant undergo
REPEATED pinch singularities (the case-1009 crash: probe on-axis at
z = -1.67 vs z_base = -1.31, dt 1e-5 -> 7.6e-10 with frozen t). Those must
not be resolved at full depth — they are byproducts, not the observable.
The jet at/above the base is smooth and carries q_jet(r_jet), so it keeps
the drilled ceiling `maxlevelLocal`. The buffer keeps the base itself (and
its flux plane) inside the full-resolution zone.
*/
int drillMLFun(double x, double y, double z) {
  if (g_zbase > -900. && x < g_zbase - DRILL_BASE_BUFFER)
    // never exceed the global ceiling: after terminal relaxation
    // maxlevelLocal can drop BELOW the focus cap, and the below-base zone
    // must relax with it rather than stay pinned at drillMaxlevelFocus
    return drillMaxlevelFocus < maxlevelLocal ? drillMaxlevelFocus : maxlevelLocal;
  return maxlevelLocal;
}

event adapt(i++) {
  scalar KAPPA[];
  curvature(f, KAPPA);

  /**
  Post-inception (and only when the focus cap + a valid base exist), the
  ceiling is REGIONAL via adapt_wavelet_limited: full `maxlevelLocal` on the
  jet, `drillMaxlevelFocus` below the base. Pre-inception the plain call is
  byte-for-byte the reference adapt — the global focus cap in drillProbe
  already regularises the collapse there.
  */
  if (drillAMR && jetFormed && drillMaxlevelFocus > 0 && g_zbase > -900.)
    adapt_wavelet_limited((scalar *){f, u.x, u.y, KAPPA},
      (double[]){fErr, VelErr, VelErr, KErr},
      drillMLFun, MINlevel);
  else
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
  ### Gas-wisp cleanup (regularisation, case-1006 lesson)

  The cavity reconnection sheds sub-grid gas wisps; at high MAXlevel the CFL
  condition chases the diverging velocity inside them and dt stalls (case
  1005: dt 2e-6 -> 7e-10 at fixed t=0.46887, then ke 5.3 -> 1100 in two
  steps). Absorb gas components smaller than drillRemoveGasSize^2 cells into
  the liquid each step. bubbles=true touches ONLY gas; shed liquid droplets
  (tip droplet, satellites) are physics and are never removed. This runs
  BEFORE the probe/getBase tagging below, so the diagnostics never see the
  wisps either.
  */
  if (drillRemoveGasSize > 0)
    remove_droplets(f, minsize = drillRemoveGasSize, bubbles = true);

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

  `jetFormed` is file-scope (restart-persistent via `drillstate`, see init).
  */
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
    // Terminal relaxation: latch once the jet is FULLY formed — the first tip
    // droplet has shed and STAYED shed. Persistence (LATCH_STEPS consecutive
    // steps with more than one liquid component) guards against a transient
    // liquid fragment near reconnection faking n > 1 for a step or two and
    // relaxing the mesh mid-jet. (tipPinched is file-scope, restart-persistent
    // via drillstate; the counter is transient and re-counts after restart.)
    if (!tipPinched && jetFormed) {
      tipPinchSteps = (n > 1) ? tipPinchSteps + 1 : 0;
      if (tipPinchSteps >= LATCH_STEPS) tipPinched = 1;
    }

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

      /**
      ### Robust arm/fire reconnection latch (case-1006 lesson)

      The classic latch above (max-|k| on axis AND lowest point off axis)
      relies on a transient window between reconnection and the first
      entrained satellite; in case 1006 that window had zero width — a
      satellite pinned the lowest point on-axis from the first
      post-reconnection step, jetFormed never fired, and the focus cap held
      the whole jet phase at level 12 (the L14 release never happened).

      The fix uses the ROBUST base observable (previous step's g_rbase,
      satellite- and droplet-proof) with a two-stage arm/fire design,
      because r_base alone is ambiguous: the pre-reconnection dimple ring
      puts the outer-surface base off-axis (up to r ~ 0.11 around t ~ 0.456
      in case 1006) long before the singular instant.

      ARM in the final singular approach — uncapped demand has reached
      MAXlevel (1/kmax collapsed; in case 1005 this first happened at
      t = 0.46659, ~2e-3 before reconnection, and never during the dimple
      wander) AND the base is pinned on-axis (r_base < ARM_BAND; the dimple
      has r_base >= 0.018 whenever it is off-axis, so the two states are
      well separated). FIRE when, after arming, the base stays off-axis
      (> AXIS_BAND) for LATCH_STEPS consecutive steps — reconnection has
      been crossed and the jet base is opening for good. Self-heal clause:
      a restart that jumps straight into the developed jet phase (base
      off-axis at restore, so arming can never trigger) latches directly
      once r_base > RBASE_JET. (drillArmed/baseOffAxis are file-scope;
      drillArmed is restart-persistent via drillstate.)
      */
      if (!jetFormed) {
        if (!drillArmed && L >= MAXlevel && g_rbase > -900. && g_rbase < ARM_BAND)
          drillArmed = 1;
        if (drillArmed) {
          baseOffAxis = (g_rbase > AXIS_BAND) ? baseOffAxis + 1 : 0;
          if (baseOffAxis >= LATCH_STEPS) jetFormed = 1;
        }
        if (g_rbase > RBASE_JET) jetFormed = 1;
      }

      /**
      Pre-inception cap (case-1006 lesson): the cavity-focus collapse is a
      genuine singularity — 1/kmax -> 0 demands MAXlevel exactly at the
      singular instant, where deeper resolution is self-defeating (dt stalls,
      reconnection blows up; see the file header Status note). Cap the focus
      regime at drillMaxlevelFocus (12 is validated to step over the topology
      change) and release the full MAXlevel only once the inception latch
      fires: the erupting slender jet is fast but SMOOTH, so deep refinement
      is safe there — and that is where the q_jet(r_jet) observable wants it.
      */
      if (!jetFormed && drillMaxlevelFocus > 0 && L > drillMaxlevelFocus)
        L = drillMaxlevelFocus;
      /**
      ### Post-inception floor (case-2331 lesson)

      The pre-inception cap above stops the ramp chasing the focus
      singularity. Nothing stopped the ramp COARSENING once the jet exists,
      and the coarsening is irreversible. The demanded level is derived from
      the tracked jet-base radius `rb = rlow`, the radius of the deepest
      interfacial point of the MAIN tagged liquid body. While the jet base
      sits on the axis, `rb ~ Delta/2` and the ramp holds MAXlevel. The
      instant the main body's deepest point relocates off-axis — a detaching
      tip fragment leaving the main tag is sufficient, and nothing about that
      is a statement that jetting has ended — `rb` jumps by two orders of
      magnitude, `L` collapses, and the hysteresis walks the ceiling down one
      level per step to `drillMaxlevelStart`. At that level the slender jet
      cannot be represented, so `rlow` can never return to the axis and the
      probe can never re-demand resolution. Case 2331 (De = 0.04) lost the
      feature at t = 0.4770 (i = 4189: rb 1.22e-3 -> 3.78e-2 in one step),
      reached level 8 by i = 4291, and ran the remaining 1.02 capillary times
      there — 4738 steps against ~33 000 for its neighbours, and no resolved
      drop. `drillMinlevelJet` is that missing floor: once `jetFormed`, never
      coarsen below it. <=0 keeps the old behaviour exactly.
      */
      if (jetFormed && drillMinlevelJet > 0 && L < drillMinlevelJet)
        L = drillMinlevelJet;

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
  // Persist the drill latch state alongside the dump so a restart resumes
  // with the correct regime instead of resetting the latches (init reads it).
  if (pid() == 0) {
    FILE *fs = fopen("drillstate", "w");
    if (fs) {
      fprintf(fs, "%d %d %d\n", jetFormed, drillArmed, tipPinched);
      fclose(fs);
    }
  }
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

    // Blow-up gate, now a knob (params.keStopMax, historical default 1e2).
    // The threshold is ad hoc: a localised transient spike can exceed it and
    // still self-recover, while a genuine divergence also stalls dt — so a
    // relaxed gate plus an external dt/progress watchdog is a legitimate way
    // to force a run through the singular instant (case-1006 protocol).
    if (ke > params.keStopMax && i > 1e1) {
      fprintf(ferr, "The kinetic energy blew up (ke = %g > keStopMax = %g). Stopping simulation\n",
              ke, params.keStopMax);
      fp = fopen("log", "a");
      fprintf(fp, "The kinetic energy blew up (ke = %g > keStopMax = %g). Stopping simulation\n",
              ke, params.keStopMax);
      fclose(fp);
      dump(file = dumpFile);
      return 1;
    }

    // Check for energy dissipation below threshold
    if (ke < params.keStopMin && i > 1e1) {
      fprintf(ferr, "kinetic energy too small now! Stopping!\n");
      dump(file = dumpFile);
      fp = fopen("log", "a");
      fprintf(fp, "kinetic energy too small now! Stopping!\n");
      fclose(fp);
      return 1;
    }
  }
}
