/**
# Bursting Bubbles in Newtonian Fluids — jet-base tracking (LOGGING ONLY)

This is `burstingBubble.c` with extra per-log-step diagnostics that track the
Worthington jet-base / cavity-focus PROBE LOCATION on the fly. It is a pure
*instrumentation* extension:

- The physics, boundary conditions, properties, and time stepping are byte-for-
  byte identical to `burstingBubble.c`.
- The adaptive-mesh `adapt()` event is UNCHANGED — no refinement criterion is
  added or modified here. This file does NOT alter the mesh in any way.

The only behavioural change is inside `logWriting(i++)`: after the kinetic
energy `ke` is computed we evaluate the geometric jet-base probe and append its
coordinates (r_b, z_b) to the log. The detection algorithm and constants are
reproduced from the post-processing helper `getJetFoot.c`, with the
time-ordered inception latch realised in-solver as a static, forward-in-time
flag (`jetFormed`). Fluxes (q_jet, q_l) and the tip height (z_jet) are
intentionally NOT computed here — they are post-processing diagnostics
(getJetFoot.c); the solver only needs the probe location, which is the eventual
hook for adaptive refinement.

Probe (over the MAIN connected liquid body; detached drops excluded),
interfacial cells f in (1e-6, 1-1e-6) with y < RCAV:
  - (z_low, r_low) : globally lowest interfacial point (min axial x).
  - (z_maxk, r_maxk): point of max |curvature| with x < ZSURF_CURV.
Inception latch (never resets once set): jetFormed = 1 when
  rmaxk in [0, R_AXIS_K) and rlow > AXIS_BAND.
Probe selection: jetFormed -> (z_b,r_b) = (z_low,r_low)  [jet base];
                 else        (z_b,r_b) = (z_maxk,r_maxk)  [cavity focus].
Sentinel -1000 when no probe exists.

Coords: x = axial (= z), y = radial (= r >= 0). Dump carries only f and u.

@file burstingBubble-adaptiveResolution.c
@author Vatsal Sanjay (vatsal.sanjay@comphy-lab.org) / CoMPhy Lab
@version 2.0
@date Jan 04, 2025 (jet-base logging added Jun 2026)
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
    fprintf(ferr, "Cannot restored from a dump file!\n");
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
}

/**
## Adaptive Mesh Refinement

Refines the mesh based on gradients of key fields:
- Volume fraction
- Velocity components
- Curvature

The wavelet error tolerances (`fErr`, `VelErr`, `KErr`) and the refinement
band (`MINlevel` to `MAXlevel`) are runtime parameters. The interface is
always resolved to `MAXlevel` through the `fErr` criterion, while `MINlevel`
sets how coarse the far field is allowed to become.

NOTE: This event is intentionally IDENTICAL to the one in `burstingBubble.c`.
The jet-base diagnostics added in this file are logging-only and must not
influence refinement.
*/
event adapt(i++) {
  scalar KAPPA[];
  curvature(f, KAPPA);

  adapt_wavelet((scalar *){f, u.x, u.y, KAPPA},
    (double[]){fErr, VelErr, VelErr, KErr},
    MAXlevel, MINlevel);
}

/**
## Output Management

Creates periodic snapshots of the simulation state.
- Dumps restart files for simulation recovery
- Saves intermediate snapshots at regular intervals defined by `tsnap`
*/
event writingFiles(t = 0; t += tsnap; t <= tmax) {
  dump(file = dumpFile);
  sprintf(nameOut, "intermediate/snapshot-%5.4f", t);
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
## Simulation Logging (+ jet-base tracking)

Records key simulation data at each timestep:
- Iteration number, timestep size, current simulation time, kinetic energy
- Jet-base probe location: r_b, z_b (fluxes/tip are post-processing; getJetFoot.c)

Also performs the original safety checks (kinetic-energy blow-up / too-small),
which are left fully intact.

The jet-base block reproduces `getJetFoot.c`: candidate detection over the
main connected liquid body restricted to interfacial cells with y < RCAV
(`foreach(serial)`, exactly as in the helper), a forward-in-time inception
latch, regime-dependent probe selection, and an MPI-safe (reduction) band
flux integral. Computed EVERY step (per-step, for correctness).
*/
event logWriting(i++) {
  // Calculate kinetic energy
  double ke = 0.;
  foreach(reduction(+:ke)) {
    ke += (2*pi*y)*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }

  /**
  ### Jet-base / cavity-focus probe (logging only)
  */
  // Main connected liquid region (exclude detached drops).
  scalar d[];
  foreach() d[] = (f[] > 1e-4);
  int n = tag(d);
  int MainPhase = 0;
  if (n > 0) {
    double *sz = calloc(n, sizeof(double));
    foreach(serial) if (d[] > 0) sz[(int)d[] - 1] += 1.;
    double sm = -1.;
    for (int j = 0; j < n; j++) if (sz[j] > sm) { sm = sz[j]; MainPhase = j + 1; }
    free(sz);
  }

  scalar kappa[];
  curvature(f, kappa);

  // Candidates (mirror getJetFoot.c: serial sweep over interfacial cells).
  double zlow = HUGE, rlow = -1.;
  double zk = -1000., rk = -1000., kmax = -1.;
  foreach(serial) {
    if (f[] <= 1e-6 || f[] >= 1. - 1e-6) continue;   // interfacial only
    if (d[] != MainPhase) continue;
    if (y > RCAV) continue;
    if (x < zlow) { zlow = x; rlow = y; }
    if (x < ZSURF_CURV && kappa[] != nodata) {
      double ak = fabs(kappa[]);
      if (ak > kmax) { kmax = ak; zk = x; rk = y; }
    }
  }
  if (rlow < 0.) { zlow = -1000.; rlow = -1000.; }

  /**
  ### Inception latch (static, forward-in-time)

  In-solver analogue of the time-ordered latch in the post-processing
  consumer. Once the cavity focus has collapsed onto the axis (rmaxk small)
  while the lowest interfacial point has lifted off the axis (rlow large),
  we declare the jet "formed" and never reset.
  */
  static int jetFormed = 0;
  if (!jetFormed && rk >= 0 && rk < R_AXIS_K && rlow > AXIS_BAND)
    jetFormed = 1;

  /**
  ### Probe selection

  Rule 1 (jetFormed): probe the jet base = lowest interfacial point.
  Rule 2 (otherwise): probe the cavity focus = max |curvature| point.
  */
  double zb, rb;
  if (jetFormed) { zb = zlow; rb = rlow; }   // rule 1: jet base
  else           { zb = zk;   rb = rk;   }   // rule 2: cavity focus

  if (rb <= 0.) { zb = -1000.; rb = -1000.; }   // no probe found this step

  if (pid() == 0) {
    static FILE *fp;
    if (i == 0) {
      fprintf(ferr, "Level %d, Oh %2.1e, Oha %2.1e, Bo %4.3f, zWall %g, Ldomain %g\n",
              MAXlevel, Oh, Oha, Bond, zWall, Ldomain);
      fprintf(ferr, "i dt t ke r_b z_b\n");
      fp = fopen("log", "w");
      fprintf(fp, "Level %d, Oh %2.1e, Oha %2.1e, Bo %4.3f, zWall %g, Ldomain %g\n",
              MAXlevel, Oh, Oha, Bond, zWall, Ldomain);
      fprintf(fp, "i dt t ke r_b z_b\n");
      fprintf(fp, "%d %g %g %g %7.6e %7.6e\n",
              i, dt, t, ke, rb, zb);
      fclose(fp);
    } else {
      fp = fopen("log", "a");
      fprintf(fp, "%d %g %g %g %7.6e %7.6e\n",
              i, dt, t, ke, rb, zb);
      fclose(fp);
    }
    fprintf(ferr, "%d %g %g %g %7.6e %7.6e\n",
            i, dt, t, ke, rb, zb);

    assert(ke > -1e-10);

    // Check for energy blowup (numerical instability)
    if (ke > 1e2 && i > 1e1) {
      if (pid() == 0) {
        fprintf(ferr, "The kinetic energy blew up. Stopping simulation\n");
        fp = fopen("log", "a");
        fprintf(fp, "The kinetic energy blew up. Stopping simulation\n");
        fclose(fp);
        dump(file = dumpFile);
        return 1;
      }
    }
    assert(ke < 1e2);

    // Check for energy dissipation below threshold
    if (ke < 1e-6 && i > 1e1) {
      if (pid() == 0) {
        fprintf(ferr, "kinetic energy too small now! Stopping!\n");
        dump(file = dumpFile);
        fp = fopen("log", "a");
        fprintf(fp, "kinetic energy too small now! Stopping!\n");
        fclose(fp);
        return 1;
      }
    }
  }
}
