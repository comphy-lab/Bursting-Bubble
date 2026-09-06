/**
# Bursting Bubbles in Viscoelastic Fluids — usual AMR

Oldroyd-B counterpart of `burstingBubble.c`. Geometry, two-stage init,
runtime `case.params`, and wavelet AMR are the same; the extra physics is
the CoMPhy log-conformation solver vendored from MultiRheoFlow
(`src-local/log-conform-viscoelastic-scalar-2D.h` + `two-phaseVE.h`).

The liquid is an Oldroyd-B fluid with solvent Ohnesorge `Oh`, relaxation
time `De`, and elastic modulus `Ec`. Gas remains Newtonian
(`lambda2 = G2 = 0`). Polymeric viscosity is the derived product
`Oh_p = Ec * De`; keep `Oh`, `De`, and `Ec` as independent controls.

`De = 0` or `Ec = 0` recovers the Newtonian limit of this solver (identity
conformation, zero polymeric stress). For a bit-identical Newtonian
baseline use `burstingBubble.c`, not this file.

Mesh: fixed-ceiling wavelet AMR on `f`, `u`, curvature, and the
conformation components `A11`, `A12`, `A22`, `AThTh`. For
feature-tracking resolution into the cavity-focus / jet singularity use
`burstingBubbleVE-drillResolution.c`.

@file burstingBubbleVE.c
@author Vatsal Sanjay (vatsal.sanjay@comphy-lab.org) / CoMPhy Lab
@version 1.0
@date Aug 27, 2026
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "log-conform-viscoelastic-scalar-2D.h"

/**
## Solver Configuration

`FILTERED` is a compile-time switch from `case.params`. The runner
adds `-DFILTERED` when `FILTERED=1`. `two-phaseVE.h` smears density
and viscosity jumps only when that macro is defined.
*/
#include "two-phaseVE.h"
#include "navier-stokes/conserving.h"
#include "tension.h"

#if !_MPI
#include "distance.h"
#endif

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
//   Oh  -> solvent Ohnesorge number (liquid)
//   Oha -> Ohnesorge number (gas) = OhRatio * Oh
//   De  -> Deborah number (liquid relaxation time)
//   Ec  -> elasto-capillary number (liquid elastic modulus)
double Oh, Oha, De, Ec, Bond, tmax;

// Domain parameters:
//   zWall   -> distance from bubble south pole to bottom wall
//   Ldomain -> computed domain size: min(zWall + 6.0, 16.0)
double zWall, Ldomain;

// Adaptive-resolution controls (set from params in main)
double fErr, VelErr, KErr, AErr;

// tsnap -> snapshot/restart dump interval. Needs a non-zero static initial
// value: Basilisk classifies event expressions (e.g. `t += tsnap`) before
// main() runs, and a zero increment would be misread as a second condition.
// main() overrides this with params.tsnap for the actual firing interval.
double tsnap = 1e-2;

char nameOut[80], dumpFile[80];

/**
## Main Function

Reads parameters, configures the domain and fluid properties, and starts
the run.

- Parses `case.params` (or legacy positional CLI) into `params`
- Validates the configuration before allocating the grid
- Sets up the physical domain with appropriate dimensions
- Configures Newtonian and Oldroyd-B properties for both phases
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
  De = params.De;
  Ec = params.Ec;
  Bond = params.Bond;
  tmax = params.tmax;
  zWall = params.zWall;
  fErr = params.fErr;
  VelErr = params.VelErr;
  KErr = params.KErr;
  AErr = params.AErr;
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
  the ceiling is therefore a safety limit and the effective step is adaptive
  (it scales with the finest cell size and the resolved physics).

  The ceiling must be written to `DT`, not to `dtmax`. `centered.h` carries
  `event set_dtmax (i++,last) dtmax = DT;`, so anything assigned to `dtmax`
  here is discarded on the first step and the knob silently does nothing.
  */
  CFL = params.CFL;
  CFL_elastic = params.CFLelastic;   // elastic-wave limit, see two-phaseVE.h
  CFL_conform = params.CFLconform;   // conformation-source limit, see two-phaseVE.h
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
  - `mu1`, `mu2`: Solvent viscosity of liquid and gas phases
  - `lambda1`, `G1`: Liquid relaxation time and elastic modulus (`De`, `Ec`)
  - gas remains Newtonian (`lambda2 = G2 = 0`)
  */
  rho1 = 1., rho2 = 1e-3;
  mu1 = Oh, mu2 = Oha;
  lambda1 = De; lambda2 = 0.;
  G1 = Ec; G2 = 0.;

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

Fresh init leaves the conformation at identity (unstressed polymers),
which is the log-conformation default.
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
}

/**
## Adaptive Mesh Refinement

Refines the mesh based on gradients of key fields:
- Volume fraction
- Velocity components
- Curvature
- Conformation-tensor components (liquid polymeric stress layers)

The wavelet error tolerances (`fErr`, `VelErr`, `KErr`, `AErr`) and the
refinement band (`MINlevel` to `MAXlevel`) are runtime parameters. The
interface is always resolved to `MAXlevel` through the `fErr` criterion,
while `MINlevel` sets how coarse the far field is allowed to become.
*/
event adapt(i++) {
  scalar KAPPA[];
  curvature(f, KAPPA);

  adapt_wavelet((scalar *){f, u.x, u.y, KAPPA, A11, A12, A22, AThTh},
    (double[]){fErr, VelErr, VelErr, KErr, AErr, AErr, AErr, AErr},
    MAXlevel, MINlevel);
}

/**
## Output Management

Creates periodic snapshots of the simulation state.
- Dumps restart files for simulation recovery
- Saves intermediate snapshots at regular intervals defined by `tsnap`

Dumps carry `f`, `u`, and the conformation / polymeric-stress fields.
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
    fprintf(ferr, "Level %d, De %2.1e, Ec %2.1e, Oh %2.1e, Oha %2.1e, Bo %4.3f, zWall %g, Ldomain %g\n",
            MAXlevel, De, Ec, Oh, Oha, Bond, zWall, Ldomain);
}

/**
## Simulation Logging

Records key simulation data at each timestep:
- Iteration number
- Timestep size
- Current simulation time
- Kinetic energy

Also performs safety checks:
- Monitors kinetic energy for stability
- Terminates simulation if energy becomes too high or too low
- Creates log files for post-processing analysis
*/
event logWriting(i++) {
  // Calculate kinetic energy
  double ke = 0.;
  foreach(reduction(+:ke)) {
    ke += (2*pi*y)*(0.5*rho(f[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
  }

  int stopBlowUp = (ke > params.keStopMax && i > 1e1);
  int stopTooSmall = (ke < params.keStopMin && i > 1e1);

  if (pid() == 0) {
    static FILE *fp;
    if (i == 0) {
      fprintf(ferr, "Level %d, De %2.1e, Ec %2.1e, Oh %2.1e, Oha %2.1e, Bo %4.3f, zWall %g, Ldomain %g\n",
              MAXlevel, De, Ec, Oh, Oha, Bond, zWall, Ldomain);
      fprintf(ferr, "i dt t ke\n");
      fp = fopen("log", "w");
      fprintf(fp, "Level %d, De %2.1e, Ec %2.1e, Oh %2.1e, Oha %2.1e, Bo %4.3f, zWall %g, Ldomain %g\n",
              MAXlevel, De, Ec, Oh, Oha, Bond, zWall, Ldomain);
      fprintf(fp, "i dt t ke\n");
      fprintf(fp, "%d %g %g %g\n", i, dt, t, ke);
      fclose(fp);
    } else {
      fp = fopen("log", "a");
      fprintf(fp, "%d %g %g %g\n", i, dt, t, ke);
      fclose(fp);
    }
    fprintf(ferr, "%d %g %g %g\n", i, dt, t, ke);

    assert(ke > -1e-10);

    if (stopBlowUp) {
      fprintf(ferr, "The kinetic energy blew up (ke = %g > keStopMax = %g). Stopping simulation\n",
              ke, params.keStopMax);
      fp = fopen("log", "a");
      fprintf(fp, "The kinetic energy blew up (ke = %g > keStopMax = %g). Stopping simulation\n",
              ke, params.keStopMax);
      fclose(fp);
    }

    if (stopTooSmall) {
      fprintf(ferr, "kinetic energy too small now! Stopping!\n");
      fp = fopen("log", "a");
      fprintf(fp, "kinetic energy too small now! Stopping!\n");
      fclose(fp);
    }
  }

  if (stopBlowUp || stopTooSmall) {
    dump(file = dumpFile);
    return 1;
  }
}
