/**
# params.h

Runtime parameter management for the bursting-bubble simulation.

Defines a single configuration structure plus helpers to set defaults,
parse a `key=value` parameter file (`case.params`), apply `key=value`
command-line overrides, fall back to the legacy positional CLI, validate,
and print the resolved configuration.

This is the C-side analogue of the shell `parse_params.sh` layer: the shell
layer orchestrates case folders and sweeps, while this header lets the
simulation read its own configuration directly from `case.params`. Adding a
new tunable knob requires only (1) a struct field, (2) a default in
`set_default_params`, and (3) a line in `apply_param_kv` — which is what
makes the adaptive time/space resolution upgrade a configuration change
rather than a code edit.

## Author
Vatsal Sanjay (vatsal.sanjay@comphy-lab.org)
CoMPhy Lab, Durham University
*/

#ifndef PARAMS_H
#define PARAMS_H

#include <ctype.h>      // isspace()
#include <math.h>       // isfinite()
#include <sys/stat.h>   // mkdir()
#include <errno.h>      // errno

/**
## SimulationParams

All runtime parameters in one structure. Grouped into physical,
geometry, adaptive-space, adaptive-time, and time-control blocks.
*/
struct SimulationParams {
  // Case identification
  int CaseNo;            /**< Case number for folder naming (4-digit: 1000-9999) */

  // Physical parameters (dimensionless numbers)
  double Oh;             /**< Ohnesorge number (liquid solvent): mu_s/sqrt(rho*sigma*R) */
  double Bond;           /**< Bond number: rho*g*R^2/sigma */
  double OhRatio;        /**< Gas/liquid Ohnesorge ratio; Oha = OhRatio*Oh */
  double De;             /**< Deborah number (liquid): lambda / t_cap. 0 = Newtonian */
  double Ec;             /**< Elasto-capillary number (liquid): G / (sigma/R). 0 = Newtonian */
  int filtered;          /**< 1 = smear density/viscosity jumps; 0 = sharp VOF
                              properties. File key remains FILTERED. The
                              runner adds `-DFILTERED` only when this is 1
                              (`two-phase.h` / `two-phaseVE.h` use `#ifdef`).
                              The member cannot be named FILTERED: that
                              token is a compile macro when smearing is on. */

  // Geometry
  double zWall;          /**< Distance from bubble south pole to bottom wall */

  // Adaptive SPACE resolution
  int MAXlevel;          /**< Maximum refinement level (2^MAXlevel cells at finest) */
  int MINlevel;          /**< Minimum refinement level (far-field coarsening floor) */
  int init_grid_level;   /**< Initial uniform grid level (2^init_grid_level cells) */
  double fErr;           /**< Wavelet error tolerance on the VOF field f */
  double VelErr;         /**< Wavelet error tolerance on velocity components */
  double KErr;           /**< Wavelet error tolerance on interface curvature */
  double AErr;           /**< Wavelet error tolerance on conformation A_ij (VE solvers) */

  // Adaptive TIME resolution
  double CFL;            /**< Advective CFL number */
  double CFLelastic;     /**< Elastic-wave CFL safety factor. The polymeric
                              stress supports a shear wave of speed
                              sqrt(Gp*tr(A)/rho); nothing else in the timestep
                              selection accounts for it, and the first VE
                              campaign lost seven runs to exactly that
                              omission. 0 disables the condition. */
  double dtmax;          /**< Ceiling on the timestep; surface tension reduces it
                              to the capillary-wave limit each step (adaptive dt) */
  double TOLERANCE;      /**< Poisson/viscous solver convergence tolerance */

  // Time control
  double tmax;           /**< Maximum simulation time (capillary units) */
  double tsnap;          /**< Snapshot/restart dump interval */

  // Safety gates (kinetic-energy sanity checks in logWriting)
  double keStopMax;      /**< Stop when ke exceeds this (blow-up guard). The
                              historical 1e2 is ad hoc: a localised, transient
                              spike (e.g. at a refinement release near the
                              singular instant) can exceed it and still
                              self-recover, while a genuine divergence also
                              stalls dt. Relax deliberately (with a dt/progress
                              watchdog) to force a run through the singularity;
                              see the case-1005/1006 notes. */
  double keStopMin;      /**< Stop when ke falls below this (dissipated /
                              nothing left to compute) */

  // Drill adaptive-resolution trigger (feature-tracking AMR + time)
  // Only consumed by burstingBubble-drillResolution.c. The plain adaptive
  // solver ignores these. See that file's header for the mechanism.
  int drillAMR;              /**< Master switch: 1 = feature-tracking ramp of the
                                  local ceiling maxlevelLocal; 0 = pin at MAXlevel
                                  (reproduces the fixed-level reference run) */
  int drillMaxlevelStart;    /**< Coarsest level the ramp is allowed to fall to
                                  (ramp floor); the far field still coarsens to
                                  MINlevel via the wavelet criterion */
  double drillNcellsK;       /**< Cells the tracked length must span before the
                                  current level is deemed sufficient — cavity-focus
                                  (curvature-radius) regime, pre-inception */
  double drillNcellsJet;     /**< Same, jet-base-radius regime, post-inception */
  int drillRelaxLevel;       /**< Level to relax to after the first tip droplet
                                  sheds; <=0 disables relaxation (hold resolution
                                  on the receding base — the safe default) */
  int drillTsnapStages;      /**< 1 = tighten the snapshot interval as the mesh
                                  refines (adaptive time-resolution of output);
                                  0 = uniform tsnap */
  double drillTsnapMinFactor;/**< Floor on the staged snapshot interval as a
                                  fraction of tsnap (guards against a snapshot
                                  storm when start is far below MAXlevel) */
  int drillMaxlevelFocus;    /**< Pre-inception cap on the demanded level. The
                                  cavity-focus collapse is a genuine singularity:
                                  chasing it beyond the level that safely steps
                                  over the topology change stalls dt (CFL chases
                                  diverging u in ever-smaller cells) and blows up
                                  ke at reconnection (case 1005, MAXlevel=14,
                                  t=0.46887). Cap the focus regime here (12 is
                                  validated) and let the full MAXlevel loose only
                                  after the inception latch, where the slender
                                  jet is fast but smooth. <=0 disables (cap =
                                  MAXlevel, the pre-case-1006 behaviour). */
  int drillRemoveGasSize;    /**< Remove gas fragments (bubbles=true) smaller
                                  than this side length in cells each step:
                                  components below drillRemoveGasSize^2 cells
                                  (2D) are absorbed into the liquid. Kills the
                                  sub-grid gas wisps shed by the cavity
                                  reconnection that drive the CFL stall. Liquid
                                  droplets are never touched (shed tip droplets
                                  are physics). 0 disables. */
  int drillAssumeJet;        /**< 1 = force the inception latch (jetFormed) on
                                  at init after a successful restore. For
                                  restarts from LEGACY post-inception dumps
                                  that predate the drillstate file — without
                                  it such a restart can never re-latch and the
                                  focus cap silently re-binds the whole jet.
                                  Only ever sets the latch, never clears it;
                                  a drillstate file, when present, is read
                                  first. Default 0. */
};

/**
### set_default_params()

Populate the parameter structure with defaults representative of a
water-like bursting bubble. These reproduce the historical hard-coded
configuration, except that `dtmax` is now a generous ceiling (surface
tension sets the actual adaptive step) and `MINlevel` is an explicit
far-field floor (previously the implicit `MAXlevel-6`).
*/
static inline void set_default_params(struct SimulationParams *p) {
  // Case identification
  p->CaseNo = 1000;

  // Physical (water-air-like)
  p->Oh = 1.0e-2;
  p->Bond = 1.0e-3;
  p->OhRatio = 2.0e-2;
  p->De = 0.0;           // Newtonian unless a VE solver is selected
  p->Ec = 0.0;
  p->filtered = 1;       // historical two-phase smear; set 0 for sharp properties

  // Geometry
  p->zWall = 0.05;

  // Adaptive space
  p->MAXlevel = 10;
  p->MINlevel = 4;
  p->init_grid_level = 5;
  p->fErr = 1.0e-3;
  p->VelErr = 1.0e-3;
  p->KErr = 1.0e-6;
  p->AErr = 1.0e-3;

  // Adaptive time
  p->CFL = 0.1;
  p->CFLelastic = 0.25;  // margin demonstrated to cross the cavity-focus instant
  p->dtmax = 1.0e-2;
  p->TOLERANCE = 1.0e-4;
  p->keStopMax = 1.0e2;    // historical blow-up gate (ad hoc; see struct note)
  p->keStopMin = 1.0e-6;

  // Time control
  p->tmax = 1.0;
  p->tsnap = 1.0e-2;

  // Drill adaptive-resolution trigger (safe defaults: track the singularity,
  // never relax, mildly stage the output cadence)
  p->drillAMR = 1;
  p->drillMaxlevelStart = 8;
  p->drillNcellsK = 5.0;
  p->drillNcellsJet = 5.0;
  p->drillRelaxLevel = -1;      // relaxation disabled by default
  p->drillTsnapStages = 1;
  p->drillTsnapMinFactor = 0.1;
  p->drillMaxlevelFocus = -1;   // no pre-inception cap by default
  p->drillRemoveGasSize = 0;    // gas-fragment cleanup off by default
  p->drillAssumeJet = 0;        // don't assume a formed jet on restore
}

/**
### apply_param_kv()

Apply a single `key`/`value` assignment to the parameter structure.

Shared by file parsing and command-line override parsing so the
key dispatch lives in exactly one place.

#### Returns
- `1` if the key was recognised
- `0` for an unknown key (caller may warn)
*/
static inline int apply_param_kv(const char *key, const char *value,
                                 struct SimulationParams *p) {
  if      (strcmp(key, "CaseNo")          == 0) p->CaseNo = atoi(value);
  else if (strcmp(key, "Solver")          == 0) return 1; /* runner-only */
  else if (strcmp(key, "Oh")              == 0) p->Oh = atof(value);
  else if (strcmp(key, "Bond")            == 0) p->Bond = atof(value);
  else if (strcmp(key, "OhRatio")         == 0) p->OhRatio = atof(value);
  else if (strcmp(key, "De")              == 0) p->De = atof(value);
  else if (strcmp(key, "Ec")              == 0) p->Ec = atof(value);
  else if (strcmp(key, "FILTERED")        == 0) p->filtered = atoi(value);
  else if (strcmp(key, "zWall")           == 0) p->zWall = atof(value);
  else if (strcmp(key, "MAXlevel")        == 0) p->MAXlevel = atoi(value);
  else if (strcmp(key, "MINlevel")        == 0) p->MINlevel = atoi(value);
  else if (strcmp(key, "init_grid_level") == 0) p->init_grid_level = atoi(value);
  else if (strcmp(key, "fErr")            == 0) p->fErr = atof(value);
  else if (strcmp(key, "VelErr")          == 0) p->VelErr = atof(value);
  else if (strcmp(key, "KErr")            == 0) p->KErr = atof(value);
  else if (strcmp(key, "AErr")            == 0) p->AErr = atof(value);
  else if (strcmp(key, "CFL")             == 0) p->CFL = atof(value);
  else if (strcmp(key, "CFLelastic")      == 0) p->CFLelastic = atof(value);
  else if (strcmp(key, "dtmax")           == 0) p->dtmax = atof(value);
  else if (strcmp(key, "TOLERANCE")       == 0) p->TOLERANCE = atof(value);
  else if (strcmp(key, "keStopMax")       == 0) p->keStopMax = atof(value);
  else if (strcmp(key, "keStopMin")       == 0) p->keStopMin = atof(value);
  else if (strcmp(key, "tmax")            == 0) p->tmax = atof(value);
  else if (strcmp(key, "tsnap")           == 0) p->tsnap = atof(value);
  else if (strcmp(key, "drillAMR")            == 0) p->drillAMR = atoi(value);
  else if (strcmp(key, "drillMaxlevelStart")  == 0) p->drillMaxlevelStart = atoi(value);
  else if (strcmp(key, "drillNcellsK")        == 0) p->drillNcellsK = atof(value);
  else if (strcmp(key, "drillNcellsJet")      == 0) p->drillNcellsJet = atof(value);
  else if (strcmp(key, "drillRelaxLevel")     == 0) p->drillRelaxLevel = atoi(value);
  else if (strcmp(key, "drillTsnapStages")    == 0) p->drillTsnapStages = atoi(value);
  else if (strcmp(key, "drillTsnapMinFactor") == 0) p->drillTsnapMinFactor = atof(value);
  else if (strcmp(key, "drillMaxlevelFocus")  == 0) p->drillMaxlevelFocus = atoi(value);
  else if (strcmp(key, "drillRemoveGasSize")  == 0) p->drillRemoveGasSize = atoi(value);
  else if (strcmp(key, "drillAssumeJet")      == 0) p->drillAssumeJet = atoi(value);
  else return 0;
  return 1;
}

/**
### trim_inplace()

Strip leading/trailing whitespace from a string in place and return a
pointer to the trimmed start.
*/
static inline char *trim_inplace(char *s) {
  while (*s && isspace((unsigned char)*s)) s++;
  if (*s == '\0') return s;
  char *end = s + strlen(s) - 1;
  while (end > s && isspace((unsigned char)*end)) *end-- = '\0';
  return s;
}

/**
### parse_params_from_file()

Parse parameters from a `key=value` configuration file. Comments (`#`),
inline comments, and blank lines are ignored. Unknown keys produce a
warning but do not abort (forward compatibility with the shell layer).

#### Returns
- `0` on success
- `-1` if the file cannot be opened
*/
static inline int parse_params_from_file(const char *filename,
                                         struct SimulationParams *p) {
  FILE *fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "ERROR: Cannot open parameter file: %s\n", filename);
    return -1;
  }

  char line[512];
  int line_num = 0;
  while (fgets(line, sizeof(line), fp)) {
    line_num++;

    // Strip inline comments
    char *hash = strchr(line, '#');
    if (hash) *hash = '\0';

    // Skip blank lines
    char *start = trim_inplace(line);
    if (*start == '\0') continue;

    // Split on first '='
    char *eq = strchr(start, '=');
    if (!eq) continue;
    *eq = '\0';
    char *key = trim_inplace(start);
    char *value = trim_inplace(eq + 1);
    if (*key == '\0' || *value == '\0') continue;

    if (!apply_param_kv(key, value, p))
      fprintf(stderr, "WARNING: Unknown parameter '%s' at %s:%d\n",
              key, filename, line_num);
  }

  fclose(fp);
  return 0;
}

/**
### apply_cli_overrides()

Apply `key=value` tokens from `argv[start..argc)` on top of the already
parsed configuration. Enables stage-specific overrides such as
`./burstingBubble case.params tmax=0.10` for the Stage 1 restart run.
*/
static inline void apply_cli_overrides(int argc, char **argv, int start,
                                       struct SimulationParams *p) {
  for (int k = start; k < argc; k++) {
    char buf[256];
    strncpy(buf, argv[k], sizeof(buf) - 1);
    buf[sizeof(buf) - 1] = '\0';
    char *eq = strchr(buf, '=');
    if (!eq) {
      fprintf(stderr, "WARNING: Ignoring malformed override '%s' (expected key=value)\n",
              argv[k]);
      continue;
    }
    *eq = '\0';
    char *key = trim_inplace(buf);
    char *value = trim_inplace(eq + 1);
    if (!apply_param_kv(key, value, p))
      fprintf(stderr, "WARNING: Unknown override key '%s'\n", key);
  }
}

/**
### parse_params_from_cli()

Legacy positional fallback so historical invocations keep working:

`./burstingBubble <MAXlevel> <Oh> <Bond> <tmax> <zWall>`

New adaptive-resolution knobs take their defaults in this mode. Prefer the
`case.params` file mode for full control.

#### Returns
- `0` on success
- `-1` if too few arguments
*/
static inline int parse_params_from_cli(int argc, char **argv,
                                        struct SimulationParams *p) {
  if (argc < 6) {
    fprintf(stderr, "ERROR: Insufficient command line arguments\n");
    fprintf(stderr, "Legacy form: %s <MAXlevel> <Oh> <Bond> <tmax> <zWall>\n", argv[0]);
    fprintf(stderr, "Preferred:   %s <case.params> [key=value ...]\n", argv[0]);
    return -1;
  }
  p->MAXlevel = atoi(argv[1]);
  p->Oh = atof(argv[2]);
  p->Bond = atof(argv[3]);
  p->tmax = atof(argv[4]);
  p->zWall = atof(argv[5]);
  return 0;
}

/**
### params_init()

Single entry point. Detects file mode versus legacy positional mode:

- File mode (preferred): `argv[1]` is a readable file or ends in `.params`.
  The file is parsed, then any further `key=value` tokens override it.
- Legacy mode: otherwise the first five positional arguments are read.

Always starts from `set_default_params`, so unspecified knobs are defaulted.

#### Returns
- `0` on success
- `-1` on error (caller should abort the run)
*/
static inline int params_init(int argc, char **argv,
                              struct SimulationParams *p) {
  set_default_params(p);

  if (argc < 2) {
    fprintf(stderr, "ERROR: Missing arguments\n");
    fprintf(stderr, "Usage: %s <case.params> [key=value ...]\n", argv[0]);
    fprintf(stderr, "   or: %s <MAXlevel> <Oh> <Bond> <tmax> <zWall>  (legacy)\n", argv[0]);
    return -1;
  }

  const char *a1 = argv[1];
  int is_file = 0;
  FILE *probe = fopen(a1, "r");
  if (probe) { fclose(probe); is_file = 1; }
  if (!is_file) {
    size_t n = strlen(a1);
    if (n >= 7 && strcmp(a1 + n - 7, ".params") == 0) is_file = 1;
  }

  if (is_file) {
    if (parse_params_from_file(a1, p) != 0) return -1;
    apply_cli_overrides(argc, argv, 2, p);
  } else {
    if (parse_params_from_cli(argc, argv, p) != 0) return -1;
  }
  return 0;
}

/**
### validate_params()

Check physical and numerical consistency.

#### Returns
- `1` if valid
- `0` if invalid
*/
static inline int validate_params(const struct SimulationParams *p) {
  int valid = 1;

  if (p->CaseNo < 1000 || p->CaseNo > 9999) {
    fprintf(stderr, "ERROR: CaseNo must be 4-digit (1000-9999), got %d\n", p->CaseNo);
    valid = 0;
  }
  if (p->Oh <= 0) {
    fprintf(stderr, "ERROR: Oh must be positive (Oh = %g)\n", p->Oh);
    valid = 0;
  }
  if (p->OhRatio <= 0) {
    fprintf(stderr, "ERROR: OhRatio must be positive (OhRatio = %g)\n", p->OhRatio);
    valid = 0;
  }
  if (p->De < 0) {
    fprintf(stderr, "ERROR: De must be non-negative (De = %g)\n", p->De);
    valid = 0;
  }
  if (p->Ec < 0) {
    fprintf(stderr, "ERROR: Ec must be non-negative (Ec = %g)\n", p->Ec);
    valid = 0;
  }
  if (!isfinite(p->De) || !isfinite(p->Ec)) {
    fprintf(stderr, "ERROR: De and Ec must be finite (De = %g, Ec = %g)\n",
            p->De, p->Ec);
    valid = 0;
  }
  if (p->Bond < 0) {
    fprintf(stderr, "ERROR: Bond must be non-negative (Bond = %g)\n", p->Bond);
    valid = 0;
  }
  if (p->filtered != 0 && p->filtered != 1) {
    fprintf(stderr, "ERROR: FILTERED must be 0 or 1 (FILTERED = %d)\n", p->filtered);
    valid = 0;
  }
  if (p->MAXlevel < p->MINlevel) {
    fprintf(stderr, "ERROR: MAXlevel (%d) must be >= MINlevel (%d)\n",
            p->MAXlevel, p->MINlevel);
    valid = 0;
  }
  if (p->MINlevel < 2) {
    fprintf(stderr, "ERROR: MINlevel (%d) must be >= 2\n", p->MINlevel);
    valid = 0;
  }
  if (p->init_grid_level < p->MINlevel || p->init_grid_level > p->MAXlevel) {
    fprintf(stderr, "WARNING: init_grid_level (%d) outside [MINlevel, MAXlevel] = [%d, %d]\n",
            p->init_grid_level, p->MINlevel, p->MAXlevel);
  }
  if (p->MAXlevel > 15) {
    fprintf(stderr, "WARNING: Very high MAXlevel (%d) may exhaust memory\n", p->MAXlevel);
  }
  if (p->fErr <= 0 || p->VelErr <= 0 || p->KErr <= 0 || p->AErr <= 0) {
    fprintf(stderr, "ERROR: Wavelet error tolerances must be positive\n");
    valid = 0;
  }
  if (p->CFL <= 0 || p->CFL > 1) {
    fprintf(stderr, "ERROR: CFL must be in (0, 1] (CFL = %g)\n", p->CFL);
    valid = 0;
  }
  if (p->CFLelastic < 0 || p->CFLelastic > 1) {
    fprintf(stderr, "ERROR: CFLelastic must be in [0, 1] (0 disables) (CFLelastic = %g)\n",
            p->CFLelastic);
    valid = 0;
  }
  if (p->dtmax <= 0) {
    fprintf(stderr, "ERROR: dtmax must be positive (dtmax = %g)\n", p->dtmax);
    valid = 0;
  }
  if (p->TOLERANCE <= 0) {
    fprintf(stderr, "ERROR: TOLERANCE must be positive (TOLERANCE = %g)\n", p->TOLERANCE);
    valid = 0;
  }
  if (p->keStopMax <= 0 || p->keStopMin < 0 || p->keStopMin >= p->keStopMax) {
    fprintf(stderr, "ERROR: need 0 <= keStopMin < keStopMax (keStopMin = %g, keStopMax = %g)\n",
            p->keStopMin, p->keStopMax);
    valid = 0;
  }
  if (p->tmax <= 0) {
    fprintf(stderr, "ERROR: tmax must be positive (tmax = %g)\n", p->tmax);
    valid = 0;
  }
  if (p->tsnap <= 0 || p->tsnap > p->tmax) {
    fprintf(stderr, "ERROR: Invalid tsnap (tsnap = %g, tmax = %g)\n", p->tsnap, p->tmax);
    valid = 0;
  }
  // Drill trigger consistency (only meaningful for the drill solver, but a
  // malformed value should still fail fast rather than silently mis-refine).
  if (p->drillAMR) {
    if (p->drillMaxlevelStart < p->MINlevel || p->drillMaxlevelStart > p->MAXlevel) {
      fprintf(stderr, "ERROR: drillMaxlevelStart (%d) must be in [MINlevel, MAXlevel] = [%d, %d]\n",
              p->drillMaxlevelStart, p->MINlevel, p->MAXlevel);
      valid = 0;
    }
    if (p->drillNcellsK <= 0 || p->drillNcellsJet <= 0) {
      fprintf(stderr, "ERROR: drillNcellsK/drillNcellsJet must be positive (%g / %g)\n",
              p->drillNcellsK, p->drillNcellsJet);
      valid = 0;
    }
    if (p->drillRelaxLevel > 0 &&
        (p->drillRelaxLevel < p->MINlevel || p->drillRelaxLevel > p->MAXlevel)) {
      fprintf(stderr, "ERROR: drillRelaxLevel (%d) must be <=0 (disabled) or in [MINlevel, MAXlevel]\n",
              p->drillRelaxLevel);
      valid = 0;
    }
    if (p->drillMaxlevelFocus > 0 &&
        (p->drillMaxlevelFocus < p->drillMaxlevelStart || p->drillMaxlevelFocus > p->MAXlevel)) {
      fprintf(stderr, "ERROR: drillMaxlevelFocus (%d) must be <=0 (disabled) or in [drillMaxlevelStart, MAXlevel] = [%d, %d]\n",
              p->drillMaxlevelFocus, p->drillMaxlevelStart, p->MAXlevel);
      valid = 0;
    }
    if (p->drillRemoveGasSize < 0) {
      fprintf(stderr, "ERROR: drillRemoveGasSize (%d) must be >= 0 (0 disables)\n",
              p->drillRemoveGasSize);
      valid = 0;
    }
    if (p->drillAssumeJet != 0 && p->drillAssumeJet != 1) {
      fprintf(stderr, "ERROR: drillAssumeJet (%d) must be 0 or 1\n",
              p->drillAssumeJet);
      valid = 0;
    }
    if (p->drillTsnapMinFactor <= 0 || p->drillTsnapMinFactor > 1) {
      fprintf(stderr, "ERROR: drillTsnapMinFactor must be in (0, 1] (%g)\n", p->drillTsnapMinFactor);
      valid = 0;
    }
  }
  return valid;
}

/**
### print_params()

Print a formatted summary for logging and reproducibility.
*/
static inline void print_params(const struct SimulationParams *p, FILE *fp) {
  fprintf(fp, "\n========================================\n");
  fprintf(fp, "Bursting Bubble Simulation Configuration\n");
  fprintf(fp, "========================================\n");
  fprintf(fp, "Case Number:              %04d\n", p->CaseNo);
  fprintf(fp, "Physical Parameters:\n");
  fprintf(fp, "  Ohnesorge (solvent):    %g\n", p->Oh);
  fprintf(fp, "  Ohnesorge (gas):        %g  (OhRatio=%g)\n", p->OhRatio * p->Oh, p->OhRatio);
  fprintf(fp, "  Bond number:            %g\n", p->Bond);
  fprintf(fp, "  Deborah (liquid):       %g\n", p->De);
  fprintf(fp, "  Elasto-capillary:       %g\n", p->Ec);
  fprintf(fp, "  Polymeric Oh (Ec*De):   %g\n", p->Ec * p->De);
  fprintf(fp, "  FILTERED (compile):     %d\n", p->filtered);
  fprintf(fp, "Geometry:\n");
  fprintf(fp, "  zWall:                  %g\n", p->zWall);
  fprintf(fp, "Adaptive Space:\n");
  fprintf(fp, "  Levels (min/max):       %d / %d\n", p->MINlevel, p->MAXlevel);
  fprintf(fp, "  Initial grid level:     %d (2^%d = %d cells)\n",
          p->init_grid_level, p->init_grid_level, 1 << p->init_grid_level);
  fprintf(fp, "  Error tol (f/Vel/K/A):  %g / %g / %g / %g\n",
          p->fErr, p->VelErr, p->KErr, p->AErr);
  fprintf(fp, "Adaptive Time:\n");
  fprintf(fp, "  CFL:                    %g\n", p->CFL);
  if (p->CFLelastic > 0)
    fprintf(fp, "  CFL (elastic wave):     %g\n", p->CFLelastic);
  else
    fprintf(fp, "  CFL (elastic wave):     DISABLED\n");
  fprintf(fp, "  dtmax (ceiling):        %g\n", p->dtmax);
  fprintf(fp, "  Solver TOLERANCE:       %g\n", p->TOLERANCE);
  fprintf(fp, "  ke stop gates (min/max): %g / %g\n", p->keStopMin, p->keStopMax);
  fprintf(fp, "Time Control:\n");
  fprintf(fp, "  tmax:                   %g\n", p->tmax);
  fprintf(fp, "  tsnap:                  %g\n", p->tsnap);
  fprintf(fp, "Drill Trigger (drill solver only):\n");
  if (p->drillAMR) {
    fprintf(fp, "  drillAMR:               ON\n");
    fprintf(fp, "  maxlevel start/floor:   %d\n", p->drillMaxlevelStart);
    fprintf(fp, "  cells/feature (K/jet):  %g / %g\n", p->drillNcellsK, p->drillNcellsJet);
    if (p->drillRelaxLevel > 0)
      fprintf(fp, "  relax level (post-pinch): %d\n", p->drillRelaxLevel);
    else
      fprintf(fp, "  relax level (post-pinch): disabled\n");
    fprintf(fp, "  staged tsnap:           %s (floor factor %g)\n",
            p->drillTsnapStages ? "ON" : "OFF", p->drillTsnapMinFactor);
    if (p->drillMaxlevelFocus > 0)
      fprintf(fp, "  focus (pre-incept) cap: %d\n", p->drillMaxlevelFocus);
    else
      fprintf(fp, "  focus (pre-incept) cap: disabled\n");
    if (p->drillRemoveGasSize > 0)
      fprintf(fp, "  gas-wisp removal:       < %d^2 cells\n", p->drillRemoveGasSize);
    else
      fprintf(fp, "  gas-wisp removal:       OFF\n");
    if (p->drillAssumeJet)
      fprintf(fp, "  assume jet on restore:  ON\n");
  } else {
    fprintf(fp, "  drillAMR:               OFF (pinned at MAXlevel = %d)\n", p->MAXlevel);
  }
  fprintf(fp, "========================================\n\n");
  fflush(fp);
}

#endif // PARAMS_H
