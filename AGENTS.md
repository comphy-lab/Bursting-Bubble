# Bursting Bubble Simulation - Developer Guide

See @README.md for the user-facing overview.

## Overview

This codebase simulates bursting bubble dynamics using Basilisk C, a finite volume solver for PDEs. The simulation uses adaptive mesh refinement (AMR) and supports both serial and parallel (MPI) execution.

## Architecture

### Two-Stage Execution Model

The simulation uses a two-stage execution model due to a Basilisk limitation:

1. **Stage 1 (Initialization)**: Generates a restart file containing the initial condition
   - Uses `distance.h` which is incompatible with MPI
   - Must run serial or with OpenMP (Linux only)
   - Short runtime (~5e-2 time units)

2. **Stage 2 (Full Simulation)**: Runs the actual simulation from the restart file
   - Can use MPI for parallelization
   - Longer runtime (configurable via `tmax`)

### Directory Structure

```
.
├── runSimulation.sh           # Main simulation runner (single case)
├── runParameterSweep.sh       # Local parameter sweep runner
├── runPostProcess-Ncases.sh   # Post-processing pipeline
├── runSweepHamilton.sbatch    # HPC sweep runner (Hamilton Stage 2)
├── runSweepHamilton-serial.sbatch   # HPC sweep runner (Hamilton Stage 1)
├── runSweepSnellius.sbatch    # HPC sweep runner (Snellius Stage 2)
├── runSweepSnellius-serial.sbatch   # HPC sweep runner (Snellius Stage 1)
├── default.params             # Default parameter file (Newtonian)
├── default-ve.params          # Oldroyd-B, usual wavelet AMR
├── default-ve-drill.params    # Oldroyd-B + drill (singularity path)
├── sweep.params               # Parameter sweep configuration
├── src-local/                 # Shared shell + C runtime libraries
│   ├── params.h               # C-side runtime parameter layer (case.params)
│   ├── parse_params.sh        # Parameter file parsing (shell layer)
│   ├── sweep_utils.sh         # Sweep generation utilities
│   ├── basilisk_version.sh    # Basilisk version pinning
│   ├── two-phaseVE.h          # Two-phase VE properties (from MultiRheoFlow)
│   ├── log-conform-viscoelastic-scalar-2D.h
│   └── log-conform-viscoelastic.h
├── simulationCases/           # Output directory
│   ├── burstingBubble.c                         # Newtonian, fixed-ceiling AMR
│   ├── burstingBubble-drillResolution.c         # Newtonian drill
│   ├── burstingBubbleVE.c                       # Oldroyd-B, usual AMR
│   ├── burstingBubbleVE-drillResolution.c       # Oldroyd-B + drill
│   ├── DataFiles/             # Initial condition data
│   └── <CaseNo>/              # Per-case output folders
└── postProcess/               # Post-processing scripts and helpers
```

### Shared Libraries (src-local/)

- **parse_params.sh**: Parse key=value parameter files, export as `PARAM_*` environment variables
- **sweep_utils.sh**: Generate Cartesian product of sweep parameters, progress tracking
- **basilisk_version.sh**: Centralized Basilisk version configuration (`BASILISK_REF`)

## Parameter Files

### Format

```
# Comments start with #
key=value
Oh=1e-2      # Inline comments allowed
Bond=1e-3
```

### Parameters

The simulation reads its configuration directly from a `case.params` file via the C-side
parameter layer in `src-local/params.h` (struct `SimulationParams`,
`parse_params_from_file`, `key=value` CLI overrides via `apply_cli_overrides`, a legacy
positional fallback `parse_params_from_cli`, plus `validate_params` / `print_params`).
Compiles must add `-I../../src-local` so `#include "params.h"` resolves.

The binary is invoked as:

```bash
./burstingBubble case.params [key=value ...]      # preferred; trailing tokens override the file
./burstingBubble <MAXlevel> <Oh> <Bond> <tmax> <zWall>   # legacy positional fallback
```

Stage 1 (restart generation) uses the override form `./burstingBubble case.params tmax=0.10`.

| Parameter | Group | Description | Default |
|-----------|-------|-------------|---------|
| `CaseNo` | case | 4-digit case identifier (1000-9999) | 1000 |
| `Solver` | case | source stem in `simulationCases/` (runner-only) | burstingBubble |
| `Oh` | physical | solvent Ohnesorge number, liquid | 1e-2 |
| `Bond` | physical | Bond number (gravity) | 1e-3 |
| `OhRatio` | physical | gas/liquid Ohnesorge ratio; `Oha = OhRatio*Oh` | 2e-2 |
| `De` | physical | Deborah number (liquid relaxation time); `0` = Newtonian | 0 |
| `Ec` | physical | elasto-capillary number (liquid modulus); `0` = Newtonian | 0 |
| `FILTERED` | numerical | compile-time smear of density/viscosity jumps (`-DFILTERED` when 1) | 1 |
| `zWall` | geometry | distance from bubble south pole to bottom wall | 0.05 |
| `MAXlevel` | space | maximum refinement level | 10 |
| `MINlevel` | space | far-field coarsening floor | 4 |
| `init_grid_level` | space | initial uniform grid level | 5 |
| `fErr` | space | wavelet tolerance on VOF `f` | 1e-3 |
| `VelErr` | space | wavelet tolerance on velocity | 1e-3 |
| `KErr` | space | wavelet tolerance on curvature | 1e-6 |
| `AErr` | space | wavelet tolerance on conformation `A_ij` (VE solvers) | 1e-3 |
| `CFL` | time | advective CFL number | 0.1 |
| `dtmax` | time | timestep ceiling (see note) | 1e-2 |
| `TOLERANCE` | time | Poisson/viscous solver tolerance | 1e-4 |
| `tmax` | time | simulation end time | 1.0 |
| `tsnap` | time | snapshot/restart dump interval | 1e-2 |

### Adaptive time & space resolution

Both mesh and timestep are adaptive. **`dtmax` is a ceiling, not a fixed step**: the
surface-tension scheme (`tension.h`) is time-explicit and reduces the step each iteration to
the capillary-wave limit `T = sqrt(rho_m * Delta_min^3 / (pi * sigma))` with
`rho_m = (rho1+rho2)/2`, so the effective step is set adaptively and scales with the finest
cell. Earlier versions hard-capped `dtmax` at `1e-5` (below that limit), pinning the step;
it is now a generous ceiling. `MINlevel` is an explicit far-field floor (previously the
implicit `MAXlevel-6` inside `adapt_wavelet`); the interface itself is always refined to
`MAXlevel` through the `fErr` criterion.

### Sweep File Format

```
BASE_CONFIG=default.params
CASE_START=1000
CASE_END=1003
SWEEP_Oh=0.01,0.02
SWEEP_Bond=0.001,0.002
```

Generates Cartesian product: 2 Oh values x 2 Bond values = 4 cases

## Workflow

### Local Development

```bash
# Single simulation (both stages)
./runSimulation.sh default.params

# Stage 1 only (generate restart)
./runSimulation.sh --stage1 default.params

# Stage 2 with MPI
./runSimulation.sh --stage2 --mpi 8 default.params

# Force overwrite existing files
./runSimulation.sh --force default.params
```

### Parameter Sweeps

```bash
# Local sweep (sequential)
./runParameterSweep.sh sweep.params

# Dry run to preview combinations
./runParameterSweep.sh --dry-run sweep.params

# Stage 1 only for all cases
./runParameterSweep.sh --stage1-only sweep.params
```

### HPC Submission (Snellius/Hamilton)

1. Run Stage 1 locally or with serial sbatch:
   ```bash
   sbatch runSweepSnellius-serial.sbatch
   ```

2. After Stage 1 completes, submit Stage 2:
   ```bash
   sbatch runSweepSnellius.sbatch
   ```

## File Preservation

The scripts preserve existing files for reruns:

- **case.params**: Preserved unless `--force` is used
- **solver source** (`Solver=` stem): Preserved unless `--force` is used
- **restart**: Stage 2 requires this file from Stage 1

This allows manual parameter/code edits between runs.

## Common Issues

### "restart file not found"
Run Stage 1 first to generate the restart file:
```bash
./runSimulation.sh --stage1 default.params
```

### "restart file is empty"
Stage 1 may have failed. Check the output logs.

### "qcc not found"
Install Basilisk first:
```bash
curl -sL https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/reset_install_basilisk-ref-locked.sh | bash -s -- --ref=v2026-01-13
```

### Stale parameters in reruns
Use `--force` to overwrite preserved files:
```bash
./runSimulation.sh --force default.params
```

## Code Style

- **Shell**: POSIX-compatible with bash extensions
- **Error handling**: `set -euo pipefail` for strict mode
- **Shellcheck**: All scripts should pass shellcheck
- **Comments**: Explain "why", not "what"
