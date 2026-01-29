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
├── default.params             # Default parameter file
├── sweep.params               # Parameter sweep configuration
├── src-local/                 # Shared shell libraries
│   ├── parse_params.sh        # Parameter file parsing
│   ├── sweep_utils.sh         # Sweep generation utilities
│   └── basilisk_version.sh    # Basilisk version pinning
├── simulationCases/           # Output directory
│   ├── burstingBubble.c       # Main source file (template)
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

### Required Parameters

| Parameter | Description | Typical Values |
|-----------|-------------|----------------|
| `CaseNo` | 4-digit case identifier (1000-9999) | 1000 |
| `Oh` | Ohnesorge number (viscosity) | 1e-3 to 1e-1 |
| `Bond` | Bond number (gravity) | 1e-3 |
| `MAXlevel` | Max refinement level | 10-12 |
| `tmax` | Simulation end time | 0.5-2.0 |
| `zWall` | Distance to bottom wall | 0.025-4.0 |

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
- **burstingBubble.c**: Preserved unless `--force` is used
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
