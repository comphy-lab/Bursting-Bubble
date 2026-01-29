# Bursting Bubble Simulations

Computational fluid dynamics simulations for bursting bubble studies using the Basilisk C framework.

## Basilisk (Required)

First-time install (or reinstall):
```bash
curl -sL https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/reset_install_basilisk-ref-locked.sh | bash -s -- --ref=v2026-01-13 --hard
```

Subsequent runs (reuses existing `basilisk/` if same ref):
```bash
curl -sL https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/reset_install_basilisk-ref-locked.sh | bash -s -- --ref=v2026-01-13
```

> **Note**: Replace `v2026-01-13` with the [latest release tag](https://github.com/comphy-lab/basilisk-C/releases).

## Overview

This repository contains axisymmetric two-phase flow simulations with adaptive mesh refinement for studying bubble bursting phenomena. The simulations use the Volume-of-Fluid (VOF) method to track the interface between the bubble and surrounding fluid, with automatic mesh refinement focused on regions of interest.

## Quick Start

### Single Simulation

```bash
# Edit parameters
vim default.params      # Set CaseNo, Oh, Bond, etc.

# Run simulation (serial)
./runSimulation.sh

# Run with MPI (4 cores)
./runSimulation.sh --mpi
```

### Parameter Sweep

```bash
# Configure sweep
vim sweep.params        # Set CASE_START, CASE_END, sweep variables

# Run sweep (serial)
./runParameterSweep.sh

# Run sweep with MPI (4 cores per case)
./runParameterSweep.sh --mpi
```

## Repository Structure

```
├── src-local/                     Modular helper files
│   ├── parse_params.sh            Parameter parsing utilities
│   ├── sweep_utils.sh             Sweep generation utilities
│   └── basilisk_version.sh        Centralized version pinning
├── postProcess/                   Post-processing tools and visualization
│   ├── getData.c                  Field extraction on structured grids
│   ├── getFacet.c                 Interface geometry extraction
│   └── Video.py                   Frame-by-frame visualization pipeline
├── simulationCases/               Case-based simulation outputs
│   ├── burstingBubble.c           Main simulation case
│   └── DataFiles/                 Input geometry data
├── runSimulation.sh               Single case runner
├── runParameterSweep.sh           Parameter sweep runner (local)
├── runSweepHamilton-serial.sbatch HPC Stage 1 runner (Durham Hamilton)
├── runSweepHamilton.sbatch        HPC sweep runner (Durham Hamilton)
├── runSweepSnellius-serial.sbatch HPC Stage 1 runner (SURF Snellius)
├── runSweepSnellius.sbatch        HPC sweep runner (SURF Snellius)
├── runPostProcess-Ncases.sh       Post-processing pipeline
├── default.params                 Single-case configuration
├── sweep.params                   Sweep configuration
```

## Key Parameters

- **Ohnesorge Number (Oh)**: `Oh = mu/sqrt(rho*sigma*R)` - ratio of viscous to inertial-capillary forces
- **Bond Number (Bo)**: `Bo = rho*g*R^2/sigma` - ratio of gravitational to surface tension forces
- **Maximum Refinement Level**: Controls mesh resolution (e.g., level 10 = 1024 cells)
- **tmax**: Maximum simulation time (dimensionless, based on capillary time scale)

## Requirements

- **Basilisk Framework**: Install via the ref-locked script above (upstream docs: [basilisk.fr](https://basilisk.fr))
- **MPI** (optional): For parallel execution
  - macOS: `brew install open-mpi`
  - Linux: `sudo apt-get install libopenmpi-dev`

## Two-Stage Execution

The simulation uses a two-stage execution model due to a Basilisk limitation (`distance.h` is incompatible with MPI):

1. **Stage 1**: Generate restart file (serial or OpenMP)
2. **Stage 2**: Run full simulation from restart (supports MPI)

```bash
# Run both stages (default)
./runSimulation.sh default.params

# Or separately:
./runSimulation.sh --stage1 default.params    # Generate restart
./runSimulation.sh --stage2 --mpi 8 default.params  # Full simulation
```

## Troubleshooting

### "restart file not found"

Stage 2 requires a restart file from Stage 1. Run Stage 1 first:

```bash
./runSimulation.sh --stage1 default.params
```

### "restart file is empty"

Stage 1 may have failed silently. Check for:
- Compilation errors in the case directory
- Memory issues (reduce MAXlevel)
- Invalid parameters

### "qcc not found"

Basilisk is not installed. Run the install script:

```bash
curl -sL https://raw.githubusercontent.com/comphy-lab/basilisk-C/main/reset_install_basilisk-ref-locked.sh | bash -s -- --ref=v2026-01-13
```

### Parameters not updating on reruns

The scripts preserve existing `case.params` and source files for reruns. Use `--force` to overwrite:

```bash
./runSimulation.sh --force default.params
```

### HPC jobs failing

1. Ensure Stage 1 completed locally before submitting Stage 2
2. Check that restart files exist in each case directory
3. Verify SLURM parameters match your allocation

## License

Copyright (C) 2026 CoMPhy Lab.

This project is licensed under the GNU General Public License v3.0. See [LICENSE](LICENSE) for details.

## Contact

For questions or collaboration inquiries, please contact the [CoMPhy Lab](https://comphy-lab.org).
