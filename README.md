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
│   ├── parse_params.sh            Parameter parsing utilities (shell layer)
│   ├── sweep_utils.sh             Sweep generation utilities
│   ├── basilisk_version.sh        Centralized version pinning
│   ├── params.h                   C-side runtime parameter layer (struct + file parser + CLI overrides + validation)
│   ├── two-phaseVE.h              Two-phase VE properties (vendored from MultiRheoFlow)
│   ├── log-conform-viscoelastic-scalar-2D.h  Axisymmetric Oldroyd-B log-conformation
│   └── log-conform-viscoelastic.h Tensor log-conformation (optional)
├── postProcess/                   Post-processing tools and visualization
│   ├── getData.c                  Field extraction on structured grids
│   ├── getFacet.c                 Interface geometry extraction
│   ├── Video.py                   Frame-by-frame visualization pipeline
│   ├── plot_comphy_fields.py      Split-axi stills (dissipation or tr A, |u|)
│   └── plot_comphy_video.py       Same layout as a mathtext MP4
├── simulationCases/               Case-based simulation outputs
│   ├── burstingBubble.c           Newtonian, fixed-ceiling wavelet AMR
│   ├── burstingBubble-drillResolution.c     Newtonian feature-tracking drill
│   ├── burstingBubbleVE.c         Oldroyd-B, usual wavelet AMR
│   ├── burstingBubbleVE-drillResolution.c   Oldroyd-B + drill
│   ├── DataFiles/                 Input geometry data (BoXXXX.dat)
│   └── initialConditions/         Young-Laplace and Bo=0 shape generators
├── runSimulation.sh               Single case runner
├── runParameterSweep.sh           Parameter sweep runner (local)
├── runSweepHamilton-serial.sbatch HPC Stage 1 runner (Durham Hamilton)
├── runSweepHamilton.sbatch        HPC sweep runner (Durham Hamilton)
├── runSweepSnellius-serial.sbatch HPC Stage 1 runner (SURF Snellius)
├── runSweepSnellius.sbatch        HPC sweep runner (SURF Snellius)
├── runPostProcess-Ncases.sh       Post-processing pipeline
├── default.params                 Newtonian single-case configuration
├── default-ve.params              Oldroyd-B usual-AMR configuration
├── default-ve-drill.params        Oldroyd-B drill configuration
├── sweep.params                   Sweep configuration
```

## Key Parameters

All knobs live in `default.params` (single case) or per-case `case.params` files and are
read at runtime by the simulation; you never edit the source to change a run. The main
physical parameters are:

- **Ohnesorge Number (Oh)**: `Oh = mu_s/sqrt(rho*sigma*R)` - solvent viscous to inertial-capillary forces (`OhRatio` sets the gas-phase value, `Oha = OhRatio*Oh`)
- **Bond Number (Bo)**: `Bo = rho*g*R^2/sigma` - ratio of gravitational to surface tension forces
- **Deborah Number (De)**: liquid relaxation time over the capillary time. `0` is Newtonian. Used only by the VE solvers.
- **Elasto-capillary Number (Ec)**: liquid elastic modulus over capillary stress. `0` is Newtonian. Polymeric viscosity is the derived product `Oh_p = Ec*De`.
- **tmax**: Maximum simulation time (dimensionless, based on the capillary time scale)
- **zWall**: Distance from the bubble south pole to the bottom wall (sets the domain size)
- **Solver**: which `simulationCases/*.c` entry point to compile. `burstingBubble` is the Newtonian baseline; `burstingBubbleVE` is usual-AMR Oldroyd-B; `burstingBubbleVE-drillResolution` is Oldroyd-B with the singularity drill.

### Adaptive resolution

Both the mesh and the timestep are adaptive, and every control is a tunable parameter:

- **Space** — `MAXlevel` (finest level, e.g. 12 = 4096 cells), `MINlevel` (far-field
  coarsening floor; the interface is always kept at `MAXlevel`), `init_grid_level`
  (initial uniform grid), and the wavelet error tolerances `fErr`, `VelErr`, `KErr`.
- **Time** — `CFL`, `TOLERANCE` (Poisson/viscous solver), and `dtmax`. **`dtmax` is a
  *ceiling*, not a fixed step.** Surface tension is time-explicit, so `tension.h` reduces
  the real step to the capillary-wave limit
  `T = sqrt(rho_m * Delta_min^3 / (pi * sigma))` every iteration. The effective timestep is
  therefore genuinely adaptive and scales with the finest cell size — coarser runs
  automatically take larger steps. (Earlier versions hard-capped `dtmax` at `1e-5`, below
  the capillary limit, which fixed the step and throttled the run.)

To change resolution, edit `default.params`; to study its effect, sweep it (e.g.
`SWEEP_MAXlevel=10,11,12` in `sweep.params`). A legacy positional CLI
(`./burstingBubble <MAXlevel> <Oh> <Bond> <tmax> <zWall>`) is still accepted as a fallback,
but the parameter-file interface is preferred.

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
