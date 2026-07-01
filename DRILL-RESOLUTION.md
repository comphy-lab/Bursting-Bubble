# Drill adaptive resolution — mesh + time

`simulationCases/burstingBubble-drillResolution.c`

This is the "drill": a feature-tracking adaptive-mesh **and** adaptive-time
trigger for the singular bursting-bubble jet. It drives the local refinement
ceiling from the jet-base / cavity-focus probe so the grid drills resolution
into the collapsing region as the singularity approaches, and relaxes it
everywhere else. The refinement idea is taken from
[`comphy-lab/ElasticPinchOff`](https://github.com/comphy-lab/ElasticPinchOff)
(`LiquidOutThinning.c`) and generalised from a single neck-radius scalar to the
two-regime bursting-bubble probe.

The drill solver is **self-contained**: the jet-base / cavity-focus probe is
inlined (the same detection algorithm developed as a logging-only diagnostic in
`burstingBubble-adaptiveResolution.c`), and here it *drives* the mesh rather
than just being logged. The A/B baseline is `burstingBubble.c` (fixed level) or
this same solver with `drillAMR=0`.

This PR also bundles the logging-only probe solver
(`burstingBubble-adaptiveResolution.c`) and the flux post-processing
(`getJetFoot.c`, `VideoFoot.py`, `footplots.py`, `conefit.py`), so the base
flux `q_jet`/`q_l`, cone-fit, and `R_j × Q_L` observable can be computed from a
run. Note the drill's `log` adds a `maxlevel` column (`i dt t ke maxlevel r_b
z_b`) vs the logging solver's `i dt t ke r_b z_b`; the post-processing parsers
were written for the latter, so drill logs need a one-column offset (a small
follow-up).

## Why the fixed-level run wastes work

With a fixed `MAXlevel = 12`, the wavelet `fErr` criterion pins the **entire**
interface near level 12 for the whole run — the slow pre-collapse buildup, the
fast focus/jet/pinch event, and the long post-pinch relaxation all pay
level-12 cost even though only the focusing region needs it. In case 1000 the
scientifically interesting dynamics (cavity focus → jet emergence → first
pinch) sit in roughly `t ∈ [0.40, 0.55]` — under 10 % of an ~11 h / 4-rank run.

## The drill (two regimes)

A **local** ceiling `maxlevelLocal ≤ MAXlevel` is handed to `adapt_wavelet`
each step instead of the fixed `MAXlevel`. It is set from the resolved length of
the active feature:

| Phase | Tracked length `s` | Meaning |
|---|---|---|
| pre-inception (cavity-focus collapse) | `s = 1/|κ|_max` | radius of curvature at the focusing cavity base; diverging curvature ⇒ shrinking `s` |
| post-inception (jet growth) | `s = r_b` | jet-base radius; on the slender/constant-flux branch it keeps shrinking as the jet lengthens, so demand keeps climbing |

The two regimes are separated by the probe's forward-in-time **inception latch**
`jetFormed` (set when the max-curvature point has collapsed onto the axis while
the lowest interfacial point has lifted off it).

**Demanded level** = smallest `L ∈ [drillMaxlevelStart, MAXlevel]` such that `s`
spans at least `nCells` cells: `s ≥ nCells · Ldomain / 2^L`. `nCells` differs by
regime (`drillNcellsK`, `drillNcellsJet`) because a curving front and a slender
column are geometrically different features.

**Hysteresis:** refine immediately to the demanded level (never under-resolve a
growing feature), coarsen at most one level per step (kills refine/coarsen
churn). Because `adapt_wavelet` itself moves the grid by at most one level per
call, the grid ramps smoothly regardless.

**Terminal relaxation (off by default):** once the first tip droplet sheds
(`jetFormed && n_components > 1`), if `drillRelaxLevel > 0` the ceiling relaxes
to that level so the long post-pinch tail runs cheaply. Disabled by default —
the target observable `R_j × Q_L` is measured at the receding base, which keeps
demanding resolution.

## Adaptive time

Two coupled effects:

1. **Solver timestep (automatic).** Surface tension is time-explicit, so
   `tension.h` limits the step to the capillary-wave period
   `dt ~ sqrt(ρ_m · Δ_min³ / (π σ))`. As `maxlevelLocal` ramps and `Δ_min`
   shrinks the timestep tightens on its own — the mesh drill drives the time
   drill for free. `dtmax` is only a ceiling. Nothing to configure.
2. **Output cadence (staged).** `drillTsnapStages=1` tightens the snapshot
   interval as the mesh refines, `tsnap = base · 2^(start − level)` floored at
   `drillTsnapMinFactor · base`, so the fast jet/pinch window is sampled densely.

## Why not just port ElasticPinchOff

Pinch-off tracks **one** scalar whose global minimum *is* the neck
(`statsf(Y).min` is safe by construction), has **one** singular event, and never
resets `broken`. Bursting-bubble is harder on every axis:

- multiple interface features (flat free surface, wall meniscus, cavity wall,
  shed drop) ⇒ needs the `MainPhase`-tagged, near-axis-restricted probe, not a
  global `statsf`;
- (at least) **two** singular events in sequence governed by **different**
  scalars (curvature-radius, then jet-base radius);
- the jet base keeps demanding resolution *after* inception rather than
  stopping, so a monotone up-only ratchet is wrong — the ceiling must be able to
  coarsen again, hence the refine-fast / coarsen-slow rule.

## Knobs

All runtime (`case.params`), parsed by `src-local/params.h`. Consumed only by
the drill solver; the plain solvers parse and ignore them.

| Key | Default | Meaning |
|---|---|---|
| `drillAMR` | `1` | master switch; `0` pins `maxlevelLocal = MAXlevel` (reproduces the fixed-level reference bit-for-bit — use for A/B) |
| `drillMaxlevelStart` | `8` | ramp floor (coarsest level the tracked band may reach) |
| `drillNcellsK` | `5` | cells per curvature-radius, pre-inception |
| `drillNcellsJet` | `5` | cells per jet-base radius, post-inception |
| `drillRelaxLevel` | `-1` | relax level after first tip pinch; `≤0` = hold resolution (safe) |
| `drillTsnapStages` | `1` | stage snapshot cadence with the mesh |
| `drillTsnapMinFactor` | `0.1` | floor on the staged `tsnap` as a fraction of base |

## Log format

`log` columns are now `i dt t ke maxlevel r_b z_b` (the reference solver had no
`maxlevel` column). `maxlevel` is `maxlevelLocal`, the live ceiling — plot it
against `t` to see the drill work. The probe and level are computed once per
step in `drillAdapt(i++)` and reused in `logWriting`, so curvature/tagging is
evaluated as in the reference solver `burstingBubble.c`.

## Build

Same as the other solvers (see `runSimulation.sh`), just a different source
file. From a case directory (`simulationCases/<CaseNo>/`):

```sh
# serial / OpenMP
qcc -O2 -Wall -disable-dimensions -fopenmp -I../../src-local \
    burstingBubble-drillResolution.c -o drill -lm

# MPI (note the -lm at the end and the glibc guard)
CC99='mpicc -std=c99 -D_GNU_SOURCE=1' qcc -Wall -O2 -D_MPI=1 \
    -disable-dimensions -I../../src-local \
    burstingBubble-drillResolution.c -o drill -lm

./drill case.params                 # serial
mpirun -np 4 ./drill case.params    # MPI
```

`distance.h` (fresh-init from `DataFiles/Bo*.dat`) is incompatible with MPI, so
a from-scratch **first** run is serial; MPI restarts from the dump.

## Recommended workflow (two-stage)

The initial condition has a thin liquid **film** whose retraction is already
sharp, so a naive from-scratch drill with a low `drillMaxlevelStart` would
under-resolve it. Mirror the existing two-stage pattern:

1. **Stage 1 — film rupture, fixed level.** Serial, `drillAMR=0`, to a restart
   point past film rupture (e.g. `tmax=0.10`). This resolves the film at full
   `MAXlevel`.
2. **Stage 2 — collapse/jet/pinch, drill on.** MPI restart from the Stage-1
   dump, `drillAMR=1`. The drill coarsens where the cavity is smooth and ramps
   to `MAXlevel` through the focus.

Alternatively, restart the drill from any existing case-1000 snapshot
(`intermediate/snapshot-*`) — dumps carry only `f` and `u`, so they are
solver-agnostic.

## Validation

Environment: `machine-ts`, Basilisk `v2026-01-13`, case 1000
(`Oh=0.01, Oha=2e-4, Bo=0.001, MAXlevel=12, zWall=4, Ldomain=10`). Reference =
the fixed-level-12 run at `.../2026-06-21-Singular-Bursting-Bubbles/simulationCases/1000`.

### 1. Compiles and runs

Serial and MPI builds both clean (only Basilisk's own header warnings). Runs
from scratch (serial, `DataFiles` init) and from a snapshot restart (MPI).

### 2. The trigger discriminates

From scratch with `MAXlevel=12`, the drill holds `maxlevelLocal = 10` through
the whole pre-focus window `t ∈ [0.05, 0.30]`: the sharpest cavity feature has
`1/|κ|_max ≈ 0.05`, which is ≈ 5 cells at level 10 (`Δ₁₀ = 10/1024 ≈ 0.0098`).
So it runs **two levels below** the reference's fixed 12 wherever the criterion
allows, and ramps toward 12 only as the cavity focuses.

### 3. Coarsening does not corrupt the physics

Drill (running at level 10) vs reference (fixed level 12), kinetic energy over
`t ∈ [0.02, 0.30]` (890 points):

| metric | value |
|---|---|
| mean relative `ke` difference | **1.41 %** |
| median | 1.08 % |
| max | 4.81 % (startup transient only) |

| `t` | `ke` drill (L10) | `ke` ref (L12) | rel |
|---|---|---|---|
| 0.05 | 0.4869 | 0.4999 | 2.61 % |
| 0.10 | 1.0855 | 1.0990 | 1.23 % |
| 0.20 | 2.4907 | 2.5171 | 1.05 % |
| 0.30 | 4.0512 | 4.0954 | 1.08 % |

A ~1 % energy match while running two levels coarser through the entire
pre-focus phase.

### 4. Full regime sweep (serial restart, t=0.29 → 0.60)

Restart from `snapshot-0.2900`, `drillMaxlevelStart=8`, `drillNcells=5`,
`drillRelaxLevel=9`, run to completion (tmax=0.60). Every code path fired, in
order, and the run completed cleanly:

| Regime | Observed |
|---|---|
| coarsen on restart | 12 → 8 over the first steps (feature broad) |
| ramp-up through focus | 8 → 9 → 10 → 11 → 12 as `1/κ` shrank |
| level 12 held | t ≈ 0.34 → 0.44 (the focus / jet-emergence window) |
| `jetFormed` latch | i = 7536; probe switched to jet-base radius |
| tip-pinch + relaxation | ceiling relaxed 12 → 9 from t ≈ 0.44, held 9 to end |
| completed | ran to tmax = 0.60, no crash |

Crude interface-band cost proxy (Σ 2^level over steps): **~62 %** of the
fixed-level-12 run.

Caveat — restart transient: the dump carries only `f` and `u` (no pressure
history), so a snapshot restart re-solves pressure from scratch. Combined with
the aggressive `start=8` coarsening, the drill's kinetic-energy peak lands at
t ≈ 0.404 (magnitude 6.43) vs the continuous reference's 6.33 at t ≈ 0.461 — a
timing shift. So the *serial restart* validates the mechanism (all regimes
fire, stable to completion) but is NOT a clean quantitative benchmark. The
clean quantitative check is the from-scratch pre-focus comparison in §3 (~1 %),
which has no restart confound. A from-scratch drill run (Stage 1 fixed-level →
Stage 2 drill) is the correct A/B for collapse timing; queued as follow-up.

### 5. MPI status — known issue (localised)

Serial and OpenMP builds are stable end-to-end (above). Under **MPI**, a
snapshot restart raises an FPE mid-run. It has been localised by three
controls:

| Run | Mode | Coarsens to | Outcome |
|---|---|---|---|
| fixed-level reference | MPI, 4 ranks | — (stays 12) | stable past the focus |
| drill `start=8` | MPI, 4 ranks | 8–10 | FPE at t≈0.337 (grid at level 10) |
| drill `start=11` | MPI, 2 ranks | 11 only | stable past t=0.345 (no FPE) |
| drill `start=8` | serial / OpenMP | 8 | stable, full sweep to tmax |

So the FPE is **not** physics (serial/OpenMP coarsen to level 8 and complete),
**not** the dynamic mechanism (MPI `start=11` is fine), and **not** the
snapshot (the fixed-level reference restarts under MPI cleanly). It appears when
the ceiling drops to a **low level (≲10) under MPI** — i.e. MPI rebalancing of
a coarsened adaptive grid, a Basilisk-infrastructure corner rather than a
solver-logic bug.

**Workaround (immediate):** for multi-rank MPI, set `drillMaxlevelStart` high
(10–11). This keeps the ramp-to-`MAXlevel` at the focus and the relaxation, but
sacrifices the deep pre-focus coarsening — so the large savings regime
(coarsen to 8) is currently serial/OpenMP-only.

**Recommended path:** run single-node serial/OpenMP, which is validated
end-to-end and sufficient for case-1000-class runs (the reference was ~11 h on
4 MPI ranks; comparable single-node on 8 OpenMP threads). Root-causing the
low-level MPI rebalance FPE (compile without FP-trapping, locate the first
NaN; test rank counts) is the follow-up for multi-node scaling. Track:
`memory/projects/singular-bursting-bubbles.md`.

## Tuning

- **`drillNcellsK` / `drillNcellsJet`** are the accuracy dial. `5` gave ~1 %
  on `ke` pre-focus; raise toward `8–12` for a tighter match at higher cost.
  Split them if the jet regime needs a different value than the curvature
  regime — the fastest calibration is to replay a finished case's `1/κ`, `r_b`
  trajectories (the `log` already carries `r_b`; `1/κ` is recoverable from the
  snapshots) through candidate `nCells` offline before spending compute.
- **`drillMaxlevelStart`** sets pre-collapse savings vs film-resolution safety;
  keep it high enough to resolve the initial film if running from scratch, or
  use the two-stage workflow.
- **`drillRelaxLevel`** reclaims the post-pinch tail. Leave disabled until the
  base measurement (`R_j × Q_L`) is confirmed complete, then set to ~9–10.

## Files

Drill (this work):
- `simulationCases/burstingBubble-drillResolution.c` — the drill solver.
- `src-local/params.h` — 7 new `drill*` knobs (struct, defaults, parse,
  validate, print).
- `default.params` — documented `drill*` block.
- `DRILL-RESOLUTION.md` — this document.

Probe + post-processing (bundled from the jet-base-tracking line):
- `simulationCases/burstingBubble-adaptiveResolution.c` — logging-only probe
  solver (fixed ceiling); the A/B reference and the origin of the inlined probe.
- `postProcess/getJetFoot.c` — base-flux `q_jet`/`q_l` and jet-foot geometry.
- `postProcess/VideoFoot.py`, `footplots.py`, `conefit.py` — overlay video,
  PRL figures, cone fit.

A/B baseline is `simulationCases/burstingBubble.c` (fixed level) or this solver
with `drillAMR=0`.
