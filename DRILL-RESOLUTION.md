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
run. The post-processing reads the snapshot **dumps** (`getJetFoot`/`getFacet`
via `restore()`), not the solver `log`, and the drill dumps the same `f`,`u`
fields — so the pipeline runs on drill snapshots unchanged. (The drill `log`
does add a `maxlevel` column, `i dt t ke maxlevel r_b z_b`, for its own
diagnostics; nothing consumes it downstream.)

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

`log` columns are now `i dt t ke maxlevel r_b z_b r_base z_base q_jet q_l`
(the reference solver had no `maxlevel` column). `maxlevel` is
`maxlevelLocal`, the live ceiling — plot it against `t` to see the drill
work. The probe and level are computed once per step in `drillProbe(i++)`
and reused in `logWriting`, so curvature/tagging is evaluated as in the
reference solver `burstingBubble.c`.

Two base-point columns coexist deliberately:

- `r_b z_b` — the AMR probe (max-|kappa| cavity focus pre-inception, lowest
  MainPhase interfacial point post-inception). This drives the drill; it can
  latch onto satellites late in the jet phase, which is fine for refinement
  (it over-refines the thin column) but wrong as an observable.
- `r_base z_base` — the robust outer-free-surface base (inlined
  `postProcess/getBase.c` logic: MainLiq + MainGas double-tag, satellite- and
  droplet-proof). `q_jet = INT_0^{r_base} u_z r dr` and
  `q_l = INT_0^{r_base} u_z dr` are the single-plane base fluxes at `z_base`
  (getJetFoot.c definitions, no 2*pi, no f-weighting). These are the
  on-the-fly science observables: `q_jet(r_jet)` / `q_l(r_jet)` and
  `R_j x Q_L` come straight off the log.

Time is logged to 8 decimals, observables to 6-decimal scientific — dense
staged-cadence data near inception is not precision-starved.

## Snapshot names

Snapshots are dumped as `intermediate/snapshot-%8.6f` (6 decimals; was 4).
The staged `tsnap` can drop below `1e-4` near inception, where 4-decimal
names would collide and silently overwrite dumps. The postProcess tools
discover snapshots by glob and parse `t` from the filename, so both old
(4-dp) and new (6-dp) cases remain readable.

## Build

Same as the other solvers (see `runSimulation.sh`), just a different source
file. From a case directory (`simulationCases/<CaseNo>/`):

```sh
# serial / OpenMP
qcc -O2 -Wall -disable-dimensions -fopenmp -I../../src-local \
    burstingBubble-drillResolution.c -o drill -lm

# MPI — build WITHOUT -D_GNU_SOURCE (see the MPI note below): this leaves
# Basilisk's FPE trap off, which the drill's dynamic coarsening needs. The
# ke blow-up/decay checks in logWriting remain the guard for a real divergence.
# -D_DEFAULT_SOURCE exposes madvise()/MADV_DONTNEED on glibc (a bare
# -std=c99 hides them and the build fails) without re-arming the trap.
CC99='mpicc -std=c99 -D_DEFAULT_SOURCE' qcc -Wall -O2 -D_MPI=1 \
    -disable-dimensions -I../../src-local \
    burstingBubble-drillResolution.c -o drill -lm

./drill case.params                 # serial
mpirun -np 4 ./drill case.params    # MPI
```

`distance.h` (fresh-init from `DataFiles/Bo*.dat`) is incompatible with MPI, so
a from-scratch **first** run is serial; MPI restarts from the dump.

> **MPI build flag.** Do NOT add `-D_GNU_SOURCE` to the MPI build. On Linux
> that flag turns on Basilisk's floating-point trap, which fires *spuriously*
> on the drill's aggressive coarsen/refine/rebalance (see the MPI note in
> Validation §5) and aborts an otherwise-correct run with SIGFPE. Without it
> the trap is off and the run is stable and bit-for-bit reproducible across
> rank counts.

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

### 5. MPI status — resolved (spurious FP trap)

An MPI snapshot restart originally raised a `SIGFPE` mid-run when the ceiling
had coarsened to a low level. It is a **spurious floating-point trap, not a
numerical failure.** Basilisk's trap (enabled by `-D_GNU_SOURCE` on Linux)
fires on its own `undefined` NaN-sentinel in a transient ghost cell left by the
coarsen/refine/rebalance; the value never enters the physics.

Root-caused by controls:

| Run | Build | Result |
|---|---|---|
| drill `start=8`, MPI np=4, trap **on** | `-D_GNU_SOURCE` | SIGFPE at t≈0.337 (10→11 refine) |
| drill `start=8`, MPI np=2, trap **off** | no `-D_GNU_SOURCE` | clean past crash pt; `ke=4.7729` |
| drill `start=8`, MPI np=4, trap **off** | no `-D_GNU_SOURCE` | clean; `ke=4.7729` — **bit-identical to np=2** |
| drill `start=8`, MPI np=4, trap **off**, full sweep | no `-D_GNU_SOURCE` | ramps 10→11→12 through the focus, ke peak ≈ ref, runs to `tmax` |
| fixed-level reference, MPI np=4 | `-D_GNU_SOURCE` | stable (never coarsens → never leaves an `undefined` cell) |
| drill `start=8`, serial / OpenMP | (OpenMP doesn't trap) | stable, full sweep |

`np=2` and `np=4` giving **bit-identical** `ke` is the clincher: real MPI data
corruption would make the two decompositions diverge. They don't — the physics
is correct; only the trap was firing.

**Fix:** build the MPI binary **without** `-D_GNU_SOURCE` (see Build). The
`ke` blow-up/decay checks in `logWriting` remain the guard for a genuine
divergence. This reclaims the deep-coarsening savings under multi-rank MPI — no
`drillMaxlevelStart` workaround needed. Deep coarsening (down to level 8) is now
validated serial, OpenMP, and MPI.

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

## Visualizing the mesh: `getView2D`

`postProcess/getView2D.c` renders one snapshot's interface + adaptive mesh
(`cells()`) to a PNG, mirrored across the r=0 axis — a 2D analogue of
`comphy-lab/Jumping-Drops`' `postProcess/getView3D_v2.c`, using the same
`view.h`/`draw.h` (bview) machinery but without the 3D camera angles/multi-axis
mirrors (axisymmetric only needs one mirror). This is the tool to *see* the
drill working — the AMR "onion layers" coarsening away from the tracked
feature.

Compilation is different from the plain solvers — link against Basilisk's
headless software-rendering framebuffer (`fb_tiny`), no display needed:

```sh
qcc -O2 -Wall -disable-dimensions postProcess/getView2D.c -o getView2D \
    -L$BASILISK/gl -lglutils -lfb_tiny -lm
```

Usage: `./getView2D <snapshot> <output.png> [fov tx ty width height]`.
`fov` (degrees) controls zoom — smaller is tighter. `tx`, `ty` pan the camera;
to center on a physical point `(z0, r0)` (e.g. the current jet-base/cavity-focus
probe from the solver `log`), set `tx = -z0/L0`, `ty = -r0/L0` (Basilisk scales
the scene by `1/L0` before translating). A `fov` of 4–6 with these
auto-centering formulas gives a tight, well-framed zoom on the singularity;
`fov=24` (default) frames roughly the whole box.

**Rotate the output 90° CCW.** Basilisk's native (x,y) are (axial, radial),
unrotated, so the raw render has axial horizontal / radial vertical — sideways
relative to the r-horizontal/z-vertical, jet-points-up convention every other
figure in this project uses. Post-process with:

```sh
magick getView2D_output.png -rotate -90 final.png
```

(ImageMagick's `-rotate` is clockwise for positive degrees, so `-90` is the
needed 90° CCW turn.) Verified against a late-time frame with a fully-formed
jet and pinched tip droplet: after rotation the free surface is horizontal,
the crater dips below it, and the jet rises correctly upward through it —
matching every other rendering in this project.

## Files

Drill (this work):
- `simulationCases/burstingBubble-drillResolution.c` — the drill solver.
- `src-local/params.h` — 7 new `drill*` knobs (struct, defaults, parse,
  validate, print).
- `default.params` — documented `drill*` block.
- `postProcess/getView2D.c` — mesh + interface PNG renderer (see above).
- `DRILL-RESOLUTION.md` — this document.

Probe + post-processing (bundled from the jet-base-tracking line):
- `simulationCases/burstingBubble-adaptiveResolution.c` — logging-only probe
  solver (fixed ceiling); the A/B reference and the origin of the inlined probe.
- `postProcess/getJetFoot.c` — base-flux `q_jet`/`q_l` and jet-foot geometry.
- `postProcess/VideoFoot.py`, `footplots.py`, `conefit.py` — overlay video,
  PRL figures, cone fit.

A/B baseline is `simulationCases/burstingBubble.c` (fixed level) or this solver
with `drillAMR=0`.
