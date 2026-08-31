/**
# getDropStats.c — per-snapshot statistics of the main body and every shed drop

The solvers log kinetic energy, the drill probe, the jet base and the base
fluxes. None of them records anything about an emitted drop, and neither did
any tool in `postProcess/`. This one does: for a single snapshot it emits the
volume, equivalent radius, surface area, centroid, extent, centre-of-mass
velocity and kinetic energy of the main liquid body and of every detached
liquid fragment, plus the jet tip and its velocity.

That is what the first-drop radius `R_d` and velocity `V_d` are computed from,
and what the emitted-drop count, total emitted surface, total emitted volume
and total emitted kinetic energy are summed from. Drop identity across
snapshots, emission ordering (`n = 1, 2, ...`) and the settled-volume
measurement rule live in the Python driver, not here — this tool is a pure
per-snapshot reduction with no memory.

## Component definition, and why it is not `getBase.c`'s

`getBase.c` tags PURE-liquid cells (`f > 1 - 1e-4`) because it needs the main
body's outer free surface. That threshold is wrong for enumerating drops: a
late satellite two or three cells across has no pure-liquid core at all and
would simply vanish. Here components are tagged on `f > 0.5`, which every
resolved drop has, and are then dilated by one cell into the surrounding
interfacial rim so the volume integral is not truncated at the `f = 0.5`
contour. The main body is the largest component by volume.

The `f > 0.5` tagging and the pure-liquid tagging agree on the main body,
which is checked and reported (`MAIN` row, `nliq` column) rather than assumed.

## Output (stderr, one line per row, whitespace separated)

```text
MAIN t nliq ndrop V Rv S Rs zc zmin zmax rmax vz vr Ek ztip rtip vtip Vtot dmin cells
DROP t id     V Rv S Rs zc zmin zmax rmax vz vr Ek dmin cells
```

- `V`   volume, axisymmetric `sum 2 pi y f Delta^2`
- `Rv`  volume-equivalent radius, `(3V/4pi)^(1/3)` — this is the reported `R_d`
- `S`   interfacial area, axisymmetric `sum 2 pi (y + p_y Delta) len Delta`
- `Rs`  area-equivalent radius, `sqrt(S/4pi)`; `Rs/Rv` is a free sphericity
        check — a value well above 1 means the fragment is a ligament, not a
        drop, and must not be reported as a drop radius
- `zc`, `vz`, `vr` mass-weighted centroid and centre-of-mass velocity
- `Ek`  kinetic energy, `sum 2 pi y (rho/2)(u.x^2 + u.y^2) Delta^2`
- `ztip`, `rtip`, `vtip` highest near-axis interfacial point of the MAIN body
        and the axial velocity there (the jet tip, distinct from any drop)
- `Vtot` total liquid volume in the domain; `Vtot - V(MAIN) - sum V(DROP)` is
        the mass left unassigned by the tagging, and should be a rounding error
- `dmin`, `cells` finest cell resolving the body, and `Rv/dmin` — the number of
        cells across the drop's radius. **This is not optional bookkeeping.** A
        drop only a few cells across is reported at whatever radius the mesh can
        represent, not its physical one, and sphericity does not catch it: a
        well-resolved sphere and a four-cell blob both score ~0.99. Any drop
        below roughly eight cells per radius is a mesh artefact until a
        refinement study says otherwise.

Lengths are in `R_0`, velocities in the inertio-capillary `V_c`, energies in
`rho V_c^2 R_0^3`. Convert to the experimental viscous-capillary velocity with
`V/V_mu = Oh * (V/V_c)`.

## Serial only

Like `getBase.c`, run this SERIALLY, one process per snapshot. The size
tallies are reduced so the logic is MPI-safe, but the bare tool is not meant
to be launched under `mpirun`. Dumps carry only `f` and `u`, so the same
binary reads Newtonian and viscoelastic snapshots.

@author CoMPhy Lab
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "fractions.h"
#include "tag.h"

/* Snapshot paths routinely exceed 80 bytes once a campaign run root is
   included, and the historical `char filename[80]` in the other postProcess
   tools overflows on those. Size for a real path and truncate safely. */
char filename[4096];

#define FDROP  0.5       // component threshold: every resolved drop exceeds it
#define FEPS   1e-6      // interfacial band
#define FPURE  1e-4      // getBase.c's pure-phase band, for the cross-check
#define R_TIP  0.25      // near-axis band defining the jet tip

int main (int a, char const *arguments[]) {
  if (a < 2) {
    fprintf (ferr, "usage: %s <snapshot>\n", arguments[0]);
    return 1;
  }
  if (snprintf (filename, sizeof(filename), "%s", arguments[1]) >= (int) sizeof(filename)) {
    fprintf (ferr, "ERROR: snapshot path longer than %zu bytes\n", sizeof(filename));
    return 1;
  }
  restore (file = filename);
#if TREE
  f.prolongation = fraction_refine;
#endif
  rho1 = 1., rho2 = 1e-3;          // match the solver, so Ek is comparable
  boundary ((scalar *){f});

  /**
  ## Components

  Tag on `f > FDROP`, then dilate one cell into the interfacial rim so that
  the volume integral captures the whole fragment rather than the part above
  the `f = 0.5` contour. */
  scalar d[];
  foreach() d[] = (f[] > FDROP);
  int n = tag (d);
  if (n < 1) {
    fprintf (ferr, "MAIN %.8f 0 0 0 0 0 0 0 0 0 0 0 0 0 -1000 -1000 -1000 0 0 -1\n", t);
    fflush (ferr);
    return 0;
  }
  boundary ((scalar *){d});

  scalar dm[];
  foreach() dm[] = d[];
  foreach() {
    if (d[] == 0 && f[] > FEPS) {
      int tg = 0;
      foreach_dimension() {
        if (!tg && d[1] > 0)  tg = (int) d[1];
        if (!tg && d[-1] > 0) tg = (int) d[-1];
      }
      if (tg) dm[] = tg;
    }
  }
  boundary ((scalar *){dm});

  /**
  ## Per-component reductions

  `foreach(serial)` for the same reason as `getBase.c`: accumulating into a
  per-component array from a threaded loop is a race. */
  double *cV = calloc (n, sizeof(double)), *cS  = calloc (n, sizeof(double));
  double *cZ = calloc (n, sizeof(double)), *cUZ = calloc (n, sizeof(double));
  double *cUR = calloc (n, sizeof(double)), *cEK = calloc (n, sizeof(double));
  double *cZMIN = calloc (n, sizeof(double)), *cZMAX = calloc (n, sizeof(double));
  double *cRMAX = calloc (n, sizeof(double));
  double *cDMIN = calloc (n, sizeof(double));
  for (int j = 0; j < n; j++) {
    cZMIN[j] = HUGE; cZMAX[j] = -HUGE; cRMAX[j] = 0.; cDMIN[j] = HUGE;
  }

  foreach(serial) {
    int j = (int) dm[] - 1;
    if (j < 0 || j >= n) continue;
    double dv = 2.*pi*y*f[]*sq(Delta);              // axisymmetric liquid volume
    cV[j]  += dv;
    cZ[j]  += dv*x;
    cUZ[j] += dv*u.x[];
    cUR[j] += dv*u.y[];
    cEK[j] += 2.*pi*y*0.5*rho(f[])*(sq(u.x[]) + sq(u.y[]))*sq(Delta);
    if (x < cZMIN[j]) cZMIN[j] = x;
    if (x > cZMAX[j]) cZMAX[j] = x;
    if (f[] > FDROP && y > cRMAX[j]) cRMAX[j] = y;
    if (Delta < cDMIN[j]) cDMIN[j] = Delta;      // finest cell resolving this body
    if (f[] > FEPS && f[] < 1. - FEPS) {            // interfacial area
      coord m = mycs (point, f);
      double alpha = plane_alpha (f[], m);
      coord p;
      double len = plane_area_center (m, alpha, &p);
      cS[j] += 2.*pi*(y + p.y*Delta)*len*Delta;
    }
  }
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, cV,  n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, cS,  n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, cZ,  n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, cUZ, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, cUR, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, cEK, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, cZMIN, n, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, cZMAX, n, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, cRMAX, n, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (MPI_IN_PLACE, cDMIN, n, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
#endif

  /**
  Total liquid volume over the whole domain. Reported so that mass closure —
  `Vtot` against `V(MAIN) + sum V(DROP)` — is checkable on every snapshot
  rather than trusted. A shortfall means the tagging or the one-cell dilation
  is losing a fragment. */
  double Vtot = 0.;
  foreach(serial) Vtot += 2.*pi*y*f[]*sq(Delta);
#if _MPI
  MPI_Allreduce (MPI_IN_PLACE, &Vtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif

  int main_id = 0;                                   // largest by volume
  double vmax = -1.;
  for (int j = 0; j < n; j++) if (cV[j] > vmax) { vmax = cV[j]; main_id = j; }

  /**
  ## Cross-check against `getBase.c`'s main-body definition

  Reported, not assumed: `nliq` is the number of PURE-liquid components at
  `getBase.c`'s threshold. If the two definitions ever disagree about which
  body is the main one, that shows up here rather than silently biasing a
  drop count. */
  scalar dl[];
  foreach() dl[] = (f[] > 1. - FPURE);
  int nliq = tag (dl);

  /**
  ## Jet tip of the main body */
  double ztip = -1000., rtip = -1000., vtip = -1000.;
  foreach(serial) {
    if ((int) dm[] - 1 != main_id) continue;
    if (f[] <= FEPS || f[] >= 1. - FEPS || y > R_TIP) continue;
    if (x > ztip) { ztip = x; rtip = y; vtip = u.x[]; }
  }

  int ndrop = n - 1;
  double Rv = pow (3.*cV[main_id]/(4.*pi), 1./3.);
  double Rs = sqrt (cS[main_id]/(4.*pi));
  fprintf (ferr,
    "MAIN %.8f %d %d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.3f\n",
    t, nliq, ndrop, cV[main_id], Rv, cS[main_id], Rs,
    cZ[main_id]/cV[main_id], cZMIN[main_id], cZMAX[main_id], cRMAX[main_id],
    cUZ[main_id]/cV[main_id], cUR[main_id]/cV[main_id], cEK[main_id],
    ztip, rtip, vtip, Vtot, cDMIN[main_id],
    cDMIN[main_id] > 0. && cDMIN[main_id] < HUGE ? Rv/cDMIN[main_id] : -1.);

  for (int j = 0; j < n; j++) {
    if (j == main_id || cV[j] <= 0.) continue;
    fprintf (ferr,
      "DROP %.8f %d %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.3f\n",
      t, j + 1, cV[j], pow (3.*cV[j]/(4.*pi), 1./3.), cS[j], sqrt (cS[j]/(4.*pi)),
      cZ[j]/cV[j], cZMIN[j], cZMAX[j], cRMAX[j], cUZ[j]/cV[j], cUR[j]/cV[j], cEK[j],
      cDMIN[j],
      cDMIN[j] > 0. && cDMIN[j] < HUGE ? pow (3.*cV[j]/(4.*pi), 1./3.)/cDMIN[j] : -1.);
  }
  fflush (ferr);

  free (cV); free (cS); free (cZ); free (cUZ); free (cUR); free (cEK);
  free (cZMIN); free (cZMAX); free (cRMAX); free (cDMIN);
  return 0;
}
