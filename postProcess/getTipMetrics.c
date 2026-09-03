/**
# getTipMetrics.c — incipient-jet tip curvature and velocity

Extract one row of tip diagnostics from a Basilisk dump. The liquid component
with the largest cell count defines the connected pool/jet and excludes shed
tip droplets. The highest VOF facet endpoint on the symmetry axis defines the
geometric tip. Its owning interfacial cell supplies the axisymmetric mean
curvature and cell-centred velocity; the helper refuses to substitute a nearby
cell when the curvature is unavailable.

For a locally spherical or paraboloidal axisymmetric apex, the two principal
curvatures are equal and the corresponding apex radius is
$R_\kappa=2/|\kappa|$. The helper emits the raw mean curvature as well as the
local grid spacing so downstream analysis can reject grid-limited values.
This operational quantity is not assumed to equal the theoretical minimum jet
radius $R_m$; that identification requires a separate convergence test.

Coordinates follow the simulation convention: `x` is axial ($z$), `y` is
radial ($r$), and `u.x` is the axial velocity.

## Output

One stderr line beginning with `TIP_METRICS`:

```
TIP_METRICS t z_tip r_tip z_cell r_cell kappa u_z u_r speed Delta level f
            n_components cs_h cs_f cs_a cs_c
```

Compile once and run one serial process per snapshot:

```
qcc -O2 -Wall -disable-dimensions getTipMetrics.c -o getTipMetrics -lm
./getTipMetrics snapshot-0.5 2> tip-metrics.txt
```
*/

#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "fractions.h"
#include "curvature.h"
#include "tag.h"

#include <float.h>
#include <limits.h>

#define AXIS_TOLERANCE 1e-10

static int largest_component (scalar labels, int count)
{
  int main = 0;
  if (count > 0) {
    double * sizes = calloc (count, sizeof(double));
    if (!sizes) {
      fprintf (stderr, "getTipMetrics: component-size allocation failed\n");
      return -1;
    }
    foreach (serial)
      if (labels[] > 0.)
        sizes[(int) labels[] - 1] += 1.;
    double largest = -1.;
    for (int index = 0; index < count; index++)
      if (sizes[index] > largest) {
        largest = sizes[index];
        main = index + 1;
      }
    free (sizes);
  }
  return main;
}

int main (int argc, char const * argv[])
{
  if (argc != 2) {
    fprintf (stderr, "usage: %s <dumpfile>\n", argv[0]);
    return 2;
  }
  if (!restore (file = argv[1])) {
    fprintf (stderr, "%s: could not restore '%s'\n", argv[0], argv[1]);
    return 1;
  }
#if TREE
  f.prolongation = fraction_refine;
#endif
  boundary ((scalar *) {f, u.x, u.y});

  /**
  ## Isolate the connected liquid body

  The largest component of `f > 1e-4` contains the pool and the connected jet.
  A detached tip droplet is therefore excluded after pinch-off. */
  scalar liquid_component[];
  foreach()
    liquid_component[] = f[] > 1e-4;
  int component_count = tag (liquid_component);
  int main_liquid = largest_component (liquid_component, component_count);
  if (main_liquid < 0)
    return 6;
  boundary ((scalar *) {liquid_component});
  if (!main_liquid) {
    fprintf (stderr, "%s: no connected liquid component in '%s'\n", argv[0], argv[1]);
    return 3;
  }

  scalar kappa[];
  cstats cs = curvature (f, kappa);

  /**
  ## Locate the main-interface intercept with the symmetry axis

  Restricting a maximum-$z$ search to a broad near-axis band can select the
  cavity shoulder before the axial jet protrudes. Instead, require the PLIC
  segment to touch $r=0$ and select the highest such endpoint. */
  double z_tip = -DBL_MAX, r_tip = DBL_MAX;
  double z_cell = -DBL_MAX, r_cell = -DBL_MAX;
  double kappa_tip = nodata, uz_tip = nodata, ur_tip = nodata;
  double delta_tip = nodata, fraction_tip = nodata;
  int level_tip = -1;
  foreach (serial)
    if (f[] > 1e-6 && f[] < 1. - 1e-6 &&
        (int) liquid_component[] == main_liquid && y <= Delta) {
      coord normal = interface_normal (point, f);
      double alpha = plane_alpha (f[], normal);
      coord segment[2];
      if (facets (normal, alpha, segment) == 2)
        for (int endpoint = 0; endpoint < 2; endpoint++) {
          double z = x + segment[endpoint].x*Delta;
          double r = y + segment[endpoint].y*Delta;
          double tolerance = AXIS_TOLERANCE*max(1., Delta);
          if (fabs(r) <= tolerance &&
              (z > z_tip || (z == z_tip && fabs(r) < fabs(r_tip)))) {
            z_tip = z;
            r_tip = r;
            z_cell = x;
            r_cell = y;
            kappa_tip = kappa[];
            uz_tip = u.x[];
            ur_tip = u.y[];
            delta_tip = Delta;
            level_tip = level;
            fraction_tip = f[];
          }
        }
    }
  if (z_tip == -DBL_MAX) {
    fprintf (stderr, "%s: no near-axis main-interface facet in '%s'\n", argv[0], argv[1]);
    return 4;
  }

  /**
  ## Require curvature in the axis-intercept cell

  `curvature()` may use height-function, fit, averaged or centroid fallbacks.
  The `cstats` values emitted below count those methods over the full
  interface; they do not identify the method used in this one cell. A nearby
  valid cell would be a different observable, so an unavailable apex value is
  a hard diagnostic failure. */
  if (kappa_tip == nodata || delta_tip <= 0.) {
    fprintf (stderr, "%s: curvature unavailable in the tip cell for '%s'\n",
             argv[0], argv[1]);
    return 5;
  }

  fprintf (stderr,
           "TIP_METRICS %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g "
           "%.17g %.17g %d %.17g %d %d %d %d %d\n",
           t, z_tip, r_tip, z_cell, r_cell, kappa_tip, uz_tip, ur_tip,
           sqrt(sq(uz_tip) + sq(ur_tip)), delta_tip, level_tip, fraction_tip,
           component_count, cs.h, cs.f, cs.a, cs.c);
  fflush (stderr);
  return 0;
}
