/**
# Polymer stretch and stretch rate from a viscoelastic dump

Answers one question: is the fluid that forms the jet tip actually stretched,
and how hard is the polymer being pulled at the cavity focus?

This decides whether finite extensibility (FENE-P) can change anything. FENE-P
caps the conformation at `tr(A) = L^2`, with `L^2 ~ 1e3-1e4` for a dilute
flexible polymer. If the tip fluid carries `tr(A) << L^2` in an Oldroyd-B run,
the cap never engages there and no choice of `L^2` alters the tip. If instead
the focus carries `tr(A) >> L^2`, Oldroyd-B is over-predicting the stress at
jet inception, and that is where FENE-P will act.

In steady uniaxial extension the FENE-P spring function satisfies
`f ~ 2*lambda*dudz`, so `2*lambda*dudz` is the value `f` that a FENE-P run
would have to carry. That number also sets the cost, because the elastic wave
speed goes as `f` and the timestep with it.

## Fields

Only `f`, `u` and the conformation components are declared. `restore()` matches
a dump's fields by name, so the fields this program does not declare are simply
not restored: no solver header stack is needed, and the same binary reads any
dump that carries these names. See `getFacet.c` for the same argument.

## Output (to stderr, one tagged row per line)

    AXIS    z f trA A11 AThTh uz dudz 2*lambda*dudz
    SUMMARY t maxTrA z_maxTrA r_maxTrA maxRate z_maxRate zTip trA_tip

`AXIS` rows are the cells adjacent to the axis, bottom to top. `zTip` is the
highest near-axis cell with `f > 0.5`, and `trA_tip` the conformation there:
that pair is the tip measurement. Maxima are taken over liquid cells only.

## Usage

    ./getStretch <dumpfile> [lambda]

`lambda` is the relaxation time in code units (equal to `De`, since time is
already scaled by the inertio-capillary time). It only scales the last two
columns; omit it and they are reported with `lambda = 1`.
*/

#include "utils.h"
#include "output.h"
#include "fractions.h"

scalar f[];
vector u[];
scalar A11[], A12[], A22[], AThTh[];

char filename[4096];

#define FLIQ 0.5   // liquid test, same threshold the drop statistics use

int main (int a, char const *arguments[]) {
  if (a < 2) {
    fprintf (ferr, "usage: %s <dumpfile> [lambda]\n", arguments[0]);
    return 1;
  }
  if (snprintf (filename, sizeof(filename), "%s", arguments[1]) >=
      (int) sizeof(filename)) {
    fprintf (ferr, "ERROR: snapshot path exceeds %zu characters\n",
             sizeof(filename) - 1);
    return 1;
  }
  double lambda = (a > 2 ? atof (arguments[2]) : 1.);

  restore (file = filename);

  double maxTrA = -1., zMaxTrA = 0., rMaxTrA = 0.;
  double maxRate = 0., zMaxRate = 0.;
  double zTip = -HUGE, trATip = 0.;

  foreach(serial) {
    if (f[] > FLIQ) {
      double trA = A11[] + A22[] + AThTh[];
      if (trA > maxTrA) { maxTrA = trA; zMaxTrA = x; rMaxTrA = y; }
      double dudz = (u.x[1] - u.x[-1])/(2.*Delta);
      if (fabs(dudz) > fabs(maxRate)) { maxRate = dudz; zMaxRate = x; }
    }
  }

  /**
  The near-axis column. `y < Delta` picks the first cell off the axis; on an
  adaptive grid that cell's size varies, which is why `Delta` is used rather
  than a fixed tolerance. */

  foreach(serial) {
    if (y < Delta) {
      double trA = A11[] + A22[] + AThTh[];
      double dudz = (u.x[1] - u.x[-1])/(2.*Delta);
      fprintf (ferr, "AXIS %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
               x, f[], trA, A11[], AThTh[], u.x[], dudz, 2.*lambda*dudz);
      if (f[] > FLIQ && x > zTip) { zTip = x; trATip = trA; }
    }
  }

  fprintf (ferr, "SUMMARY %.8f %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
           t, maxTrA, zMaxTrA, rMaxTrA, 2.*lambda*maxRate, zMaxRate,
           zTip, trATip);
  fflush (ferr);
  return 0;
}
