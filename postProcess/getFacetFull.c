/**
# getFacetFull.c — interface facets from a two-phase axisymmetric dump

Same job as `getFacet.c` (extract the interface facets of `f` from a dump and
write them to stderr as `x y` segments), but compiled against the **full solver
header stack** so it can `restore()` the drill / bursting-bubble solver's dumps.

## Why this exists

The bare `getFacet.c` declares only `scalar f[]` and includes the minimal
`utils.h / output.h / fractions.h`. That is fine for a plain VOF dump, but the
bursting-bubble solver (`burstingBubble-drillResolution.c`) writes dumps that
carry the **axisymmetric metric** (`cm`, `fm`) and the full
`navier-stokes/centered` + `two-phase` + `conserving` + `tension` field set
(`u`, `g`, `rhov`, `sf`, ...). `restore()` into a program that has not
allocated those fields (nor `axi`'s metric) either segfaults or aborts.

Including the same headers as the solver, in the same order, makes the restored
field layout match the dump exactly, so `restore()` succeeds. `getBase.c`
already includes this stack and works as-is; only facet extraction needed a
matching build.

## Usage

    qcc -O2 -disable-dimensions getFacetFull.c -o getFacetFull -lm
    ./getFacetFull <dumpfile> 2> facets.dat      # facets go to stderr

Output: `x y` interface segments (axi: x = axial, y = radial >= 0), one blank
line between segments — the same format `getFacet.c` emits and that
`plotJetMetricsTheory.py --facet` reads for the inception cone fit.

@author Vatsal Sanjay (vatsal.sanjay@comphy-lab.org) / CoMPhy Lab
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "output.h"

int main(int argc, char const *arguments[]) {
  if (argc < 2) {
    fprintf(ferr, "usage: %s <dumpfile>   (facets -> stderr)\n", arguments[0]);
    return 1;
  }
  restore(file = arguments[1]);
  output_facets(f, ferr);
  fflush(ferr);
  return 0;
}
