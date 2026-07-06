/**
# getFacetMain.c — MAIN-body interface facets (entrained bubbles/droplets excluded)

Like getFacetFull.c, but outputs only the facets of the MAIN interface — the
outer free surface + collapsing cavity + Worthington jet — excluding the
entrapped bubble at the cavity bottom and any shed satellite droplets/bubbles.

Uses the same tag.h isolation as getBase.c: an interfacial cell is on the main
outer surface iff a face-neighbour is pure gas tagged MainGas (the outer
atmosphere + connected cavity, NOT an entrained bubble) AND a face-neighbour is
pure liquid tagged MainLiq (the main body, NOT a detached droplet). Facets are
emitted (to stderr, `x y` segments) only for those cells.

This is what the self-similar profile figure needs so that min(z) / the shift
reference and the drawn interface are not polluted by the entrapped bubble.

    qcc -O2 -disable-dimensions getFacetMain.c -o getFacetMain -lm
    ./getFacetMain <dumpfile> 2> facets_main.dat

@author Vatsal Sanjay (vatsal.sanjay@comphy-lab.org) / CoMPhy Lab
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "navier-stokes/conserving.h"
#include "tension.h"
#include "fractions.h"
#include "curvature.h"
#include "tag.h"

char filename[256];

static int largest_component(scalar d, int n) {
  int main = 0;
  if (n > 0) {
    double *sz = calloc(n, sizeof(double));
    foreach(serial) if (d[] > 0) sz[(int)d[] - 1] += 1.;
#if _MPI
    MPI_Allreduce(MPI_IN_PLACE, sz, n, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
    double sm = -1.;
    for (int j = 0; j < n; j++) if (sz[j] > sm) { sm = sz[j]; main = j + 1; }
    free(sz);
  }
  return main;
}

int main(int a, char const *arguments[]) {
  if (a < 2) { fprintf(ferr, "usage: %s <dumpfile>\n", arguments[0]); return 1; }
  sprintf(filename, "%s", arguments[1]);
  restore(file = filename);
#if TREE
  f.prolongation = fraction_refine;
#endif
  boundary((scalar *){f});

  scalar dl[]; foreach() dl[] = (f[] > 1. - 1e-4);
  int MainLiq = largest_component(dl, tag(dl));
  scalar dg[]; foreach() dg[] = (f[] < 1e-4);
  int MainGas = largest_component(dg, tag(dg));
  boundary((scalar *){dl, dg});

  FILE *fp = ferr;
  foreach(serial)
    if (f[] > 1e-6 && f[] < 1. - 1e-6) {
      bool touchGas = ((int)dg[1,0] == MainGas) || ((int)dg[-1,0] == MainGas) ||
                      ((int)dg[0,1] == MainGas) || ((int)dg[0,-1] == MainGas);
      bool touchLiq = ((int)dl[1,0] == MainLiq) || ((int)dl[-1,0] == MainLiq) ||
                      ((int)dl[0,1] == MainLiq) || ((int)dl[0,-1] == MainLiq);
      if (!(touchGas && touchLiq)) continue;
      coord n = interface_normal(point, f);
      double alpha = plane_alpha(f[], n);
      coord segment[2];
      if (facets(n, alpha, segment) == 2)
        fprintf(fp, "%g %g\n%g %g\n\n",
                x + segment[0].x*Delta, y + segment[0].y*Delta,
                x + segment[1].x*Delta, y + segment[1].y*Delta);
    }
  fflush(fp);
  return 0;
}
