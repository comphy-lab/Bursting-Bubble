/**
# getBase.c — robust jet-base / cavity-focus probe (per-snapshot)

Emits the base of the Worthington jet, robust to satellite bubbles /
pinched fragments that the older "globally lowest interfacial point" probe
(getJetFoot) latches onto.

Protocol (the fix Vatsal asked for -- tag.h to isolate the main body,
then find the base):

  1. tag the PURE-LIQUID cells (f > 1-eps) -> MainLiq = largest component.
     This drops DETACHED liquid droplets (the shed tip droplet, liquid
     satellites).
  2. tag the PURE-GAS cells (f < eps) -> MainGas = largest component
     (the outer atmosphere + the connected cavity). Entrained gas bubbles
     in the liquid are SMALL, SEPARATE gas components -- NOT MainGas.
  3. An interfacial cell is on the OUTER free surface iff a face-neighbour is
     pure gas tagged MainGas AND a face-neighbour is pure liquid tagged
     MainLiq. Entrained-bubble surfaces fail the first test (their gas side is
     a satellite gas region); detached-droplet surfaces fail the second.
  4. base = the lowest such outer-surface cell (min axial x), over y < RCAV.
     Pre-inception this is the collapsing cavity floor (on axis); post-
     inception it becomes the off-axis shoulder where the jet meets the
     receding rim -- one continuous, latch-free definition.

Also emits the jet tip (max axial x near the axis) for reference.

STDERR, one line:  t  z_base r_base  z_tip r_tip  n_out
(sentinel -1000 if no outer-surface cell found)

Coords: x = axial (= z), y = radial (= r >= 0). Newtonian dump: f, u.

@author Vatsal Sanjay (vatsal.sanjay@comphy-lab.org) / CoMPhy Lab
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "fractions.h"
#include "tag.h"

char filename[80];

#define RCAV   1.20      // exclude the flat outer free surface (r spans to L0)
#define R_TIP  0.25      // near-axis band for the jet tip (max z at r->0)

int main(int a, char const *arguments[]) {
  sprintf(filename, "%s", arguments[1]);
  restore(file = filename);
#if TREE
  f.prolongation = fraction_refine;
#endif
  boundary((scalar *){f});

  // --- MainLiq = largest pure-liquid component ---
  scalar dl[];
  foreach() dl[] = (f[] > 1. - 1e-4);
  int nl = tag(dl);
  int MainLiq = 0;
  if (nl > 0) {
    double *sz = calloc(nl, sizeof(double));
    foreach(serial) if (dl[] > 0) sz[(int)dl[] - 1] += 1.;
    double sm = -1.;
    for (int j = 0; j < nl; j++) if (sz[j] > sm) { sm = sz[j]; MainLiq = j + 1; }
    free(sz);
  }

  // --- MainGas = largest pure-gas component (outer atmosphere + cavity) ---
  scalar dg[];
  foreach() dg[] = (f[] < 1e-4);
  int ng = tag(dg);
  int MainGas = 0;
  if (ng > 0) {
    double *sz = calloc(ng, sizeof(double));
    foreach(serial) if (dg[] > 0) sz[(int)dg[] - 1] += 1.;
    double sm = -1.;
    for (int j = 0; j < ng; j++) if (sz[j] > sm) { sm = sz[j]; MainGas = j + 1; }
    free(sz);
  }
  boundary((scalar *){dl, dg});

  // --- lowest OUTER-surface interfacial cell ---
  double zbase = HUGE, rbase = -1.;
  double ztip = -1000., rtip = -1000.;
  int nout = 0;
  foreach(serial) {
    if (f[] <= 1e-6 || f[] >= 1. - 1e-6) continue;   // interfacial only
    if (y > RCAV) continue;
    // tip: highest near-axis interfacial point
    if (y < R_TIP && x > ztip) { ztip = x; rtip = y; }
    // outer-surface test: a face-neighbour is MainGas, another is MainLiq
    bool touchGas = ((int)dg[1,0] == MainGas) || ((int)dg[-1,0] == MainGas) ||
                    ((int)dg[0,1] == MainGas) || ((int)dg[0,-1] == MainGas);
    bool touchLiq = ((int)dl[1,0] == MainLiq) || ((int)dl[-1,0] == MainLiq) ||
                    ((int)dl[0,1] == MainLiq) || ((int)dl[0,-1] == MainLiq);
    if (touchGas && touchLiq) {
      nout++;
      if (x < zbase) { zbase = x; rbase = y; }
    }
  }
  if (rbase < 0.) { zbase = -1000.; rbase = -1000.; }

  fprintf(ferr, "%f %7.6e %7.6e %7.6e %7.6e %d\n",
          t, zbase, rbase, ztip, rtip, nout);
  fflush(ferr); fclose(ferr);
  return 0;
}
