/**
# getJetFoot.c — Worthington jet base/focus probe + base flux (per-snapshot)

Emits, for ONE Basilisk dump, the geometric probe candidates AND the
integrated axial flux through each, so the time-ordered consumer can pick
the regime-appropriate base and build grid-robust observables (the
jet-base flow rate is the validation metric; tip velocity does not
grid-converge — see project notes).

Candidates (over the MAIN connected liquid body; detached drops excluded):
  - (z_low, r_low) : globally lowest interfacial point (min axial x).
  - (z_maxk, r_maxk): point of maximum |curvature| below the free surface.

Base flux, evaluated on a thin axial band  z_b - HB < x < z_b + HB  and
0 < r < r_b, for EACH candidate (r_b = that candidate's radius):
  - q_jet = INT_0^{r_b} u_z r dr    [L^3/T]   (band sum of u_z*y*dA / 2HB)
            physical meaning: the flow rate feeding into the jet,
            q_jet ~ r_jet^2 v_jet ~ r_jet^((3*alpha-1)/alpha).
  - q_l   = INT_0^{r_b} u_z dr      [L^2/T]   (band sum of u_z*dA  / 2HB)
            flow rate per unit length.
where u_z = u.x (axial), r = y (radial), dA = sq(Delta). HB ~ 2 finest cells.
(No 2*pi factor and no f-weighting, per the requested definition; the base
region is essentially all liquid so f-weighting is ~identical.)

Coords: x = axial (= z), y = radial (= r >= 0). Newtonian dump: only f, u.

STDERR, one line:
  t  z_low r_low  z_maxk r_maxk  qjet_low ql_low  qjet_maxk ql_maxk  z_jet
where z_jet = max axial position of the interface near the axis (jet tip).
Sentinel -1000 for any candidate (and its fluxes) that does not exist.

Author: Vatsal Sanjay (vatsal.sanjay@comphy-lab.org) / CoMPhy Lab
*/
#include "axi.h"
#include "navier-stokes/centered.h"
#include "two-phase.h"
#include "fractions.h"
#include "tag.h"
#include "curvature.h"

char filename[80];

// geometry-tuned for case 1000: origin(-6,0), L0=10, free surface near z=0
#define RCAV       1.20    // exclude the flat outer free surface (r spans to L0)
#define ZSURF_CURV 0.0     // max|k| search restricted below this axial level
#define R_TIP      0.25    // near-axis band for the jet tip height z_jet (max z at r->0)

int main(int a, char const *arguments[]) {
  sprintf(filename, "%s", arguments[1]);
  restore(file = filename);
#if TREE
  f.prolongation = fraction_refine;
#endif
  boundary((scalar *){f, u.x, u.y});

  // --- main connected liquid region (exclude detached drops) ---
  scalar d[];
  foreach() d[] = (f[] > 1e-4);
  int n = tag(d);
  int MainPhase = 0;
  if (n > 0) {
    double *sz = calloc(n, sizeof(double));
    foreach(serial) if (d[] > 0) sz[(int)d[] - 1] += 1.;
    double sm = -1.;
    for (int j = 0; j < n; j++) if (sz[j] > sm) { sm = sz[j]; MainPhase = j + 1; }
    free(sz);
  }

  scalar kappa[];
  curvature(f, kappa);

  // --- candidates ---
  double zlow = HUGE, rlow = -1.;
  double zk = -1000., rk = -1000., kmax = -1.;
  double zjet = -1000.;            // jet tip = max axial position near the axis
  double dmin = HUGE;
  foreach(serial) {
    if (Delta < dmin) dmin = Delta;
    if (f[] <= 1e-6 || f[] >= 1. - 1e-6) continue;   // interfacial only
    if (d[] != MainPhase) continue;
    if (y > RCAV) continue;
    if (x < zlow) { zlow = x; rlow = y; }
    if (y < R_TIP && x > zjet) zjet = x;             // tip height (max z at r->0)
    if (x < ZSURF_CURV && kappa[] != nodata) {
      double ak = fabs(kappa[]);
      if (ak > kmax) { kmax = ak; zk = x; rk = y; }
    }
  }
  if (rlow < 0.) { zlow = -1000.; rlow = -1000.; }

  // --- base flux through each candidate (thin axial band, 0<r<r_b) ---
  double HB = 2.0 * dmin;          // band half-width ~ 2 finest cells
  double q1 = 0., ql1 = 0., q2 = 0., ql2 = 0.;
  int have1 = (rlow > 0.), have2 = (rk > 0.);
  if (have1 || have2) {
    foreach(reduction(+:q1) reduction(+:ql1) reduction(+:q2) reduction(+:ql2)) {
      double dA = sq(Delta);
      if (have1 && fabs(x - zlow) < HB && y > 0. && y < rlow) {
        q1  += u.x[] * y * dA;
        ql1 += u.x[] * dA;
      }
      if (have2 && fabs(x - zk) < HB && y > 0. && y < rk) {
        q2  += u.x[] * y * dA;
        ql2 += u.x[] * dA;
      }
    }
    double norm = 2.0 * HB;
    q1 /= norm; ql1 /= norm; q2 /= norm; ql2 /= norm;
  }
  if (!have1) { q1 = -1000.; ql1 = -1000.; }
  if (!have2) { q2 = -1000.; ql2 = -1000.; }

  fprintf(ferr, "%f %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e %7.6e\n",
          t, zlow, rlow, zk, rk, q1, ql1, q2, ql2, zjet);
  fflush(ferr); fclose(ferr);
  return 0;
}
