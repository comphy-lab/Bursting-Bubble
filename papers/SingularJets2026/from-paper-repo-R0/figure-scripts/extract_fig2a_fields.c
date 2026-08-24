/**
 * Extract regular-grid fields for the Fig. 2(a) streamline diagnostic.
 *
 * Usage:
 *   ./extract_fig2a_fields <snapshot> <zmin> <rmin> <zmax> <rmax> <nr>
 *
 * Output columns on stdout:
 *   z r f u_z u_r |u|
 *
 * Basilisk coordinates for these axisymmetric runs are x=z and y=r. The
 * velocity components are therefore u.x=u_z and u.y=u_r.
 */

#include "utils.h"
#include "output.h"

scalar f[];
vector u[];

int main (int argc, char const * argv[])
{
  if (argc != 7) {
    fprintf (stderr,
             "usage: %s <snapshot> <zmin> <rmin> <zmax> <rmax> <nr>\n",
             argv[0]);
    return 1;
  }

  char filename[4096];
  snprintf (filename, sizeof(filename), "%s", argv[1]);

  double zmin = atof (argv[2]);
  double rmin = atof (argv[3]);
  double zmax = atof (argv[4]);
  double rmax = atof (argv[5]);
  int nr = atoi (argv[6]);

  if (zmax <= zmin || rmax <= rmin || nr <= 0) {
    fprintf (stderr, "invalid bounds or grid size\n");
    return 1;
  }

  restore (file = filename);

  double dr = (rmax - rmin)/((double) nr);
  int nz = (int) ((zmax - zmin)/dr);
  if (nz <= 0) {
    fprintf (stderr, "computed nz <= 0\n");
    return 1;
  }

  for (int iz = 0; iz < nz; iz++) {
    double z = zmin + (iz + 0.5)*dr;
    for (int ir = 0; ir < nr; ir++) {
      double r = rmin + (ir + 0.5)*dr;
      double c = interpolate (f, z, r);
      double uz = interpolate (u.x, z, r);
      double ur = interpolate (u.y, z, r);
      fprintf (stdout, "%.12g %.12g %.12g %.12g %.12g %.12g\n",
               z, r, c, uz, ur, sqrt(sq(uz) + sq(ur)));
    }
  }

  return 0;
}
