/**
# getView2D - render one snapshot to a PNG with the adaptive mesh

Restores one `intermediate/snapshot-*` dump from the bursting-bubble
solver (axisymmetric: fields f, u.x, u.y) and renders the interface plus
the adaptive mesh (cells), mirrored across the symmetry axis to reconstruct
the full "bowtie" bubble/jet cross-section. Optionally zoomed to a region
of interest (e.g. the collapsing cavity / jet base) via fov/tx/ty.

2D analogue of comphy-lab/Jumping-Drops postProcess/getView3D_v2.c, adapted
for an axisymmetric quadtree grid (no octree, no 3D camera angles/mirrors —
one mirror across the r=0 axis is enough to reconstruct the full picture).

## Usage

```bash
./getView2D <snapshot> <output.png> [fov tx ty width height]
```

`fov` (field of view, degrees) controls zoom — smaller is more zoomed in.
`tx`, `ty` pan the camera to center on a region of interest; `ty` in
particular should be set to (minus) the axial coordinate of the jet-base /
cavity-focus probe to center the zoom there. Camera defaults frame the whole
box; override for a zoomed region-of-interest render.

@author Vatsal Sanjay (vatsal.sanjay@comphy-lab.org) / CoMPhy Lab
*/

#include <string.h>
#include <stdlib.h>

#include "utils.h"
#include "fractions.h"
#include "view.h"

scalar f[];
char filename[512], Imagename[512];

static void draw_scene()
{
  draw_vof ("f", lw = 3, lc = {0.0, 0.50, 0.0});   // green, matches the marker-view panel
  cells (lw = 1);
}

/**
Mirror across the axisymmetry axis (y = r = 0) to reconstruct the full
bubble/jet cross-section from the half-domain that Basilisk actually solves.
*/
static void draw_mirrored_scene()
{
  draw_scene();
  mirror (n = {0, 1})
    draw_scene();
}

static void draw_time_label()
{
  double time = 0.;
  char label[80];
  char * dash = strrchr (filename, '-');

  if (dash && sscanf (dash + 1, "%lf", &time) == 1)
    snprintf (label, sizeof(label), "t = %.4f", time);
  else
    snprintf (label, sizeof(label), "t = %.4f", t);

  draw_string (label, pos = 2, size = 40, lc = {0., 0., 0.}, lw = 1.2);
}

int main (int a, char const * arguments[])
{
  if (a < 3) {
    fprintf (stderr, "usage: %s snapshot output.png [fov tx ty width height]\n",
             arguments[0]);
    return 1;
  }

  /* A truncated output path used to make save() fail silently (rc = 0, no
     file written) — fail loudly instead of rendering into the void. */
  if (snprintf (filename, sizeof(filename), "%s", arguments[1]) >= (int) sizeof(filename) ||
      snprintf (Imagename, sizeof(Imagename), "%s", arguments[2]) >= (int) sizeof(Imagename)) {
    fprintf (stderr, "%s: path longer than %zu characters\n",
             arguments[0], sizeof(filename) - 1);
    return 1;
  }
  restore (file = filename);

  double fov = 24.;   // degrees; smaller = more zoomed in
  double tx = 0., ty = 0.;
  int width = 1000, height = 1000, label = 1;
  if (a > 3) fov    = atof (arguments[3]);
  if (a > 4) tx      = atof (arguments[4]);
  if (a > 5) ty      = atof (arguments[5]);
  if (a > 6) width   = atoi (arguments[6]);
  if (a > 7) height  = atoi (arguments[7]);
  if (a > 8) label   = atoi (arguments[8]);   // 0 to suppress the t=... label

  view (fov = fov, tx = tx, ty = ty,
        width = width, height = height, samples = 4,
        bg = {1.0, 1.0, 1.0});

  draw_mirrored_scene();
  if (label) draw_time_label();
  save (Imagename);
}
