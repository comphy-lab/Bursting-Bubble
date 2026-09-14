/**
# Getting Facets

A utility for extracting interface facets from fluid simulation data.

## Description
This program extracts and outputs the facets representing the interface
between different phases in a multiphase flow simulation. The facets
define the boundary between fluid phases, useful for geometric analysis
and visualization of the interface morphology.

## Physics Background
In multiphase fluid simulations, interfaces between different fluids are
critical features that determine many physical phenomena like surface tension
effects, droplet formation, and coalescence events. This utility identifies
these interfaces by extracting facets from volume fraction data, allowing
for quantitative analysis of interfacial dynamics.

## Usage

```
./getFacets input_file
```

- Author: Vatsal Sanjay
vatsalsanjay@gmail.com
Physics of Fluids Department
University of Twente

*/

#include "utils.h"
#include "output.h"
#include "fractions.h"

scalar f[];  // Volume fraction field

/**
Facets need the volume fraction and nothing else. `restore()` matches the
dump's fields by name, so a program declaring only `f` reads a dump from any
of this repository's solvers -- Newtonian or viscoelastic -- and the extra
fields it does not declare are simply not restored. Verified on a
`burstingBubbleVE-drillResolution` dump: the facets are byte-identical to
those from a program built with the full viscoelastic header stack. Do not
add `two-phase.h`, `two-phaseVE.h` or the log-conformation header here to
"make restore work"; they are not needed, and pinning the header stack to one
solver is what makes a post-processing tool solver-specific.

The path is held in a buffer large enough for a real campaign run path. The
historical `char filename[80]` here silently overflowed on any path longer
than 79 characters and aborted with "buffer overflow detected" before reading
anything, which is how a campaign case 140 characters deep looked like a
corrupt dump. `getDropStats.c` documents the same hazard.
*/
char filename[4096];

/**
### Main Function

Loads a simulation snapshot and extracts the interface facets.

- Input parameters:
  - `arguments[1]`: Filename of the simulation snapshot to process

- Process:
  1. Restores the simulation state from the specified file
  2. Extracts interface facets from the volume fraction field
  3. Outputs facet data to standard error

- Return value:
  - Returns 0 on successful completion

- Note:
  The facet extraction algorithm identifies where the volume fraction
  field crosses a threshold value (typically 0.5) between adjacent cells.
*/
int main(int a, char const *arguments[]) {
  if (a < 2) {
    fprintf (ferr, "usage: %s <dumpfile>   (facets -> stderr)\n", arguments[0]);
    return 1;
  }
  if (snprintf (filename, sizeof(filename), "%s", arguments[1]) >=
      (int) sizeof(filename)) {
    fprintf (ferr, "ERROR: snapshot path exceeds %zu characters\n",
             sizeof(filename) - 1);
    return 1;
  }
  restore (file = filename);

  output_facets (f, ferr);
  fflush (ferr);
  /* ferr is stderr: flush it, never close it. */

  return 0;
}
