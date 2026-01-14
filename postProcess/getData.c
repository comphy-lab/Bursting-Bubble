/**
# Getting Data from Simulation Snapshot

A utility for extracting field data from Basilisk simulation snapshots onto
a structured Cartesian grid for post-processing and visualization.

## Description

This program samples strain-rate and velocity fields from simulation snapshots
and outputs them on a regular grid suitable for visualization tools. It supports
both axisymmetric and 2D Cartesian geometries.

## Usage

```
./getData <filename> <xmin> <ymin> <xmax> <ymax> <ny>
```

Where:
- `filename`: Path to the Basilisk snapshot file
- `xmin`, `ymin`: Lower bounds of the sampling domain
- `xmax`, `ymax`: Upper bounds of the sampling domain
- `ny`: Number of grid points in y-direction (nx computed automatically)

## Geometry Configuration

Set `AXI=1` for axisymmetric (default) or `AXI=0` for 2D Cartesian:
- **Axisymmetric**: x=radial, y=axial (includes azimuthal D22 term)
- **2D Cartesian**: x=x-coordinate, y=y-coordinate (no D22 term)

To change geometry:
- Method 1: Edit `#define AXI 1` below to `#define AXI 0`
- Method 2: Compile with flag: `qcc -DAXI=0 ...`

## Workflow

1. Parse CLI bounds/grid spacing into `extraction_config`
2. Restore the snapshot via `restore(file=...)`
3. Register each derived scalar in `field_list`
4. Compute fields and interpolate onto regular grid
5. Stream `x y field0 field1 ...` rows to stderr

## Adding New Fields

To add a new derived quantity (e.g., `Aij`):
1. Declare scalar: `scalar Aij[];`
2. Register in `register_fields()`: `field_list = list_add(field_list, Aij);`
3. Compute in `compute_fields()`: `compute_Aij_field(Aij);`
4. Write compute function: `static void compute_Aij_field(scalar target) { ... }`

Author: Vatsal Sanjay (vatsal.sanjay@comphy-lab.org)
Affiliation: CoMPhy Lab, Durham University
*/

#include "utils.h"
#include "output.h"

#ifndef AXI
#define AXI 1
#endif

scalar f[];
vector u[];

/**
## Data Structures
*/
typedef struct {
  char filename[4096];
  double xmin, ymin, xmax, ymax;
  double Deltax, Deltay;
  int nx, ny;
} extraction_config;

scalar D2c[], vel[];
scalar * field_list = NULL;

static int parse_arguments(int argc, char const *argv[],
                           extraction_config *cfg);
static int configure_grid(extraction_config *cfg);
static void register_fields(void);
static void compute_fields(void);
static double ** allocate_field_buffer(const extraction_config *cfg,
                                       int field_count);
static void sample_fields(const extraction_config *cfg, double **field_buffer,
                          int field_count);
static void write_fields(const extraction_config *cfg, double **field_buffer,
                         int field_count, FILE *fp);
static void cleanup_output(FILE *fp, double **field_buffer);
static void compute_D2c_field(scalar target);
static void compute_velocity_field(scalar target);

/**
## Main Function

Entry point for simulation snapshot extraction and processing.
Validates command-line arguments and orchestrates snapshot restoration,
field computation, and grid interpolation.
*/
int main(int a, char const *arguments[])
{
  extraction_config cfg;
  if (!parse_arguments(a, arguments, &cfg))
    return 1;

  if (!configure_grid(&cfg))
    return 1;

  register_fields();
  restore (file = cfg.filename);
  compute_fields();

  int registered_fields = list_len(field_list);
  double ** field =
    allocate_field_buffer(&cfg, registered_fields);
  sample_fields(&cfg, field, registered_fields);

  FILE * fp = ferr;
  write_fields(&cfg, field, registered_fields, fp);
  cleanup_output(fp, field);
}

/**
## Argument Parsing

Read CLI arguments and guard against invalid bounds/grid sizes.
*/
static int parse_arguments(int argc, char const *argv[],
                           extraction_config *cfg)
{
  if (argc != 7) {
    fprintf(stderr, "Error: Expected 6 arguments\n");
    fprintf(stderr,
            "Usage: %s <filename> <xmin> <ymin> "
            "<xmax> <ymax> <ny>\n", argv[0]);
    return 0;
  }

  snprintf(cfg->filename, sizeof(cfg->filename), "%s", argv[1]);
  cfg->xmin = atof(argv[2]);
  cfg->ymin = atof(argv[3]);
  cfg->xmax = atof(argv[4]);
  cfg->ymax = atof(argv[5]);
  cfg->ny = atoi(argv[6]);

  if (cfg->ny <= 0) {
    fprintf(stderr, "Error: ny must be positive.\n");
    return 0;
  }

  if (cfg->xmax <= cfg->xmin || cfg->ymax <= cfg->ymin) {
    fprintf(stderr, "Error: Bounds must satisfy xmax>xmin "
                    "and ymax>ymin.\n");
    return 0;
  }

  return 1;
}

/**
## Grid Configuration

Translate bounds and ny into nx, Δx, Δy for regular sampling.
*/
static int configure_grid(extraction_config *cfg)
{
  cfg->Deltay = (cfg->ymax - cfg->ymin)/((double) cfg->ny);
  cfg->nx = (int) ((cfg->xmax - cfg->xmin)/cfg->Deltay);

  if (cfg->nx <= 0) {
    fprintf(stderr, "Error: Computed nx <= 0. "
                    "Check the provided bounds.\n");
    return 0;
  }

  cfg->Deltax = (cfg->xmax - cfg->xmin)/((double) cfg->nx);
  return 1;
}

/**
## Field Registration

Populate Basilisk list with each scalar field.
To add a new field, declare the scalar at the top and add it here.
*/
static void register_fields(void)
{
  field_list = list_add(field_list, D2c);
  field_list = list_add(field_list, vel);
}

/**
## Field Computation

Dispatch compute callbacks for each registered field.
*/
static void compute_fields(void)
{
  compute_D2c_field(D2c);
  compute_velocity_field(vel);
}

static double ** allocate_field_buffer(const extraction_config *cfg,
                                       int registered_fields)
{
  return (double **) matrix_new (cfg->nx, cfg->ny + 1,
                                 registered_fields*sizeof(double));
}

/**
## Field Sampling

Interpolate every registered scalar on the regular grid.
The matrix layout follows Basilisk's `matrix_new`: row-major on i (x),
with contiguous blocks of `registered_fields` entries per (i, j).
*/
static void sample_fields(const extraction_config *cfg, double **field_buffer,
                          int registered_fields)
{
  for (int i = 0; i < cfg->nx; i++) {
    double x = cfg->Deltax*(i + 1./2) + cfg->xmin;
    for (int j = 0; j < cfg->ny; j++) {
      double y = cfg->Deltay*(j + 1./2) + cfg->ymin;
      int k = 0;
      for (scalar s in field_list)
        field_buffer[i][registered_fields*j + k++] =
          interpolate (s, x, y);
    }
  }
}

/**
## Output Writing

Stream rows in the format: `x y field0 field1 ...` to the output stream.
*/
static void write_fields(const extraction_config *cfg, double **field_buffer,
                         int registered_fields, FILE *fp)
{
  for (int i = 0; i < cfg->nx; i++) {
    double x = cfg->Deltax*(i + 1./2) + cfg->xmin;
    for (int j = 0; j < cfg->ny; j++) {
      double y = cfg->Deltay*(j + 1./2) + cfg->ymin;
      fprintf (fp, "%g %g", x, y);
      int k = 0;
      for (scalar s in field_list)
        fprintf (fp, " %g",
                 field_buffer[i][registered_fields*j + k++]);
      fputc ('\n', fp);
    }
  }
}

static void cleanup_output(FILE *fp, double **field_buffer)
{
  fflush (fp);
  fclose (fp);
  matrix_free (field_buffer);
}

/**
## Strain-Rate Field (D²)

Compute log₁₀(D²) where D² is the second invariant of the strain-rate tensor.

### Geometry-Dependent Formulation

**Axisymmetric (AXI=1, x=radial, y=axial):**
$$D_{11} = \partial u_y/\partial y \quad \text{(axial velocity gradient)}$$
$$D_{22} = u_y/y \quad \text{(azimuthal component)}$$
$$D_{33} = \partial u_x/\partial x \quad \text{(radial velocity gradient)}$$
$$D_{13} = (\partial u_y/\partial x + \partial u_x/\partial y)/2$$
$$D^2 = D_{11}^2 + D_{22}^2 + D_{33}^2 + 2D_{13}^2$$

**2D Cartesian (AXI=0):**
Same as above but without the D₂₂ term.

Returns log₁₀(μᵣ·D²) where μᵣ is the viscosity ratio (1 in liquid, 0.02 in gas).
Floor value of -10 for non-positive values.
*/
static void compute_D2c_field(scalar target)
{
  foreach() {
    double D11 = (u.y[0,1] - u.y[0,-1])/(2*Delta);
#if AXI
    double D22 = (y > 1e-10) ? u.y[]/y : 0.0;  // Epsilon guard for axis
#endif
    double D33 = (u.x[1,0] - u.x[-1,0])/(2*Delta);
    double D13 =
      0.5*((u.y[1,0] - u.y[-1,0] + u.x[0,1] - u.x[0,-1])/(2*Delta));
#if AXI
    double D2 = sq(D11) + sq(D22) + sq(D33) + 2.*sq(D13);
#else
    double D2 = sq(D11) + sq(D33) + 2.*sq(D13);
#endif
    double mu_r = f[] + (1. - f[])*2e-2;  // viscosity ratio: 1 in liquid, 0.02 in gas
    target[] = mu_r * D2;
    if (target[] > 0.)
      target[] = log(target[])/log(10);
    else
      target[] = -10;
  }
}

/**
## Velocity Magnitude Field

Compute velocity magnitude: $|\mathbf{u}| = \sqrt{u_x^2 + u_y^2}$

Geometry-independent calculation:
- Axisymmetric (AXI=1): u.x=radial, u.y=axial
- 2D Cartesian (AXI=0): u.x=x-component, u.y=y-component
*/
static void compute_velocity_field(scalar target)
{
  foreach()
    target[] = sqrt(sq(u.x[]) + sq(u.y[]));
}
