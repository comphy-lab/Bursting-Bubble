# Initial conditions

Bond-number-dependent cavity shapes for Stage 1. The solvers read

```
simulationCases/DataFiles/Bo%5.4f.dat
```

column 1 axial, column 2 radial, cavity in −axial, far-field free
surface at axial = 0. `distance.h` is MPI-incompatible, so this
polyline is converted to a restart dump in serial and then reused.

Lengths are scaled by the equivalent-sphere radius $R_0$ (gas volume
$4\pi/3$). The Bond number is $\mathcal{B}o=\rho g R_0^2/\gamma$.
That is not the same as a Bond number based on the bottom curvature
radius.

## Young–Laplace, $\mathcal{B}o>0$

`generate_bond_shape.py` ports the MATLAB `InitialCondition.m` driver
(Lhuissier & Villermaux 2012): submerged profile, spherical cap, and
outer meniscus, with nested shooting on the bottom curvature $R_b$
(volume) and the contact angle $\varphi_c$ (meniscus). A cold start
walks Bond from $10^{-3}$ to the target (ratio $\le 2$), using each
accepted $(R_b,\varphi_c)$ as the next guess. A comma-separated list
continues from the previous requested value. Pass `--no-continue` to
force a one-shot solve.

```bash
./.venv/bin/python generate_bond_shape.py --bond 0.001
./.venv/bin/python generate_bond_shape.py --bond 0.01,0.04 --out-dir ../DataFiles
./.venv/bin/python generate_bond_shape.py --bond 10
```

The experimental window $\mathcal{B}o=0.01$–$0.04$ needs those files
generated here. The repository ships `Bo0.0010.dat` as the current
default.

## Zero Bond

`generate_zero_bond.py` is the strictly $\mathcal{B}o=0$ geometry: a
unit circle meeting a line (axisymmetric sphere–plane), regularised
by a fillet of scale `delta`. Same Basilisk convention as the
Young–Laplace files.

```bash
./.venv/bin/python generate_zero_bond.py --delta 0.01 --out ../DataFiles/Bo0.0000.dat
```

## Opening-angle check

`plot_opening_angle.py` rebuilds the $2\alpha_c/\pi$ versus
$(R_c/R_0)\sqrt{\mathcal{B}o}$ comparison, including the digitised
Lhuissier & Villermaux points in `reference/Villermaux.csv`.

```bash
./.venv/bin/python plot_opening_angle.py
```

## Environment

```bash
uv venv --python 3.12 .venv
uv pip install --python .venv/bin/python -r requirements.txt
./.venv/bin/python test_shapes.py
```
