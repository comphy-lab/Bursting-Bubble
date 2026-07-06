#!/usr/bin/env python3
"""
generate_Bo0_IC.py — Bo=0 bursting-bubble initial condition (sphere-plane
contact), for the drill solver's DataFiles/Bo0.0000.dat.

At Bo=0 the bubble is an UNDEFORMED sphere of radius 1 intersecting the flat
free surface — the physically correct shape of a ~10 micron bubble, vs the
gravity-deformed Bo>0 shapes. A point contact cannot be resolved, so the
contact line is regularised by a fillet of radius ~delta (a small capillary
bridge), leaving an initial rupture "hole" of radius ~2*delta at the top.

Geometry from comphy-lab/Circle-Contacts-Line (GetCircles): a radius-1 cavity
arc from the bottom (on axis) up to the neck, a fillet, and the flat free
surface. That code emits (x_axial-into-liquid, y_radial) with the cavity in
+x and the surface at x=0. The drill solver's convention is the OPPOSITE
axial sign (cavity in -x, south pole near -2, free surface at x=0, then the
solver's init negates the distance field), so we flip x here. The polyline
orientation (cavity bottom -> surface far field) then matches the existing
DataFiles/Bo0.0010.dat exactly, so the same init() code fills liquid outside
the cavity and gas inside + above the surface.

Usage:
    python3 generate_Bo0_IC.py [--delta 0.01] [--rmax 32] [--out DataFiles/Bo0.0000.dat]
"""
import argparse
import numpy as np


def get_circles(delta, rmax, n=2000):
    X1c = -(1 + delta)
    phic1 = np.arcsin(2 * delta)
    phi1 = np.linspace(np.pi, phic1, n)
    X1 = X1c + np.cos(phi1)
    Y1 = np.sin(phi1)

    Yfc = (1 + delta) * np.tan(phic1)
    Rf = (1 + delta) / np.cos(phic1) - 1
    phif = np.linspace(np.pi / 2 - phic1, -np.pi / 2, n)
    Xf = -Rf * np.sin(phif)
    Yf = Yfc - Rf * np.cos(phif)

    # shift so the neck (fillet end) sits at x=0, and mirror
    X1 = -(X1 - Xf[-1])
    Xf = -(Xf - Xf[-1])

    X2 = Xf[-1] * np.ones(n)          # flat free surface at x=0
    Y2 = np.linspace(Yf[-1], rmax, n)  # extend radially past the domain

    X = np.concatenate([X1, Xf, X2])
    Y = np.concatenate([Y1, Yf, Y2])
    return X, Y, Rf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--delta", type=float, default=0.01,
                    help="fillet/neck regularisation scale (initial hole radius ~2*delta)")
    ap.add_argument("--rmax", type=float, default=32.0,
                    help="radial extent of the flat free surface (>= domain size)")
    ap.add_argument("--out", default="simulationCases/DataFiles/Bo0.0000.dat")
    a = ap.parse_args()

    X, Y, Rf = get_circles(a.delta, a.rmax)
    # flip axial sign to the drill-solver convention (cavity in -x, south pole
    # near -2, surface at x=0); keeps the cavity-bottom -> surface ordering.
    Xbb = -X
    Ybb = np.maximum(Y, 1e-8)   # avoid exact 0 on the axis (matches Bo*.dat)

    with open(a.out, "w") as fh:
        for x, y in zip(Xbb, Ybb):
            fh.write("% .7e   % .7e\n" % (x, y))
    print("wrote %s  (delta=%g, Rf=%.5f, %d points)" % (a.out, a.delta, Rf, len(Xbb)))
    print("  x_axial [%.4f, %.4f]  (south pole %.4f, surface 0)" % (Xbb.min(), Xbb.max(), Xbb.min()))
    print("  y_radial[%.4f, %.4f]" % (Ybb.min(), Ybb.max()))


if __name__ == "__main__":
    main()
