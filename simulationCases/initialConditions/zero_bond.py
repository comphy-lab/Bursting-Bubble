"""Strictly zero-Bond cavity: a unit sphere meeting a flat free surface.

The 2-D generating curve is a circle and a line; the solver runs it as
axisymmetric. A point contact is not resolvable, so a fillet of scale
``delta`` regularises the neck (initial hole radius ~ 2 δ).

Geometry follows ``comphy-lab/Circle-Contacts-Line`` (GetCircles) and
the Bursting-Bubble drill convention: cavity in −axial, south pole near
−2, free surface at axial = 0. The same Stage-1 ``init()`` that reads
``Bo0.0010.dat`` then fills liquid outside the cavity.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["ZeroBondShape", "sphere_plane"]


@dataclass
class ZeroBondShape:
    delta: float
    fillet_radius: float
    axial: np.ndarray
    radial: np.ndarray

    def south_pole_axial(self) -> float:
        return float(self.axial.min())


def sphere_plane(delta: float = 0.01, rmax: float = 32.0, n: int = 2000) -> ZeroBondShape:
    if not np.isfinite(delta):
        raise ValueError("delta must be finite")
    if not np.isfinite(rmax):
        raise ValueError("rmax must be finite")
    if not isinstance(n, (int, np.integer)) or int(n) < 2:
        raise ValueError("n must be an integer >= 2")
    if delta <= 0.0:
        raise ValueError("delta must be positive")
    if delta >= 0.5:
        raise ValueError("delta < 0.5 is required for arcsin(2 delta)")
    if rmax <= 0.0:
        raise ValueError("rmax must be positive")
    n = int(n)

    X1c = -(1.0 + delta)
    phic1 = np.arcsin(2.0 * delta)
    phi1 = np.linspace(np.pi, phic1, n)
    X1 = X1c + np.cos(phi1)
    Y1 = np.sin(phi1)

    Yfc = (1.0 + delta) * np.tan(phic1)
    Rf = (1.0 + delta) / np.cos(phic1) - 1.0
    phif = np.linspace(np.pi / 2.0 - phic1, -np.pi / 2.0, n)
    Xf = -Rf * np.sin(phif)
    Yf = Yfc - Rf * np.cos(phif)

    # Shift so the neck sits at x = 0, then mirror into +x (cavity).
    X1 = -(X1 - Xf[-1])
    Xf = -(Xf - Xf[-1])

    X2 = np.full(n, Xf[-1])
    Y2 = np.linspace(Yf[-1], rmax, n)

    axial = -np.concatenate([X1, Xf, X2])
    radial = np.maximum(np.concatenate([Y1, Yf, Y2]), 1e-8)
    return ZeroBondShape(
        delta=float(delta),
        fillet_radius=float(Rf),
        axial=axial,
        radial=radial,
    )


def write_basilisk_dat(shape: ZeroBondShape, path) -> None:
    with open(path, "w") as fh:
        for x, y in zip(shape.axial, shape.radial):
            fh.write("% .7e   % .7e\n" % (x, y))
