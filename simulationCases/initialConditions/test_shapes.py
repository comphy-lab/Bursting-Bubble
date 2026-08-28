#!/usr/bin/env python3
"""Geometry checks against the shipped DataFiles polylines."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

_DIR = Path(__file__).resolve().parent
if str(_DIR) not in sys.path:
    sys.path.insert(0, str(_DIR))

from generate_bond_shape import _parse_bonds  # noqa: E402
from young_laplace import continuation_ladder, solve_equilibrium  # noqa: E402
from zero_bond import sphere_plane  # noqa: E402

DATA = _DIR.parent / "DataFiles"


def test_bo001_matches_datafile():
    ref = np.loadtxt(DATA / "Bo0.0010.dat")
    shape = solve_equilibrium(1e-3)
    south_delta = abs(shape.south_pole_axial() - ref[0, 0])
    assert south_delta < 0.01, south_delta
    assert abs(shape.volume_err) < 1e-6
    theory = shape.capillary_metric / (np.pi * np.sqrt(3.0))
    assert abs(shape.opening_metric - theory) / theory < 0.08


def test_continuation_ladder():
    assert continuation_ladder(1e-3) == [1e-3]
    ladder = continuation_ladder(10.0, seed=1e-3, max_ratio=2.0)
    assert ladder[0] == 1e-3
    assert ladder[-1] == 10.0
    ratios = [ladder[i + 1] / ladder[i] for i in range(len(ladder) - 1)]
    assert all(r <= 2.01 for r in ratios)


def test_cold_start_continues_in_bond():
    shape = solve_equilibrium(1.0)
    assert abs(shape.volume_err) < 1e-5
    walked = shape.notes.get("continuation")
    assert walked, "expected a Bond continuation ladder"
    assert walked[-1] == 1.0
    assert walked[0] <= 1e-3 * 1.01
    # Fillet must start at the neck, not walk back from the apex: after
    # leaving the south pole the written polyline must not return to the axis.
    left = int(np.argmax(shape.radial > 0.30))
    assert left > 0
    assert shape.radial[left:].min() > 0.15


def test_parse_bonds_rejects_nonfinite_and_negative():
    import argparse

    for bad in ("nan", "inf", "-1", "-0.01"):
        try:
            _parse_bonds(bad)
        except argparse.ArgumentTypeError:
            continue
        raise AssertionError(bad)


def test_zero_bond_sphere_plane():
    ref = np.loadtxt(DATA / "Bo0.0000.dat")
    shape = sphere_plane(delta=0.01, rmax=32.0)
    assert abs(shape.south_pole_axial() - ref[0, 0]) < 1e-6
    assert abs(shape.axial[-1]) < 1e-12
    assert abs(shape.radial[-1] - 32.0) < 1e-8
    # unit sphere: south pole near -2
    assert -2.05 < shape.south_pole_axial() < -2.00
    try:
        sphere_plane(delta=0.5)
    except ValueError:
        pass
    else:
        raise AssertionError("delta=0.5 must be rejected")
    for kwargs in (
        {"delta": float("nan")},
        {"rmax": float("inf")},
        {"n": 0},
    ):
        try:
            sphere_plane(**kwargs)
        except ValueError:
            continue
        raise AssertionError(kwargs)


if __name__ == "__main__":
    test_continuation_ladder()
    test_parse_bonds_rejects_nonfinite_and_negative()
    test_bo001_matches_datafile()
    test_cold_start_continues_in_bond()
    test_zero_bond_sphere_plane()
    print("test_shapes: ok")
