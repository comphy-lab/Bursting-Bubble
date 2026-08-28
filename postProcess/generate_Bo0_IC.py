#!/usr/bin/env python3
"""Compatibility wrapper. The generator lives in simulationCases/initialConditions."""

from __future__ import annotations

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
_IC = _ROOT / "simulationCases" / "initialConditions"
if str(_IC) not in sys.path:
    sys.path.insert(0, str(_IC))

from generate_zero_bond import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
