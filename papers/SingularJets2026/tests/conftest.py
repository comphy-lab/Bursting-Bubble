from __future__ import annotations

import sys
from pathlib import Path


CAPSULE = Path(__file__).resolve().parents[1]
FIGURE_SCRIPTS = CAPSULE / "figure-scripts"
sys.path.insert(0, str(FIGURE_SCRIPTS))
