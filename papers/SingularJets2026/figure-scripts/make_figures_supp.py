#!/usr/bin/env python3
"""Write the exact supplemental figure PDF renders used by the paper."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "output" / "supplement")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    for name in ("figure-s1.pdf", "figure-s2.pdf", "figure-s3.pdf"):
        source = ROOT / "assets" / "supplement" / name
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copyfile(source, output / name)
    print(f"Wrote supplemental figures to {output}")


if __name__ == "__main__":
    main()
