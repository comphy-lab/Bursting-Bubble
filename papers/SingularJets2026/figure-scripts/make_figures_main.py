#!/usr/bin/env python3
"""Write exact paper PDFs and regenerate Figure 1 and Figure 4(a,b) panels."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "output")
    parser.add_argument("--no-tex", action="store_true")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)

    command = [sys.executable, str(ROOT / "schematic" / "make_figure1.py"),
               "--output-dir", str(output / "regenerated")]
    if args.no_tex:
        command.append("--no-tex")
    subprocess.run(command, check=True)
    shutil.copyfile(ROOT / "assets" / "main" / "figure1.pdf", output / "figure1.pdf")
    shutil.copyfile(output / "regenerated" / "sketch_paper_worthington.pdf",
                    output / "figure1-redrawn.pdf")
    for script, name in (("make_figure2.py", "figure2.pdf"),
                         ("make_figure4a.py", "figure4a.pdf"),
                         ("make_figure4b.py", "figure4b.pdf")):
        command = [sys.executable, str(ROOT / "figure-scripts" / script),
                   "--output", str(output / name)]
        if args.no_tex and script != "make_figure2.py":
            command.append("--no-tex")
        subprocess.run(command, check=True)
    for name in ("figure3.pdf", "figure4.pdf", "figure5.pdf"):
        source = ROOT / "assets" / "main" / name
        if not source.is_file():
            raise FileNotFoundError(source)
        shutil.copyfile(source, output / name)
    print(f"Wrote main figures to {output}")


if __name__ == "__main__":
    main()
