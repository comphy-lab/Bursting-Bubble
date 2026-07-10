#!/usr/bin/env python3
"""Shared metadata and atomic-output helpers for the paper capsule."""

from __future__ import annotations

import json
import os
import tempfile
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
METADATA_PATH = ROOT / "metadata.json"


def load_metadata(path: Path = METADATA_PATH) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        metadata = json.load(handle)
    if metadata.get("schema_version") != 1:
        raise ValueError(f"Unsupported capsule metadata schema in {path}")
    return metadata


def require_nonempty(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Expected a non-empty output: {path}")


def _render_figure_temp(fig, output: Path, **savefig_kwargs: Any) -> Path:
    output.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        prefix=f".{output.stem}.", suffix=output.suffix, dir=output.parent
    )
    os.close(descriptor)
    temporary = Path(temporary_name)
    try:
        fig.savefig(temporary, format=output.suffix.lstrip("."), **savefig_kwargs)
        require_nonempty(temporary)
    except Exception:
        temporary.unlink(missing_ok=True)
        raise
    return temporary


def atomic_savefig(
    fig,
    output: Path,
    *,
    companion_png: bool = True,
    **savefig_kwargs: Any,
) -> tuple[Path, ...]:
    """Save a requested figure and its PNG companion via sibling temp files."""
    output = output.resolve()
    outputs = [output]
    if companion_png and output.suffix.lower() != ".png":
        outputs.append(output.with_suffix(".png"))

    temporaries: list[Path] = []
    try:
        for target in outputs:
            temporaries.append(_render_figure_temp(fig, target, **savefig_kwargs))
        for temporary, target in zip(temporaries, outputs):
            os.replace(temporary, target)
    finally:
        for temporary in temporaries:
            temporary.unlink(missing_ok=True)
    return tuple(outputs)
