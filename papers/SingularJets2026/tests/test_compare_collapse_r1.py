from __future__ import annotations

from pathlib import Path

import numpy as np

import compare_collapse_r1 as collapse


def _write_facets(path: Path, points: np.ndarray) -> None:
    lines = []
    for first, second in zip(points[:-1], points[1:]):
        lines.extend((f"{first[0]} {first[1]}",
                      f"{second[0]} {second[1]}", ""))
    path.write_text("\n".join(lines), encoding="utf-8")


def _synthetic_facets(root: Path, *, alpha: float = 0.63, t0: float = 0.5) -> None:
    root.mkdir()
    template_z = np.linspace(0.0, 4.0, 31)
    template_r = 0.7+0.18*template_z+0.025*template_z**2
    for phase, prefix, times in (
        ("pre", "facetpremain", [0.490, 0.492, 0.494, 0.496, 0.498, 0.499]),
        ("post", "facetmain", [0.501, 0.502, 0.504, 0.506, 0.508, 0.510]),
    ):
        index = [f"t0 {t0}"]
        for time in times:
            tau = abs(time-t0)
            scale = tau**alpha
            radius = 0.7*scale
            index.append(f"{time:.6f} {radius:.12g} 0.0")
            points = np.column_stack((template_z*scale, template_r*scale))
            _write_facets(root/f"{prefix}_{time:.6f}.txt", points)
        (root/("index_pre.txt" if phase == "pre" else "index.txt")).write_text(
            "\n".join(index)+"\n", encoding="utf-8")


def test_true_exponent_has_lower_shape_collapse_error(tmp_path):
    facets = tmp_path/"facets"
    _synthetic_facets(facets)
    arguments = collapse.build_parser().parse_args([
        "--facets", str(facets),
        "--candidate", "true", "0.63",
        "--candidate", "wrong", str(2/3),
        "--t0-fit-window", "0.008", "0.06",
        "--t0-offset", "0",
        "--tau-window", "0.0005", "0.02",
        "--z-window", "-0.1", "4.1",
        "--samples", "80",
        "--output-json", str(tmp_path/"result.json"),
    ])
    report = collapse.analyse(arguments)
    rows = {row["candidate"]: row for row in report["scenarios"]}

    assert rows["true"]["pre_profile_count"] == rows["wrong"]["pre_profile_count"]
    assert rows["true"]["post_profile_count"] == rows["wrong"]["post_profile_count"]
    assert rows["true"]["post_normalised_chamfer_rms"] < 1e-8
    assert (rows["true"]["post_normalised_chamfer_rms"] <
            rows["wrong"]["post_normalised_chamfer_rms"])
    assert (tmp_path/"result.json").is_file()


def test_equal_snapshot_metric_ignores_duplicate_point_density():
    first = np.column_stack((np.linspace(0, 1, 20), np.linspace(0, 2, 20)))
    second = first+np.array([0.05, 0.0])
    duplicated = np.repeat(second, 3, axis=0)
    # Duplicating a whole profile does not change the order of magnitude and
    # cannot make a displaced profile look perfectly collapsed.
    assert collapse.normalised_chamfer(first, duplicated) > 0
    np.testing.assert_allclose(
        collapse.normalised_chamfer(first, duplicated),
        collapse.normalised_chamfer(first, second), rtol=0.2)


def test_profile_window_is_applied_after_similarity_rescaling(monkeypatch, tmp_path):
    raw = np.array([[0.020, 0.010], [0.025, 0.012], [0.030, 0.014]])
    monkeypatch.setattr(collapse, "chained_file", lambda _: (raw,))

    points = collapse.profile(tmp_path/"unused.txt", z_base=0.0, scale=0.01,
                              z_window=(1.5, 3.1), samples=12)

    assert len(points) == 12
    assert points[:, 1].min() >= 1.5
    assert points[:, 1].max() <= 3.1


def test_profile_sample_budget_represents_every_component(monkeypatch, tmp_path):
    paths = tuple(np.array([[0.0, float(index)], [1.0, float(index)]])
                  for index in range(3))
    allocations = []
    monkeypatch.setattr(collapse, "chained_file", lambda _: paths)

    def record_allocation(path, count):
        allocations.append(count)
        return np.repeat(path[:1], count, axis=0)

    monkeypatch.setattr(collapse, "resample_path", record_allocation)
    points = collapse.profile(tmp_path/"unused.txt", z_base=0.0, scale=1.0,
                              z_window=(-0.1, 1.1), samples=10)

    assert allocations == [4, 3, 3]
    assert len(points) == 10
