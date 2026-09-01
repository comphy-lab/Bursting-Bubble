from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

import fit_alpha_grid_r1 as fit


def _write_log(path: Path, rows: list[list[float]]) -> None:
    path.write_text(
        "header\n" + " ".join(fit.COLUMNS) + "\n"
        + "\n".join(" ".join(str(value) for value in row) for row in rows)
        + "\n",
        encoding="utf-8",
    )


def _row(i: int, t: float, level: int, r: float, q_jet: float) -> list[float]:
    return [i, 1e-5, t, 0.0, level, r, 0.0, r, 0.0, q_jet, 0.0]


def _synthetic_runs(slope: float, *, count: int = 80):
    radii = np.geomspace(0.025, 0.005, count)
    return [
        fit.RunSeries("case-a", 15, Path("a"), radii, 1.7 * radii**slope),
        fit.RunSeries("case-b", 15, Path("b"), radii, 3.1 * radii**slope),
    ]


def test_restart_rows_are_deduplicated_with_last_value_winning(tmp_path):
    path = tmp_path / "case_L15_log.txt"
    _write_log(
        path,
        [
            _row(1, 0.5, 15, 0.01, 1.0),
            _row(2, 0.500000001, 15, 0.02, 2.0),
            _row(3, 0.6, 15, 0.03, 3.0),
        ],
    )
    parsed = fit.read_log(path)
    assert parsed.level == 15
    assert parsed.values["t"].size == 2
    assert parsed.values["r_base"][0] == pytest.approx(0.02)


def test_log_binning_includes_both_window_boundaries():
    radii = np.array([0.01, np.sqrt(0.0002), 0.02])
    points = fit.log_radius_bins(
        radii, 4.0 * radii**1.25, lower=0.01, upper=0.02, bins=2,
        label="case", level=15
    )
    assert sum(point.raw_count for point in points) == 3
    assert {point.bin_index for point in points} == {0, 1}


def test_same_level_cases_remain_isolated_and_recover_slope():
    expected_slope = 1.41
    binned = fit.bin_runs(_synthetic_runs(expected_slope), lower=0.005,
                          upper=0.025, bins=18, min_occupied_bins=10)
    pooled = fit.pooled_run_intercepts(binned)
    assert pooled.slope == pytest.approx(expected_slope, abs=2e-14)
    assert fit.alpha_from_flux_slope(pooled.slope) == pytest.approx(
        1.0 / (3.0 - expected_slope)
    )
    assert set(pooled.intercepts) == {"case-a", "case-b"}
    assert pooled.intercepts["case-a"] != pooled.intercepts["case-b"]


def test_equal_run_aggregation_is_invariant_to_duplicate_sampling():
    def points(label, slope, count):
        x = np.linspace(-5.0, -3.0, count)
        return [fit.BinnedPoint(label, 15, index, np.exp(value),
                                np.exp(2.0+slope*value), 1)
                for index, value in enumerate(x)]

    baseline = {"coarse": points("coarse", 1.0, 5),
                "fine": points("fine", 2.0, 20)}
    duplicated = {"coarse": baseline["coarse"],
                  "fine": [point for point in baseline["fine"]
                           for _ in range(3)]}
    equal_baseline = fit.pooled_run_intercepts(baseline, aggregation="equal-run")
    equal_duplicated = fit.pooled_run_intercepts(duplicated, aggregation="equal-run")
    point_baseline = fit.pooled_run_intercepts(baseline, aggregation="point-weighted")
    point_duplicated = fit.pooled_run_intercepts(duplicated, aggregation="point-weighted")

    assert equal_duplicated.slope == pytest.approx(equal_baseline.slope)
    assert point_duplicated.slope != pytest.approx(point_baseline.slope)


def test_invalid_diagnostic_aggregation_does_not_erase_selected_fit(monkeypatch):
    original = fit.pooled_run_intercepts

    def fail_equal_run(points_by_run, *, fixed_slope=None,
                       aggregation="point-weighted"):
        if aggregation == "equal-run":
            raise ValueError("diagnostic aggregation failed")
        return original(points_by_run, fixed_slope=fixed_slope,
                        aggregation=aggregation)

    monkeypatch.setattr(fit, "pooled_run_intercepts", fail_equal_run)
    summary, _ = fit.fit_window(
        _synthetic_runs(1.4), lower=0.005, upper=0.025, bins=12,
        min_occupied_bins=6, aggregation="point-weighted")

    assert summary["pooled"]["alpha_flux"] > 0
    alternative = summary["aggregation_comparison"]["equal-run"]
    assert alternative == {"valid": False,
                           "reason": "diagnostic aggregation failed"}


def test_temporal_bootstrap_is_seed_deterministic():
    kwargs = dict(lower=0.005, upper=0.025, bins=12, min_occupied_bins=6,
                  replicates=40, block_size=5, seed=17)
    first = fit.temporal_block_inference(_synthetic_runs(1.4), **kwargs)
    second = fit.temporal_block_inference(_synthetic_runs(1.4), **kwargs)
    assert first == second


def test_insufficient_occupied_bins_is_rejected():
    with pytest.raises(ValueError, match="Insufficient occupied bins"):
        fit.bin_runs([fit.RunSeries("case", 15, Path("case"),
                     np.array([0.01, 0.011]), np.array([1.0, 1.1]))],
                     lower=0.01, upper=0.02, bins=10, min_occupied_bins=3)


def test_we_slope_is_definition_linked_to_flux_slope():
    slopes = fit.derived_slopes(1.37)
    assert slopes["s_We"] == pytest.approx(2.0 * slopes["s_Q"] - 3.0)
    assert slopes["s_q"] == pytest.approx(slopes["s_Q"] - 1.0)


def test_null_bootstrap_accepts_exact_two_thirds_and_rejects_clear_alternative():
    kwargs = dict(lower=0.005, upper=0.025, bins=12, min_occupied_bins=6,
                  replicates=99, block_size=5, seed=81)
    null = fit.temporal_block_inference(_synthetic_runs(1.5), **kwargs)
    alternative = fit.temporal_block_inference(_synthetic_runs(1.1), **kwargs)
    assert null["null_test_p_two_sided"] > 0.5
    assert alternative["null_test_p_two_sided"] < 0.05


def test_nonphysical_flux_slope_is_rejected():
    with pytest.raises(ValueError, match="Nonphysical"):
        fit.alpha_from_flux_slope(3.0)


def test_report_includes_case_isolation_and_per_run_intervals(tmp_path):
    specifications = []
    for index, run in enumerate(_synthetic_runs(1.4), start=1):
        path = tmp_path / f"case-{index}_L15_log.txt"
        rows = [
            _row(row, 0.5 + row * 1e-4, 15, radius, flux / (2.0 * np.pi))
            for row, (radius, flux) in enumerate(zip(run.r_j, run.Q_j), start=1)
        ]
        _write_log(path, rows)
        specifications.append([run.label, "15", str(path)])

    arguments = fit.build_parser().parse_args(
        [
            "--physics-group",
            "synthetic-matched-grid",
            *(item for specification in specifications for item in ["--series", *specification]),
            "--window",
            "0.005",
            "0.025",
            "--bins",
            "12",
            "--min-occupied-bins",
            "6",
            "--bootstrap",
            "20",
            "--block-size",
            "5",
            "--output-json",
            str(tmp_path / "report.json"),
            "--output-csv",
            str(tmp_path / "report.csv"),
        ]
    )
    report = fit.analyse(arguments)

    assert report["method"]["case_isolation"] == "each --series is a separate trajectory"
    assert set(report["fit"]["per_run"]) == {"case-a", "case-b"}
    for result in report["fit"]["per_run"].values():
        inference = result["temporal_block_inference"]
        assert len(inference["alpha_flux_ci95"]) == 2
        assert 0.0 <= inference["null_test_p_two_sided"] <= 1.0
    assert (tmp_path / "report.json").is_file()
    assert (tmp_path / "report.csv").is_file()
