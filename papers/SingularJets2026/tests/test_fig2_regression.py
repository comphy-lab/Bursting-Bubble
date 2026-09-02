from __future__ import annotations

import json
import hashlib
from pathlib import Path

import numpy as np
import pytest

import make_fig2_flux_scalings as flux
import make_fig2_v2 as figure_v2
import make_fig2a_streamlines as figure_2a


CAPSULE = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def series_by_run():
    return [
        (
            run,
            flux.processed_series(
                flux.read_log(CAPSULE / "data-Oh-0.03" / run.filename)
            ),
        )
        for run in flux.RUNS
    ]


def test_canonical_window_is_metadata_backed():
    with (CAPSULE / "metadata.json").open(encoding="utf-8") as handle:
        metadata = json.load(handle)["figure_2"]
    window = metadata["canonical_fit_window"]
    assert (window["minimum"], window["maximum"]) == (0.005, 0.023952)
    assert flux.CONE_FIT_WINDOW == (0.005, 0.023952)


def test_figure_legend_matches_archived_focus_caps():
    assert figure_v2.SHORT_LEGEND_LABELS[r"L13"] == r"Level 13, focus 12"
    assert figure_v2.SHORT_LEGEND_LABELS[r"L14"] == r"Level 14, focus 12"


def test_panel_a_metadata_matches_archive_schema():
    assert figure_2a.FIG2A_METADATA["regular_grid"]["columns"] == list(
        figure_2a.Field.__dataclass_fields__
    )


def test_raw_snapshot_manifest_has_verified_public_downloads():
    with figure_2a.DEFAULT_SNAPSHOT_MANIFEST.open(encoding="utf-8") as handle:
        manifest = json.load(handle)
    assert manifest["link_status"] == "public-downloads-verified"
    assert [entry["snapshot"] for entry in manifest["snapshots"]] == list(
        figure_2a.DEFAULT_SNAPSHOTS
    )
    for entry in manifest["snapshots"]:
        assert entry["url"].startswith("https://www.dropbox.com/scl/fi/")
        assert entry["url"].endswith("&dl=1")


def test_reference_prefactors_and_counts(series_by_run):
    expected = flux.FIG2_METADATA["expected_prefactors"]
    counts = flux.FIG2_METADATA["expected_point_counts"]
    q_slope = (3.0 * flux.ALPHA - 1.0) / flux.ALPHA
    q_line_slope = (2.0 * flux.ALPHA - 1.0) / flux.ALPHA
    we_slope = (3.0 * flux.ALPHA - 2.0) / flux.ALPHA

    assert flux.reference_run_normalisation(
        series_by_run, "Q_j", q_slope, flux.CONE_FIT_WINDOW
    ) == pytest.approx(expected["Q_j"], rel=1e-12)
    assert flux.reference_run_normalisation(
        series_by_run, "q_j", q_line_slope, flux.CONE_FIT_WINDOW
    ) == pytest.approx(expected["q_j"], rel=1e-12)
    assert flux.reference_run_normalisation(
        series_by_run, "We_j", we_slope, flux.CONE_FIT_WINDOW
    ) == pytest.approx(expected["We_j_raw"], rel=1e-12)

    reference = next(
        series
        for run, series in series_by_run
        if run.level == 15 and run.focus == 15
    )
    in_window = (
        (reference["r_j"] >= flux.CONE_FIT_WINDOW[0])
        & (reference["r_j"] <= flux.CONE_FIT_WINDOW[1])
    )
    assert np.count_nonzero(in_window) == counts["reference_window_raw"]


def test_interpolated_panel_c_prefactor(series_by_run):
    q_slope = (3.0 * flux.ALPHA - 1.0) / flux.ALPHA
    we_slope = (3.0 * flux.ALPHA - 2.0) / flux.ALPHA
    interpolated, _ = figure_v2.build_interpolated_we_series(
        series_by_run,
        q_slope,
        flux.CONE_FIT_WINDOW,
        anchor_r=0.1,
        blend_start_r=0.005,
        marker_target=52,
    )
    expected = flux.FIG2_METADATA["expected_prefactors"]
    counts = flux.FIG2_METADATA["expected_point_counts"]
    assert flux.reference_run_normalisation(
        interpolated, "We_j", we_slope, flux.CONE_FIT_WINDOW
    ) == pytest.approx(expected["We_j_interpolated_panel_c"], rel=1e-12)

    reference = next(
        series
        for run, series in interpolated
        if run.level == 15 and run.focus == 15
    )
    in_window = (
        (reference["r_j"] >= flux.CONE_FIT_WINDOW[0])
        & (reference["r_j"] <= flux.CONE_FIT_WINDOW[1])
    )
    assert np.count_nonzero(in_window) == counts["reference_window_interpolated_panel_c"]


def test_default_panel_a_inputs_are_complete_and_offline(monkeypatch):
    monkeypatch.setattr(
        figure_2a.urllib.request,
        "urlopen",
        lambda *_args, **_kwargs: pytest.fail("offline load attempted network access"),
    )
    monkeypatch.setattr(
        figure_2a.shutil,
        "which",
        lambda *_args, **_kwargs: pytest.fail("offline load attempted compiler lookup"),
    )
    fields, segments = figure_2a.load_archived_inputs(
        figure_2a.DEFAULT_DATA_DIR, figure_2a.DEFAULT_SNAPSHOTS
    )
    assert len(fields) == len(segments) == 4
    assert all(field.z.shape[1] == 190 for field in fields)
    assert all(len(snapshot_segments) > 0 for snapshot_segments in segments)


def test_panel_a_manifest_matches_committed_extracts():
    with (figure_2a.DEFAULT_DATA_DIR / "manifest.json").open(encoding="utf-8") as handle:
        manifest = json.load(handle)
    assert manifest["grid"]["shape"] == [294, 190]
    for entry in manifest["files"]:
        path = figure_2a.DEFAULT_DATA_DIR / entry["path"]
        assert path.stat().st_size == entry["bytes"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == entry["sha256"]
