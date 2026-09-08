"""
layer4/test_sigma_cycles.py — tests for the alignment/disorder cycle diagnostic.

Written by the coordinator after three successive subagent failures left
layer4/sigma_cycles.py on disk completely untested (and an untested module
immediately cost two rounds of wrong-key debugging in manual use).

Design note on the synthetic fixtures: `_find_local_extrema` only inspects
interior indices, so every fixture below PADS the series with a higher value
before the first minimum and after the last minimum, otherwise the boundary
minima are never detected. Fixtures also pass smoothing_window=1 (which is a
no-op copy) so that extremum positions are exactly the constructed ones and
the tests are not testing the moving average by accident; smoothing and the
prominence filter get their own dedicated test.
"""

from __future__ import annotations

import math
import os
import sqlite3
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))

from layer4.sigma_cycles import (  # noqa: E402
    analyze_logged_experiment,
    build_cycle_log_row,
    cascade_verdict,
    cycle_ratios,
    detect_cycles,
    log_cycle_analysis,
)


# --- Fixtures: synthetic sigma_star(t) series with known cycle structure ---

def _cascade_series() -> tuple[np.ndarray, np.ndarray]:
    """3 cycles, periods 4 -> 2 -> 1 (ratio 0.5), minima 1.0 -> 0.8 -> 0.6 -> 0.4."""
    times = np.array([-1.0, 0.0, 2.0, 4.0, 5.0, 6.0, 6.5, 7.0, 8.0])
    sigma = np.array([3.0, 1.0, 5.0, 0.8, 5.0, 0.6, 5.0, 0.4, 3.0])
    return times, sigma


def _persistent_series() -> tuple[np.ndarray, np.ndarray]:
    """3 cycles, constant period 2, constant minima 1.0."""
    times = np.array([-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0])
    sigma = np.array([3.0, 1.0, 5.0, 1.0, 5.0, 1.0, 5.0, 1.0, 3.0])
    return times, sigma


def _rippled_series() -> tuple[np.ndarray, np.ndarray]:
    """One genuine large cycle with small sub-prominence ripples near the peak."""
    times = np.arange(11, dtype=float)
    sigma = np.array([3.0, 1.0, 4.0, 3.9, 4.1, 5.0, 4.1, 3.9, 4.0, 1.0, 3.0])
    return times, sigma


# --- detect_cycles ---

def test_detect_cycles_finds_three_shrinking_cycles() -> None:
    times, sigma = _cascade_series()
    cycles = detect_cycles(times, sigma, smoothing_window=1)
    assert len(cycles) == 3
    assert [round(c["period"], 6) for c in cycles] == [4.0, 2.0, 1.0]


def test_detect_cycles_reports_values_from_the_original_series() -> None:
    times, sigma = _cascade_series()
    cycles = detect_cycles(times, sigma, smoothing_window=1)
    assert [c["sigma_min_start"] for c in cycles] == [1.0, 0.8, 0.6]
    assert cycles[-1]["sigma_min_end"] == pytest.approx(0.4)
    assert all(c["sigma_max"] == pytest.approx(5.0) for c in cycles)


def test_detect_cycles_shares_boundary_minimum_between_consecutive_cycles() -> None:
    times, sigma = _cascade_series()
    cycles = detect_cycles(times, sigma, smoothing_window=1)
    for earlier, later in zip(cycles, cycles[1:]):
        assert earlier["t_min_end"] == pytest.approx(later["t_min_start"])


def test_detect_cycles_amplitude_is_peak_minus_mean_of_bounding_minima() -> None:
    times, sigma = _cascade_series()
    first = detect_cycles(times, sigma, smoothing_window=1)[0]
    assert first["amplitude"] == pytest.approx(5.0 - 0.5 * (1.0 + 0.8))


def test_detect_cycles_monotonic_ramp_has_no_cycles() -> None:
    times = np.arange(12, dtype=float)
    sigma = np.linspace(0.5, 9.0, 12)
    assert detect_cycles(times, sigma, smoothing_window=1) == []


def test_prominence_filter_suppresses_subthreshold_ripples() -> None:
    """With the default prominence the ripples collapse to one cycle; with the
    filter disabled the same series fragments into three."""
    times, sigma = _rippled_series()
    filtered = detect_cycles(times, sigma, smoothing_window=1, min_prominence_ratio=0.25)
    unfiltered = detect_cycles(times, sigma, smoothing_window=1, min_prominence_ratio=0.0)
    assert len(filtered) == 1
    assert filtered[0]["period"] == pytest.approx(8.0)
    assert len(unfiltered) == 3


def test_detect_cycles_does_not_mutate_inputs() -> None:
    times, sigma = _cascade_series()
    t_ref, s_ref = times.copy(), sigma.copy()
    detect_cycles(times, sigma, smoothing_window=1)
    np.testing.assert_array_equal(times, t_ref)
    np.testing.assert_array_equal(sigma, s_ref)


@pytest.mark.parametrize(
    "times, sigma",
    [
        (np.array([]), np.array([])),
        (np.array([0.0, 1.0, 2.0]), np.array([1.0, 2.0])),
        (np.array([0.0]), np.array([1.0])),
    ],
    ids=["empty", "length-mismatch", "single-point"],
)
def test_detect_cycles_rejects_degenerate_input(times, sigma) -> None:
    with pytest.raises(ValueError):
        detect_cycles(times, sigma)


# --- cycle_ratios ---

def test_cycle_ratios_recovers_geometric_shrinkage_and_decreasing_minima() -> None:
    times, sigma = _cascade_series()
    ratios = cycle_ratios(detect_cycles(times, sigma, smoothing_window=1))
    assert ratios["period_ratios"] == pytest.approx([0.5, 0.5])
    assert ratios["geom_mean_period_ratio"] == pytest.approx(0.5)
    assert ratios["sigma_min_values"] == pytest.approx([1.0, 0.8, 0.6, 0.4])
    assert ratios["sigma_min_trend_slope"] < 0.0
    assert ratios["sigma_min_decreasing"] is True


def test_cycle_ratios_flags_constant_cycling_as_not_decreasing() -> None:
    times, sigma = _persistent_series()
    ratios = cycle_ratios(detect_cycles(times, sigma, smoothing_window=1))
    assert ratios["geom_mean_period_ratio"] == pytest.approx(1.0)
    assert ratios["sigma_min_trend_slope"] == pytest.approx(0.0, abs=1e-12)
    assert ratios["sigma_min_decreasing"] is False


def test_cycle_ratios_is_nan_guarded_below_two_cycles() -> None:
    single = [{
        "t_min_start": 0.0, "t_max": 1.0, "t_min_end": 2.0, "period": 2.0,
        "sigma_min_start": 1.0, "sigma_max": 5.0, "sigma_min_end": 1.5,
        "amplitude": 3.75,
    }]
    ratios = cycle_ratios(single)
    assert ratios["n_cycles"] == 1
    assert ratios["period_ratios"] == []
    assert math.isnan(ratios["geom_mean_period_ratio"])
    assert math.isnan(ratios["sigma_min_trend_slope"])
    assert ratios["sigma_min_decreasing"] is False


def test_cycle_ratios_empty_is_well_defined() -> None:
    ratios = cycle_ratios([])
    assert ratios["n_cycles"] == 0
    assert ratios["sigma_min_values"] == []


def test_cycle_ratios_rejects_non_list() -> None:
    with pytest.raises(ValueError):
        cycle_ratios({"not": "a list"})


# --- cascade_verdict ---

def test_cascade_verdict_cascade_requires_shrinking_periods_and_falling_minima() -> None:
    times, sigma = _cascade_series()
    cycles = detect_cycles(times, sigma, smoothing_window=1)
    assert cascade_verdict(cycles, cycle_ratios(cycles)) == "CASCADE-LEANING"


def test_cascade_verdict_persistent_for_constant_cycling() -> None:
    times, sigma = _persistent_series()
    cycles = detect_cycles(times, sigma, smoothing_window=1)
    assert cascade_verdict(cycles, cycle_ratios(cycles)) == "PERSISTENT"


def test_cascade_verdict_insufficient_below_three_cycles() -> None:
    """Two cycles must never be stretched into a trend, even a cascade-shaped one."""
    times = np.array([-1.0, 0.0, 2.0, 4.0, 5.0, 6.0, 7.0])
    sigma = np.array([3.0, 1.0, 5.0, 0.8, 5.0, 0.6, 3.0])
    cycles = detect_cycles(times, sigma, smoothing_window=1)
    assert len(cycles) == 2
    assert cascade_verdict(cycles, cycle_ratios(cycles)) == "INSUFFICIENT-CYCLES"


def test_cascade_verdict_ambiguous_on_mixed_signal() -> None:
    """Shrinking periods but RISING minima is not a cascade."""
    cycles = [{"period": p, "sigma_min_start": s, "sigma_min_end": s + 0.1,
               "amplitude": 1.0, "t_min_start": 0.0, "t_max": 0.5, "t_min_end": 1.0,
               "sigma_max": 5.0}
              for p, s in zip([4.0, 2.0, 1.0], [0.4, 0.6, 0.8])]
    assert cascade_verdict(cycles, cycle_ratios(cycles)) == "AMBIGUOUS"


# --- growth-phase split and DB round-trip ---

def _seed_timeseries_db(path: str, exp_id: str, times, M_vals, sigma_vals) -> None:
    conn = sqlite3.connect(path)
    try:
        conn.execute(
            """
            CREATE TABLE IF NOT EXISTS layer4_bridge_timeseries
            (exp_id TEXT, step INTEGER, t REAL, M REAL, sigma_star REAL)
            """
        )
        conn.executemany(
            "INSERT INTO layer4_bridge_timeseries (exp_id, step, t, M, sigma_star)"
            " VALUES (?, ?, ?, ?, ?)",
            [(exp_id, i, float(t), float(m), float(s))
             for i, (t, m, s) in enumerate(zip(times, M_vals, sigma_vals))],
        )
        conn.commit()
    finally:
        conn.close()


def test_growth_phase_is_restricted_to_data_up_to_and_including_the_M_peak(tmp_path) -> None:
    db = str(tmp_path / "r.db")
    times, sigma = _cascade_series()
    M_vals = np.array([1.0, 2.0, 3.0, 9.0, 4.0, 3.0, 2.0, 1.5, 1.0])  # peak at index 3
    _seed_timeseries_db(db, "EXP-TEST-001", times, M_vals, sigma)

    report = analyze_logged_experiment("EXP-TEST-001", db_path=db)

    assert report["m_peak_index"] == 3
    assert report["m_peak_time"] == pytest.approx(4.0)
    assert report["n_points"] == 9
    assert report["growth_phase"]["n_points"] == 4
    assert report["full_trajectory"]["n_points"] == 9
    # the two later cycles live past the M peak, so the growth phase cannot see them
    assert len(report["growth_phase"]["cycles"]) < len(report["full_trajectory"]["cycles"])


def test_analyze_logged_experiment_rejects_unknown_table_and_missing_exp(tmp_path) -> None:
    db = str(tmp_path / "r.db")
    times, sigma = _persistent_series()
    _seed_timeseries_db(db, "EXP-TEST-002", times, np.ones_like(times), sigma)
    with pytest.raises(ValueError):
        analyze_logged_experiment("EXP-TEST-002", db_path=db, table="not_a_table")
    with pytest.raises(ValueError):
        analyze_logged_experiment("EXP-DOES-NOT-EXIST", db_path=db)


def test_cycle_log_row_round_trips_through_the_database(tmp_path) -> None:
    db = str(tmp_path / "r.db")
    times, sigma = _cascade_series()
    cycles = detect_cycles(times, sigma, smoothing_window=1)
    phase_report = {
        "cycles": cycles,
        "ratios": cycle_ratios(cycles),
        "verdict": cascade_verdict(cycles, cycle_ratios(cycles)),
        "n_points": int(times.size),
    }
    row = build_cycle_log_row("EXP-L4-CYCLE-TEST-001", "EXP-SRC-001",
                              "growth_phase", phase_report, notes="unit test")
    log_cycle_analysis(db, row)

    conn = sqlite3.connect(db)
    try:
        got = conn.execute(
            "SELECT exp_id, source_exp_id, phase, n_cycles, verdict FROM layer4_sigma_cycles"
        ).fetchall()
    finally:
        conn.close()

    assert got == [("EXP-L4-CYCLE-TEST-001", "EXP-SRC-001", "growth_phase", 3,
                    "CASCADE-LEANING")]
