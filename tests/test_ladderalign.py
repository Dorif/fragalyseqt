# Tests for ladderalign.py — ILS peak-to-ladder alignment.
#
# All tests use synthetic data with a known ground-truth linear polynomial:
#   bp = 0.05 * dp + 10
# which maps data-point positions to base-pair sizes.

import numpy as np
import pytest
from fragalyseqt.ladderalign import (
    align_ils_peaks,
    _fit_poly,
    _compute_rss,
    _pattern_match_align,
    _ror,
    _violation_count,
    _fix_monotonicity,
)


# ---------------------------------------------------------------------------
# Shared synthetic ladder generator
# ---------------------------------------------------------------------------

def _make_ladder(sizes, slope=0.05, intercept=10.0, noise_dp=0.0, rng=None):
    if rng is None:
        rng = np.random.default_rng(42)
    sizes = np.array(sizes, dtype=float)
    dp = (sizes - intercept) / slope
    if noise_dp > 0:
        dp = dp + rng.normal(0, noise_dp, size=dp.shape)
    heights = np.full(len(sizes), 1000.0)
    return dp, sizes, heights


# ---------------------------------------------------------------------------
# _fit_poly
# ---------------------------------------------------------------------------

class TestFitPoly:
    def test_linear_exact(self):
        dp = [100.0, 900.0]
        bp = [15.0,  55.0]
        f = _fit_poly(dp, bp, degree=1)
        assert abs(f(100.0) - 15.0) < 1e-9
        assert abs(f(900.0) - 55.0) < 1e-9

    def test_degree_clamped(self):
        dp = [100.0, 200.0]
        bp = [15.0,  20.0]
        f = _fit_poly(dp, bp, degree=3)
        assert f(150.0) == pytest.approx(17.5, abs=0.1)

    def test_quadratic(self):
        dp = [100.0, 500.0, 900.0]
        bp = [15.0,  35.0,  55.0]
        f = _fit_poly(dp, bp, degree=2)
        assert abs(f(100.0) - 15.0) < 0.1


# ---------------------------------------------------------------------------
# _compute_rss
# ---------------------------------------------------------------------------

class TestComputeRss:
    def test_perfect_fit_zero_rss(self):
        dp = [100.0, 500.0, 900.0]
        bp = [15.0,  35.0,  55.0]
        f = _fit_poly(dp, bp, 1)
        assert _compute_rss(f, dp, bp) < 1e-10

    def test_nonzero_rss(self):
        f = _fit_poly([0.0, 1000.0], [0.0, 100.0], 1)
        assert _compute_rss(f, [500.0], [60.0]) > 0


# ---------------------------------------------------------------------------
# _pattern_match_align
# ---------------------------------------------------------------------------

class TestPatternMatchAlign:
    def test_perfect_match_no_noise(self):
        sizes = [75.0, 100.0, 139.0, 150.0, 160.0, 200.0]
        dp, bp, _ = _make_ladder(sizes)
        result, score = _pattern_match_align(dp, sizes)
        assert result is not None
        assert score < 1e-6
        np.testing.assert_allclose(result, dp, rtol=1e-9)

    def test_extra_peaks_at_end(self):
        # Short definition (60–200) with full run reaching 500 bp.
        sizes_def = [75.0, 100.0, 139.0, 150.0, 160.0, 200.0]
        dp_true, _, _ = _make_ladder(sizes_def)
        extra = np.array([dp_true[-1] + 400.0,
                          dp_true[-1] + 900.0,
                          dp_true[-1] + 1500.0])
        dp_all = np.sort(np.concatenate([dp_true, extra]))
        result, score = _pattern_match_align(dp_all, sizes_def)
        assert result is not None
        assert score < 0.01
        np.testing.assert_allclose(result, dp_true, rtol=1e-6)

    def test_noise_peaks_before_ladder(self):
        sizes = [75.0, 100.0, 139.0, 150.0, 160.0, 200.0]
        dp, _, _ = _make_ladder(sizes)
        noise = np.array([dp[0] - 1200.0, dp[0] - 700.0, dp[0] - 300.0])
        dp_all = np.sort(np.concatenate([noise, dp]))
        result, score = _pattern_match_align(dp_all, sizes)
        assert result is not None
        assert score < 0.01
        np.testing.assert_allclose(result, dp, rtol=1e-6)

    def test_noise_before_and_after(self):
        # Noise at both ends — correct window is in the middle.
        sizes = [75.0, 100.0, 139.0, 150.0, 160.0, 200.0]
        dp, _, _ = _make_ladder(sizes)
        noise_before = np.array([dp[0] - 800.0, dp[0] - 300.0])
        noise_after  = np.array([dp[-1] + 400.0, dp[-1] + 900.0])
        dp_all = np.sort(np.concatenate([noise_before, dp, noise_after]))
        result, score = _pattern_match_align(dp_all, sizes)
        assert result is not None
        assert score < 0.01
        np.testing.assert_allclose(result, dp, rtol=1e-6)

    def test_too_few_peaks(self):
        sizes = [75.0, 100.0, 139.0, 150.0, 160.0, 200.0]
        dp = np.array([1000.0, 2000.0, 3000.0])  # only 3 peaks for 6 sizes
        result, score = _pattern_match_align(dp, sizes)
        assert result is None

    def test_score_is_zero_for_perfect_linear_ladder(self):
        # All sizes equally spaced → uniform pattern → score = 0.
        sizes = [100.0, 200.0, 300.0, 400.0, 500.0]
        dp, _, _ = _make_ladder(sizes)
        _, score = _pattern_match_align(dp, sizes)
        assert score < 1e-10


# ---------------------------------------------------------------------------
# align_ils_peaks (integration)
# ---------------------------------------------------------------------------

class TestRor:
    def test_linear_all_ones(self):
        # Equal bp and dp gaps everywhere → ror = 1.
        sizes  = np.array([60, 80, 100, 120, 140], dtype=float)
        window = np.array([600, 800, 1000, 1200, 1400], dtype=float)
        r = _ror(window, sizes)
        np.testing.assert_allclose(r, 1.0, atol=1e-9)

    def test_noise_peak_produces_paired_violation(self):
        # Noise peak at dp=1210 between 120bp (dp=1200) and 140bp (dp=1400).
        # Creates ror << 1 at 120→140 gap, then ror >> 1 at 140→160.
        sizes  = np.array([60, 80, 100, 120, 140, 160], dtype=float)
        window = np.array([600, 800, 1000, 1200, 1210, 1600], dtype=float)
        r = _ror(window, sizes)
        assert r[2] < 0.35   # 120→140 (noise): gap too small
        assert r[3] > 3.0    # 140→160: gap too large (compensates)


class TestViolationCount:
    def test_clean_window(self):
        sizes  = np.array([60, 80, 100, 120, 140], dtype=float)
        window = np.array([600, 800, 1000, 1200, 1400], dtype=float)
        n, ok  = _violation_count(window, sizes)
        assert ok and n == 0

    def test_non_monotone_dp_invalid(self):
        sizes  = np.array([60, 80, 100], dtype=float)
        window = np.array([600, 500, 1000], dtype=float)  # dp[1] < dp[0]
        _, ok  = _violation_count(window, sizes)
        assert not ok

    def test_noise_peak_gives_violations(self):
        sizes  = np.array([60, 80, 100, 120, 140, 160], dtype=float)
        window = np.array([600, 800, 1000, 1200, 1210, 1600], dtype=float)
        n, ok  = _violation_count(window, sizes)
        assert ok and n >= 2


class TestFixMonotonicity:
    def test_clean_window_unchanged(self):
        sizes  = np.array([60, 80, 100, 120, 140], dtype=float)
        dp     = np.array([600, 800, 1000, 1200, 1400], dtype=float)
        result, _ = _fix_monotonicity(dp, sizes, dp.copy())
        np.testing.assert_array_equal(result, dp)

    def test_noise_peak_removed_when_extra_available(self):
        # 7 detected peaks: 6 real ILS + 1 noise (dp=1210) between 120 and 140.
        # Pattern matcher picks [600,800,1000,1200,1210,1600]:
        #   120bp→dp=1200 ✓, 140bp→dp=1210 (noise!) → RoR violation.
        # After removing dp=1210, pattern match on 6 remaining peaks
        # finds [600,800,1000,1200,1400,1600] with 0 violations.
        sizes      = np.array([60, 80, 100, 120, 140, 160], dtype=float)
        dp_all     = np.array([600, 800, 1000, 1200, 1210, 1400, 1600], dtype=float)
        bad_window = np.array([600, 800, 1000, 1200, 1210, 1600], dtype=float)
        n_before, _ = _violation_count(bad_window, sizes)
        fixed, _ = _fix_monotonicity(dp_all, sizes, bad_window)
        n_after, ok = _violation_count(fixed, sizes)
        assert ok and n_after < n_before

    def test_no_left_shift_accepted(self):
        # If removing a peak causes the new window to start earlier (left-shift),
        # the result must be rejected and the original window returned.
        sizes      = np.array([60, 80, 100, 120, 140, 160], dtype=float)
        # dp_all has an extra peak before the ladder (dp=500)
        # and a noise peak (dp=1210) inside it.
        dp_all     = np.array([500, 600, 800, 1000, 1200, 1210, 1400, 1600], dtype=float)
        bad_window = np.array([600, 800, 1000, 1200, 1210, 1600], dtype=float)
        fixed, _ = _fix_monotonicity(dp_all, sizes, bad_window)
        # The fixed window must not start before dp=600 (original anchor)
        assert float(fixed[0]) >= 600.0 - 1.0

    def test_two_satellites_both_removed(self):
        # Satellite peaks on BOTH sides of the genuine 120bp ILS peak.
        # dp_all: real ladder + blob at 1205 (right of 1200) + blob at 1395 (left of 1400).
        # Pattern window: [600,800,1000,1200,1205,1395,1600] → two violations.
        # After two iterations both satellites should be removed.
        sizes      = np.array([60, 80, 100, 120, 140, 160, 180], dtype=float)
        dp_all     = np.array([600, 800, 1000, 1200, 1205, 1395, 1400, 1600, 1800],
                              dtype=float)
        bad_window = np.array([600, 800, 1000, 1200, 1205, 1395, 1600], dtype=float)
        n_before, _ = _violation_count(bad_window, sizes)
        fixed, _ = _fix_monotonicity(dp_all, sizes, bad_window)
        n_after, ok = _violation_count(fixed, sizes)
        assert ok and n_after < n_before

    def test_unchanged_when_no_alternatives(self):
        # Removing any window peak leaves fewer than len(sizes) → no fix possible.
        sizes      = np.array([60, 80, 100, 120, 140, 160], dtype=float)
        bad_window = np.array([600, 800, 1000, 1200, 1210, 1600], dtype=float)
        result, _ = _fix_monotonicity(bad_window.copy(), sizes, bad_window.copy())
        np.testing.assert_array_equal(result, bad_window)


class TestAlignIlsPeaks:
    def test_perfect_ladder(self):
        sizes = [75.0, 100.0, 139.0, 150.0, 160.0,
                 200.0, 250.0, 300.0, 340.0, 350.0]
        dp, bp, heights = _make_ladder(sizes)
        result = align_ils_peaks(dp, sizes, heights)
        assert not np.any(np.isnan(result))
        np.testing.assert_allclose(result, dp, rtol=1e-9)

    def test_extra_peaks_at_end(self):
        # Short definition with longer actual run.
        sizes_full  = [75, 100, 139, 150, 160, 200, 250,
                       300, 340, 350, 400, 450, 490, 500]
        sizes_short = sizes_full[:10]
        dp_full, _, heights_full = _make_ladder(sizes_full)
        result = align_ils_peaks(dp_full, sizes_short, heights_full)
        assert not np.any(np.isnan(result))
        assert np.all(np.diff(result) > 0)

    def test_noise_peaks_around_ladder(self):
        sizes = [80.0, 100.0, 120.0, 140.0, 160.0, 180.0]
        dp, bp, heights = _make_ladder(sizes, noise_dp=0.5)
        rng = np.random.default_rng(7)
        noise_dp = dp[:-1] + rng.uniform(3, 7, len(dp) - 1)
        all_dp = np.sort(np.concatenate([dp, noise_dp]))
        order = np.argsort(np.concatenate([dp, noise_dp]))
        all_h = np.concatenate([heights, np.full(len(noise_dp), 200.0)])[order]
        result = align_ils_peaks(all_dp, sizes, all_h)
        matched = int(np.sum(~np.isnan(result)))
        assert matched >= len(sizes) - 1

    def test_insufficient_peaks_returns_nan(self):
        result = align_ils_peaks([1000.0], [80.0, 100.0], [500.0])
        assert np.all(np.isnan(result))

    def test_result_length_matches_sizes(self):
        sizes = [75.0, 100.0, 139.0, 160.0, 200.0]
        dp, bp, heights = _make_ladder(sizes)
        result = align_ils_peaks(dp, sizes, heights)
        assert len(result) == len(sizes)

    def test_result_monotone(self):
        sizes = [75.0, 100.0, 139.0, 150.0, 160.0, 200.0]
        dp, bp, heights = _make_ladder(sizes)
        result = align_ils_peaks(dp, sizes, heights)
        non_nan = result[~np.isnan(result)]
        assert np.all(np.diff(non_nan) > 0)
