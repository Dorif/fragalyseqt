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
