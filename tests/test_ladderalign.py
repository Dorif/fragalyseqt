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
        noise_after = np.array([dp[-1] + 400.0, dp[-1] + 900.0])
        dp_all = np.sort(np.concatenate([noise_before, dp, noise_after]))
        result, score = _pattern_match_align(dp_all, sizes)
        assert result is not None
        assert score < 0.01
        np.testing.assert_allclose(result, dp, rtol=1e-6)

    def test_too_few_peaks(self):
        sizes = [75.0, 100.0, 139.0, 150.0, 160.0, 200.0]
        dp = np.array([1000.0, 2000.0, 3000.0])
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
        sizes = np.array([60, 80, 100, 120, 140], dtype=float)
        window = np.array([600, 800, 1000, 1200, 1400], dtype=float)
        r = _ror(window, sizes)
        np.testing.assert_allclose(r, 1.0, atol=1e-9)

    def test_noise_peak_produces_paired_violation(self):
        # Noise peak at dp=1210 between 120bp (dp=1200) and 140bp (dp=1400).
        # Creates ror << 1 at 120→140 gap, then ror >> 1 at 140→160.
        sizes = np.array([60, 80, 100, 120, 140, 160], dtype=float)
        window = np.array([600, 800, 1000, 1200, 1210, 1600], dtype=float)
        r = _ror(window, sizes)
        assert r[2] < 0.35
        assert r[3] > 3.0


class TestViolationCount:
    def test_clean_window(self):
        sizes = np.array([60, 80, 100, 120, 140], dtype=float)
        window = np.array([600, 800, 1000, 1200, 1400], dtype=float)
        n, ok = _violation_count(window, sizes)
        assert ok and n == 0

    def test_non_monotone_dp_invalid(self):
        sizes = np.array([60, 80, 100], dtype=float)
        window = np.array([600, 500, 1000], dtype=float)
        _, ok = _violation_count(window, sizes)
        assert not ok

    def test_noise_peak_gives_violations(self):
        sizes = np.array([60, 80, 100, 120, 140, 160], dtype=float)
        window = np.array([600, 800, 1000, 1200, 1210, 1600], dtype=float)
        n, ok = _violation_count(window, sizes)
        assert ok and n >= 2


class TestFixMonotonicity:
    def test_clean_window_unchanged(self):
        sizes = np.array([60, 80, 100, 120, 140], dtype=float)
        dp = np.array([600, 800, 1000, 1200, 1400], dtype=float)
        result, _ = _fix_monotonicity(dp, sizes, dp.copy())
        np.testing.assert_array_equal(result, dp)

    def test_noise_peak_removed_when_extra_available(self):
        # 7 detected peaks: 6 real ILS + 1 noise (dp=1210) between 120 and 140.
        # Pattern matcher picks [600,800,1000,1200,1210,1600]:
        #   120bp→dp=1200 ✓, 140bp→dp=1210 (noise!) → RoR violation.
        # After removing dp=1210, pattern match on 6 remaining peaks
        # finds [600,800,1000,1200,1400,1600] with 0 violations.
        sizes = np.array([60, 80, 100, 120, 140, 160], dtype=float)
        dp_all = np.array([600, 800, 1000, 1200, 1210, 1400, 1600], dtype=float)
        bad_window = np.array([600, 800, 1000, 1200, 1210, 1600], dtype=float)
        n_before, _ = _violation_count(bad_window, sizes)
        fixed, _ = _fix_monotonicity(dp_all, sizes, bad_window)
        n_after, ok = _violation_count(fixed, sizes)
        assert ok and n_after < n_before

    def test_no_left_shift_accepted(self):
        # If removing a peak causes the new window to start earlier (left-shift),
        # the result must be rejected and the original window returned.
        sizes = np.array([60, 80, 100, 120, 140, 160], dtype=float)
        # dp_all has an extra peak before the ladder (dp=500)
        # and a noise peak (dp=1210) inside it.
        dp_all = np.array([500, 600, 800, 1000, 1200, 1210, 1400, 1600], dtype=float)
        bad_window = np.array([600, 800, 1000, 1200, 1210, 1600], dtype=float)
        fixed, _ = _fix_monotonicity(dp_all, sizes, bad_window)
        # The fixed window must not start before dp=600 (original anchor)
        assert float(fixed[0]) >= 600.0 - 1.0

    def test_two_satellites_both_removed(self):
        # Satellite peaks on BOTH sides of the genuine 120bp ILS peak.
        # dp_all: real ladder + blob at 1205 (right of 1200) + blob at 1395 (left of 1400).
        # Pattern window: [600,800,1000,1200,1205,1395,1600] → two violations.
        # After two iterations both satellites should be removed.
        sizes = np.array([60, 80, 100, 120, 140, 160, 180], dtype=float)
        dp_all = np.array([600, 800, 1000, 1200, 1205, 1395, 1400, 1600, 1800],
                          dtype=float)
        bad_window = np.array([600, 800, 1000, 1200, 1205, 1395, 1600], dtype=float)
        n_before, _ = _violation_count(bad_window, sizes)
        fixed, _ = _fix_monotonicity(dp_all, sizes, bad_window)
        n_after, ok = _violation_count(fixed, sizes)
        assert ok and n_after < n_before

    def test_unchanged_when_no_alternatives(self):
        # Removing any window peak leaves fewer than len(sizes) → no fix possible.
        sizes = np.array([60, 80, 100, 120, 140, 160], dtype=float)
        bad_window = np.array([600, 800, 1000, 1200, 1210, 1600], dtype=float)
        result, _ = _fix_monotonicity(bad_window.copy(), sizes, bad_window.copy())
        np.testing.assert_array_equal(result, bad_window)

    def test_fewest_violations_wins_over_pattern_score(self):
        # Regression guard for the candidate-selection criterion.
        #
        # Real GS500 ladder with three satellites hugging the 2580, 5800 and
        # 8800 dp peaks. Several removals are viable here and they differ in
        # how many RoR violations they leave behind, while their pattern
        # scores are nearly identical (all ~0.001x). Selecting purely on the
        # pattern score therefore settles on a window that still carries
        # violations; preferring the candidate with the FEWEST violations
        # (pattern score only breaking ties) clears them completely.
        sizes = np.array([75, 100, 139, 150, 160, 200, 250, 300,
                          340, 350, 400, 450, 490, 500], dtype=float)
        true_dp = (sizes - 10.0) / 0.05
        satellites = np.array([2589.091, 5796.455, 8804.855])
        dp_all = np.sort(np.concatenate([true_dp, satellites]))
        window, _ = _pattern_match_align(dp_all, sizes)
        n_before, ok_before = _violation_count(window, sizes)
        assert ok_before and n_before > 0
        fixed, _ = _fix_monotonicity(dp_all, sizes, window.copy())
        n_after, ok_after = _violation_count(fixed, sizes)
        assert ok_after
        # All violations must be gone, not merely reduced.
        assert n_after == 0, f'{n_before} violations -> {n_after}, expected 0'
        # And the surviving window must be built from genuine ladder peaks.
        kept_satellites = [float(x) for x in fixed
                           if np.min(np.abs(satellites - float(x))) < 1e-6]
        assert len(kept_satellites) <= 1


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
        sizes_full = [75, 100, 139, 150, 160, 200, 250,
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

    def test_pre_ladder_dye_blobs_filtered(self):
        # Mirror the production capillary-electrophoresis case: dye blobs
        # (bleedthrough, unincorporated primers) form a tight cluster of
        # TALL peaks before the actual ILS ladder begins.  These blobs
        # exceed the median*3 strong-peak threshold AND survive the
        # iterative refinement because their heights span the same range
        # as the ladder peaks (sometimes higher).  Only the gap-ratio
        # transition reliably separates them from the ladder.
        sizes_full = [20, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214,
                      220, 240, 250, 260, 280, 300, 314, 320, 340, 360, 380,
                      400, 414, 420, 440, 460, 480, 500, 514, 520, 540, 560,
                      580, 600]
        slope, intercept = 8.5, 0.0
        # Ladder peaks for sizes 60..600 (34 peaks at strong height).
        ladder_dp = np.array([slope * sz + intercept for sz in sizes_full[2:]],
                             dtype=float)
        ladder_h = np.full(len(ladder_dp), 3000.0)
        # Pre-ladder dye blob cluster: 8 tall peaks tightly spaced at low dp,
        # well below the first ladder peak.  Spacing within the cluster
        # must be much smaller than the first ladder gap (~170 dp) so the
        # gap-ratio filter triggers.
        first_ladder_dp = float(ladder_dp[0])
        blob_dp = np.array([200, 215, 230, 280, 295, 310, 325, 345],
                           dtype=float)
        # ensure they sit before the ladder
        assert blob_dp[-1] < first_ladder_dp - 100
        blob_h = np.array([5000, 6000, 12000, 3000, 12000, 7000, 1500, 9000],
                          dtype=float)
        # Add abundant weak noise so the all-peak median sits well below
        # the ladder height — matching real instrument data where most
        # detected peaks are sub-threshold baseline noise.
        rng = np.random.default_rng(2024)
        weak_dp = np.sort(rng.uniform(blob_dp[0],
                                       ladder_dp[-1] + 50, 80))
        weak_h = np.full(80, 150.0)
        all_dp = np.concatenate([blob_dp, ladder_dp, weak_dp])
        all_h = np.concatenate([blob_h, ladder_h, weak_h])
        order = np.argsort(all_dp)
        all_dp = all_dp[order]
        all_h = all_h[order]
        result = align_ils_peaks(all_dp, sizes_full, all_h)
        # Sizes 20 and 40 must be NaN; sizes 60..600 must match ladder
        # peaks even though the blob cluster has TALLER peaks at lower dp.
        assert np.isnan(result[0]) and np.isnan(result[1])
        for i, sz in enumerate(sizes_full[2:], start=2):
            expected = slope * sz + intercept
            assert not np.isnan(result[i]), f'size {sz} should match'
            assert abs(result[i] - expected) < 0.5, (
                f'size {sz}: expected ~{expected}, got {result[i]}')

    def test_boundary_trim_with_middle_noise_spike(self):
        # Synthetic mirror of the SeqStudio GS600 case where the lowest two
        # ladder fragments (20, 40) were not injected, and a single weak
        # noise peak in the middle of the ladder happened to exceed
        # median*3 but sits at <30% of the ladder-peak median.  The
        # boundary-trim subset finder must skip the missing sizes (using
        # strong peaks against a trimmed size_std subset), and the iterative
        # strong-peak refinement must drop the noise spike before alignment.
        sizes_full = [20, 40, 60, 80, 100, 114, 120, 140, 160, 180, 200, 214,
                      220, 240, 250, 260, 280, 300, 314, 320, 340, 360, 380,
                      400, 414, 420, 440, 460, 480, 500, 514, 520, 540, 560,
                      580, 600]
        slope, intercept = 8.5, 0.0
        # Ladder peaks for 60..600 (34 fragments) at strong height.
        ladder_dp = np.array([slope * sz + intercept for sz in sizes_full[2:]],
                             dtype=float)
        ladder_h = np.full(len(ladder_dp), 3000.0)
        # Inject a noise spike between sizes 300 and 314 at moderate height
        # (passes median*3 but well below ladder-peak median).
        dp300 = slope * 300 + intercept
        dp314 = slope * 314 + intercept
        noise_dp = np.array([(dp300 + dp314) / 2])
        noise_h = np.array([700.0])
        # Surround with ~60 weak (h~150) noise peaks scattered throughout to
        # stress the strong-peak iterative refinement.
        rng = np.random.default_rng(2024)
        weak_dp = np.sort(rng.uniform(ladder_dp[0] - 50,
                                       ladder_dp[-1] + 50, 60))
        weak_h = np.full(60, 150.0)
        all_dp = np.concatenate([ladder_dp, noise_dp, weak_dp])
        all_h = np.concatenate([ladder_h, noise_h, weak_h])
        order = np.argsort(all_dp)
        all_dp = all_dp[order]
        all_h = all_h[order]
        result = align_ils_peaks(all_dp, sizes_full, all_h)
        # Sizes 20 and 40 must be NaN (no ladder peak); sizes 60..600 must
        # match the ladder dp positions to within sub-bp precision.
        assert np.isnan(result[0]) and np.isnan(result[1])
        for i, sz in enumerate(sizes_full[2:], start=2):
            expected = slope * sz + intercept
            assert not np.isnan(result[i]), f'size {sz} should match'
            assert abs(result[i] - expected) < 0.5, (
                f'size {sz}: expected ~{expected}, got {result[i]}')


class TestRealLadderGF2800M:
    # Regression for the GS600 LIZ(60-460) sub-ladder mis-alignment on the
    # NIST GF_2800M.hid orange (LIZ) channel.  The lane is hard because:
    #   * the ladder peaks dominate the channel and are uniformly tall, so the
    #     median*3 strong-peak gate selects nothing;
    #   * tall pre-ladder dye blobs (~18000, taller than the ~10-14k ladder)
    #     precede the ladder;
    #   * the selected 60-460 sub-ladder is an interior slice of the full
    #     20-600 standard, whose extra 20/40 and 480-600 fragments are all
    #     physically present, creating registration ambiguity for a normalised-
    #     spacing matcher;
    #   * the 320 fragment is essentially absent (height 99) and a spurious
    #     peak sits at dp 6167.
    # The normalised-spacing objective preferred a mis-registered window; the
    # RoR-minimising consensus aligner recovers the correct, smooth assignment.
    # Ground truth (size -> dp) was confirmed against the smooth ~11 dp/bp
    # migration curve (every size lands on a tall, evenly-spaced real peak).
    GF_PEAKS = [
        (2991, 294), (3011, 18056), (3027, 17307), (3041, 103), (3104, 148),
        (3121, 9253), (3139, 163), (3157, 17118), (3168, 80), (3180, 3486),
        (3199, 363), (3216, 17389), (3230, 219), (3257, 299), (3269, 12448),
        (3288, 18108), (3348, 18347), (3355, 18312), (3399, 18203),
        (3481, 18133), (3496, 18505), (3539, 114), (3594, 7318), (3675, 555),
        (3699, 18446), (3750, 120), (3865, 12045), (4133, 9048), (4393, 13021),
        (4469, 108), (4544, 255), (4570, 14273), (4648, 12112), (4899, 10023),
        (5136, 10455), (5373, 13309), (5602, 11591), (5765, 11478),
        (5832, 9973), (6057, 14301), (6167, 6175), (6280, 12390), (6394, 155),
        (6496, 13281), (6712, 12893), (6862, 9022), (6921, 10174), (6983, 99),
        (7127, 13404), (7228, 132), (7330, 11106), (7532, 8280), (7731, 9744),
        (7871, 10354), (7929, 9070), (8123, 9824), (8312, 8461), (8502, 8090),
        (8684, 3868), (8810, 4144), (8865, 8403), (9048, 8517), (9224, 5933),
        (9396, 7734), (9565, 5428),
    ]
    SIZES = [60, 80, 100, 114, 120, 140, 160, 180, 200, 214, 220, 240, 250,
             260, 280, 300, 314, 320, 340, 360, 380, 400, 414, 420, 440, 460]
    TRUTH = {60: 3865, 80: 4133, 100: 4393, 114: 4570, 120: 4648, 140: 4899,
             160: 5136, 180: 5373, 200: 5602, 214: 5765, 220: 5832, 240: 6057,
             250: 6167, 260: 6280, 280: 6496, 300: 6712, 314: 6862, 320: 6921,
             340: 7127, 360: 7330, 380: 7532, 400: 7731, 414: 7871, 420: 7929,
             440: 8123, 460: 8312}

    def test_gs600_60_460_subladder_aligns_smoothly(self):
        dp = np.array([p for p, _ in self.GF_PEAKS], dtype=float)
        ht = np.array([h for _, h in self.GF_PEAKS], dtype=float)
        result = align_ils_peaks(dp, self.SIZES, ht)
        for sz, got in zip(self.SIZES, result):
            assert not np.isnan(got), f'size {sz} unmatched (NaN)'
            assert abs(got - self.TRUTH[sz]) < 1.0, (
                f'size {sz}: expected dp {self.TRUTH[sz]}, got {got}')
        assert np.all(np.diff(result) > 0)

    def test_extra_band_peaks_flags_subladder(self):
        # Root-cause "D" advisory: the GF lane runs the full 20-600 standard
        # but only 60-460 was selected, so many extra ladder-band peaks remain.
        from fragalyseqt.ladderalign import extra_ladder_band_peaks
        dp = np.array([p for p, _ in self.GF_PEAKS], dtype=float)
        ht = np.array([h for _, h in self.GF_PEAKS], dtype=float)
        assert extra_ladder_band_peaks(dp, ht, self.SIZES) >= 4

    def test_extra_band_peaks_zero_for_exact_ladder(self):
        # A clean ladder whose definition matches the run must not flag.
        from fragalyseqt.ladderalign import extra_ladder_band_peaks
        sizes = [80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 220.0]
        dp = np.array([8.5 * s for s in sizes], dtype=float)
        ht = np.full(len(sizes), 3000.0)
        assert extra_ladder_band_peaks(dp, ht, sizes) == 0
