# Tests for peakrefine.py — sub-datapoint peak-position refinement.

import numpy as np
import pytest
from fragalyseqt import peakrefine
from fragalyseqt.peakrefine import (refine_peak_positions, recover_clipped_peak,
                                    recover_saturated_peaks)


def _gaussian(x, mu, sigma, amp=1000.0):
    return amp * np.exp(-((x - mu) ** 2) / (2.0 * sigma ** 2))


def _clipped(x, mu, sigma, amp, ceiling):
    return np.minimum(_gaussian(x, mu, sigma, amp), ceiling)


class TestRefinePeakPositions:
    def test_empty(self):
        out = refine_peak_positions(np.zeros(10), [])
        assert len(out) == 0

    def test_edge_peaks_unchanged(self):
        sig = np.array([5.0, 1.0, 1.0, 1.0, 5.0])
        out = refine_peak_positions(sig, [0, 4])
        np.testing.assert_array_equal(out, [0.0, 4.0])

    def test_symmetric_apex_no_shift(self):
        # Symmetric neighbours → vertex exactly on the sample.
        sig = np.array([0.0, 1.0, 2.0, 1.0, 0.0])
        out = refine_peak_positions(sig, [2])
        assert out[0] == pytest.approx(2.0, abs=1e-12)

    def test_log_parabolic_recovers_gaussian_apex(self):
        # A Gaussian is exactly parabolic in log space, so log-parabolic
        # interpolation recovers the true sub-sample centre to float precision.
        x = np.arange(0, 21, dtype=float)
        for mu in (10.3, 9.8, 12.49, 7.5):
            sig = _gaussian(x, mu, sigma=3.0)
            p = int(np.argmax(sig))
            out = refine_peak_positions(sig, [p])
            assert out[0] == pytest.approx(mu, abs=1e-6), f'mu={mu}'

    def test_log_parabolic_beats_raw_parabola_bias(self):
        # Raw-value parabola is biased on a Gaussian; log-parabolic is not.
        x = np.arange(0, 21, dtype=float)
        mu = 10.35
        sig = _gaussian(x, mu, sigma=2.5)
        p = int(np.argmax(sig))
        refined = refine_peak_positions(sig, [p])[0]
        # raw 3-point parabola on the same samples
        yl, yc, yr = sig[p - 1], sig[p], sig[p + 1]
        raw = p + 0.5 * (yl - yr) / (yl - 2 * yc + yr)
        assert abs(refined - mu) < abs(raw - mu)

    def test_flat_top_uses_plateau_centre(self):
        # Clipped/saturated peak: find_peaks reports the first tied sample;
        # refinement must return the plateau centre instead.
        sig = np.array([0, 5, 10, 32767, 32767, 32767, 10, 5, 0], dtype=float)
        p = int(np.argmax(sig))            # = 3 (first tied ceiling sample)
        out = refine_peak_positions(sig, [p])
        assert out[0] == pytest.approx(4.0, abs=1e-9)   # centre of [3,4,5]

    def test_flat_top_independent_of_ceiling_value(self):
        # Works for a 65535-style ceiling too (no hard-coded 32767).
        sig = np.array([0, 100, 65535, 65535, 65535, 65535, 100, 0],
                       dtype=float)
        p = int(np.argmax(sig))            # = 2
        out = refine_peak_positions(sig, [p])
        assert out[0] == pytest.approx(3.5, abs=1e-9)   # centre of [2,3,4,5]

    def test_noise_shoulder_kept_integer(self):
        # A near-tie on one side makes the parabola vertex leave the interval;
        # the integer position must be kept rather than thrown outside it.
        sig = np.array([0.0, 999.0, 1000.0, 100.0, 0.0])
        out = refine_peak_positions(sig, [2])
        assert -0.5 <= (out[0] - 2.0) <= 0.5

    def test_non_positive_samples_fall_back_to_parabola(self):
        # Baseline-subtracted data with non-positive flanks must not crash on
        # log() and should still produce a bounded sub-sample shift.
        sig = np.array([-5.0, 0.0, 10.0, 2.0, -3.0])
        out = refine_peak_positions(sig, [2])
        assert -0.5 <= (out[0] - 2.0) <= 0.5

    def test_toggle_off_returns_integers(self):
        sig = _gaussian(np.arange(21, dtype=float), 10.3, 3.0)
        peakrefine.REFINE_PEAKS = False
        try:
            out = refine_peak_positions(sig, [int(np.argmax(sig))])
            assert out[0] == 10.0
        finally:
            peakrefine.REFINE_PEAKS = True

    def test_multiple_peaks_mixed(self):
        x = np.arange(0, 40, dtype=float)
        sig = _gaussian(x, 10.4, 2.0) + _gaussian(x, 28.7, 2.0)
        peaks = [int(np.argmax(sig[:20])), 20 + int(np.argmax(sig[20:]))]
        out = refine_peak_positions(sig, peaks)
        assert out[0] == pytest.approx(10.4, abs=1e-3)
        assert out[1] == pytest.approx(28.7, abs=1e-3)


class TestRecoverClippedPeak:
    def test_mild_clipping_recovers_height_pos_width(self):
        x = np.arange(0, 101, dtype=float)
        mu, sigma, amp, ceiling = 50.3, 4.0, 2000.0, 1000.0
        sig = _clipped(x, mu, sigma, amp, ceiling)
        p = int(np.argmax(sig))
        fit = recover_clipped_peak(sig, p)
        assert fit is not None
        assert fit.height == pytest.approx(amp, rel=0.05)
        assert fit.position == pytest.approx(mu, abs=0.3)
        assert fit.sigma == pytest.approx(sigma, rel=0.05)
        assert fit.fwhm == pytest.approx(2.35482 * sigma, rel=0.05)

    def test_moderate_clipping_recovers_height(self):
        x = np.arange(0, 121, dtype=float)
        mu, sigma, amp, ceiling = 60.0, 5.0, 5000.0, 1000.0
        sig = _clipped(x, mu, sigma, amp, ceiling)
        fit = recover_clipped_peak(sig, int(np.argmax(sig)))
        assert fit is not None
        assert fit.height == pytest.approx(amp, rel=0.10)

    def test_not_clipped_returns_none(self):
        x = np.arange(0, 101, dtype=float)
        sig = _gaussian(x, 50.3, 4.0, amp=800.0)  # below any ceiling
        assert recover_clipped_peak(sig, int(np.argmax(sig))) is None

    def test_severe_clipping_refuses(self):
        # True apex 50x the ceiling: only far-tail flanks survive → refuse.
        x = np.arange(0, 201, dtype=float)
        sig = _clipped(x, 100.0, 4.0, 50000.0, 1000.0)
        fit = recover_clipped_peak(sig, int(np.argmax(sig)))
        # Either refused outright, or — if attempted — never wildly low.
        assert fit is None or fit.height >= 1000.0

    def test_height_always_exceeds_ceiling_when_returned(self):
        x = np.arange(0, 121, dtype=float)
        for amp in (1500.0, 3000.0, 8000.0):
            sig = _clipped(x, 60.4, 4.5, amp, 1000.0)
            fit = recover_clipped_peak(sig, int(np.argmax(sig)))
            if fit is not None:
                assert fit.height > 1000.0
                assert fit.height == pytest.approx(amp, rel=0.15)


class TestRecoverSaturatedPeaks:
    def test_flags_and_recovers_clipped_keeps_others(self):
        x = np.arange(0, 160, dtype=float)
        # peak 0: clipped; peak 1: normal unclipped
        sig = (_clipped(x, 40.3, 4.0, 2500.0, 1000.0)
               + _gaussian(x, 120.6, 3.0, 700.0))
        int_pos = [int(np.argmax(sig[:80])), 80 + int(np.argmax(sig[80:]))]
        refined = refine_peak_positions(sig, int_pos)
        heights = np.array([sig[int_pos[0]], sig[int_pos[1]]])
        fwhms = np.array([9.0, 7.0])
        pos, H, F, sat = recover_saturated_peaks(
            sig, int_pos, refined, heights, fwhms)
        assert list(sat) == [True, False]
        assert H[0] == pytest.approx(2500.0, rel=0.08)   # recovered
        assert H[1] == heights[1]                        # untouched
        assert F[1] == fwhms[1]                          # untouched
        assert pos[1] == refined[1]

    def test_no_saturation_all_false(self):
        x = np.arange(0, 60, dtype=float)
        sig = _gaussian(x, 30.2, 4.0, 500.0)
        ip = [int(np.argmax(sig))]
        ref = refine_peak_positions(sig, ip)
        pos, H, F, sat = recover_saturated_peaks(
            sig, ip, ref, [sig[ip[0]]], [9.0])
        assert not sat.any()
