# This file is part of FragalyseQt.
#
# FragalyseQt is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
#
# FragalyseQt is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License
# for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with FragalyseQt. If not, see <https://www.gnu.org/licenses/>.

# Sub-datapoint peak-position refinement.
#
# scipy.signal.find_peaks only ever reports the INTEGER index of a local
# maximum, so every peak position carries up to +-0.5 dp of quantisation
# error.  Refining to sub-datapoint precision reduces quantisation noise in
# the ILS sizing curve and in sample-peak sizes — a second-order but free
# precision/reproducibility gain that matters most for microvariant edges and
# for run-to-run sizing SD.
#
# This implementation is "smart" in three ways over a naive 3-point parabola:
#   1. Log-parabolic interpolation.  CE peaks are approximately Gaussian, and
#      a Gaussian is exactly parabolic in log space — interpolating log(y)
#      removes the small systematic bias a raw-value parabola has on Gaussian
#      peaks.  Falls back to a plain parabola when any of the three samples is
#      non-positive (e.g. after baseline subtraction).
#   2. Flat-top / saturation handling.  When a peak clips at the ADC ceiling
#      its apex samples are tied, and find_peaks places the "position" on the
#      FIRST tied sample — biasing it left by up to half the plateau width
#      (often several dp).  The apex of a clipped/flat-topped peak is estimated
#      as the CENTRE of the plateau instead.  This is instrument-agnostic (it
#      detects a flat top rather than hard-coding a 32767 vs 65767 ceiling).
#      Note: this recovers POSITION only — the true height/area of a clipped
#      peak is unrecoverable from the apex samples and is left unchanged.
#   3. Robustness clamp.  The true apex of a genuine local maximum lies within
#      +-0.5 sample of the discrete max; a parabola vertex beyond that interval
#      signals noise or an asymmetric shoulder, so the integer position is kept
#      rather than letting noise throw the estimate outside the sample step.

from collections import namedtuple

from numpy import array, log, exp, sqrt, polyfit, polyval, zeros, float64

# Feature toggle (A/B): set False to keep integer find_peaks positions.
REFINE_PEAKS = True

# FWHM = 2*sqrt(2*ln2) * sigma for a Gaussian.
_FWHM_PER_SIGMA = 2.3548200450309493

# A sample within this fraction of the apex value counts as part of a flat top.
_FLAT_REL = 1e-3
# Number of consecutive at-ceiling samples that defines a clipped/flat peak.
_MIN_PLATEAU_SAMPLES = 3


def refine_peak_positions(signal, positions):
    # Return sub-datapoint-refined float positions for each integer peak index
    # in *positions*, given the channel *signal*.  See module docstring.
    sig = array(signal, dtype=float64)
    pos = array(positions)
    refined = pos.astype(float64)
    if not REFINE_PEAKS or len(pos) == 0:
        return refined
    n = len(sig)
    for k in range(len(pos)):
        p = int(pos[k])
        if p <= 0 or p >= n - 1:
            continue  # edge peak — cannot form a 3-point window
        yc = sig[p]
        # ── flat-top / saturated peak → plateau centre ──────────────────────
        thr = yc - max(abs(yc) * _FLAT_REL, 1.0)
        lo = p
        while lo - 1 >= 0 and sig[lo - 1] >= thr:
            lo -= 1
        hi = p
        while hi + 1 < n and sig[hi + 1] >= thr:
            hi += 1
        if hi - lo + 1 >= _MIN_PLATEAU_SAMPLES:
            refined[k] = (lo + hi) / 2.0
            continue
        # ── log-parabolic (Gaussian-exact) else plain parabola ──────────────
        yl, yr = sig[p - 1], sig[p + 1]
        if yl > 0.0 and yc > 0.0 and yr > 0.0:
            a, b, c = log(yl), log(yc), log(yr)
        else:
            a, b, c = yl, yc, yr
        denom = a - 2.0 * b + c
        if denom == 0.0:
            continue  # degenerate (flat) — keep integer position
        delta = 0.5 * (a - c) / denom
        if -0.5 <= delta <= 0.5:
            refined[k] = p + delta
        # else: vertex outside the sample interval → noise/shoulder, keep int
    return refined


# ── Saturated / clipped-peak height recovery ────────────────────────────────
#
# When a peak rails at the ADC ceiling the apex samples are flattened, so the
# measured height (the ceiling value) and the find_peaks half-max width are
# both wrong.  The UNCLIPPED flank samples still lie on the underlying
# Gaussian, and ln(y) of a Gaussian is exactly a parabola in x — so a
# log-quadratic least-squares fit to those flanks recovers the true apex
# position (mu), width (sigma → FWHM/area) and height (A = the vertex value).
#
# This is an EXTRAPOLATION whose error grows with clipping severity, so the
# estimator REFUSES (returns None) when it cannot be trusted: too few flank
# samples, a non-concave fit, an apex outside the plateau, an implausibly large
# extrapolation, or a poor log-space fit.  Saturated peaks rise far above the
# baseline, so baseline offset is negligible and ignored.

ClippedPeakFit = namedtuple(
    "ClippedPeakFit", "position height fwhm sigma n_flank rms")


def recover_clipped_peak(signal, peak_index, flat_rel=_FLAT_REL,
                         min_plateau=_MIN_PLATEAU_SAMPLES,
                         min_flank_per_side=2, max_flank_per_side=8,
                         max_extrapolation=30.0, max_rms=0.15):
    # Recover the true apex of a clipped/saturated peak from its unclipped
    # flanks, or return None when the peak is not clipped or the estimate is
    # not trustworthy.  See section docstring.
    sig = array(signal, dtype=float64)
    n = len(sig)
    p = int(peak_index)
    if p <= 0 or p >= n - 1:
        return None
    ceiling = float(sig[p])
    if ceiling <= 0.0:
        return None
    # Plateau extent (the clipped region).
    thr = ceiling - max(abs(ceiling) * flat_rel, 1.0)
    lo = p
    while lo - 1 >= 0 and sig[lo - 1] >= thr:
        lo -= 1
    hi = p
    while hi + 1 < n and sig[hi + 1] >= thr:
        hi += 1
    if hi - lo + 1 < min_plateau:
        return None  # not a flat-topped/clipped peak

    # Collect strictly-decreasing, positive, below-ceiling flank samples.
    def _flank(start, step):
        out = []
        i = start
        prev = thr
        while 0 <= i < n and len(out) < max_flank_per_side:
            y = float(sig[i])
            if y <= 0.0 or y >= prev:
                break  # baseline reached, or rising into a neighbouring peak
            out.append((i, y))
            prev = y
            i += step
        return out

    left = _flank(lo - 1, -1)
    right = _flank(hi + 1, 1)
    if len(left) < min_flank_per_side or len(right) < min_flank_per_side:
        return None

    xs = array([x for x, _ in left] + [x for x, _ in right], dtype=float64)
    ys = array([y for _, y in left] + [y for _, y in right], dtype=float64)
    # Centre x on the plateau to keep the quadratic well-conditioned.
    x0 = float(p)
    xc = xs - x0
    ln_y = log(ys)
    a, b, c = polyfit(xc, ln_y, 2)
    if a >= 0.0:
        return None  # not concave — not a real peak shape
    mu = x0 - b / (2.0 * a)
    sigma = sqrt(-1.0 / (2.0 * a))
    height = float(exp(c - b * b / (4.0 * a)))
    # Sanity gates.
    if not (lo - 1.0 <= mu <= hi + 1.0):
        return None
    if height <= ceiling or height > ceiling * max_extrapolation:
        return None
    rms = float(sqrt(((polyval([a, b, c], xc) - ln_y) ** 2).mean()))
    if rms > max_rms:
        return None
    return ClippedPeakFit(mu, height, _FWHM_PER_SIGMA * sigma, sigma,
                          len(left) + len(right), rms)


def recover_saturated_peaks(signal, int_positions, refined_positions,
                            heights, fwhms):
    # For each peak, recover the true apex of clipped peaks.  Returns
    # (positions, heights, fwhms, saturated_mask) where non-clipped peaks keep
    # their refined position / measured height / measured width, and clipped
    # peaks recoverable from their flanks get the Gaussian-fit apex (better
    # position), true height and true FWHM.  Clipped peaks that cannot be
    # recovered are still flagged saturated but keep their measured values.
    pos = array(refined_positions, dtype=float64).copy()
    H = array(heights, dtype=float64).copy()
    F = array(fwhms, dtype=float64).copy()
    sat = zeros(len(int_positions), dtype=bool)
    for k in range(len(int_positions)):
        p = int(int_positions[k])
        # Cheap clip test: is there a flat top at this peak?
        if p <= 0 or p >= len(signal) - 1:
            continue
        yc = float(signal[p])
        thr = yc - max(abs(yc) * _FLAT_REL, 1.0)
        run = 1
        i = p - 1
        while i >= 0 and signal[i] >= thr:
            run += 1
            i -= 1
        i = p + 1
        while i < len(signal) and signal[i] >= thr:
            run += 1
            i += 1
        if run < _MIN_PLATEAU_SAMPLES:
            continue
        sat[k] = True  # clipped (flagged even if not recoverable)
        fit = recover_clipped_peak(signal, p)
        if fit is not None:
            pos[k] = fit.position
            H[k] = fit.height
            F[k] = fit.fwhm
    return pos, H, F, sat
