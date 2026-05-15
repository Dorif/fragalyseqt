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

# ILS (Internal Lane Standard) peak-to-ladder alignment for FragalyseQt.
#
# Problem: given a list of known ladder sizes (bp) and a list of detected
# peak positions (data points / scans) in the ILS channel, find the optimal
# one-to-one assignment of each detected peak to its expected ladder size.
#
# The naive approach (take the last N peaks) fails when the actual CE run
# used a longer protocol than the selected size standard definition, or when
# noise produces extra peaks around the ladder.
#
# Solution: relative spacing pattern matching.
#   For electrophoresis the absolute dp→bp mapping is non-linear, but the
#   RELATIVE pattern of spacings between consecutive ladder fragments is
#   preserved.  Two consecutive fragments that are close in bp space will
#   also be close in dp space relative to fragments farther apart.
#
#   We search for the window of exactly N consecutive detected peaks whose
#   normalised inter-peak spacing pattern best matches the known size standard
#   spacing pattern.  The comparison is fully scale-invariant, so it works
#   regardless of run length, capillary length, or voltage.
#
# All algorithms are independent re-implementations based on algorithmic
# principles; no source code from fatoolsng (LGPLv3) is incorporated.
# This file is part of FragalyseQt and is distributed under AGPLv3.

from numpy import array, polyfit, poly1d, inf, nan, float64


# ---------------------------------------------------------------------------
# Utility helpers (used by tests and for polynomial sanity checks)
# ---------------------------------------------------------------------------

def _fit_poly(dp_pts, bp_pts, degree):
    # Fit the polynomial bp = poly(dp) of the given degree.
    # Automatically clamps degree so it never exceeds len(points)-1.
    actual_degree = min(degree, len(dp_pts) - 1)
    z = polyfit(array(dp_pts, dtype=float64),
                array(bp_pts, dtype=float64),
                actual_degree)
    return poly1d(z)


def _compute_rss(f, dp_pts, bp_pts):
    # Residual sum of squares for polynomial f evaluated at dp_pts vs bp_pts.
    dp = array(dp_pts, dtype=float64)
    bp = array(bp_pts, dtype=float64)
    return float(((f(dp) - bp) ** 2).sum())


# ---------------------------------------------------------------------------
# Primary alignment method: relative spacing pattern matching
# ---------------------------------------------------------------------------

def _pattern_match_align(dp_all, size_std):
    # Find the window of exactly len(size_std) CONSECUTIVE detected peaks
    # whose relative inter-peak spacing pattern best matches the spacing
    # pattern of the known size standard.
    #
    # Algorithm:
    #   1. Compute the normalised bp spacing pattern:
    #        bp_diffs[i] = size_std[i+1] - size_std[i]
    #        bp_pattern  = bp_diffs / (size_std[-1] - size_std[0])
    #
    #   2. For each candidate window of n_sizes consecutive detected peaks:
    #        dp_diffs[i] = dp_all[k+i+1] - dp_all[k+i]
    #        dp_pattern  = dp_diffs / (dp_all[k+n-1] - dp_all[k])
    #        score       = mean((dp_pattern - bp_pattern)^2)
    #
    #   3. Return the window with the lowest score.
    #
    # The normalisation makes the comparison scale-invariant: it works for any
    # run length or capillary, and for any size standard definition (60–460,
    # 60–500, 60–600 …) against the same actual electropherogram.
    #
    # Returns (matched_dp, best_score):
    #   matched_dp : 1-D ndarray of shape (n_sizes,) — dp positions of the
    #                best-matching window, directly usable as ladder_peaks
    #   best_score : mean squared deviation of the two normalised patterns
    #                (0 = perfect match, > 0.1 indicates a poor match)
    # Returns (None, inf) if the peak count is less than n_sizes.

    n_sz = len(size_std)
    n_dp = len(dp_all)

    if n_dp < n_sz:
        return None, inf

    sizes = array(sorted(size_std), dtype=float64)

    # normalised bp spacing pattern
    bp_diffs  = sizes[1:] - sizes[:-1]
    total_bp  = sizes[-1] - sizes[0]
    if total_bp <= 0:
        return None, inf
    bp_pattern = bp_diffs / total_bp

    best_score = inf
    best_k     = 0

    for k in range(n_dp - n_sz + 1):
        window   = dp_all[k:k + n_sz]
        total_dp = window[-1] - window[0]
        # skip degenerate windows (zero span)
        if total_dp <= 0:
            continue
        dp_diffs   = window[1:] - window[:-1]
        dp_pattern = dp_diffs / total_dp

        score = float(((dp_pattern - bp_pattern) ** 2).mean())

        if score < best_score:
            best_score = score
            best_k     = k

    return dp_all[best_k:best_k + n_sz].copy(), best_score


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------

def align_ils_peaks(dp_positions, size_std, heights):
    # Align detected ILS peak positions (data points) to known ladder sizes (bp).
    #
    # Parameters
    # ----------
    # dp_positions : 1-D array-like — detected ILS peak positions (data points)
    # size_std     : 1-D array-like — known ladder sizes in bp
    # heights      : 1-D array-like — peak heights (kept for API compatibility)
    #
    # Returns
    # -------
    # ndarray of shape (len(size_std),) — dp position matched to each
    # size_std[i], in ascending size order.
    # Returns an all-NaN array when fewer peaks are detected than expected sizes.

    dp_arr   = array(sorted(dp_positions), dtype=float64)
    size_arr = array(sorted(size_std),     dtype=float64)

    n_sizes = len(size_arr)
    n_peaks = len(dp_arr)

    if n_peaks < 2 or n_sizes < 2:
        return array([nan] * n_sizes)

    matched_dp, _ = _pattern_match_align(dp_arr, size_arr)

    if matched_dp is None:
        # fewer detected peaks than expected sizes — nothing we can do
        return array([nan] * n_sizes)

    return matched_dp
