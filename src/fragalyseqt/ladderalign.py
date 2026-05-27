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
# Problem: given a list of known ladder sizes (bp) and a list of detected peak
# positions (data points/scans) in the ILS channel, find the optimal one-to-one
# assignment of each detected peak to its expected ladder size.
#
# The naive approach (take the last N peaks) fails when the actual CE run used
# a longer protocol than the selected size standard definition, or when noise
# produces extra peaks around the ladder.
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

from numpy import array, polyfit, poly1d, inf, nan, float64, zeros, exp, int32
# ---------------------------------------------------------------------------
# Utility helpers (used by tests and for polynomial sanity checks)
# ---------------------------------------------------------------------------


def _fit_poly(dp_pts, bp_pts, degree):
    # Fit the polynomial bp = poly(dp) of the given degree.
    # Automatically clamps degree so it never exceeds len(points)-1.
    actual_degree = min(degree, len(dp_pts) - 1)
    z = polyfit(array(dp_pts, dtype=float64), array(bp_pts, dtype=float64),
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


def _detect_pre_ladder_boundary(dp_sorted, ratio_threshold=4.0,
                                max_prev_gap=100.0):
    # Scan consecutive-gap ratios in a SORTED dp array and locate the index
    # of the gap that marks the pre-ladder → ladder transition.  The
    # transition is characterised by:
    #   - gaps[i] / gaps[i-1] >= ratio_threshold (sudden widening)
    #   - gaps[i-1] <= max_prev_gap (previous gap is "pre-ladder-tight";
    #     once we're past the cluster, gaps are wider and this filter
    #     prevents the function from firing inside the ladder zone)
    #
    # Returns -1 if no qualifying transition is found, otherwise the index
    # i such that dp_sorted[i+1:] starts at the first ladder peak.  Within
    # the ladder zone the maximum gap-ratio for GS-600 sizes is ~3.3
    # (6 bp → 20 bp); 4.0 leaves a safety margin.  max_prev_gap caps the
    # previous gap to a magnitude characteristic of pre-ladder clusters
    # (dye blobs sit < 50 dp apart for typical CE conditions).
    n = len(dp_sorted)
    if n < 4:
        return -1
    gaps = [dp_sorted[i + 1] - dp_sorted[i] for i in range(n - 1)]
    for i in range(1, len(gaps)):
        if gaps[i - 1] <= 0 or gaps[i - 1] > max_prev_gap:
            continue
        if gaps[i] / gaps[i - 1] >= ratio_threshold:
            return i
    return -1


def _strip_pre_ladder_cluster(strong_dp, ratio_threshold=5.0):
    # Pre-ladder dye blobs (bleedthrough, primer dimers, unincorporated
    # dye terminators) often produce a tight cluster of TALL peaks before
    # the actual ILS ladder begins.  They appear in the strong-peak set
    # but are spaced much more densely than the ladder.  The transition
    # from the cluster to the first ladder peak is marked by a gap that
    # is several times larger than the within-cluster gaps.
    #
    # Detect that boundary: walk through consecutive strong-peak gaps and
    # find the first index where gap[i] / gap[i-1] >= ratio_threshold.
    # Drop strong peaks at indices 0..i and keep i+1 onwards.
    #
    # Within the ladder, the largest gap-ratio transitions are bp 6→20
    # (~3.3x in dp) — well below 5.0 — so a threshold of 5.0 catches the
    # pre-ladder boundary without ever firing inside the ladder zone.
    #
    # Returns the filtered strong_dp; if no qualifying jump is found,
    # returns the input unchanged.
    n = len(strong_dp)
    if n < 4:
        return strong_dp
    sd = sorted(strong_dp)
    gaps = [sd[i + 1] - sd[i] for i in range(n - 1)]
    big_gap_idx = -1
    for i in range(1, len(gaps)):
        if gaps[i - 1] <= 0:
            continue
        if gaps[i] / gaps[i - 1] >= ratio_threshold:
            big_gap_idx = i
            break
    if big_gap_idx < 0:
        return strong_dp
    # gaps[i] is the gap from sd[i] to sd[i+1]; keep sd[i+1:] (drop the
    # pre-ladder peaks up to and including sd[i]).
    keep_from = big_gap_idx + 1
    remaining = n - keep_from
    if remaining < 4:
        return strong_dp
    return sd[keep_from:]


def _find_best_alignment_subset(strong_dp, size_std, max_trim_start=8,
                                max_trim_end=4, min_k=4,
                                score_threshold=1e-3):
    # Find the smallest boundary-trim of `size_std` that fits the given
    # strong-peak positions with acceptable quality.  Tries combinations
    # (trim_s, trim_e) where trim_s ladder sizes are dropped from the start
    # and trim_e from the end.  For each (trim_s, trim_e), runs
    # `_pattern_match_align` against the corresponding sub-list of size_std
    # and accepts the result if its score is below `score_threshold`.
    #
    # Rationale: in real ILS runs the smallest ladder sizes (20/40/60) and
    # sometimes the largest (580/600) are absent or poorly resolved.  The
    # 36-anchor pattern matcher then forces those sizes onto pre-ladder
    # noise or end-of-run noise, shifting the whole window by 1-3 anchor
    # positions.  Running pattern match against strong peaks only AND
    # allowing the size_std to be boundary-trimmed lets the algorithm
    # choose the right subset (e.g. 80-600 if 20/40/60 are absent) and
    # leaves the trimmed sizes as NaN at the caller level.
    #
    # Among acceptable candidates we prefer the SMALLEST trim count.
    # Pattern-match scores tend to decrease as the window shrinks (smaller
    # mean residual with fewer intervals), so a pure minimum-score search
    # would over-trim and discard genuine boundary peaks.
    #
    # Returns (trim_start, trim_end, best_score) or None when no subset
    # with len >= min_k passes `score_threshold`.
    from numpy import asarray, float64
    n_sizes = len(size_std)
    n_strong = len(strong_dp)
    if n_strong < min_k or n_sizes < min_k:
        return None
    sd_arr = asarray(strong_dp, dtype=float64)
    best_score = inf
    best_result = None
    best_trim_count = inf
    for trim_s in range(0, min(max_trim_start, n_sizes - min_k) + 1):
        for trim_e in range(0, min(max_trim_end, n_sizes - trim_s - min_k) + 1):
            k_target = n_sizes - trim_s - trim_e
            if k_target > n_strong:
                continue
            sub_sizes = size_std[trim_s:n_sizes - trim_e]
            matched, score = _pattern_match_align(sd_arr, sub_sizes)
            if matched is None or score > score_threshold:
                continue
            trim_count = trim_s + trim_e
            # Prefer smallest trim count; break ties with lowest score.
            if (trim_count < best_trim_count
                    or (trim_count == best_trim_count and score < best_score)):
                best_trim_count = trim_count
                best_score = score
                best_result = (trim_s, trim_e)
    if best_result is None:
        return None
    return (best_result[0], best_result[1], best_score)


def _pattern_match_align(dp_all, size_std, heights=None):
    # Find the window of exactly len(size_std) CONSECUTIVE detected peaks whose
    # relative inter-peak spacing pattern best matches the spacing pattern of
    # the known size standard.  The `heights` argument is accepted but
    # currently unused — kept on the signature for future height-aware
    # scoring experiments.
    #
    # Algorithm:
    #   1. Compute the normalised bp spacing pattern:
    #        bp_diffs[i] = size_std[i+1] - size_std[i]
    #        bp_pattern = bp_diffs/(size_std[-1] - size_std[0])
    #
    #   2. For each candidate window of n_sizes consecutive detected peaks:
    #        dp_diffs[i] = dp_all[k+i+1] - dp_all[k+i]
    #        dp_pattern = dp_diffs/(dp_all[k+n-1] - dp_all[k])
    #        score = mean((dp_pattern - bp_pattern)^2)
    #
    #   3. Return the window with the lowest score.
    #
    # The normalisation makes the comparison scale-invariant: it works for any
    # run length or capillary, and for any size standard definition against
    # the same actual electropherogram.
    #
    # Returns (matched_dp, best_score):
    #   matched_dp : 1-D ndarray of shape (n_sizes,) — dp positions of the
    #                best-matching window, directly usable as ladder_peaks
    #   best_score : mean squared deviation of the two normalised patterns
    #                (0 = perfect match, > 0.1 indicates a poor match)
    # Returns (None, inf) if the peak count is less than n_sizes.
    _ = heights  # currently unused
    n_sz = len(size_std)
    n_dp = len(dp_all)
    if n_dp < n_sz:
        return None, inf
    sizes = array(sorted(size_std), dtype=float64)
    # normalised bp spacing pattern
    bp_diffs = sizes[1:] - sizes[:-1]
    total_bp = sizes[-1] - sizes[0]
    if total_bp <= 0:
        return None, inf
    bp_pattern = bp_diffs / total_bp
    best_score = inf
    best_k = 0
    for k in range(n_dp - n_sz + 1):
        window = dp_all[k:k + n_sz]
        total_dp = window[-1] - window[0]
        # skip degenerate windows (zero span)
        if total_dp <= 0:
            continue
        dp_diffs = window[1:] - window[:-1]
        dp_pattern = dp_diffs / total_dp
        score = float(((dp_pattern - bp_pattern) ** 2).mean())
        if score < best_score:
            best_score = score
            best_k = k
    return dp_all[best_k:best_k + n_sz].copy(), best_score
# ---------------------------------------------------------------------------
# Monotonicity check and correction
# ---------------------------------------------------------------------------


def _ror(window, sizes):
    # Ratio-of-ratios for each consecutive triple in the matched window:
    #
    #   ror[i] = (Δdp[i+1] / Δdp[i]) / (Δbp[i+1] / Δbp[i])
    #          = (Δdp[i+1] × Δbp[i]) / (Δdp[i] × Δbp[i+1])
    #
    # For a valid electrophoresis run this value changes slowly and stays
    # roughly in [0.35, 3.0]. A non-monotone peak (blob, undetected ILS
    # fragment replaced by a neighbour) produces a pair of values outside
    # this range: one anomalously small, the next anomalously large.
    d_dp = window[1:] - window[:-1]
    d_bp = sizes[1:] - sizes[:-1]
    return (d_dp[1:] * d_bp[:-1]) / (d_dp[:-1] * d_bp[1:])


def _violation_count(window, sizes, lo=0.35, hi=2.0):
    # Return the number of RoR values outside [lo, hi].
    # Returns (count, valid): valid=False when any Δdp ≤ 0 (non-monotone dp).
    d_dp = window[1:] - window[:-1]
    if float(d_dp.min()) <= 0:
        return 1000, False
    r = _ror(window, sizes)
    return int(((r < lo) | (r > hi)).sum()), True


def _fix_monotonicity(dp_all, sizes, window, lo=0.35, hi=2.0, max_iter=16):
    # Fix RoR violations by removing the offending (satellite/blob) peak from
    # dp_all and re-running pattern matching.
    #
    # For each RoR violation at triple (k, k+1, k+2), both window[k+1] and
    # window[k+2] are candidates for removal (skip anchors 0 and n-1).
    # All candidates are tried in every iteration; the one whose removal gives
    # the BEST pattern-match score (most faithful to the size-standard pattern)
    # among those that do not increase violations is accepted.
    #
    # Using pattern score as a tiebreaker resolves ambiguity when two adjacent
    # peaks are nearly equidistant from the expected position — the removal
    # that produces a perfectly spaced window (score≈0) wins over one that
    # leaves a slight residual spacing error.
    #
    # Guard: the new window must start at the same dp as the original (±1 dp)
    # to prevent left-shift where a pre-ladder noise peak is pulled in as the
    # size-range anchor, shifting every size assignment by one step.
    dp_arr = array(sorted(dp_all), dtype=float64)
    size_arr = array(sorted(sizes), dtype=float64)
    anchor_start = float(window[0])
    anchor_end = float(window[-1])
    visited = {tuple(window.tolist())}
    for _ in range(max_iter):
        n_viol, _ = _violation_count(window, size_arr, lo, hi)
        if n_viol == 0:
            break
        r = _ror(window, size_arr)
        n = len(window)
        # Collect suspect positions from all violated triples.
        # ror[k] involves positions k, k+1, k+2 — both inner positions are
        # candidates; anchors (0 and n-1) are never removed.
        suspects = set()
        for k, rv in enumerate(r):
            if lo <= float(rv) <= hi:
                continue
            if 0 < k + 1 < n - 1:
                suspects.add(k + 1)
            if 0 < k + 2 < n - 1:
                suspects.add(k + 2)
        # Try every suspect; keep the removal that gives the lowest pattern
        # score (best alignment) while not increasing violations.
        best_win = None
        best_red = None
        best_pscore = inf
        best_viol = n_viol
        for pos in sorted(suspects):
            suspect_dp = float(window[pos])
            reduced = dp_arr[dp_arr != suspect_dp]
            if len(reduced) < len(size_arr):
                continue
            new_win, pat_score = _pattern_match_align(reduced, size_arr)
            if new_win is None:
                continue
            n_new, ok = _violation_count(new_win, size_arr, lo, hi)
            if not ok or n_new > n_viol:
                continue
            # Start-anchor guard: never shift start to the left (would pull in
            # pre-ladder noise peaks that precede the genuine ladder region).
            if abs(float(new_win[0]) - anchor_start) > 1.0:
                continue
            # End-anchor guard: allow the window end to shift only when this
            # step strictly reduces violations (n_new < n_viol). A sideways
            # move that WORSENS violations is never allowed to shift the end,
            # preventing the "exactly-full" ILS channel from silently pulling
            # in an OL peak. Sideways moves (n_new == n_viol) with improved
            # pattern score are allowed: they occur when two independent blobs
            # happen to produce the same violation count before and after
            # removal (removing one blob exposes the other). anchor_end is
            # updated after each accepted step so that subsequent iterations
            # can build on the progress.
            if abs(float(new_win[-1]) - anchor_end) > 1.0:
                if n_new > n_viol:
                    continue
            new_key = tuple(new_win.tolist())
            if new_key in visited:
                continue
            # Among all valid candidates prefer the best pattern score.
            if pat_score < best_pscore:
                best_pscore = pat_score
                best_viol = n_new
                best_win = new_win
                best_red = reduced
        if best_win is None:
            # no useful removal found
            break
        visited.add(tuple(best_win.tolist()))
        window = best_win
        dp_arr = best_red
        # update so next step can build on it
        anchor_end = float(best_win[-1])
    # Return both the matched window and the cleaned dp_arr (no satellites).
    # The caller can use cleaned_dp_arr as the DP input so that already-removed
    # satellites do not reappear in the DP alignment.
    return window, dp_arr
# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def _backextrapolate_window_start(clean_dp, size_arr, window, lo=0.35, hi=2.0):
    # When the window starts too late (genuine ILS peaks before window[0] were
    # missed), RoR violations cluster at the beginning while a longer interior
    # segment is clean.  This function tries to recover the correct start.
    #
    # Strategy: find the longest RoR-monotone interior segment (best_s >= 1),
    # use its last peak as a second anchor, then run NW alignment with each
    # peak in clean_dp that lies BEFORE window[0] as the tentative first-size
    # anchor.  The NW polynomial is initialised with (candidate→size_arr[0],
    # segment_end→segment_end_bp) so it sees genuine ILS peaks as well as
    # pre-window satellites and naturally skips the satellites (extra-peak
    # moves are free).  Accept the result with the fewest violations.
    n = len(window)
    if n < 4:
        return window
    r = _ror(window, size_arr)
    n_viol_cur = int(sum(1 for rv in r if not (lo <= float(rv) <= hi)))
    if n_viol_cur == 0:
        return window
    # Longest monotone segment that does NOT start at position 0.
    best_s, best_l, cur_s, cur_l = 1, 1, 1, 1
    for i in range(1, len(r)):
        if lo <= float(r[i]) <= hi:
            cur_l += 1
        else:
            if cur_l > best_l:
                best_l, best_s = cur_l, cur_s
            cur_s, cur_l = i + 1, 1
    if cur_l > best_l:
        best_l, best_s = cur_l, cur_s
    if best_s == 0 or best_l < 3:
        return window
    # Second anchor: the last peak of the monotone segment.  It has the same
    # size assignment in the wrong window and the correct window (violations
    # only occur at the start, not the end).
    anchor2_dp = float(window[best_s + best_l - 1])
    anchor2_bp = float(size_arr[best_s + best_l - 1])
    bp0 = float(size_arr[0])

    best_result = window
    best_viol = n_viol_cur
    # Try every peak that comes before the current window start as a candidate
    # for the first ILS size (size_arr[0]).  For each candidate build a linear
    # polynomial through (candidate_dp, bp0) and (anchor2_dp, anchor2_bp),
    # run NW on the full clean_dp, and keep the result with fewest violations.
    w_start = float(window[0])
    for cand_dp_v in clean_dp:
        cand_dp = float(cand_dp_v)
        if cand_dp >= w_start:
            # clean_dp is sorted; no need to go further
            break
        if cand_dp >= anchor2_dp:
            continue
        slope = (anchor2_bp - bp0) / (anchor2_dp - cand_dp)
        init_poly = poly1d([slope, bp0 - slope * cand_dp])
        nw_matches, _ = _dp_align(clean_dp, size_arr, init_poly,
                                  gap_penalty=-0.5, tol_bp=15.0, max_iter=6)
        if len(nw_matches) != n:
            continue
        nw_win = array([float(clean_dp[j]) for _, j in nw_matches],
                       dtype=float64)
        n_new, ok = _violation_count(nw_win, size_arr)
        if ok and n_new < best_viol:
            best_viol = n_new
            best_result = nw_win
            if best_viol == 0:
                break
    return best_result


def _global_linear_init(dp_all, sizes):
    # Rough global linear bp=f(dp) polynomial for NW initialisation.
    #
    # The pre-ladder zone (bleedthrough, unincorporated primers) creates a
    # dense cluster of low-dp peaks that must not anchor the polynomial.
    # Heuristic: skip peaks until the first inter-peak gap that is at least
    # 130 % of the minimum expected ILS gap (min_bp_step × rough_rate).
    # This typically places dp_lo at the first genuine ILS peak.
    # dp_hi uses the 90th-percentile peak to avoid post-ILS noise.
    n = len(dp_all)
    size_arr = array(sorted(sizes), dtype=float64)
    bp_lo = float(size_arr[0])
    bp_hi = float(size_arr[-1])
    bp_range = bp_hi - bp_lo
    dp_range = float(dp_all[-1] - dp_all[0])
    if n < 2 or bp_range <= 0 or dp_range <= 0:
        return poly1d([0.0, bp_lo])
    rough_rate = dp_range / bp_range
    min_bp_step = float((size_arr[1:] - size_arr[:-1]).min())
    # 130 % of min expected gap
    gap_threshold = min_bp_step * rough_rate * 1.3
    dp_lo_idx = 0
    for i in range(n - 1):
        if float(dp_all[i + 1] - dp_all[i]) >= gap_threshold:
            dp_lo_idx = i
            break
    dp_hi_idx = min(n - 1, int(n * 0.90))
    dp_lo = float(dp_all[dp_lo_idx])
    dp_hi = float(dp_all[dp_hi_idx])
    if dp_hi <= dp_lo:
        return poly1d([0.0, bp_lo])
    slope = (bp_hi - bp_lo) / (dp_hi - dp_lo)
    return poly1d([slope, bp_lo - slope * dp_lo])


def _nw_align(dp_all, sizes, poly_func, gap_penalty=-0.5, tol_bp=8.0):
    # Needleman-Wunsch global alignment of sizes to detected peaks.
    #
    # Sizes run along rows, detected peaks along columns.
    # Moving diagonally = matching size[i] to peak[j]: score = S[i,j]
    # Moving up = size[i] has no matching peak: gap_penalty (< 0)
    # Moving left = peak[j] is an extra peak: free (0)
    #
    # This allows:
    #   - Missing sizes (e.g. 60bp below detection threshold) → unmatched,
    #     penalised
    #   - Extra peaks (satellites, dye blobs) → unmatched, free
    #
    # Returns list of (size_idx, peak_idx) matches in ascending order.
    n_sz = len(sizes)
    n_dp = len(dp_all)
    # Score S[i,j]: Gaussian similarity between predicted bp for peak j and
    # the expected size i.
    predicted_bp = poly_func(dp_all)
    S = zeros((n_sz, n_dp))
    for i, sz in enumerate(sizes):
        diff = (predicted_bp - sz) / tol_bp
        S[i] = exp(-diff * diff / 2.0)
    # DP matrix and traceback
    D = zeros((n_sz + 1, n_dp + 1))
    trace = zeros((n_sz + 1, n_dp + 1), dtype=int32)
    for i in range(1, n_sz + 1):
        D[i, 0] = D[i - 1, 0] + gap_penalty
        # forced up
        trace[i, 0] = 1
    for j in range(1, n_dp + 1):
        # forced left (free)
        trace[0, j] = 2
    # stop
    trace[0, 0] = 3
    for i in range(1, n_sz + 1):
        for j in range(1, n_dp + 1):
            diag = D[i - 1, j - 1] + S[i - 1, j - 1]
            up = D[i - 1, j] + gap_penalty
            left = D[i, j - 1] + 0.0
            if diag >= up and diag >= left:
                D[i, j] = diag
                trace[i, j] = 0
            elif up >= left:
                D[i, j] = up
                trace[i, j] = 1
            else:
                D[i, j] = left
                trace[i, j] = 2
    # Traceback
    matches = []
    i, j = n_sz, n_dp
    while trace[i, j] != 3:
        t = int(trace[i, j])
        if t == 0:
            i -= 1
            j -= 1
            matches.append((i, j))
        elif t == 1:
            i -= 1
        else:
            j -= 1
    matches.reverse()
    return matches, float(D[n_sz, n_dp])


def _dp_align(dp_all, sizes, initial_poly, gap_penalty=-0.5, tol_bp=8.0,
              max_iter=10):
    # Iterative DP: alternate between NW alignment and polynomial refit until
    # the DP score stops improving.
    #
    # Starts from the polynomial provided, then refines it through successive
    # alignments. After convergence returns (matched_pairs, final_polynomial).
    dp_arr = array(dp_all, dtype=float64)
    size_arr = array(sizes, dtype=float64)
    f = initial_poly
    prev_score = -inf
    best_matches = []
    for _ in range(max_iter):
        matches, score = _nw_align(dp_arr, size_arr, f, gap_penalty, tol_bp)
        if score <= prev_score:
            break
        best_matches = matches
        prev_score = score
        if not matches:
            break
        matched_dp = array([dp_arr[j] for _, j in matches], dtype=float64)
        matched_bp = array([size_arr[i] for i, _ in matches], dtype=float64)
        degree = min(3, len(matches) - 1)
        if degree < 1:
            break
        f = _fit_poly(matched_dp, matched_bp, degree)
    return best_matches, f


def align_ils_peaks(dp_positions, size_std, heights, saturated_threshold=32000):
    # Align detected ILS peak positions (dp) to known ladder sizes (bp).
    #
    # Parameters
    # ----------
    # dp_positions : 1-D array-like — detected ILS peak positions (data points)
    # size_std : 1-D array-like — known ladder sizes in bp
    # heights : 1-D array-like — peak heights in the same order as dp_positions
    # saturated_threshold : peaks with height > this value are oversaturated
    #                       and are excluded before alignment. Default 32000
    #                       matches the ABI 16-bit instrument ceiling.
    #
    # Returns
    # -------
    # ndarray of shape (len(size_std),) — dp position matched to each
    # size_std[i], in ascending size order.  Individual entries are NaN when
    # the algorithm cannot confidently match that size to a detected peak
    # (an anchor outside the longest RoR-monotone stretch); callers should
    # filter NaN out before using the result for sizing.  All-NaN is
    # returned when too few peaks are detected or no monotone stretch
    # exists at all.

    # Drop saturated peaks before doing anything else.
    dp_pos = array(dp_positions, dtype=float64)
    ht = array(heights, dtype=float64)
    if len(dp_pos) == len(ht):
        mask = ht <= saturated_threshold
        dp_pos = dp_pos[mask]
        ht = ht[mask]
    # Keep dp_arr and ht_arr aligned by dp ordering so the height-aware
    # snap step at the end can map an anchor's dp back to its height.
    # We also keep an UNTOUCHED copy of the full peak set (dp_arr_full,
    # ht_arr_full) because Step 1 deletes peaks it considers satellites,
    # and the snap step needs the original list to discover that the
    # "satellite" was actually a stronger, real peak.
    if len(ht) == len(dp_pos):
        from numpy import argsort
        _order = argsort(dp_pos)
        dp_arr = array(dp_pos[_order], dtype=float64)
        ht_arr = array(ht[_order], dtype=float64)
        dp_arr_full = dp_arr.copy()
        ht_arr_full = ht_arr.copy()
    else:
        dp_arr = array(sorted(dp_pos), dtype=float64)
        ht_arr = None
        dp_arr_full = None
        ht_arr_full = None
    full_size_arr = array(sorted(size_std), dtype=float64)
    n_full_sizes = len(full_size_arr)
    # Phase 1: strong-peak subset selection.  When the run is missing the
    # smallest (20/40/60) or largest (580/600) ladder fragments — because
    # the injection underloaded the low-bp end, or the run was cut short —
    # forcing the pattern matcher to anchor every size against weak
    # pre-ladder/post-ladder noise shifts the whole window.  By matching
    # only the strong peaks against boundary-trimmed subsets of size_std,
    # the algorithm picks the right alignment range; trimmed sizes are
    # filled with NaN below.
    trim_start = 0
    trim_end = 0
    strong_dp_final = None
    strong_score = inf
    if (ht_arr_full is not None and len(ht_arr_full) >= 4
            and n_full_sizes >= 6):
        from numpy import median as npmedian
        median_h = float(npmedian(ht_arr_full))
        if median_h > 0:
            strong_thresh = median_h * 3.0
            strong_mask = ht_arr_full >= strong_thresh
            strong_dp = dp_arr_full[strong_mask]
            strong_ht = ht_arr_full[strong_mask]
            # Iterative refinement: noise spikes occasionally exceed the
            # initial median*3 cut but sit far below the actual ladder peaks.
            # Re-median from within the strong subset and drop anything less
            # than 30% of that. Stabilises within one or two passes.
            for _ in range(3):
                if len(strong_ht) < 4:
                    break
                med_s = float(npmedian(strong_ht))
                if med_s <= 0:
                    break
                refined_mask = strong_ht >= med_s * 0.3
                if int(refined_mask.sum()) == len(strong_ht):
                    break
                strong_dp = strong_dp[refined_mask]
                strong_ht = strong_ht[refined_mask]
            # Strip pre-ladder dye-blob cluster: those blobs are often
            # TALLER than ladder peaks (so iterative refinement keeps them)
            # but sit much closer together.  A gap-ratio jump separates the
            # cluster from the first ladder peak.
            strong_dp = _strip_pre_ladder_cluster(strong_dp)
            if len(strong_dp) >= 4:
                subset = _find_best_alignment_subset(
                    strong_dp, list(sorted(size_std)))
                if subset is not None:
                    trim_start, trim_end, strong_score = subset
                    strong_dp_final = strong_dp
    size_arr = full_size_arr[trim_start:n_full_sizes - trim_end] \
        if (trim_start or trim_end) else full_size_arr
    n_sizes = len(size_arr)
    n_peaks = len(dp_arr)

    def _expand(sub_result):
        # Expand a sub-result aligned to size_arr (length n_sizes) back to
        # the full size_std order (length n_full_sizes), with NaN for the
        # boundary-trimmed sizes.
        if trim_start == 0 and trim_end == 0:
            return sub_result
        from numpy import full as npfull
        out = npfull(n_full_sizes, nan, dtype=float64)
        out[trim_start:trim_start + len(sub_result)] = sub_result
        return out

    # Strong-peak fast path: when the strong subset matched the trimmed
    # size_std with one-to-one correspondence and a low residual, that IS
    # the answer.  The full-peak pipeline below is for messy data with
    # satellites and noise — using it when we already have a clean fit
    # from the strong subset risks pulling the alignment off into nearby
    # noise.  Accept the strong-peak alignment when:
    #   - n_strong equals n_sub (every trimmed size matched a strong peak)
    #   - pattern-match score < 1e-3 (excellent fit, ~1% spacing residual)
    #   - resulting window has valid monotonic ratios
    if (strong_dp_final is not None and len(strong_dp_final) == n_sizes
            and strong_score < 1e-3):
        from numpy import asarray as npasarray
        candidate = npasarray(sorted(strong_dp_final), dtype=float64)
        viol, ok = _violation_count(candidate, size_arr)
        if ok and viol == 0:
            return _expand(candidate)

    if n_peaks < n_sizes or n_sizes < 2:
        return _expand(array([nan] * n_sizes))
    # Step 1: Iterative unconditional removal of distance-based satellites.
    #
    # After each pattern match, check every inner window position with:
    #   threshold = 2 * (window[i+1] - window[i-1]) / (sizes[i] - sizes[i-1])
    # A peak is a satellite when window[i] - window[i-1] < threshold, meaning
    # it sits much closer to its predecessor than the local dp/bp rate predicts.
    #
    # Identified satellites are removed from dp_all UNCONDITIONALLY (no anchor
    # guards, no violation-count comparison) and pattern matching is re-run.
    # This repeats until the window is stable (no satellites found).
    max_rounds = min(n_peaks - n_sizes, n_sizes // 2) + 1
    matched_dp = None
    for _ in range(max_rounds + 1):
        matched_dp, _ = _pattern_match_align(dp_arr, size_arr, ht_arr)
        if matched_dp is None:
            return _expand(array([nan] * n_sizes))
        satellites = []
        for i in range(1, len(matched_dp) - 1):
            # gap i→next
            d_bp_next = float(size_arr[i + 1]) - float(size_arr[i])
            if d_bp_next <= 0:
                continue
            # dp prev→next
            span = float(matched_dp[i + 1]) - float(matched_dp[i - 1])
            actual = float(matched_dp[i]) - float(matched_dp[i - 1])
            # Formula: threshold = 2*span/(sizes[i+1] - sizes[i]) A peak is a
            # satellite when its gap to the predecessor is less than
            # threshold, meaning it's closer than the local dp/bp rate allows.
            if actual < 2.0 * span / d_bp_next:
                satellites.append(float(matched_dp[i]))
        if not satellites:
            # stable — no satellites in current window
            break
        remove_set = set(satellites)
        keep_mask = array([float(x) not in remove_set for x in dp_arr],
                          dtype=bool)
        dp_arr = dp_arr[keep_mask]
        if ht_arr is not None:
            ht_arr = ht_arr[keep_mask]
        if len(dp_arr) < n_sizes:
            matched_dp, _ = _pattern_match_align(dp_arr, size_arr, ht_arr)
            if matched_dp is None:
                return _expand(array([nan] * n_sizes))
            break
    # Step 2: RoR-based fix for any remaining anomalies (conditional).
    matched_dp, clean_dp = _fix_monotonicity(dp_arr, size_arr, matched_dp)
    # Step 3: Left-shift exploration — if violations remain, try shifting the
    # window start leftward (exhaustive search, pick minimum violations).
    n_viol_cur, _ = _violation_count(matched_dp, size_arr)
    max_shifts = min(len(clean_dp) - n_sizes, n_sizes // 2)
    if n_viol_cur > 0 and max_shifts > 0:
        start_idx = int((clean_dp == matched_dp[0]).argmax())
        best_viol = n_viol_cur
        best_dp = matched_dp
        for shift in range(1, max_shifts + 1):
            new_start = start_idx - shift
            if new_start < 0:
                break
            candidate = clean_dp[new_start: new_start + n_sizes]
            if len(candidate) < n_sizes:
                break
            n_viol_new, ok = _violation_count(candidate, size_arr)
            if ok and n_viol_new < best_viol:
                best_viol = n_viol_new
                best_dp = candidate
                if best_viol == 0:
                    break
        matched_dp = best_dp
    # Step 4: Monotone-segment back-extrapolation.
    #
    # When the window starts at the wrong peak, the first few peaks may
    # carry RoR violations while a longer interior segment is clean.
    # Fit a polynomial to that clean segment, evaluate it at the first
    # ILS size, and look for a peak in clean_dp at the predicted position.
    # If found, try a shifted window from that peak; accept if violations
    # strictly decrease.
    matched_dp = _backextrapolate_window_start(clean_dp, size_arr, matched_dp)
    # Step 5: NW alignment with global linear initialisation.
    #
    # The pattern-matching window can be anchored at the wrong starting peak
    # when satellite blobs cause the correct window (with more satellites)
    # to score worse than a shifted wrong window (with fewer).  NW avoids
    # this by treating extra peaks as free (no penalty) and searching
    # globally with a polynomial that is not anchored to any single window.
    #
    # Acceptance criterion: the NW result must be complete (every size
    # matched to some peak) and have strictly fewer RoR violations than the
    # current best.  This keeps the change conservative — if NW finds
    # nothing better it is silently ignored.
    init_poly = _global_linear_init(clean_dp, size_arr)
    nw_matches, _ = _dp_align(clean_dp, size_arr, init_poly, gap_penalty=-0.5,
                              tol_bp=20.0, max_iter=8)
    if len(nw_matches) == n_sizes:
        nw_window = array([float(clean_dp[j]) for _, j in nw_matches],
                          dtype=float64)
        nw_viol, nw_ok = _violation_count(nw_window, size_arr)
        cur_viol, _ = _violation_count(matched_dp, size_arr)
        if nw_ok and nw_viol < cur_viol:
            matched_dp = nw_window
    # Step 6: confidence masking.  Replace anchors that fall outside the
    # longest RoR-monotone stretch with NaN.  Downstream sizing must drop
    # NaN entries before fitting — see fragalyseqt.py:align_ils_peaks call.
    matched_dp = _mask_low_confidence_anchors(matched_dp, size_arr)
    # Step 7: height-aware snap.  The pattern matcher ignores peak heights,
    # so it can anchor a size at a weak satellite when a much stronger real
    # peak sits a few dp away.  For each non-NaN anchor, look ±40 dp in the
    # ORIGINAL peak set (before Step 1 satellite removal — that step may
    # have stripped the real stronger peak) for a peak >=3× taller and
    # snap to it, but only if RoR remains valid in all triples involving
    # this position.  40 dp tolerance picks up satellites within typical
    # +A/-A artifact range and adjacent-noise peaks; RoR keeps us from
    # crossing into the next ladder anchor's slot.
    if ht_arr_full is not None:
        matched_dp = _snap_to_stronger_peaks(matched_dp, dp_arr_full,
                                             ht_arr_full, size_arr)
    return _expand(matched_dp)


def _snap_to_stronger_peaks(window, dp_arr, ht_arr, sizes,
                            tolerance=40.0, height_factor=3.0,
                            lo=0.35, hi=2.0):
    # For each non-NaN anchor in `window`, search the (dp_arr, ht_arr)
    # input for a peak within ±tolerance dp that is at least
    # `height_factor` times taller than the current anchor's peak.  If
    # found AND moving the anchor there does not introduce a RoR
    # violation in any triple involving that position, snap.
    #
    # The pattern matcher uses only normalised spacings, so when two
    # peaks sit close in dp it can pick the one that gives a slightly
    # better global score even if a much stronger peak (a real ladder
    # anchor) is 5–15 dp away.  This step is the final cleanup that
    # corrects those sub-pixel mis-anchors.
    from numpy import isnan, asarray, argmin, abs as npabs, where as npwhere
    n = len(window)
    result = window.copy()
    if n < 3:
        return result
    sz_arr = asarray(sizes, dtype=float64)
    dp_a = asarray(dp_arr, dtype=float64)
    ht_a = asarray(ht_arr, dtype=float64)
    for i in range(n):
        d = result[i]
        if isnan(d):
            continue
        # Locate this anchor's peak in dp_arr.
        idx = int(argmin(npabs(dp_a - d)))
        if npabs(dp_a[idx] - d) > 0.5:
            continue  # anchor not in input — should not happen
        cur_h = float(ht_a[idx])
        # dps already used by other anchors (don't steal them)
        used_dps = set()
        for k in range(n):
            if k != i and not isnan(result[k]):
                used_dps.add(float(result[k]))
        # Find the strongest nearby peak above the snap threshold.
        threshold_h = cur_h * height_factor
        best_dp = d
        best_h = threshold_h
        in_range = npwhere((dp_a >= d - tolerance) &
                           (dp_a <= d + tolerance))[0]
        for j in in_range:
            cand_dp = float(dp_a[j])
            cand_h = float(ht_a[j])
            if cand_dp == d or cand_dp in used_dps:
                continue
            if cand_h > best_h:
                best_h = cand_h
                best_dp = cand_dp
        if best_dp == d:
            continue
        # Verify RoR remains valid in every triple involving position i.
        test = result.copy()
        test[i] = best_dp
        ok = True
        for k in range(max(0, i - 2), min(n - 2, i) + 1):
            if k + 2 >= n:
                continue
            if isnan(test[k]) or isnan(test[k+1]) or isnan(test[k+2]):
                continue
            d_dp_a = test[k+1] - test[k]
            d_dp_b = test[k+2] - test[k+1]
            d_bp_a = sz_arr[k+1] - sz_arr[k]
            d_bp_b = sz_arr[k+2] - sz_arr[k+1]
            if d_dp_a <= 0:
                ok = False
                break
            ror = (d_dp_b * d_bp_a) / (d_dp_a * d_bp_b)
            if not (lo <= ror <= hi):
                ok = False
                break
        if ok:
            result[i] = best_dp
    return result


def _mask_low_confidence_anchors(window, sizes, lo=0.35, hi=2.0, min_run=2):
    # Return a copy of `window` (float64) with NaN at every position that
    # is not part of the longest contiguous run of valid RoR triples.
    #
    # A RoR is valid when it lies in [lo, hi].  A run of L valid RoRs
    # implies L+2 mutually consistent anchors; that range is kept, every
    # other anchor is NaN-ed.  When no run reaches `min_run` valid RoRs in
    # a row, every anchor is NaN-ed — the caller treats that as "no
    # confident ladder" and falls back.
    #
    # This is the gate that lets `align_ils_peaks` report partial
    # alignments: sizes that could not be confidently matched to a peak
    # come back as NaN instead of forced onto a satellite.
    from numpy import nan
    n = len(window)
    result = array(window, dtype=float64).copy()
    if n < 3:
        return result
    r = _ror(array(window, dtype=float64), array(sizes, dtype=float64))
    valid = (r >= lo) & (r <= hi)
    best_start, best_len = 0, 0
    cur_start, cur_len = 0, 0
    for i in range(len(valid)):
        if bool(valid[i]):
            cur_len += 1
            if cur_len > best_len:
                best_len = cur_len
                best_start = cur_start
        else:
            cur_start = i + 1
            cur_len = 0
    if best_len < min_run:
        result[:] = nan
        return result
    keep_lo = best_start
    keep_hi = best_start + best_len + 1
    for i in range(n):
        if i < keep_lo or i > keep_hi:
            result[i] = nan
    return result
