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

"""
Standalone implementation of the Joint Baseline Correction and Denoising (jBCD)
algorithm (Liu et al., Applied Spectroscopy 2015, 69(9), 1013-1022).

Ported from pybaselines.morphological (Donald Erb, 2021) with all framework
scaffolding removed.
Original pybaselines package is distributed under BSD 3-clause license.
"""

from numpy import full as npfull
from numpy import ones as npones
from numpy import finfo as npfinfo
from numpy import asarray as npasarray
from numpy import minimum as npminimum
from numpy.linalg import norm
from scipy.linalg import solveh_banded
from scipy.ndimage import grey_opening, grey_dilation, grey_erosion

_MIN_FLOAT = npfinfo(float).tiny


# ---------------------------------------------------------------------------
# Banded penalty matrix construction  (D^T D in lower-only banded format)
# ---------------------------------------------------------------------------

def _diff_penalty(n, diff_order):
    """
    Lower-only banded form of the diff-order finite-difference penalty matrix.
    Returns array of shape (diff_order+1, n) suitable for solveh_banded with
    lower=True.  Hard-coded for orders 1–3; sparse fallback for higher orders.
    """
    if diff_order == 1:
        out = npfull((2, n), -1.)
        out[0, 0] = out[0, -1] = 1.
        out[0, 1:-1] = 2.
        out[1, -1] = 0.
        return out

    if diff_order == 2:
        out = npones((3, n))
        out[0, 1] = out[0, -2] = 5.
        out[0, 2:-2] = 6.
        out[1, 0] = out[1, -2] = -2.
        out[1, 1:-2] = -4.
        out[1, -1] = 0.
        out[2, -1] = out[2, -2] = 0.
        return out

    if diff_order == 3:
        out = npfull((4, n), -1.)
        out[0, 0] = out[0, -1] = 1.
        out[0, 1] = out[0, -2] = 10.
        out[0, 2] = out[0, -3] = 19.
        out[0, 3:-3] = 20.
        out[1, 0] = out[1, -2] = -3.
        out[1, 1] = out[1, -3] = -12.
        out[1, 2:-3] = -15.
        out[1, -1] = 0.
        out[2, 0] = out[2, -3] = 3.
        out[2, 1:-3] = 6.
        out[2, -1] = out[2, -2] = 0.
        out[3, -1] = out[3, -2] = out[3, -3] = 0.
        return out

    # General case via sparse
    from scipy.sparse import eye as speye
    D = speye(n - 1, n, 1) - speye(n - 1, n)
    for _ in range(diff_order - 1):
        m = D.shape[0] - 1
        D = (speye(m, m + 1, 1) - speye(m, m + 1)) @ D
    P = (D.T @ D).todia()
    diags = P.data[::-1]          # full banded, main diagonal first
    return diags[diff_order:].copy()   # keep lower half only


# ---------------------------------------------------------------------------
# Morphological helper
# ---------------------------------------------------------------------------

def _avg_opening(y, half_window, opening):
    # Element-wise average of dilation and erosion of the morphological opening
    w = 2 * half_window + 1
    return 0.5 * (grey_dilation(opening, w) + grey_erosion(opening, w))


# ---------------------------------------------------------------------------
# Public function
# ---------------------------------------------------------------------------

def jbcd(data, half_window=25, alpha=0.1, beta=10., gamma=1.,
         beta_mult=1.1, gamma_mult=0.909, diff_order=1,
         max_iter=20, tol=1e-2, tol_2=1e-3, robust_opening=True):
    """
    Joint Baseline Correction and Denoising.

    Parameters
    ----------
    data : array-like, shape (N,)
        Raw signal data.
    half_window : int, optional
        Half-size of the morphological structuring element.  Default is 25.
    alpha : float, optional
        Regularisation weight tying the baseline to the morphological opening.
        Default is 0.1.
    beta : float, optional
        Smoothness penalty for the baseline (grows by beta_mult each iteration).
        Default is 10.
    gamma : float, optional
        Smoothness penalty for the signal (shrinks by gamma_mult each iteration).
        Default is 1.
    beta_mult : float, optional
        Multiplicative update for beta per iteration.  Default is 1.1.
    gamma_mult : float, optional
        Multiplicative update for gamma per iteration.  Default is 0.909.
    diff_order : int, optional
        Order of the finite-difference penalty matrix.  Default is 1.
    max_iter : int, optional
        Maximum number of iterations.  Default is 20.
    tol : float, optional
        Convergence threshold for the signal.  Default is 1e-2.
    tol_2 : float, optional
        Convergence threshold for the baseline.  Default is 1e-3.
    robust_opening : bool, optional
        If True (default), use the element-wise minimum of the plain opening
        and the average of dilation/erosion of the opening as the initial
        baseline estimate.

    Returns
    -------
    baseline : numpy.ndarray, shape (N,)
        Estimated baseline.
    params : dict
        ``{'signal': ndarray}`` — denoised signal (data minus baseline).

    References
    ----------
    Liu, H., et al. Joint Baseline-Correction and Denoising for Raman Spectra.
    Applied Spectroscopy, 2015, 69(9), 1013-1022.
    """
    y = npasarray(data, dtype=float)
    n = y.size
    half_window = int(half_window)
    w = 2 * half_window + 1

    penalty = _diff_penalty(n, diff_order)
    # shape (diff_order+1, n), lower-only

    opening = grey_opening(y, w)
    if robust_opening:
        opening = npminimum(opening, _avg_opening(y, half_window, opening))

    baseline_old = opening.copy()
    signal_old = y.copy()
    partial_rhs_2 = 2. * alpha * opening

    for _ in range(max_iter + 1):
        lhs_1 = gamma * penalty.copy()
        lhs_1[0] += 1.

        lhs_2 = 2. * beta * penalty.copy()
        lhs_2[0] += 1. + 2. * alpha

        signal = solveh_banded(lhs_1, y - baseline_old,
                               lower=True, overwrite_ab=True, overwrite_b=True)
        baseline = solveh_banded(lhs_2, y - signal + partial_rhs_2,
                                 lower=True, overwrite_ab=True,
                                 overwrite_b=True)

        d_signal = norm(signal - signal_old) / max(norm(signal_old),
                                                   _MIN_FLOAT)
        d_baseline = norm(baseline - baseline_old) / max(norm(baseline_old),
                                                         _MIN_FLOAT)

        if d_signal < tol and d_baseline < tol_2:
            break

        signal_old = signal
        baseline_old = baseline
        gamma *= gamma_mult
        beta *= beta_mult

    return baseline, {'signal': signal}
