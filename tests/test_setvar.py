import numpy as np
import pytest
from fragalyseqt.setvar import (
    southern_fit_local, southern_fit_global,
    set_spl_dgr, set_lsq_ord, set_knots, chk_key_valid,
)


# ---------------------------------------------------------------------------
# southern_fit_local
# ---------------------------------------------------------------------------

def _ladder():
    """Synthetic ladder following the Southern hyperbolic model exactly:
    bp = c / (dp - m0) + L0  with c=1e6, m0=100, L0=50.
    """
    c, m0, L0 = 1e6, 100.0, 50.0
    ladder_peaks = np.array([500, 800, 1200, 1800, 2500, 3500], dtype=float)
    size_std = c / (ladder_peaks - m0) + L0
    return ladder_peaks, size_std


def test_southern_local_at_ladder_points():
    lp, ss = _ladder()
    result = southern_fit_local(lp, ss, lp)
    np.testing.assert_allclose(result, ss, atol=1.0)


def test_southern_local_interpolation():
    lp, ss = _ladder()
    # Query between first two points — result should lie between their sizes
    query = np.array([(lp[0] + lp[1]) / 2])
    result = southern_fit_local(lp, ss, query)
    assert ss[1] < result[0] < ss[0]


def test_southern_local_output_shape():
    lp, ss = _ladder()
    query = np.linspace(lp[0], lp[-1], 50)
    result = southern_fit_local(lp, ss, query)
    assert result.shape == query.shape


# ---------------------------------------------------------------------------
# southern_fit_global
# ---------------------------------------------------------------------------

def test_southern_global_at_ladder_points():
    lp, ss = _ladder()
    result = southern_fit_global(lp, ss, lp)
    np.testing.assert_allclose(result, ss, atol=2.0)


def test_southern_global_output_shape():
    lp, ss = _ladder()
    query = np.linspace(lp[0], lp[-1], 100)
    result = southern_fit_global(lp, ss, query)
    assert result.shape == query.shape


# ---------------------------------------------------------------------------
# set_spl_dgr
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('method,expected', [
    ('cubic spline',          3),
    ('Cubic Spline',          3),
    ('linear spline',         1),
    ('5th degree spline',     5),
    ('LSQ weighted cubic',    3),
    ('spline (no keyword)',   3),
])
def test_set_spl_dgr(method, expected):
    assert set_spl_dgr(method) == expected


# ---------------------------------------------------------------------------
# set_lsq_ord
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('method,expected', [
    ('2nd order polynomial', 2),
    ('3rd order polynomial', 3),
    ('5th order polynomial', 5),
    ('something else',       5),
])
def test_set_lsq_ord(method, expected):
    assert set_lsq_ord(method) == expected


# ---------------------------------------------------------------------------
# set_knots
# ---------------------------------------------------------------------------

def test_set_knots_weighted_returns_subarray():
    arr = np.arange(20, dtype=float)
    knots = set_knots('LSQ weighted cubic', arr, spldeg=3)
    assert knots is not None
    assert len(knots) > 0
    assert knots[0] >= arr[0]
    assert knots[-1] <= arr[-1]


def test_set_knots_non_weighted_returns_none():
    arr = np.arange(20, dtype=float)
    assert set_knots('cubic spline', arr, spldeg=3) is None


# ---------------------------------------------------------------------------
# chk_key_valid
# ---------------------------------------------------------------------------

def test_chk_key_valid_present():
    assert chk_key_valid('a', {'a': 1}) is True


def test_chk_key_valid_none_value():
    assert chk_key_valid('a', {'a': None}) is False


def test_chk_key_valid_missing():
    assert chk_key_valid('x', {'a': 1}) is False
