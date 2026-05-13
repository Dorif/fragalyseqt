import numpy as np
import pytest
from fragalyseqt.jbcd import jbcd


def _flat_signal(n=500):
    rng = np.random.default_rng(42)
    baseline = np.linspace(100, 300, n)
    peaks = np.zeros(n)
    for pos in [100, 200, 300, 400]:
        x = np.arange(n)
        peaks += 2000 * np.exp(-0.5 * ((x - pos) / 8) ** 2)
    noise = rng.normal(0, 5, n)
    return baseline + peaks + noise


def test_output_shapes():
    signal = _flat_signal()
    baseline, params = jbcd(signal)
    assert baseline.shape == signal.shape
    assert params['signal'].shape == signal.shape


def test_baseline_does_not_exceed_signal():
    signal = np.abs(_flat_signal()) + 10
    baseline, _ = jbcd(signal)
    assert np.all(baseline <= signal + 1e-9)


def test_denoised_signal_plus_baseline_equals_input():
    signal = _flat_signal()
    baseline, params = jbcd(signal)
    np.testing.assert_allclose(params['signal'] + baseline, signal, rtol=1e-10)


def test_baseline_below_peaks():
    n = 300
    x = np.arange(n, dtype=float)
    baseline_true = 50 + 0.1 * x
    peak = 3000 * np.exp(-0.5 * ((x - 150) / 10) ** 2)
    signal = baseline_true + peak
    baseline_est, _ = jbcd(signal, half_window=30)
    # At the peak apex the baseline should be well below the peak
    assert baseline_est[150] < signal[150] * 0.5


def test_flat_signal_baseline_near_signal():
    # Flat signal with no peaks — baseline ≈ signal
    signal = np.ones(200) * 500.0
    baseline, params = jbcd(signal, half_window=10)
    np.testing.assert_allclose(baseline, signal, rtol=0.01)


def test_accepts_list_input():
    signal = list(range(1, 201))
    baseline, params = jbcd(signal)
    assert baseline.shape == (200,)


@pytest.mark.parametrize('diff_order', [1, 2, 3])
def test_diff_orders(diff_order):
    signal = _flat_signal(300)
    baseline, params = jbcd(signal, diff_order=diff_order)
    assert baseline.shape == signal.shape
    assert np.all(baseline <= signal + 1e-9)
