import numpy as np
import pytest
from fragalyseqt.stutterfilter import apply_stutter_filter, _try_allele


# ---------------------------------------------------------------------------
# _try_allele
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('label,expected', [
    ('14',            14.0),
    ('D3S1358:14',    14.0),
    ('9.3',           9.3),
    ('SE33:9.3',      9.3),
    ('D3S1358:14*',   14.0),
    ('D3S1358:14.0',  14.0),
])
def test_try_allele_valid(label, expected):
    assert _try_allele(label) == expected


@pytest.mark.parametrize('label', [
    'D3S1358:OL', 'ILS', '', 'OL', 'D3S1358:?',
])
def test_try_allele_invalid(label):
    assert _try_allele(label) is None


# ---------------------------------------------------------------------------
# apply_stutter_filter
# ---------------------------------------------------------------------------

def _panel(minus=0.15, plus=None):
    return {
        'D3S1358': {
            'dye': 'blue',
            'min_size': 90.0,
            'max_size': 160.0,
            'stutter': {'minus': minus, 'plus': plus},
            'alleles': [],
        }
    }


def test_stutter_marked_as_empty_string():
    # allele 13 height is 5% of allele 14 → stutter (below 15% threshold)
    result = apply_stutter_filter(
        peaksizes=np.array([120.0, 116.0]),
        peakheights=np.array([1000.0, 50.0]),
        peakchannels=np.array([1, 1]),
        peakalleles=['D3S1358:14', 'D3S1358:13'],
        panel=_panel(minus=0.15),
        dye_names=['blue'],
    )
    assert result[0] == 'D3S1358:14'
    assert result[1] == ''


def test_balanced_heterozygote_not_stutter():
    # Both alleles roughly equal in height → neither is stutter
    result = apply_stutter_filter(
        peaksizes=np.array([120.0, 116.0]),
        peakheights=np.array([1000.0, 900.0]),
        peakchannels=np.array([1, 1]),
        peakalleles=['D3S1358:14', 'D3S1358:13'],
        panel=_panel(minus=0.15),
        dye_names=['blue'],
    )
    assert result[0] == 'D3S1358:14'
    assert result[1] == 'D3S1358:13'


def test_no_stutter_data_skips_marker():
    result = apply_stutter_filter(
        peaksizes=np.array([120.0, 116.0]),
        peakheights=np.array([1000.0, 50.0]),
        peakchannels=np.array([1, 1]),
        peakalleles=['D3S1358:14', 'D3S1358:13'],
        panel=_panel(minus=None, plus=None),
        dye_names=['blue'],
    )
    # No stutter data → no filtering
    assert result == ['D3S1358:14', 'D3S1358:13']


def test_zero_stutter_threshold_filters_nothing():
    # Threshold 0: h / parent < 0 is never true → nothing filtered
    result = apply_stutter_filter(
        peaksizes=np.array([120.0, 116.0]),
        peakheights=np.array([1000.0, 1.0]),
        peakchannels=np.array([1, 1]),
        peakalleles=['D3S1358:14', 'D3S1358:13'],
        panel=_panel(minus=0.0),
        dye_names=['blue'],
    )
    assert result[0] == 'D3S1358:14'
    assert result[1] == 'D3S1358:13'


def test_plus_stutter():
    # allele 15 height is 2% of allele 14 → n+1 stutter (below 5% threshold)
    result = apply_stutter_filter(
        peaksizes=np.array([120.0, 124.0]),
        peakheights=np.array([1000.0, 20.0]),
        peakchannels=np.array([1, 1]),
        peakalleles=['D3S1358:14', 'D3S1358:15'],
        panel=_panel(minus=None, plus=0.05),
        dye_names=['blue'],
    )
    assert result[0] == 'D3S1358:14'
    assert result[1] == ''


def test_ils_and_ol_not_touched():
    result = apply_stutter_filter(
        peaksizes=np.array([120.0, 116.0, 200.0]),
        peakheights=np.array([1000.0, 50.0, 5000.0]),
        peakchannels=np.array([1, 1, 5]),
        peakalleles=['D3S1358:14', 'D3S1358:13', 'ILS'],
        panel=_panel(minus=0.15),
        dye_names=['blue', 'green', 'yellow', 'red', 'orange'],
    )
    assert result[2] == 'ILS'


def test_ol_stutter_below_threshold_is_suppressed():
    result = apply_stutter_filter(
        peaksizes=np.array([120.0, 114.5, 116.0]),
        peakheights=np.array([1000.0, 60.0, 80.0]),
        peakchannels=np.array([1, 1, 1]),
        peakalleles=['D3S1358:14', 'D3S1358:OL', 'D3S1358:13'],
        panel=_panel(minus=0.15),
        dye_names=['blue'],
    )
    assert result[1] == '', "OL peak below stutter threshold should be suppressed"


def test_ol_stutter_above_threshold_is_kept():
    result = apply_stutter_filter(
        peaksizes=np.array([120.0, 114.5, 116.0]),
        peakheights=np.array([1000.0, 300.0, 80.0]),
        peakchannels=np.array([1, 1, 1]),
        peakalleles=['D3S1358:14', 'D3S1358:OL', 'D3S1358:13'],
        panel=_panel(minus=0.15),
        dye_names=['blue'],
    )
    assert result[1] == 'D3S1358:OL', "OL peak above stutter threshold should remain"


def test_input_not_mutated():
    alleles = ['D3S1358:14', 'D3S1358:13']
    apply_stutter_filter(
        peaksizes=np.array([120.0, 116.0]),
        peakheights=np.array([1000.0, 50.0]),
        peakchannels=np.array([1, 1]),
        peakalleles=alleles,
        panel=_panel(minus=0.15),
        dye_names=['blue'],
    )
    assert alleles == ['D3S1358:14', 'D3S1358:13']
