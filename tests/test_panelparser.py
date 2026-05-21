import os
import pytest
import numpy as np
from fragalyseqt.panelparser import (
    parse_genemapper, parse_osiris,
    assign_alleles, has_bin_data, has_stutter_data,
)

PANELS_DIR = os.path.join(os.path.dirname(__file__), '..', 'docs', 'NIST_DATASET_PANELS')
TFS_DIR = os.path.join(PANELS_DIR, 'Thermo_Fisher_Scientific')
PPG_DIR = os.path.join(PANELS_DIR, 'Promega')
OSIRIS_DIR = os.path.join(os.path.dirname(__file__), '..', 'docs', 'OSIRIS_PANELS')

AMPFLSTR_PANELS = os.path.join(TFS_DIR, 'AmpFLSTR_Panels_v5X_Panels.txt')
AMPFLSTR_BINS = os.path.join(TFS_DIR, 'AmpFLSTR_Panels_v5X_AmpFLSTR_Bins_v5X_bins.txt')
AMPFLSTR_STUTTER = os.path.join(TFS_DIR, 'AmpFLSTR_Panels_v5X_stutter.txt')
PPF6C_PANELS = os.path.join(PPG_DIR, 'PowerPlex_Fusion_6C_Panels_IDX_v1.2_Panels.txt')
PPF6C_BINS = os.path.join(PPG_DIR, 'PowerPlex_Fusion_6C_Panels_IDX_v1.2_PowerPlex_Fusion_6C_Bins_IDX_v1.2_bins.txt')
PPY23_PANELS = os.path.join(PPG_DIR, 'PowerPlexY23_Panels_IDX_v1.2_Panels.txt')
PPY23_BINS = os.path.join(PPG_DIR, 'PowerPlexY23_Panels_IDX_v1.2_PowerPlexY23_Bins_IDX_v1.2_bins.txt')
PPY23_STUTTER = os.path.join(PPG_DIR, 'PowerPlexY23_Panels_IDX_v1.2_stutter.txt')
GLOBALFILER_XML = os.path.join(OSIRIS_DIR, 'GlobalFiler_LadderInfo.xml')


# ---------------------------------------------------------------------------
# GeneMapper — panels file only
# ---------------------------------------------------------------------------

def test_genemapper_panels_parses():
    panels = parse_genemapper(AMPFLSTR_PANELS)
    assert len(panels) > 0


def test_genemapper_panels_have_markers():
    panels = parse_genemapper(AMPFLSTR_PANELS)
    for panel in panels.values():
        assert len(panel) > 0, "Panel must have at least one marker"


def test_genemapper_panels_no_bin_data_without_bins_file():
    panels = parse_genemapper(AMPFLSTR_PANELS)
    assert not has_bin_data(panels), "No bin data expected without bins file"


def test_genemapper_panels_no_stutter_without_stutter_file():
    panels = parse_genemapper(AMPFLSTR_PANELS)
    assert not has_stutter_data(panels)


# ---------------------------------------------------------------------------
# GeneMapper — with bins
# ---------------------------------------------------------------------------

def test_genemapper_bins_loaded():
    panels = parse_genemapper(AMPFLSTR_PANELS, AMPFLSTR_BINS)
    assert has_bin_data(panels), "Bin data expected after loading bins file"


def test_genemapper_bins_sizes_are_floats():
    panels = parse_genemapper(AMPFLSTR_PANELS, AMPFLSTR_BINS)
    for panel in panels.values():
        for marker in panel.values():
            for allele in marker['alleles']:
                assert isinstance(allele['size'], float)
                assert isinstance(allele['left_bin'], float)
                assert isinstance(allele['right_bin'], float)


# ---------------------------------------------------------------------------
# GeneMapper — with stutter
# ---------------------------------------------------------------------------

def test_genemapper_stutter_loaded():
    panels = parse_genemapper(AMPFLSTR_PANELS, AMPFLSTR_BINS, AMPFLSTR_STUTTER)
    assert has_stutter_data(panels)


def test_genemapper_stutter_ratios_positive():
    panels = parse_genemapper(PPY23_PANELS, PPY23_BINS, PPY23_STUTTER)
    for panel in panels.values():
        for marker in panel.values():
            for key in ('minus', 'plus'):
                v = marker['stutter'][key]
                if v is not None:
                    assert v >= 0, f"Stutter ratio must be >= 0, got {v}"


# ---------------------------------------------------------------------------
# OSIRIS
# ---------------------------------------------------------------------------

def test_osiris_parses():
    panels = parse_osiris(GLOBALFILER_XML)
    assert len(panels) > 0


def test_osiris_markers_have_min_max():
    panels = parse_osiris(GLOBALFILER_XML)
    for panel in panels.values():
        for marker in panel.values():
            assert marker['min_size'] is not None
            assert marker['max_size'] is not None
            assert marker['min_size'] < marker['max_size']


def test_osiris_default_bins_are_positive():
    panels = parse_osiris(GLOBALFILER_XML)
    for panel in panels.values():
        for marker in panel.values():
            for allele in marker['alleles']:
                assert allele['left_bin'] > 0
                assert allele['right_bin'] > 0


# ---------------------------------------------------------------------------
# assign_alleles
# ---------------------------------------------------------------------------

def _simple_panel():
    return {
        'D3S1358': {
            'dye': 'blue',
            'min_size': 100.0,
            'max_size': 160.0,
            'stutter': {'minus': None, 'plus': None},
            'alleles': [
                {'label': '14', 'size': 120.0, 'left_bin': 0.5, 'right_bin': 0.5, 'virtual': False},
                {'label': '15', 'size': 124.0, 'left_bin': 0.5, 'right_bin': 0.5, 'virtual': False},
                {'label': '16', 'size': 128.0, 'left_bin': 0.5, 'right_bin': 0.5, 'virtual': True},
            ]
        }
    }


def test_assign_alleles_match():
    result = assign_alleles([120.1, 124.3], [1, 1], _simple_panel())
    assert result[0] == 'D3S1358:14'
    assert result[1] == 'D3S1358:15'


def test_assign_alleles_virtual_suffix():
    result = assign_alleles([128.0], [1], _simple_panel())
    assert result[0] == 'D3S1358:16*'


def test_assign_alleles_outside_range():
    result = assign_alleles([80.0], [1], _simple_panel())
    assert result[0] == ''


def test_assign_alleles_in_range_no_bin_match():
    # In range but outside all bin windows → OL
    result = assign_alleles([115.0], [1], _simple_panel())
    assert result[0] == 'D3S1358:OL'


def test_assign_alleles_wrong_channel():
    # Channel 2 = green, panel expects blue (ch 1) → no match
    result = assign_alleles([120.0], [2], _simple_panel())
    assert result[0] == ''


def test_assign_alleles_no_bins():
    panel = {
        'D3S1358': {
            'dye': 'blue', 'min_size': 100.0, 'max_size': 160.0,
            'stutter': {'minus': None, 'plus': None},
            'alleles': [{'label': '14', 'size': None,
                         'left_bin': None, 'right_bin': None, 'virtual': False}]
        }
    }
    result = assign_alleles([120.0], [1], panel)
    assert result[0] == 'D3S1358:?'


def test_assign_alleles_empty():
    result = assign_alleles([], [], _simple_panel())
    assert result == []
