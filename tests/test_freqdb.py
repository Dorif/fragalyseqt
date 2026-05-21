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

import os
import pytest
from fragalyseqt.freqdb import (
    FrequencyTable,
    import_freq_csv,
    import_freq_fam,
    save_freq_table,
    load_freq_table,
    list_freq_tables,
    get_allele_freq,
)

NIST_CSV = os.path.join(os.path.dirname(__file__), '..', 'docs',
           'SPECS_AND_REFERENCES', 'NIST_1036_Combined_1036.csv')


def test_import_csv_returns_frequency_table():
    t = import_freq_csv(NIST_CSV, 'test', '', 'Combined')
    assert isinstance(t, FrequencyTable)


def test_import_csv_has_markers():
    t = import_freq_csv(NIST_CSV, 'test', '', 'Combined')
    assert len(t.markers) > 0
    for marker, alleles in t.markers.items():
        assert isinstance(marker, str) and marker
        assert len(alleles) > 0


def test_import_csv_freq_sum_per_marker():
    t = import_freq_csv(NIST_CSV, 'test', '', 'Combined')
    for marker, alleles in t.markers.items():
        assert abs(sum(alleles.values()) - 1.0) < 0.001, marker


def test_import_csv_microvariant_key():
    t = import_freq_csv(NIST_CSV, 'test', '', 'Combined')
    assert '15.2' in t.markers['D3S1358']


def test_import_csv_integer_key_no_dot_zero():
    t = import_freq_csv(NIST_CSV, 'test', '', 'Combined')
    assert '12' in t.markers['CSF1PO']
    assert '12.0' not in t.markers['CSF1PO']


def test_import_csv_tab_delimiter(tmp_path):
    tsv = tmp_path / 'test.tsv'
    tsv.write_text('marker\tallele\tfrequency\nD3S1358\t15\t0.5\nD3S1358\t16\t0.5\n')
    t = import_freq_csv(str(tsv), 'test', '', 'pop')
    assert '15' in t.markers['D3S1358']


def test_import_csv_bad_frequency_raises(tmp_path):
    bad = tmp_path / 'bad.csv'
    bad.write_text('marker,allele,frequency\nCSF1PO,12,1.5\n')
    with pytest.raises(ValueError):
        import_freq_csv(str(bad), 'test', '', 'pop')


def test_import_csv_missing_column_raises(tmp_path):
    bad = tmp_path / 'bad.csv'
    bad.write_text('marker,allele\nCSF1PO,12\n')
    with pytest.raises(ValueError):
        import_freq_csv(str(bad), 'test', '', 'pop')


def test_save_load_roundtrip(tmp_path):
    t = import_freq_csv(NIST_CSV, 'NIST Combined', 'GlobalFiler', 'Combined')
    path = str(tmp_path / 'table.json')
    save_freq_table(t, path)
    t2 = load_freq_table(path)
    assert t2.name == t.name
    assert t2.markers == t.markers
    assert t2.theta == t.theta


def test_list_freq_tables(tmp_path):
    t = import_freq_csv(NIST_CSV, 'test', '', 'Combined')
    save_freq_table(t, str(tmp_path / 'freqtables' / 'a.json'))
    save_freq_table(t, str(tmp_path / 'freqtables' / 'b.json'))
    assert len(list_freq_tables(str(tmp_path))) == 2


def test_get_allele_freq_known():
    t = import_freq_csv(NIST_CSV, 'test', '', 'Combined')
    assert get_allele_freq(t, 'CSF1PO', '12') > 0.3


def test_get_allele_freq_missing_allele_returns_min_freq():
    t = import_freq_csv(NIST_CSV, 'test', '', 'Combined')
    assert get_allele_freq(t, 'CSF1PO', '99') == t.min_freq


def test_get_allele_freq_missing_marker_returns_min_freq():
    t = import_freq_csv(NIST_CSV, 'test', '', 'Combined')
    assert get_allele_freq(t, 'NONEXISTENT', '12') == t.min_freq


def test_get_allele_freq_normalizes_input():
    t = import_freq_csv(NIST_CSV, 'test', '', 'Combined')
    assert get_allele_freq(t, 'CSF1PO', '12.0') == get_allele_freq(t, 'CSF1PO', '12')


NORWEGIAN_FAM = os.path.join(os.path.dirname(__file__), '..', 'Norwegian_DB.fam')
SOUTHAM_FAM = os.path.join(os.path.dirname(__file__), '..', 'popSTR_South_America_DB.fam')


def test_import_fam_returns_frequency_table():
    t = import_freq_fam(NORWEGIAN_FAM, 'Norwegian')
    assert isinstance(t, FrequencyTable)


def test_import_fam_population_from_file():
    t = import_freq_fam(NORWEGIAN_FAM, 'Norwegian')
    assert t.population == 'Norwegian'


def test_import_fam_marker_count_norwegian():
    t = import_freq_fam(NORWEGIAN_FAM, 'Norwegian')
    assert len(t.markers) == 35


def test_import_fam_marker_count_southam():
    t = import_freq_fam(SOUTHAM_FAM, 'SouthAm')
    assert len(t.markers) == 70


def test_import_fam_freq_sum_per_marker():
    t = import_freq_fam(NORWEGIAN_FAM, 'Norwegian')
    for marker, alleles in t.markers.items():
        assert abs(sum(alleles.values()) - 1.0) < 0.001, marker


def test_import_fam_microvariant_key():
    t = import_freq_fam(NORWEGIAN_FAM, 'Norwegian')
    assert '8.3' in t.markers['TH01']


def test_import_fam_integer_key_no_dot_zero():
    t = import_freq_fam(NORWEGIAN_FAM, 'Norwegian')
    assert '14' in t.markers['D3S1358']
    assert '14.0' not in t.markers['D3S1358']


def test_import_fam_population_override():
    t = import_freq_fam(NORWEGIAN_FAM, 'Norwegian', population='Nordic')
    assert t.population == 'Nordic'


def test_import_fam_empty_file_raises(tmp_path):
    bad = tmp_path / 'empty.fam'
    bad.write_text('')
    with pytest.raises(ValueError):
        import_freq_fam(str(bad), 'test')
