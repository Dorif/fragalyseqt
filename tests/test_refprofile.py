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
from fragalyseqt.database import SQLiteBackend
from fragalyseqt.forensicstats import AlleleCall
from fragalyseqt.refprofile import (
    ReferenceProfile,
    store_profile,
    update_profile,
    delete_profile,
    get_profile,
    list_profiles,
)

_CALLS = [
    AlleleCall('D3S1358', '15', '17'),
    AlleleCall('vWA', '18', None),
    AlleleCall('FGA', '22', '24'),
]


@pytest.fixture
def db(tmp_path):
    backend = SQLiteBackend(str(tmp_path / 'test.db'))
    yield backend
    backend.close()


def _profile(**kwargs):
    defaults = dict(name='Test', role='Suspect, Known',
                    notes=None, calls=list(_CALLS))
    defaults.update(kwargs)
    return ReferenceProfile(**defaults)


def test_store_returns_int(db):
    pid = store_profile(db, _profile())
    assert isinstance(pid, int) and pid > 0


def test_get_profile_name(db):
    pid = store_profile(db, _profile(name='Victim A'))
    assert get_profile(db, pid).name == 'Victim A'


def test_get_profile_role(db):
    pid = store_profile(db, _profile(role='Biological Father'))
    assert get_profile(db, pid).role == 'Biological Father'


def test_get_profile_notes(db):
    pid = store_profile(db, _profile(notes='Case 2024-001'))
    assert get_profile(db, pid).notes == 'Case 2024-001'


def test_alleles_roundtrip(db):
    pid = store_profile(db, _profile())
    calls = get_profile(db, pid).calls
    assert len(calls) == 3
    markers = {c.marker for c in calls}
    assert markers == {'D3S1358', 'vWA', 'FGA'}


def test_homo_allele2_is_none(db):
    pid = store_profile(db, _profile(calls=[AlleleCall('vWA', '18', None)]))
    calls = get_profile(db, pid).calls
    assert calls[0].allele2 is None


def test_list_shows_stored(db):
    store_profile(db, _profile(name='A'))
    store_profile(db, _profile(name='B'))
    names = {p['name'] for p in list_profiles(db)}
    assert names == {'A', 'B'}


def test_update_returns_new_id(db):
    pid = store_profile(db, _profile())
    new_pid = update_profile(db, pid, name='Updated')
    assert new_pid != pid


def test_update_name(db):
    pid = store_profile(db, _profile(name='Old'))
    new_pid = update_profile(db, pid, name='New')
    assert get_profile(db, new_pid).name == 'New'


def test_update_copies_alleles(db):
    pid = store_profile(db, _profile())
    new_pid = update_profile(db, pid, role='Victim, Known')
    assert len(get_profile(db, new_pid).calls) == 3


def test_update_replaces_alleles(db):
    pid = store_profile(db, _profile())
    new_calls = [AlleleCall('CSF1PO', '11', '12')]
    new_pid = update_profile(db, pid, calls=new_calls)
    calls = get_profile(db, new_pid).calls
    assert len(calls) == 1
    assert calls[0].marker == 'CSF1PO'


def test_update_hides_old_in_list(db):
    pid = store_profile(db, _profile(name='Old'))
    update_profile(db, pid, name='New')
    names = [p['name'] for p in list_profiles(db)]
    assert 'New' in names
    assert 'Old' not in names


def test_update_old_still_retrievable(db):
    pid = store_profile(db, _profile(name='Old'))
    update_profile(db, pid, name='New')
    assert get_profile(db, pid).name == 'Old'


def test_delete_hides_from_list(db):
    pid = store_profile(db, _profile())
    delete_profile(db, pid)
    assert list_profiles(db) == []


def test_delete_row_still_retrievable(db):
    pid = store_profile(db, _profile(name='Preserved'))
    delete_profile(db, pid)
    assert get_profile(db, pid).name == 'Preserved'


def test_list_newest_first(db):
    pid1 = store_profile(db, _profile(name='First'))
    pid2 = store_profile(db, _profile(name='Second'))
    profiles = list_profiles(db)
    assert profiles[0]['id'] == pid2
    assert profiles[1]['id'] == pid1


def test_multiple_updates_only_latest_in_list(db):
    pid = store_profile(db, _profile(name='v1'))
    pid2 = update_profile(db, pid, name='v2')
    update_profile(db, pid2, name='v3')
    names = [p['name'] for p in list_profiles(db)]
    assert names == ['v3']


def test_profiles_from_codis_xml(tmp_path):
    from fragalyseqt.refprofile import profiles_from_codis_xml
    NS = 'urn:CODISImportFile-schema'
    xml = f"""<?xml version="1.0"?>
<CODISImportFile xmlns="{NS}">
  <SPECIMEN>
    <SPECIMENID>CASE-001</SPECIMENID>
    <SPECIMENCATEGORY>Victim, Known</SPECIMENCATEGORY>
    <LOCUS>
      <LOCUSNAME>D3S1358</LOCUSNAME>
      <ALLELE><ALLELEVALUE>15</ALLELEVALUE></ALLELE>
      <ALLELE><ALLELEVALUE>17</ALLELEVALUE></ALLELE>
    </LOCUS>
    <LOCUS>
      <LOCUSNAME>vWA</LOCUSNAME>
      <ALLELE><ALLELEVALUE>18</ALLELEVALUE></ALLELE>
    </LOCUS>
  </SPECIMEN>
  <SPECIMEN>
    <SPECIMENID>CASE-002</SPECIMENID>
    <LOCUS>
      <LOCUSNAME>FGA</LOCUSNAME>
      <ALLELE><ALLELEVALUE>22</ALLELEVALUE></ALLELE>
      <ALLELE><ALLELEVALUE>24</ALLELEVALUE></ALLELE>
    </LOCUS>
  </SPECIMEN>
</CODISImportFile>"""
    path = tmp_path / 'test.xml'
    path.write_text(xml)
    profiles = profiles_from_codis_xml(str(path))
    assert len(profiles) == 2
    assert profiles[0].name == 'CASE-001'
    assert profiles[0].role == 'Victim, Known'
    assert len(profiles[0].calls) == 2
    d3 = next(c for c in profiles[0].calls if c.marker == 'D3S1358')
    assert d3.allele1 == '15' and d3.allele2 == '17'
    vwa = next(c for c in profiles[0].calls if c.marker == 'vWA')
    assert vwa.allele2 is None
    assert profiles[1].name == 'CASE-002'
    assert profiles[1].role is None
