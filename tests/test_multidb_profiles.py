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

# Reference profiles live in TWO databases: the casework database and the
# dedicated reference-profile one. These tests cover storing profiles in
# either, listing and searching across both, and telling apart two profiles
# that share an id because they sit in different databases.

import pytest

from fragalyseqt.database import RefProfileBackend, SQLiteBackend
from fragalyseqt.comparison import search_profiles_multi
from fragalyseqt.forensicstats import AlleleCall
from fragalyseqt.refprofile import (
    ReferenceProfile,
    get_profile,
    list_profiles,
    list_profiles_multi,
    store_profile,
)

_CALLS = [
    AlleleCall('D3S1358', '15', '17'),
    AlleleCall('vWA', '18', None),
    AlleleCall('FGA', '22', '24'),
]


@pytest.fixture
def casework(tmp_path):
    backend = SQLiteBackend(str(tmp_path / 'casework.db'))
    yield backend
    backend.close()


@pytest.fixture
def refdb(tmp_path):
    backend = RefProfileBackend(str(tmp_path / 'ref.db'))
    yield backend
    backend.close()


@pytest.fixture
def sources(casework, refdb):
    return [('casework', 'Casework', casework),
            ('refdb', 'Reference profiles', refdb)]


def _store(db, name, role='Suspect, Known', calls=None):
    return store_profile(db, ReferenceProfile(
        name=name, role=role, notes=None,
        calls=list(_CALLS if calls is None else calls)))


class TestCaseworkDatabaseStoresProfiles:
    # Saving to the casework database used to raise AttributeError: the
    # backend simply had no reference-profile methods.

    def test_store_and_read_back(self, casework):
        pid = _store(casework, 'Suspect A')
        assert isinstance(pid, int) and pid > 0
        got = get_profile(casework, pid)
        assert got.name == 'Suspect A'
        assert [c.marker for c in got.calls] == [c.marker for c in _CALLS]

    def test_listed_after_store(self, casework):
        _store(casework, 'Suspect A')
        assert [p['name'] for p in list_profiles(casework)] == ['Suspect A']

    def test_casework_features_still_work(self, casework):
        # The schema gained tables; it must not have lost any.
        assert casework.get_session_list() == []
        assert hasattr(casework, 'store_file')

    def test_both_backends_share_the_same_api(self, casework, refdb):
        for name in ('store_reference_profile', 'get_reference_profile',
                     'list_reference_profiles', 'get_reference_alleles',
                     'store_reference_alleles', 'delete_reference_profile'):
            assert hasattr(casework, name), name
            assert hasattr(refdb, name), name


class TestListAcrossDatabases:
    def test_profiles_from_both_databases(self, sources, casework, refdb):
        _store(casework, 'From a case')
        _store(refdb, 'Reference one')
        names = {p['name'] for p in list_profiles_multi(sources)}
        assert names == {'From a case', 'Reference one'}

    def test_each_row_knows_its_source(self, sources, casework, refdb):
        _store(casework, 'From a case')
        _store(refdb, 'Reference one')
        by_name = {p['name']: p for p in list_profiles_multi(sources)}
        assert by_name['From a case']['source'] == 'casework'
        assert by_name['Reference one']['source'] == 'refdb'
        assert by_name['From a case']['source_label'] == 'Casework'

    def test_same_id_in_two_databases_stays_distinct(self, sources,
                                                     casework, refdb):
        # Both are id 1: only the source tells them apart.
        assert _store(casework, 'Case one') == _store(refdb, 'Ref one')
        rows = list_profiles_multi(sources)
        assert len(rows) == 2
        assert len({(r['source'], r['id']) for r in rows}) == 2

    def test_empty_databases_give_no_rows(self, sources):
        assert list_profiles_multi(sources) == []

    def test_missing_database_is_skipped(self, casework):
        _store(casework, 'From a case')
        rows = list_profiles_multi([('casework', 'Casework', casework),
                                    ('refdb', 'Reference profiles', None)])
        assert [r['name'] for r in rows] == ['From a case']


class TestSearchAcrossDatabases:
    def test_finds_a_profile_stored_in_casework(self, sources, casework):
        # The search used to look in the reference database only.
        _store(casework, 'From a case')
        hits = search_profiles_multi(sources, _CALLS)
        assert [h['name'] for h in hits] == ['From a case']
        assert hits[0]['source'] == 'casework'

    def test_ranks_globally_across_databases(self, sources, casework, refdb):
        _store(refdb, 'Partial', calls=[AlleleCall('D3S1358', '15', '17'),
                                        AlleleCall('vWA', '19', '20')])
        _store(casework, 'Exact')
        hits = search_profiles_multi(sources, _CALLS)
        # The exact match wins even though it lives in the other database.
        assert hits[0]['name'] == 'Exact'
        assert hits[0]['matched'] == hits[0]['common']

    def test_hits_carry_their_database(self, sources, casework, refdb):
        _store(casework, 'From a case')
        _store(refdb, 'Reference one')
        for hit in search_profiles_multi(sources, _CALLS):
            assert hit['source'] in ('casework', 'refdb')
            assert hit['db'] is not None


class TestAppendOnly:
    # The casework database gained the reference_* tables, so the append-only
    # guarantee has to hold there too: updates supersede, deletes tombstone,
    # and nothing is ever physically overwritten or removed.

    def _count(self, backend, table):
        cur = backend._conn.execute(f'SELECT COUNT(*) FROM {table}')
        return cur.fetchone()[0]

    @pytest.mark.parametrize('which', ['casework', 'refdb'])
    def test_update_supersedes_instead_of_overwriting(self, which, casework,
                                                      refdb):
        db = casework if which == 'casework' else refdb
        from fragalyseqt.refprofile import update_profile
        old_id = _store(db, 'Original')
        new_id = update_profile(db, old_id, name='Corrected')
        assert new_id != old_id
        assert self._count(db, 'reference_profile') == 2
        # The superseded row is hidden from listings but still readable.
        assert [p['name'] for p in list_profiles(db)] == ['Corrected']
        assert get_profile(db, old_id).name == 'Original'

    @pytest.mark.parametrize('which', ['casework', 'refdb'])
    def test_delete_writes_a_tombstone_and_keeps_the_row(self, which, casework,
                                                         refdb):
        db = casework if which == 'casework' else refdb
        from fragalyseqt.refprofile import delete_profile
        pid = _store(db, 'To be deleted')
        delete_profile(db, pid)
        assert self._count(db, 'reference_profile') == 1
        assert self._count(db, 'reference_profile_deletion') == 1
        assert list_profiles(db) == []
        # Still physically present for audit purposes.
        assert get_profile(db, pid).name == 'To be deleted'

    def test_module_issues_no_destructive_sql(self):
        # The SQL this module runs lives in string literals, so inspect those
        # directly. Docstrings are excluded: prose such as "no UPDATE or
        # DELETE involved" describes the guarantee rather than breaking it.
        import ast
        import pathlib
        import re

        import fragalyseqt.database as dbmod

        tree = ast.parse(
            pathlib.Path(dbmod.__file__).read_text(encoding='utf-8'))

        docstrings = set()
        for node in ast.walk(tree):
            if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef,
                                 ast.AsyncFunctionDef)):
                body = getattr(node, 'body', None)
                if (body and isinstance(body[0], ast.Expr)
                        and isinstance(body[0].value, ast.Constant)
                        and isinstance(body[0].value.value, str)):
                    docstrings.add(id(body[0].value))

        sql = ' '.join(
            node.value for node in ast.walk(tree)
            if isinstance(node, ast.Constant)
            and isinstance(node.value, str)
            and id(node) not in docstrings)

        for bad in (r'\bUPDATE\s+\w', r'\bDELETE\s+FROM\b',
                    r'\bREPLACE\s+INTO\b', r'\bDROP\s+TABLE\b',
                    r'\bTRUNCATE\b'):
            assert not re.search(bad, sql, re.IGNORECASE), bad

        # Sanity: the extracted text really is the module's SQL, otherwise
        # the assertions above would pass against an empty string.
        assert 'INSERT INTO reference_profile' in sql
