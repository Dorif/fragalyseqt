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

# Tests for the profile comparison dialog: the preselected frequency table and
# comparing profiles straight from the reference database with no files open.

import os

import pytest

pytest.importorskip('pyqtgraph')

os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')

from fragalyseqt.database import RefProfileBackend  # noqa: E402
from fragalyseqt.forensicstats import AlleleCall  # noqa: E402
from fragalyseqt.refprofile import ReferenceProfile, store_profile  # noqa: E402

_FREQTABLES_DIR = os.path.join(os.path.dirname(__file__), '..', 'src',
                               'fragalyseqt', 'freqtables')

_CALLS = [
    AlleleCall('D3S1358', '15', '17'),
    AlleleCall('vWA', '18', None),
    AlleleCall('FGA', '22', '24'),
]


@pytest.fixture(scope='module')
def qapp():
    from pyqtgraph.Qt.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


@pytest.fixture
def iface():
    # localizefq() picks the language from $LANG, so the dialog strings differ
    # per machine. The dialog only needs the keys to exist, and a plain dict
    # subclass keeps the test independent of the ambient locale.
    class _Msg(dict):
        def __missing__(self, key):
            return key

    from fragalyseqt.localize import localizefq
    msg = _Msg()
    localizefq(msg)
    return msg


@pytest.fixture
def db(tmp_path):
    backend = RefProfileBackend(str(tmp_path / 'ref.db'))
    yield backend
    backend.close()


def _store(db, name, role):
    return store_profile(db, ReferenceProfile(name=name, role=role,
                                              notes=None, calls=list(_CALLS)))


def _dialog(iface, tab_names=(), states=(), db=None, sources=None):
    from fragalyseqt.comparisondialog import ComparisonDialog
    return ComparisonDialog(list(states), list(tab_names), _FREQTABLES_DIR,
                            iface, db=db, sources=sources)


def _selected(combo):
    return combo.itemData(combo.currentIndex())


class TestDefaultFrequencyTable:
    def test_strider_entire_database_is_preselected(self, qapp, iface):
        from fragalyseqt.comparisondialog import _DEFAULT_TABLE_NAME
        dlg = _dialog(iface)
        assert dlg._table_combo.count() > 1
        assert dlg._table_combo.currentText() == _DEFAULT_TABLE_NAME

    def test_default_table_name_matches_a_shipped_table(self, qapp, iface):
        # Guards against the constant drifting away from the JSON "name".
        from fragalyseqt.comparisondialog import _DEFAULT_TABLE_NAME
        dlg = _dialog(iface)
        assert _DEFAULT_TABLE_NAME in dlg._tables

    def test_default_table_is_usable_for_calculation(self, qapp, iface):
        dlg = _dialog(iface)
        table = dlg._tables.get(dlg._table_combo.currentText())
        assert table is not None
        assert table.markers


class TestCompareWithoutOpenFiles:
    def test_saved_profiles_offered_with_no_files_open(self, qapp, iface, db):
        _store(db, 'Suspect A', 'Suspect, Known')
        _store(db, 'Victim B', 'Victim, Known')
        dlg = _dialog(iface, db=db)
        assert _selected(dlg._prof1_combo) is not None
        assert _selected(dlg._prof2_combo) is not None

    def test_two_combos_do_not_preselect_the_same_profile(self, qapp, iface,
                                                          db):
        # Otherwise the dialog opens comparing a profile against itself.
        _store(db, 'Suspect A', 'Suspect, Known')
        _store(db, 'Victim B', 'Victim, Known')
        dlg = _dialog(iface, db=db)
        assert _selected(dlg._prof1_combo) != _selected(dlg._prof2_combo)

    def test_no_leading_separator_when_no_files_open(self, qapp, iface, db):
        # A separator as the first row would leave itemData() == None and no
        # profile selected at all.
        _store(db, 'Suspect A', 'Suspect, Known')
        dlg = _dialog(iface, db=db)
        assert dlg._prof1_combo.itemData(0) is not None

    def test_calls_resolve_from_database_profiles(self, qapp, iface, db):
        _store(db, 'Suspect A', 'Suspect, Known')
        _store(db, 'Victim B', 'Victim, Known')
        dlg = _dialog(iface, db=db)
        calls = dlg._get_calls_for_combo(dlg._prof1_combo)
        assert [c.marker for c in calls] == [c.marker for c in _CALLS]

    def test_identity_comparison_runs_without_files(self, qapp, iface, db):
        from fragalyseqt.comparison import compare_identity
        _store(db, 'Suspect A', 'Suspect, Known')
        _store(db, 'Victim B', 'Victim, Known')
        dlg = _dialog(iface, db=db)
        table = dlg._tables[dlg._table_combo.currentText()]
        result = compare_identity(
            dlg._get_calls_for_combo(dlg._prof1_combo),
            dlg._get_calls_for_combo(dlg._prof2_combo),
            table, 0.01)
        # Both profiles carry identical calls, so loci must actually be
        # compared and the identity LR must favour the same-source hypothesis.
        assert result.n_loci > 0
        assert result.combined_stat > 1.0

    def test_open_files_still_take_precedence(self, qapp, iface, db):
        # Regression guard: with files open the first two tabs stay the
        # default, exactly as before saved profiles could be preselected.
        _store(db, 'Suspect A', 'Suspect, Known')
        dlg = _dialog(iface, tab_names=('FileA', 'FileB'),
                      states=({}, {}), db=db)
        assert _selected(dlg._prof1_combo) == ('tab', 0)
        assert _selected(dlg._prof2_combo) == ('tab', 1)

    def test_separator_present_when_both_groups_exist(self, qapp, iface, db):
        _store(db, 'Suspect A', 'Suspect, Known')
        dlg = _dialog(iface, tab_names=('FileA',), states=({},), db=db)
        datas = [dlg._prof1_combo.itemData(i)
                 for i in range(dlg._prof1_combo.count())]
        assert ('tab', 0) in datas
        assert None in datas  # the separator between files and profiles

    def test_dialog_builds_with_neither_files_nor_database(self, qapp, iface):
        dlg = _dialog(iface)
        assert dlg._prof1_combo.count() == 0
        assert dlg._table_combo.count() > 0


class TestBothDatabasesInDialog:
    # Profiles live in two databases; the dialog must offer both and resolve
    # each entry against the database it actually came from.

    @pytest.fixture
    def casework(self, tmp_path):
        from fragalyseqt.database import SQLiteBackend
        backend = SQLiteBackend(str(tmp_path / 'casework.db'))
        yield backend
        backend.close()

    @pytest.fixture
    def sources(self, casework, db):
        return [('casework', 'Casework', casework),
                ('refdb', 'Reference profiles', db)]

    def test_profiles_from_both_databases_are_offered(self, qapp, iface,
                                                      casework, db, sources):
        _store(casework, 'From case', 'Suspect, Known')
        _store(db, 'From reference', 'Victim, Known')
        dlg = _dialog(iface, sources=sources)
        texts = [dlg._prof1_combo.itemText(i)
                 for i in range(dlg._prof1_combo.count())]
        assert any('From case' in t for t in texts)
        assert any('From reference' in t for t in texts)

    def test_entries_are_labelled_with_their_database(self, qapp, iface,
                                                      casework, db, sources):
        _store(casework, 'From case', 'Suspect, Known')
        _store(db, 'From reference', 'Victim, Known')
        dlg = _dialog(iface, sources=sources)
        texts = [dlg._prof1_combo.itemText(i)
                 for i in range(dlg._prof1_combo.count())]
        assert any('Casework' in t for t in texts)
        assert any('Reference profiles' in t for t in texts)

    def test_item_data_carries_the_source_key(self, qapp, iface,
                                              casework, db, sources):
        _store(casework, 'From case', 'Suspect, Known')
        _store(db, 'From reference', 'Victim, Known')
        dlg = _dialog(iface, sources=sources)
        datas = [dlg._prof1_combo.itemData(i)
                 for i in range(dlg._prof1_combo.count())]
        assert ('profile', ('casework', 1)) in datas
        assert ('profile', ('refdb', 1)) in datas

    def test_same_id_in_both_databases_resolves_to_the_right_profile(
            self, qapp, iface, casework, db, sources):
        # Ids are unique only WITHIN a database: both profiles below are id 1.
        # Losing the source key would silently compare the wrong person.
        other = [AlleleCall('D3S1358', '14', '16'),
                 AlleleCall('vWA', '19', '20')]
        store_profile(casework, ReferenceProfile(
            name='Case person', role='Suspect, Known', notes=None,
            calls=list(other)))
        _store(db, 'Reference person', 'Victim, Known')
        dlg = _dialog(iface, sources=sources)

        by_source = {}
        for i in range(dlg._prof1_combo.count()):
            data = dlg._prof1_combo.itemData(i)
            if data is not None:
                by_source[data[1][0]] = i

        dlg._prof1_combo.setCurrentIndex(by_source['casework'])
        case_calls = dlg._get_calls_for_combo(dlg._prof1_combo)
        dlg._prof1_combo.setCurrentIndex(by_source['refdb'])
        ref_calls = dlg._get_calls_for_combo(dlg._prof1_combo)

        assert [c.allele1 for c in case_calls] == ['14', '19']
        assert [c.allele1 for c in ref_calls] == ['15', '18', '22']

    def test_no_database_label_when_only_one_source(self, qapp, iface, db):
        # A single source needs no disambiguating suffix.
        _store(db, 'Solo', 'Suspect, Known')
        dlg = _dialog(iface, sources=[('refdb', 'Reference profiles', db)])
        assert '·' not in dlg._prof1_combo.itemText(0)

    def test_plain_db_argument_still_works(self, qapp, iface, db):
        # Backward compatibility: callers passing db= keep working.
        _store(db, 'Solo', 'Suspect, Known')
        dlg = _dialog(iface, db=db)
        assert dlg._prof1_combo.count() == 1
        assert dlg._get_calls_for_combo(dlg._prof1_combo)


class TestLaunchWithoutOpenFiles:
    # compare_profiles() used to `return` silently unless files were open, so
    # profiles saved in the database could never be compared on their own.

    @pytest.fixture
    def window(self, qapp, tmp_path, monkeypatch):
        import fragalyseqt.fragalyseqt as fq
        from fragalyseqt.comparisondialog import ComparisonDialog
        monkeypatch.setattr(fq, '_USER_DATA', str(tmp_path))
        monkeypatch.setattr(fq, '_REFPROFILE_DB_PATH',
                            str(tmp_path / 'refprofiles.db'))
        seen = []
        monkeypatch.setattr(fq, 'msgbox',
                            lambda title, text, kind=0: seen.append(text))
        opened = []
        monkeypatch.setattr(ComparisonDialog, 'exec',
                            lambda self: (opened.append(self), 0)[1])
        from fragalyseqt.main import FragalyseApp
        win = FragalyseApp()
        yield win, seen, opened
        win.close()

    def test_refuses_when_nothing_to_compare(self, window):
        win, seen, opened = window
        assert win.file_states == []
        win.compare_profiles()
        # No files and an empty database: explain, do not silently do nothing.
        assert opened == []
        assert seen and seen[0]

    def test_opens_with_two_saved_profiles_and_no_files(self, window):
        win, seen, opened = window
        db = win._get_refdb()
        _store(db, 'Suspect A', 'Suspect, Known')
        _store(db, 'Victim B', 'Victim, Known')
        win.compare_profiles()
        assert seen == []
        assert len(opened) == 1
        dlg = opened[0]
        first = _selected(dlg._prof1_combo)
        second = _selected(dlg._prof2_combo)
        assert first is not None and second is not None
        assert first != second


class TestSearchWithoutOpenFiles:
    # search_profile_in_db() used to `return` silently unless a file was open,
    # so material from previously analysed cases could never be looked up.

    @pytest.fixture
    def window(self, qapp, tmp_path, monkeypatch):
        import fragalyseqt.fragalyseqt as fq
        import fragalyseqt.boxes as boxes
        monkeypatch.setattr(fq, '_USER_DATA', str(tmp_path))
        monkeypatch.setattr(fq, '_CASEWORK_DB_PATH',
                            str(tmp_path / 'casework.db'))
        monkeypatch.setattr(fq, '_REFPROFILE_DB_PATH',
                            str(tmp_path / 'refprofiles.db'))
        seen = []
        # The handler imports msgbox from .boxes, so patch it at the source.
        monkeypatch.setattr(boxes, 'msgbox',
                            lambda title, text, kind=0: seen.append(text))
        from fragalyseqt.main import FragalyseApp
        win = FragalyseApp()
        yield win, seen, monkeypatch
        win.close()

    def _accept_picking(self, monkeypatch, index=0):
        import fragalyseqt.fragalyseqt as fq
        from pyqtgraph.Qt.QtWidgets import QComboBox

        def _exec(dlg):
            combo = dlg.findChild(QComboBox)
            if combo is not None:
                combo.setCurrentIndex(index)
            return 1
        monkeypatch.setattr(fq.QDialog, 'exec', _exec)

    def test_says_so_when_there_is_nothing_to_search_with(self, window):
        win, seen, _ = window
        assert win.file_states == []
        assert win._pick_query_profile(win._profile_sources()) is None
        assert seen, 'no files and empty databases must not fail silently'

    def test_query_profile_can_come_from_a_database(self, window):
        win, _seen, monkeypatch = window
        _store(win._get_db(), 'From an old case', 'Evidence')
        self._accept_picking(monkeypatch)
        picked = win._pick_query_profile(win._profile_sources())
        assert picked is not None
        calls, origin = picked
        assert len(calls) == len(_CALLS)
        assert origin == ('casework', 1)

    def test_cancelling_the_picker_returns_none(self, window):
        import fragalyseqt.fragalyseqt as fq
        win, _seen, monkeypatch = window
        _store(win._get_db(), 'From an old case', 'Evidence')
        monkeypatch.setattr(fq.QDialog, 'exec', lambda self: 0)
        assert win._pick_query_profile(win._profile_sources()) is None

    def test_search_finds_material_across_both_databases(self, window):
        win, _seen, monkeypatch = window
        _store(win._get_db(), 'Trace from case 12/2024', 'Evidence')
        _store(win._get_refdb(), 'Known suspect', 'Suspect, Known')
        from fragalyseqt.comparison import search_profiles_multi
        results = search_profiles_multi(win._profile_sources(), list(_CALLS))
        found = {(r['source'], r['name']) for r in results}
        assert ('casework', 'Trace from case 12/2024') in found
        assert ('refdb', 'Known suspect') in found

    def test_query_profile_is_not_reported_as_its_own_match(self, window):
        win, _seen, monkeypatch = window
        _store(win._get_db(), 'Trace from case 12/2024', 'Evidence')
        _store(win._get_refdb(), 'Known suspect', 'Suspect, Known')
        self._accept_picking(monkeypatch)
        calls, origin = win._pick_query_profile(win._profile_sources())
        from fragalyseqt.comparison import search_profiles_multi
        results = search_profiles_multi(win._profile_sources(), calls)
        kept = [r for r in results if (r['source'], r['id']) != origin]
        assert origin in {(r['source'], r['id']) for r in results}
        assert origin not in {(r['source'], r['id']) for r in kept}
        assert kept, 'the other database still has a genuine match'
