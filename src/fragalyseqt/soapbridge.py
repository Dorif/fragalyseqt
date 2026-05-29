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

from threading import RLock, Event
from base64 import b64decode
from csv import writer as csv_writer
from io import StringIO
from tempfile import NamedTemporaryFile
from os import unlink
from os.path import splitext
from pyqtgraph.Qt.QtCore import QObject, Signal
from .setvar import set_graph_name


class SOAPBridge(QObject):
    # Thread-safe facade exposing FileState data to the SOAP server thread.

    # Read operations access file_states directly under an RLock (safe after
    # reanalyse() completes under Python's GIL).  Write operations are dispatched
    # to the Qt main thread via a queued signal and block until completion.

    _dispatch = Signal(object)

    def __init__(self, main_window, timeout=30):
        super().__init__()
        self._w = main_window
        self._timeout = timeout
        self._lock = RLock()
        self._dispatch.connect(self._execute)

    def _execute(self, func):
        func()

    def set_timeout(self, seconds):
        self._timeout = seconds

    def _run_in_main_thread(self, func):
        result = [None]
        error = [None]
        done = Event()

        def wrapped():
            try:
                result[0] = func()
            except Exception as exc:
                error[0] = exc
            finally:
                done.set()

        self._dispatch.emit(wrapped)
        if not done.wait(self._timeout):
            raise TimeoutError('Qt main thread did not respond within timeout')
        if error[0]:
            raise error[0]
        return result[0]

    # ------------------------------------------------------------------
    # Read operations — called directly from SOAP thread
    # ------------------------------------------------------------------

    def get_session_list(self):
        with self._lock:
            result = []
            for i, s in enumerate(self._w.file_states):
                from .localize import localizefq
                _msg = {}
                localizefq(_msg)
                no_panel = _msg.get('nopanel', 'No panel')
                panel_name = s.panel_combo.currentText()
                result.append({
                    'session_id': i,
                    'file_name': self._w.file_tab.tabText(i),
                    'instrument': set_graph_name(s.abif_raw),
                    'dye_count': len(s.dyerange),
                    'has_sizing': len(s.peaksizes) > 0,
                    'has_panel': bool(s.panel_data and panel_name != no_panel),
                    'panel_name': panel_name,
                })
            return result

    def get_peak_table(self, session_id):
        with self._lock:
            s = self._w.file_states[session_id]
            peaks = []
            for i in range(len(s.peakchannels)):
                ch = int(s.peakchannels[i])
                dye = s.Dye[ch - 1] if 0 < ch <= len(s.Dye) else str(ch)
                size_bp = float(s.peaksizes[i]) if len(s.peaksizes) > 0 else float('nan')
                allele = s.peakalleles[i] if s.peakalleles else ''
                peaks.append({
                    'channel': ch,
                    'dye_name': dye,
                    'position_dp': float(s.peakpositions[i]),
                    'size_bp': size_bp,
                    'height': float(s.peakheights[i]),
                    'fwhm': float(s.peakfwhms[i]),
                    'area': float(s.peakareas[i]),
                    'allele': allele,
                })
            return peaks

    def get_analysis_params(self, session_id):
        with self._lock:
            s = self._w.file_states[session_id]
            return {
                'min_height': s.getheight.value(),
                'min_width': s.getwidth.value(),
                'min_prominence': s.getprominence.value(),
                'halfwindow_width': s.gethalfwin.value(),
                'peak_window': s.getpkwlen.value(),
                'baseline_correction': s.bcd.isChecked(),
                'size_standard': s.ILS.currentText(),
                'sizing_method': s.SM.currentText(),
                'panel': s.panel_combo.currentText(),
                'allele_min_height': s.allele_min_height.value(),
            }

    def export_csv(self, session_id):
        with self._lock:
            s = self._w.file_states[session_id]
            header, rows = self._w._build_csv_data(s)
        buf = StringIO()
        w = csv_writer(buf)
        w.writerow(header)
        w.writerows(rows)
        return buf.getvalue()

    def get_raw_signal(self, session_id, channel):
        with self._lock:
            s = self._w.file_states[session_id]
            if not s.ch or channel < 1 or channel > len(s.ch):
                return []
            return list(s.ch[channel - 1])

    # ------------------------------------------------------------------
    # Write operations — dispatched to Qt main thread
    # ------------------------------------------------------------------

    def submit_file(self, file_name, content_b64):
        content = b64decode(content_b64)
        suffix = splitext(file_name)[1] or '.fsa'

        def _do():
            with NamedTemporaryFile(suffix=suffix, delete=False) as f:
                f.write(content)
                tmp_path = f.name
            try:
                return self._w._load_file(tmp_path, tab_name=file_name)
            finally:
                try:
                    unlink(tmp_path)
                except OSError:
                    pass

        return self._run_in_main_thread(_do)

    def set_analysis_params(self, session_id, params):
        def _do():
            s = self._w.file_states[session_id]
            widgets = [(s.getheight, 'min_height', int),
                       (s.getwidth, 'min_width', int),
                       (s.getprominence, 'min_prominence', int),
                       (s.gethalfwin, 'halfwindow_width', int),
                       (s.getpkwlen, 'peak_window', int),
                       (s.allele_min_height, 'allele_min_height', int),]
            for widget, key, cast in widgets:
                if key in params:
                    widget.blockSignals(True)
                    widget.setValue(cast(params[key]))
                    widget.blockSignals(False)
            if 'baseline_correction' in params:
                checked = params['baseline_correction'].lower() == 'true'
                s.bcd.blockSignals(True)
                s.bcd.setChecked(checked)
                s.bcd.blockSignals(False)
                s.do_BCD = checked
            for widget, key in ((s.ILS, 'size_standard'),
                                (s.SM, 'sizing_method'),
                                (s.panel_combo, 'panel')):
                if key in params:
                    widget.blockSignals(True)
                    widget.setCurrentText(params[key])
                    widget.blockSignals(False)
            self._w.reanalyse(s)
            return True

        return self._run_in_main_thread(_do)

    def trigger_reanalysis(self, session_id):
        def _do():
            s = self._w.file_states[session_id]
            self._w.reanalyse(s)
            return True
        return self._run_in_main_thread(_do)

    def trigger_batch_process(self, session_id):
        def _do():
            orig = self._w.file_tab.currentIndex()
            self._w.file_tab.setCurrentIndex(session_id)
            self._w.process_whole_batch()
            self._w.file_tab.setCurrentIndex(orig)
            return True
        return self._run_in_main_thread(_do)

    def close_session(self, session_id):
        def _do():
            self._w.file_tab.setCurrentIndex(session_id)
            self._w._close_tab_action()
            return True
        return self._run_in_main_thread(_do)

    # ------------------------------------------------------------------
    # Panel operations
    # ------------------------------------------------------------------

    def list_panels(self) -> list:
        # Return names of all panels currently in the library.
        from .panelparser import load_panel_library
        from .fragalyseqt import _PANEL_LIBRARY
        with self._lock:
            return list(load_panel_library(_PANEL_LIBRARY).keys())

    def import_panel(self, panels_name: str, panels_b64: str,
                     bins_name: str = '', bins_b64: str = '',
                     stutter_name: str = '', stutter_b64: str = '') -> list:
        # Import a panel into the library. Returns list of imported panel names.
        def _do():
            tmpfiles = []
            try:
                def _write_tmp(name, b64):
                    suffix = splitext(name)[1] or '.txt'
                    with NamedTemporaryFile(suffix=suffix, delete=False) as f:
                        f.write(b64decode(b64))
                        tmpfiles.append(f.name)
                        return f.name
                panels_path = _write_tmp(panels_name, panels_b64)
                bins_path = _write_tmp(bins_name, bins_b64) if bins_b64 else ''
                stutter_path = _write_tmp(stutter_name, stutter_b64) if stutter_b64 else ''
                return self._w._do_import_panel(panels_path, bins_path, stutter_path)
            finally:
                for f in tmpfiles:
                    try:
                        unlink(f)
                    except OSError:
                        pass

        return self._run_in_main_thread(_do)

    # ------------------------------------------------------------------
    # Database operations
    # ------------------------------------------------------------------

    def save_session(self, name: str) -> int:
        def _do():
            return self._w._do_save_session(name)
        return self._run_in_main_thread(_do)

    def list_sessions(self) -> list:
        with self._lock:
            return self._w._get_db().get_session_list()

    def open_session(self, session_id: int) -> dict:
        from .database import verify_session
        db = self._w._get_db()
        statuses = verify_session(db, session_id)
        readonly = not all(st.status == 'ok' for st in statuses)

        def _do():
            return self._w._do_open_session(session_id, readonly=readonly)
        return self._run_in_main_thread(_do)

    def delete_session(self, session_id: int) -> bool:
        def _do():
            self._w._get_db().hide_session(session_id)
            return True
        return self._run_in_main_thread(_do)

    # ------------------------------------------------------------------
    # Frequency table and reference profile operations
    # ------------------------------------------------------------------

    def _load_freq_table_by_name(self, name: str):
        from os import listdir
        from os.path import isdir, join
        from .fragalyseqt import _FREQTABLES_DIR
        from .freqdb import load_freq_table
        if isdir(_FREQTABLES_DIR):
            for fname in listdir(_FREQTABLES_DIR):
                if fname.endswith('.json'):
                    t = load_freq_table(join(_FREQTABLES_DIR, fname))
                    if t.name == name:
                        return t
        raise ValueError(f'Frequency table not found: {name!r}')

    def _get_calls_soap(self, session_id=None, profile_id=None):
        from .comparison import allele_calls_from_state
        from .refprofile import get_profile
        if session_id is not None:
            return allele_calls_from_state(self._w.file_states[int(session_id)])
        if profile_id is not None:
            return get_profile(self._w._get_refdb(), int(profile_id)).calls
        raise ValueError('Provide session_id or profile_id')

    def list_freq_tables(self) -> list:
        from os import listdir
        from os.path import isdir, join
        from .fragalyseqt import _FREQTABLES_DIR
        from .freqdb import load_freq_table
        names = []
        if not isdir(_FREQTABLES_DIR):
            return names
        for fname in sorted(listdir(_FREQTABLES_DIR)):
            if fname.endswith('.json'):
                try:
                    t = load_freq_table(join(_FREQTABLES_DIR, fname))
                    names.append(t.name)
                except Exception:
                    pass
        return names

    def list_ref_profiles(self) -> list:
        with self._lock:
            db = self._w._get_refdb()
            result = []
            for p in db.list_reference_profiles():
                alleles = db.get_reference_alleles(p['id'])
                result.append({
                    'id': p['id'],
                    'name': p['name'],
                    'role': p['role'] or '',
                    'n_loci': len(alleles),
                })
            return result

    def import_freq_table(self, file_name: str, content_b64: str,
                          table_name: str) -> str:
        from os.path import splitext, join
        from .freqdb import import_freq_csv, import_freq_fam, save_freq_table
        from .fragalyseqt import _FREQTABLES_DIR

        content = b64decode(content_b64)
        suffix = splitext(file_name)[1] or '.csv'
        with NamedTemporaryFile(suffix=suffix, delete=False) as f:
            f.write(content)
            tmp = f.name
        try:
            if suffix.lower() == '.fam':
                t = import_freq_fam(tmp, table_name.strip())
            else:
                t = import_freq_csv(tmp, table_name.strip(), '', '')
        finally:
            try:
                unlink(tmp)
            except OSError:
                pass

        dest = join(_FREQTABLES_DIR,
                    table_name.strip().replace(' ', '_') + '.json')
        save_freq_table(t, dest)
        return t.name

    _VALID_DATABASES = ('casework', 'refprofiles')

    def _resolve_db(self, database: str):
        # Resolve the mandatory "database" parameter to a backend instance.
        # Raises ValueError when the value is missing or not recognised so
        # callers receive a clear error rather than silently writing to the
        # wrong store — important in forensic / medical contexts.
        if not database:
            raise ValueError(
                'database is required; must be one of: '
                + ', '.join(f'"{v}"' for v in self._VALID_DATABASES))
        db_key = database.strip().lower()
        if db_key == 'casework':
            return self._w._get_db()
        if db_key == 'refprofiles':
            return self._w._get_refdb()
        raise ValueError(
            f'Unknown database {database!r}; '
            'must be "casework" or "refprofiles".')

    def import_ref_profile(self, xml_b64: str, database: str,
                           role: str = '') -> list:
        db = self._resolve_db(database)
        content = b64decode(xml_b64)
        with NamedTemporaryFile(suffix='.xml', delete=False) as f:
            f.write(content)
            tmp = f.name
        try:
            from .refprofile import profiles_from_codis_xml, store_profile
            profiles = profiles_from_codis_xml(tmp)
        finally:
            try:
                unlink(tmp)
            except OSError:
                pass

        def _do():
            ids = []
            for p in profiles:
                if role:
                    p.role = role
                ids.append(store_profile(db, p))
            return ids

        return self._run_in_main_thread(_do)

    def _result_to_pdf_b64(self, result) -> str:
        from base64 import b64encode
        from tempfile import NamedTemporaryFile
        from .pdfreport import export_comparison_pdf

        def _do():
            with NamedTemporaryFile(suffix='.pdf', delete=False) as f:
                tmp = f.name
            try:
                export_comparison_pdf(result, tmp)
                with open(tmp, 'rb') as f:
                    return b64encode(f.read()).decode('ascii')
            finally:
                try:
                    unlink(tmp)
                except OSError:
                    pass

        return self._run_in_main_thread(_do)

    def compare_identity(self, params: dict):
        from .comparison import compare_identity
        with self._lock:
            calls_q = self._get_calls_soap(
                params.get('session_id_q'), params.get('profile_id_q'))
            calls_r = self._get_calls_soap(
                params.get('session_id_r'), params.get('profile_id_r'))
            table = self._load_freq_table_by_name(params['table_name'])
            theta = float(params.get('theta', 0.01))
            return compare_identity(calls_q, calls_r, table, theta)

    def export_identity_pdf(self, params: dict) -> str:
        return self._result_to_pdf_b64(self.compare_identity(params))

    def export_kinship_pdf(self, params: dict) -> str:
        return self._result_to_pdf_b64(self.compare_kinship(params))

    def compare_kinship(self, params: dict):
        from .comparison import compare_kinship
        from .forensicstats import RELATIONSHIPS
        with self._lock:
            calls1 = self._get_calls_soap(
                params.get('session_id_q'), params.get('profile_id_q'))
            calls2 = self._get_calls_soap(
                params.get('session_id_r'), params.get('profile_id_r'))
            table = self._load_freq_table_by_name(params['table_name'])
            theta = float(params.get('theta', 0.01))
            rel_name = params.get('relationship', 'Unrelated')
            if rel_name not in RELATIONSHIPS:
                raise ValueError(f'Unknown relationship: {rel_name!r}')
            return compare_kinship(calls1, calls2, table, RELATIONSHIPS[rel_name], theta)

    def store_ref_profile(self, params: dict) -> int:
        from .refprofile import ReferenceProfile, store_profile
        db = self._resolve_db(params.get('database', ''))
        calls = self._get_calls_soap(
            params.get('session_id'), params.get('profile_id_q'))
        name = params.get('name', '').strip()
        if not name:
            raise ValueError('name is required')
        role = params.get('role') or None
        notes = params.get('notes') or None
        profile = ReferenceProfile(name=name, role=role, notes=notes, calls=calls)

        def _do():
            return store_profile(db, profile)

        return self._run_in_main_thread(_do)

    def search_profiles_soap(self, params: dict) -> list:
        from .comparison import search_profiles
        with self._lock:
            calls = self._get_calls_soap(
                params.get('session_id'), params.get('profile_id_q'))
            return search_profiles(self._w._get_refdb(), calls)

    def get_ref_profile(self, profile_id: int) -> dict:
        from .refprofile import get_profile
        with self._lock:
            p = get_profile(self._w._get_refdb(), profile_id)
            return {
                'id': p.id,
                'name': p.name,
                'role': p.role or '',
                'notes': p.notes or '',
                'created_at': p.created_at or '',
                'calls': [
                    {'marker': c.marker,
                     'allele1': c.allele1,
                     'allele2': c.allele2 or ''}
                    for c in p.calls
                ],
            }

    def update_ref_profile(self, profile_id: int, params: dict) -> int:
        from .refprofile import update_profile
        kwargs = {k: v for k, v in params.items()
                  if k in ('name', 'role', 'notes') and v is not None}

        def _do():
            return update_profile(self._w._get_refdb(), profile_id, **kwargs)

        return self._run_in_main_thread(_do)

    def delete_ref_profile(self, profile_id: int) -> bool:
        from .refprofile import delete_profile

        def _do():
            delete_profile(self._w._get_refdb(), profile_id)
            return True

        return self._run_in_main_thread(_do)
