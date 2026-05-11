import threading
import base64
import csv
import io
import tempfile
import os
from os.path import splitext
from pyqtgraph.Qt import QtCore
from .setvar import set_graph_name


class SOAPBridge(QtCore.QObject):
    """Thread-safe facade exposing FileState data to the SOAP server thread.

    Read operations access file_states directly under an RLock (safe after
    reanalyse() completes under Python's GIL).  Write operations are dispatched
    to the Qt main thread via a queued signal and block until completion.
    """

    _dispatch = QtCore.Signal(object)

    def __init__(self, main_window, timeout=30):
        super().__init__()
        self._w = main_window
        self._timeout = timeout
        self._lock = threading.RLock()
        self._dispatch.connect(self._execute)

    def _execute(self, func):
        func()

    def set_timeout(self, seconds):
        self._timeout = seconds

    def _run_in_main_thread(self, func):
        result = [None]
        error = [None]
        done = threading.Event()

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
                'window_width': s.getwinwidth.value(),
                'baseline_correction': s.bcd.isChecked(),
                'size_standard': s.ILS.currentText(),
                'sizing_method': s.SM.currentText(),
                'panel': s.panel_combo.currentText(),
            }

    def export_csv(self, session_id):
        with self._lock:
            s = self._w.file_states[session_id]
            header, rows = self._w._build_csv_data(s)
        buf = io.StringIO()
        w = csv.writer(buf)
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
        content = base64.b64decode(content_b64)
        suffix = splitext(file_name)[1] or '.fsa'

        def _do():
            with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as f:
                f.write(content)
                tmp_path = f.name
            try:
                return self._w._load_file(tmp_path, tab_name=file_name)
            finally:
                try:
                    os.unlink(tmp_path)
                except OSError:
                    pass

        return self._run_in_main_thread(_do)

    def set_analysis_params(self, session_id, params):
        def _do():
            s = self._w.file_states[session_id]
            widgets = [
                (s.getheight,     'min_height',          int),
                (s.getwidth,      'min_width',           int),
                (s.getprominence, 'min_prominence',      int),
                (s.getwinwidth,   'window_width',        int),
            ]
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
        """Return names of all panels currently in the library."""
        from .panelparser import load_panel_library
        from .fragalyseqt import _PANEL_LIBRARY
        with self._lock:
            return list(load_panel_library(_PANEL_LIBRARY).keys())

    def import_panel(self, panels_name: str, panels_b64: str,
                     bins_name: str = '', bins_b64: str = '',
                     stutter_name: str = '', stutter_b64: str = '') -> list:
        """Import a panel into the library. Returns list of imported panel names."""
        import base64, tempfile, os
        from os.path import splitext

        def _do():
            tmpfiles = []
            try:
                def _write_tmp(name, b64):
                    suffix = splitext(name)[1] or '.txt'
                    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as f:
                        f.write(base64.b64decode(b64))
                        tmpfiles.append(f.name)
                        return f.name

                panels_path  = _write_tmp(panels_name, panels_b64)
                bins_path    = _write_tmp(bins_name, bins_b64)    if bins_b64    else ''
                stutter_path = _write_tmp(stutter_name, stutter_b64) if stutter_b64 else ''

                return self._w._do_import_panel(panels_path, bins_path, stutter_path)
            finally:
                for f in tmpfiles:
                    try:
                        os.unlink(f)
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
