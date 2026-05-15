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

from .boxes import msgbox
from .localize import localizefq
from .soapsettings import load_soap_settings, save_soap_settings, SOAPSettingsDialog
from .database import (DatabaseBackend, SQLiteBackend, open_backend,
                       compute_hashes, verify_session, compress_signal,
                       decompress_signal, InstrumentFileRecord, PeakCallRecord,
                       DyeChannelRecord, AlleleCallRecord, ChannelSignalRecord,
                       AnalysisRunRecord, SavedSessionRecord, SessionTabRecord)
from .session_dialog import SaveSessionDialog, OpenSessionDialog, VerificationDialog
from .codisexport import CODISExportDialog
from .stutterfilter import apply_stutter_filter
from os import makedirs
from os.path import expanduser, dirname, basename, join, isfile, splitext
from shutil import copy2
from csv import writer as csvwriter
from concurrent.futures import ThreadPoolExecutor
from Bio.SeqIO import read as fsaread
from numpy import around, multiply, array, concatenate, transpose, where, abs as npabs, any as npany, isnan as npisnan
from scipy.signal import find_peaks
from scipy.interpolate import splrep, splev
from numpy.polynomial.polynomial import Polynomial
from .jbcd import jbcd
from .ladderalign import align_ils_peaks
# Using FileDialog and SpinBox from pyqtgraph to prevent some possible problems
# for macOS users and to allow more fine variable setting.
from pyqtgraph import PlotWidget, FileDialog, SpinBox, ComboBox, TableWidget
# Using pyqtgraph widgets to make program independent of Qt for Python binding.
from pyqtgraph.Qt.QtWidgets import QCheckBox
from . import fillhid
from xml.etree.ElementTree import parse as xmlparse, Element, SubElement
from .sizestdeditor import SizeStandardEditor
from .setvar import (set_dye_array, set_graph_name, set_spl_dgr, set_lsq_ord,
                     set_knots, chk_key_valid, southern_fit_local,
                     southern_fit_global)
from .panelparser import (parse_genemapper, parse_genemarker, parse_osiris,
                          assign_alleles, _xml_root_tag, load_panel_library,
                          save_panel_library)
from platformdirs import user_data_dir
ftype = "ABI fragment analysis files (*.fsa *.hid);;"
ftype += "Native Nanophore files (*.frf)"
ifacemsg = {}
localizefq(ifacemsg)
homedir = expanduser('~')
_PEN_COLORS = ('b', 'g', 'y', 'r', 'orange', 'c', 'm', 'k')
_USER_DATA = user_data_dir('fragalyseqt', appauthor=False)
_PANEL_LIBRARY = join(_USER_DATA, 'panels.json')
_SIZESTANDARDS = join(_USER_DATA, 'sizestandards.xml')
_SIZESTANDARDS_DEFAULT = join(dirname(__file__), 'sizestandards.xml')
_DB_PATH = join(_USER_DATA, 'sessions.db')
if not isfile(_SIZESTANDARDS):
    makedirs(_USER_DATA, exist_ok=True)
    copy2(_SIZESTANDARDS_DEFAULT, _SIZESTANDARDS)
size_standards = {
    e.get('name'): {
        'channel': e.get('channel'),
        'sizes': [int(x) for x in e.find('Sizes').text.split()]
    }
    for e in xmlparse(_SIZESTANDARDS).getroot()
}


def _refine_peak_positions(signal, positions):
    signal = array(signal, dtype=float)
    positions = array(positions)
    if len(positions) == 0:
        return positions.astype(float)
    refined = positions.astype(float)
    mask = (positions > 0) & (positions < len(signal) - 1)
    pos = positions[mask]
    y_left = signal[pos - 1]
    y_center = signal[pos]
    y_right = signal[pos + 1]
    denom = y_left - 2.0 * y_center + y_right
    valid = denom != 0.0
    delta = where(valid, 0.5*(y_left - y_right)/where(valid, denom, 1.0), 0.0)
    refined[mask] = pos + delta
    return refined


class FileState:
    # Holds all per-file/per-tab analysis state and widget references.
    def __init__(self):
        # Analysis data
        self.abif_raw = None
        self.dyerange = range(0)
        self.Dye = []
        self.udatac = []
        self.ch = []
        self.x_plot = []
        self.show_channels = [1] * 8
        self.do_BCD = False
        self.should_sizecall = False
        self.winwidth = 51
        self.issouthern = False
        self.lsq_order = 0
        self.farr = []
        self.size_std = []
        self.peakpositions = array([])
        self.peakheights = array([])
        self.peakfwhms = array([])
        self.peakchannels = array([])
        self.peaksizes = array([])
        self.peakareas = array([])
        # allele labels from panel binning
        self.peakalleles = []
        # loaded panel for this tab
        self.panel_data = {}
        # Session / database state
        # absolute path of source file
        self.file_path = None
        # True when loaded from DB without source file
        self.readonly = False
        # analysis_run.id in DB (set after save/restore)
        self.db_run_id = None
        # Per-tab widget references (set by _create_tab_content)
        self.plot_widget = None
        self.table_widget = None
        self.getheight = None
        self.getwidth = None
        self.getprominence = None
        self.getwinwidth = None
        self.ILS = None
        self.SM = None
        self.sizecall = None
        self.bcd = None
        self.hidech = []
        self.panel_combo = None
        self.batch_btn = None


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        from pyqtgraph.Qt.QtWidgets import (
            QWidget, QVBoxLayout, QTabWidget,)
        from pyqtgraph.Qt.QtGui import QIcon
        try:
            from pyqtgraph.Qt.QtGui import QAction
        except ImportError:
            from pyqtgraph.Qt.QtWidgets import QAction
        MainWindow.setWindowTitle("FragalyseQt")
        MainWindow.setWindowIcon(QIcon("FragalyseQt.png"))
        MainWindow.resize(1024, 768)
        self.centralwidget = QWidget(MainWindow)
        MainWindow.setCentralWidget(self.centralwidget)

        self.file_states = []
        self._soap_server = None
        self._soap_bridge = None
        self._db = None

        menubar = MainWindow.menuBar()

        file_menu = menubar.addMenu(ifacemsg['menu_file'])
        act_open = QAction(ifacemsg['openfragmentfile'], MainWindow)
        act_open.setShortcut("Ctrl+O")
        act_open.triggered.connect(self.open_and_plot)
        file_menu.addAction(act_open)
        act_close = QAction(ifacemsg['closetab'], MainWindow)
        act_close.setShortcut("Ctrl+W")
        act_close.triggered.connect(self._close_tab_action)
        file_menu.addAction(act_close)
        file_menu.addSeparator()
        act_save_ses = QAction(ifacemsg['savesession'], MainWindow)
        act_save_ses.setShortcut('Ctrl+Shift+S')
        act_save_ses.triggered.connect(self._save_session)
        file_menu.addAction(act_save_ses)
        act_open_ses = QAction(ifacemsg['opensession'], MainWindow)
        act_open_ses.setShortcut('Ctrl+Shift+O')
        act_open_ses.triggered.connect(self._open_session)
        file_menu.addAction(act_open_ses)
        file_menu.addSeparator()
        act_csv = QAction(ifacemsg['csvexport'], MainWindow)
        act_csv.setShortcut("Ctrl+E")
        act_csv.triggered.connect(self.export_csv)
        file_menu.addAction(act_csv)
        act_session = QAction(ifacemsg['exportwholesession'], MainWindow)
        act_session.setShortcut("Ctrl+Shift+E")
        act_session.triggered.connect(self.export_session)
        file_menu.addAction(act_session)
        act_internal = QAction(ifacemsg['exportinternal'], MainWindow)
        act_internal.setShortcut("Ctrl+I")
        act_internal.triggered.connect(self.export_internal)
        file_menu.addAction(act_internal)
        act_codis = QAction(ifacemsg['codisexport'], MainWindow)
        act_codis.setShortcut("Ctrl+Shift+C")
        act_codis.triggered.connect(self.export_codis)
        file_menu.addAction(act_codis)

        settings_menu = menubar.addMenu(ifacemsg['menu_settings'])
        act_soap = QAction(ifacemsg['soapapimenu'], MainWindow)
        act_soap.triggered.connect(self._soap_settings)
        settings_menu.addAction(act_soap)
        settings_menu.addSeparator()
        act_panel = QAction(ifacemsg['importpanel'], MainWindow)
        act_panel.setShortcut("Ctrl+Shift+P")
        act_panel.triggered.connect(self.import_panel_to_library)
        settings_menu.addAction(act_panel)
        act_size_std = QAction(ifacemsg['addsizestd'], MainWindow)
        act_size_std.setShortcut("Ctrl+Shift+A")
        act_size_std.triggered.connect(self.add_size_standard)
        settings_menu.addAction(act_size_std)

        help_menu = menubar.addMenu(ifacemsg['menu_help'])
        act_about = QAction(ifacemsg['aboutbtn'], MainWindow)
        act_about.setShortcut("F1")
        act_about.triggered.connect(self.about)
        help_menu.addAction(act_about)

        root_layout = QVBoxLayout(self.centralwidget)
        root_layout.setContentsMargins(8, 8, 8, 8)
        self.file_tab = QTabWidget(self.centralwidget)
        root_layout.addWidget(self.file_tab)

        self._apply_soap_settings(load_soap_settings())

    @property
    def _state(self):
        idx = self.file_tab.currentIndex()
        if 0 <= idx < len(self.file_states):
            return self.file_states[idx]
        return None

    def _create_tab_content(self, state):
    # Creates the per-tab widget (plot, table, controls) and stores widget
    # references in the given FileState.
        from pyqtgraph.Qt.QtWidgets import (QWidget, QPushButton, QVBoxLayout,
                                            QLabel, QHBoxLayout, QGridLayout,
                                            QSizePolicy,)
        tab_widget = QWidget()
        tab_layout = QHBoxLayout(tab_widget)
        tab_layout.setContentsMargins(4, 4, 4, 4)
        tab_layout.setSpacing(6)

        # Left side: plot + table stacked vertically
        left_widget = QWidget()
        left_layout = QVBoxLayout(left_widget)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.setSpacing(4)

        try:
            expanding_policy = QSizePolicy.Expanding
        except AttributeError:
            expanding_policy = QSizePolicy.Policy.Expanding

        plot = PlotWidget()
        plot.setBackground('#cacaca')
        plot.showGrid(x=True, y=True)
        plot.setLabel("left", "Signal intensity, RFU", color='k')
        plot.setSizePolicy(expanding_policy, expanding_policy)
        left_layout.addWidget(plot, stretch=3)

        table = TableWidget(sortable=False)
        table.setSizePolicy(expanding_policy, expanding_policy)
        left_layout.addWidget(table, stretch=2)

        tab_layout.addWidget(left_widget, stretch=3)

        # Right side: controls panel
        controls_widget = QWidget()
        controls_layout = QGridLayout(controls_widget)
        controls_layout.setContentsMargins(0, 0, 0, 0)
        controls_layout.setHorizontalSpacing(6)
        controls_layout.setVerticalSpacing(6)

        getheightlabel = QLabel()
        getheightlabel.setText(ifacemsg["minph"])
        getheightlabel.setStyleSheet(''' font-size: 10pt; ''')
        controls_layout.addWidget(getheightlabel, 0, 0)

        getheight = SpinBox(minStep=1, dec=True)
        getheight.setRange(1, 64000)
        getheight.setValue(175)
        getheight.setMinimumHeight(20)
        getheight.setStyleSheet(''' font-size: 8pt; ''')
        getheight.valueChanged.connect(self.reanalyse)
        controls_layout.addWidget(getheight, 0, 1)

        getwidthlabel = QLabel()
        getwidthlabel.setText(ifacemsg["minpw"])
        getwidthlabel.setStyleSheet(''' font-size: 10pt; ''')
        controls_layout.addWidget(getwidthlabel, 1, 0)

        getwidth = SpinBox(dec=True)
        getwidth.setRange(1, 16000)
        getwidth.setValue(4)
        getwidth.setMinimumHeight(20)
        getwidth.setStyleSheet(''' font-size: 8pt; ''')
        getwidth.valueChanged.connect(self.reanalyse)
        controls_layout.addWidget(getwidth, 1, 1)

        getprominencelabel = QLabel()
        getprominencelabel.setText(ifacemsg["minpp"])
        getprominencelabel.setStyleSheet(''' font-size: 10pt; ''')
        controls_layout.addWidget(getprominencelabel, 2, 0)

        getprominence = SpinBox(minStep=1, dec=True)
        getprominence.setRange(1, 64000)
        getprominence.setValue(175)
        getprominence.setMinimumHeight(20)
        getprominence.setStyleSheet(''' font-size: 8pt; ''')
        getprominence.valueChanged.connect(self.reanalyse)
        controls_layout.addWidget(getprominence, 2, 1)

        getwinwidthlabel = QLabel()
        getwinwidthlabel.setText(ifacemsg["minww"])
        getwinwidthlabel.setStyleSheet(''' font-size: 10pt; ''')
        controls_layout.addWidget(getwinwidthlabel, 3, 0)

        getwinwidth = SpinBox(minStep=1, dec=True)
        getwinwidth.setRange(1, 1000)
        getwinwidth.setValue(51)
        getwinwidth.setMinimumHeight(20)
        getwinwidth.setStyleSheet(''' font-size: 8pt; ''')
        getwinwidth.valueChanged.connect(self.reanalyse)
        controls_layout.addWidget(getwinwidth, 3, 1)

        hidech = []
        i = 0
        while i < 8:
            cb = QCheckBox()
            cb.setText(ifacemsg['ch_inact_msg'])
            cb.setStyleSheet(''' font-size: 10pt; ''')
            cb.toggled.connect(self.hide_ch)
            cb.number = i
            controls_layout.addWidget(cb, 4 + (i // 2), i % 2)
            hidech.append(cb)
            i += 1

        bcd = QCheckBox()
        bcd.setText(ifacemsg["bcd"])
        bcd.toggled.connect(self.setbcd)
        bcd.setStyleSheet(''' font-size: 10pt; ''')
        controls_layout.addWidget(bcd, 8, 0, 1, 2)

        ILS_combo = ComboBox()
        ILS_combo.setItems(list(size_standards.keys()))
        ILS_combo.setStyleSheet(''' font-size: 10pt; ''')
        controls_layout.addWidget(ILS_combo, 9, 0, 1, 2)

        SM_combo = ComboBox()
        SM_combo.setItems(["Local Southern", "Global Southern",
                           "Linear spline sizing", "Cubic spline sizing",
                           "5th degree spline sizing",
                           "LSQ weighted linear spline sizing",
                           "LSQ weighted cubic spline sizing",
                           "LSQ weighted 5th degree spline sizing",
                           "LSQ 2nd order", "LSQ 3rd order", "LSQ 5th order"])
        SM_combo.setStyleSheet(''' font-size: 10pt; ''')
        controls_layout.addWidget(SM_combo, 10, 0)

        sizecall = QPushButton()
        sizecall.setCheckable(True)
        sizecall.setText("SizeCall")
        sizecall.setStyleSheet(''' font-size: 10pt; ''')
        sizecall.clicked.connect(self.reanalyse)
        controls_layout.addWidget(sizecall, 10, 1)

        panel_combo = ComboBox()
        _library = load_panel_library(_PANEL_LIBRARY)
        panel_combo.setItems([ifacemsg["nopanel"]] + list(_library.keys()))
        panel_combo.setEnabled(True)
        state.panel_data = _library
        panel_combo.setStyleSheet(''' font-size: 10pt; ''')
        panel_combo.currentIndexChanged.connect(self.reanalyse)
        controls_layout.addWidget(panel_combo, 11, 0, 1, 2)

        batch_btn = QPushButton(ifacemsg['processbatch'])
        batch_btn.setStyleSheet(''' font-size: 10pt; ''')
        batch_btn.clicked.connect(self.process_whole_batch)
        controls_layout.addWidget(batch_btn, 12, 0, 1, 2)

        state.batch_btn = batch_btn

        controls_layout.setColumnStretch(0, 1)
        controls_layout.setColumnStretch(1, 0)

        tab_layout.addWidget(controls_widget, stretch=1)

        # Store widget references in the state object
        state.plot_widget = plot
        state.table_widget = table
        state.getheight = getheight
        state.getwidth = getwidth
        state.getprominence = getprominence
        state.getwinwidth = getwinwidth
        state.hidech = hidech
        state.bcd = bcd
        state.ILS = ILS_combo
        state.SM = SM_combo
        state.sizecall = sizecall
        state.panel_combo = panel_combo

        return tab_widget

    def _close_tab_action(self):
        idx = self.file_tab.currentIndex()
        if idx < 0:
            return
        self.file_tab.removeTab(idx)
        self.file_states.pop(idx)

    # Checkboxes w/o designations or with designations of nonexistent channels
    # would look weird, so let's inactivate them correctly.
    def inactivatechkboxes(self):
        s = self._state
        if s is None:
            return
        for cb in s.hidech:
            cb.setText(ifacemsg['ch_inact_msg'])

    def _parse_panel_files(self):
    # Open file dialog(s) and parse a panel file. Returns the parsed panel
    # dict, or None if the user cancelled or an error occurred.
        global homedir
        path, _ = FileDialog.getOpenFileName(
            self, ifacemsg['loadpaneldlg'], homedir,
            "Panel files (*.txt *.xml)")
        if not path:
            return None
        try:
            if path.lower().endswith('.xml'):
                if _xml_root_tag(path) == 'KitData':
                    data = parse_osiris(path)
                else:
                    data = parse_genemarker(path)
            else:
                data = parse_genemapper(path)
                bins_path, _ = FileDialog.getOpenFileName(
                    self, ifacemsg['loadbinsdlg'], dirname(path),
                    "Bins files (*.txt)")
                if bins_path:
                    data = parse_genemapper(path, bins_path)
                else:
                    msgbox("", ifacemsg['nobinsmsg'], 0)
                stutter_path, _ = FileDialog.getOpenFileName(
                    self, ifacemsg.get('loadstutterdlg'),
                    dirname(path), "Stutter files (*.txt)")
                if stutter_path:
                    data = parse_genemapper(path, bins_path, stutter_path,)
        except Exception as exc:
            msgbox("", str(exc), 2)
            return None
        if not data:
            msgbox("", ifacemsg['nodatamsg'], 1)
            return None
        return data

    def _do_import_panel(self, panels_path: str, bins_path: str = '',
                         stutter_path: str = '') -> list:
        # Parse panel files and save to library. Returns list of panel names.
        # Raises ValueError on parse failure. Must run on Qt main thread.
        if panels_path.lower().endswith('.xml'):
            if _xml_root_tag(panels_path) == 'KitData':
                data = parse_osiris(panels_path)
            else:
                data = parse_genemarker(panels_path)
        else:
            data = parse_genemapper(panels_path, bins_path, stutter_path,)
        if not data:
            raise ValueError('No panel data found in file')
        save_panel_library(data, _PANEL_LIBRARY)
        library = load_panel_library(_PANEL_LIBRARY)
        for s in self.file_states:
            current = s.panel_combo.currentText()
            s.panel_data = library
            s.panel_combo.setItems([ifacemsg["nopanel"]] + list(library.keys()))
            s.panel_combo.setEnabled(True)
            if current in library:
                s.panel_combo.setValue(current)
        return list(data.keys())

    def import_panel_to_library(self):
        # Parse a panel file and save it permanently to panels.json.
        data = self._parse_panel_files()
        if data is None:
            return
        try:
            # _parse_panel_files already parsed; save and refresh
            save_panel_library(data, _PANEL_LIBRARY)
            library = load_panel_library(_PANEL_LIBRARY)
            for s in self.file_states:
                current = s.panel_combo.currentText()
                s.panel_data = library
                s.panel_combo.setItems([ifacemsg["nopanel"]] + list(library.keys()))
                s.panel_combo.setEnabled(True)
                if current in library:
                    s.panel_combo.setValue(current)
        except Exception as exc:
            msgbox("", str(exc), 2)
        msgbox("", ifacemsg.get('panelimported'), 0)

    def add_size_standard(self):
        # Add a new size standard to the library.
        dlg = SizeStandardEditor(self, ifacemsg)
        if dlg.exec():
            data = dlg.get_data()
            if not data['name'] or not data['sizes']:
                return
            # Update XML
            tree = xmlparse(_SIZESTANDARDS)
            root = tree.getroot()
            ladder = SubElement(root, 'Ladder', {
                'name': data['name'], 'channel': data['channel']
            })
            sizes_elem = SubElement(ladder, 'Sizes')
            sizes_elem.text = data['sizes']
            tree.write(_SIZESTANDARDS, encoding='UTF-8', xml_declaration=True)
            # Refresh global dictionary
            global size_standards
            size_standards[data['name']] = {
                'channel': data['channel'],
                'sizes': [int(x) for x in data['sizes'].split()]
            }
            # Update open tabs
            for s in self.file_states:
                current = s.ladder_combo.currentText()
                s.ladder_combo.setItems(list(size_standards.keys()))
                if current in size_standards:
                    s.ladder_combo.setValue(current)
            msgbox("", ifacemsg.get('stdadded'), 0)

    def _load_file(self, fname, tab_name=None):
        # Load a single file into a new tab. Returns new session_id or -1.
        global homedir
        udatac = fillhid.UDATAC
        label = tab_name or basename(fname)
        if fname.lower().endswith('.frf'):
            from . import fillfrf
            try:
                abif_result = fillfrf.parse_frf(fname)
            except Exception:
                msgbox(ifacemsg['dmgdfile'], ifacemsg['nodatamsg'], 2)
                return -1
            homedir = dirname(fname)
        else:
            FAfile = open(fname, "rb")
            try:
                tmprecord = fsaread(FAfile, "abi")
            except AssertionError:
                class record():
                    annotations = {"abif_raw": {
                        "DATA1": None, "DATA2": None, "DATA3": None,
                        "DATA4": None, "Dye#1": None, "DyeN1": None,
                        "DyeN2": None, "DyeN3": None, "DyeN4": None,
                        "MODL1": None}}
                tmprecord = record()
# Preventing data corruption in a case if target file is corrupted.
            FAfile.close()
# Closing file to save memory and avoid unexpected things.
            tmpabif = tmprecord.annotations["abif_raw"]
            if tmpabif["DATA1"] is None:
                try:
                    fillhid.parse_hid(fname, tmpabif, ifacemsg)
                except Exception:
                    msgbox(ifacemsg['dmgdfile'], ifacemsg['nodatamsg'], 2)
                    return -1
            abif_result = tmpabif
# We need raw ABIF data only, no need in whole data structure, created by
# BioPython's AbiIO. This way multiple brackets constructions are evaded.
            homedir = dirname(fname)
        state = FileState()
        state.abif_raw = abif_result
        state.file_path = fname
        state.udatac = udatac
        state.Dye = set_dye_array(abif_result)
        state.dyerange = range(abif_result["Dye#1"])
        tab_widget = self._create_tab_content(state)
        self.file_states.append(state)
        self.file_tab.addTab(tab_widget, label)
        self.file_tab.setCurrentIndex(len(self.file_states) - 1)
        self.reanalyse()
        return len(self.file_states) - 1

    def open_and_plot(self):
        fnames, _ = FileDialog.getOpenFileNames(self,
                                               'Open files for analysis',
                                               homedir, ftype)
        if not fnames:
            return
        for fname in fnames:
            self._load_file(fname)

    # ------------------------------------------------------------------
    # Database / session helpers
    # ------------------------------------------------------------------

    def _get_db(self) -> DatabaseBackend:
        if self._db is None:
            makedirs(_USER_DATA, exist_ok=True)
            self._db = open_backend({'backend': 'sqlite', 'path': _DB_PATH})
        return self._db

    def _save_session(self):
        if not self.file_states:
            return
        dlg = SaveSessionDialog(parent=self)
        if not dlg.exec():
            return
        name = dlg.get_name()
        if not name:
            return
        self._do_save_session(name)

    def _do_save_session(self, name: str) -> int:
        # Save current open tabs to DB. Returns new session_id.
        import getpass
        created_by = getpass.getuser()
        db = self._get_db()
        existing = db.find_session_by_name(name)
        from .setvar import set_graph_name
        with db.transaction():
            session_id = db.store_session(SavedSessionRecord(
                created_by=created_by, name=name,
                supersedes_id=existing['id'] if existing else None,))
            for tab_order, s in enumerate(self.file_states):
                if s.readonly and s.db_run_id is not None:
                    # Reuse existing run record — no need to re-hash
                    db.store_session_tab(SessionTabRecord(
                        session_id=session_id, tab_order=tab_order,
                        run_id=s.db_run_id))
                    continue
                if not s.file_path:
                    continue
                hashes = compute_hashes(s.file_path)
                file_id = db.store_file(InstrumentFileRecord(
                    created_by=created_by, file_name=basename(s.file_path),
                    file_path=s.file_path, file_size=hashes['size'],
                    hash_md5=hashes['hash_md5'], hash_sha1=hashes['hash_sha1'],
                    hash_sha256=hashes['hash_sha256'],
                    hash_sha3_256=hashes['hash_sha3_256'],
                    instrument=set_graph_name(s.abif_raw),))
                # Dye channels (only if not already stored for this file_id)
                existing_dyes = db.get_dye_channels(file_id)
                if not existing_dyes:
                    for i, dye_name in enumerate(s.Dye):
                        db.store_dye_channel(DyeChannelRecord(
                            file_id=file_id, channel=i + 1, dye_name=dye_name))
                # Channel signals (only if not already stored)
                if not db.get_channel_signals(file_id) and s.ch:
                    for i, signal in enumerate(s.ch):
                        if len(signal):
                            db.store_channel_signal(ChannelSignalRecord(
                                file_id=file_id,
                                channel=i + 1,
                                signal=compress_signal(signal)))
                run_id = db.store_analysis_run(AnalysisRunRecord(
                    file_id=file_id, created_by=created_by,
                    min_height=s.getheight.value(),
                    min_prominence=s.getprominence.value(),
                    min_width=s.getwidth.value(),
                    window_width=s.getwinwidth.value(),
                    baseline_correction=s.do_BCD,
                    sizing_method=s.SM.currentText(),
                    size_standard=s.ILS.currentText(),
                    panel=s.panel_combo.currentText(),
                    supersedes_id=s.db_run_id,
                ))
                s.db_run_id = run_id
                # Peaks
                n = len(s.peakchannels)
                for i in range(n):
                    ch = int(s.peakchannels[i])
                    dye = s.Dye[ch - 1] if 0 < ch <= len(s.Dye) else str(ch)
                    bp = float(s.peaksizes[i]) if len(s.peaksizes) > 0 else None
                    allele_full = s.peakalleles[i] if s.peakalleles else ''
                    is_ladder = allele_full == 'ILS'
                    peak_id = db.store_peak_call(PeakCallRecord(
                        run_id=run_id, created_by=created_by, channel=ch,
                        dye_name=dye, position_dp=float(s.peakpositions[i]),
                        position_bp=bp, height=float(s.peakheights[i]),
                        area=float(s.peakareas[i]) if len(s.peakareas) > 0 else None,
                        fwhm=float(s.peakfwhms[i]) if len(s.peakfwhms) > 0 else None,
                        is_ladder=is_ladder,))
                    # Parse "Marker:Allele" or "Marker:Allele:Stutter" format
                    if allele_full and not is_ladder and allele_full != 'OL':
                        parts = allele_full.split(':', 2)
                        marker_name  = parts[0] if len(parts) > 1 else ''
                        allele_label = parts[1] if len(parts) > 1 else allele_full
                        is_stutter   = len(parts) == 3 and parts[2].lower() == 'stutter'
                        db.store_allele_call(AlleleCallRecord(
                            peak_id=peak_id, created_by=created_by,
                            allele=allele_label, marker=marker_name,
                            is_stutter=is_stutter,))
                db.store_session_tab(SessionTabRecord(
                    session_id=session_id, tab_order=tab_order, run_id=run_id,))
        return session_id

    def _open_session(self):
        db = self._get_db()
        sessions = db.get_session_list()
        if not sessions:
            msgbox('', ifacemsg.get('nosessions'), 1)
            return
        dlg = OpenSessionDialog(sessions, parent=self)
        if not dlg.exec():
            return
        session_id = dlg.get_selected_session_id()
        if session_id is None:
            return
        statuses = verify_session(db, session_id)
        all_ok = all(st.status == 'ok' for st in statuses)
        if not all_ok:
            vdlg = VerificationDialog(statuses, parent=self)
            if not vdlg.exec():
                return
        self._do_open_session(session_id, readonly=not all_ok)

    def _do_open_session(self, session_id: int, readonly: bool=False) -> dict:
        # Restore a session. Returns {'readonly': bool, 'tabs': int}.
        db = self._get_db()
        tabs = db.get_session_tabs(session_id)
        for tab in tabs:
            run = db.get_run_info(tab['run_id'])
            file_info = db.get_file_info(run['file_id'])
            if readonly:
                self._restore_readonly_tab(run, file_info)
            else:
                self._restore_full_tab(run, file_info)
        return {'readonly': readonly, 'tabs': len(tabs)}

    def _restore_full_tab(self, run: dict, file_info: dict):
        # Load file from disk and apply stored analysis parameters.
        sid = self._load_file(file_info['file_path'])
        if sid < 0:
            return
        s = self.file_states[sid]
        s.db_run_id = run['id']
        widgets = [(s.getheight, 'min_height', int),
                   (s.getwidth, 'min_width', int),
                   (s.getprominence, 'min_prominence', int),
                   (s.getwinwidth, 'window_width', int),]
        for widget, key, cast in widgets:
            if run.get(key) is not None:
                widget.blockSignals(True)
                widget.setValue(cast(run[key]))
                widget.blockSignals(False)
        bcd = bool(run.get('baseline_correction', 0))
        s.bcd.blockSignals(True)
        s.bcd.setChecked(bcd)
        s.bcd.blockSignals(False)
        s.do_BCD = bcd
        for widget, key in ((s.ILS, 'size_standard'), (s.SM, 'sizing_method'),
                            (s.panel_combo, 'panel')):
            if run.get(key):
                widget.blockSignals(True)
                widget.setCurrentText(run[key])
                widget.blockSignals(False)
        self.reanalyse(s)

    def _restore_readonly_tab(self, run: dict, file_info: dict):
        # Reconstruct a tab entirely from database data (no source file).
        from numpy import array as nparray, nan, isnan

        db = self._get_db()
        peak_rows = db.get_peak_calls_for_run(run['id'])
        allele_rows = db.get_allele_calls_for_run(run['id'])
        dye_rows = db.get_dye_channels(file_info['id'])
        allele_map = {r['peak_id']: r for r in allele_rows}
        s = FileState()
        s.readonly = True
        s.abif_raw = {}
        s.db_run_id = run['id']
        s.file_path = file_info['file_path']
        s.Dye = [r['dye_name'] for r in dye_rows]
        s.dyerange = range(len(s.Dye))
        s.udatac, channels, areas, alleles = [], [], [], []
        pos_dp, pos_bp, heights, fwhms = [], [], [], []
        for p in peak_rows:
            pos_dp.append(p['position_dp'])
            pos_bp.append(p['position_bp'] if p['position_bp'] is not None
                          else nan)
            heights.append(p['height'])
            fwhms.append(p['fwhm'])
            channels.append(p['channel'])
            areas.append(p['area'])
            if p['is_ladder']:
                alleles.append('ILS')
            elif p['id'] in allele_map:
                ac = allele_map[p['id']]
                marker = ac.get('marker')
                label  = ac.get('allele')
                suffix = ':Stutter' if ac.get('is_stutter') else ''
                alleles.append(f'{marker}:{label}{suffix}' if marker else label)
            else:
                alleles.append('')
        s.peakpositions = nparray(pos_dp)
        s.peakheights = nparray(heights)
        s.peakfwhms = nparray(fwhms)
        s.peakchannels = nparray(channels)
        s.peakareas = nparray(areas)
        s.peakalleles = alleles
        bp_arr = nparray(pos_bp)
        s.peaksizes = bp_arr if not all(isnan(bp_arr)) else nparray([])
        # Raw channel signals → electropherogram in read-only mode
        signal_rows = db.get_channel_signals(file_info['id'])
        if signal_rows:
            s.ch = [decompress_signal(r['signal']) for r in signal_rows]
            s.x_plot = list(range(1, len(s.ch[0]) + 1))
            ladder_pts = sorted(
                (p['position_dp'], p['position_bp'])
                for p in peak_rows
                if p['is_ladder'] and p['position_bp'] is not None)
            if len(ladder_pts) >= 2:
                dp_pts = [pt[0] for pt in ladder_pts]
                bp_pts = [pt[1] for pt in ladder_pts]
                spl = splrep(dp_pts, bp_pts, k=min(3, len(ladder_pts) - 1))
                s.x_plot = list(around(splev(s.x_plot, spl), 3))
                s.should_sizecall = True
        tab_widget = self._create_tab_content(s)
        self._disable_tab_controls(s)
        tab_label = (basename(file_info['file_path']) + ' [RO]')
        self.file_states.append(s)
        self.file_tab.addTab(tab_widget, tab_label)
        self.file_tab.setCurrentIndex(len(self.file_states) - 1)
        # Fill table with stored data (skips findpeaks + allele binning)
        self.retab(s)
        self.replot(s)

    def _disable_tab_controls(self, s: 'FileState'):
        # Grey out all analysis controls for a read-only tab.
        for w in ([s.getheight, s.getwidth, s.getprominence, s.getwinwidth,
                   s.bcd, s.ILS, s.SM, s.sizecall, s.panel_combo, s.batch_btn]
                  + s.hidech):
            if w is not None:
                w.setEnabled(False)

    def _soap_settings(self):
        settings = load_soap_settings()
        dlg = SOAPSettingsDialog(settings, parent=self)
        if dlg.exec():
            new_settings = dlg.get_settings()
            save_soap_settings(new_settings)
            self._apply_soap_settings(new_settings)

    def _apply_soap_settings(self, settings):
        if self._soap_server is not None:
            self._soap_server.stop()
            self._soap_server = None
        if not settings.get('enabled'):
            return
        from .soapbridge import SOAPBridge
        from .soapserver import SOAPServerThread
        if self._soap_bridge is None:
            self._soap_bridge = SOAPBridge(self)
        self._soap_bridge.set_timeout(settings.get('timeout', 30))
        self._soap_server = SOAPServerThread(
            bridge=self._soap_bridge, host=settings.get('host', '127.0.0.1'),
            port=settings.get('port', 8742),
            token=settings.get('token') or None,)
        self._soap_server.start()

    def about(self):
        msgbox(ifacemsg['aboutbtn'], ifacemsg['infoboxtxt'], 0)

    def findpeaks(self, state=None):
        s = state if state is not None else self._state
        if s is None or s.readonly:
            return

        def _sizingerror():
            msgbox("", ifacemsg['wrongsizing'], 1)
            s.sizecall.setChecked(False)
            s.should_sizecall = False
            # Guard: tell reanalyse() not to re-trigger the auto-sizing that
            # just failed, otherwise panel_data being set would cause a loop.
            self._sizing_recovery = True
            self.reanalyse(s)
            self._sizing_recovery = False

        # Detecting peaks and calculating peaks data.
        h = s.getheight.value()
        w = s.getwidth.value()
        p = s.getprominence.value()
        s.winwidth = s.getwinwidth.value()
        _positions = []
        _heights = []
        _fwhms = []
        _channels = []
        _sizes = []
        s.ch = []
        s.farr = []
        s.lsq_order = 0
        s.issouthern = False
        spline = None
        spline_degree = 0
        func = None
        southern_func = None
        _bcd_cache = {}
        s.x_plot = list(range(1, len(s.abif_raw["DATA1"]) + 1))
        if s.should_sizecall:
            ILS_Name = s.ILS.currentText()
            Sizing_Method = s.SM.currentText()
            try:
                ils_channel = size_standards[ILS_Name]['channel']
                ils_data = s.abif_raw[ils_channel]
                if s.do_BCD:
                    _, params = jbcd(ils_data,
                                     half_window=(s.winwidth-1)//2)
                    ils_data = params['signal']
                    _bcd_cache[ils_channel] = ils_data
                s.size_std = size_standards[ILS_Name]['sizes']
                n_expected = len(s.size_std)
                for _attempt in range(4):
                    ILSP = find_peaks(ils_data, height=h, width=w,
                                      prominence=p, wlen=s.winwidth,
                                      rel_height=0.5)
                    if len(ILSP[0]) >= n_expected:
                        break
                    h = max(h // 2, 1)
                    p = max(p // 2, 1)
                # Refine all detected ILS peak positions before alignment.
                all_refined = _refine_peak_positions(ils_data, ILSP[0])
                all_heights = ILSP[1]['peak_heights']

                # Smart alignment: iterative polynomial refinement + DP.
                # Correctly handles the common mismatch where the actual CE run
                # used a longer protocol than the selected size standard
                # (e.g. full LIZ-600 run with a 60–460 definition selected).
                aligned_dp = align_ils_peaks(all_refined, s.size_std,
                                             all_heights)

                if not npisnan(aligned_dp).any():
                    # pattern match succeeded — use directly
                    ladder_peaks   = aligned_dp
                    _size_std_used = array(s.size_std)
                else:
                    # fewer detected peaks than expected — fall back
                    beginning_index = len(ILSP[0]) - n_expected
                    ladder_peaks   = _refine_peak_positions(
                        ils_data, ILSP[0][beginning_index:])
                    _size_std_used = array(s.size_std)

                if 'spline' in Sizing_Method:
                    spline_degree = set_spl_dgr(Sizing_Method)
                    knots = set_knots(Sizing_Method, ladder_peaks,
                                      spline_degree)
                if spline_degree != 0:
                    spline = splrep(ladder_peaks, _size_std_used,
                                    k=spline_degree, t=knots)
                    s.x_plot = around(splev(s.x_plot, spline), 3)
                elif 'order' in Sizing_Method:
                    s.lsq_order = set_lsq_ord(Sizing_Method)
                    func = Polynomial.fit(ladder_peaks, _size_std_used,
                                          s.lsq_order)
                    s.x_plot = around(func(array(s.x_plot)), 3)
                elif 'Southern' in Sizing_Method:
                    s.issouthern = True
                    southern_func = (southern_fit_local
                                     if 'Local' in Sizing_Method
                                     else southern_fit_global)
                    s.x_plot = around(southern_func(
                        ladder_peaks, _size_std_used, s.x_plot), 3)
            except (ValueError, TypeError, KeyError):
                _sizingerror()
                return
    # By default, find_peaks function measures width at half maximum of height.
    # But explicit is always better, then implicit, so it is specified clearly.
        if s.do_BCD:
            half_win = (s.winwidth-1)//2

            def _bcd_channel(chnum):
                key = s.udatac[chnum]
                if key in _bcd_cache:
                    return _bcd_cache[key]
                _, params = jbcd(s.abif_raw[key], half_window=half_win)
                return params['signal']
            with ThreadPoolExecutor() as executor:
                s.ch = list(executor.map(_bcd_channel, s.dyerange))
        else:
            for chnum in s.dyerange:
                s.ch.append(array(s.abif_raw[s.udatac[chnum]], dtype=float))

        def _detect_peaks(chnum):
            return find_peaks(s.ch[chnum], height=h, width=w, prominence=p,
                              wlen=s.winwidth, rel_height=0.5)
        with ThreadPoolExecutor() as executor:
            chP = list(executor.map(_detect_peaks, s.dyerange))
        for chnum in s.dyerange:
            refined_pos = _refine_peak_positions(s.ch[chnum], chP[chnum][0])
            _positions.append(refined_pos)
            _heights.append(chP[chnum][1]['peak_heights'])
            _fwhms.append(chP[chnum][1]['widths'])
            if s.should_sizecall and len(refined_pos) != 0:
                if spline_degree != 0:
                    _sizes.append(splev(refined_pos, spline))
                elif s.issouthern:
                    _sizes.append(southern_func(
                        ladder_peaks, s.size_std, refined_pos))
                else:
                    _sizes.append(func(refined_pos))
            _channels.append([chnum + 1]*len(chP[chnum][0]))  # 1-based index
        s.peakpositions = concatenate(_positions) if _positions else array([])
        s.peakheights = concatenate(_heights) if _heights else array([])
        s.peakfwhms = concatenate(_fwhms) if _fwhms else array([])
        s.peakchannels = concatenate(_channels) if _channels else array([])
        s.peaksizes = concatenate(_sizes) if _sizes else array([])
        if len(s.peaksizes) > 0:
            valid = s.peaksizes >= 0
            s.peakpositions = s.peakpositions[valid]
            s.peakheights = s.peakheights[valid]
            s.peakfwhms = s.peakfwhms[valid]
            s.peakchannels = s.peakchannels[valid]
            s.peaksizes = s.peaksizes[valid]
# Calculate areas from full-precision values, then round everything.
        s.peakareas = around(multiply(s.peakheights, s.peakfwhms)*1.0645, 2)
        s.peaksizes = around(s.peaksizes, 2)
        s.peakheights = around(s.peakheights, 2)
        s.peakfwhms = around(s.peakfwhms, 2)

# Peak areas are calculated using formula for Gaussian peak area
# (https://www.physicsforums.com/threads/area-under-gaussian-peak-by-easy-measurements.419285/):
# A = FWHM*H/(2sqrt(2ln(2))/sqrt(2pi)) = 1.0645*FWHM*H, where FWHM is Full
# Width at Half Maximum. Real area may differ for non-Gaussian peaks, but at
# least majority of them are of Gaussian shape. If peaks are well separated -
# just calculate their area, but if your peaks are crowded (e.g. in TP-PCR or
# allelic ladders), oversaturated or you have noisy data - you MUST use
# baseline correction and denoising prior peak area calculation.

    def replot(self, state=None):
        s = state if state is not None else self._state
        if s is None:
            return
        if s.readonly:
            s.plot_widget.clear()
            s.plot_widget.plotItem.setLimits(xMin=None, xMax=None, yMin=None,
                                             yMax=None)
            for i in s.dyerange:
                s.hidech[i].setText(ifacemsg['hidechannel'] + s.Dye[i])
            ro_title = (basename(s.file_path) if s.file_path else 'Unknown') \
                       + ' [RO]'
            s.plot_widget.setTitle(ro_title, color='r', size='10pt')
            x_label = 'Size, bases' if s.should_sizecall else 'Size, data points'
            if not s.ch or not s.x_plot:
                s.plot_widget.setLabel('bottom', x_label, color='k')
                return
            s.plot_widget.setLabel('bottom', x_label, color='k')
            max_y = 0
            for i in s.dyerange:
                if s.show_channels[i] and i < len(s.ch) and len(s.ch[i]):
                    ch_max = s.ch[i].max()
                    if ch_max > max_y:
                        max_y = ch_max
                    s.plot_widget.plot(s.x_plot, s.ch[i], pen=_PEN_COLORS[i])
            if max_y == 0:
                max_y = 64000
            s.plot_widget.plotItem.setLimits(xMin=0, xMax=s.x_plot[-1], yMin=0,
                                             yMax=max_y)
            return
        s.plot_widget.clear()
        s.plot_widget.plotItem.setLimits(xMin=None, xMax=None, yMin=None,
                                         yMax=None)
        for i in s.dyerange:
            s.hidech[i].setText(ifacemsg['hidechannel'] + s.Dye[i])
        s.plot_widget.setTitle(set_graph_name(s.abif_raw), color="k",
                               size="10pt")
        max_x = len(s.x_plot)
        if s.should_sizecall or len(s.peaksizes) > 0:
            # In the most normal case if you have good overall CE data quality,
            # the last member of x_plot array should have the biggest size.
            x_max = s.x_plot[len(s.x_plot)-1]
            s.plot_widget.setLabel('bottom', 'Size, bases', color='k')
            max_ladder = max(s.size_std)
            if max_ladder+200 < x_max:
                max_x = max_ladder+200
            # EXTREMELY weird situation, but it sometimes happens, e.g. for low
            # CE data quality with a lot of noise or poorly injected ladder.
            elif max_ladder > x_max:
                max_x = max_ladder
            else:
                max_x = x_max
        else:
            s.plot_widget.setLabel('bottom', 'Size, data points', color='k')
        max_y = 0
        for i in s.dyerange:
            if s.show_channels[i]:
                if s.should_sizecall or len(s.peaksizes) > 0:
                    valid = array(s.x_plot) >= 0
                    ch_max = float(s.ch[i][valid].max()) if valid.any() else 0
                else:
                    ch_max = s.ch[i].max()
                if ch_max > max_y:
                    max_y = ch_max
                s.plot_widget.plot(s.x_plot, s.ch[i], pen=_PEN_COLORS[i])
        if max_y == 0:
            max_y = 64000
        s.plot_widget.plotItem.setLimits(xMin=0, xMax=max_x, yMin=0,
                                         yMax=max_y)

    def _build_csv_data(self, s):
        # Return (header, rows) for a FileState, ready to write as CSV.
        header = ['Peak Channel', 'Peak Position (Datapoints)', 'Peak Height',
                  'Peak FWHM', 'Peak Area (Datapoints)']
        ch_names = [s.Dye[int(ch) - 1] if 0 < int(ch) <= len(s.Dye)
                    else str(ch) for ch in s.peakchannels]
        pdarray = [ch_names, s.peakpositions, s.peakheights,
                   s.peakfwhms, s.peakareas]
        if len(s.peaksizes) > 0:
            pdarray.append(s.peaksizes)
            header += ['Peak Size (Bases)']
        if any(s.peakalleles):
            pdarray.append(s.peakalleles)
            header += ['Allele']
        return header, transpose(pdarray)

    def export_csv(self):
        s = self._state
        if s is None or s.abif_raw is None:
            return
        header, peak_data = self._build_csv_data(s)
        csvname, _ = FileDialog.getSaveFileName(self, ifacemsg['savecsv'],
                                                homedir, 'CSV(*.csv)')
        if not csvname:
            return
        with open(csvname, 'w', encoding='UTF8', newline='') as f:
            w = csvwriter(f)
            w.writerow(header)
            w.writerows(peak_data)

    def export_session(self):
        if not self.file_states:
            return
        folder = FileDialog.getExistingDirectory(self,
                                                 ifacemsg['choosefolder'],
                                                 homedir)
        if not folder:
            return
        seen = {}
        for i, s in enumerate(self.file_states):
            stem = splitext(self.file_tab.tabText(i))[0]
            count = seen.get(stem, 0) + 1
            seen[stem] = count
            fname = stem if count == 1 else f'{stem}_{count}'
            csvpath = join(folder, fname + '.csv')
            header, peak_data = self._build_csv_data(s)
            with open(csvpath, 'w', encoding='UTF8', newline='') as f:
                w = csvwriter(f)
                w.writerow(header)
                w.writerows(peak_data)

    def export_internal(self):
        # Exporting internal analysis data for ABI 3500 / SeqStudio files.
        s = self._state
        if s is None or s.abif_raw is None:
            return
        if not chk_key_valid("Peak1", s.abif_raw):
            msgbox(ifacemsg['unsuppeq'], ifacemsg['unsuppeqmsg'], 1)
            return
        header = ['Peak Channel', 'Peak Position (Datapoints)', 'Peak Height',
                  'Peak FWHM', 'Peak Area (Datapoints)',
                  'Peak Size (Bases)', 'Peak Area (Bases)']
        peak_chn = []
        for channel in s.abif_raw["Peak1"]:
            idx = min(max(channel - 1, 0), len(s.Dye) - 1)
            peak_chn.append(s.Dye[idx])
        pdarray = [peak_chn, s.abif_raw["Peak2"], s.abif_raw["Peak7"],
                   s.abif_raw["Peak5"], s.abif_raw["Peak10"],
                   s.abif_raw["Peak12"], s.abif_raw["Peak17"]]
        peak_data = transpose(pdarray)
        csvname, _ = FileDialog.getSaveFileName(self, ifacemsg['savecsv'],
                                                homedir, 'CSV(*.csv)')
        if not csvname:
            return
        with open(csvname, 'w', encoding='UTF8', newline='') as f:
            w = csvwriter(f)
            w.writerow(header)
            w.writerows(peak_data)

    def export_codis(self):
        if not self.file_states:
            return
        tab_names = [self.file_tab.tabText(i)
                     for i in range(self.file_tab.count())]
        dlg = CODISExportDialog(self.file_states, tab_names, ifacemsg,
                                parent=self)
        dlg.exec()

    def hide_ch(self):
        s = self._state
        if s is None:
            return
        checkBox = self.sender()
        s.show_channels[checkBox.number] = 0 if checkBox.isChecked() else 1
        if s.abif_raw is None:
            return
        self.replot()

    def retab(self, state=None):
        self.findpeaks(state)
        s = state if state is not None else self._state
        if s is None:
            return
        if not s.readonly:
    # Identify the ILS channel (1-based index) so its peaks are labelled
    # "ILS" rather than going through allele binning.
            try:
                ils_channel = s.udatac.index(
                    size_standards[s.ILS.currentText()]['channel']) + 1
            except (ValueError, KeyError):
                ils_channel = None
    # Allele binning — only meaningful if sizes exist and a panel is selected.
    # Unmatched peaks get 'OL' (off ladder); if no panel is loaded the column
    # is left blank so it doesn't clutter the view.
            panel_name = s.panel_combo.currentText()
            if (s.panel_data and panel_name in s.panel_data
                    and len(s.peaksizes) > 0):
                s.peakalleles = assign_alleles(
                    s.peaksizes, s.peakchannels, s.panel_data[panel_name])
            else:
                s.peakalleles = [''] * len(s.peakchannels)
            # Stamp ILS peaks after binning.
            if ils_channel is not None and len(s.peaksizes) > 0:
                ils_sizes_arr = array(s.size_std, dtype=float)
                new_alleles = []
                for ch, sz, a in zip(s.peakchannels, s.peaksizes, s.peakalleles):
                    if int(ch) != ils_channel:
                        new_alleles.append(a)
                    elif npany(npabs(float(sz) - ils_sizes_arr) <= 0.5):
                        new_alleles.append('ILS')
                    else:
                        new_alleles.append('OL')
                s.peakalleles = new_alleles
            # Stutter filtering — runs after binning and ILS stamping.
            if (s.panel_data and panel_name in s.panel_data
                    and len(s.peaksizes) > 0):
                s.peakalleles = apply_stutter_filter(
                    s.peaksizes, s.peakheights, s.peakchannels, s.peakalleles,
                    s.panel_data[panel_name], s.Dye,)
        # Convert 1-based channel indices to dye names for display.
        ch_names = [s.Dye[int(ch) - 1] if 0 < int(ch) <= len(s.Dye)
                    else str(ch) for ch in s.peakchannels]
        rowcount = len(s.peakchannels)
        s.table_widget.setRowCount(rowcount)
        basic_data = [ch_names, s.peakpositions, s.peakheights,
                      s.peakfwhms, s.peakareas]
        if len(s.peaksizes) <= 0:
            basic_data.append(["NaN"]*len(s.peakchannels))
        else:
            basic_data.append(s.peaksizes)
        basic_data.append(s.peakalleles)
        s.table_widget.setData(transpose(basic_data))
        s.table_widget.setHorizontalHeaderLabels(['Peak Channel',
                                                  'Peak Position\n(Datapoints)',
                                                  'Peak Height', 'Peak FWHM',
                                                  'Peak Area\n(Datapoints)',
                                                  'Peak Size', 'Allele'])
        s.table_widget.resizeColumnsToContents()

    def process_whole_batch(self):
        source = self._state
        if source is None:
            return
        params = {
            'height': source.getheight.value(),
            'width': source.getwidth.value(),
            'prominence': source.getprominence.value(),
            'winwidth': source.getwinwidth.value(),
            'bcd': source.bcd.isChecked(),
            'ils': source.ILS.currentText(),
            'sm': source.SM.currentText(),
            'sizecall': source.sizecall.isChecked(),
            'panel': source.panel_combo.currentText(),
            'channels': list(source.show_channels),
        }
        current_idx = self.file_tab.currentIndex()
        for i, t in enumerate(self.file_states):
            if i == current_idx:
                continue
            for w in (t.getheight, t.getwidth, t.getprominence, t.getwinwidth,
                      t.bcd, t.ILS, t.SM, t.sizecall, t.panel_combo):
                w.blockSignals(True)
            for cb in t.hidech:
                cb.blockSignals(True)
            t.getheight.setValue(params['height'])
            t.getwidth.setValue(params['width'])
            t.getprominence.setValue(params['prominence'])
            t.getwinwidth.setValue(params['winwidth'])
            t.bcd.setChecked(params['bcd'])
            t.do_BCD = params['bcd']
            t.ILS.setCurrentText(params['ils'])
            t.SM.setCurrentText(params['sm'])
            t.sizecall.setChecked(params['sizecall'])
            t.panel_combo.setCurrentText(params['panel'])
            for ch_idx, cb in enumerate(t.hidech):
                hidden = not params['channels'][ch_idx]
                cb.setChecked(hidden)
                t.show_channels[ch_idx] = params['channels'][ch_idx]
            for w in (t.getheight, t.getwidth, t.getprominence, t.getwinwidth,
                      t.bcd, t.ILS, t.SM, t.sizecall, t.panel_combo):
                w.blockSignals(False)
            for cb in t.hidech:
                cb.blockSignals(False)
            self.reanalyse(t)

    def setbcd(self):
        s = self._state
        if s is None:
            return
        s.do_BCD = self.sender().isChecked()
        self.reanalyse()

    def reanalyse(self, state=None):
        if not isinstance(state, FileState):
            state = None
        s = state if state is not None else self._state
        if s is None or s.readonly or s.abif_raw is None:
            return
        s.should_sizecall = False
        if s.sizecall.isChecked():
            s.should_sizecall = True
            s.sizecall.setChecked(False)
    # If a panel is loaded, sizing must run so allele assignment has fragment
    # sizes to work with. Skip this if we are already in a recovery pass after
    # a sizing error (_sizing_recovery flag set by _sizingerror()), which would
    # otherwise cause infinite recursion.
        if s.panel_data:
            if (s.panel_combo.currentText() != ifacemsg["nopanel"]
                and not s.should_sizecall
                and not getattr(self, '_sizing_recovery', False)):
                s.should_sizecall = True
        self.retab(s)
        self.replot(s)
        s.should_sizecall = False
