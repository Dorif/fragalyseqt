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
from .codisexport import CODISExportDialog
from .stutterfilter import apply_stutter_filter
from os.path import expanduser, dirname, basename
from csv import writer as csvwriter
from concurrent.futures import ThreadPoolExecutor
from Bio.SeqIO import read as fsaread
from numpy import around, multiply, array, concatenate, transpose, where
from scipy.signal import find_peaks
from scipy.interpolate import splrep, splev
from numpy.polynomial.polynomial import Polynomial
from pybaselines.morphological import jbcd
# Using FileDialog and SpinBox from pyqtgraph to prevent some possible problems
# for macOS users and to allow more fine variable setting.
from pyqtgraph import PlotWidget, FileDialog, SpinBox, ComboBox, TableWidget
# Using pyqtgraph widgets to make program independent from
# Qt for Python implementation.
from pyqtgraph.Qt.QtWidgets import QCheckBox
from .sizestandards import size_standards
from . import fillhid
from .setvar import (set_dye_array, set_graph_name,
                     set_spl_dgr, set_knots, set_lsq_ord, chk_key_valid,
                     southern_fit_local, southern_fit_global)
from .panelparser import (parse_genemapper, parse_genemarker, parse_osiris,
                          assign_alleles, _xml_root_tag)
ftype = "ABI fragment analysis files (*.fsa *.hid);;"
ftype += "Native Nanophore files (*.frf)"
ifacemsg = {}
localizefq(ifacemsg)
homedir = expanduser('~')
_PEN_COLORS = ('b', 'g', 'y', 'r', 'orange', 'c', 'm', 'k')



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
    delta = where(valid, 0.5 * (y_left - y_right) / where(valid, denom, 1.0),
                  0.0)
    refined[mask] = pos + delta
    return refined


class FileState:
    """Holds all per-file/per-tab analysis state and widget references."""
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
        self.peakalleles = []       # allele labels from panel binning
        self.panel_data = {}        # loaded panel for this tab
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
        self.load_panel_btn = None


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        from pyqtgraph.Qt.QtWidgets import (
            QWidget,
            QPushButton,
            QVBoxLayout,
            QHBoxLayout,
            QSizePolicy,
            QTabWidget,
        )
        from pyqtgraph.Qt.QtGui import QIcon
        MainWindow.setWindowTitle("FragalyseQt")
        MainWindow.setWindowIcon(QIcon("FragalyseQt.png"))
        MainWindow.resize(960, 640)
        self.centralwidget = QWidget(MainWindow)
        MainWindow.setCentralWidget(self.centralwidget)

        self.file_states = []

        root_layout = QVBoxLayout(self.centralwidget)
        root_layout.setContentsMargins(8, 8, 8, 8)
        root_layout.setSpacing(6)

        top_bar = QHBoxLayout()
        top_bar.setSpacing(6)
        root_layout.addLayout(top_bar)

        self.openFSA = QPushButton(self.centralwidget)
        self.openFSA.setCheckable(True)
        self.openFSA.setText(ifacemsg["openfragmentfile"])
        self.openFSA.setShortcut("Ctrl+O")
        self.openFSA.clicked.connect(self.open_and_plot)
        self.openFSA.setMinimumWidth(120)
        top_bar.addWidget(self.openFSA)

        self.closeTab = QPushButton(self.centralwidget)
        self.closeTab.setText(ifacemsg["closetab"])
        self.closeTab.setShortcut("Ctrl+W")
        self.closeTab.clicked.connect(self._close_tab_action)
        self.closeTab.setMinimumWidth(90)
        top_bar.addWidget(self.closeTab)

        self.aboutInfo = QPushButton(self.centralwidget)
        self.aboutInfo.setCheckable(True)
        self.aboutInfo.setText(ifacemsg["aboutbtn"])
        self.aboutInfo.setShortcut("F1")
        self.aboutInfo.clicked.connect(self.about)
        self.aboutInfo.setMinimumWidth(100)
        top_bar.addWidget(self.aboutInfo)

        self.exportInternalAnalysisData = QPushButton(self.centralwidget)
        self.exportInternalAnalysisData.setText(ifacemsg["exportinternal"])
        self.exportInternalAnalysisData.setShortcut("Ctrl+I")
        self.exportInternalAnalysisData.setObjectName("IA")
        self.exportInternalAnalysisData.clicked.connect(self.export_csv)
        self.exportInternalAnalysisData.setMinimumWidth(260)
        top_bar.addWidget(self.exportInternalAnalysisData)

        self.exportCSV = QPushButton(self.centralwidget)
        self.exportCSV.setText(ifacemsg["csvexport"])
        self.exportCSV.setShortcut("Ctrl+E")
        self.exportCSV.setObjectName("CSV")
        self.exportCSV.clicked.connect(self.export_csv)
        self.exportCSV.setMinimumWidth(120)
        top_bar.addWidget(self.exportCSV)

        self.exportCODIS = QPushButton(self.centralwidget)
        self.exportCODIS.setText(ifacemsg["codisexport"])
        self.exportCODIS.setShortcut("Ctrl+Shift+C")
        self.exportCODIS.clicked.connect(self.export_codis)
        self.exportCODIS.setMinimumWidth(160)
        top_bar.addWidget(self.exportCODIS)

        top_bar.addStretch(1)

        self.file_tab = QTabWidget(self.centralwidget)
        root_layout.addWidget(self.file_tab, stretch=1)

    @property
    def _state(self):
        idx = self.file_tab.currentIndex()
        if 0 <= idx < len(self.file_states):
            return self.file_states[idx]
        return None

    def _create_tab_content(self, state):
        """Creates the per-tab widget (plot + table + controls) and stores
        widget references in the given FileState."""
        from pyqtgraph.Qt.QtWidgets import (
            QWidget,
            QPushButton,
            QLabel,
            QVBoxLayout,
            QHBoxLayout,
            QGridLayout,
            QSizePolicy,
            QDoubleSpinBox,
        )
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
        plot.setBackground(None)
        plot.showGrid(x=True, y=True)
        plot.setLabel("left", "Signal intensity, RFU")
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
                           "Cubic spline sizing", "Linear spline sizing",
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

        load_panel_btn = QPushButton()
        load_panel_btn.setText(ifacemsg["loadpanel"])
        load_panel_btn.setStyleSheet(''' font-size: 10pt; ''')
        load_panel_btn.clicked.connect(self.load_panel_action)
        controls_layout.addWidget(load_panel_btn, 11, 0, 1, 2)

        panel_combo = ComboBox()
        panel_combo.setItems([ifacemsg["nopanel"]])
        panel_combo.setEnabled(False)
        panel_combo.setStyleSheet(''' font-size: 10pt; ''')
        panel_combo.currentIndexChanged.connect(self.reanalyse)
        controls_layout.addWidget(panel_combo, 12, 0, 1, 2)

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
        state.load_panel_btn = load_panel_btn
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

    def load_panel_action(self):
        """Open a file dialog to load a GeneMapper or GeneMarker panel file
        for the current tab.

        GeneMapper (.txt): after selecting the Panels file, immediately prompts
        for the companion Bins file.  Bins are optional — the user can decline
        and only marker-range annotation will be available.
        GeneMarker (.xml): self-contained; bins are embedded, no extra dialog.
        """
        s = self._state
        global homedir
        path, _ = FileDialog.getOpenFileName(
            self, ifacemsg['loadpaneldlg'], homedir,
            "Panel files (*.txt *.xml)"
        )
        if not path:
            return
        try:
            if path.lower().endswith('.xml'):
                if _xml_root_tag(path) == 'KitData':
                    data = parse_osiris(path)
                else:
                    data = parse_genemarker(path)
            else:
                # Load panels; immediately ask for the bins file.
                data = parse_genemapper(path, '')
                bins_path, _ = FileDialog.getOpenFileName(
                    self, ifacemsg['loadbinsdlg'], dirname(path),
                    "Bins files (*.txt)"
                )
                if bins_path:
                    data = parse_genemapper(path, bins_path)
                else:
                    msgbox("", ifacemsg['nobinsmsg'], 0)
                # Ask for the optional stutter file.
                stutter_path, _ = FileDialog.getOpenFileName(
                    self,
                    ifacemsg.get('loadstutterdlg',
                                 'Select stutter ratios file (optional)'),
                    dirname(path), "Stutter files (*.txt)"
                )
                if stutter_path:
                    data = parse_genemapper(
                        path,
                        bins_path if bins_path else '',
                        stutter_path,
                    )
        except Exception as exc:
            msgbox("", str(exc), 2)
            return
        if not data:
            msgbox("", ifacemsg['nodatamsg'], 1)
            return
        if s is not None:
            s.panel_data = data
            s.panel_combo.setItems(list(data.keys()))
            s.panel_combo.setEnabled(True)

    def open_and_plot(self):
        openBtn = self.sender()
        if openBtn.isChecked():
            openBtn.setChecked(False)
            global homedir
            udatac = fillhid.UDATAC
            fnames, _ = FileDialog.getOpenFileNames(self,
                                                   'Open files for analysis',
                                                   homedir, ftype)
            if not fnames: return
# If file open is cancelled, no error rises.
            for fname in fnames:
                if fname.lower().endswith('.frf'):
                    from . import fillfrf
                    try:
                        abif_result = fillfrf.parse_frf(fname)
                    except Exception:
                        msgbox(ifacemsg['dmgdfile'], ifacemsg['nodatamsg'], 2)
                        continue
                    homedir = dirname(fname)
                    state = FileState()
                    state.abif_raw = abif_result
                    state.udatac = udatac
                    state.Dye = set_dye_array(abif_result)
                    state.dyerange = range(abif_result["Dye#1"])
                    tab_widget = self._create_tab_content(state)
                    self.file_states.append(state)
                    self.file_tab.addTab(tab_widget, basename(fname))
                    self.file_tab.setCurrentIndex(len(self.file_states) - 1)
                    self.reanalyse()
                    continue
                FAfile = open(fname, "rb")
                try:
                    tmprecord = fsaread(FAfile, "abi")
                except AssertionError:
                    class record():
                        annotations = {"abif_raw": {
                            "DATA1": None,
                            "DATA2": None,
                            "DATA3": None,
                            "DATA4": None,
                            "Dye#1": None,
                            "DyeN1": None,
                            "DyeN2": None,
                            "DyeN3": None,
                            "DyeN4": None,
                            "MODL1": None}}
                    tmprecord = record()
# Preventing data corruption in a case if target file is corrupted.
                FAfile.close()
# Closing file to save memory and avoid unexpected things.
                tmpabif = tmprecord.annotations["abif_raw"]
                abif_result = None
                if tmpabif["DATA1"] is None:
                    # Assuming what it may be HID file.
                    try:
                        fillhid.parse_hid(fname, tmpabif, ifacemsg)
                    except Exception:
                        msgbox(ifacemsg['dmgdfile'], ifacemsg['nodatamsg'], 2)
                        continue
                abif_result = tmpabif
# We need raw data from ABIF file only, no need in entire data structure,
# created by BioPython's AbiIO. This way multiple brackets constructions
# are evaded.
                homedir = dirname(fname)
                state = FileState()
                state.abif_raw = abif_result
                state.udatac = udatac
                state.Dye = set_dye_array(abif_result)
                state.dyerange = range(abif_result["Dye#1"])
                tab_widget = self._create_tab_content(state)
                self.file_states.append(state)
                self.file_tab.addTab(tab_widget, basename(fname))
                self.file_tab.setCurrentIndex(len(self.file_states) - 1)
                self.reanalyse()

    def about(self):
        msgbox(ifacemsg['aboutbtn'], ifacemsg['infoboxtxt'], 0)
        self.aboutInfo.setChecked(False)

    def findpeaks(self):
        s = self._state
        if s is None:
            return

        def _sizingerror():
            msgbox("", ifacemsg['wrongsizing'], 1)
            s.sizecall.setChecked(False)
            s.should_sizecall = False
            # Guard: tell reanalyse() not to re-trigger the auto-sizing that
            # just failed, otherwise panel_data being set would cause a loop.
            self._sizing_recovery = True
            self.reanalyse()
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
        s.x_plot = list(dict(enumerate(s.abif_raw["DATA1"], start=1)))
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
                beginning_index = len(ILSP[0]) - n_expected
                ladder_peaks = _refine_peak_positions(
                    ils_data, ILSP[0][beginning_index:])
                if 'spline' in Sizing_Method:
                    spline_degree = set_spl_dgr(Sizing_Method)
                    knots = set_knots(Sizing_Method, ladder_peaks,
                                      spline_degree)
                if spline_degree != 0:
                    spline = splrep(ladder_peaks, s.size_std, k=spline_degree,
                                    t=knots)
                    s.x_plot = around(splev(s.x_plot, spline), 3)
                elif 'order' in Sizing_Method:
                    s.lsq_order = set_lsq_ord(Sizing_Method)
                    func = Polynomial.fit(ladder_peaks, s.size_std,
                                          s.lsq_order)
                    s.x_plot = around(func(array(s.x_plot)), 3)
                elif 'Southern' in Sizing_Method:
                    s.issouthern = True
                    southern_func = (southern_fit_local
                                     if 'Local' in Sizing_Method
                                     else southern_fit_global)
                    s.x_plot = around(southern_func(
                        ladder_peaks, s.size_std, s.x_plot), 3)
            except (ValueError, TypeError, KeyError):
                _sizingerror()
                return
        # By default, find_peaks function measures width at
        # half maximum of height (rel_height=0.5). But
        # explicit is always better, then implicit, so
        # rel_height is specified clearly.
        if s.do_BCD:
            half_win = (s.winwidth-1)//2
            def _bcd_channel(chnum):
                _, params = jbcd(s.abif_raw[s.udatac[chnum]],
                                 half_window=half_win)
                return list(params['signal'])
            with ThreadPoolExecutor() as executor:
                s.ch = list(executor.map(_bcd_channel, s.dyerange))
        else:
            for chnum in s.dyerange:
                s.ch.append(list(s.abif_raw[s.udatac[chnum]]))
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

    def replot(self):
        s = self._state
        if s is None:
            return
        s.plot_widget.clear()
        s.plot_widget.plotItem.setLimits(xMin=None, xMax=None, yMin=None,
                                         yMax=None)
        for i in s.dyerange:
            s.hidech[i].setText(ifacemsg['hidechannel'] + s.Dye[i])
        s.plot_widget.setTitle(set_graph_name(s.abif_raw), color="c",
                               size="10pt")
        max_x = len(s.x_plot)
        if s.should_sizecall or len(s.peaksizes) > 0:
            # In the most normal case if you have good overall CE data quality,
            # the last member of x_plot array should have the biggest size.
            x_max = s.x_plot[len(s.x_plot)-1]
            s.plot_widget.setLabel('bottom', 'Size, bases')
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
            s.plot_widget.setLabel('bottom', 'Size, data points')
        max_y = 0
        for i in s.dyerange:
            if s.show_channels[i]:
                if s.should_sizecall or len(s.peaksizes) > 0:
                    ch_arr = array(s.ch[i])
                    valid = array(s.x_plot) >= 0
                    ch_max = float(ch_arr[valid].max()) if valid.any() else 0
                else:
                    ch_max = max(s.ch[i])
                if ch_max > max_y:
                    max_y = ch_max
                s.plot_widget.plot(s.x_plot, s.ch[i], pen=_PEN_COLORS[i])
        if max_y == 0:
            max_y = 64000
        s.plot_widget.plotItem.setLimits(xMin=0, xMax=max_x, yMin=0,
                                         yMax=max_y)

    def export_csv(self):
        # Exporting CSV with data generated by findpeaks().
        s = self._state
        if s is None or s.abif_raw is None:
            return
        expbox = self.sender()
        header = ['Peak Channel', 'Peak Position (Datapoints)', 'Peak Height',
                  'Peak FWHM', 'Peak Area (Datapoints)']
        do_export = False
        if expbox.focusWidget().objectName() == "CSV":
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
            do_export = True
        elif (expbox.focusWidget().objectName() == "IA" and
              chk_key_valid("Peak1", s.abif_raw)):
            # Exporting internal analysis data, but first checking if file has
            # them, assuming if Peak1 field is valid, other fields are too.
            peak_chn = []
            for channel in s.abif_raw["Peak1"]:
                idx = min(max(channel - 1, 0), len(s.Dye) - 1)
                peak_chn.append(s.Dye[idx])
            pdarray = [peak_chn, s.abif_raw["Peak2"], s.abif_raw["Peak7"],
                       s.abif_raw["Peak5"], s.abif_raw["Peak10"],
                       s.abif_raw["Peak12"], s.abif_raw["Peak17"]]
            header += ['Peak Size (Bases)', 'Peak Area (Bases)']
            do_export = True
        else:
            msgbox(ifacemsg['unsuppeq'], ifacemsg['unsuppeqmsg'], 1)
        if do_export:
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
                                parent=self.exportCODIS)
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

    def retab(self):
        self.findpeaks()
        s = self._state
        if s is None:
            return

        # Identify the ILS channel (1-based index) so its peaks are labelled
        # "ILS" rather than going through allele binning.
        ils_channel = None
        if len(s.peaksizes) > 0:
            try:
                ils_channel = s.udatac.index(
                    size_standards[s.ILS.currentText()]['channel']) + 1
            except ValueError:
                pass

        # Allele binning — only meaningful when sizes exist and a panel is
        # selected.  Unmatched peaks get 'OL' (out of ladder); when no panel
        # is loaded the column is left blank so it doesn't clutter the view.
        panel_name = s.panel_combo.currentText()
        if (s.panel_data and panel_name in s.panel_data
                and len(s.peaksizes) > 0):
            s.peakalleles = assign_alleles(
                s.peaksizes, s.peakchannels, s.panel_data[panel_name])
        else:
            s.peakalleles = [''] * len(s.peakchannels)

        # Stamp ILS peaks after binning.
        # Only peaks at the ILS channel whose sized value matches a known
        # ladder fragment (within 0.5 bp) receive "ILS".  Peaks at the ILS
        # channel that do not match any ladder size are labelled "OL".
        if ils_channel is not None and len(s.size_std) > 0:
            ils_sizes = s.size_std
            new_alleles = []
            for ch, sz, a in zip(s.peakchannels, s.peaksizes, s.peakalleles):
                if int(ch) != ils_channel:
                    new_alleles.append(a)
                elif any(abs(float(sz) - ref) <= 0.5 for ref in ils_sizes):
                    new_alleles.append('ILS')
                else:
                    new_alleles.append('OL')
            s.peakalleles = new_alleles

        # Stutter filtering — runs after binning and ILS stamping so that
        # only properly labelled integer alleles are tested.
        if (s.panel_data and panel_name in s.panel_data
                and len(s.peaksizes) > 0):
            s.peakalleles = apply_stutter_filter(
                s.peaksizes, s.peakheights, s.peakchannels, s.peakalleles,
                s.panel_data[panel_name], s.Dye,
            )

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
                                                  'Peak Size',
                                                  'Allele'])
        s.table_widget.resizeColumnsToContents()

    def setbcd(self):
        s = self._state
        if s is None:
            return
        s.do_BCD = self.sender().isChecked()
        self.reanalyse()

    def reanalyse(self):
        s = self._state
        if s is None or s.abif_raw is None:
            return
        s.should_sizecall = False
        if s.sizecall.isChecked():
            s.should_sizecall = True
            s.sizecall.setChecked(False)
        # When a panel is loaded, sizing must run so allele assignment has
        # fragment sizes to work with.  Skip this if we are already in a
        # recovery pass after a sizing error (_sizing_recovery flag set by
        # _sizingerror()), which would otherwise cause infinite recursion.
        if (s.panel_data and not s.should_sizecall
                and not getattr(self, '_sizing_recovery', False)):
            s.should_sizecall = True
        self.retab()
        self.replot()
        s.should_sizecall = False
