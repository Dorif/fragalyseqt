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
from os.path import expanduser, dirname
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
from . import fillarray
from .fillarray import _safe_find
from .setvar import (set_dye_array, set_graph_name, set_ILS_channel,
                     set_spl_dgr, set_knots, set_lsq_ord, chk_key_valid,
                     southern_fit_local, southern_fit_global)
ftype = "ABI fragment analysis files (*.fsa *.hid)"
global show_channels, ifacemsg, do_BCD
do_BCD = False
ifacemsg = {}
localizefq(ifacemsg)
show_channels = [1] * 8
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


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        from pyqtgraph.Qt.QtWidgets import (
            QWidget,
            QPushButton,
            QLabel,
            QVBoxLayout,
            QHBoxLayout,
            QGridLayout,
            QSizePolicy,
        )
        from pyqtgraph.Qt.QtGui import QIcon
        MainWindow.setWindowTitle("FragalyseQt")
        MainWindow.setWindowIcon(QIcon("FragalyseQt.png"))
        MainWindow.resize(960, 640)
        self.centralwidget = QWidget(MainWindow)
        MainWindow.setCentralWidget(self.centralwidget)

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
        top_bar.addStretch(1)

        self.graphWidget = PlotWidget(self.centralwidget)
        self.graphWidget.setBackground(None)
        self.graphWidget.showGrid(x=True, y=True)
        self.graphWidget.setLabel("left", "Signal intensity, RFU")
        try:
            expanding_policy = QSizePolicy.Expanding
        except AttributeError:
            expanding_policy = QSizePolicy.Policy.Expanding
        self.graphWidget.setSizePolicy(expanding_policy, expanding_policy)
        root_layout.addWidget(self.graphWidget, stretch=3)

        bottom_layout = QHBoxLayout()
        bottom_layout.setSpacing(8)
        root_layout.addLayout(bottom_layout, stretch=2)

        self.fsatab = TableWidget(self.centralwidget, sortable=False)
        self.fsatab.setSizePolicy(expanding_policy, expanding_policy)
        bottom_layout.addWidget(self.fsatab, stretch=3)

        controls_widget = QWidget(self.centralwidget)
        controls_layout = QGridLayout(controls_widget)
        controls_layout.setContentsMargins(0, 0, 0, 0)
        controls_layout.setHorizontalSpacing(6)
        controls_layout.setVerticalSpacing(6)
        bottom_layout.addWidget(controls_widget, stretch=1)

        self.getheightlabel = QLabel(self)
        self.getheightlabel.setText(ifacemsg["minph"])
        self.getheightlabel.setStyleSheet(''' font-size: 10pt; ''')
        controls_layout.addWidget(self.getheightlabel, 0, 0)

        self.getheight = SpinBox(self, minStep=1, dec=True)
        self.getheight.setRange(1, 64000)
        self.getheight.setValue(175)
        self.getheight.setMinimumHeight(20)
        self.getheight.setStyleSheet(''' font-size: 8pt; ''')
        self.getheight.valueChanged.connect(self.reanalyse)
        controls_layout.addWidget(self.getheight, 0, 1)

        self.getwidthlabel = QLabel(self)
        self.getwidthlabel.setText(ifacemsg["minpw"])
        self.getwidthlabel.setStyleSheet(''' font-size: 10pt; ''')
        controls_layout.addWidget(self.getwidthlabel, 1, 0)

        self.getwidth = SpinBox(self, dec=True)
        self.getwidth.setRange(1, 16000)
        self.getwidth.setValue(4)
        self.getwidth.setMinimumHeight(20)
        self.getwidth.setStyleSheet(''' font-size: 8pt; ''')
        self.getwidth.valueChanged.connect(self.reanalyse)
        controls_layout.addWidget(self.getwidth, 1, 1)

        self.getprominencelabel = QLabel(self)
        self.getprominencelabel.setText(ifacemsg["minpp"])
        self.getprominencelabel.setStyleSheet(''' font-size: 10pt; ''')
        controls_layout.addWidget(self.getprominencelabel, 2, 0)

        self.getprominence = SpinBox(self, minStep=1, dec=True)
        self.getprominence.setRange(1, 64000)
        self.getprominence.setValue(175)
        self.getprominence.setMinimumHeight(20)
        self.getprominence.setStyleSheet(''' font-size: 8pt; ''')
        self.getprominence.valueChanged.connect(self.reanalyse)
        controls_layout.addWidget(self.getprominence, 2, 1)

        self.getwinwidthlabel = QLabel(self)
        self.getwinwidthlabel.setText(ifacemsg["minww"])
        self.getwinwidthlabel.setStyleSheet(''' font-size: 10pt; ''')
        controls_layout.addWidget(self.getwinwidthlabel, 3, 0)

        self.getwinwidth = SpinBox(self, minStep=1, dec=True)
        self.getwinwidth.setRange(1, 1000)
        self.getwinwidth.setValue(51)
        self.getwinwidth.setMinimumHeight(20)
        self.getwinwidth.setStyleSheet(''' font-size: 8pt; ''')
        self.getwinwidth.valueChanged.connect(self.reanalyse)
        controls_layout.addWidget(self.getwinwidth, 3, 1)

        self.hidech = []
        i = 0
        while i < 8:
            self.hidech.append(QCheckBox(self.centralwidget))
            self.hidech[i].setStyleSheet(''' font-size: 10pt; ''')
            self.hidech[i].toggled.connect(self.hide_ch)
            self.hidech[i].number = i
            controls_layout.addWidget(self.hidech[i], 4 + (i // 2), i % 2)
            i += 1
        self.bcd = QCheckBox(self.centralwidget)
        self.bcd.setText(ifacemsg["bcd"])
        self.bcd.toggled.connect(self.setbcd)
        self.bcd.setStyleSheet(''' font-size: 10pt; ''')
        controls_layout.addWidget(self.bcd, 8, 0, 1, 2)
        self.ILS = ComboBox(self.centralwidget)
        self.ILS.setItems(size_standards)
        self.ILS.setStyleSheet(''' font-size: 10pt; ''')
        controls_layout.addWidget(self.ILS, 9, 0, 1, 2)
        self.SM = ComboBox(self.centralwidget)
        self.SM.setItems(["Local Southern", "Global Southern",
                          "Cubic spline sizing", "Linear spline sizing",
                          "5th degree spline sizing",
                          "LSQ weighted linear spline sizing",
                          "LSQ weighted cubic spline sizing",
                          "LSQ weighted 5th degree spline sizing",
                          "LSQ 2nd order", "LSQ 3rd order", "LSQ 5th order"])
        self.SM.setStyleSheet(''' font-size: 10pt; ''')
        controls_layout.addWidget(self.SM, 10, 0)
        self.sizecall = QPushButton(self.centralwidget)
        self.sizecall.setCheckable(True)
        self.sizecall.setText("SizeCall")
        self.sizecall.setStyleSheet(''' font-size: 10pt; ''')
        self.sizecall.clicked.connect(self.reanalyse)
        controls_layout.addWidget(self.sizecall, 10, 1)

        controls_layout.setColumnStretch(0, 1)
        controls_layout.setColumnStretch(1, 0)
        bottom_layout.setStretch(0, 3)
        bottom_layout.setStretch(1, 1)
        self.inactivatechkboxes()

# Checkboxes w/o designations or with designations of nonexistent channels
# would look weird, so let's inactivate them correctly.
    def inactivatechkboxes(self):
        i = 0
        while i < 8:
            self.hidech[i].setText(ifacemsg['ch_inact_msg'])
            i += 1

    def open_and_plot(self):
        openBtn = self.sender()
        if openBtn.isChecked():
            openBtn.setChecked(False)
            global homedir, dyerange, Dye, abif_raw, udatac
            udatac = ["DATA1", "DATA2", "DATA3", "DATA4", "DATA105",
                      "DATA106", "DATA107", "DATA108"]
            wavelng = ["DyeW1", "DyeW2", "DyeW3", "DyeW4", "DyeW5", "DyeW6",
                       "DyeW7", "DyeW8"]
            dyen = ["DyeN1", "DyeN2", "DyeN3", "DyeN4", "DyeN5", "DyeN6",
                    "DyeN7", "DyeN8"]
            fname, _ = FileDialog.getOpenFileName(self,
                                                  'Open file for analysis',
                                                  homedir, ftype)
            if not fname: return
# If file open is cancelled, no error rises.
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
            if tmpabif["DATA1"] is None:
                # Assuming what it may be HID file.
                try:
                    HIDfile = open(fname, "rb")
                    s = HIDfile.read()
                    tmpkeys = tmpabif.keys()
                    MODL1_hex = b'\x4d\x4f\x44\x4c\x00\x00\x00\x01'
                    if "MODL1" in tmpkeys and tmpabif["MODL1"] is None:
                        modl1_pos = _safe_find(s, MODL1_hex)
                        HIDfile.seek(modl1_pos + 12, 0)
                        namesize = int.from_bytes(HIDfile.read(4), 'big')
                        # If we have file generated by old ABI 310 software -
                        # it may have different format of records.
                        if namesize < 4:
                            HIDfile.seek(modl1_pos + 10, 0)
                            namesize = int.from_bytes(HIDfile.read(2), 'big')
                            HIDfile.seek(8, 1)
                        else:
                            HIDfile.seek(4, 1)
                        tmpabif["MODL1"] = HIDfile.read(namesize)
                    if "Peak1" in tmpkeys and tmpabif["Peak1"] is None:
                        pshorthexarray = [b'\x50\x65\x61\x6b\x00\x00\x00\x01',
                                          b'\x50\x65\x61\x6b\x00\x00\x00\x05']
                        pinthexarray = [b'\x50\x65\x61\x6b\x00\x00\x00\x02',
                                        b'\x50\x65\x61\x6b\x00\x00\x00\x03',
                                        b'\x50\x65\x61\x6b\x00\x00\x00\x04',
                                        b'\x50\x65\x61\x6b\x00\x00\x00\x07',
                                        b'\x50\x65\x61\x6b\x00\x00\x00\x08',
                                        b'\x50\x65\x61\x6b\x00\x00\x00\x09',
                                        b'\x50\x65\x61\x6b\x00\x00\x00\x0a']
                        pdoublehexarray = [b'\x50\x65\x61\x6b\x00\x00\x00\x06',
                                           b'\x50\x65\x61\x6b\x00\x00\x00\x0b',
                                           b'\x50\x65\x61\x6b\x00\x00\x00\x0c',
                                           b'\x50\x65\x61\x6b\x00\x00\x00\x0d',
                                           b'\x50\x65\x61\x6b\x00\x00\x00\x0e',
                                           b'\x50\x65\x61\x6b\x00\x00\x00\x0f',
                                           b'\x50\x65\x61\x6b\x00\x00\x00\x10',
                                           b'\x50\x65\x61\x6b\x00\x00\x00\x11',
                                           b'\x50\x65\x61\x6b\x00\x00\x00\x12',
                                           b'\x50\x65\x61\x6b\x00\x00\x00\x15']
                        pshortname = ["Peak1", "Peak5"]
                        pintname = ["Peak2", "Peak3", "Peak4", "Peak7",
                                    "Peak8", "Peak9", "Peak10"]
                        pdoublename = ["Peak6", "Peak11", "Peak12", "Peak13",
                                       "Peak14", "Peak15", "Peak16", "Peak17",
                                       "Peak18", "Peak21"]
                        HIDfile.seek(_safe_find(s, pshorthexarray[0]) + 12, 0)
                        peakarraylen = int.from_bytes(HIDfile.read(4), 'big')
                        fillarray.fill_num_array(tmpabif, HIDfile, s,
                                                 pshortname, pshorthexarray,
                                                 peakarraylen, '>H')
                        fillarray.fill_num_array(tmpabif, HIDfile, s, pintname,
                                                 pinthexarray, peakarraylen,
                                                 '>I')
                        fillarray.fill_num_array(tmpabif, HIDfile, s,
                                                 pdoublename, pdoublehexarray,
                                                 peakarraylen, '>d')
                    if tmprecord.annotations["abif_raw"]["Dye#1"] is None:
                        dyenum = b'\x44\x79\x65\x23\x00\x00\x00\x01'
                        HIDfile.seek(_safe_find(s, dyenum) + 20, 0)
                        tmpabif["Dye#1"] = min(int.from_bytes(
                            HIDfile.read(2), 'big'), 8)
                        for i in range(tmpabif["Dye#1"]):
                            tmpabif[udatac[i]] = None
                    fillarray.fill_char_array(tmpabif, HIDfile, s, udatac,
                                              dyen, wavelng)
                    HIDfile.close()
                    abif_raw = tmpabif
                except Exception:
                    HIDfile.close()
                    # If it's not HID file and we can't obtain data - we
                    # should tell about this.
                    msgbox(ifacemsg['dmgdfile'], ifacemsg['nodatamsg'], 2)
                    try:
                        record
                    except NameError:
                        self.open_and_plot()
            else:
                abif_raw = tmprecord.annotations["abif_raw"]
# We need raw data from ABIF file only, no need in entire data structure,
# created by BioPython's AbiIO. This way multiple brackets constructions
# are evaded.
            homedir = dirname(fname)
            self.inactivatechkboxes()
            Dye = set_dye_array(abif_raw)
            dyerange = range(abif_raw["Dye#1"])
# Assuming no more than 8 dyes are met at once.
            self.reanalyse()

    def about(self):
        msgbox(ifacemsg['aboutbtn'], ifacemsg['infoboxtxt'], 0)
        self.aboutInfo.setChecked(False)

    def findpeaks(self):

        def _sizingerror():
            msgbox("", ifacemsg['wrongsizing'], 1)
            self.sizecall.setChecked(False)
            self.reanalyse()

        # Detecting peaks and calculating peaks data.
        global winwidth, peakpositions, peakheights, peakfwhms, peakchannels
        global peakareas, peaksizes, ch, x_plot, size_std, lsq_order, farr
        global issouthern
        h = self.getheight.value()
        w = self.getwidth.value()
        p = self.getprominence.value()
        winwidth = self.getwinwidth.value()
        _positions = []
        _heights = []
        _fwhms = []
        _channels = []
        _sizes = []
        ch = []
        chP = []
        farr = []
        lsq_order = 0
        issouthern = False
        x_plot = list(dict(enumerate(abif_raw["DATA1"], start=1)))
        if should_sizecall:
            spline_degree = 0
            ILS_Name = self.ILS.currentText()
            Sizing_Method = self.SM.currentText()
            try:
                ils_data = set_ILS_channel(abif_raw, ILS_Name)
                if do_BCD:
                    _, params = jbcd(ils_data,
                                     half_window=(winwidth-1)//2)
                    ils_data = params['signal']
                size_std = size_standards[ILS_Name]
                n_expected = len(size_std)
                ils_h, ils_p = h, p
                for _attempt in range(4):
                    ILSP = find_peaks(ils_data, height=ils_h, width=w,
                                      prominence=ils_p, wlen=winwidth,
                                      rel_height=0.5)
                    if len(ILSP[0]) >= n_expected:
                        break
                    ils_h = max(ils_h // 2, 1)
                    ils_p = max(ils_p // 2, 1)
                beginning_index = len(ILSP[0]) - n_expected
                ladder_peaks = _refine_peak_positions(
                    ils_data, ILSP[0][beginning_index:])
                if 'spline' in Sizing_Method:
                    spline_degree = set_spl_dgr(Sizing_Method)
                    knots = set_knots(Sizing_Method, ladder_peaks,
                                      spline_degree)
                if spline_degree != 0:
                    spline = splrep(ladder_peaks, size_std, k=spline_degree,
                                    t=knots)
                    x_plot = around(splev(x_plot, spline), 3)
                elif 'order' in Sizing_Method:
                    lsq_order = set_lsq_ord(Sizing_Method)
                    func = Polynomial.fit(ladder_peaks, size_std,
                                          lsq_order)
                    x_plot = around(func(array(x_plot)), 3)
                elif 'Southern' in Sizing_Method:
                    issouthern = True
                    southern_func = (southern_fit_local
                                     if 'Local' in Sizing_Method
                                     else southern_fit_global)
                    x_plot = around(southern_func(
                        ladder_peaks, size_std, x_plot), 3)
            except (ValueError, TypeError, KeyError):
                _sizingerror()
        # By default, find_peaks function measures width at
        # half maximum of height (rel_height=0.5). But
        # explicit is always better, then implicit, so
        # rel_height is specified clearly.
        if do_BCD:
            half_win = (winwidth-1)//2
            def _bcd_channel(chnum):
                _, params = jbcd(abif_raw[udatac[chnum]],
                                 half_window=half_win)
                return list(params['signal'])
            with ThreadPoolExecutor() as executor:
                ch = list(executor.map(_bcd_channel, dyerange))
        else:
            for chnum in dyerange:
                ch.append(list(abif_raw[udatac[chnum]]))
        def _detect_peaks(chnum):
            return find_peaks(ch[chnum], height=h, width=w, prominence=p,
                              wlen=winwidth, rel_height=0.5)
        with ThreadPoolExecutor() as executor:
            chP = list(executor.map(_detect_peaks, dyerange))
        for chnum in dyerange:
            refined_pos = _refine_peak_positions(ch[chnum], chP[chnum][0])
            _positions.append(refined_pos)
            _heights.append(chP[chnum][1]['peak_heights'])
            _fwhms.append(chP[chnum][1]['widths'])
            if should_sizecall and len(refined_pos) != 0:
                if spline_degree != 0:
                    _sizes.append(splev(refined_pos, spline))
                elif issouthern:
                    _sizes.append(southern_func(
                        ladder_peaks, size_std, refined_pos))
                else:
                    _sizes.append(func(refined_pos))
            _channels.append([Dye[chnum]]*len(chP[chnum][0]))
        peakpositions = concatenate(_positions) if _positions else array([])
        peakheights = concatenate(_heights) if _heights else array([])
        peakfwhms = concatenate(_fwhms) if _fwhms else array([])
        peakchannels = concatenate(_channels) if _channels else array([])
        peaksizes = concatenate(_sizes) if _sizes else array([])
        if len(peaksizes) > 0:
            valid = peaksizes >= 0
            peakpositions = peakpositions[valid]
            peakheights = peakheights[valid]
            peakfwhms = peakfwhms[valid]
            peakchannels = peakchannels[valid]
            peaksizes = peaksizes[valid]
# Calculate areas from full-precision values, then round everything.
        peakareas = around(multiply(peakheights, peakfwhms)*1.0645, 2)
        peaksizes = around(peaksizes, 2)
        peakheights = around(peakheights, 2)
        peakfwhms = around(peakfwhms, 2)

# Peak areas are calculated using formula for Gaussian peak area
# (https://www.physicsforums.com/threads/area-under-gaussian-peak-by-easy-measurements.419285/):
# A = FWHM*H/(2sqrt(2ln(2))/sqrt(2pi)) = 1.0645*FWHM*H, where FWHM is Full
# Width at Half Maximum. Real area may differ for non-Gaussian peaks, but at
# least majority of them are of Gaussian shape. If peaks are well separated -
# just calculate their area, but if your peaks are crowded (e.g. in TP-PCR or
# allelic ladders), oversaturated or you have noisy data - you MUST use
# baseline correction and denoising prior peak area calculation.

    def replot(self):
        self.graphWidget.clear()
        self.graphWidget.plotItem.setLimits(xMin=None, xMax=None, yMin=None, yMax=None)
        for i in dyerange:
            self.hidech[i].setText(ifacemsg['hidechannel'] + Dye[i])
        self.graphWidget.setTitle(set_graph_name(abif_raw), color="c",
                                  size="10pt")
        max_x = len(x_plot)
        if should_sizecall or len(peaksizes) > 0:
            # In the most normal case if you have good overall CE data quality,
            # the last member of x_plot array should have the biggest size.
            x_max = x_plot[len(x_plot)-1]
            self.graphWidget.setLabel('bottom', 'Size, bases')
            max_ladder = max(size_std)
            if max_ladder+200 < x_max:
                max_x = max_ladder+200
            # EXTREMELY weird situation, but it sometimes happens, e.g. for low
            # CE data quality with a lot of noise or poorly injected ladder.
            elif max_ladder > x_max:
                max_x = max_ladder
            else:
                max_x = x_max
        else:
            self.graphWidget.setLabel('bottom', 'Size, data points')
        max_y = 0
        for i in dyerange:
            if show_channels[i]:
                if should_sizecall or len(peaksizes) > 0:
                    ch_arr = array(ch[i])
                    valid = array(x_plot) >= 0
                    ch_max = float(ch_arr[valid].max()) if valid.any() else 0
                else:
                    ch_max = max(ch[i])
                if ch_max > max_y:
                    max_y = ch_max
                self.graphWidget.plot(x_plot, ch[i], pen=_PEN_COLORS[i])
        if max_y == 0:
            max_y = 64000
        self.graphWidget.plotItem.setLimits(xMin=0, xMax=max_x, yMin=0, yMax=max_y)

    def export_csv(self):
        # Exporting CSV with data generated by findpeaks().
        try:
            peakchannels
        except NameError:
            return
        expbox = self.sender()
        header = ['Peak Channel', 'Peak Position (Datapoints)', 'Peak Height',
                  'Peak FWHM', 'Peak Area (Datapoints)']
        do_export = False
        if expbox.focusWidget().objectName() == "CSV":
            pdarray = [peakchannels, peakpositions, peakheights, peakfwhms,
                       peakareas]
            if len(peaksizes) > 0:
                pdarray.append(peaksizes)
                header += ['Peak Size (Bases)']
            do_export = True
        elif (expbox.focusWidget().objectName() == "IA" and
              chk_key_valid("Peak1", abif_raw)):
            # Exporting internal analysis data, but first checking if file has
            # them, assuming if Peak1 field is valid, other fields are too.
            peak_chn = []
            for channel in abif_raw["Peak1"]:
                idx = min(max(channel - 1, 0), len(Dye) - 1)
                peak_chn.append(Dye[idx])
            pdarray = [peak_chn, abif_raw["Peak2"], abif_raw["Peak7"],
                       abif_raw["Peak5"], abif_raw["Peak10"],
                       abif_raw["Peak12"], abif_raw["Peak17"]]
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

    def hide_ch(self):
        checkBox = self.sender()
        if checkBox.isChecked():
            show_channels[checkBox.number] = 0
        else:
            show_channels[checkBox.number] = 1
        try:
            abif_raw
        except NameError:
            return
        self.replot()

    def retab(self):
        self.findpeaks()
        rowcount = len(peakchannels)
        self.fsatab.setRowCount(rowcount)
        basic_data = [peakchannels, peakpositions, peakheights,
                      peakfwhms, peakareas]
        if len(peaksizes) <= 0:
            basic_data.append(["NaN"]*len(peakchannels))
        else:
            basic_data.append(peaksizes)
        self.fsatab.setData(transpose(basic_data))
        self.fsatab.setHorizontalHeaderLabels(['Peak Channel',
                                               'Peak Position\n(Datapoints)',
                                               'Peak Height', 'Peak FWHM',
                                               'Peak Area\n(Datapoints)',
                                               'Peak Size'])
        self.fsatab.resizeColumnsToContents()

    def setbcd(self):
        checkBox = self.sender()
        global do_BCD
        if checkBox.isChecked():
            do_BCD = True
        else:
            do_BCD = False
        self.reanalyse()

    def reanalyse(self):
        try:
            abif_raw
        except NameError:
            return
        global should_sizecall
        should_sizecall = False
        if self.sizecall.isChecked():
            should_sizecall = True
            self.sizecall.setChecked(False)
        self.retab()
        self.replot()
        should_sizecall = False
