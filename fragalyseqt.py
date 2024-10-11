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

from boxes import msgbox
from localize import localizefq
from os.path import expanduser, dirname
# Using FileDialog and SpinBox from pyqtgraph to prevent some possible problems
# for macOS users and to allow more fine variable setting.
from pyqtgraph import PlotWidget, FileDialog, SpinBox, ComboBox, TableWidget
# Using pyqtgraph widgets to make program independent from
# Qt for Python implementation.
from pyqtgraph.Qt.QtWidgets import QCheckBox
from sizestandards import size_standards
ftype = "ABI fragment analysis files (*.fsa *.hid)"
global show_channels, ifacemsg, do_BCD
do_BCD = False
ifacemsg = {}
localizefq(ifacemsg)
show_channels = [1] * 8
homedir = expanduser('~')


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        from pyqtgraph.Qt.QtWidgets import QWidget, QPushButton, QLabel
        from pyqtgraph.Qt.QtGui import QIcon
        MainWindow.setWindowTitle("FragalyseQt")
        MainWindow.setWindowIcon(QIcon('FragalyseQt.png'))
        MainWindow.resize(1024, 640)
        self.centralwidget = QWidget(MainWindow)
        MainWindow.setCentralWidget(self.centralwidget)
        self.openFSA = QPushButton(self.centralwidget)
        self.openFSA.setGeometry(0, 0, 120, 20)
        self.openFSA.setCheckable(True)
        self.openFSA.setText(ifacemsg['openfragmentfile'])
        self.openFSA.clicked.connect(self.open_and_plot)
        self.aboutInfo = QPushButton(self.centralwidget)
        self.aboutInfo.setGeometry(120, 0, 100, 20)
        self.aboutInfo.setCheckable(True)
        self.aboutInfo.setText(ifacemsg['aboutbtn'])
        self.aboutInfo.clicked.connect(self.about)
        self.exportInternalAnalysisData = QPushButton(self.centralwidget)
        self.exportInternalAnalysisData.setGeometry(220, 0, 260, 20)
        self.exportInternalAnalysisData.setText(ifacemsg['exportinternal'])
        self.exportInternalAnalysisData.setObjectName("IA")
        self.exportInternalAnalysisData.clicked.connect(self.export_csv)
        self.exportCSV = QPushButton(self.centralwidget)
        self.exportCSV.setGeometry(480, 0, 120, 20)
        self.exportCSV.setText(ifacemsg['csvexport'])
        self.exportCSV.setObjectName("CSV")
        self.exportCSV.clicked.connect(self.export_csv)
        self.graphWidget = PlotWidget(self.centralwidget)
        self.graphWidget.setGeometry(0, 20, 1024, 360)
        self.graphWidget.setBackground(None)
        self.graphWidget.showGrid(x=True, y=True)
        self.graphWidget.setLabel('left', 'Signal intensity, RFU')
        self.fsatab = TableWidget(self.centralwidget, sortable=False)
        self.fsatab.setGeometry(0, 380, 724, 260)
        self.getheightlabel = QLabel(self)
        self.getheightlabel.setGeometry(724, 380, 236, 20)
        self.getheightlabel.setText(ifacemsg['minph'])
        self.getheightlabel.setStyleSheet(''' font-size: 10px; ''')
        self.getheight = SpinBox(self, minStep=1, dec=True)
        self.getheight.setGeometry(960, 380, 64, 20)
        self.getheight.setRange(1, 64000)
        self.getheight.setValue(175)
        self.getheight.setStyleSheet(''' font-size: 10px; ''')
        self.getheight.valueChanged.connect(self.reanalyse)
        self.getwidthlabel = QLabel(self)
        self.getwidthlabel.setGeometry(724, 400, 236, 20)
        self.getwidthlabel.setText(ifacemsg['minpw'])
        self.getwidthlabel.setStyleSheet(''' font-size: 10px; ''')
        self.getwidth = SpinBox(self, dec=True)
        self.getwidth.setGeometry(960, 400, 64, 20)
        self.getwidth.setRange(1, 16000)
        self.getwidth.setValue(2)
        self.getwidth.setStyleSheet(''' font-size: 10px; ''')
        self.getwidth.valueChanged.connect(self.reanalyse)
        self.getprominencelabel = QLabel(self)
        self.getprominencelabel.setGeometry(724, 420, 236, 20)
        self.getprominencelabel.setText(ifacemsg['minpp'])
        self.getprominencelabel.setStyleSheet(''' font-size: 10px; ''')
        self.getprominence = SpinBox(self, minStep=1, dec=True)
        self.getprominence.setGeometry(960, 420, 64, 20)
        self.getprominence.setRange(1, 64000)
        self.getprominence.setValue(175)
        self.getprominence.setStyleSheet(''' font-size: 10px; ''')
        self.getprominence.valueChanged.connect(self.reanalyse)
        self.getwinwidthlabel = QLabel(self)
        self.getwinwidthlabel.setGeometry(724, 440, 236, 20)
        self.getwinwidthlabel.setText(ifacemsg['minww'])
        self.getwinwidthlabel.setStyleSheet(''' font-size: 10px; ''')
        self.getwinwidth = SpinBox(self, minStep=1, dec=True)
        self.getwinwidth.setGeometry(960, 440, 64, 20)
        self.getwinwidth.setRange(1, 1000)
        self.getwinwidth.setValue(15)
        self.getwinwidth.setStyleSheet(''' font-size: 10px; ''')
        self.getwinwidth.valueChanged.connect(self.reanalyse)
        self.hidech = []
        i = 0
        while i < 8:
            self.hidech.append(QCheckBox(self.centralwidget))
            if i % 2 == 0:
                self.hidech[i].setGeometry(724, 460+20*(i//2), 150, 20)
            else:
                self.hidech[i].setGeometry(874, 460+20*(i//2), 150, 20)
            self.hidech[i].setStyleSheet(''' font-size: 10px; ''')
            self.hidech[i].toggled.connect(self.hide_ch)
            self.hidech[i].number = i
            i += 1
        self.bcd = QCheckBox(self.centralwidget)
        self.bcd.setGeometry(724, 540, 300, 20)
        self.bcd.setText(ifacemsg['bcd'])
        self.bcd.toggled.connect(self.setbcd)
        self.bcd.setStyleSheet(''' font-size: 10px; ''')
        self.ILS = ComboBox(self.centralwidget)
        self.ILS.setGeometry(724, 560, 300, 20)
        self.ILS.setItems(size_standards)
        self.ILS.setStyleSheet(''' font-size: 10px; ''')
        self.SM = ComboBox(self.centralwidget)
        self.SM.setGeometry(724, 580, 216, 20)
        self.SM.setItems(['Cubic spline sizing', 'Linear spline sizing',
                          '5th degree spline sizing',
                          'LSQ weighted linear spline sizing',
                          'LSQ weighted cubic spline sizing',
                          'LSQ weighted 5th degree spline sizing',
                          'LSQ 2nd order', 'LSQ 3rd order', 'LSQ 5th order'])
        self.SM.setStyleSheet(''' font-size: 10px; ''')
        self.sizecall = QPushButton(self.centralwidget)
        self.sizecall.setGeometry(940, 580, 84, 20)
        self.sizecall.setCheckable(True)
        self.sizecall.setText('SizeCall')
        self.sizecall.setStyleSheet(''' font-size: 10px; ''')
        self.sizecall.clicked.connect(self.reanalyse)
        self.inactivatechkboxes()

# Checkboxes w/o designations or with designations of nonexistent channels
# would look weird, so let's inactivate them correctly.
    def inactivatechkboxes(self):
        i = 0
        while i < 8:
            self.hidech[i].setText(ifacemsg['ch_inact_msg'])
            i += 1

    def open_and_plot(self):
        from Bio.SeqIO import read as fsaread
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
                    import fillarray
                    HIDfile = open(fname, "rb")
                    s = HIDfile.read()
                    tmpkeys = tmpabif.keys()
                    MODL1_hex = b'\x4d\x4f\x44\x4c\x00\x00\x00\x01'
                    if "MODL1" in tmpkeys and tmpabif["MODL1"] is None:
                        HIDfile.seek(s.find(MODL1_hex) + 12, 0)
                        namesize = int.from_bytes(HIDfile.read(4), 'big')
                        # If we have file generated by old ABI 310 software -
                        # it may have different format of records.
                        if namesize < 4:
                            HIDfile.seek(s.find(MODL1_hex) + 10, 0)
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
                        HIDfile.seek(s.find(pshorthexarray[0]) + 12, 0)
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
                        HIDfile.seek(s.find(dyenum) + 20, 0)
                        tmpabif["Dye#1"] = int.from_bytes(HIDfile.read(2),
                                                          'big')
                        for i in range(tmpabif["Dye#1"]):
                            tmpabif[udatac[i]] = None
                    fillarray.fill_char_array(tmpabif, HIDfile, s, udatac,
                                              dyen, wavelng)
                    HIDfile.close()
                    abif_raw = tmpabif
                except:
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
            from setvar import set_dye_array
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
        from pybaselines.morphological import jbcd
        from numpy import around, multiply, array
        from scipy.signal import find_peaks
        global winwidth, peakpositions, peakheights, peakfwhms, peakchannels
        global peakareas, peaksizes, ch, x_plot, size_std
        h = self.getheight.value()
        w = self.getwidth.value()
        p = self.getprominence.value()
        winwidth = self.getwinwidth.value()
        peakpositions = []
        peakheights = []
        peakfwhms = []
        peakchannels = []
        peaksizes = []
        ch = []
        chN = []
        chP = []
        x_plot = list(dict(enumerate(abif_raw["DATA1"], start=1)))
        if should_sizecall:
            from scipy.interpolate import splrep, splev
            from numpy.polynomial.polynomial import Polynomial
            spline_degree = 0
            ILS_Name = self.ILS.currentText()
            Sizing_Method = self.SM.currentText()
            try:
                from setvar import set_ILS_channel
                ILSP = find_peaks(set_ILS_channel(abif_raw, ILS_Name),
                                  height=h, width=w, prominence=p,
                                  wlen=winwidth, rel_height=0.5)
                size_std = size_standards[ILS_Name]
                beginning_index = len(ILSP[0]) - len(size_std)
                ladder_peaks = ILSP[0][beginning_index:]
                if 'spline' in Sizing_Method:
                    from setvar import set_spl_dgr, set_knots
                    spline_degree = set_spl_dgr(Sizing_Method)
                    knots = set_knots(Sizing_Method, ladder_peaks,
                                      spline_degree)
                if spline_degree != 0:
                    spline = splrep(ladder_peaks, size_std, k=spline_degree,
                                    t=knots)
                    x_plot = around(splev(x_plot, spline), 3)
                elif 'order' in Sizing_Method:
                    from setvar import set_lsq_ord
                    func = Polynomial.fit(ladder_peaks, size_std,
                                          set_lsq_ord(Sizing_Method))
                    x_plot = around(func(array(x_plot)), 3)
            except ValueError:
                _sizingerror()
            except TypeError:
                _sizingerror()
            except KeyError:
                _sizingerror()
        # By default, find_peaks function measures width at
        # half maximum of height (rel_height=0.5). But
        # explicit is always better, then implicit, so
        # rel_height is specified clearly.
        for chnum in dyerange:
            if do_BCD:
                _, params = jbcd(abif_raw[udatac[chnum]], half_window=winwidth)
                ch.append(list(params['signal']))
            else:
                ch.append(list(abif_raw[udatac[chnum]]))
            chP.append(find_peaks(ch[chnum], height=h, width=w, prominence=p,
                                  wlen=winwidth, rel_height=0.5))
            chN.append([Dye[chnum]]*len(chP[chnum][0]))
            peakpositions += chP[chnum][0].tolist()
            peakheights += chP[chnum][1]['peak_heights'].tolist()
            peakfwhms += chP[chnum][1]['widths'].tolist()
            if should_sizecall and len(chP[chnum][0]) != 0:
                if spline_degree != 0:
                    peaksizes += list(splev(chP[chnum][0], spline))
                else:
                    peaksizes += list(func(chP[chnum][0]))
            peakchannels += chN[chnum]
# Well, we don't need all the digits after the point.
        peaksizes = around(peaksizes, 2)
        peakheights = around(peakheights, 2)
        peakfwhms = around(peakfwhms, 2)
        peakareas = around(multiply(peakheights, peakfwhms)*1.0645, 2)
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
        for i in dyerange:
            self.hidech[i].setText(ifacemsg['hidechannel'] + Dye[i])
        from setvar import set_graph_name
        self.graphWidget.setTitle(set_graph_name(abif_raw), color="c",
                                  size="10pt")
        max_x = len(x_plot)
        tmppen = ['b', 'g', 'y', 'r', 'orange', 'c', 'm', 'k']
        pen = []
        for i in dyerange:
            pen.append(tmppen[i])
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
        self.graphWidget.plotItem.setLimits(xMin=0, xMax=max_x, yMax=64000)
# Maximum peak height in files generated by new ABI 3500
# and SeqStudio family sequencers is 64000 arbitrary units.
        for i in dyerange:
            if show_channels[i]:
                self.graphWidget.plot(x_plot, ch[i], pen=pen[i])

    def export_csv(self):
        # Exporting CSV with data generated by findpeaks().
        expbox = self.sender()
        header = ['Peak Channel', 'Peak Position (Datapoints)', 'Peak Height',
                  'Peak FWHM', 'Peak Area (Datapoints)']
        do_export = False
        from setvar import chk_key_valid
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
                peak_chn.append(Dye[channel-1])
            pdarray = [peak_chn, abif_raw["Peak2"], abif_raw["Peak7"],
                       abif_raw["Peak5"], abif_raw["Peak10"],
                       abif_raw["Peak12"], abif_raw["Peak17"]]
            header += ['Peak Size (Bases)', 'Peak Area (Bases)']
            do_export = True
        else:
            msgbox(ifacemsg['unsuppeq'], ifacemsg['unsuppeqmsg'], 1)
        if do_export:
            from csv import writer
            from numpy import transpose
            peak_data = transpose(pdarray)
            csvname, _ = FileDialog.getSaveFileName(self, ifacemsg['savecsv'],
                                                    homedir, 'CSV(*.csv)')
            f = open(csvname, 'w', encoding='UTF8', newline='')
            writer(f).writerow(header)
            for row in peak_data:
                writer(f).writerow(row)
            f.close()

    def hide_ch(self):
        checkBox = self.sender()
        if checkBox.isChecked():
            show_channels[checkBox.number] = 0
        else:
            show_channels[checkBox.number] = 1
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
        from numpy import transpose
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
        global should_sizecall
        should_sizecall = False
        if self.sizecall.isChecked():
            should_sizecall = True
            self.sizecall.setChecked(False)
        self.retab()
        self.replot()
        should_sizecall = False
