# This file is part of FragalyseQt.
#
# FragalyseQt is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published 
# by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# FragalyseQt is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License along with FragalyseQt. If not, see <https://www.gnu.org/licenses/>.
import boxes, localize, fillarray
from os.path import expanduser, dirname
#Using FileDialog and SpinBox from pyqtgraph to prevent some possible problems for macOS users and to allow more fine variable setting.
from pyqtgraph import PlotWidget, FileDialog, SpinBox, ComboBox, TableWidget
#Using widgets from pyqtgraph to make program independent from Qt for Python implementation.
from pyqtgraph.Qt.QtWidgets import QCheckBox
from sizestandards import size_standards
ftype = "ABI fragment analysis files (*.fsa *.hid)"
global show_channels, ifacemsg, do_BCD
do_BCD = False
ifacemsg = {}
localize.localizefq(ifacemsg)
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
        self.graphWidget.setLabel('left', 'Signal intensity, relative fluorescent units')
        self.fsatab = TableWidget(self.centralwidget, sortable = False)
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
            if i%2 == 0:
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
        self.SM.setItems(['Cubic spline sizing','Linear spline sizing','5th degree spline sizing', 'LSQ weighted linear spline sizing', 'LSQ weighted cubic spline sizing', 'LSQ weighted 5th degree spline sizing', 'LSQ 2nd order', 'LSQ 3rd order', 'LSQ 5th order'])
        self.SM.setStyleSheet(''' font-size: 10px; ''')
        self.sizecall = QPushButton(self.centralwidget)
        self.sizecall.setGeometry(940, 580, 84, 20)
        self.sizecall.setCheckable(True)
        self.sizecall.setText('SizeCall')
        self.sizecall.setStyleSheet(''' font-size: 10px; ''')
        self.sizecall.clicked.connect(self.reanalyse)
        self.inactivatechkboxes()
    def inactivatechkboxes(self):
#Checkboxes without designations or with designations of nonexistent channels would look weird, so let's inactivate them correctly.
        i = 0
        while i < 8:
            self.hidech[i].setText(ifacemsg['ch_inact_msg'])
            i += 1
    def open_and_plot(self):
        from Bio.SeqIO import read as fsaread
        from charset_normalizer import from_bytes
        openBtn = self.sender()
        if openBtn.isChecked():
            openBtn.setChecked(False)
            global homedir, DN, x, Dye, graph_name, pen, keysarray, updatachnls, abif_raw
            updatachnls = ["DATA1","DATA2","DATA3","DATA4","DATA105","DATA106","DATA107","DATA108"]
            wavelng = ["DyeW1","DyeW2","DyeW3","DyeW4","DyeW5","DyeW6","DyeW7","DyeW8"]
            dyen = ["DyeN1","DyeN2","DyeN3","DyeN4","DyeN5","DyeN6","DyeN7","DyeN8"]
            fname, _ = FileDialog.getOpenFileName(self, 'Open file for analysis', homedir, ftype)
            FAfile = open(fname, "rb")
            try:
                tmprecord = fsaread(FAfile, "abi")
            except AssertionError:
                class record():
                    annotations = {"abif_raw":{"DATA1":None,"DATA2":None,"DATA3":None,"DATA4":None,"Dye#1": None,"DyeN1":None,"DyeN2":None,"DyeN3":None,"DyeN4":None,"MODL1":None}}
                tmprecord = record()
#Preventing data corruption in a case if target file is corrupted.
            FAfile.close()
#Closing file to save memory and avoid unexpected things.
            if tmprecord.annotations["abif_raw"]["DATA1"] == None:
#Assuming what it may be HID file.
                try:
                    HIDfile = open(fname, "rb")
                    s = HIDfile.read()
                    tmpkeys = tmprecord.annotations["abif_raw"].keys()
                    if "MODL1" in tmpkeys and tmprecord.annotations["abif_raw"]["MODL1"] == None:
                        HIDfile.seek(s.find(b'\x4d\x4f\x44\x4c\x00\x00\x00\x01') + 12, 0)
                        namesize = int.from_bytes(HIDfile.read(4), 'big')
#If we have file generated by old ABI 310 software - it may have different format of records.
                        if namesize < 4:
                            HIDfile.seek(s.find(b'\x4d\x4f\x44\x4c\x00\x00\x00\x01') + 10, 0)
                            namesize = int.from_bytes(HIDfile.read(2), 'big')
                            HIDfile.seek(8, 1)
                        else:
                            HIDfile.seek(4, 1)
                        tmprecord.annotations["abif_raw"]["MODL1"] = HIDfile.read(namesize)
                    if "Peak1" in tmpkeys and tmprecord.annotations["abif_raw"]["Peak1"] == None:
                        pshorthexarray = [b'\x50\x65\x61\x6b\x00\x00\x00\x01',b'\x50\x65\x61\x6b\x00\x00\x00\x05']
                        pinthexarray = [b'\x50\x65\x61\x6b\x00\x00\x00\x02',b'\x50\x65\x61\x6b\x00\x00\x00\x03',b'\x50\x65\x61\x6b\x00\x00\x00\x04',b'\x50\x65\x61\x6b\x00\x00\x00\x07',
                                        b'\x50\x65\x61\x6b\x00\x00\x00\x08',b'\x50\x65\x61\x6b\x00\x00\x00\x09',b'\x50\x65\x61\x6b\x00\x00\x00\x0a']
                        pdoublehexarray = [b'\x50\x65\x61\x6b\x00\x00\x00\x06',b'\x50\x65\x61\x6b\x00\x00\x00\x0b',b'\x50\x65\x61\x6b\x00\x00\x00\x0c',b'\x50\x65\x61\x6b\x00\x00\x00\x0d',
                                        b'\x50\x65\x61\x6b\x00\x00\x00\x0e',b'\x50\x65\x61\x6b\x00\x00\x00\x0f',b'\x50\x65\x61\x6b\x00\x00\x00\x10',b'\x50\x65\x61\x6b\x00\x00\x00\x11',
                                        b'\x50\x65\x61\x6b\x00\x00\x00\x12',b'\x50\x65\x61\x6b\x00\x00\x00\x15']
                        pshortname = ["Peak1","Peak5"]
                        pintname = ["Peak2","Peak3","Peak4","Peak7","Peak8","Peak9","Peak10"]
                        pdoublename = ["Peak6","Peak11","Peak12","Peak13","Peak14","Peak15","Peak16","Peak17","Peak18","Peak21"]
                        HIDfile.seek(s.find(pshorthexarray[0]) + 12, 0)
                        peakarraylength = int.from_bytes(HIDfile.read(4), 'big')
                        fillarray.fill_num_array(tmprecord.annotations["abif_raw"], HIDfile, s, pshortname, pshorthexarray, peakarraylength, '>H')
                        fillarray.fill_num_array(tmprecord.annotations["abif_raw"], HIDfile, s, pintname, pinthexarray, peakarraylength, '>I')
                        fillarray.fill_num_array(tmprecord.annotations["abif_raw"], HIDfile, s, pdoublename, pdoublehexarray, peakarraylength, '>d')
                    if tmprecord.annotations["abif_raw"]["Dye#1"] == None:
                        HIDfile.seek(s.find(b'\x44\x79\x65\x23\x00\x00\x00\x01') + 20, 0)
                        tmprecord.annotations["abif_raw"]["Dye#1"] = int.from_bytes(HIDfile.read(2), 'big')
                        if tmprecord.annotations["abif_raw"]["Dye#1"] > 4:
                            tmprecord.annotations["abif_raw"]["DATA105"] = None
                        if tmprecord.annotations["abif_raw"]["Dye#1"] > 5:
                            tmprecord.annotations["abif_raw"]["DATA106"] = None
                        if tmprecord.annotations["abif_raw"]["Dye#1"] > 6:
                            tmprecord.annotations["abif_raw"]["DATA107"] = None
                        if tmprecord.annotations["abif_raw"]["Dye#1"] > 7:
                            tmprecord.annotations["abif_raw"]["DATA108"] = None
                    fillarray.fill_char_array(tmprecord.annotations["abif_raw"], HIDfile, s, updatachnls, dyen, wavelng)
                    HIDfile.close()
                    abif_raw = tmprecord.annotations["abif_raw"]
                except:
#If it is not HID file and we fail to obtain any data - we should tell about this.
                    boxes.msgbox(ifacemsg['dmgdfile'], ifacemsg['nodatamsg'], 2)
                    try:
                        record
                    except NameError:
                        self.open_and_plot()
            else:
                abif_raw = tmprecord.annotations["abif_raw"]
#We need raw data from ABIF file only, no need in entire data structure, created by BioPython's AbiIO. This way multiple brackets constructions are evaded.
            x = list(dict(enumerate(abif_raw["DATA1"], start=1)))
            keysarray = abif_raw.keys()
            homedir = dirname(fname)
            self.inactivatechkboxes()
            graph_name = size_standard = equipment = ""
            DN = abif_raw["Dye#1"]
            Dye = ['']*DN
#First four channels are always present.
            pen = ['b', 'g', 'y', 'r']
            if DN >= 5:
                pen.append('orange')
            if DN >= 6:
                pen.append('c')
            if DN >= 7:
                pen.append('m')
            if DN == 8:
                pen.append('k')
            if "DyeN1" not in keysarray or abif_raw["DyeN1"] == None:
                Dye = ["FAM", "VIC", "TAMRA", "ROX"]
                if DN >= 5:
                    Dye.append("LIZ")
                if DN >= 6:
                    Dye.append("SID")
                if DN >= 7:
                    Dye.append("Channel 7")
                if DN == 8:
                    Dye.append("Channel 8")
            else:
                iteration = 0
                while iteration < DN:
                    Dye[iteration] = str(abif_raw[dyen[iteration]], 'UTF-8')
#Checking if dye names are present, but no emission wavelengths or wavelengths equal to zero. Assuming if DyeW1 is present and nonzero, others are present and non-zero too.
                    if "DyeW1" in keysarray and abif_raw["DyeW1"] != (None or 0):
                        Dye[iteration] += " " + str(abif_raw[wavelng[iteration]]) + " nm"
                    iteration += 1
            i = 0
            while i < DN:
                self.hidech[i].setText(ifacemsg['hidechannel'] + Dye[i])
                i += 1
#Assuming no more than 8 dyes are met at once.
            if "StdF1" in keysarray and abif_raw["StdF1"]!=b'':
                size_standard = str(abif_raw["StdF1"], 'UTF-8') + " size standard"
            else:
                size_standard = "Unknown size standard"
            if ("DySN1" and "MODF1") not in keysarray and abif_raw["MODL1"] != b'310 ':
                equipment = "RapidHIT ID v1.X"
#RapidHIT ID v1.X *.FSA files lack DySN1 and MODF1 keys, because there are only one dye set and only one run module.
            elif ("RunN1" and "DySN1" and "HCFG3") in keysarray and abif_raw["DySN1"] != None and (b'\xd1\xca' in abif_raw["DySN1"] or b'.avt' in abif_raw["RunN1"]) and abif_raw["HCFG3"] == b'3130xl':
                equipment = "Nanophore-05"
            elif abif_raw["MODL1"] == b'3200':
                equipment = "SeqStudio"
            elif "HCFG3" not in keysarray and ("NLNE1" and "DyeW1") in keysarray and abif_raw["DyeW1"] == 0:
                equipment = "Superyears Honor "
                if "DATA108" in keysarray:
                    if abif_raw["NLNE1"] == 16:
                        equipment += "1816"
                    else:
                        equipment += "1824"
                else:
                    if abif_raw["NLNE1"] == 16:
                        equipment += "1616"
                    elif abif_raw["NLNE1"] == 24:
                        equipment += "1624"
                    else:
                        equipment += "1696"
            elif "HCFG3" in keysarray and abif_raw["HCFG3"] != None:
                equipment = str(abif_raw["HCFG3"], 'UTF-8')
            elif abif_raw["MODL1"] != None:
                equipment = str(abif_raw["MODL1"], 'UTF-8')
            else:
                equipment = "Unknown equipment"
            graph_name = size_standard + ", " + equipment
            if "SpNm1" in keysarray:
                graph_name = str(from_bytes(abif_raw["SpNm1"]).best()) + ", " + graph_name
            elif "CTNM1" in keysarray:
                graph_name = str(from_bytes(abif_raw["CTNM1"]).best()) + ", " + graph_name
            self.reanalyse()
    def about(self):
        boxes.msgbox(ifacemsg['aboutbtn'], ifacemsg['infoboxtxt'], 0)
        self.aboutInfo.setChecked(False)
    def findpeaks(self):
#Detecting peaks and calculating peaks data.
        from numpy.polynomial.polynomial import Polynomial
        from numpy import around, multiply, array
        from scipy.signal import find_peaks
        global winwidth, peakpositions, peakheights, peakfwhms, peakchannels, peakareas, peaksizes, ch, x_plot, ILS_Name
        dgr = 0
        h = self.getheight.value()
        w = self.getwidth.value()
        p = self.getprominence.value()
        winwidth = self.getwinwidth.value()
        peakpositions = []
        peakheights = []
        peakfwhms = []
        peakchannels = []
        peaksizes = []
        ch = [0]*DN
        chN = ['']*DN
        chP = [dict]*DN
        iterator = 0
        while iterator < DN:
            ch[iterator] = list(abif_raw[updatachnls[iterator]])
            iterator += 1
        if do_BCD == False:
            pass
        else:
            from pybaselines.morphological import jbcd
            iterator = 0
            while iterator < DN:
                _, params = jbcd(ch[iterator], half_window=winwidth)
                ch[iterator] = list(params['signal'])
                iterator += 1
        x_plot = x
        if should_sizecall == True:
            try:
                from scipy.interpolate import splrep, splev
                spline_degree = 0
                knots = None
                ILS_Name = self.ILS.currentText()
                Sizing_Method = self.SM.currentText()
                if ILS_Name.find('ROX') != -1 or ILS_Name.find('CXR') != -1:
                    ILSchannel = ch[3]
                elif ILS_Name.find('LIZ') != -1 or ILS_Name.find('CC5') != -1 or ILS_Name.find('WEN') != -1 or ILS_Name.find('BTO') != -1 or ILS_Name.find('GDZ') != -1:
                    ILSchannel = ch[4]
                elif ILS_Name.find('CC0') != -1:
                    ILSchannel = ch[7]
                ILSP = find_peaks(ILSchannel, height=h, width=w, prominence=p, wlen=winwidth, rel_height=0.5)
                beginning_index = len(ILSP[0]) - len(size_standards[ILS_Name])
                if Sizing_Method.find('5th degree') != -1 and Sizing_Method.find('LSQ') == -1:
                    spline_degree = 5
                elif Sizing_Method.find('Cubic') != -1 and Sizing_Method.find('LSQ') == -1:
                    spline_degree = 3
                elif Sizing_Method.find('Linear') != -1 and Sizing_Method.find('LSQ') == -1:
                    spline_degree = 1
                elif Sizing_Method.find('order') == -1:
                    s_len = len(ILSP[0])
                    k1 = beginning_index + s_len//2 - s_len//3
                    k2 = s_len//2 + s_len//3
                    if Sizing_Method.find('5') != -1:
                        spline_degree = 5
                        if len(ILSP[0])-len(ILSP[0][k1:k2]) > 4:
                            knots = ILSP[0][k1:k2]
                        else:
                            #Making it work with GS120LIZ ladder too.
                            knots = ILSP[0][k1+3:k2]
                    elif Sizing_Method.find('cubic') != -1:
                        spline_degree = 3
                        if len(ILSP[0])-len(ILSP[0][k1:k2]) > 4:
                            knots = ILSP[0][k1:k2]
                        else:
                            #Making it work with GS120LIZ ladder too.
                            knots = ILSP[0][k1+1:k2]
                    elif Sizing_Method.find('linear') != -1:
                        spline_degree = 1
                        knots = ILSP[0][k1:k2]
                else:
                    if Sizing_Method.find('2') != -1:
                         dgr = 2
                    elif Sizing_Method.find('3') != -1:
                         dgr = 3
                    elif Sizing_Method.find('5') != -1:
                         dgr = 5
                if spline_degree != 0:
                    spline = splrep(ILSP[0][beginning_index:], size_standards[ILS_Name], k=spline_degree, t=knots)
                else:
                    func = Polynomial.fit(ILSP[0][beginning_index:], size_standards[ILS_Name], dgr)
                if dgr == 0:
                    x_plot = list(around(splev(x, spline), 3))
                else:
                    x_plot = around(func(array(x)), 3)
            except:
                boxes.msgbox("", ifacemsg['wrongsizing'], 1)
                self.sizecall.setChecked(False)
                self.reanalyse()
#By default, find_peaks function measures width at half maximum of height (rel_height=0.5).
#But explicit is always better, then implicit, so rel_height is specified clearly.
        channumber = 0
        while channumber < DN:
            chP[channumber] = find_peaks(ch[channumber], height=h, width=w, prominence=p, wlen=winwidth, rel_height=0.5)
            chN[channumber] = [Dye[channumber]]*len(chP[channumber][0].tolist())
            peakpositions += chP[channumber][0].tolist()
            peakheights += chP[channumber][1]['peak_heights'].tolist()
            peakfwhms += chP[channumber][1]['widths'].tolist()
            if should_sizecall == True and len(chP[channumber][0]) != 0:
                if dgr == 0:
                    peaksizes += list(around(splev(chP[channumber][0], spline), 2))
                else:
                    peaksizes += list(around(func(array(chP[channumber][0])), 2))
            channumber += 1
        for channel in chN:
            peakchannels += list(channel)
#Well, we don't need all the digits after the point.
        peakheights = around(peakheights, 2)
        peakfwhms = around(peakfwhms, 2)
        peakareas = around(multiply(multiply(peakheights, peakfwhms), 1.0645), 2)
#Peaks areas are calculated using formula for Gaussian peaks area ( https://www.physicsforums.com/threads/area-under-gaussian-peak-by-easy-measurements.419285/ ):
#A = FWHM*H/(2sqrt(2ln(2))/sqrt(2pi)) = 1.0645*FWHM*H, where FWHM is Full Width at Half Maximum. Real area may differ if peak is non-Gaussian, but at least 
#majority of peaks are of Gaussian shape.
#If peaks are well separated, you can directly calculate peak area, but if your peaks are crowded (e.g. in TP-PCR), oversaturated or you have noisy data - you MUST 
#use baseline correction and denoising prior peak area calculation.
    def replot(self):
        self.graphWidget.clear()
        self.graphWidget.setTitle(graph_name, color="c", size="10pt")
        if should_sizecall == True and len(peaksizes) > 0:
            #In the most normal case if you have good overall CE data quality, the last member of x_plot array should be the biggest size of all.
            x_max = x_plot[len(x_plot)-1]
            self.graphWidget.setLabel('bottom', 'Size, bases')
            max_ladder = max(size_standards[ILS_Name])
            if max_ladder+200 < x_max:
                self.graphWidget.plotItem.setLimits(xMin=0, xMax=max_ladder+200, yMax=64000)
            #EXTREMELY weird situation, but it sometimes happens, e.g. for low CE data quality with a lot of noise or poorly injected ladder.
            elif max_ladder > x_max:
                self.graphWidget.plotItem.setLimits(xMin=0, xMax=max_ladder, yMax=64000)
            else:
                self.graphWidget.plotItem.setLimits(xMin=0, xMax=x_max, yMax=64000)
        else:
            self.graphWidget.setLabel('bottom', 'Size, data points')
            self.graphWidget.plotItem.setLimits(xMin=0, xMax=len(x_plot), yMax=64000)
#Maximum peak height in files generated by new ABI 3500 and SeqStudio family sequencers is 64000 arbitrary units.
        i = 0
        while i < DN:
            if show_channels[i]:
                self.graphWidget.plot(x_plot, ch[i], pen=pen[i])
            i += 1
    def export_csv(self):
#Exporting CSV with data generated by findpeaks().
        expbox = self.sender()
        header = ['Peak Channel', 'Peak Position (Datapoints)', 'Peak Height', 'Peak FWHM', 'Peak Area (Datapoints)']
        do_export = False
        from numpy import transpose
        if expbox.focusWidget().objectName() == "CSV":
            peak_data_array = [peakchannels, peakpositions, peakheights, peakfwhms, peakareas]
            if len(peaksizes) > 0:
                peak_data_array.append(peaksizes)
                header += ['Peak Size (Bases)']
            peak_data = transpose(peak_data_array)
            do_export = True
        elif expbox.focusWidget().objectName() == "IA":
#Exporting internal analysis data.
            if "Peak1" in keysarray and abif_raw["Peak1"] != None:
#Checking if file has internal analysis data, assuming if Peak1 field is present, other fields are too.
                peak_channel = list(abif_raw["Peak1"])
                i = 0
                while i < len(peak_channel):
                    peak_channel[i] = Dye[peak_channel[i]-1]
                    i += 1
                peak_data = transpose([peak_channel, list(abif_raw["Peak2"]), list(abif_raw["Peak7"]), list(abif_raw["Peak5"]),
                list(abif_raw["Peak10"]), list(abif_raw["Peak12"]), list(abif_raw["Peak17"])])
                header += ['Peak Size (Bases)', 'Peak Area (Bases)']
                do_export = True
            else:
                boxes.msgbox(ifacemsg['unsupportedeq'], ifacemsg['unsupportedeqmsg'], 1)
        if do_export == True:
            from csv import writer
            csvname, _ = FileDialog.getSaveFileName(self, ifacemsg['savecsv'], homedir, 'CSV(*.csv)')
            f = open(csvname, 'w', encoding='UTF8', newline ='')
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
        basic_data = [peakchannels, peakpositions, peakheights, peakfwhms, peakareas]
        if len(peaksizes) <= 0:
            basic_data.append(["NaN"]*len(peakchannels))
        else:
            basic_data.append(peaksizes)
        from numpy import transpose
        self.fsatab.setData(transpose(basic_data))
        self.fsatab.setHorizontalHeaderLabels(['Peak Channel', 'Peak Position\n(Datapoints)', 'Peak Height', 'Peak FWHM', 'Peak Area\n(Datapoints)', 'Peak Size'])
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
