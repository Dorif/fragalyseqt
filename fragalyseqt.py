# This file is part of FragalyseQt.
#
# FragalyseQt is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# FragalyseQt is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with FragalyseQt. If not, see <https://www.gnu.org/licenses/>.
import boxes, localize, fillarray
from os import path
#Using FileDialog and SpinBox from pyqtgraph to prevent some possible problems for macOS users and to allow more fine variable setting.
from pyqtgraph import PlotWidget, FileDialog, SpinBox
#Using widgets from pyqtgraph to make program independent from Qt for Python implementation.
from pyqtgraph.Qt.QtWidgets import QWidget, QCheckBox, QTableWidget, QTableWidgetItem
ftype = "ABI fragment analysis files (*.fsa *.hid)"
global show_channels, ifacemsg
do_BCD = False
ifacemsg = {}
localize.localizefq(ifacemsg)
show_channels = [1] * 8
homedir = path.expanduser('~')
class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        from pyqtgraph.Qt.QtWidgets import QPushButton, QLabel
        from pyqtgraph.Qt.QtGui import QIcon
        MainWindow.setWindowTitle("FragalyseQt")
        MainWindow.setWindowIcon(QIcon('FragalyseQt.png'))
        MainWindow.resize(1280, 720)
        self.centralwidget = QWidget(MainWindow)
        MainWindow.setCentralWidget(self.centralwidget)
        self.openFSA = QPushButton(self.centralwidget)
        self.openFSA.setGeometry(0, 0, 120, 20)
        self.openFSA.setCheckable(True)
        self.openFSA.setText(ifacemsg['openfragmentfile'])
        self.openFSA.clicked.connect(self.open_and_plot)
        self.aboutInfo = QPushButton(self.centralwidget)
        self.aboutInfo.setGeometry(121, 0, 100, 20)
        self.aboutInfo.setCheckable(True)
        self.aboutInfo.setText(ifacemsg['aboutbtn'])
        self.aboutInfo.clicked.connect(self.about)
        self.exportInternalAnalysisData = QPushButton(self.centralwidget)
        self.exportInternalAnalysisData.setGeometry(221, 0, 260, 20)
        self.exportInternalAnalysisData.setText(ifacemsg['exportinternal'])
        self.exportInternalAnalysisData.setObjectName("IA")
        self.exportInternalAnalysisData.clicked.connect(self.export_csv)
        self.exportCSV = QPushButton(self.centralwidget)
        self.exportCSV.setGeometry(481, 0, 120, 20)
        self.exportCSV.setText(ifacemsg['csvexport'])
        self.exportCSV.setObjectName("CSV")
        self.exportCSV.clicked.connect(self.export_csv)
        self.graphWidget = PlotWidget(self.centralwidget)
        self.graphWidget.setGeometry(0, 21, 1280, 360)
        self.graphWidget.setBackground(None)
        self.graphWidget.showGrid(x=True, y=True)
        self.graphWidget.setLabel('left', 'Signal intensity, relative fluorescent units')
        self.graphWidget.setLabel('bottom', 'Size, data points')
        self.fsatab = QTableWidget(self.centralwidget)
        self.fsatab.setGeometry(0, 380, 920, 340)
        self.fsatab.setColumnCount(5)
        self.fsatab.setHorizontalHeaderLabels(['Peak Channel', 'Peak Position in Datapoints', 'Peak Height', 'Peak FWHM', 'Peak Area in Datapoints'])
        self.fsatab.resizeColumnsToContents()
        self.getheightlabel = QLabel(self)
        self.getheightlabel.setGeometry(921, 381, 280, 20)
        self.getheightlabel.setText(ifacemsg['minph'])
        self.getheight = SpinBox(self)
        self.getheight.setGeometry(1201, 381, 80, 20)
        self.getheight.setRange(1, 64000)
        self.getheight.setValue(175)
        self.getheight.setOpts(minStep=1, dec=True)
        self.getheight.valueChanged.connect(self.reanalyse)
        self.getwidthlabel = QLabel(self)
        self.getwidthlabel.setGeometry(921, 401, 280, 20)
        self.getwidthlabel.setText(ifacemsg['minpw'])
        self.getwidth = SpinBox(self)
        self.getwidth.setGeometry(1201, 401, 80, 20)
        self.getwidth.setRange(1, 16000)
        self.getwidth.setValue(2)
        self.getwidth.setOpts(dec=True)
        self.getwidth.valueChanged.connect(self.reanalyse)
        self.getprominencelabel = QLabel(self)
        self.getprominencelabel.setGeometry(921, 421, 280, 20)
        self.getprominencelabel.setText(ifacemsg['minpp'])
        self.getprominence = SpinBox(self)
        self.getprominence.setGeometry(1201, 421, 80, 20)
        self.getprominence.setRange(1, 64000)
        self.getprominence.setValue(175)
        self.getprominence.setOpts(minStep=1, dec=True)
        self.getprominence.valueChanged.connect(self.reanalyse)
        self.getwinwidthlabel = QLabel(self)
        self.getwinwidthlabel.setGeometry(921, 441, 280, 20)
        self.getwinwidthlabel.setText(ifacemsg['minww'])
        self.getwinwidth = SpinBox(self)
        self.getwinwidth.setGeometry(1201, 441, 80, 20)
        self.getwinwidth.setRange(1, 1000)
        self.getwinwidth.setValue(15)
        self.getwinwidth.setOpts(minStep=1, dec=True)
        self.getwinwidth.valueChanged.connect(self.reanalyse)
        self.hidech1 = QCheckBox(self.centralwidget)
        self.hidech1.setGeometry(921, 461, 360, 20)
        self.hidech1.number = 0
        self.hidech1.toggled.connect(self.hide_ch)
        self.hidech2 = QCheckBox(self.centralwidget)
        self.hidech2.setGeometry(921, 481, 360, 20)
        self.hidech2.number = 1
        self.hidech2.toggled.connect(self.hide_ch)
        self.hidech3 = QCheckBox(self.centralwidget)
        self.hidech3.setGeometry(921, 501, 360, 20)
        self.hidech3.number = 2
        self.hidech3.toggled.connect(self.hide_ch)
        self.hidech4 = QCheckBox(self.centralwidget)
        self.hidech4.setGeometry(921, 521, 360, 20)
        self.hidech4.number = 3
        self.hidech4.toggled.connect(self.hide_ch)
        self.hidech5 = QCheckBox(self.centralwidget)
        self.hidech5.setGeometry(921, 541, 360, 20)
        self.hidech5.number = 4
        self.hidech5.toggled.connect(self.hide_ch)
        self.hidech6 = QCheckBox(self.centralwidget)
        self.hidech6.setGeometry(921, 561, 360, 20)
        self.hidech6.number = 5
        self.hidech6.toggled.connect(self.hide_ch)
        self.hidech7 = QCheckBox(self.centralwidget)
        self.hidech7.setGeometry(921, 581, 360, 20)
        self.hidech7.number = 6
        self.hidech7.toggled.connect(self.hide_ch)
        self.hidech8 = QCheckBox(self.centralwidget)
        self.hidech8.setGeometry(921, 601, 360, 20)
        self.hidech8.number = 7
        self.hidech8.toggled.connect(self.hide_ch)
        self.bcd = QCheckBox(self.centralwidget)
        self.bcd.setGeometry(921, 621, 360, 20)
        self.bcd.setText(ifacemsg['bcd'])
        self.bcd.toggled.connect(self.setbcd)
        self.inactivatechkboxes()
    def inactivatechkboxes(self):
#Checkboxes without designations or with designations of nonexistent channels would look weird, so let's inactivate them correctly.
        self.hidech1.setText(ifacemsg['ch_inact_msg'])
        self.hidech2.setText(ifacemsg['ch_inact_msg'])
        self.hidech3.setText(ifacemsg['ch_inact_msg'])
        self.hidech4.setText(ifacemsg['ch_inact_msg'])
        self.hidech5.setText(ifacemsg['ch_inact_msg'])
        self.hidech6.setText(ifacemsg['ch_inact_msg'])
        self.hidech7.setText(ifacemsg['ch_inact_msg'])
        self.hidech8.setText(ifacemsg['ch_inact_msg'])
    def open_and_plot(self):
        from Bio.SeqIO import read as fsaread
        from charset_normalizer import from_bytes
        openBtn = self.sender()
        if openBtn.isChecked():
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
            x = list(dict(enumerate(abif_raw["DATA1"])))
            keysarray = abif_raw.keys()
            homedir = path.dirname(fname)
            Dye = ['']*8
            self.inactivatechkboxes()
            graph_name = size_standard = equipment = ""
            DN = abif_raw["Dye#1"]
            pen = [''] * DN
            pen[0] = 'b'
            pen[1] = 'g'
            pen[2] = 'y'
            pen[3] = 'r'
            if DN >= 5:
                pen[4] = 'orange'
            if DN >= 6:
                pen[5] = 'c'
            if DN >= 7:
                pen[6] = 'm'
            if DN == 8:
                pen[7] = 'k'
            if "DyeN1" not in keysarray or abif_raw["DyeN1"] == None:
                Dye[0] = "FAM"
                Dye[1] = "VIC"
                Dye[2] = "TAMRA"
                Dye[3] = "ROX"
                if DN >= 5:
                    Dye[4] = "LIZ"
                if DN >= 6:
                    Dye[5] = "SID"
                if DN >= 7:
                    Dye[5] = "Channel 7"
                if DN == 8:
                    Dye[5] = "Channel 8"
            else:
                iteration = 0
                while iteration < DN:
                    Dye[iteration] = str(abif_raw[dyen[iteration]], 'UTF-8')
#Checking if dye names are present, but no emission wavelengths or wavelengths equal to zero. Assuming if DyeW1 is present and nonzero, others are present and non-zero too.
                    if "DyeW1" in keysarray and abif_raw["DyeW1"] != (None or 0):
                        Dye[iteration] += " " + str(abif_raw[wavelng[iteration]]) + " nm"
                    iteration += 1
            self.hidech1.setText(ifacemsg['hidechannel'] + Dye[0])
            self.hidech2.setText(ifacemsg['hidechannel'] + Dye[1])
            self.hidech3.setText(ifacemsg['hidechannel'] + Dye[2])
            self.hidech4.setText(ifacemsg['hidechannel'] + Dye[3])
            if DN >= 5:
                self.hidech5.setText(ifacemsg['hidechannel'] + Dye[4])
            if DN >= 6:
                self.hidech6.setText(ifacemsg['hidechannel'] + Dye[5])
            if DN >= 7:
                self.hidech7.setText(ifacemsg['hidechannel'] + Dye[6])
            if DN == 8:
                self.hidech8.setText(ifacemsg['hidechannel'] + Dye[7])
#Assuming no more than 8 dyes are met at once.
            if "StdF1" in keysarray and abif_raw["StdF1"]!=b'':
                size_standard = str(abif_raw["StdF1"], 'UTF-8') + " size standard"
            else:
                size_standard = "Unknown size standard "
                if "DyeB1" in keysarray:
#Most usual channels for size standards are 4th (ROX) and 5th (LIZ), so let's begin from them.
                    if abif_raw["DyeB4"]==b'S':
                        size_standard += "at channel 4"
                    elif DN > 4 and abif_raw["DyeB5"]==b'S':
                        size_standard += "at channel 5"
#Less usual, but possible case - size standard at 3rd (TAMRA) channel.
                    elif abif_raw["DyeB3"]==b'S':
                        size_standard += "at channel 3"
#Exotic cases.
                    elif abif_raw["DyeB1"]==b'S':
                        size_standard += "at channel 1"
                    elif abif_raw["DyeB2"]==b'S':
                        size_standard += "at channel 2"
                    elif DN > 5 and abif_raw["DyeB6"]==b'S':
                        size_standard += "at channel 6"
                    elif DN > 6 and abif_raw["DyeB7"]==b'S':
                        size_standard += "at channel 7"
                    elif DN > 7 and abif_raw["DyeB8"]==b'S':
                        size_standard += "at channel 8"
#If file contains no info about size standard and channel used for it...
                    else:
                        size_standard += "at unknown channel"
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
    def findpeaks(self):
#Detecting peaks and calculating peaks data.
        from scipy.signal import find_peaks
        global winwidth
        h = self.getheight.value()
        w = self.getwidth.value()
        p = self.getprominence.value()
        winwidth = self.getwinwidth.value()
        global peakpositions, peakheights, peakfwhms, peakchannels, peakareas, ch
        ch = [0]*DN
        chN = ['']*DN
        iterator = 0
        while iterator < DN:
            ch[iterator] = list(abif_raw[updatachnls[iterator]])
            iterator += 1
        if do_BCD == False:
            ch1data = find_peaks(ch[0], height=h, width=w, prominence=p, wlen=winwidth, rel_height=0.5)
            ch2data = find_peaks(ch[1], height=h, width=w, prominence=p, wlen=winwidth, rel_height=0.5)
            ch3data = find_peaks(ch[2], height=h, width=w, prominence=p, wlen=winwidth, rel_height=0.5)
            ch4data = find_peaks(ch[3], height=h, width=w, prominence=p, wlen=winwidth, rel_height=0.5)
            if DN>=5:
                ch5data = find_peaks(ch[4], height=h, width=w, prominence=p, wlen=winwidth, rel_height=0.5)
            if DN>=6:
                ch6data = find_peaks(ch[4], height=h, width=w, prominence=p, wlen=winwidth, rel_height=0.5)
            if DN>=7:
                ch7data = find_peaks(ch[4], height=h, width=w, prominence=p, wlen=winwidth, rel_height=0.5)
            if DN==8:
                ch8data = find_peaks(ch[4], height=h, width=w, prominence=p, wlen=winwidth, rel_height=0.5)
        else:
            from pybaselines.morphological import jbcd
            _, params = jbcd(ch[0], half_window=winwidth)
            ch1data = find_peaks(params['signal'], height=h, width=w, prominence=p, wlen=winwidth, rel_height=0.5)
            _, params = jbcd(ch[1], half_window=winwidth)
            ch2data = find_peaks(params['signal'], height=h, width=w, prominence=p, wlen=winwidth, rel_height=0.5)
            _, params = jbcd(ch[2], half_window=winwidth)
            ch3data = find_peaks(params['signal'], height=h, width=w, prominence=p, wlen=winwidth, rel_height=0.5)
            _, params = jbcd(ch[3], half_window=winwidth)
            ch4data = find_peaks(params['signal'], height=h, width=w, prominence=p, wlen=winwidth, rel_height=0.5)
            if DN>=5:
                _, params = jbcd(ch[4], half_window=winwidth)
                ch5data = find_peaks(params['signal'], height=h, width=w, prominence=p, wlen=winwidth, rel_height=0.5)
            if DN>=6:
                _, params = jbcd(ch[5], half_window=winwidth)
                ch6data = find_peaks(params['signal'], height=h, width=w, prominence=p, wlen=winwidth, rel_height=0.5)
            if DN>=7:
                _, params = jbcd(ch[6], half_window=winwidth)
                ch7data = find_peaks(params['signal'], height=h, width=w, prominence=p, wlen=winwidth, rel_height=0.5)
            if DN==8:
                _, params = jbcd(ch[7], half_window=winwidth)
                ch8data = find_peaks(params['signal'], height=h, width=w, prominence=p, wlen=winwidth, rel_height=0.5)
        chN[0] = [Dye[0]]*len(ch1data[0])
        chN[1] = [Dye[1]]*len(ch2data[0])
        chN[2] = [Dye[2]]*len(ch3data[0])
        chN[3] = [Dye[3]]*len(ch4data[0])
        peakchannels = list(chN[0]) + list(chN[1]) + list(chN[2]) + list(chN[3])
        peakpositions = ch1data[0].tolist() + ch2data[0].tolist() + ch3data[0].tolist() + ch4data[0].tolist()
        peakheights = ch1data[1]['peak_heights'].tolist() + ch2data[1]['peak_heights'].tolist() + ch3data[1]['peak_heights'].tolist() + ch4data[1]['peak_heights'].tolist()
        peakfwhms = ch1data[1]['widths'].tolist() + ch2data[1]['widths'].tolist() + ch3data[1]['widths'].tolist() + ch4data[1]['widths'].tolist()
        if DN>=5:
            peakpositions += ch5data[0].tolist()
            peakheights += ch5data[1]['peak_heights'].tolist()
            peakfwhms += ch5data[1]['widths'].tolist()
            chN[4] = [Dye[4]]*len(ch5data[0])
            peakchannels += list(chN[4])
        if DN>=6:
            peakpositions += ch6data[0].tolist()
            peakheights += ch6data[1]['peak_heights'].tolist()
            peakfwhms += ch6data[1]['widths'].tolist()
            chN[5] = [Dye[5]]*len(ch6data[0])
            peakchannels += list(chN[5])
        if DN>=7:
            peakpositions += ch7data[0].tolist()
            peakheights += ch7data[1]['peak_heights'].tolist()
            peakfwhms += ch7data[1]['widths'].tolist()
            chN[6] = [Dye[6]]*len(ch7data[0])
            peakchannels += list(chN[6])
        if DN==8:
            peakpositions += ch8data[0].tolist()
            peakheights += ch8data[1]['peak_heights'].tolist()
            peakfwhms += ch8data[1]['widths'].tolist()
            chN[7] = [Dye[7]]*len(ch8data[0])
            peakchannels += list(chN[7])
#Well, we don't need all the digits after the point.
        peakareas = list(peakheights)
        i = 0
        while i < len(peakfwhms):
                peakfwhms[i] = round(peakfwhms[i], 2)
                peakareas[i] = peakareas[i]*peakfwhms[i]/0.94
                i += 1
#Calculating peaks areas using formula for Gaussian peaks: A = FWHM*H/(2sqrt(2ln(2))/sqrt(2*pi)) = FWHM*H/0.94.
#FWHM is Full Width at Half Maximum.
#https://www.physicsforums.com/threads/area-under-gaussian-peak-by-easy-measurements.419285/
#Real area may be different if peak is non-Gaussian, but at least majority of them are.
#If peaks are well separated, peak prominence roughly equals peak height and either of them may be used to calculate peak area.
#If peaks are crowded (e.g. in TP-PCR) - you MUST use baseline correction and denoising prior peak area calculation.
#By default, find_peaks function measures width at half maximum of height (rel_height=0.5).
#But explicit is either way better, then implicit, so rel_height is specified clearly.
    def replot(self):
        from pybaselines.morphological import jbcd
        self.graphWidget.clear()
        self.graphWidget.setTitle(graph_name, color="b", size="12pt")
        self.graphWidget.plotItem.setLimits(xMin=0, xMax=len(x), yMax=64000)
#Maximum peak height in files generated by new ABI 3500 and SeqStudio family sequencers is 64000 arbitrary units.
        i = 0
        while i < DN:
            if show_channels[i]:
                if do_BCD == False:
                    self.graphWidget.plot(x, ch[i], pen=pen[i])
                else:
                    _, params = jbcd(ch[i], half_window=winwidth)
                    self.graphWidget.plot(x, params['signal'], pen=pen[i])
            i += 1
    def export_csv(self):
#Exporting CSV with data generated by findpeaks().
        expbox = self.sender()
        header = ['Peak Channel', 'Peak Position in Datapoints', 'Peak Height', 'Peak FWHM', 'Peak Area in Datapoints']
        do_export = False
        if expbox.focusWidget().objectName() == "CSV":
            peak_data = zip(peakchannels, peakpositions, peakheights, peakfwhms, peakareas)
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
                peak_data = zip(peak_channel,
                    list(abif_raw["Peak2"]),
                    list(abif_raw["Peak7"]),
                    list(abif_raw["Peak5"]),
                    list(abif_raw["Peak10"]),
                    list(abif_raw["Peak12"]),
                    list(abif_raw["Peak17"]))
                header.extend(['Peak Position in Bases', 'Peak Area in Bases'])
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
        count = 0
        while count < rowcount:
            self.fsatab.setItem(count, 0, QTableWidgetItem(str(peakchannels[count])))
            self.fsatab.setItem(count, 1, QTableWidgetItem(str(peakpositions[count])))
            self.fsatab.setItem(count, 2, QTableWidgetItem(str(peakheights[count])))
            self.fsatab.setItem(count, 3, QTableWidgetItem(str(peakfwhms[count])))
            self.fsatab.setItem(count, 4, QTableWidgetItem(str(peakareas[count])))
            count += 1
    def setbcd(self):
        checkBox = self.sender()
        global do_BCD
        if checkBox.isChecked():
            do_BCD = True
        else:
            do_BCD = False
        self.reanalyse()
    def reanalyse(self):
        self.retab()
        self.replot()
