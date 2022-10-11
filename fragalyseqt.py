# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

import sys, os, pyqtgraph, Bio.SeqIO, chardet
from collections import defaultdict
from PyQt5 import QtCore, QtGui, QtWidgets
ftype = "ABI fragment analysis files (*.fsa)"
global show_channels, homedir
show_channels = [1] * 8
homedir = os.path.expanduser('~')
class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1024, 640)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.pushButton = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton.setGeometry(QtCore.QRect(0, 0, 100, 32))
        self.pushButton.setObjectName("pushButton")
        self.pushButton.setCheckable(True)
        self.pushButton.clicked.connect(self.open_and_plot)
        self.pushButton_2 = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_2.setGeometry(QtCore.QRect(101, 0, 100, 32))
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_2.setCheckable(True)
        self.pushButton_2.clicked.connect(self.about)
        self.pushButton_3 = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_3.setGeometry(QtCore.QRect(201, 0, 100, 32))
        self.pushButton_3.setObjectName("pushButton_3")
        self.pushButton_3.setCheckable(True)
        self.pushButton_3.clicked.connect(self.export_ia)
        self.pushButton_4 = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_4.setGeometry(QtCore.QRect(301, 0, 100, 32))
        self.pushButton_4.setObjectName("pushButton_3")
        self.pushButton_4.setCheckable(True)
        self.pushButton_4.clicked.connect(self.export_csv)
        self.graphWidget = pyqtgraph.PlotWidget(self.centralwidget)
        self.graphWidget.setGeometry(QtCore.QRect(0, 30, 1024, 400))
        self.graphWidget.setObjectName("graphicsView")
        self.graphWidget.setBackground('w')
        self.graphWidget.showGrid(x=True, y=True)
        MainWindow.setCentralWidget(self.centralwidget)
        self.checkBox = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox.setGeometry(QtCore.QRect(20, 440, 150, 24))
        self.checkBox.setObjectName("checkBox")
        self.checkBox.number = 0
        self.checkBox.toggled.connect(self.hide_ch)
        self.checkBox_2 = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_2.setGeometry(QtCore.QRect(20, 480, 150, 24))
        self.checkBox_2.setObjectName("checkBox_2")
        self.checkBox_2.number = 1
        self.checkBox_2.toggled.connect(self.hide_ch)
        self.checkBox_3 = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_3.setGeometry(QtCore.QRect(20, 520, 150, 24))
        self.checkBox_3.setObjectName("checkBox_3")
        self.checkBox_3.number = 2
        self.checkBox_3.toggled.connect(self.hide_ch)
        self.checkBox_4 = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_4.setGeometry(QtCore.QRect(20, 560, 150, 24))
        self.checkBox_4.setObjectName("checkBox_4")
        self.checkBox_4.number = 3
        self.checkBox_4.toggled.connect(self.hide_ch)
        self.checkBox_5 = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_5.setGeometry(QtCore.QRect(220, 440, 150, 24))
        self.checkBox_5.setObjectName("checkBox_5")
        self.checkBox_5.number = 4
        self.checkBox_5.toggled.connect(self.hide_ch)
        self.checkBox_6 = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_6.setGeometry(QtCore.QRect(220, 480, 150, 24))
        self.checkBox_6.setObjectName("checkBox_6")
        self.checkBox_6.number = 5
        self.checkBox_6.toggled.connect(self.hide_ch)
        self.checkBox_7 = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_7.setGeometry(QtCore.QRect(220, 520, 150, 24))
        self.checkBox_7.setObjectName("checkBox_7")
        self.checkBox_7.number = 6
        self.checkBox_7.toggled.connect(self.hide_ch)
        self.checkBox_8 = QtWidgets.QCheckBox(self.centralwidget)
        self.checkBox_8.setGeometry(QtCore.QRect(220, 560, 150, 24))
        self.checkBox_8.setObjectName("checkBox_8")
        self.checkBox_8.number = 7
        self.checkBox_8.toggled.connect(self.hide_ch)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 30))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
    def inactivatechkboxes(self):
#Checkboxes without designations or with designations of nonexistent channels would look weird, so let's inactivate them correctly.
        self.checkBox.setText(QtCore.QCoreApplication.translate("MainWindow", "Inactive channel"))
        self.checkBox_2.setText(QtCore.QCoreApplication.translate("MainWindow", "Inactive channel"))
        self.checkBox_3.setText(QtCore.QCoreApplication.translate("MainWindow", "Inactive channel"))
        self.checkBox_4.setText(QtCore.QCoreApplication.translate("MainWindow", "Inactive channel"))
        self.checkBox_5.setText(QtCore.QCoreApplication.translate("MainWindow", "Inactive channel"))
        self.checkBox_6.setText(QtCore.QCoreApplication.translate("MainWindow", "Inactive channel"))
        self.checkBox_7.setText(QtCore.QCoreApplication.translate("MainWindow", "Inactive channel"))
        self.checkBox_8.setText(QtCore.QCoreApplication.translate("MainWindow", "Inactive channel"))
    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtCore.QCoreApplication.translate("MainWindow", "FragalyseQt"))
        self.pushButton.setText(QtCore.QCoreApplication.translate("MainWindow", "Open FSA file"))
        self.pushButton_2.setText(QtCore.QCoreApplication.translate("MainWindow", "About"))
        self.pushButton_3.setText(QtCore.QCoreApplication.translate("MainWindow", "Export IA"))
        self.pushButton_4.setText(QtCore.QCoreApplication.translate("MainWindow", "Export CSV"))
        self.inactivatechkboxes()
    def open_and_plot(self):
        openBtn = self.sender()
        if openBtn.isChecked():
            fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open FSA file for analysis', homedir, ftype)
            global record, DN, x, Dye, graph_name, pen1, pen2, pen3, pen4, pen5, pen6, pen7, pen8
            tmprecord = Bio.SeqIO.read(open(fname, "rb"), "abi")
#Preventing data corruption in a case if target file is corrupted.
            if tmprecord.annotations["abif_raw"].get("DATA1") == None or tmprecord.annotations["abif_raw"]["DATA1"] == None:
#Assuming what if no data in first channel - no data would be in other channels too.
                ferrbox = QtWidgets.QMessageBox()
                ferrbox.setIcon(QtWidgets.QMessageBox.Critical)
                ferrbox.setWindowTitle("Damaged file!")
                ferrbox.setText("There are no readable data in file!")
                ferrbox.exec_()
            else:
                record = tmprecord
            x = []
            Dye = ['']*8
            self.inactivatechkboxes()
            graph_name = ""
            size_standard = equipment = ""
            pen1 = pyqtgraph.mkPen(color = (0, 0, 255))
            pen2 = pyqtgraph.mkPen(color = (0, 255, 0))
            pen3 = pyqtgraph.mkPen(color = (255, 240, 0))
            pen4 = pyqtgraph.mkPen(color = (255, 0, 0))
            DN = 4
#Assuming no less than 4 dyes are present.
            DN = record.annotations["abif_raw"]["Dye#1"]
            if record.annotations["abif_raw"]["DySN1"] == b'J6' or record.annotations["abif_raw"]["DySN1"] == b'J6-T' or record.annotations["abif_raw"]["DySN1"] == b'D':
                Dye[0] = "6-FAM 522 nm"
                Dye[1] = "VIC 554 nm"
                Dye[2] = "NED 575 nm"
                Dye[3] = "PET 595 nm"
                Dye[4] = "LIZ 655 nm"
                Dye[5] = "SID 620 nm"
                self.checkBox.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide 6-FAM channel"))
                self.checkBox_2.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide VIC channel"))
                self.checkBox_3.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide NED channel"))
                self.checkBox_4.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide PET channel"))
                self.checkBox_5.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide LIZ channel"))
                self.checkBox_6.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide SID channel"))
                pen5 = pyqtgraph.mkPen(color = (255, 165, 0))
                pen6 = pyqtgraph.mkPen(color = (0, 255, 255))
            elif record.annotations["abif_raw"]["DySN1"] == b'G5':
                Dye[0] = "6-FAM 522 nm"
                Dye[1] = "VIC 554 nm"
                Dye[2] = "NED 575 nm"
                Dye[3] = "PET 595 nm"
                Dye[4] = "LIZ 655 nm"
                self.checkBox.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide 6-FAM channel"))
                self.checkBox_2.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide VIC channel"))
                self.checkBox_3.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide NED channel"))
                self.checkBox_4.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide PET channel"))
                self.checkBox_5.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide LIZ channel"))
                pen5 = pyqtgraph.mkPen(color = (255, 165, 0))
            elif record.annotations["abif_raw"]["DySN1"] == b'E5':
                Dye[0] = "dR110"
                Dye[1] = "dR6G"
                Dye[2] = "dTAMRA"
                Dye[3] = "dROX"
                Dye[4] = "LIZ"
                self.checkBox.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide dR110 channel"))
                self.checkBox_2.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide dR6G channel"))
                self.checkBox_3.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide dTAMRA channel"))
                self.checkBox_4.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide dROX channel"))
                self.checkBox_5.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide LIZ channel"))
                pen5 = pyqtgraph.mkPen(color = (255, 165, 0))
            elif record.annotations["abif_raw"]["DySN1"] == b'D' or record.annotations["abif_raw"]["DySN1"] == b'F':
                Dye[0] = "6-FAM 522 nm"
                Dye[1] = "VIC 554 nm"
                Dye[2] = "NED 575 nm"
                Dye[3] = "PET 595 nm"
                self.checkBox.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide 6-FAM channel"))
                self.checkBox_2.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide VIC channel"))
                self.checkBox_3.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide NED channel"))
                self.checkBox_4.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide PET channel"))
            else:
                if record.annotations["abif_raw"].get("DyeN1") == None or record.annotations["abif_raw"]["DyeN1"] == None:
#If dye names are not indicated... Well, it is absolutely sure, wavelengths are not indicated too.
                    Dye[0] = "FAM"
                    Dye[1] = "VIC"
                    Dye[2] = "TAMRA"
                    Dye[3] = "ROX"
                    self.checkBox.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide FAM channel"))
                    self.checkBox_2.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide VIC channel"))
                    self.checkBox_3.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide TAMRA channel"))
                    self.checkBox_4.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide ROX channel"))
                    if DN >= 5:
                        Dye[4] = "LIZ"
                        self.checkBox_5.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide LIZ channel"))
                        pen5 = pyqtgraph.mkPen(color = (255, 165, 0))
                    if DN >= 6:
                        Dye[5] = "SID"
                        self.checkBox_5.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide SID channel"))
                        pen6 = pyqtgraph.mkPen(color = (0, 255, 255))
                    if DN >= 7:
                        Dye[5] = "Channel 7"
                        self.checkBox_7.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide channel 7"))
                        pen7 = pyqtgraph.mkPen(color = (255, 0, 255))
                    if DN == 8:
                        Dye[5] = "Channel 8"
                        self.checkBox_8.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide channel 8"))
                        pen8 = pyqtgraph.mkPen(color = (0, 0, 0))
                else:
                    DyeN = ['']*8
                    DyeN[0] = str(record.annotations["abif_raw"]["DyeN1"], 'UTF-8')
                    DyeN[1] = str(record.annotations["abif_raw"]["DyeN2"], 'UTF-8')
                    DyeN[2] = str(record.annotations["abif_raw"]["DyeN3"], 'UTF-8')
                    DyeN[3] = str(record.annotations["abif_raw"]["DyeN4"], 'UTF-8')
                    if record.annotations["abif_raw"].get("DyeW1") == None or record.annotations["abif_raw"]["DyeW1"] == None:
#A little bit stranger situation, but also possible - dye names are present, but without emission wavelengths.
                        Dye[0] = DyeN[0]
                        self.checkBox.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN[0] + " channel"))
                        Dye[1] = DyeN[1]
                        self.checkBox_2.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN[1] + " channel"))
                        Dye[2] = DyeN[2]
                        self.checkBox_3.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN[2] + " channel"))
                        Dye[3] = DyeN[3]
                        self.checkBox_4.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN[3] + " channel"))
                        if DN >= 5:
                            DyeN[4] = str(record.annotations["abif_raw"]["DyeN5"], 'UTF-8')
                            self.checkBox_5.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN[4] + " channel"))
                            Dye[4] = DyeN[4]
                            pen5 = pyqtgraph.mkPen(color = (255, 165, 0))
                        if DN >= 6:
                            DyeN[5] = str(record.annotations["abif_raw"]["DyeN6"], 'UTF-8')
                            self.checkBox_6.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN[5] + " channel"))
                            Dye[5] = DyeN[5]
                            pen6 = pyqtgraph.mkPen(color = (0, 255, 255))
                        if DN >= 7:
                            DyeN[6] = str(record.annotations["abif_raw"]["DyeN7"], 'UTF-8')
                            self.checkBox_7.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN[6] + " channel"))
                            Dye[6] = DyeN[6]
                            pen7 = pyqtgraph.mkPen(color = (255, 0, 255))
                        if DN == 8:
                            DyeN[7] = str(record.annotations["abif_raw"]["DyeN8"], 'UTF-8')
                            self.checkBox_8.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN[7] + " channel"))
                            Dye[7] = DyeN[7]
                            pen8 = pyqtgraph.mkPen(color = (0, 0, 0))
                    else:
                        Dye[0] = DyeN[0] + " " + str(record.annotations["abif_raw"]["DyeW1"]) + " nm"
                        Dye[1] = DyeN[1] + " " + str(record.annotations["abif_raw"]["DyeW2"]) + " nm"
                        Dye[2] = DyeN[2] + " " + str(record.annotations["abif_raw"]["DyeW3"]) + " nm"
                        Dye[3] = DyeN[3] + " " + str(record.annotations["abif_raw"]["DyeW4"]) + " nm"
                        self.checkBox.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN[0] + " channel"))
                        self.checkBox_2.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN[1] + " channel"))
                        self.checkBox_3.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN[2] + " channel"))
                        self.checkBox_4.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN[3] + " channel"))
                        if DN >= 5:
                            DyeN[4] = str(record.annotations["abif_raw"]["DyeN5"], 'UTF-8')
                            self.checkBox_5.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN[4] + " channel"))
                            Dye[4] = DyeN[4] + " " + str(record.annotations["abif_raw"]["DyeW5"]) + " nm"
                            pen5 = pyqtgraph.mkPen(color = (255, 165, 0))
                        if DN >= 6:
                            DyeN6 = str(record.annotations["abif_raw"]["DyeN6"], 'UTF-8')
                            self.checkBox_6.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN[5] + " channel"))
                            Dye[5] = DyeN[5] + " " + str(record.annotations["abif_raw"]["DyeW6"]) + " nm"
                            pen6 = pyqtgraph.mkPen(color = (0, 255, 255))
                        if DN >= 7:
                            DyeN7 = str(record.annotations["abif_raw"]["DyeN7"], 'UTF-8')
                            self.checkBox_7.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN[6] + " channel"))
                            Dye[6] = DyeN[6] + " " + str(record.annotations["abif_raw"]["DyeW7"]) + " nm"
                            pen7 = pyqtgraph.mkPen(color = (255, 0, 255))
                        if DN == 8:
                            DyeN8 = str(record.annotations["abif_raw"]["DyeN8"], 'UTF-8')
                            self.checkBox_8.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN[7] + " channel"))
                            Dye[7] = DyeN[7] + " " + str(record.annotations["abif_raw"]["DyeW8"]) + " nm"
                            pen8 = pyqtgraph.mkPen(color = (0, 0, 0))
#Assuming no more than 8 dyes are met at once.
            x = list(dict(enumerate(record.annotations["abif_raw"]["DATA1"])))
#Assuming traces for different dyes have equal number of data points.
            if record.annotations["abif_raw"].get("StdF1") == None or record.annotations["abif_raw"]["StdF1"] == None:
#More checks for the God of Checks! Error capture for the Throne!
                size_standard = "Unknown size standard"
            else:
                size_standard = str(record.annotations["abif_raw"]["StdF1"], 'UTF-8') + " size standard"
            if record.annotations["abif_raw"].get("HCFG3") == None or record.annotations["abif_raw"]["HCFG3"] == None:
                equipment = "Unknown equipment"
            else:
                equipment = "ABI " + str(record.annotations["abif_raw"]["HCFG3"], 'UTF-8')
            graph_name = size_standard + ", " + equipment
            if record.annotations["abif_raw"].get("RunN1") != None and record.annotations["abif_raw"]["RunN1"] != None:
                encoding = defaultdict(list)
                encoding = chardet.detect(record.annotations["abif_raw"]["RunN1"])
                graph_name = str(record.annotations["abif_raw"]["RunN1"], encoding['encoding']) + ", " + graph_name
            self.replot()
    def replot(self):
        self.graphWidget.clear()
        self.graphWidget.setTitle(graph_name, color="b", size="16pt")
        self.graphWidget.addLegend()
        if show_channels[0] == 1:
            self.graphWidget.plot(x, record.annotations["abif_raw"]["DATA1"], name=Dye[0], pen=pen1)
        if show_channels[1] == 1:
            self.graphWidget.plot(x, record.annotations["abif_raw"]["DATA2"], name=Dye[1], pen=pen2)
        if show_channels[2] == 1:
            self.graphWidget.plot(x, record.annotations["abif_raw"]["DATA3"], name=Dye[2], pen=pen3)
        if show_channels[3] == 1:
            self.graphWidget.plot(x, record.annotations["abif_raw"]["DATA4"], name=Dye[3], pen=pen4)
        if show_channels[4] == 1 and DN >= 5:
            self.graphWidget.plot(x, record.annotations["abif_raw"]["DATA105"], name=Dye[4], pen=pen5)
        if show_channels[5] == 1 and DN >= 6:
            self.graphWidget.plot(x, record.annotations["abif_raw"]["DATA106"], name=Dye[5], pen=pen6)
        if show_channels[6] == 1 and DN >= 7:
            self.graphWidget.plot(x, record.annotations["abif_raw"]["DATA107"], name=Dye[6], pen=pen7)
        if show_channels[7] == 1 and DN == 8:
            self.graphWidget.plot(x, record.annotations["abif_raw"]["DATA108"], name=Dye[7], pen=pen8)
    def about(self):
        infobox = QtWidgets.QMessageBox()
        infobox.setIcon(QtWidgets.QMessageBox.Information)
        infobox.setWindowTitle("Program Info")
        infobox.setText("FragalyseQt version 0.2, codename \"Friedreich\".\n\nThis program version supports" +
                        " assays with up to 8 different dyes used simultaneously, selective channel hiding" +
                        ", non-Latin run names and can correctly handle damaged files, exporting peaks " +
                        "locations, areas, FWHM's and channel names in CSV for any *.FSA files and export " +
                        "in CSV internal analysis data for *.FSA files generated by ABI 3500 and SeqStudio " +
                        "series equipment.\n\nPeaks areas are calculated assuming they are Gaussian peaks." +
                        "\n\nLicensed under GNU GPL version 3.\n\n" +
                        "If you wish to contact author for any reason - please write at dorif11@gmail.com")
        infobox.exec_()
    def findpeaks(self):
#Detecting peaks and calculating peaks data.
        import numpy
        from scipy.signal import find_peaks
        ch1 = ch2 = ch3 = ch4 = ch5 = ch6 = ch7 = ch8 = []
        h, _ = QtWidgets.QInputDialog.getInt(self, 'Input Minimal Peak Height', 'Input minimal peak height:', value = 250)
        w, _ = QtWidgets.QInputDialog.getInt(self, 'Input Minimal Peak Width', 'Input minimal peak width:', value = 5)
        p, _ = QtWidgets.QInputDialog.getInt(self, 'Input Minimal Peak Prominence', 'Input minimal peak prominence:', value = 100)
        ch1 = list(record.annotations["abif_raw"]["DATA1"])
        ch2 = list(record.annotations["abif_raw"]["DATA2"])
        ch3 = list(record.annotations["abif_raw"]["DATA3"])
        ch4 = list(record.annotations["abif_raw"]["DATA4"])
        ch1data = find_peaks(ch1, height=h, width=w, rel_height=0.5, prominence=p)
        ch2data = find_peaks(ch2, height=h, width=w, rel_height=0.5, prominence=p)
        ch3data = find_peaks(ch3, height=h, width=w, rel_height=0.5, prominence=p)
        ch4data = find_peaks(ch4, height=h, width=w, rel_height=0.5, prominence=p)
        global peakpositions, peakprominences, peakheights, peakfwhms, peakchannels, peakareas, peakareaspro
        chN1 = [Dye[0]]*len(ch1data[0])
        chN2 = [Dye[1]]*len(ch2data[0])
        chN3 = [Dye[2]]*len(ch3data[0])
        chN4 = [Dye[3]]*len(ch4data[0])
        peakchannels = []
        peakchannels = list(chN1) + list(chN2) + list(chN3) + list(chN4)
        peakpositions = ch1data[0].tolist() + ch2data[0].tolist() + ch3data[0].tolist() + ch4data[0].tolist()
        peakheights = ch1data[1]['peak_heights'].tolist() + ch2data[1]['peak_heights'].tolist() + ch3data[1]['peak_heights'].tolist() + ch4data[1]['peak_heights'].tolist()
        peakprominences = ch1data[1]['prominences'].tolist() + ch2data[1]['prominences'].tolist() + ch3data[1]['prominences'].tolist() + ch4data[1]['prominences'].tolist()
        peakfwhms = ch1data[1]['widths'].tolist() + ch2data[1]['widths'].tolist() + ch3data[1]['widths'].tolist() + ch4data[1]['widths'].tolist()
        if DN>=5:
            ch5 = list(record.annotations["abif_raw"]["DATA105"])
            ch5data = find_peaks(ch5, height=h, width=w, rel_height=0.5, prominence=p)
            peakpositions = peakpositions + ch5data[0].tolist()
            peakheights = peakheights + ch5data[1]['peak_heights'].tolist()
            peakprominences = peakprominences + ch5data[1]['prominences'].tolist()
            peakfwhms = peakfwhms + ch5data[1]['widths'].tolist()
            chN5 = [Dye[4]]*len(ch5data[0])
            peakchannels = peakchannels + list(chN5)
        if DN>=6:
            ch6 = list(record.annotations["abif_raw"]["DATA106"])
            ch6data = find_peaks(ch6, height=h, width=w, rel_height=0.5, prominence=p)
            peakpositions = peakpositions + ch6_peaks[0].tolist()
            peakheights = peakheights + ch6data[1]['peak_heights'].tolist()
            peakprominences = peakprominences + ch6data[1]['prominences'].tolist()
            peakfwhms = peakfwhms + ch6data[1]['widths'].tolist()
            chN6 = [Dye[5]]*len(ch6data[0])
            peakchannels = peakchannels + list(chN6)
        if DN>=7:
            ch7 = list(record.annotations["abif_raw"]["DATA107"])
            ch7data = find_peaks(ch7, height=h, width=w, rel_height=0.5, prominence=p)
            peakpositions = peakpositions + ch7_peaks[0].tolist()
            peakheights = peakheights + ch7data[1]['peak_heights'].tolist()
            peakprominences = peakprominences + ch7data[1]['prominences'].tolist()
            peakfwhms = peakfwhms + ch7data[1]['widths'].tolist()
            chN7 = [Dye[6]]*len(ch7data[0])
            peakchannels = peakchannels + list(chN7)
        if DN==8:
            ch8 = list(record.annotations["abif_raw"]["DATA108"])
            ch8data = find_peaks(ch8, height=h, width=w, rel_height=0.5, prominence=p)
            peakpositions = peakpositions + ch8_peaks[0].tolist()
            peakheights = peakheights + ch8data[1]['peak_heights'].tolist()
            peakprominences = peakprominences + ch8data[1]['prominences'].tolist()
            peakfwhms = peakfwhms + ch8data[1]['widths'].tolist()
            chN7 = [Dye[7]]*len(ch8data[0])
            peakchannels = peakchannels + list(chN8)
#Well, we don't need all the digits after the point.
        for i, n in enumerate(peakfwhms):
                peakfwhms[i] = round(peakfwhms[i], 2)
#Calculating peaks areas using formula for Gaussian peaks: A = FWHM*H/(2sqrt(2ln(2))/sqrt(2*pi)) = FWHM*H/0.94.
#FWHM is Full Width at Half Maximum.
#https://www.physicsforums.com/threads/area-under-gaussian-peak-by-easy-measurements.419285/
#Real area may be different if peak is non-Gaussian, but at least majority of them are.
        peakareas = list(peakheights)
        for i, n in enumerate(peakheights):
                peakareas[i] = round((peakareas[i]*peakfwhms[i]/0.94), 2)
    def export_csv(self):
#Exporting CSV with data generated by findpeaks().
        import csv
        self.findpeaks()
        peak_data = zip(
                peakchannels,
                peakpositions,
                peakheights,
                peakprominences,
                peakfwhms,
                peakareas,
                )
        csvname, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save CSV', homedir, 'CSV(*.csv)')
        with open(csvname, 'w', encoding='UTF8', newline ='') as f:
            writer = csv.writer(f)
            writer.writerow(['Peak Channel', 'Peak Position in Datapoints', 'Peak Height', 'Peak Prominence', 'Peak FWHM',
                             'Peak Area in Datapoints'])
            for row in peak_data:
                writer.writerow(row)
    def export_ia(self):
#Exporting internal analysis data.
        import csv
        if (record.annotations["abif_raw"].get("HCFG2") != None and record.annotations["abif_raw"]["HCFG3"] != None) and (record.annotations["abif_raw"]["HCFG2"] == b'35XX' or record.annotations["abif_raw"]["HCFG2"] == b'SEQSTUDIO'):
#Checking if file has record about equipment used and it is ABI 3500 or SeqStudio.
            peak_channel = list(record.annotations["abif_raw"]["Peak1"])
            for i, n in enumerate(peak_channel):
                peak_channel[i] = Dye[n-1]
#Using ABI3500 and SeqStudio abilities for primary data analysis.
            ia_data = zip(
                peak_channel,
                list(record.annotations["abif_raw"]["Peak2"]),
                list(record.annotations["abif_raw"]["Peak5"]),
                list(record.annotations["abif_raw"]["Peak7"]),
                list(record.annotations["abif_raw"]["Peak10"]),
                list(record.annotations["abif_raw"]["Peak12"]),
                list(record.annotations["abif_raw"]["Peak17"])
                )
            ianame, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save CSV', homedir, 'CSV(*.csv)')
            with open(ianame, 'w', encoding='UTF8', newline ='') as f:
                writer = csv.writer(f)
                writer.writerow(['Peak Channel', 'Peak Position in Datapoints', 'Peak FWHM', 'Peak Height', 'Peak Area in Datapoints',
                                 'Peak Position in Bases', 'Peak Area in Bases'])
                for row in ia_data:
                    writer.writerow(row)
        else:
            msgbox = QtWidgets.QMessageBox()
            msgbox.setIcon(QtWidgets.QMessageBox.Warning)
            msgbox.setWindowTitle("Unsupported equipment!")
            msgbox.setText("Internal analysis data could be exported only from files generated by ABI 3500 and" + 
                           " SeqStudio family sequencers!")
            msgbox.exec_()
    def hide_ch(self):
        checkBox = self.sender()
        if checkBox.isChecked():
            show_channels[checkBox.number] = 0
            self.replot()
        else:
            show_channels[checkBox.number] = 1
            self.replot()
