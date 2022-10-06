# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
# You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>. 

import sys, os, pyqtgraph, Bio.SeqIO, chardet, csv
from collections import defaultdict
from PyQt5 import QtCore, QtGui, QtWidgets
ftype = "ABI fragment analysis files (*.fsa)"
global show_channels
show_channels = [1] * 8
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
        self.pushButton.clicked.connect(self.fsa_open_and_draw)
        self.pushButton_2 = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_2.setGeometry(QtCore.QRect(101, 0, 100, 32))
        self.pushButton_2.setObjectName("pushButton_2")
        self.pushButton_2.setCheckable(True)
        self.pushButton_2.clicked.connect(self.about)
        self.pushButton_3 = QtWidgets.QPushButton(self.centralwidget)
        self.pushButton_3.setGeometry(QtCore.QRect(201, 0, 100, 32))
        self.pushButton_3.setObjectName("pushButton_3")
        self.pushButton_3.setCheckable(True)
        self.pushButton_3.clicked.connect(self.export_csv)
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
        self.pushButton_3.setText(QtCore.QCoreApplication.translate("MainWindow", "Export CSV"))
        self.inactivatechkboxes()
    def fsa_open_and_draw(self, pressed):
        source = self.sender()
        if pressed:
            fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open FSA file for analysis', os.path.expanduser('~'), ftype)
            global record, DN, x, Dye, graph_name, pen1, pen2, pen3, pen4, pen5, pen6, pen7, pen8
            record = Bio.SeqIO.read(open(fname, "rb"), "abi")
            encoding = defaultdict(list)
            x = []
            Dye = ['']*8
            self.inactivatechkboxes()
            graph_name = ""
            DyeN1 = DyeN2 = DyeN3 = DyeN4 = DyeN5 = DyeN6 = DyeN7 = DyeN8 = size_standard = equipment = ""
            DyeN1=str(record.annotations["abif_raw"]["DyeN1"], 'UTF-8')
            DyeN2=str(record.annotations["abif_raw"]["DyeN2"], 'UTF-8')
            DyeN3=str(record.annotations["abif_raw"]["DyeN3"], 'UTF-8')
            DyeN4=str(record.annotations["abif_raw"]["DyeN4"], 'UTF-8')
            Dye[0] = DyeN1 + " " + str(record.annotations["abif_raw"]["DyeW1"]) + " nm"
            pen1 = pyqtgraph.mkPen(color = (0, 0, 255))
            self.checkBox.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN1 + " channel"))
            Dye[1] = DyeN2 + " " + str(record.annotations["abif_raw"]["DyeW2"]) + " nm"
            pen2 = pyqtgraph.mkPen(color = (0, 255, 0))
            self.checkBox_2.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN2 + " channel"))
            Dye[2] = DyeN3 + " " + str(record.annotations["abif_raw"]["DyeW3"]) + " nm"
            pen3 = pyqtgraph.mkPen(color = (255, 240, 0))
            self.checkBox_3.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN3 + " channel"))
            Dye[3] = DyeN4 + " " + str(record.annotations["abif_raw"]["DyeW4"]) + " nm"
            pen4 = pyqtgraph.mkPen(color = (255, 0, 0))
            self.checkBox_4.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN4 + " channel"))
            DN = 4
#Assuming no less than 4 dyes are present.
            DN = record.annotations["abif_raw"]["Dye#1"]
            if DN >= 5:
                DyeN5 = str(record.annotations["abif_raw"]["DyeN5"], 'UTF-8')
                self.checkBox_5.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN5 + " channel"))
                Dye[4] = DyeN5 + " " + str(record.annotations["abif_raw"]["DyeW5"]) + " nm"
                pen5 = pyqtgraph.mkPen(color = (255, 165, 0))
            if DN >= 6:
                DyeN6 = str(record.annotations["abif_raw"]["DyeN6"], 'UTF-8')
                self.checkBox_6.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN6 + " channel"))
                Dye[5] = DyeN6 + " " + str(record.annotations["abif_raw"]["DyeW6"]) + " nm"
                pen6 = pyqtgraph.mkPen(color = (0, 255, 255))
            if DN >= 7:
                DyeN7 = str(record.annotations["abif_raw"]["DyeN7"], 'UTF-8')
                self.checkBox_7.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN7 + " channel"))
                Dye[6] = DyeN7 + " " + str(record.annotations["abif_raw"]["DyeW7"]) + " nm"
                pen7 = pyqtgraph.mkPen(color = (255, 0, 255))
            if DN == 8:
                DyeN8 = str(record.annotations["abif_raw"]["DyeN8"], 'UTF-8')
                self.checkBox_8.setText(QtCore.QCoreApplication.translate("MainWindow", "Hide " + DyeN8 + " channel"))
                Dye[7] = DyeN8 + " " + str(record.annotations["abif_raw"]["DyeW8"]) + " nm"
                pen8 = pyqtgraph.mkPen(color = (0, 0, 0))
#Assuming no more than 8 dyes are met at once.
            x = list(dict(enumerate(record.annotations["abif_raw"]["DATA1"])))
#Assuming traces for different dyes have equal number of data points.
            size_standard = str(record.annotations["abif_raw"]["StdF1"], 'UTF-8') + " size standard"
            equipment = "ABI " + str(record.annotations["abif_raw"]["HCFG3"], 'UTF-8')
            encoding = chardet.detect(record.annotations["abif_raw"]["RunN1"])
            graph_name = str(record.annotations["abif_raw"]["RunN1"], encoding['encoding']) + ", " + size_standard + ", " + equipment
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
        infobox.setText("FragalyseQt version 0.1, codename \"Huntington\".\n\n" +
                        "This program supports assays with up to 8 different dyes used simultaneously, selective dye channels hiding" +
                        " and non-Latin run names, but currently it only supports CSV export for FSA files generated by ABI 3500 and" +
                        " SeqStudio series equipment.\n\nLicensed under GNU GPL version 3.\n\n" +
                        "If you wish to contact author for any reason - please write at dorif11@gmail.com")
        infobox.exec_()
    def export_csv(self):
        if record.annotations["abif_raw"]["HCFG2"] == b'35XX' or record.annotations["abif_raw"]["HCFG2"] == b'SEQSTUDIO':
#Using ABI3500 and SeqStudio abilities for primary data analysis.
            peak_channel = list(record.annotations["abif_raw"]["Peak1"])
            for i, n in enumerate(peak_channel):
                peak_channel[i] = Dye[n-1]
            peak_data = zip(
                peak_channel,
                list(record.annotations["abif_raw"]["Peak2"]),
                list(record.annotations["abif_raw"]["Peak7"]),
                list(record.annotations["abif_raw"]["Peak10"]),
                list(record.annotations["abif_raw"]["Peak12"]),
                list(record.annotations["abif_raw"]["Peak17"])
                )
            csvname, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save CSV', os.path.expanduser('~'), 'CSV(*.csv)')
            with open(csvname, 'w', encoding='UTF8', newline ='') as f:
                writer = csv.writer(f)
                writer.writerow(['Peak Channel', 'Peak Position in Datapoints', 'Peak Hight', 'Peak Area in Datapoints', 'Peak Position in Bases', 'Peak Area in Bases'])
                for row in peak_data:
                    writer.writerow(row)
        else:
            msgbox = QtWidgets.QMessageBox()
            msgbox.setIcon(QtWidgets.QMessageBox.Warning)
            msgbox.setWindowTitle("Unsupported equipment!")
            msgbox.setText("This feature is currently supported only for ABI 3500 and SeqStudio family sequencers!")
            msgbox.exec_()
    def hide_ch(self):
        checkBox = self.sender()
        if checkBox.isChecked():
            show_channels[checkBox.number] = 0
            self.replot()
        else:
            show_channels[checkBox.number] = 1
            self.replot()