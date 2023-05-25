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
        global ch_inact_msg, aboutbtn, infoboxtxt, openfragmentfile, exportinternal, csvexport, unsupportedeq, unsupportedeqmsg, dmgdfile, nodatamsg, hidechannel, minph, minpw, minpp, savecsv
        if os.name == 'posix':
            langvar = os.getenv('LANG')
        else:
            try:
                import locale, ctypes
                langvar = locale.windows_locale[ ctypes.windll.kernel32.GetUserDefaultUILanguage() ]
            except:
                langvar = "en"
        if "en" in langvar:
            ch_inact_msg = "Inactive channel"
            aboutbtn = "About"
            infoboxtxt = "FragalyseQt version 0.2.3, codename \"Friedreich\".\n\nThis program version supports analysis of up to 8 different fluorescense channels simultaneously, selective channel hiding, non-Latin run names and can correctly handle damaged files, exporting peaks locations, areas, FWHM's and channel names in CSV for any *.FSA files and export in CSV internal analysis data for *.FSA files generated by ABI 3500 and SeqStudio series equipment.\n\nPeaks areas are calculated assuming they are Gaussian peaks.\n\nLicensed under GNU GPL version 3.\n\nIf you wish to contact author for any reason - please write at dorif11@gmail.com"
            openfragmentfile = "Open FSA file"
            exportinternal = "Export Internal Analysis Data"
            csvexport = "Export in CSV"
            unsupportedeq = "Unsupported equipment!"
            unsupportedeqmsg = "Internal analysis data could be exported only from files generated by ABI 3500 and SeqStudio family sequencers!"
            dmgdfile = "Damaged file!"
            nodatamsg = "There are no data in file!"
            hidechannel = "Hide "
            minph = "Select minimal peak height"
            minpw = "Select minimal peak width"
            minpp = "Select minimal peak prominence"
            savecsv = "Save CSV"
        elif "ru" in langvar:
            ch_inact_msg = "Неактивный канал"
            aboutbtn = "О программе"
            infoboxtxt = "FragalyseQt версия 0.2.3, кодовое имя \"Фридрейх\".\n\nЭта версия программы поддерживает одновременный анализ до 8 каналов флуоресценции, избирательное сокрытие каналов, не-латинские имена запусков, может правильно работать с повреждёнными файлами, поддерживает экспорт имён каналов, положений, площадей и FWHM пиков в формат CSV для любого файла *.FSA и экспорт данных внутреннего анализа данных для файлов *.FSA, полученных на инструментах ABI 3500 и SeqStudio.\n\nПлощади пиков рассчитываются исходя из того, что пики являются гауссовыми.\n\nЛицензировано под GNU GPL версии 3.\n\nЕсли по какой-либо причине вы хотите проконтактировать с автором - пишите, пожалуйста, на dorif11@gmail.com"
            openfragmentfile = "Открыть файл FSA"
            exportinternal = "Экспорт Данных Внутреннего Анализа"
            csvexport = "Экспорт в CSV"
            unsupportedeq = "Неподдерживаемое оборудование!"
            unsupportedeqmsg = "Данные внутреннего анализа могут быть экспортированы только из файлов, полученных на ABI 3500 или семействе SeqStudio!"
            dmgdfile = "Повреждённый файл!"
            nodatamsg = "В файле отсутствуют данные!"
            hidechannel = "Скрыть "
            minph = "Выберите минимальную высоту пика"
            minpw = "Выберите минимальную ширину пика"
            minpp = "Выберите минимальное возвышение пика"
            savecsv = "Сохранить CSV"
        elif "ua" in langvar:
            ch_inact_msg = "Неактивний канал"
            aboutbtn = "О програмі"
            infoboxtxt = "FragalyseQt версія 0.2.3, кодове ім'я \"Фридрейх\".\n\n Версія програми підтримує одночасний аналіз до 8 каналів флуоресценции, вибіркове приховування каналів, не-латинські імена запусків, може правильно працювати з пошкодженими файлами, підтримує експорт імен каналів, положень, площадей і FWHM піков у формат CSV для будь-якого файла *.FSA та експорт даних внутрішнього аналізу даних для файлів *.FSA, отриманих на інструментах ABI 3500 и SeqStudio.\n\nПлощаді піков розраховуються виходячи з того, що піки є гауссовими.\n\Ліцензовано під GNU GPL версії 3.\n\nЯкщо з будь-якої причини ви хочете проконтактувати з автором - пишіть, будь ласка, на dorif11@gmail.com"
            openfragmentfile = "Відкрити  файл FSA"
            exportinternal = "Експорт Даних Внутренього Аналіза"
            csvexport = "Експорт в CSV"
            unsupportedeq = "Непідтримуване оболаднення!"
            unsupportedeqmsg = "Дані внутренього аналіза можуть бути експортовані тільки із файлів, отриманих на ABI 3500 або сімействі SeqStudio!"
            dmgdfile = "Пошкоджений файл!"
            nodatamsg = "У файлі відсутні дані!"
            hidechannel = "Сховати "
            minph = "Виберіть мінімальну висоту піку"
            minpw = "Виберіть мінімальну ширину піку"
            minpp = "Виберіть мінімальне підвищення піку"
            savecsv = "Зберегти CSV"
        elif "ro" in langvar:
            ch_inact_msg = "Canal inactiv"
            aboutbtn = "Despre program"
            infoboxtxt = "FragalyseQt versiune 0.2.3, codename \"Friedreich\".\n\nAceastă versiune permite analiza simultană a pînă la 8 canale de fluorescență, ascunderea selectivă a canalelor, nume non-latine a runurilor, poate funcționa corect cu fișiere deteriorate, permite exportarea numelor de canale, locațiile, suprafețele și FWHMurile a vîrfurilor în formatul CSV, pentru orice fișierul de tip *.FSA și exportarea în formatul CSV a datelor analizei interne pentru fișiere de tip *.FSA generate de echipamentele de model ABI 3500 și SeqStudio.\n\nSuprafețele vîrfurilor, sunt calculate presupunînd că vîrfurile sunt gausiene.\n\nLicențiat sub GNU GPL versiunea 3.\n\nDacă doriți să contactați autorul, pentru orice întrebare - scrieți, vă rog, la adresa dorif11@gmail.com"
            openfragmentfile = "Deschideți fișierul FSA"
            exportinternal = "Exportul Datelor Analizei Interne"
            csvexport = "Exportul în CSV"
            unsupportedeq = "Echipament nesusținut!"
            unsupportedeqmsg = "Datele analizei interne pot fi exportate numai pentru fișierile generate de echipamentele de model ABI 3500 și SeqStudio!"
            dmgdfile = "Fișier deteriorat!"
            nodatamsg = "In fișier lipsesc datele!"
            hidechannel = "Ascunde "
            minph = "Selectaţi înălţimea minimală a vîrfului"
            minpw = "Selectaţi lăţimea minimală a vîrfului"
            minpp = "Selectaţi proeminenţa minimală a vîrfului"
            savecsv = "Salvează CSV"
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.openFSA = QtWidgets.QPushButton(self.centralwidget)
        self.openFSA.setGeometry(QtCore.QRect(0, 0, 160, 32))
        self.openFSA.setObjectName("pushButton")
        self.openFSA.setCheckable(True)
        self.openFSA.clicked.connect(self.open_and_plot)
        self.aboutInfo = QtWidgets.QPushButton(self.centralwidget)
        self.aboutInfo.setGeometry(QtCore.QRect(161, 0, 100, 32))
        self.aboutInfo.setObjectName("aboutInfo")
        self.aboutInfo.setCheckable(True)
        self.aboutInfo.clicked.connect(self.about)
        self.exportInternalAnalysisData = QtWidgets.QPushButton(self.centralwidget)
        self.exportInternalAnalysisData.setGeometry(QtCore.QRect(261, 0, 280, 32))
        self.exportInternalAnalysisData.setObjectName("IA")
        self.exportInternalAnalysisData.setCheckable(True)
        self.exportInternalAnalysisData.clicked.connect(self.export_csv)
        self.exportCSV = QtWidgets.QPushButton(self.centralwidget)
        self.exportCSV.setGeometry(QtCore.QRect(541, 0, 120, 32))
        self.exportCSV.setObjectName("CSV")
        self.exportCSV.setCheckable(True)
        self.exportCSV.clicked.connect(self.export_csv)
        self.graphWidget = pyqtgraph.PlotWidget(self.centralwidget)
        self.graphWidget.setGeometry(QtCore.QRect(0, 30, 1024, 400))
        self.graphWidget.setObjectName("graphicsView")
        self.graphWidget.setBackground('w')
        self.graphWidget.showGrid(x=True, y=True)
        MainWindow.setCentralWidget(self.centralwidget)
        self.hidech1 = QtWidgets.QCheckBox(self.centralwidget)
        self.hidech1.setGeometry(QtCore.QRect(20, 440, 200, 24))
        self.hidech1.setObjectName(ch_inact_msg)
        self.hidech1.number = 0
        self.hidech1.toggled.connect(self.hide_ch)
        self.hidech2 = QtWidgets.QCheckBox(self.centralwidget)
        self.hidech2.setGeometry(QtCore.QRect(20, 480, 200, 24))
        self.hidech2.setObjectName(ch_inact_msg)
        self.hidech2.number = 1
        self.hidech2.toggled.connect(self.hide_ch)
        self.hidech3 = QtWidgets.QCheckBox(self.centralwidget)
        self.hidech3.setGeometry(QtCore.QRect(20, 520, 200, 24))
        self.hidech3.setObjectName(ch_inact_msg)
        self.hidech3.number = 2
        self.hidech3.toggled.connect(self.hide_ch)
        self.hidech4 = QtWidgets.QCheckBox(self.centralwidget)
        self.hidech4.setGeometry(QtCore.QRect(20, 560, 200, 24))
        self.hidech4.setObjectName(ch_inact_msg)
        self.hidech4.number = 3
        self.hidech4.toggled.connect(self.hide_ch)
        self.hidech5 = QtWidgets.QCheckBox(self.centralwidget)
        self.hidech5.setGeometry(QtCore.QRect(270, 440, 200, 24))
        self.hidech5.setObjectName(ch_inact_msg)
        self.hidech5.number = 4
        self.hidech5.toggled.connect(self.hide_ch)
        self.hidech6 = QtWidgets.QCheckBox(self.centralwidget)
        self.hidech6.setGeometry(QtCore.QRect(270, 480, 200, 24))
        self.hidech6.setObjectName(ch_inact_msg)
        self.hidech6.number = 5
        self.hidech6.toggled.connect(self.hide_ch)
        self.hidech7 = QtWidgets.QCheckBox(self.centralwidget)
        self.hidech7.setGeometry(QtCore.QRect(270, 520, 200, 24))
        self.hidech7.setObjectName(ch_inact_msg)
        self.hidech7.number = 6
        self.hidech7.toggled.connect(self.hide_ch)
        self.hidech8 = QtWidgets.QCheckBox(self.centralwidget)
        self.hidech8.setGeometry(QtCore.QRect(270, 560, 200, 24))
        self.hidech8.setObjectName(ch_inact_msg)
        self.hidech8.number = 7
        self.hidech8.toggled.connect(self.hide_ch)
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
        self.hidech1.setText(QtCore.QCoreApplication.translate("MainWindow", ch_inact_msg))
        self.hidech2.setText(QtCore.QCoreApplication.translate("MainWindow", ch_inact_msg))
        self.hidech3.setText(QtCore.QCoreApplication.translate("MainWindow", ch_inact_msg))
        self.hidech4.setText(QtCore.QCoreApplication.translate("MainWindow", ch_inact_msg))
        self.hidech5.setText(QtCore.QCoreApplication.translate("MainWindow", ch_inact_msg))
        self.hidech6.setText(QtCore.QCoreApplication.translate("MainWindow", ch_inact_msg))
        self.hidech7.setText(QtCore.QCoreApplication.translate("MainWindow", ch_inact_msg))
        self.hidech8.setText(QtCore.QCoreApplication.translate("MainWindow", ch_inact_msg))
    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtCore.QCoreApplication.translate("MainWindow", "FragalyseQt"))
        self.openFSA.setText(QtCore.QCoreApplication.translate("MainWindow", openfragmentfile))
        self.aboutInfo.setText(QtCore.QCoreApplication.translate("MainWindow", aboutbtn))
        self.exportInternalAnalysisData.setText(QtCore.QCoreApplication.translate("MainWindow", exportinternal))
        self.exportCSV.setText(QtCore.QCoreApplication.translate("MainWindow", csvexport))
        self.inactivatechkboxes()
    def open_and_plot(self):
        openBtn = self.sender()
        if openBtn.isChecked():
            fname, _ = QtWidgets.QFileDialog.getOpenFileName(self, 'Open FSA file for analysis', homedir, ftype)
            global record, DN, x, Dye, graph_name, pen
            FAfile = open(fname, "rb")
            tmprecord = Bio.SeqIO.read(FAfile, "abi")
#Preventing data corruption in a case if target file is corrupted.
            FAfile.close()
#Closing file to save memory and avoid unexpected things.
            if tmprecord.annotations["abif_raw"].get("DATA1") == None or tmprecord.annotations["abif_raw"]["DATA1"] == None:
#Assuming what if no data in first channel - no data would be in other channels too.
                ferrbox = QtWidgets.QMessageBox()
                ferrbox.setIcon(QtWidgets.QMessageBox.Critical)
                ferrbox.setWindowTitle(dmgdfile)
                ferrbox.setText(nodatamsg)
                ferrbox.exec_()
            else:
                record = tmprecord
            x = []
            Dye = ['']*8
            self.inactivatechkboxes()
            graph_name = ""
            size_standard = equipment = ""
            DN = 4
#Assuming no less than 4 dyes are present.
            DN = record.annotations["abif_raw"]["Dye#1"]
            pen = [''] * DN
            pen[0] = pyqtgraph.mkPen(color = (0, 0, 255))
            pen[1] = pyqtgraph.mkPen(color = (0, 255, 0))
            pen[2] = pyqtgraph.mkPen(color = (255, 240, 0))
            pen[3] = pyqtgraph.mkPen(color = (255, 0, 0))
            if DN >= 5:
                pen[4] = pyqtgraph.mkPen(color = (255, 165, 0))
                y5 = record.annotations["abif_raw"]["DATA105"]
            if DN >= 6:
                pen[5] = pyqtgraph.mkPen(color = (0, 255, 255))
                y6 = record.annotations["abif_raw"]["DATA106"]
            if DN >= 7:
                pen[6] = pyqtgraph.mkPen(color = (255, 0, 255))
                y7 = record.annotations["abif_raw"]["DATA107"]
            if DN == 8:
                pen[7] = pyqtgraph.mkPen(color = (0, 0, 0))
                y8 = record.annotations["abif_raw"]["DATA108"]
            if record.annotations["abif_raw"].get("DyeN1") == None or record.annotations["abif_raw"]["DyeN1"] == None:
#If dye names are not indicated... Well, it is absolutely sure, wavelengths are not indicated too.
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
                iteration = 1
                while iteration <= DN:
                    DyeCH = "DyeN" + str(iteration)
                    Dye[iteration-1] = str(record.annotations["abif_raw"][DyeCH], 'UTF-8')
#A little strange, but possible situation - dye names are present, but without emission wavelengths or with wavelengths equal to zero.
                    DyeWL = "DyeW" + str(iteration)
                    if DyeWL in record.annotations["abif_raw"].keys():
                        Dye[iteration-1] += " " + str(record.annotations["abif_raw"][DyeWL]) + " nm"
                    iteration += 1
            self.hidech1.setText(QtCore.QCoreApplication.translate("MainWindow", hidechannel + Dye[0]))
            self.hidech2.setText(QtCore.QCoreApplication.translate("MainWindow", hidechannel + Dye[1]))
            self.hidech3.setText(QtCore.QCoreApplication.translate("MainWindow", hidechannel + Dye[2]))
            self.hidech4.setText(QtCore.QCoreApplication.translate("MainWindow", hidechannel + Dye[3]))
            if DN >= 5:
                self.hidech5.setText(QtCore.QCoreApplication.translate("MainWindow", hidechannel + Dye[4]))
            if DN >= 6:
                self.hidech6.setText(QtCore.QCoreApplication.translate("MainWindow", hidechannel + Dye[5]))
            if DN >= 7:
                self.hidech7.setText(QtCore.QCoreApplication.translate("MainWindow", hidechannel + Dye[6]))
            if DN == 8:
                self.hidech8.setText(QtCore.QCoreApplication.translate("MainWindow", hidechannel + Dye[7]))
#Assuming no more than 8 dyes are met at once.
            x = list(dict(enumerate(record.annotations["abif_raw"]["DATA1"])))
#Traces for different dyes have equal number of data points.
            if "StdF1" in record.annotations["abif_raw"].keys():
                size_standard = str(record.annotations["abif_raw"]["StdF1"], 'UTF-8') + " size standard"
            else:
                size_standard = "Unknown size marker"
            if ("DySN1" and "DyeW1") not in record.annotations["abif_raw"].keys() and "HCFG3" not in record.annotations["abif_raw"].keys():
                equipment = "RapidHIT ID v1.X"
#RapidHIT ID v1.X use very limited subset of ABIF file specification - legacy from IntegeneX times, when they were concurrents with ABI.
            elif (b'\xd1\xca' in record.annotations["abif_raw"]["DySN1"] or b'.avt' in record.annotations["abif_raw"]["RunN1"]) and record.annotations["abif_raw"]["HCFG3"] == b'3130xl':
                equipment = "Nanophore-05"
            elif "HCFG3" in record.annotations["abif_raw"].keys() and record.annotations["abif_raw"]["HCFG3"] != None:
                if record.annotations["abif_raw"]["HCFG3"] == b'3200':
                    equipment = "SeqStudio"
                else:
                    equipment = "ABI " + str(record.annotations["abif_raw"]["HCFG3"], 'UTF-8')
            else:
                equipment = "Unknown equipment"
            graph_name = size_standard + ", " + equipment
            if "RunN1" in record.annotations["abif_raw"].keys():
                encoding = defaultdict(list)
                encoding = chardet.detect(record.annotations["abif_raw"]["RunN1"])
                graph_name = str(record.annotations["abif_raw"]["RunN1"], encoding['encoding']) + ", " + graph_name
            else:
                pass
            self.replot()
    def replot(self):
        self.graphWidget.clear()
        self.graphWidget.setTitle(graph_name, color="b", size="12pt")
        self.graphWidget.addLegend()
        i = 1
        while i <= DN:
            chd = ""
            if i <= 4:
                chd = "DATA" + str(i)
            else:
                chd = "DATA10" + str(i)
            if show_channels[i-1]:
                self.graphWidget.plot(x, record.annotations["abif_raw"][chd], pen=pen[i-1])
            i += 1
    def about(self):
        infobox = QtWidgets.QMessageBox()
        infobox.setIcon(QtWidgets.QMessageBox.Information)
        infobox.setWindowTitle("Program Info")
        infobox.setText(infoboxtxt)
        infobox.exec_()
    def findpeaks(self):
#Detecting peaks and calculating peaks data.
        import numpy
        from scipy.signal import find_peaks
        h, _ = QtWidgets.QInputDialog.getInt(self, minph, minph + ':', value = 250)
        w, _ = QtWidgets.QInputDialog.getInt(self, minpw, minpw + ':', value = 5)
        p, _ = QtWidgets.QInputDialog.getInt(self, minpp, minpp + ':', value = 100)
        global peakpositions, peakprominences, peakheights, peakfwhms, peakchannels, peakareas, peakareaspro
        ch = [0]*DN
        chN = ['']*DN
        iterator = 0
        while iterator < DN:
            chd = ''
            if iterator < 4:
                chd = "DATA" + str(iterator + 1)
            else:
                chd = "DATA10" + str(iterator + 1)
            ch[iterator] = list(record.annotations["abif_raw"][chd])
            iterator += 1
        ch1data = find_peaks(ch[0], height=h, width=w, prominence=p)
        ch2data = find_peaks(ch[1], height=h, width=w, prominence=p)
        ch3data = find_peaks(ch[2], height=h, width=w, prominence=p)
        ch4data = find_peaks(ch[3], height=h, width=w, prominence=p)
        chN[0] = [Dye[0]]*len(ch1data[0])
        chN[1] = [Dye[1]]*len(ch2data[0])
        chN[2] = [Dye[2]]*len(ch3data[0])
        chN[3] = [Dye[3]]*len(ch4data[0])
        peakchannels = list(chN[0]) + list(chN[1]) + list(chN[2]) + list(chN[3])
        peakpositions = ch1data[0].tolist() + ch2data[0].tolist() + ch3data[0].tolist() + ch4data[0].tolist()
        peakprominences = ch1data[1]['prominences'].tolist() + ch2data[1]['prominences'].tolist() + ch3data[1]['prominences'].tolist() + ch4data[1]['prominences'].tolist()
        peakfwhms = ch1data[1]['widths'].tolist() + ch2data[1]['widths'].tolist() + ch3data[1]['widths'].tolist() + ch4data[1]['widths'].tolist()
        if DN>=5:
            ch5data = find_peaks(ch[4], height=h, width=w, prominence=p)
            peakpositions += ch5data[0].tolist()
            peakprominences += ch5data[1]['prominences'].tolist()
            peakfwhms += ch5data[1]['widths'].tolist()
            chN[4] = [Dye[4]]*len(ch5data[0])
            peakchannels += list(chN[4])
        if DN>=6:
            ch6data = find_peaks(ch[5], height=h, width=w, prominence=p)
            peakpositions += ch6data[0].tolist()
            peakprominences += ch6data[1]['prominences'].tolist()
            peakfwhms += ch6data[1]['widths'].tolist()
            chN[5] = [Dye[5]]*len(ch6data[0])
            peakchannels += list(chN[5])
        if DN>=7:
            ch7data = find_peaks(ch[6], height=h, width=w, prominence=p)
            peakpositions += ch7data[0].tolist()
            peakprominences += ch7data[1]['prominences'].tolist()
            peakfwhms += ch7data[1]['widths'].tolist()
            chN[6] = [Dye[6]]*len(ch7data[0])
            peakchannels += list(chN[6])
        if DN==8:
            ch8data = find_peaks(ch[7], height=h, width=w, prominence=p)
            peakpositions += ch8data[0].tolist()
            peakprominences += ch8data[1]['prominences'].tolist()
            peakfwhms += ch8data[1]['widths'].tolist()
            chN[7] = [Dye[7]]*len(ch8data[0])
            peakchannels += list(chN[7])
#Well, we don't need all the digits after the point.
        for i, n in enumerate(peakfwhms):
                peakfwhms[i] = round(peakfwhms[i], 2)
#Calculating peaks areas using formula for Gaussian peaks: A = FWHM*H/(2sqrt(2ln(2))/sqrt(2*pi)) = FWHM*H/0.94.
#FWHM is Full Width at Half Maximum.
#https://www.physicsforums.com/threads/area-under-gaussian-peak-by-easy-measurements.419285/
#Real area may be different if peak is non-Gaussian, but at least majority of them are.
#For real, peak prominence is used as a height value, because only that part of peak has meaning.
#By default, find_peaks function measures width at half maximum of prominence.
        peakareas = list(peakprominences)
        for i, n in enumerate(peakprominences):
                peakareas[i] = round((peakareas[i]*peakfwhms[i]/0.94), 2)
    def export_csv(self):
#Exporting CSV with data generated by findpeaks().
        expbox = self.sender()
        header = []
        export = False
        if expbox.focusWidget().objectName() == "CSV":
            self.findpeaks()
            peak_data = zip(
                    peakchannels,
                    peakpositions,
                    peakprominences,
                    peakfwhms,
                    peakareas,
                    )
            header = ['Peak Channel', 'Peak Position in Datapoints', 'Peak Height', 'Peak FWHM', 'Peak Area in Datapoints']
            do_export = True
        elif expbox.focusWidget().objectName() == "IA":
#Exporting internal analysis data.            
            if "Peak1" in record.annotations["abif_raw"].keys():
#Checking if file has internal analysis data, assuming if Peak1 field is present, other fields are too.
                peak_channel = list(record.annotations["abif_raw"]["Peak1"])
                for i, n in enumerate(peak_channel):
                    peak_channel[i] = Dye[n-1]
#Using ABI3500 and SeqStudio abilities for primary data analysis.
                peak_data = zip(
                    peak_channel,
                    list(record.annotations["abif_raw"]["Peak2"]),
                    list(record.annotations["abif_raw"]["Peak5"]),
                    list(record.annotations["abif_raw"]["Peak7"]),
                    list(record.annotations["abif_raw"]["Peak10"]),
                    list(record.annotations["abif_raw"]["Peak12"]),
                    list(record.annotations["abif_raw"]["Peak17"])
                    )
                header = ['Peak Channel', 'Peak Position in Datapoints', 'Peak FWHM', 'Peak Height', 'Peak Area in Datapoints',
                          'Peak Position in Bases', 'Peak Area in Bases']
                do_export = True
            else:
                msgbox = QtWidgets.QMessageBox()
                msgbox.setIcon(QtWidgets.QMessageBox.Warning)
                msgbox.setWindowTitle(unsupportedeq)
                msgbox.setText(unsupportedeqmsg)
                msgbox.exec_()
        if do_export == True:
            import csv
            csvname, _ = QtWidgets.QFileDialog.getSaveFileName(self, savecsv, homedir, 'CSV(*.csv)')
            f = open(csvname, 'w', encoding='UTF8', newline ='')
            writer = csv.writer(f)
            writer.writerow(header)
            for row in peak_data:
                writer.writerow(row)
            f.close()
    def hide_ch(self):
        checkBox = self.sender()
        if checkBox.isChecked():
            show_channels[checkBox.number] = 0
            self.replot()
        else:
            show_channels[checkBox.number] = 1
            self.replot()