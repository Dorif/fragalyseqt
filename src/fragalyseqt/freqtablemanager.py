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

from os import listdir
from os.path import isdir, join, splitext, basename
from pyqtgraph.Qt.QtWidgets import (QDialog, QHBoxLayout, QVBoxLayout, QLabel,
                                    QTableWidget, QInputDialog, QFileDialog,
                                    QHeaderView, QDialogButtonBox,
                                    QTableWidgetItem,)
from pyqtgraph.Qt.QtCore import Qt
from .freqdb import (load_freq_table, import_freq_csv, import_freq_fam,
                     import_freq_gm, _is_gm_format, save_freq_table)
from .boxes import msgbox


class FreqTableManagerDialog(QDialog):
    def __init__(self, freqtables_dir: str, iface: dict, parent=None):
        super().__init__(parent)
        self._ftdir = freqtables_dir
        self._iface = iface
        self._tables = []
        self.setWindowTitle(iface['freq_tables_menu'])
        self.setMinimumSize(800, 450)
        self._build_ui()
        self._load_tables()

    def _build_ui(self):
        root = QVBoxLayout(self)

        split = QHBoxLayout()

        self._list_tbl = QTableWidget(0, 4)
        self._list_tbl.setHorizontalHeaderLabels(['Name', 'Population',
                                                  'Markers', 'θ'])
        self._list_tbl.setSelectionBehavior(
            QTableWidget.SelectRows if hasattr(QTableWidget, 'SelectRows')
            else QTableWidget.SelectionBehavior.SelectRows)
        self._list_tbl.setEditTriggers(
            QTableWidget.NoEditTriggers if hasattr(QTableWidget,
                                                   'NoEditTriggers')
            else QTableWidget.EditTrigger.NoEditTriggers)
        self._list_tbl.verticalHeader().setVisible(False)
        hdr = self._list_tbl.horizontalHeader()
        try:
            hdr.setSectionResizeMode(0, QHeaderView.Stretch)
            hdr.setSectionResizeMode(1, QHeaderView.ResizeToContents)
            hdr.setSectionResizeMode(2, QHeaderView.ResizeToContents)
            hdr.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        except AttributeError:
            hdr.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
            hdr.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
            hdr.setSectionResizeMode(2, QHeaderView.ResizeMode.ResizeToContents)
            hdr.setSectionResizeMode(3, QHeaderView.ResizeMode.ResizeToContents)
        self._list_tbl.currentCellChanged.connect(
            lambda row, col, prow, pcol: self._on_select(row))

        right = QVBoxLayout()
        self._detail_label = QLabel('Select a table to view its markers.')
        self._detail_label.setWordWrap(True)
        self._detail_tbl = QTableWidget(0, 2)
        self._detail_tbl.setHorizontalHeaderLabels(['Marker', 'Alleles'])
        self._detail_tbl.setEditTriggers(
            QTableWidget.NoEditTriggers if hasattr(QTableWidget,
                                                   'NoEditTriggers')
            else QTableWidget.EditTrigger.NoEditTriggers)
        self._detail_tbl.verticalHeader().setVisible(False)
        dhdr = self._detail_tbl.horizontalHeader()
        QHV = QHeaderView
        try:
            dhdr.setSectionResizeMode(0, QHV.Stretch)
            dhdr.setSectionResizeMode(1, QHV.ResizeToContents)
        except AttributeError:
            dhdr.setSectionResizeMode(0, QHV.ResizeMode.Stretch)
            dhdr.setSectionResizeMode(1, QHV.ResizeMode.ResizeToContents)
        right.addWidget(self._detail_label)
        right.addWidget(self._detail_tbl)

        split.addWidget(self._list_tbl, 2)
        split.addLayout(right, 3)
        root.addLayout(split)

        try:
            btns = QDialogButtonBox(QDialogButtonBox.Close)
        except AttributeError:
            btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        import_btn = btns.addButton(
            self._iface['importfreqtable'],
            QDialogButtonBox.ActionRole if hasattr(QDialogButtonBox,
                                                   'ActionRole')
            else QDialogButtonBox.ButtonRole.ActionRole)
        import_btn.clicked.connect(self._import_table)
        btns.rejected.connect(self.reject)
        root.addWidget(btns)

    def _load_tables(self):
        self._tables.clear()
        if not isdir(self._ftdir):
            return
        for fname in sorted(listdir(self._ftdir)):
            if not fname.endswith('.json'):
                continue
            try:
                t = load_freq_table(join(self._ftdir, fname))
                self._tables.append(t)
            except Exception:
                pass

        self._list_tbl.setRowCount(len(self._tables))
        try:
            no_edit = Qt.ItemIsEditable
        except AttributeError:
            no_edit = Qt.ItemFlag.ItemIsEditable

        for row, t in enumerate(self._tables):
            for col, text in enumerate([t.name, t.population,
                                        str(len(t.markers)),
                                        f'{t.theta:.3f}',]):
                item = QTableWidgetItem(text)
                item.setFlags(item.flags() & ~no_edit)
                self._list_tbl.setItem(row, col, item)

    def _import_table(self):
        path, _ = QFileDialog.getOpenFileName(
            self, self._iface['importfreqtable'], '',
            'Frequency tables (*.csv *.tsv *.fam *.txt);;'
            'CSV / TSV (*.csv *.tsv);;'
            'Familias (*.fam);;'
            'GeneMarker (*.txt);;'
            'All files (*)')
        if not path:
            return
        default_name = splitext(basename(path))[0]
        dlg = QInputDialog(self)
        dlg.setWindowTitle(self._iface['importfreqtable'])
        dlg.setLabelText('Table name:')
        dlg.setTextValue(default_name)
        if not dlg.exec():
            return
        name = dlg.textValue().strip()
        if not name:
            return
        try:
            ext = splitext(path)[1].lower()
            if ext == '.fam':
                t = import_freq_fam(path, name.strip())
            elif _is_gm_format(path):
                t = import_freq_gm(path, name.strip())
            else:
                t = import_freq_csv(path, name.strip(), '', '')
            dest = join(self._ftdir, name.strip().replace(' ', '_') + '.json')
            save_freq_table(t, dest)
            self._load_tables()
            msgbox(self._iface['importfreqtable'],
                   f'{t.name}: {len(t.markers)} markers imported.', 0)
        except Exception as exc:
            msgbox(self._iface['importfreqtable'], str(exc), 2)

    def _on_select(self, row: int):
        if row < 0 or row >= len(self._tables):
            return
        t = self._tables[row]
        self._detail_label.setText(
            f'{t.name}  ·  {t.population}  ·  θ={t.theta}  ·  min_freq={t.min_freq}')
        markers = sorted(t.markers.items())
        self._detail_tbl.setRowCount(len(markers))
        try:
            no_edit = Qt.ItemIsEditable
        except AttributeError:
            no_edit = Qt.ItemFlag.ItemIsEditable
        for i, (marker, alleles) in enumerate(markers):
            for col, text in enumerate([marker, str(len(alleles))]):
                item = QTableWidgetItem(text)
                item.setFlags(item.flags() & ~no_edit)
                self._detail_tbl.setItem(i, col, item)
