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

from pyqtgraph.Qt.QtWidgets import (
    QDialog, QHBoxLayout, QVBoxLayout,
    QTableWidget, QTableWidgetItem, QHeaderView, QLabel,
    QDialogButtonBox,
)
from pyqtgraph.Qt.QtCore import Qt

from .refprofile import list_profiles, get_profile


class RefProfileManagerDialog(QDialog):
    def __init__(self, db, iface: dict, parent=None):
        super().__init__(parent)
        self._db = db
        self._profiles = []
        self.setWindowTitle(iface['ref_profiles_menu'])
        self.setMinimumSize(800, 450)
        self._build_ui()
        self._load_profiles()

    def _build_ui(self):
        root = QVBoxLayout(self)

        split = QHBoxLayout()

        self._list_tbl = QTableWidget(0, 4)
        self._list_tbl.setHorizontalHeaderLabels(
            ['Name', 'Role', 'Created', 'Loci'])
        self._list_tbl.setSelectionBehavior(
            QTableWidget.SelectRows if hasattr(QTableWidget, 'SelectRows')
            else QTableWidget.SelectionBehavior.SelectRows)
        self._list_tbl.setEditTriggers(
            QTableWidget.NoEditTriggers if hasattr(QTableWidget, 'NoEditTriggers')
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
        self._detail_label = QLabel('Select a profile to view its alleles.')
        self._detail_tbl = QTableWidget(0, 3)
        self._detail_tbl.setHorizontalHeaderLabels(['Marker', 'Allele 1', 'Allele 2'])
        self._detail_tbl.setEditTriggers(
            QTableWidget.NoEditTriggers if hasattr(QTableWidget, 'NoEditTriggers')
            else QTableWidget.EditTrigger.NoEditTriggers)
        self._detail_tbl.verticalHeader().setVisible(False)
        dhdr = self._detail_tbl.horizontalHeader()
        try:
            dhdr.setSectionResizeMode(QHeaderView.Stretch)
        except AttributeError:
            dhdr.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        right.addWidget(self._detail_label)
        right.addWidget(self._detail_tbl)

        split.addWidget(self._list_tbl, 2)
        split.addLayout(right, 3)
        root.addLayout(split)

        try:
            btns = QDialogButtonBox(QDialogButtonBox.Close)
        except AttributeError:
            btns = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        btns.rejected.connect(self.reject)
        root.addWidget(btns)

    def _load_profiles(self):
        self._profiles = list_profiles(self._db)
        self._list_tbl.setRowCount(len(self._profiles))
        try:
            no_edit = Qt.ItemIsEditable
        except AttributeError:
            no_edit = Qt.ItemFlag.ItemIsEditable

        for row, p in enumerate(self._profiles):
            created = (p.get('created_at') or '')[:10]
            alleles = self._db.get_reference_alleles(p['id'])
            for col, text in enumerate([
                p.get('name', ''),
                p.get('role') or '',
                created,
                str(len(alleles)),
            ]):
                item = QTableWidgetItem(text)
                item.setFlags(item.flags() & ~no_edit)
                self._list_tbl.setItem(row, col, item)

    def _on_select(self, row: int):
        if row < 0 or row >= len(self._profiles):
            return
        p_info = self._profiles[row]
        try:
            profile = get_profile(self._db, p_info['id'])
        except Exception:
            return

        label_parts = [profile.name]
        if profile.role:
            label_parts.append(profile.role)
        if profile.notes:
            label_parts.append(profile.notes)
        self._detail_label.setText('  ·  '.join(label_parts))

        calls = sorted(profile.calls, key=lambda c: c.marker)
        self._detail_tbl.setRowCount(len(calls))
        try:
            no_edit = Qt.ItemIsEditable
        except AttributeError:
            no_edit = Qt.ItemFlag.ItemIsEditable

        for i, c in enumerate(calls):
            for col, text in enumerate([
                c.marker,
                c.allele1,
                c.allele2 or '—',
            ]):
                item = QTableWidgetItem(text)
                item.setFlags(item.flags() & ~no_edit)
                self._detail_tbl.setItem(i, col, item)
