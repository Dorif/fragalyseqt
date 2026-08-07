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

from math import isfinite, floor, log10
from os.path import isdir, join, splitext, basename
from os import listdir
from pyqtgraph import FileDialog
from pyqtgraph.Qt.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel,
                                    QGridLayout, QComboBox, QPushButton,
                                    QRadioButton, QDoubleSpinBox, QTableWidget,
                                    QTableWidgetItem, QHeaderView, QWidget,
                                    QFileDialog, QInputDialog,)
from pyqtgraph.Qt.QtGui import QColor
from pyqtgraph.Qt.QtCore import Qt
from .boxes import msgbox
from .forensicstats import RELATIONSHIPS
from .freqdb import load_freq_table, save_freq_table
from .comparison import (allele_calls_from_state, compare_identity,
                         compare_kinship, export_comparison_csv,)
from .pdfreport import export_comparison_pdf
from .refprofile import (list_profiles, get_profile, store_profile,
                         profiles_from_codis_xml)

_SUP_TABLE = str.maketrans('0123456789-', '⁰¹²³⁴⁵⁶⁷⁸⁹⁻')

# Frequency table preselected when the comparison dialog opens. Must match the
# "name" field of src/fragalyseqt/freqtables/STRidER_2025_Combined.json.
_DEFAULT_TABLE_NAME = 'STRidER 2025 — Entire Database'


def _fmt_lr(val: float) -> str:
    if not isfinite(val):
        return '∞' if val > 0 else '-∞'
    if val == 0:
        return '0'
    exp = int(floor(log10(abs(val))))
    if abs(exp) < 4:
        return f'{val:.4g}'
    mantissa = val / 10 ** exp
    sup = str(exp).translate(_SUP_TABLE)
    return f'{mantissa:.2f}×10{sup}'


class ComparisonDialog(QDialog):
    def __init__(self, file_states, tab_names, freqtables_dir, iface,
                 db=None, parent=None):
        super().__init__(parent)
        self._states = file_states
        self._tab_names = tab_names
        self._ftdir = freqtables_dir
        self._msg = iface
        self._db = db
        self._tables: dict = {}
        self._last_result = None
        self.setWindowTitle(iface['cmp_title'])
        self.setMinimumWidth(720)
        self._build_ui()
        self._refresh_tables()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setSpacing(8)

        grid = QGridLayout()
        grid.setVerticalSpacing(6)
        grid.setHorizontalSpacing(8)

        self._prof1_combo = QComboBox()
        self._prof2_combo = QComboBox()
        self._populate_profile_combo(self._prof1_combo, default_second=False)
        self._populate_profile_combo(self._prof2_combo, default_second=True)

        self._table_combo = QComboBox()

        ft_row = QWidget()
        ft_layout = QHBoxLayout(ft_row)
        ft_layout.setContentsMargins(0, 0, 0, 0)
        ft_layout.addWidget(self._table_combo, 1)
        if self._db is not None:
            import_codis_btn = QPushButton(self._msg['cmp_import_codis'])
            import_codis_btn.clicked.connect(self._import_codis)
            ft_layout.addWidget(import_codis_btn)

        self._theta_spin = QDoubleSpinBox()
        self._theta_spin.setRange(0.0, 0.05)
        self._theta_spin.setDecimals(3)
        self._theta_spin.setSingleStep(0.001)
        self._theta_spin.setValue(0.01)

        grid.addWidget(QLabel(self._msg['cmp_prof1']), 0, 0)
        grid.addWidget(self._prof1_combo, 0, 1)
        grid.addWidget(QLabel(self._msg['cmp_prof2']), 1, 0)
        grid.addWidget(self._prof2_combo, 1, 1)
        grid.addWidget(QLabel(self._msg['cmp_freqtable']), 2, 0)
        grid.addWidget(ft_row, 2, 1)
        grid.addWidget(QLabel(self._msg['cmp_theta']), 3, 0)
        grid.addWidget(self._theta_spin, 3, 1)
        layout.addLayout(grid)

        mode_row = QHBoxLayout()
        self._mode_id = QRadioButton(self._msg['cmp_identity'])
        self._mode_ki = QRadioButton(self._msg['cmp_kinship'])
        self._mode_id.setChecked(True)
        mode_row.addWidget(self._mode_id)
        mode_row.addWidget(self._mode_ki)
        mode_row.addStretch()
        layout.addLayout(mode_row)

        rel_row = QHBoxLayout()
        self._rel_label = QLabel(self._msg['cmp_relationship'])
        self._rel_combo = QComboBox()
        for name in RELATIONSHIPS:
            self._rel_combo.addItem(name)
        self._rel_combo.setCurrentText('Half Siblings')
        rel_row.addWidget(self._rel_label)
        rel_row.addWidget(self._rel_combo)
        rel_row.addStretch()
        layout.addLayout(rel_row)

        self._mode_id.toggled.connect(self._on_mode_changed)
        self._on_mode_changed()

        calc_btn = QPushButton(self._msg['cmp_calculate'])
        calc_btn.clicked.connect(self._calculate)
        layout.addWidget(calc_btn)

        self._results_widget = QWidget()
        res_layout = QVBoxLayout(self._results_widget)
        res_layout.setContentsMargins(0, 8, 0, 0)

        summary = QGridLayout()
        self._lr_val = QLabel()
        self._log10_val = QLabel()
        self._concl_val = QLabel()
        self._loci_val = QLabel()
        summary.addWidget(QLabel(self._msg['cmp_combined_lr']), 0, 0)
        summary.addWidget(self._lr_val, 0, 1)
        summary.addWidget(QLabel(self._msg['cmp_log10_lr']), 0, 2)
        summary.addWidget(self._log10_val, 0, 3)
        summary.addWidget(QLabel(self._msg['cmp_conclusion']), 1, 0)
        summary.addWidget(self._concl_val, 1, 1, 1, 3)
        summary.addWidget(self._loci_val, 2, 0, 1, 4)
        res_layout.addLayout(summary)

        self._tbl = QTableWidget()
        self._tbl.setColumnCount(7)
        self._tbl.setHorizontalHeaderLabels(['Marker', 'Profile 1',
                                             'Profile 2', 'p(a1)', 'p(a2)',
                                             'LR / KI', 'Note',])
        self._tbl.setMinimumHeight(180)
        hdr = self._tbl.horizontalHeader()
        try:
            hdr.setSectionResizeMode(QHeaderView.Stretch)
        except AttributeError:
            hdr.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self._tbl.verticalHeader().setVisible(False)
        res_layout.addWidget(self._tbl)

        btn_row = QHBoxLayout()
        export_btn = QPushButton(self._msg['cmp_export_csv'])
        export_btn.clicked.connect(self._export_csv)
        pdf_btn = QPushButton(self._msg['cmp_export_pdf'])
        pdf_btn.clicked.connect(self._export_pdf)
        close_btn = QPushButton('Close')
        close_btn.clicked.connect(self.reject)
        btn_row.addWidget(export_btn)
        btn_row.addWidget(pdf_btn)
        btn_row.addStretch()
        btn_row.addWidget(close_btn)
        res_layout.addLayout(btn_row)

        layout.addWidget(self._results_widget)
        self._results_widget.setVisible(False)

    def _populate_profile_combo(self, combo, default_second=False):
        combo.clear()
        for i, name in enumerate(self._tab_names):
            combo.addItem(name, ('tab', i))
        if self._db is not None:
            saved = list_profiles(self._db)
            if saved:
                # A separator only belongs BETWEEN two non-empty groups. With
                # no files open the saved profiles must start at index 0,
                # otherwise the combo opens on a blank separator row whose
                # itemData is None and no profile is selected at all.
                if combo.count():
                    combo.insertSeparator(combo.count())
                for p in saved:
                    label = f"[{p['role'] or '—'}] {p['name']}"
                    combo.addItem(label, ('profile', p['id']))
        # Preselect a real entry, skipping separators: the second combo takes
        # the second selectable item so the dialog never opens with one and
        # the same profile compared against itself.
        selectable = [i for i in range(combo.count())
                      if combo.itemData(i) is not None]
        if not selectable:
            return
        if default_second and len(selectable) > 1:
            combo.setCurrentIndex(selectable[1])
        else:
            combo.setCurrentIndex(selectable[0])

    def _get_calls_for_combo(self, combo):
        data = combo.itemData(combo.currentIndex())
        if data is None:
            return []
        source, ref = data
        if source == 'tab':
            return allele_calls_from_state(self._states[ref])
        return get_profile(self._db, ref).calls

    def _import_codis(self):
        path, _ = QFileDialog.getOpenFileName(
            self, self._msg['cmp_import_codis'], '',
            'XML (*.xml);;All files (*)')
        if not path:
            return
        try:
            profiles = profiles_from_codis_xml(path)
            if not profiles:
                msgbox(self._msg['cmp_title'], 'No specimens found in XML.', 1)
                return
            for p in profiles:
                store_profile(self._db, p)
            msgbox(self._msg['cmp_title'],
                   self._msg['cmp_codis_imported'].format(n=len(profiles)), 0)
            self._populate_profile_combo(self._prof1_combo,
                                         default_second=False)
            self._populate_profile_combo(self._prof2_combo, default_second=True)
        except Exception as exc:
            msgbox(self._msg['cmp_title'], str(exc), 2)

    def _on_mode_changed(self):
        ki = self._mode_ki.isChecked()
        self._rel_label.setEnabled(ki)
        self._rel_combo.setEnabled(ki)

    def _refresh_tables(self):
        self._table_combo.clear()
        self._tables.clear()
        if not isdir(self._ftdir):
            return
        for fname in sorted(listdir(self._ftdir)):
            if not fname.endswith('.json'):
                continue
            try:
                t = load_freq_table(join(self._ftdir, fname))
                self._tables[t.name] = t
                self._table_combo.addItem(t.name)
            except Exception:
                pass
        self._select_default_table()

    def _select_default_table(self):
        # Preselect the STRidER 2025 whole-database table: it is the broadest
        # reference set shipped with FragalyseQt and the sensible starting
        # point for casework, whereas plain alphabetical order would land on
        # whichever population happens to sort first. Matched on the table's
        # declared name (not the file name) so renaming a file cannot silently
        # change the default; if that table is absent the first entry stays
        # selected.
        for i in range(self._table_combo.count()):
            if self._table_combo.itemText(i) == _DEFAULT_TABLE_NAME:
                self._table_combo.setCurrentIndex(i)
                return

    def _calculate(self):
        table = self._tables.get(self._table_combo.currentText())
        if table is None:
            msgbox(self._msg['cmp_title'], self._msg['cmp_no_table'], 1)
            return
        try:
            calls1 = self._get_calls_for_combo(self._prof1_combo)
            calls2 = self._get_calls_for_combo(self._prof2_combo)
        except Exception as exc:
            msgbox(self._msg['cmp_title'], str(exc), 2)
            return

        theta = self._theta_spin.value()
        if self._mode_id.isChecked():
            result = compare_identity(calls1, calls2, table, theta)
        else:
            rel = RELATIONSHIPS[self._rel_combo.currentText()]
            result = compare_kinship(calls1, calls2, table, rel, theta)

        self._last_result = result
        self._show_result(result)

    def _show_result(self, result):
        self._lr_val.setText(_fmt_lr(result.combined_stat))
        self._log10_val.setText(
            f'{result.log10_stat:.2f}' if isfinite(result.log10_stat) else '—')
        self._concl_val.setText(result.verbal_scale)
        self._loci_val.setText(
            self._msg['cmp_loci_stat'].format(
                n=result.n_loci, m=result.n_excluded)
            + (f'\nOnly in Profile 1: {result.n_only_q}   '
               f'Only in Profile 2: {result.n_only_r}'
               if result.n_only_q or result.n_only_r else ''))

        grey = QColor(150, 150, 150)
        try:
            no_edit = Qt.ItemIsEditable
        except AttributeError:
            no_edit = Qt.ItemFlag.ItemIsEditable

        self._tbl.setRowCount(len(result.loci))
        for row, l in enumerate(result.loci):
            q_str = '/'.join(str(a) for a in l.alleles_q if a is not None)
            r_str = '/'.join(str(a) for a in l.alleles_r if a is not None)
            p2_str = f'{l.freq_a2:.4f}' if l.freq_a2 is not None else '—'
            stat_str = f'{l.locus_stat:.4f}' if l.included else 'excl.'
            for col, text in enumerate(
                    [l.marker, q_str, r_str,
                     f'{l.freq_a1:.4f}', p2_str, stat_str, l.note]):
                item = QTableWidgetItem(text)
                item.setFlags(item.flags() & ~no_edit)
                if not l.included:
                    item.setForeground(grey)
                self._tbl.setItem(row, col, item)

        self._results_widget.setVisible(True)
        self.adjustSize()

    def _export_csv(self):
        if self._last_result is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, self._msg['cmp_save_csv'],
            'comparison.csv', 'CSV (*.csv)')
        if not path:
            return
        with open(path, 'w', encoding='utf-8', newline='') as f:
            f.write(export_comparison_csv(self._last_result))

    def _export_pdf(self):
        if self._last_result is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, self._msg['cmp_save_pdf'],
            'comparison.pdf', 'PDF (*.pdf)')
        if not path:
            return
        try:
            export_comparison_pdf(self._last_result, path)
        except Exception as exc:
            msgbox(self._msg['cmp_title'], str(exc), 2)
