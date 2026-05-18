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

# Dialogs for session save, open, and file verification.
from pyqtgraph.Qt import QtWidgets, QtCore
from .localize import localizefq

_msg: dict = {}
localizefq(_msg)


def _std_buttons(*kinds):
    # Return QDialogButtonBox.StandardButton flags compatible with all
    # Qt bindings.
    try:
        bb = QtWidgets.QDialogButtonBox
        flags = bb.Ok  # PyQt5 / PySide6 old-style
        _ = flags | bb.Cancel
    except AttributeError:
        bb = QtWidgets.QDialogButtonBox
        flags = None
    result = None
    for kind in kinds:
        try:
            flag = getattr(bb, kind)
        except AttributeError:
            flag = getattr(bb.StandardButton, kind)
        result = flag if result is None else result | flag
    return result


class SaveSessionDialog(QtWidgets.QDialog):
    # Ask the user for a session name before saving.

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(_msg.get('savesessiontitle', 'Save Session'))
        self.setMinimumWidth(340)
        layout = QtWidgets.QFormLayout(self)

        self._name = QtWidgets.QLineEdit()
        from datetime import datetime
        self._name.setText(datetime.now().strftime('%Y-%m-%d %H:%M'))
        layout.addRow(_msg.get('sessionnamelabel', 'Session name:'),
                      self._name)

        buttons = QtWidgets.QDialogButtonBox(_std_buttons('Ok', 'Cancel'))
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addRow(buttons)

    def get_name(self) -> str:
        return self._name.text().strip()


class OpenSessionDialog(QtWidgets.QDialog):
    # Let the user pick a previously saved session.

    def __init__(self, sessions: list[dict], parent=None):
        super().__init__(parent)
        self.setWindowTitle(_msg.get('opensessiontitle', 'Open Session'))
        self.setMinimumSize(420, 260)
        layout = QtWidgets.QVBoxLayout(self)

        self._list = QtWidgets.QListWidget()
        self._sessions = sessions
        for s in sessions:
            ts = s['created_at'][:16].replace('T', ' ')
            self._list.addItem(f"{s['name']}  —  {ts}")
        self._list.setCurrentRow(0)
        self._list.itemDoubleClicked.connect(self.accept)
        layout.addWidget(self._list)

        buttons = QtWidgets.QDialogButtonBox(_std_buttons('Ok', 'Cancel'))
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def get_selected_session_id(self) -> int | None:
        row = self._list.currentRow()
        if row < 0 or row >= len(self._sessions):
            return None
        return self._sessions[row]['id']


class VerificationDialog(QtWidgets.QDialog):
    # Show per-file verification results and offer to open read-only.

    def __init__(self, statuses: list, parent=None):
        super().__init__(parent)
        self.setWindowTitle(_msg.get('sessionverifytitle',
                                     'File Verification'))
        self.setMinimumSize(560, 300)
        layout = QtWidgets.QVBoxLayout(self)

        table = QtWidgets.QTableWidget(len(statuses), 3)
        table.setHorizontalHeaderLabels(['File', 'Status', 'Detail'])
        table.horizontalHeader().setStretchLastSection(True)
        table.verticalHeader().setVisible(False)
        table.setEditTriggers(
            QtWidgets.QAbstractItemView.NoEditTriggers
            if hasattr(QtWidgets.QAbstractItemView, 'NoEditTriggers')
            else QtWidgets.QAbstractItemView.EditTrigger.NoEditTriggers)
        table.setSelectionMode(
            QtWidgets.QAbstractItemView.NoSelection
            if hasattr(QtWidgets.QAbstractItemView, 'NoSelection')
            else QtWidgets.QAbstractItemView.SelectionMode.NoSelection)

        for row, st in enumerate(statuses):
            table.setItem(row, 0, QtWidgets.QTableWidgetItem(st.file_name))
            icon = '✓' if st.status == 'ok' else '✗'
            status_item = QtWidgets.QTableWidgetItem(f'{icon}  {st.status}')
            table.setItem(row, 1, status_item)
            table.setItem(row, 2, QtWidgets.QTableWidgetItem(st.detail))

        table.resizeColumnsToContents()
        layout.addWidget(table)

        msg = QtWidgets.QLabel(
            _msg.get('sessionverifymsg',
                     'One or more source files are missing or have been modified.\n'
                     'The session will open in read-only mode using data stored '
                     'in the database.'))
        msg.setWordWrap(True)
        layout.addWidget(msg)

        btn_layout = QtWidgets.QHBoxLayout()
        btn_layout.addStretch()
        self._open_btn = QtWidgets.QPushButton(
            _msg.get('openreadonly', 'Open read-only'))
        cancel_btn = QtWidgets.QPushButton('Cancel')
        self._open_btn.clicked.connect(self.accept)
        cancel_btn.clicked.connect(self.reject)
        btn_layout.addWidget(self._open_btn)
        btn_layout.addWidget(cancel_btn)
        layout.addLayout(btn_layout)
