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

from pyqtgraph.Qt.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel,
                                    QLineEdit, QComboBox, QPlainTextEdit,
                                    QPushButton, QDialogButtonBox)


class SizeStandardEditor(QDialog):
    def __init__(self, parent=None, ifacemsg=None, dyes=None):
        super().__init__(parent)
        self.ifacemsg = ifacemsg
        self.dyes = dyes or []
        self.setWindowTitle(ifacemsg.get('addsizestd', 'Add size standard'))
        self.setMinimumWidth(400)
        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout(self)
        self.setStyleSheet("QLineEdit, QComboBox, QPlainTextEdit { background-color: #cacaca; color: black; }")

        # Name
        name_layout = QHBoxLayout()
        name_layout.addWidget(QLabel(self.ifacemsg.get('stdname',
                                                       'Standard Name:')))
        self.name_edit = QLineEdit()
        name_layout.addWidget(self.name_edit)
        layout.addLayout(name_layout)

        # Channel/Dye
        channel_layout = QHBoxLayout()
        channel_layout.addWidget(QLabel(self.ifacemsg.get('stdchannel',
                                                          'ILS Channel (ABIF key):')))
        self.channel_combo = QComboBox()

        self.channel_map = {
            "Blue": "DATA1",
            "Green": "DATA2",
            "Yellow": "DATA3",
            "Red": "DATA4",
            "Orange": "DATA105",
            "Purple": "DATA106",
            "Aqua": "DATA107",
            "Brown": "DATA108"
        }

        self.channel_combo.addItems(list(self.channel_map.keys()))
        self.channel_combo.setEditable(True)
        channel_layout.addWidget(self.channel_combo)
        layout.addLayout(channel_layout)

        # Sizes
        layout.addWidget(QLabel(self.ifacemsg.get('stdsizes',
                                                  'Sizes (space separated):')))
        self.sizes_edit = QPlainTextEdit()
        self.sizes_edit.setPlaceholderText("20 40 60 80 ...")
        layout.addWidget(self.sizes_edit)

        # Buttons
        if hasattr(QDialogButtonBox, 'StandardButton'):
            btn = QDialogButtonBox.StandardButton
        else:
            btn = QDialogButtonBox
        btns = btn.Save | btn.Cancel

        self.button_box = QDialogButtonBox(btns)
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)
        layout.addWidget(self.button_box)

    def get_data(self):
        text = self.channel_combo.currentText()
        channel = self.channel_map.get(text, text)
        return {
            'name': self.name_edit.text().strip(),
            'channel': channel,
            'sizes': self.sizes_edit.toPlainText().strip()
        }
