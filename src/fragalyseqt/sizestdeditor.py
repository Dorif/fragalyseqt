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

import xml.etree.ElementTree as ET
from .boxes import msgbox
from pyqtgraph.Qt.QtWidgets import (QDialog, QVBoxLayout, QHBoxLayout, QLabel,
                                    QComboBox, QPlainTextEdit, QFileDialog,
                                    QDialogButtonBox, QPushButton, QLineEdit,)

# dyeIndex values found in SizeStandardContainer XML files (GeneMapper / HID)
# map to the ABIF channel key used internally.
_DYE_INDEX_TO_CHANNEL = {
    1: 'DATA1',
    2: 'DATA2',
    3: 'DATA3',
    4: 'DATA4',
    5: 'DATA105',
    6: 'DATA106',
    7: 'DATA107',
    8: 'DATA108',
}


def _parse_sizestd_xml(path):
    # Parse a size-standard XML file and return (name, channel, sizes_str)
    # or raise ValueError with a human-readable message on failure.
    #
    # Two formats are supported:
    #
    # 1. External GeneMapper / HID SizeStandardContainer format:
    #      <SizeStandardContainer>
    #        <xmlSizeStandard>
    #          <sizeStandardName>...</sizeStandardName>
    #          <dyeIndex>5</dyeIndex>
    #          <sizeStdDefinition>80.0</sizeStdDefinition>
    #          ...
    #        </xmlSizeStandard>
    #      </SizeStandardContainer>
    #
    # 2. FragalyseQt internal sizestandards.xml format:
    #      <Ladder name="GS600LIZ" channel="DATA105">
    #        <Sizes>20 40 60 80 ...</Sizes>
    #      </Ladder>
    tree = ET.parse(path)
    root = tree.getroot()
    tag = root.tag.split('}')[-1]  # strip namespace if present

    # --- Format 1: SizeStandardContainer ---
    if tag == 'SizeStandardContainer':
        ns_prefix = ''
        if '}' in root.tag:
            ns_prefix = root.tag.split('}')[0] + '}'
        std = root.find(f'{ns_prefix}xmlSizeStandard')
        if std is None:
            std = root.find('xmlSizeStandard')
        if std is None:
            raise ValueError("No <xmlSizeStandard> element found.")
        name_el = std.find('sizeStandardName')
        dye_el = std.find('dyeIndex')
        size_els = std.findall('sizeStdDefinition')
        if name_el is None or not size_els:
            raise ValueError("Missing name or size definitions.")
        name = (name_el.text or '').strip()
        try:
            dye_idx = int((dye_el.text or '5').strip())
        except ValueError:
            dye_idx = 5
        channel = _DYE_INDEX_TO_CHANNEL.get(dye_idx, f'DATA{dye_idx}')
        sizes = []
        for el in size_els:
            try:
                sizes.append(str(int(float(el.text))))
            except (TypeError, ValueError):
                pass
        if not sizes:
            raise ValueError("No valid size values found.")
        return name, channel, ' '.join(sizes)

    # --- Format 2: single <Ladder> element (internal format) ---
    if tag == 'Ladder':
        ladder = root
    else:
        ladder = root.find('Ladder')
    if ladder is not None:
        name = ladder.get('name', '')
        channel = ladder.get('channel', 'DATA105')
        sizes_el = ladder.find('Sizes')
        sizes_str = (sizes_el.text or '').strip() if sizes_el is not None else ''
        if not sizes_str:
            raise ValueError("No <Sizes> content in <Ladder> element.")
        return name, channel, sizes_str

    raise ValueError("Unrecognised size-standard XML format.")


class SizeStandardEditor(QDialog):
    def __init__(self, parent=None, ifacemsg=None, dyes=None):
        super().__init__(parent)
        self.ifacemsg = ifacemsg or {}
        self.dyes = dyes or []
        self.setWindowTitle(self.ifacemsg.get('addsizestd'))
        self.setMinimumWidth(400)
        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout(self)
        self.setStyleSheet(
            "QLineEdit, QComboBox, QPlainTextEdit { "
            "background-color: #cacaca; color: black; }")

        # Name
        name_layout = QHBoxLayout()
        name_layout.addWidget(QLabel(self.ifacemsg.get('stdname')))
        self.name_edit = QLineEdit()
        name_layout.addWidget(self.name_edit)
        layout.addLayout(name_layout)

        # Channel/Dye
        channel_layout = QHBoxLayout()
        channel_layout.addWidget(QLabel(self.ifacemsg.get('stdchannel')))
        self.channel_combo = QComboBox()

        self.channel_map = {
            "Blue":   "DATA1",
            "Green":  "DATA2",
            "Yellow": "DATA3",
            "Red":    "DATA4",
            "Orange": "DATA105",
            "Purple": "DATA106",
            "Aqua":   "DATA107",
            "Brown":  "DATA108",
        }
        self._rev_channel_map = {v: k for k, v in self.channel_map.items()}

        self.channel_combo.addItems(list(self.channel_map.keys()))
        self.channel_combo.setEditable(True)
        channel_layout.addWidget(self.channel_combo)
        layout.addLayout(channel_layout)

        # Sizes
        layout.addWidget(QLabel(self.ifacemsg.get('stdsizes')))
        self.sizes_edit = QPlainTextEdit()
        self.sizes_edit.setPlaceholderText("20 40 60 80 ...")
        layout.addWidget(self.sizes_edit)

        # Save / Cancel
        if hasattr(QDialogButtonBox, 'StandardButton'):
            btn = QDialogButtonBox.StandardButton
        else:
            btn = QDialogButtonBox
        btns = btn.Save | btn.Cancel

        self.button_box = QDialogButtonBox(btns)
        self.button_box.accepted.connect(self.accept)
        self.button_box.rejected.connect(self.reject)

        import_btn = QPushButton(self.ifacemsg.get('stdimportxml'))
        import_btn.clicked.connect(self._import_from_xml)

        bottom_row = QHBoxLayout()
        bottom_row.addWidget(import_btn)
        bottom_row.addWidget(self.button_box)
        layout.addLayout(bottom_row)

    def _import_from_xml(self):
        path, _ = QFileDialog.getOpenFileName(
            self, self.ifacemsg.get('stdimportxml'),
            '', 'XML files (*.xml);;All files (*)')
        if not path:
            return
        try:
            name, channel, sizes_str = _parse_sizestd_xml(path)
        except Exception:
            msgbox(self.ifacemsg.get('stdimportxml'),
                   self.ifacemsg.get('stdimporterr'), 1)
            return
        self.name_edit.setText(name)
        # Select matching colour name if known, otherwise set raw channel key
        colour = self._rev_channel_map.get(channel)
        if colour:
            idx = self.channel_combo.findText(colour)
            if idx >= 0:
                self.channel_combo.setCurrentIndex(idx)
            else:
                self.channel_combo.setEditText(colour)
        else:
            self.channel_combo.setEditText(channel)
        self.sizes_edit.setPlainText(sizes_str)

    def get_data(self):
        text = self.channel_combo.currentText()
        channel = self.channel_map.get(text, text)
        return {
            'name': self.name_edit.text().strip(),
            'channel': channel,
            'sizes': self.sizes_edit.toPlainText().strip()
        }
