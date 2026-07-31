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

from json import load as json_load, dump as json_dump
from os import makedirs
from os.path import join, isfile
from platformdirs import user_data_dir
from pyqtgraph.Qt.QtWidgets import (QDialog, QFormLayout, QSpinBox, QCheckBox,
                                    QLineEdit, QDialogButtonBox)

_USER_DATA = user_data_dir('fragalyseqt', appauthor=False)
_SETTINGS_FILE = join(_USER_DATA, 'soap_settings.json')

_DEFAULTS = {
    'enabled': False,
    'host': '127.0.0.1',
    'port': 8742,
    'token': '',
    'timeout': 30,
}


def load_soap_settings():
    if isfile(_SETTINGS_FILE):
        try:
            with open(_SETTINGS_FILE) as f:
                return {**_DEFAULTS, **json_load(f)}
        except Exception:
            pass
    return dict(_DEFAULTS)


def save_soap_settings(settings):
    makedirs(_USER_DATA, exist_ok=True)
    with open(_SETTINGS_FILE, 'w') as f:
        json_dump(settings, f, indent=2)


class SOAPSettingsDialog(QDialog):
    def __init__(self, settings, parent=None):
        super().__init__(parent)
        self.setWindowTitle('SOAP API Settings')
        layout = QFormLayout(self)

        self._enabled = QCheckBox()
        self._enabled.setChecked(settings['enabled'])
        layout.addRow('Enable SOAP API:', self._enabled)

        self._host = QLineEdit(settings['host'])
        layout.addRow('Bind address:', self._host)

        self._port = QSpinBox()
        self._port.setRange(1024, 65535)
        self._port.setValue(settings['port'])
        layout.addRow('Port:', self._port)

        self._token = QLineEdit(settings['token'])
        self._token.setPlaceholderText('Leave empty for local-only access')
        layout.addRow('API token:', self._token)

        self._timeout = QSpinBox()
        self._timeout.setRange(5, 300)
        self._timeout.setValue(settings['timeout'])
        self._timeout.setSuffix(' s')
        layout.addRow('Request timeout:', self._timeout)

        try:
            std = QDialogButtonBox.Ok | QDialogButtonBox.Cancel
        except AttributeError:
            std = (QDialogButtonBox.StandardButton.Ok |
                   QDialogButtonBox.StandardButton.Cancel)
        buttons = QDialogButtonBox(std)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addRow(buttons)

    def get_settings(self):
        return {
            'enabled': self._enabled.isChecked(),
            'host': self._host.text().strip() or '127.0.0.1',
            'port': self._port.value(),
            'token': self._token.text().strip(),
            'timeout': self._timeout.value(),
        }
