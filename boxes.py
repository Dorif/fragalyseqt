# This file is part of FragalyseQt.
#
# FragalyseQt is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published 
# by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# FragalyseQt is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License along with FragalyseQt. If not, see <https://www.gnu.org/licenses/>.
from pyqtgraph.Qt.QtWidgets import QMessageBox
def msgbox(headerstr, msgstr, msgtype):
    mbox = QMessageBox()
    if msgtype == 1:
        mbox.setIcon(QMessageBox.Warning)
    elif msgtype == 2:
        mbox.setIcon(QMessageBox.Critical)
    else:
        mbox.setIcon(QMessageBox.Information)
    mbox.setWindowTitle(headerstr)
    mbox.setText(msgstr)
#Yes, this button is shown by default and is deafult button by default, but better I'll do it in explicit manner.
    mbox.setStandardButtons(QMessageBox.Ok)
    mbox.setDefaultButton(QMessageBox.Ok)
    mbox.exec()
