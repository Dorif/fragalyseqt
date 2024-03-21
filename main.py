# This file is part of FragalyseQt.
#
# FragalyseQt is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# FragalyseQt is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with FragalyseQt. If not, see <https://www.gnu.org/licenses/>.
from pyqtgraph.Qt.QtWidgets import QMainWindow
from pyqtgraph import mkQApp
import fragalyseqt
class FragalyseApp(QMainWindow, fragalyseqt.Ui_MainWindow):
    def __init__(self, parent=None):
        super(FragalyseApp, self).__init__(parent)
        self.setupUi(self)
def main():
    app = mkQApp()
    form = FragalyseApp()
    form.show()
    app.exec()
if __name__ == '__main__':
    main()
