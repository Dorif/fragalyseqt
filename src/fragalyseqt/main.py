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

from pyqtgraph.Qt.QtWidgets import QMainWindow, QApplication
from . import fragalyseqt as fragalyseqt_ui
from . import __version__
from .macos import set_application_name
from sys import argv


class FragalyseApp(QMainWindow, fragalyseqt_ui.Ui_MainWindow):
    def __init__(self, parent=None):
        super(FragalyseApp, self).__init__(parent)
        self.setupUi(self)


def configure_application(app):
    """Identify the application to the desktop environment.

    Must run *before* any window is created.  Without it the process is
    only known by its interpreter, so macOS shows "Python" in the global
    menu bar and the generic Python icon in the Dock, and GNOME/KDE cannot
    match the window to its .desktop entry.  The icon is set on the
    application, not merely on the window: the Dock and the taskbar read
    the application icon, while setWindowIcon() on a window only decorates
    that window's own title bar and switcher entry.
    """
    app.setApplicationName("FragalyseQt")
    app.setApplicationDisplayName("FragalyseQt")
    app.setApplicationVersion(__version__)
    app.setOrganizationName("FragalyseQt")
# Matches packaging/fragalyseqt.desktop, so Wayland/GNOME associates the
# window with the installed launcher and shows its name and icon.
    app.setDesktopFileName("fragalyseqt")
    app.setWindowIcon(fragalyseqt_ui.app_icon())
    return app


def main():
    # macOS reads the name next to the apple menu from the bundle that owns
    # the process, not from Qt, so an unbundled interpreter shows "Python"
    # there no matter what setApplicationDisplayName() says.  This has to
    # happen before the QApplication exists, because Qt builds the Cocoa
    # menu bar from the bundle information while starting up.  No-op
    # elsewhere.
    set_application_name()
    app = QApplication(argv)
    configure_application(app)
    form = FragalyseApp()
    form.show()
    app.exec()


if __name__ == '__main__':
    main()
