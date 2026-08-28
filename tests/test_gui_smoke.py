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

# GUI smoke tests for the main window.
#
# Strategy: drive FragalyseApp exactly like a real user would through
# pytest-qt's `qtbot`, but bypass the native file-picker dialog by calling
# the same internal entry point open_and_plot() delegates to
# (FragalyseApp._load_file(path)) directly with a path to one of the real
# NIST reference .hid files already checked into docs/TEST_FILES. Everything
# downstream of that call — tab creation, peak finding, sizing, table
# population, widget wiring — runs exactly as in production.
#
# Requires a Qt binding (PyQt6/PyQt5/PySide6) and pytest-qt to be installed,
# see pyproject.toml's `test` extra. Runs headless via the Qt "offscreen"
# platform plugin, so no X server / display is needed (works in CI).

import os
import pytest

pytest.importorskip("pytestqt")

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from pyqtgraph.Qt.QtCore import Qt  # noqa: E402
from pyqtgraph.Qt.QtGui import QAction  # noqa: E402
from fragalyseqt.main import FragalyseApp, configure_application  # noqa: E402
import fragalyseqt  # noqa: E402
from fragalyseqt import fragalyseqt as fq
from fragalyseqt.fragalyseqt import ifacemsg

_KIT_DIR = os.path.join(
    os.path.dirname(__file__), "..", "docs", "TEST_FILES", "NIST",
    "Single_Source", "Promega PowerPlex Fusion 6C",
)

# 4-dye ABI310 MLPA runs: they carry no DATA105, so a LIZ size standard
# cannot be sized against them -- used to exercise the batch failure path.
_MLPA_DIR = os.path.join(
    os.path.dirname(__file__), "..", "docs", "TEST_FILES", "ABI310", "MLPA",
)

_SAMPLE_HID = os.path.join(_KIT_DIR, "PPF6C_2800M.hid")
# Two more real runs from the same kit, used by the batch-processing tests.
_SAMPLE_HID_2 = os.path.join(_KIT_DIR, "PPF6C_NTD01.hid")
_SAMPLE_HID_3 = os.path.join(_KIT_DIR, "PPF6C_NTD02.hid")


@pytest.fixture
def main_window(qtbot):
    window = FragalyseApp()
    qtbot.addWidget(window)
    window.show()
    return window


def test_main_window_starts_with_no_tabs(main_window):
    assert main_window.file_tab.count() == 0
    assert main_window.file_states == []


def test_load_real_hid_file_creates_populated_tab(main_window):
    idx = main_window._load_file(_SAMPLE_HID)

    assert idx == 0
    assert main_window.file_tab.count() == 1
    assert main_window.file_tab.tabText(0) == "PPF6C_2800M.hid"

    state = main_window.file_states[idx]
    assert len(state.peakpositions) > 0
    assert len(state.Dye) == 6  # 6-dye Fusion 6C kit
    assert state.table_widget.rowCount() == len(state.peakpositions)


def test_toggling_channel_checkbox_hides_channel(qtbot, main_window):
    idx = main_window._load_file(_SAMPLE_HID)
    state = main_window.file_states[idx]
    checkbox = state.hidech[0]

    assert state.show_channels[0] == 1
    qtbot.mouseClick(checkbox, Qt.MouseButton.LeftButton)
    assert state.show_channels[0] == 0


def test_min_height_spinbox_triggers_reanalysis(qtbot, main_window):
    idx = main_window._load_file(_SAMPLE_HID)
    state = main_window.file_states[idx]
    before = len(state.peakpositions)

    # Raising the minimum peak height should only ever reduce (or keep
    # equal) the number of detected peaks -- never increase it.
    state.getheight.setValue(state.getheight.value() + 500)
    qtbot.wait(50)

    after = len(main_window.file_states[idx].peakpositions)
    assert after <= before


def test_batch_propagates_sizing_without_a_panel(main_window):
    # Regression: batch processing used to be a no-op unless a panel had been
    # selected.  SizeCall is a momentary trigger -- reanalyse() consumes it
    # and unchecks the button -- so process_whole_batch() read it back as
    # False and sent sizecall=False to every tab; the only thing left that
    # could still switch sizing on was reanalyse()'s "a panel is selected"
    # branch.  With no panel the other tabs got no sizing at all, which is
    # exactly the "just measure fragment sizes" workflow.
    src_idx = main_window._load_file(_SAMPLE_HID)
    other_idx = main_window._load_file(_SAMPLE_HID_2)
    assert other_idx != src_idx

    source = main_window.file_states[src_idx]
    other = main_window.file_states[other_idx]
    main_window.file_tab.setCurrentIndex(src_idx)

    # No panel selected anywhere.
    assert source.panel_combo.currentText() == ifacemsg["nopanel"]
    assert other.panel_combo.currentText() == ifacemsg["nopanel"]

    # User sizes the current tab; the button unlatches itself afterwards.
    source.sizecall.setChecked(True)
    main_window.reanalyse(source)
    assert len(source.peaksizes) > 0
    assert source.sizecall.isChecked() is False

    assert len(other.peaksizes) == 0

    main_window.process_whole_batch()

    # The other tab must now be sized too, with the source's settings.
    assert len(other.peaksizes) > 0
    assert other.getheight.value() == source.getheight.value()
    assert other.ILS.currentText() == source.ILS.currentText()


def test_batch_propagates_settings_to_every_tab(main_window):
    # A batch run must not stop at the first non-current tab, and must copy
    # the detection parameters (not just sizing) onto all of them.
    src_idx = main_window._load_file(_SAMPLE_HID)
    main_window._load_file(_SAMPLE_HID_2)
    main_window._load_file(_SAMPLE_HID_3)
    assert main_window.file_tab.count() == 3

    source = main_window.file_states[src_idx]
    main_window.file_tab.setCurrentIndex(src_idx)
    source.getheight.setValue(300)
    source.bcd.setChecked(True)
    source.do_BCD = True

    main_window.process_whole_batch()

    for i, state in enumerate(main_window.file_states):
        if i == src_idx:
            continue
        assert state.getheight.value() == 300
        assert state.do_BCD is True
        assert state.table_widget.rowCount() == len(state.peakpositions)


def test_batch_reports_failures_once_and_keeps_going(main_window, monkeypatch):
    # A tab whose sizing fails must not abort the batch, and must not raise
    # one modal dialog per file: failures are collected and reported once.
    # Source is a 6-dye run where GS600LIZ (DATA105) exists, so sizing works;
    # the other tabs are 4-dye MLPA runs with no DATA105, so sizing raises.
    dialogs = []
    monkeypatch.setattr(fq, "msgbox",
                        lambda header, msg, kind: dialogs.append(msg))

    src_idx = main_window._load_file(_SAMPLE_HID)
    for name in ("C10-Control01-105A.fsa", "D10-Control01-105B.fsa",
                 "E10-Control02-105A.fsa"):
        main_window._load_file(os.path.join(_MLPA_DIR, name))
    assert main_window.file_tab.count() == 4

    source = main_window.file_states[src_idx]
    main_window.file_tab.setCurrentIndex(src_idx)
    source.ILS.setCurrentText("GS600LIZ")
    source.sizecall.setChecked(True)
    main_window.reanalyse(source)
    assert len(source.peaksizes) > 0

    dialogs.clear()
    main_window.process_whole_batch()

    # Exactly one summary dialog, naming every tab that failed.
    assert len(dialogs) == 1
    for name in ("C10-Control01-105A.fsa", "D10-Control01-105B.fsa",
                 "E10-Control02-105A.fsa"):
        assert name in dialogs[0]

    # The run went through all tabs and left no batch state behind.
    assert main_window.file_tab.count() == 4
    assert main_window._batch_mode is False
    assert main_window._batch_failures == []


def test_about_entry_stays_in_the_help_menu(main_window):
    # Regression (macOS): Qt guesses a menu role from the action text, and
    # "About" is detected as AboutRole, which moves the item into the
    # application menu.  Help then held nothing, and macOS does not draw an
    # empty menu -- so the entry disappeared from the menu bar.  The guess
    # only matches the English label, so the seven other translations were
    # unaffected, which is what made it look platform-specific rather than
    # language-specific.
    help_menu = None
    for action in main_window.menuBar().actions():
        if action.menu() is not None and ifacemsg["aboutbtn"] in [
                a.text() for a in action.menu().actions()]:
            help_menu = action.menu()
            break

    assert help_menu is not None, "the About entry is not in any menu"

    about = next(a for a in help_menu.actions()
                 if a.text() == ifacemsg["aboutbtn"])
    menu_role = getattr(QAction, 'MenuRole', QAction)
    assert about.menuRole() == menu_role.NoRole, (
        "About must keep NoRole, otherwise macOS relocates it and leaves "
        "the Help menu empty")


def test_no_menu_is_left_empty(main_window):
    # An empty top-level menu is invisible on macOS, so it is always a bug:
    # either the items were never added, or Qt's text heuristic relocated
    # every one of them into the application menu.
    for action in main_window.menuBar().actions():
        menu = action.menu()
        if menu is None:
            continue
        items = [a for a in menu.actions() if not a.isSeparator()]
        assert items, f"menu {action.text()!r} has no items"


def test_icon_is_found_regardless_of_working_directory(tmp_path, monkeypatch):
    # Regression: the icon used to be loaded from the relative path
    # "FragalyseQt.png", which only resolved when the process happened to be
    # started from the source checkout.  Anywhere else -- and after a plain
    # pip install -- the icon silently came out null, so the desktop fell
    # back to the generic interpreter icon.
    monkeypatch.chdir(tmp_path)
    path = fq.icon_path()
    assert os.path.isabs(path), "icon path must not depend on the cwd"
    assert os.path.isfile(path), f"bundled icon missing: {path}"
    assert not fq.app_icon().isNull(), "icon failed to load"


def test_icon_ships_inside_the_package():
    # The application resolves the icon next to its own __file__, so the file
    # has to live in the package itself -- a copy in the repository root is
    # not installed by pip and would leave installed users without an icon.
    packaged = os.path.join(os.path.dirname(fq.__file__), "FragalyseQt.png")
    assert os.path.isfile(packaged), (
        "FragalyseQt.png must ship inside the package, not only in the "
        "repository root")


def test_application_identifies_itself_to_the_desktop(qapp):
    # Regression (macOS Dock / global menu bar): without an application name
    # and an application-level icon the process is only known by its
    # interpreter, so the Dock showed the generic Python icon and the global
    # menu bar read "Python".  setWindowIcon() on the window is not enough --
    # the Dock and the taskbar read the icon set on the application.
    saved = (qapp.applicationName(), qapp.applicationDisplayName(),
             qapp.applicationVersion(), qapp.organizationName(),
             qapp.desktopFileName(), qapp.windowIcon())
    try:
        configure_application(qapp)
        assert qapp.applicationName() == "FragalyseQt"
        assert qapp.applicationDisplayName() == "FragalyseQt"
        assert qapp.applicationVersion() == fragalyseqt.__version__
        # Must match packaging/fragalyseqt.desktop for Wayland/GNOME to
        # associate the window with its installed launcher.
        assert qapp.desktopFileName() == "fragalyseqt"
        assert not qapp.windowIcon().isNull(), (
            "the Dock reads the application icon, which is unset")
    finally:
        # QApplication is process-wide in pytest-qt: restore it so this test
        # cannot leak its state into the others.
        (qapp.setApplicationName(saved[0]),
         qapp.setApplicationDisplayName(saved[1]),
         qapp.setApplicationVersion(saved[2]),
         qapp.setOrganizationName(saved[3]),
         qapp.setDesktopFileName(saved[4]),
         qapp.setWindowIcon(saved[5]))
