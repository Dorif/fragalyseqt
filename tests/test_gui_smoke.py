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
from fragalyseqt.main import FragalyseApp  # noqa: E402

_SAMPLE_HID = os.path.join(
    os.path.dirname(__file__), "..", "docs", "TEST_FILES", "NIST",
    "Single_Source", "Promega PowerPlex Fusion 6C", "PPF6C_2800M.hid",
)


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
