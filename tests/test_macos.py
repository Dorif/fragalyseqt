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

# The macOS menu-bar name comes from the bundle owning the process, not from
# Qt.  These tests run everywhere: the platform-specific assertions are
# skipped off macOS, while the no-op and never-raise guarantees -- the parts
# that protect Linux, Windows and the BSDs -- are checked on every platform.

import sys

import pytest

from fragalyseqt import macos


def test_is_macos_follows_the_platform(monkeypatch):
    for name in ('linux', 'win32', 'freebsd13', 'openbsd7'):
        monkeypatch.setattr(macos, 'platform', name)
        assert macos.is_macos() is False, name
    monkeypatch.setattr(macos, 'platform', 'darwin')
    assert macos.is_macos() is True


@pytest.mark.parametrize('other', ['linux', 'win32', 'freebsd13'])
def test_no_op_on_other_platforms(other, monkeypatch):
    # Off macOS the function must do nothing at all: no Objective-C lookup,
    # no exception, and a falsey result so callers can tell it did not act.
    monkeypatch.setattr(macos, 'platform', other)

    def must_not_be_called():
        raise AssertionError('the Objective-C runtime was touched on ' + other)

    monkeypatch.setattr(macos, '_objc', must_not_be_called)
    assert macos.set_application_name() is False


def test_a_broken_runtime_never_raises(monkeypatch):
    # A cosmetic menu title must never be able to stop the application from
    # starting, so every failure path is swallowed and reported as False.
    monkeypatch.setattr(macos, 'platform', 'darwin')
    monkeypatch.setattr(macos, '_OBJC', None)

    def boom(*args, **kwargs):
        raise OSError('libobjc unavailable')

    monkeypatch.setattr(macos, '_objc', boom)
    assert macos.set_application_name() is False
    assert macos.bundle_name() is None


@pytest.mark.skipif(sys.platform != 'darwin', reason='macOS only')
def test_bundle_name_is_patched_on_macos():
    # Regression: the interpreter's own bundle is called "Python", which is
    # what the menu bar next to the apple showed.  Qt cannot override it --
    # the name has to be changed in the loaded bundle's info dictionary.
    assert macos.set_application_name() is True
    assert macos.bundle_name() == macos.APP_NAME
    # Both keys matter: the menu bar reads CFBundleName, Finder prefers
    # CFBundleDisplayName.
    assert macos._info_value('CFBundleDisplayName') == macos.APP_NAME


@pytest.mark.skipif(sys.platform != 'darwin', reason='macOS only')
def test_patching_is_idempotent_and_restorable():
    # main() may be re-entered in tests and embedded uses; applying the name
    # twice must stay correct rather than accumulate anything.
    assert macos.set_application_name('ZZZ-Probe') is True
    assert macos.bundle_name() == 'ZZZ-Probe'
    assert macos.set_application_name() is True
    assert macos.set_application_name() is True
    assert macos.bundle_name() == macos.APP_NAME
