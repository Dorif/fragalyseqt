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

"""macOS-only desktop integration.

The application name shown next to the apple menu does not come from Qt.
macOS reads CFBundleName from the Info.plist of the bundle that owns the
process, and an interpreter that is not wrapped in an .app bundle is owned
by the framework's own Python.app -- whose CFBundleName is literally
"Python".  So the menu bar says "Python" no matter what
QApplication.setApplicationDisplayName() is set to.

Patching the already-loaded bundle's info dictionary before the Cocoa menu
bar is built makes macOS pick up the real name.  This is done through the
Objective-C runtime with ctypes from the standard library, so it adds no
dependency and costs nothing on other platforms.

Everything here is best-effort: a wrong title in the menu bar is cosmetic,
so no failure in this module may ever stop the application from starting.
"""

from sys import platform

APP_NAME = "FragalyseQt"

# objc_msgSend must be re-declared for every signature: it is a variadic
# trampoline, and on 64-bit a wrong restype/argtypes silently truncates
# pointers instead of failing loudly.
_OBJC = None

# Both keys are patched: CFBundleName is what the menu bar uses, while
# CFBundleDisplayName is what Finder and some system dialogs prefer.
_NAME_KEYS = ("CFBundleName", "CFBundleDisplayName")


def is_macos():
    """True when running on macOS."""
    return platform == "darwin"


def _objc():
    """The Objective-C runtime, loaded once and cached."""
    global _OBJC
    if _OBJC is None:
        from ctypes import cdll, c_char_p, c_void_p, util

        lib = cdll.LoadLibrary(util.find_library("objc"))
        lib.objc_getClass.restype = c_void_p
        lib.objc_getClass.argtypes = [c_char_p]
        lib.sel_registerName.restype = c_void_p
        lib.sel_registerName.argtypes = [c_char_p]
        _OBJC = lib
    return _OBJC


def _send(receiver, selector, *args, restype=None, argtypes=None):
    """Send an Objective-C message, declaring the signature every time."""
    from ctypes import c_void_p

    lib = _objc()
    lib.objc_msgSend.restype = c_void_p if restype is None else restype
    lib.objc_msgSend.argtypes = [c_void_p, c_void_p] + list(argtypes or ())
    return lib.objc_msgSend(receiver, lib.sel_registerName(selector), *args)


def _nsstring(text):
    """A Python str as an NSString."""
    from ctypes import c_char_p

    cls = _objc().objc_getClass(b"NSString")
    return _send(cls, b"stringWithUTF8String:", text.encode("utf-8"),
                 argtypes=[c_char_p])


def _pystring(pointer):
    """An NSString back as a Python str (None for nil)."""
    from ctypes import c_char_p

    if not pointer:
        return None
    raw = _send(pointer, b"UTF8String", restype=c_char_p)
    return raw.decode("utf-8") if raw else None


def _info_dictionary():
    """The mutable info dictionary of the bundle owning this process."""
    cls = _objc().objc_getClass(b"NSBundle")
    return _send(_send(cls, b"mainBundle"), b"infoDictionary")


def _info_value(key):
    """Current value of an Info.plist key, as a Python str or None."""
    from ctypes import c_void_p

    return _pystring(_send(_info_dictionary(), b"objectForKey:",
                           _nsstring(key), argtypes=[c_void_p]))


def _set_info_value(key, value):
    """Overwrite an Info.plist key in the loaded bundle."""
    from ctypes import c_void_p

    _send(_info_dictionary(), b"setObject:forKey:",
          _nsstring(value), _nsstring(key),
          argtypes=[c_void_p, c_void_p])


def bundle_name():
    """The CFBundleName macOS currently sees, or None if unavailable.

    Used to verify the patch and to report the state in tests; returns None
    on other platforms and whenever the Objective-C runtime cannot be
    reached.
    """
    try:
        return _info_value(_NAME_KEYS[0])
    except Exception:
        return None


def set_application_name(name=APP_NAME):
    """Make macOS show `name` instead of the interpreter's own name.

    Must be called *before* the QApplication is constructed, because Qt
    builds the Cocoa menu bar from the bundle information at start-up.

    Returns True when the name was applied, False when there was nothing to
    do (another platform) or the attempt failed.  Never raises.
    """
    if not is_macos():
        return False
    try:
        for key in _NAME_KEYS:
            _set_info_value(key, name)
        return _info_value(_NAME_KEYS[0]) == name
    except Exception:
        # A cosmetic menu title is never worth failing a launch over.
        return False
