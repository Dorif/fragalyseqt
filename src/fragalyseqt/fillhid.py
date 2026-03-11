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

from . import fillarray
from .fillarray import _safe_find

UDATAC  = ["DATA1", "DATA2", "DATA3", "DATA4",
           "DATA105", "DATA106", "DATA107", "DATA108"]
_DYEN   = ["DyeN1", "DyeN2", "DyeN3", "DyeN4",
           "DyeN5", "DyeN6", "DyeN7", "DyeN8"]
_WAVELNG = ["DyeW1", "DyeW2", "DyeW3", "DyeW4",
            "DyeW5", "DyeW6", "DyeW7", "DyeW8"]


def parse_hid(fname, tmpabif, ifacemsg):
    """Parse a HID (ABIF) file whose channel data was not decoded by BioPython.
    Fills tmpabif in-place with MODL1, Peak arrays, Dye#1 count, and channel
    data. Shows an error dialog and re-raises on any failure."""
    try:
        HIDfile = open(fname, "rb")
        s = HIDfile.read()
        tmpkeys = tmpabif.keys()
        MODL1_hex = b'\x4d\x4f\x44\x4c\x00\x00\x00\x01'
        if "MODL1" in tmpkeys and tmpabif["MODL1"] is None:
            modl1_pos = _safe_find(s, MODL1_hex)
            HIDfile.seek(modl1_pos + 12, 0)
            namesize = int.from_bytes(HIDfile.read(4), 'big')
            # Old ABI 310 software may use a different record format.
            if namesize < 4:
                HIDfile.seek(modl1_pos + 10, 0)
                namesize = int.from_bytes(HIDfile.read(2), 'big')
                HIDfile.seek(8, 1)
            else:
                HIDfile.seek(4, 1)
            tmpabif["MODL1"] = HIDfile.read(namesize)
        if "Peak1" in tmpkeys and tmpabif["Peak1"] is None:
            pshorthexarray = [
                b'\x50\x65\x61\x6b\x00\x00\x00\x01',
                b'\x50\x65\x61\x6b\x00\x00\x00\x05']
            pinthexarray = [
                b'\x50\x65\x61\x6b\x00\x00\x00\x02',
                b'\x50\x65\x61\x6b\x00\x00\x00\x03',
                b'\x50\x65\x61\x6b\x00\x00\x00\x04',
                b'\x50\x65\x61\x6b\x00\x00\x00\x07',
                b'\x50\x65\x61\x6b\x00\x00\x00\x08',
                b'\x50\x65\x61\x6b\x00\x00\x00\x09',
                b'\x50\x65\x61\x6b\x00\x00\x00\x0a']
            pdoublehexarray = [
                b'\x50\x65\x61\x6b\x00\x00\x00\x06',
                b'\x50\x65\x61\x6b\x00\x00\x00\x0b',
                b'\x50\x65\x61\x6b\x00\x00\x00\x0c',
                b'\x50\x65\x61\x6b\x00\x00\x00\x0d',
                b'\x50\x65\x61\x6b\x00\x00\x00\x0e',
                b'\x50\x65\x61\x6b\x00\x00\x00\x0f',
                b'\x50\x65\x61\x6b\x00\x00\x00\x10',
                b'\x50\x65\x61\x6b\x00\x00\x00\x11',
                b'\x50\x65\x61\x6b\x00\x00\x00\x12',
                b'\x50\x65\x61\x6b\x00\x00\x00\x15']
            pshortname  = ["Peak1", "Peak5"]
            pintname    = ["Peak2", "Peak3", "Peak4", "Peak7",
                           "Peak8", "Peak9", "Peak10"]
            pdoublename = ["Peak6",  "Peak11", "Peak12", "Peak13",
                           "Peak14", "Peak15", "Peak16", "Peak17",
                           "Peak18", "Peak21"]
            HIDfile.seek(_safe_find(s, pshorthexarray[0]) + 12, 0)
            peakarraylen = int.from_bytes(HIDfile.read(4), 'big')
            fillarray.fill_num_array(tmpabif, HIDfile, s,
                                     pshortname, pshorthexarray,
                                     peakarraylen, '>H')
            fillarray.fill_num_array(tmpabif, HIDfile, s,
                                     pintname, pinthexarray,
                                     peakarraylen, '>I')
            fillarray.fill_num_array(tmpabif, HIDfile, s,
                                     pdoublename, pdoublehexarray,
                                     peakarraylen, '>d')
        if tmpabif["Dye#1"] is None:
            dyenum = b'\x44\x79\x65\x23\x00\x00\x00\x01'
            HIDfile.seek(_safe_find(s, dyenum) + 20, 0)
            tmpabif["Dye#1"] = min(int.from_bytes(HIDfile.read(2), 'big'), 8)
            for i in range(tmpabif["Dye#1"]):
                tmpabif[UDATAC[i]] = None
        fillarray.fill_char_array(tmpabif, HIDfile, s, UDATAC, _DYEN, _WAVELNG)
        HIDfile.close()
    except Exception:
        HIDfile.close()
        raise
