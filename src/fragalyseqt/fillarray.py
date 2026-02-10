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

from struct import unpack


def _safe_find(data, pattern):
    pos = data.find(pattern)
    if pos == -1:
        raise ValueError(f"HID pattern {pattern!r} not found")
    return pos


def fill_num_array(FA_dict, FA_file, FA_str, namearray, hexarray, plen, fstr):
    for item in namearray:
        FA_dict[item] = ()
    index = 0
    while index < len(namearray):
        FA_file.seek(_safe_find(FA_str, hexarray[index]) + 20, 0)
        FA_file.seek(int.from_bytes(FA_file.read(4), 'big'), 0)
        if fstr == '>H':
            readnum = 2
        elif fstr == '>I':
            readnum = 4
        elif fstr == '>d':
            readnum = 8
        FA_dict[namearray[index]] = unpack(f'>{plen}{fstr[-1]}',
                                           FA_file.read(readnum * plen))
        index += 1


def fill_char_array(FA_dict, FA_file, FA_str, namearray, dyestr, wavestr):
    array = [b'\x44\x41\x54\x41\x00\x00\x00\x01',
             b'\x44\x41\x54\x41\x00\x00\x00\x02',
             b'\x44\x41\x54\x41\x00\x00\x00\x03',
             b'\x44\x41\x54\x41\x00\x00\x00\x04',
             b'\x44\x41\x54\x41\x00\x00\x00\x69',
             b'\x44\x41\x54\x41\x00\x00\x00\x6a',
             b'\x44\x41\x54\x41\x00\x00\x00\x6b',
             b'\x44\x41\x54\x41\x00\x00\x00\x6c']
    carray = [b'\x44\x79\x65\x57\x00\x00\x00\x01',
              b'\x44\x79\x65\x57\x00\x00\x00\x02',
              b'\x44\x79\x65\x57\x00\x00\x00\x03',
              b'\x44\x79\x65\x57\x00\x00\x00\x04',
              b'\x44\x79\x65\x57\x00\x00\x00\x05',
              b'\x44\x79\x65\x57\x00\x00\x00\x06',
              b'\x44\x79\x65\x57\x00\x00\x00\x07',
              b'\x44\x79\x65\x57\x00\x00\x00\x08']
    darray = [b'\x44\x79\x65\x4e\x00\x00\x00\x01',
              b'\x44\x79\x65\x4e\x00\x00\x00\x02',
              b'\x44\x79\x65\x4e\x00\x00\x00\x03',
              b'\x44\x79\x65\x4e\x00\x00\x00\x04',
              b'\x44\x79\x65\x4e\x00\x00\x00\x05',
              b'\x44\x79\x65\x4e\x00\x00\x00\x06',
              b'\x44\x79\x65\x4e\x00\x00\x00\x07',
              b'\x44\x79\x65\x4e\x00\x00\x00\x08']
    sarray = [b'\x44\x79\x65\x53\x00\x00\x00\x01',
              b'\x44\x79\x65\x53\x00\x00\x00\x02',
              b'\x44\x79\x65\x53\x00\x00\x00\x03',
              b'\x44\x79\x65\x53\x00\x00\x00\x04',
              b'\x44\x79\x65\x53\x00\x00\x00\x05',
              b'\x44\x79\x65\x53\x00\x00\x00\x06',
              b'\x44\x79\x65\x53\x00\x00\x00\x07',
              b'\x44\x79\x65\x53\x00\x00\x00\x08']
    FA_file.seek(_safe_find(FA_str, array[0]) + 12, 0)
    datalength = int.from_bytes(FA_file.read(4), 'big')
    FA_keys = FA_dict.keys()
    for item in namearray:
        if item in FA_keys:
            FA_dict[item] = ()
    index = 0
    while index < FA_dict["Dye#1"]:
        FA_file.seek(_safe_find(FA_str, array[index]) + 20, 0)
        FA_file.seek(int.from_bytes(FA_file.read(4), 'big'), 0)
        FA_dict[namearray[index]] = unpack(f'>{datalength}h',
                                           FA_file.read(2 * datalength))
        if carray[index] in FA_keys:
            FA_file.seek(_safe_find(FA_str, carray[index]) + 20, 0)
            FA_dict[wavestr[index]] = int.from_bytes(FA_file.read(2), 'big')
        FA_file.seek(_safe_find(FA_str, darray[index]) + 16, 0)
        nlen = int.from_bytes(FA_file.read(4), 'big') - 1
        FA_file.seek(int.from_bytes(FA_file.read(4), 'big') + 1, 0)
        FA_dict[dyestr[index]] = FA_file.read(nlen) if nlen >= 0 else b''
        if len(FA_dict[dyestr[index]]) == 0:
            FA_file.seek(_safe_find(FA_str, sarray[index]) + 20, 0)
            slen = int.from_bytes(FA_file.read(1), 'big')
            FA_dict[dyestr[index]] = FA_file.read(slen)
            if FA_dict[dyestr[index]] == b'':
                FA_dict[dyestr[index]] = bytes(str(index+1), 'utf-8')
        index += 1
