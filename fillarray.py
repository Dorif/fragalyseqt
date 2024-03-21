# This file is part of FragalyseQt.
#
# FragalyseQt is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# FragalyseQt is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with FragalyseQt. If not, see <https://www.gnu.org/licenses/>.
from struct import unpack
def fill_num_array(FA_dict, FA_file, FA_str, namearray, hexarray, plength, format_str):
    for item in namearray:
        FA_dict[item] = ()
    index = 0
    while index < len(namearray):
        FA_file.seek(FA_str.find(hexarray[index]) + 20, 0)
        FA_file.seek(int.from_bytes(FA_file.read(4), 'big'), 0)
        if format_str == '>H':
            readnum = 2
        elif format_str == '>I':
            readnum = 4
        elif format_str == '>d':
            readnum = 8
        iterator = 0
        while iterator < plength:
            FA_dict[namearray[index]] += unpack(format_str, FA_file.read(readnum))
            iterator += 1
        index += 1
def fill_char_array(FA_dict, FA_file, FA_str, namearray, dyestr, wavestr):
    array = [b'\x44\x41\x54\x41\x00\x00\x00\x01',b'\x44\x41\x54\x41\x00\x00\x00\x02',b'\x44\x41\x54\x41\x00\x00\x00\x03',b'\x44\x41\x54\x41\x00\x00\x00\x04',b'\x44\x41\x54\x41\x00\x00\x00\x69',
    b'\x44\x41\x54\x41\x00\x00\x00\x6a',b'\x44\x41\x54\x41\x00\x00\x00\x6b',b'\x44\x41\x54\x41\x00\x00\x00\x6c']
    carray = [b'\x44\x79\x65\x57\x00\x00\x00\x01',b'\x44\x79\x65\x57\x00\x00\x00\x02',b'\x44\x79\x65\x57\x00\x00\x00\x03',b'\x44\x79\x65\x57\x00\x00\x00\x04',b'\x44\x79\x65\x57\x00\x00\x00\x05',
    b'\x44\x79\x65\x57\x00\x00\x00\x06',b'\x44\x79\x65\x57\x00\x00\x00\x07',b'\x44\x79\x65\x57\x00\x00\x00\x08']
    darray = [b'\x44\x79\x65\x4e\x00\x00\x00\x01',b'\x44\x79\x65\x4e\x00\x00\x00\x02',b'\x44\x79\x65\x4e\x00\x00\x00\x03',b'\x44\x79\x65\x4e\x00\x00\x00\x04',b'\x44\x79\x65\x4e\x00\x00\x00\x05',
    b'\x44\x79\x65\x4e\x00\x00\x00\x06',b'\x44\x79\x65\x4e\x00\x00\x00\x07',b'\x44\x79\x65\x4e\x00\x00\x00\x08']
    FA_file.seek(FA_str.find(array[0]) + 12, 0)
    datalength = int.from_bytes(FA_file.read(4), 'big')
    for item in namearray:
        if item in FA_dict.keys():
            FA_dict[item] = ()
    index = 0
    while index < FA_dict["Dye#1"]:
        FA_file.seek(FA_str.find(array[index]) + 20, 0)
        FA_file.seek(int.from_bytes(FA_file.read(4), 'big'), 0)
        iterator = 0
        while iterator < datalength:
            FA_dict[namearray[index]] += unpack('>h', FA_file.read(2))
            iterator += 1
        FA_file.seek(FA_str.find(carray[index]) + 20, 0)
        FA_dict[wavestr[index]] = int.from_bytes(FA_file.read(2), 'big')
        FA_file.seek(FA_str.find(darray[index]) + 16, 0)
        slen = int.from_bytes(FA_file.read(4), 'big') - 1
        FA_file.seek(int.from_bytes(FA_file.read(4), 'big') + 1, 0)
        FA_dict[dyestr[index]] = FA_file.read(slen)
        index += 1
