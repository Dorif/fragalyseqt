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

def set_graph_name(fdata):
    from charset_normalizer import from_bytes
    size_std = "Unknown size standard, "
    equipment = "Unknown equipment"
    graph_name = ""
    k_arr = fdata.keys()
    if "StdF1" in k_arr and fdata["StdF1"] != b'':
        size_std = str(from_bytes(fdata["StdF1"]).best()) + " size standard, "
    if ("DySN1" and "MODF1") not in k_arr and fdata["MODL1"] != b'310 ':
        # RapidHIT ID v1.X *.FSA files lack DySN1 and MODF1 keys,
        # because there are only one dye set and only one run module.
        equipment = "RapidHIT ID v1.X"
    elif (("RunN1" and "HCFG3") in k_arr and fdata["HCFG3"] == b'3130xl' and
          chk_key_valid("DySN1", fdata) and ((b'.avt' in fdata["RunN1"]) or
          (b'\xd1\xca' in fdata["DySN1"]))):
        equipment = "Nanophore-05"
    elif chk_key_valid("MODL1", fdata) and fdata["MODL1"] == b'3200':
        equipment = "SeqStudio"
    elif (("NLNE1" and "DyeW1") in k_arr and
          "HCFG3" not in k_arr and fdata["DyeW1"] == 0):
        equipment = "Superyears Honor "
        if "DATA108" in k_arr:
            if fdata["NLNE1"] == 16:
                equipment += "1816"
            else:
                equipment += "1824"
        else:
            if fdata["NLNE1"] == 16:
                equipment += "1616"
            elif fdata["NLNE1"] == 24:
                equipment += "1624"
            else:
                equipment += "1696"
    elif chk_key_valid("HCFG3", fdata):
        equipment = str(fdata["HCFG3"], 'UTF-8')
    elif chk_key_valid("MODL1", fdata):
        equipment = str(fdata["MODL1"], 'UTF-8')
    if chk_key_valid("SpNm1", fdata):
        graph_name = str(from_bytes(fdata["SpNm1"]).best()) + ", "
    elif chk_key_valid("CTNM1", fdata):
        graph_name = str(from_bytes(fdata["CTNM1"]).best()) + ", "
    return (graph_name + size_standard + equipment)


def set_dye_array(fdata):
    drange = range(fdata["Dye#1"])
    darr = []
    dname = ["DyeN1", "DyeN2", "DyeN3", "DyeN4", "DyeN5", "DyeN6", "DyeN7",
             "DyeN8"]
    dwave = ["DyeW1", "DyeW2", "DyeW3", "DyeW4", "DyeW5", "DyeW6", "DyeW7",
             "DyeW8"]
    if not chk_key_valid(dname[0], fdata):
        tmpd = ["FAM", "VIC", "TAMRA", "ROX", "LIZ", "SID", "Channel 7",
                "Channel 8"]
        for i in drange:
            darr.append(tmpd[i])
    else:
        for i in drange:
            darr.append(str(fdata[dname[i]], 'UTF-8'))
        # Checking if no emission wavelengths values are present or
        # their wavelengths are equal to 0. Assuming if DyeW1 is
        # present and nonzero, others are present and nonzero too.
        if chk_key_valid(dwave[0], fdata) and fdata[dwave[0]] != 0:
            for i in drange:
                darr[i] += (" " + str(fdata[dwave[i]]) + " nm")
    return darr


def set_spl_dgr(alg):
    alg = alg.lower()
    if '5' in alg:
        return 5
    elif 'linear' in alg:
        return 1
    else:
        return 3


def set_knots(alg, arr, spldeg):
    alg = alg.lower()
    if 'weighted' in alg:
        s_len = len(arr)
        k1 = s_len//6
        k2 = 5*s_len//6
        if spldeg == 3:
            # Making LSQ weighted cubic spline work with
            # GS120LIZ ladder too.
            k1 += 1
        elif spldeg == 5:
            # Making 5th degree LSQ weighted spline work with
            # GS120LIZ ladder too.
            k1 += 3
        return arr[k1:k2]
    else:
        return None


def set_lsq_ord(alg):
    if '2' in alg:
        return 2
    elif '3' in alg:
        return 3
    else:
        return 5


def set_ILS_channel(fdata, ILS):
    ILS = ILS.upper()
    if ('ROX' or 'CXR') in ILS:
        return fdata["DATA4"]
    elif 'CC0' in ILS:
        return fdata["DATA108"]
    else:
        return fdata["DATA105"]


def chk_key_valid(key, fdata):
    if key in fdata.keys() and fdata[key] is not None:
        return True
    else:
        return False
