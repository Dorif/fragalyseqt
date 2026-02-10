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

from charset_normalizer import from_bytes
from numpy import (array, searchsorted, clip, where, zeros, full,
                   ones, column_stack)
from numpy.linalg import lstsq
from scipy.optimize import least_squares


def _decode_bytes(data):
    result = from_bytes(data).best()
    if result is not None:
        return str(result)
    return data.decode('UTF-8', errors='replace')


def set_graph_name(fdata):
    size_std = "Unknown size standard, "
    equipment = "Unknown equipment"
    graph_name = ""
    k_arr = fdata.keys()
    if "StdF1" in k_arr and fdata["StdF1"] != b'':
        size_std = _decode_bytes(fdata["StdF1"]) + " size standard, "
    if "DySN1" not in k_arr and "MODF1" not in k_arr and fdata["MODL1"] != b'310 ':
        # RapidHIT ID v1.X *.FSA files lack DySN1 and MODF1 keys,
        # because there are only one dye set and only one run module.
        equipment = "RapidHIT ID v1.X"
    elif ("RunN1" in k_arr and "HCFG3" in k_arr and fdata["HCFG3"] == b'3130xl' and
          chk_key_valid("DySN1", fdata) and ((b'.avt' in fdata["RunN1"]) or
          (b'\xd1\xca' in fdata["DySN1"]))):
        equipment = "Nanophore-05"
    elif chk_key_valid("MODL1", fdata) and fdata["MODL1"] == b'3200':
        equipment = "SeqStudio"
    elif ("NLNE1" in k_arr and "DyeW1" in k_arr and
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
    elif chk_key_valid("MCHN1", fdata) and "Promega" in _decode_bytes(
                                                            fdata["MCHN1"]):
        equipment = _decode_bytes(fdata["MCHN1"])
    elif chk_key_valid("HCFG3", fdata):
        equipment = _decode_bytes(fdata["HCFG3"])
    elif chk_key_valid("MODL1", fdata):
        equipment = _decode_bytes(fdata["MODL1"])
    if chk_key_valid("SpNm1", fdata):
        graph_name = _decode_bytes(fdata["SpNm1"]) + ", "
    elif chk_key_valid("CTNM1", fdata):
        graph_name = _decode_bytes(fdata["CTNM1"]) + ", "
    return (graph_name + size_std + equipment)


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
            darr.append(_decode_bytes(fdata[dname[i]]))
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
    if 'ROX' in ILS or 'CXR' in ILS:
        return fdata["DATA4"]
    elif 'CC0' in ILS:
        return fdata["DATA108"]
    else:
        return fdata["DATA105"]


def southern_m0(L1, m1, L2, m2, L3, m3):
    denom_L = L2 - L3
    denom_m = m2 - m1
    if denom_L == 0.0 or denom_m == 0.0:
        return 0.0
    A = (L1 - L2) / denom_L * (m3 - m2) / denom_m
    if A == 1.0:
        return 0.0
    return (m3 - m1 * A) / (1.0 - A)


def _southern_3pt_size(L1, m1, L2, m2, L3, m3, m_query):
    m0 = southern_m0(L1, m1, L2, m2, L3, m3)
    inv1 = 1.0 / (m1 - m0)
    inv2 = 1.0 / (m2 - m0)
    if inv1 == inv2:
        return L1
    c = (L1 - L2) / (inv1 - inv2)
    L0 = L1 - c * inv1
    return c / (m_query - m0) + L0


def _vec_southern_3pt(size_std, ladder_peaks, i, m_query):
    L1, m1 = size_std[i], ladder_peaks[i]
    L2, m2 = size_std[i + 1], ladder_peaks[i + 1]
    L3, m3 = size_std[i + 2], ladder_peaks[i + 2]
    denom_L = L2 - L3
    denom_m = m2 - m1
    zero_denom = (denom_L == 0.0) | (denom_m == 0.0)
    safe_dL = where(zero_denom, 1.0, denom_L)
    safe_dm = where(zero_denom, 1.0, denom_m)
    A = (L1 - L2) / safe_dL * (m3 - m2) / safe_dm
    bad_m0 = zero_denom | (A == 1.0)
    m0 = where(bad_m0, 0.0, (m3 - m1 * A) / where(bad_m0, 1.0, 1.0 - A))
    inv1 = 1.0 / (m1 - m0)
    inv2 = 1.0 / (m2 - m0)
    same_inv = inv1 == inv2
    c = where(same_inv, 0.0, (L1 - L2) / where(same_inv, 1.0, inv1 - inv2))
    L0 = L1 - c * inv1
    return where(same_inv, L1, c / (m_query - m0) + L0)


def southern_fit_local(ladder_peaks, size_std, query_points):
    ladder_peaks = array(ladder_peaks, dtype=float)
    size_std = array(size_std, dtype=float)
    query_points = array(query_points, dtype=float)
    n = len(ladder_peaks)
    idx = searchsorted(ladder_peaks, query_points).astype(int)
    mask_a = (idx >= 2) & (idx < n)
    ia = clip(idx - 2, 0, n - 3)
    mask_b = (idx >= 1) & (idx + 1 < n)
    ib = clip(idx - 1, 0, n - 3)
    est_a = _vec_southern_3pt(size_std, ladder_peaks, ia, query_points)
    est_b = _vec_southern_3pt(size_std, ladder_peaks, ib, query_points)
    both = mask_a & mask_b
    result = where(both, (est_a + est_b) / 2.0,
                   where(mask_a, est_a, where(mask_b, est_b, 0.0)))
    mask_none = ~(mask_a | mask_b)
    if mask_none.any():
        low = mask_none & (query_points <= ladder_peaks[0])
        high = mask_none & ~low
        if low.any():
            result[low] = _vec_southern_3pt(
                size_std, ladder_peaks,
                zeros(low.sum(), dtype=int), query_points[low])
        if high.any():
            result[high] = _vec_southern_3pt(
                size_std, ladder_peaks,
                full(high.sum(), n - 3, dtype=int), query_points[high])
    return result


def southern_fit_global(ladder_peaks, size_std, query_points):
    ladder_peaks = array(ladder_peaks, dtype=float)
    size_std = array(size_std, dtype=float)
    query_points = array(query_points, dtype=float)
    n = len(ladder_peaks)
    mid = n // 2
    m0_init = southern_m0(size_std[0], ladder_peaks[0],
                          size_std[mid], ladder_peaks[mid],
                          size_std[-1], ladder_peaks[-1])
    inv_m = 1.0 / (ladder_peaks - m0_init)
    A_mat = column_stack([inv_m, ones(n)])
    sol = lstsq(A_mat, size_std, rcond=None)
    c_init, L0_init = sol[0]

    def residuals(params):
        c, m0, L0 = params
        return size_std - (c / (ladder_peaks - m0) + L0)

    res = least_squares(residuals, [c_init, m0_init, L0_init])
    c, m0, L0 = res.x
    return c / (query_points - m0) + L0


def chk_key_valid(key, fdata):
    return key in fdata and fdata[key] is not None
