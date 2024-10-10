def set_graph_name(fdata):
    from charset_normalizer import from_bytes
    size_standard = "Unknown size standard, "
    equipment = "Unknown equipment"
    keysarr = fdata.keys()
    if "StdF1" in keysarr and fdata["StdF1"] != b'':
        size_standard = str(fdata["StdF1"], 'UTF-8') + " size standard, "
    if ("DySN1" and "MODF1") not in keysarr and fdata["MODL1"] != b'310 ':
        # RapidHIT ID v1.X *.FSA files lack DySN1 and MODF1 keys,
        # because there are only one dye set and only one run module.
        equipment = "RapidHIT ID v1.X"
    elif (("RunN1" and "DySN1" and "HCFG3") in keysarr and
          fdata["DySN1"] is not None and ((b'.avt' in fdata["RunN1"]) or
          (b'\xd1\xca' in fdata["DySN1"])) and fdata["HCFG3"] == b'3130xl'):
        equipment = "Nanophore-05"
    elif fdata["MODL1"] == b'3200':
        equipment = "SeqStudio"
    elif (("NLNE1" and "DyeW1") in keysarr and
          "HCFG3" not in keysarr and fdata["DyeW1"] == 0):
        equipment = "Superyears Honor "
        if "DATA108" in keysarr:
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
    elif "HCFG3" in keysarr and fdata["HCFG3"] is not None:
        equipment = str(fdata["HCFG3"], 'UTF-8')
    elif "MODL1" in keysarr and fdata["MODL1"] is not None:
        equipment = str(fdata["MODL1"], 'UTF-8')
    if "SpNm1" in keysarr:
        graph_name = str(from_bytes(fdata["SpNm1"]).best()) + ", "
    elif "CTNM1" in keysarr:
        graph_name = str(from_bytes(fdata["CTNM1"]).best()) + ", "
    return (graph_name + size_standard + equipment)


def set_dye_array(fdata):
    k_arr = fdata.keys()
    drange = range(fdata["Dye#1"])
    darr = []
    dname = ["DyeN1", "DyeN2", "DyeN3", "DyeN4", "DyeN5", "DyeN6", "DyeN7",
             "DyeN8"]
    dwave = ["DyeW1", "DyeW2", "DyeW3", "DyeW4", "DyeW5", "DyeW6", "DyeW7",
             "DyeW8"]
    if dname[0] not in k_arr or fdata[dname[0]] is None:
        tmpd = ["FAM", "VIC", "TAMRA", "ROX", "LIZ", "SID", "Channel 7",
                "Channel 8"]
        for i in drange:
            darr.append(tmpd[i])
    else:
        # Checking if no emission wavelengths values are present or
        # wavelengths if they are equal to 0. Assuming if DyeW1 is
        # present and nonzero, others are present and non-zero too.
        for i in drange:
            darr.append(str(fdata[dname[i]], 'UTF-8'))
        if dwave[0] in k_arr and (fdata[dwave[0]] is not None and
                                  fdata[dwave[0]] != 0):
            for i in drange:
                darr[i] += (" " + str(fdata[dwave[i]]) + " nm")
    return darr


def set_IA_data(fdata):
    if "Peak1" in fdata.keys() and fdata["Peak1"] is not None:
        return True
    else:
        return False


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