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

from xml.etree.ElementTree import iterparse, parse

_UDATAC = ["DATA1", "DATA2", "DATA3", "DATA4",
           "DATA105", "DATA106", "DATA107", "DATA108"]

_STD_TITLE_MAP = {
    "Стандарт длины S550": b"GDZ_S550",
    "Стандарт длины S450": b"GDZ_S450",
    "Стандарт длины S400": b"GDZ_S400",
}


def parse_frf(filepath):
    """Parse a Nanophore-05 .frf (XML) file.
    Returns an abif_raw-compatible dict for use with fragalyseqt's data model.
    Uses two passes: a full parse for compact metadata, then iterparse for the
    large <Data> block so each <Point> element is freed after reading."""

    # Pass 1 — metadata (small, parse fully)
    root = parse(filepath).getroot()
    n_channels = 8

    sample   = root.findtext("SampleName") or ""
    title    = root.findtext("Title") or ""
    instr    = root.findtext("InstrumentName") or "Нанофор 5"
    std_name = root.findtext("SizeStandard/Title") or ""

    wl_el = root.find("DyesWavelength")
    wavelengths = [int(e.text) for e in wl_el] if wl_el is not None else []

    # Pass 2 — channel arrays via iterparse (frees each <Point> after reading)
    channels = [[] for _ in range(n_channels)]
    for _, elem in iterparse(filepath, events=("end",)):
        if elem.tag == "Point":
            data_el = elem.find("Data")
            if data_el is not None:
                for i, v in enumerate(data_el):
                    if i < n_channels:
                        channels[i].append(int(v.text))
            elem.clear()

    if not channels[0]:
        raise ValueError("No data points found in FRF file")

    # FRF stores raw ADC counts with a large hardware DC offset (~45 000–
    # 121 000 per channel). FSA export from the same instrument normalises
    # these to near-zero baseline. Subtract the per-channel minimum so the
    # data looks comparable and the plot Y-axis is not dominated by the offset.
    for i in range(n_channels):
        mn = min(channels[i])
        channels[i] = [v - mn for v in channels[i]]

    abif_raw = {
        "Dye#1":  n_channels,
        # Keys matched to the Nanophore-05 branch in set_graph_name
        "HCFG3":  b"3130xl",
        "DySN1":  b"\xd1\xca",
        "RunN1":  b"run.avt",
        # Sample / standard labels for graph title
        "SpNm1":  (title or sample).encode("utf-8"),
        "StdF1":  _STD_TITLE_MAP.get(std_name, std_name.encode("utf-8")),
    }

    for i in range(n_channels):
        abif_raw[_UDATAC[i]] = channels[i]
        # Dye names are not stored in FRF; use "Ch{N}" + emission wavelength
        abif_raw[f"DyeN{i+1}"] = f"Ch{i+1}".encode()
        abif_raw[f"DyeW{i+1}"] = wavelengths[i] if i < len(wavelengths) else 0

    return abif_raw
