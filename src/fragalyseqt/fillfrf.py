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

from numpy import array as np_array
from numpy.linalg import solve as np_solve, LinAlgError
from xml.etree.ElementTree import iterparse, parse
from .fillhid import UDATAC


def parse_frf(filepath):
    """Parse a Nanophore-05 .frf (XML) file.
    Returns an abif_raw-compatible dict for use with fragalyseqt's data model.
    Uses two passes: a full parse for compact metadata, then iterparse for the
    large <Data> block so each <Point> element is freed after reading."""

    # Pass 1 — metadata (small, parse fully)
    root = parse(filepath).getroot()

    sample   = root.findtext("SampleName") or ""
    title    = root.findtext("Title") or ""
    std_name = root.findtext("SizeStandard/Title") or ""

    wl_el = root.find("DyesWavelength")
    wavelengths = [int(e.text) for e in wl_el] if wl_el is not None else []

    # <StandardChannel> is 1-based index of the ILS/size-standard channel.
    std_ch_text = root.findtext("StandardChannel")
    std_channel = int(std_ch_text) if std_ch_text is not None else None
    # Top-level <Matrix> is the spectral crosstalk calibration matrix.
    # Each row is an <ArrayOfDouble>; diagonal values are 1.
    # <Parameters><UseMatrix> controls whether the operator enabled the
    # correction; respect that flag so we do not corrupt data that was
    # collected with matrix correction disabled.
    use_matrix = root.findtext("Parameters/UseMatrix") == "true"
    matrix_el = root.find("Matrix")
    spectral_matrix = None
    if use_matrix and matrix_el is not None:
        rows = [[float(v.text) for v in row_el] for row_el in matrix_el]
        if rows:
            spectral_matrix = np_array(rows)

    # Pass 2 — channel arrays via iterparse (frees each <Point> after reading).
    # Allocate the maximum (8) slots; actual count is resolved afterwards.
    channels = [[] for _ in range(8)]
    for _, elem in iterparse(filepath, events=("end",)):
        if elem.tag == "Point":
            data_el = elem.find("Data")
            if data_el is not None:
                for i, v in enumerate(data_el):
                    if i < 8:
                        channels[i].append(int(v.text))
            elem.clear()

    if not channels[0]:
        raise ValueError("No data points found in FRF file")

    # Determine actual channel count: prefer the wavelength list (authoritative
    # metadata); fall back to counting channels that received data.
    if wavelengths:
        n_channels = min(len(wavelengths), 8)
    else:
        n_channels = sum(1 for ch in channels if ch)

    # Apply spectral crosstalk correction: solve M·x = raw for each time point.
    # Use the n_channels×n_channels top-left submatrix in case the stored
    # matrix is larger than the actual channel count.
    if spectral_matrix is not None and spectral_matrix.shape[0] >= n_channels:
        m = spectral_matrix[:n_channels, :n_channels]
        try:
            raw = np_array([channels[i] for i in range(n_channels)], dtype=float)
            corrected = np_solve(m, raw)
            for i in range(n_channels):
                channels[i] = corrected[i].tolist()
        except LinAlgError:
            pass  # singular matrix — skip correction, use raw counts

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
        "StdF1":  std_name.encode("utf-8"),
    }

    # Store 1-based ILS channel index for set_ILS_channel routing
    if std_channel is not None:
        abif_raw["STDC1"] = std_channel

    for i in range(n_channels):
        abif_raw[UDATAC[i]] = channels[i]
        # Dye names are not stored in FRF; use "Ch{N}" + emission wavelength
        abif_raw[f"DyeN{i+1}"] = f"Ch{i+1}".encode()
        abif_raw[f"DyeW{i+1}"] = wavelengths[i] if i < len(wavelengths) else 0

    return abif_raw
