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

"""Panel bin file parsers and allele binning for FragalyseQt.

Supports two panel formats:
  * GeneMapper  — tab-delimited *_Panels_*.txt (GMID-X v1.x and GM v3.0),
                  with an optional companion *_Bins_*.txt for precise allele
                  sizes and acceptance windows, and an optional
                  *_Stutter_*.txt for per-marker stutter ratios.
  * GeneMarker  — self-contained *.xml files that carry allele sizes,
                  bin widths, and stutter ratios in a single document.

Both parsers return the same unified internal structure so that the binning
engine (assign_alleles) is format-agnostic:

    {
        "<panel_name>": {
            "<marker_name>": {
                "dye":      str | None,   # GeneMapper colour word
                "min_size": float | None, # coarse range lower bound
                "max_size": float | None, # coarse range upper bound
                "stutter": {              # per-marker stutter thresholds
                    "minus": float | None,  # n-1 ratio (e.g. 0.13)
                    "plus":  float | None,  # n+1 ratio (e.g. 0.01)
                },
                "alleles": [
                    {
                        "label":     str,
                        "size":      float | None,  # None when no bins file
                        "left_bin":  float | None,
                        "right_bin": float | None,
                        "virtual":   bool,          # ladder / OL allele
                    },
                    ...
                ],
            },
            ...
        },
        ...
    }

One GeneMapper Panels file can contain several Panel sections (e.g.
SeqStudio_Panels_v7X.txt has seven); each becomes a separate top-level key.
"""

import os
from xml.etree.ElementTree import parse as _xmlparse

# ---------------------------------------------------------------------------
# Channel index → GeneMapper colour word mapping
# ---------------------------------------------------------------------------

# CE instruments always place dyes in the same channel order regardless of
# the dye chemistry used: channel 1 = blue, 2 = green, 3 = yellow, 4 = red,
# 5 = purple, 6 = orange.  Using the 1-based channel index is therefore more
# reliable than trying to map dye trade names (which vary across kits and
# instruments) to colour words.
_CHANNEL_INDEX_TO_COLOR = {
    1: 'blue', 2: 'green', 3: 'yellow',
    4: 'red',  5: 'orange', 6: 'purple',
}


# ---------------------------------------------------------------------------
# GeneMapper panel file — internal helpers
# ---------------------------------------------------------------------------

def _extract_allele_list(parts):
    """Locate the allele-list column in a GeneMapper panel marker row.

    Column layout (0-based):
      marker(0), dye(1), min_size(2), max_size(3), control_alleles(4),
      bit_precision(5), reserved(6),
      [optional: indel_flag(7), variant_flag(8), ...]  allele_list(last).

    The allele list is a comma-separated string (e.g. "12, 13, 14," or
    "A1, A2,").  Boolean/reserved tokens like 'false', 'true', 'none', '-'
    never contain commas, so scanning right-to-left for the first comma is
    a clean discriminator.  Columns 0–6 are intentionally skipped to avoid
    picking up the control_alleles field (col 4) which also uses commas.
    """
    for col in range(len(parts) - 1, 6, -1):
        if ',' in parts[col]:
            return parts[col]
    return ''


def _parse_genemapper_panels(path):
    """Parse a GeneMapper *_Panels_* / *_Panel_* text file.

    Handles both GMID-X (v1.0–1.6) and GM v3.0 variants, and multi-panel
    files (multiple "Panel  <name>  null" sections).

    Returns dict[panel_name -> dict[marker_name -> marker_entry]].
    Allele sizes are None until enriched by _parse_genemapper_bins().
    Stutter thresholds are None until enriched by _parse_genemapper_stutter().
    """
    panels = {}
    current_panel = None

    with open(path, encoding='utf-8', errors='replace') as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            if not line.strip():
                continue
            parts = line.split('\t')
            key = parts[0].strip().lower()

            # Skip metadata header rows common to all GeneMapper files
            if key in ('version', 'kit type:', 'chemistry kit', 'binset name'):
                continue

            # Start of a new panel section
            if key == 'panel':
                panel_name = parts[1].strip() if len(parts) > 1 else 'Unknown'
                current_panel = {}
                panels[panel_name] = current_panel
                continue

            # Marker data rows need at least: name, dye, min_size, max_size
            if current_panel is None or len(parts) < 4:
                continue

            try:
                min_size = float(parts[2])
                max_size = float(parts[3])
            except ValueError:
                continue  # Non-numeric sizes → column header or junk line

            marker_name = parts[0].strip()
            dye = parts[1].strip().lower()

            allele_str = _extract_allele_list(parts)
            allele_labels = ([a.strip() for a in allele_str.split(',')
                              if a.strip()]
                             if allele_str else [])

            current_panel[marker_name] = {
                'dye': dye,
                'min_size': min_size,
                'max_size': max_size,
                'stutter': {'minus': None, 'plus': None},
                # Allele sizes left as None; filled in if a Bins file is loaded
                'alleles': [{'label': lbl, 'size': None,
                             'left_bin': None, 'right_bin': None,
                             'virtual': False}
                            for lbl in allele_labels],
            }

    return panels


def _parse_genemapper_bins(path):
    """Parse a GeneMapper *_Bins_* text file (single or multi-panel).

    Bins files may contain multiple 'Panel Name' sections, each covering one
    kit panel.  Returns dict[panel_name -> dict[marker_name -> allele list]].
    """
    all_bins = {}
    current_panel_bins = None
    current_marker = None

    with open(path, encoding='utf-8', errors='replace') as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            if not line.strip():
                continue
            parts = line.split('\t')
            key = parts[0].strip().lower()

            if key == 'panel name':
                panel_name = parts[1].strip() if len(parts) > 1 else None
                if panel_name:
                    current_panel_bins = {}
                    all_bins[panel_name] = current_panel_bins
                current_marker = None
                continue
            if key == 'marker name':
                current_marker = parts[1].strip() if len(parts) > 1 else None
                if current_marker and current_panel_bins is not None:
                    current_panel_bins[current_marker] = []
                continue
            if key in ('version', 'chemistry kit', 'binset name'):
                continue

            # Allele row: label, size, left_bin, right_bin [, 'virtual']
            if current_marker is None or current_panel_bins is None or len(parts) < 4:
                continue
            try:
                label = parts[0].strip()
                size = float(parts[1])
                left_bin = float(parts[2])
                right_bin = float(parts[3])
            except ValueError:
                continue
            virtual = len(parts) > 4 and parts[4].strip().lower() == 'virtual'
            current_panel_bins[current_marker].append({
                'label': label, 'size': size,
                'left_bin': left_bin, 'right_bin': right_bin,
                'virtual': virtual,
            })

    return all_bins


def _parse_genemapper_stutter(path):
    """Parse a GeneMapper *_Stutter_*.txt file.

    Stutter files use a multi-panel / multi-marker layout:

        Panel Name  <name>
        Marker Name <name>
        <ratio>  <window_lo>  <window_hi>  <type>  ...

    *type* is "Minus" (n-1) or "Plus" (n+1).  A marker may have several
    rows; entries with window_hi < 3.0 bp represent half-repeat or n-2
    stutters and are skipped — only standard single-unit stutter
    (tetranucleotide: ~3.25–4.75, trinucleotide: ~2.25–3.75) is imported.
    When multiple qualifying rows exist for the same type, the largest ratio
    is kept.

    Returns dict[panel_name -> dict[marker_name ->
                                    {"minus": float|None, "plus": float|None}]]
    """
    all_stutter = {}
    current_panel_stutter = None
    current_marker = None

    with open(path, encoding='utf-8', errors='replace') as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            if not line.strip():
                continue
            parts = line.split('\t')
            key = parts[0].strip().lower()

            if key in ('version', 'chemistry kit'):
                continue
            if key == 'panel name':
                panel_name = parts[1].strip() if len(parts) > 1 else None
                if panel_name:
                    current_panel_stutter = {}
                    all_stutter[panel_name] = current_panel_stutter
                current_marker = None
                continue
            if key == 'marker name':
                current_marker = parts[1].strip() if len(parts) > 1 else None
                if current_marker and current_panel_stutter is not None:
                    if current_marker not in current_panel_stutter:
                        current_panel_stutter[current_marker] = {
                            'minus': None, 'plus': None}
                continue

            # Stutter data row: ratio  window_lo  window_hi  type
            if current_marker is None or current_panel_stutter is None:
                continue
            if len(parts) < 4:
                continue
            try:
                ratio = float(parts[0])
                window_hi = float(parts[2])
                stype = parts[3].strip().lower()
            except ValueError:
                continue

            # Skip sub-unit stutter (window_hi < 3.0 bp = half-repeat or n-2)
            if window_hi < 3.0:
                continue

            entry = current_panel_stutter.setdefault(
                current_marker, {'minus': None, 'plus': None})
            if stype == 'minus':
                if entry['minus'] is None or ratio > entry['minus']:
                    entry['minus'] = ratio
            elif stype == 'plus':
                if entry['plus'] is None or ratio > entry['plus']:
                    entry['plus'] = ratio

    return all_stutter


# ---------------------------------------------------------------------------
# Public GeneMapper entry point
# ---------------------------------------------------------------------------

def parse_genemapper(panels_path, bins_path='', stutter_path=''):
    """Load a GeneMapper Panels file, optionally enriched by Bins and Stutter files.

    Parameters
    ----------
    panels_path : str
        Path to a ``*_Panel*.txt`` or ``*_Panel_*.txt`` file.
    bins_path : str
        Path to a companion ``*_Bins*.txt`` file.  Pass an empty string
        (default) to load panels only without bin data.
    stutter_path : str
        Path to a companion ``*_Stutter*.txt`` file.  Pass an empty string
        (default) to skip stutter data import.

    Returns
    -------
    dict[str, dict]
        Unified panel dict.  Multi-panel files yield multiple top-level keys.
        When no bins file is given, allele sizes stay None and only coarse
        marker-range binning is available (see assign_alleles).
        When no stutter file is given, stutter thresholds stay None and the
        stutter filter falls back to its global spinbox values.
    """
    panels = _parse_genemapper_panels(panels_path)

    if bins_path and os.path.isfile(bins_path):
        all_bins = _parse_genemapper_bins(bins_path)
        if all_bins:
            for bins_panel_name, bins_data in all_bins.items():
                if bins_panel_name in panels:
                    target_panel = panels[bins_panel_name]
                elif len(panels) == 1:
                    target_panel = next(iter(panels.values()))
                else:
                    continue
                for marker_name, allele_list in bins_data.items():
                    target_panel[marker_name]['alleles'] = allele_list

    if stutter_path and os.path.isfile(stutter_path):
        all_stutter = _parse_genemapper_stutter(stutter_path)
        if all_stutter:
            for stutter_panel_name, stutter_data in all_stutter.items():
                if stutter_panel_name in panels:
                    target_panel = panels[stutter_panel_name]
                elif len(panels) == 1:
                    target_panel = next(iter(panels.values()))
                else:
                    continue
                for marker_name, stutter_entry in stutter_data.items():
                    target_panel[marker_name]['stutter'] = stutter_entry

    return panels


# ---------------------------------------------------------------------------
# GeneMarker XML parser
# ---------------------------------------------------------------------------

# DyeIndex values in GeneMarker XML correspond to the same colour words
# used by GeneMapper panel files.
_GENEMARKER_DYE_INDEX = {
    '0': None,
    '1': 'blue', '2': 'green', '3': 'yellow',
    '4': 'red',  '5': 'purple', '6': 'orange',
}


def parse_genemarker(xml_path):
    """Parse a GeneMarker XML panel file.

    GeneMarker files are self-contained: they carry allele sizes, acceptance
    windows, and stutter ratios, so no companion files are needed.

    Stutter thresholds are read from each locus's ``<LocusFilter>`` element:
      * ``StutterPer_N_L4`` + ``DecimalStutterPer_N_L4`` → n-1 (minus) ratio
      * ``StutterPer_N_R4`` + ``DecimalStutterPer_N_R4`` → n+1 (plus) ratio

    The combined percentage is computed as
    ``(integer_part + decimal_part / 10) / 100``
    and stored as a fraction in the unified ``"stutter"`` dict.

    ``Control='1'`` alleles are allelic ladder reference peaks; they are
    flagged as virtual so the analyst can identify them in the output.

    Returns the same unified dict as parse_genemapper().
    """
    tree = _xmlparse(xml_path)
    root = tree.getroot()
    panel_name = (root.findtext('PanelName')
                  or os.path.splitext(os.path.basename(xml_path))[0])

    markers = {}
    loci_node = root.find('Loci')
    if loci_node is None:
        return {panel_name: markers}

    for locus in loci_node.findall('Locus'):
        marker_name = (locus.findtext('MarkerTitle') or '').strip()
        if not marker_name:
            continue

        dye_index = (locus.findtext('DyeIndex') or '1').strip()
        dye = _GENEMARKER_DYE_INDEX.get(dye_index, 'blue')

        try:
            min_size = float(locus.findtext('LowerBoundary') or 0)
            max_size = float(locus.findtext('UpperBoundary') or 0)
        except ValueError:
            min_size = max_size = 0.0

        # Stutter thresholds from <LocusFilter> element
        stutter = {'minus': None, 'plus': None}
        lf = locus.find('LocusFilter')
        if lf is not None:
            def _stutter_ratio(per_attr, dec_attr):
                try:
                    per = float(lf.get(per_attr, '0'))
                    dec = float(lf.get(dec_attr, '0'))
                    ratio = (per + dec / 10.0) / 100.0
                    return ratio if ratio > 0 else None
                except (ValueError, TypeError):
                    return None
            stutter['minus'] = _stutter_ratio(
                'StutterPer_N_L4', 'DecimalStutterPer_N_L4')
            stutter['plus'] = _stutter_ratio(
                'StutterPer_N_R4', 'DecimalStutterPer_N_R4')

        alleles = []
        for allele_el in locus.findall('Allele'):
            label = allele_el.get('Label', '').strip()
            # GeneMarker stores 'DefSize' (theoretical) and 'Size' (measured
            # from an allelic ladder run).  Prefer the measured value.
            raw_size = (allele_el.get('Size') or allele_el.get('DefSize')
                        or '0')
            try:
                size = float(raw_size)
                left_bin = float(allele_el.get('Left_Binning', '0.5'))
                right_bin = float(allele_el.get('Right_Binning', '0.5'))
            except ValueError:
                continue
            virtual = allele_el.get('Control', '0') == '1'
            alleles.append({
                'label': label, 'size': size,
                'left_bin': left_bin, 'right_bin': right_bin,
                'virtual': virtual,
            })

        markers[marker_name] = {
            'dye': dye, 'min_size': min_size, 'max_size': max_size,
            'stutter': stutter,
            'alleles': alleles,
        }

    return {panel_name: markers}


# ---------------------------------------------------------------------------
# OSIRIS LadderInfo XML parser
# ---------------------------------------------------------------------------

def _xml_root_tag(path):
    """Return the root element tag of an XML file, or None on any error."""
    try:
        return _xmlparse(path).getroot().tag
    except Exception:
        return None


def parse_osiris(xml_path, default_bin=1.5):
    """Parse an OSIRIS LadderInfo XML file (NIST OSIRIS MarkerSet schema).

    Both v2.0 (MarkerSet.xsd) and v2.7 (MarkerSetV4.xsd) are handled
    identically — ILS search-region fields are not used; only marker names,
    channel/colour mapping, size ranges, and ladder allele BPs are extracted.

    Allele bin widths default to ±default_bin bp.  OSIRIS stores theoretical
    integer BP values that can differ from instrument-measured sizes by
    several bp, so the default is wider than the ±0.5 used for GeneMapper.

    Stutter thresholds are not stored in this format and are left as None.

    Returns the same unified dict as parse_genemapper() / parse_genemarker().
    Multiple <Set> elements per file each become a separate top-level key.
    """
    tree = _xmlparse(xml_path)
    root = tree.getroot()
    result = {}

    for kit_set in root.findall('.//Kits/Set'):
        panel_name = (kit_set.findtext('Name') or '').strip()
        if not panel_name:
            continue

        # Kit channel number → lowercase colour word (blue/green/yellow/…)
        ch_to_color = {}
        for ch_el in kit_set.findall('FsaChannelMap/Channel'):
            try:
                num = int(ch_el.findtext('KitChannelNumber') or '0')
                color = (ch_el.findtext('Color') or '').strip().lower()
                if num and color:
                    ch_to_color[num] = color
            except ValueError:
                pass

        markers = {}
        for locus in kit_set.findall('Locus'):
            marker_name = (locus.findtext('Name') or '').strip()
            if not marker_name:
                continue
            try:
                ch_num = int(locus.findtext('Channel') or '0')
                min_bp = float(locus.findtext('MinBP') or '0')
                max_bp = float(locus.findtext('MaxBP') or '0')
            except ValueError:
                continue

            alleles = []
            for a in locus.findall('LadderAlleles/Allele'):
                label = (a.findtext('Name') or '').strip()
                bp_txt = a.findtext('BP')
                if not label or bp_txt is None:
                    continue
                try:
                    size = float(bp_txt)
                except ValueError:
                    continue
                alleles.append({
                    'label': label, 'size': size,
                    'left_bin': default_bin, 'right_bin': default_bin,
                    'virtual': False,
                })

            markers[marker_name] = {
                'dye': ch_to_color.get(ch_num),
                'min_size': min_bp, 'max_size': max_bp,
                'stutter': {'minus': None, 'plus': None},
                'alleles': alleles,
            }

        if markers:
            result[panel_name] = markers

    return result


# ---------------------------------------------------------------------------
# Auto-detect loader
# ---------------------------------------------------------------------------

def has_bin_data(panel_data):
    """Return True if at least one allele in the panel has a precise size.

    Used to detect whether a companion Bins file was successfully loaded for
    a GeneMapper panel, so the UI can offer manual selection when auto-
    detection fails.
    """
    return any(a['size'] is not None
               for panel in panel_data.values()
               for marker in panel.values()
               for a in marker['alleles'])


def has_stutter_data(panel_data):
    """Return True if at least one marker in the panel has a stutter ratio."""
    return any(marker.get('stutter', {}).get('minus') is not None
               or marker.get('stutter', {}).get('plus') is not None
               for panel in panel_data.values()
               for marker in panel.values())


def load_panel(path):
    """Detect format from file extension and delegate to the right parser.

    ``.xml`` → GeneMarker  (parse_genemarker)
    ``.txt`` → GeneMapper  (parse_genemapper, with auto-detected bins)

    Returns the unified panel dict.
    """
    if path.lower().endswith('.xml'):
        return parse_genemarker(path)
    return parse_genemapper(path)


# ---------------------------------------------------------------------------
# Allele binning engine
# ---------------------------------------------------------------------------

def assign_alleles(peak_sizes, peak_channel_indices, panel_markers):
    """Assign allele labels to sized peaks using one panel's marker data.

    Algorithm (per peak):
    1. Map the peak's 1-based channel index to a colour word via
       _CHANNEL_INDEX_TO_COLOR (1=blue, 2=green, 3=yellow, 4=red,
       5=purple, 6=orange).  Using channel index rather than dye name is
       instrument-agnostic: dye trade names vary across kits, but channel
       order is fixed by the CE instrument.
    2. Skip markers whose dye colour differs (when both are known).
    3. Skip markers whose [min_size, max_size] range does not contain the peak.
    4. If allele-level bin data exist, assign the first allele whose acceptance
       window [size − left_bin, size + right_bin] covers the peak and stop.
    5. If no bin data are present, annotate with the marker name and '?'.

    Parameters
    ----------
    peak_sizes           : sequence of float
    peak_channel_indices : sequence of int  — 1-based channel numbers
    panel_markers        : dict[marker_name -> marker_entry]  — one panel's data

    Returns
    -------
    list of str, same length as peak_sizes.
    Format: "MarkerName:Allele"  e.g. "D3S1358:14"
            "MarkerName:?"       for range-only hits (no bin data)
            "MarkerName:14*"     for virtual / allelic-ladder alleles
            "OL"                 for no match (out of ladder)
    """
    results = [''] * len(peak_sizes)

    for i, (size, ch_idx) in enumerate(zip(peak_sizes, peak_channel_indices)):
        peak_color = _CHANNEL_INDEX_TO_COLOR.get(int(ch_idx))

        for marker_name, marker in panel_markers.items():
            if peak_color and marker['dye'] and peak_color != marker['dye']:
                continue

            mn, mx = marker['min_size'], marker['max_size']
            if mn and mx and (size < mn or size > mx):
                continue

            has_bins = any(a['size'] is not None for a in marker['alleles'])

            if has_bins:
                for allele in marker['alleles']:
                    if allele['size'] is None:
                        continue
                    if (allele['size'] - allele['left_bin']
                            <= size <=
                            allele['size'] + allele['right_bin']):
                        suffix = '*' if allele['virtual'] else ''
                        results[i] = f"{marker_name}:{allele['label']}{suffix}"
                        break
                else:
                    # Dye and size range matched but no allele bin covers peak
                    results[i] = f"{marker_name}:OL"
            else:
                results[i] = f"{marker_name}:?"
            break  # stop after first matching marker

    return results
