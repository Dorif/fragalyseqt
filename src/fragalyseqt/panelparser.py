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

# Panel bin file parsers and allele binning for FragalyseQt.

# Supports two panel formats:
#   * GeneMapper  — tab-delimited *_Panels_*.txt with an optional companion
#                   *_Bins_*.txt for allele sizes and acceptance windows, and
#                   an optional *_Stutter_*.txt for per-marker stutter ratios.
#   * GeneMarker  — self-contained *.xml files that carry allele sizes, bin
#                   widths, and stutter ratios in a single document.

# Both parsers return the same unified internal structure so that the binning
# engine (assign_alleles) is format-agnostic:

#     {
#         "<panel_name>": {
#             "<marker_name>": {
#                 "dye": str | None, # GeneMapper colour word
#                 "min_size": float | None, # coarse range lower bound
#                 "max_size": float | None, # coarse range upper bound
#                 "stutter": { # per-marker stutter thresholds
#                     "minus": float | None, # n-1 ratio (e.g. 0.13)
#                     "plus": float | None, # n+1 ratio (e.g. 0.01)
#                 },
#                 "alleles": [
#                     {
#                         "label": str,
#                         "size": float | None, # None when no bins file
#                         "left_bin": float | None,
#                         "right_bin": float | None,
#                         "virtual": bool, # ladder / OL allele
#                     },
#                     ...
#                 ],
#             },
#             ...
#         },
#         ...
#     }

# One GeneMapper Panels file can contain several Panel sections; each becomes a
# separate top-level key.

from os import makedirs
from os.path import isfile, splitext, basename, dirname
from json import load as json_load, dump as json_dump
from .safexml import parse as _xmlparse
from .setvar import CHANNEL_COLOR

_LIBRARY_VERSION = 1


# ---------------------------------------------------------------------------
# GeneMapper panel file — internal helpers
# ---------------------------------------------------------------------------

def _extract_allele_list(parts):
    # Locate the allele-list column in a GeneMapper panel marker row.
    # Column layout (0-based):
    # marker(0), dye(1), min_size(2), max_size(3), control_alleles(4),
    # bit_precision(5), reserved(6),
    # [optional: indel_flag(7), variant_flag(8), ...]  allele_list(last).

    # The allele list is a comma-separated string (e.g. "12, 13, 14," or "A1,
    # A2,"). Boolean/reserved tokens like 'false', 'true', 'none', '-' never
    # contain commas, so scanning right-to-left for the first comma is a clean
    # discriminator. Columns 0–6 are intentionally skipped to avoid picking up
    # the control_alleles field (col 4) which also uses commas.
    for col in range(len(parts) - 1, 6, -1):
        if ',' in parts[col]:
            return parts[col]
    return ''


def _parse_genemapper_panels(path):
    # Parse a GeneMapper *_Panels_* / *_Panel_* text file.
    # Handles both single-panel and multi-panel files (multiple
    # "Panel name null"sections).
    # Returns dict[panel_name -> dict[marker_name -> marker_entry]].
    # Allele sizes are None until enriched by _parse_genemapper_bins().
    # Stutter thresholds are None until enriched by _parse_genemapper_stutter().
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
                # Non-numeric sizes → column header or junk line
                continue
            marker_name = parts[0].strip()
            dye = parts[1].strip().lower()
            allele_str = _extract_allele_list(parts)
            allele_labels = ([a.strip() for a in allele_str.split(',')
                              if a.strip()]
                             if allele_str else [])
            current_panel[marker_name] = {
                'dye': dye, 'min_size': min_size, 'max_size': max_size,
                'stutter': {'minus': None, 'plus': None},
                # Allele sizes left as None; filled in if a Bins file is loaded
                'alleles': [{'label': lbl, 'size': None, 'left_bin': None,
                             'right_bin': None, 'virtual': False}
                            for lbl in allele_labels],
            }
    return panels


def _parse_genemapper_bins(path):
    # Parse a GeneMapper *_Bins_* text file (single or multi-panel).
    # Bins files may contain multiple 'Panel Name' sections, each covering one
    # kit panel. Returns dict[panel_name -> dict[marker_name -> allele list]].
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
                current_marker = None
                if len(parts) < 2 or not parts[1].strip():
                    continue
                panel_name = parts[1].strip()
                current_panel_bins = {}
                all_bins[panel_name] = current_panel_bins
                continue
            if key == 'marker name':
                current_marker = None
                if len(parts) > 1 and parts[1].strip():
                    current_marker = parts[1].strip()
                    if current_panel_bins is not None:
                        current_panel_bins[current_marker] = []
                continue
            if key in ('version', 'chemistry kit', 'binset name'):
                continue
            # Allele row: label, size, left_bin, right_bin [, 'virtual']
            if current_marker is None or current_panel_bins is None:
                continue
            if len(parts) < 4:
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
                'label': label, 'size': size, 'left_bin': left_bin,
                'right_bin': right_bin, 'virtual': virtual,
            })
    return all_bins


def _parse_genemapper_stutter(path):
    # Parse a GeneMapper *_Stutter_*.txt file.
    # Stutter files use a multi-panel / multi-marker layout:
    #     Panel Name  <name>
    #     Marker Name <name>
    #     <ratio>  <window_lo>  <window_hi>  <type>  ...
    # *type* is "Minus" (n-1) or "Plus" (n+1).  A marker may have several
    # rows for the same type; the largest ratio is kept.
    # Returns dict[panel_name -> dict[marker_name ->
    #                                 {"minus": float|None, "plus": float|None}]]
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
                current_marker = None
                if len(parts) < 2 or not parts[1].strip():
                    continue
                panel_name = parts[1].strip()
                current_panel_stutter = {}
                all_stutter[panel_name] = current_panel_stutter
                continue
            if key == 'marker name':
                current_marker = None
                if len(parts) > 1 and parts[1].strip():
                    current_marker = parts[1].strip()
                    if current_panel_stutter is not None:
                        current_panel_stutter.setdefault(
                            current_marker, {'minus': None, 'plus': None})
                continue
            # Stutter data row: ratio  window_lo  window_hi  type
            if current_marker is None or current_panel_stutter is None:
                continue
            if len(parts) < 4:
                continue
            try:
                ratio = float(parts[0])
                stype = parts[3].strip().lower()
            except ValueError:
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

def _apply_companion(path, parser, panels, key):
    # Enrich *panels* in-place from a companion Bins or Stutter file.
    if not path or not isfile(path):
        return
    data = parser(path)
    if not data:
        return
    for panel_name, marker_data in data.items():
        if panel_name in panels:
            target = panels[panel_name]
        elif len(panels) == 1:
            target = next(iter(panels.values()))
        else:
            continue
        for marker_name, value in marker_data.items():
            target[marker_name][key] = value


def parse_genemapper(panels_path, bins_path='', stutter_path=''):
    # Load a GeneMapper Panels file, optionally with Bins and Stutter files.

    # Parameters
    # ----------
    # panels_path : str
    # Path to a ``*_Panel*.txt`` or ``*_Panel_*.txt`` file.
    # bins_path : str
    # Path to a companion ``*_Bins*.txt`` file. Pass an empty string to load
    # panels only.
    # stutter_path : str
    # Path to a companion ``*_Stutter*.txt`` file. Pass an empty string to skip
    # stutter data import.

    # Returns
    # -------
    # dict[str, dict]
    # Unified panel dict. Multi-panel files yield multiple top-level keys. Without
    # bins file, allele sizes are None and only coarse marker-range binning is
    # available (see assign_alleles). Without stutter file, stutter thresholds are
    # None and no stutter filtering is applied for the affected markers.
    panels = _parse_genemapper_panels(panels_path)
    _apply_companion(bins_path, _parse_genemapper_bins, panels, 'alleles')
    _apply_companion(stutter_path, _parse_genemapper_stutter, panels, 'stutter')
    return panels


# ---------------------------------------------------------------------------
# GeneMarker XML parser
# ---------------------------------------------------------------------------


def parse_genemarker(xml_path):
    # Parse a GeneMarker XML panel file.

    # GeneMarker files are self-contained: they carry allele sizes, acceptance
    # windows, and stutter ratios, so no companion files are needed.

    # Stutter thresholds are read from each locus's ``<LocusFilter>`` element:
    #   * ``StutterPer_N_L4`` + ``DecimalStutterPer_N_L4`` → n-1 (minus) ratio
    #   * ``StutterPer_N_R4`` + ``DecimalStutterPer_N_R4`` → n+1 (plus) ratio

    # The combined percentage is computed as
    # ``(integer_part + decimal_part / 10) / 100``
    # and stored as a fraction in the unified ``"stutter"`` dict.

    # ``Control='1'`` alleles are allelic ladder reference peaks; they are
    # flagged as virtual so the analyst can identify them in the output.

    # Returns the same unified dict as parse_genemapper().

    tree = _xmlparse(xml_path)
    root = tree.getroot()
    panel_name = (root.findtext('PanelName')
                  or splitext(basename(xml_path))[0])
    markers = {}
    loci_node = root.find('Loci')
    if loci_node is None:
        return {panel_name: markers}
    for locus in loci_node.findall('Locus'):
        marker_name = locus.findtext('MarkerTitle').strip()
        dye_index = locus.findtext('DyeIndex').strip()
        dye = CHANNEL_COLOR.get(int(dye_index))
        min_size = float(locus.findtext('LowerBoundary'))
        max_size = float(locus.findtext('UpperBoundary'))
        # Stutter thresholds from <LocusFilter> element
        stutter = {'minus': None, 'plus': None}
        lf = locus.find('LocusFilter')
        if lf is not None:

            def _stutter_ratio(per_attr, dec_attr):
                try:
                    per = float(lf.get(per_attr))
                    dec = float(lf.get(dec_attr))
                    ratio = (per + dec / 10.0) / 100.0
                    return ratio
                except (ValueError, TypeError):
                    return None

            stutter['minus'] = _stutter_ratio('StutterPer_N_L4',
                                              'DecimalStutterPer_N_L4')
            stutter['plus'] = _stutter_ratio('StutterPer_N_R4',
                                             'DecimalStutterPer_N_R4')
        alleles = []
        for allele_el in locus.findall('Allele'):
            label = allele_el.get('Label', '').strip()
            # GeneMarker stores 'DefSize' (theoretical) and 'Size' (measured
            # from an allelic ladder run). Prefer the measured value.
            raw_size = (allele_el.get('Size') or allele_el.get('DefSize'))
            try:
                size = float(raw_size)
                left_bin = float(allele_el.get('Left_Binning'))
                right_bin = float(allele_el.get('Right_Binning'))
            except ValueError:
                continue
            virtual = allele_el.get('Control') == '1'
            alleles.append({
                'label': label, 'size': size, 'left_bin': left_bin,
                'right_bin': right_bin, 'virtual': virtual,
            })
        markers[marker_name] = {
            'dye': dye, 'min_size': min_size, 'max_size': max_size,
            'stutter': stutter, 'alleles': alleles,
        }
    return {panel_name: markers}


# ---------------------------------------------------------------------------
# OSIRIS LadderInfo XML parser
# ---------------------------------------------------------------------------

def _xml_root_tag(path):
    # Return the root element tag of an XML file, or None on any error.
    try:
        return _xmlparse(path).getroot().tag
    except Exception:
        return None


def parse_osiris(xml_path, default_bin=0.5):
    # Parse an OSIRIS LadderInfo XML file (NIST OSIRIS MarkerSet schema).
    # Both v2.0 (MarkerSet.xsd) and v2.7 (MarkerSetV4.xsd) are handled
    # identically — ILS search-region fields are not used; only marker names,
    # channel mapping, size ranges, and ladder allele BPs are extracted. Allele
    # bin widths default to ±default_bin bp. Stutter thresholds are not stored
    # in this format. Returns the same unified dict as other parsers.
    # Multiple <Set> elements per file each become a separate top-level key.
    tree = _xmlparse(xml_path)
    root = tree.getroot()
    result = {}

    for kit_set in root.findall('.//Kits/Set'):
        panel_name = (kit_set.findtext('Name')).strip()
        # Kit channel number → lowercase colour word (blue/green/yellow/…)
        ch_to_color = {
            int(ch_el.findtext('KitChannelNumber')): ch_el.findtext('Color').strip().lower()
            for ch_el in kit_set.findall('FsaChannelMap/Channel')
        }
        markers = {}
        for locus in kit_set.findall('Locus'):
            marker_name = (locus.findtext('Name')).strip()
            ch_num = int(locus.findtext('Channel'))
            min_bp = float(locus.findtext('MinBP'))
            max_bp = float(locus.findtext('MaxBP'))
            alleles = []
            for a in locus.findall('LadderAlleles/Allele'):
                label = a.findtext('Name').strip()
                size = float(a.findtext('BP'))
                alleles.append({
                    'label': label, 'size': size, 'left_bin': default_bin,
                    'right_bin': default_bin, 'virtual': False,
                })
            markers[marker_name] = {
                'dye': ch_to_color.get(ch_num), 'min_size': min_bp,
                'max_size': max_bp, 'stutter': {'minus': None, 'plus': None},
                'alleles': alleles,
            }
        if markers:
            result[panel_name] = markers
    return result


# ---------------------------------------------------------------------------
# Auto-detect loader
# ---------------------------------------------------------------------------

def has_bin_data(panel_data):
    # Return True if at least one allele in the panel has a precise size.
    # Detects whether a companion Bins file was loaded for a GeneMapper panel.
    return any(a['size'] is not None
               for panel in panel_data.values()
               for marker in panel.values()
               for a in marker['alleles'])


def has_stutter_data(panel_data):
    # Return True if at least one marker in the panel has a stutter ratio.
    return any(marker['stutter']['minus'] is not None
               or marker['stutter']['plus'] is not None
               for panel in panel_data.values() for marker in panel.values())


def load_panel(path):
    # Detect format from file extension and delegate to the right parser.
    # ``.xml`` → GeneMarker (parse_genemarker)
    # ``.txt`` → GeneMapper (parse_genemapper, with auto-detected bins)
    # Returns the unified panel dict.
    if path.lower().endswith('.xml'):
        return parse_genemarker(path)
    return parse_genemapper(path)


# ---------------------------------------------------------------------------
# Allele binning engine
# ---------------------------------------------------------------------------

def assign_alleles(peak_sizes, peak_channel_indices, panel_markers):
    # Assign allele labels to sized peaks using one panel's marker data.

    # Algorithm (per peak):
    # 1. Map the peak's 1-based channel index to a colour word via CHANNEL_COLOR
    #    (1=blue, 2=green, 3=yellow, 4=red, 5=orange, 6=purple, 7=gray). Using
    #    channel index instead of dye name is instrument-agnostic: dye trade names
    #    vary across kits, but channel order is fixed by the CE instrument.
    # 2. Skip markers whose dye colour differs (when both are known).
    # 3. Skip markers whose [min_size, max_size] range does not contain the peak.
    # 4. If allele-level bin data exist, assign the first allele whose acceptance
    #    window [size − left_bin, size + right_bin] covers the peak and stop.
    # 5. If no bin data are present, annotate with the marker name and '?'.

    # Parameters
    # ----------
    # peak_sizes : sequence of float
    # peak_channel_indices : sequence of int  — 1-based channel numbers
    # panel_markers : dict[marker_name -> marker_entry] — one panel's data

    # Returns
    # -------
    # list of str, same length as peak_sizes.
    # Format: "MarkerName:Allele"  e.g. "D3S1358:14"
    #         "MarkerName:?"       for range-only hits (no bin data)
    #         "MarkerName:14*"     for virtual / allelic-ladder alleles
    #         "OL"                 for no match (out of ladder)

    _has_bins = {
        name: any(a['size'] is not None for a in m['alleles'])
        for name, m in panel_markers.items()
    }
    results = [''] * len(peak_sizes)
    for i, (size, ch_idx) in enumerate(zip(peak_sizes, peak_channel_indices)):
        peak_color = CHANNEL_COLOR.get(int(ch_idx))
        for marker_name, marker in panel_markers.items():
            if peak_color and marker['dye'] and peak_color != marker['dye']:
                continue
            mn, mx = marker['min_size'], marker['max_size']
            if mn and mx and (size < mn or size > mx):
                continue
            has_bins = _has_bins[marker_name]
            if has_bins:
                for allele in marker['alleles']:
                    if allele['size'] is None:
                        continue
                    if (allele['size'] - allele['left_bin']
                        <= size <= allele['size'] + allele['right_bin']):
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


# ---------------------------------------------------------------------------
# Panel library (persistent JSON store)
# ---------------------------------------------------------------------------

def load_panel_library(path):
    # Return the full panel dict from the library JSON file.
    # Returns an empty dict if the file does not exist or has an incompatible
    # version number.
    if not isfile(path):
        return {}
    with open(path, encoding='utf-8') as fh:
        data = json_load(fh)
    if data.get('_version') != _LIBRARY_VERSION:
        return {}
    return data.get('panels', {})


def save_panel_library(panels, path):
    # Merge *panels* into the library at *path* and write it back.
    # Creates the directory if necessary.  Existing panels with the same name
    # are overwritten; all others are preserved.
    library = load_panel_library(path)
    library.update(panels)
    makedirs(dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as fh:
        json_dump({'_version': _LIBRARY_VERSION, 'panels': library},
                  fh, indent=2, ensure_ascii=False)
