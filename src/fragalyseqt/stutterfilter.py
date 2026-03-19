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

"""Stutter peak detection and labelling for STR data.

Algorithm (allele-label space)
─────────────────────────────
For each locus, after allele binning:
  • A peak labelled MARKER:A whose height is below the panel n-1 stutter
    ratio × the height of MARKER:A+1 is relabelled "MARKER:A:Stutter"
    (n-1 stutter of A+1).
  • A peak labelled MARKER:A whose height is below the panel n+1 stutter
    ratio × the height of MARKER:A-1 is relabelled "MARKER:A:Stutter"
    (n+1 stutter of A-1).

Thresholds come exclusively from the panel's per-marker stutter data.
Markers without stutter data are skipped (no filtering applied).

Both integer alleles ("14") and microvariants ("9.3") are candidates.
OL, ILS, and unlabelled peaks are never modified.
Virtual-allele peaks (labelled "MARKER:14*") are also eligible candidates.
Neighbour lookup uses allele ± 1, so a microvariant 9.3 checks for parents
at 10.3 (n-1) and 8.3 (n+1).  Allele numbers are rounded to one decimal
place to guarantee consistent float dict keys.

Working in allele-number space avoids the need to know the repeat-unit size
in bp.  True adjacent-allele heterozygotes are not falsely flagged because
typical inter-allele height balance (> 60 %) greatly exceeds normal stutter
thresholds (≤ 15 %).
"""


from .panelparser import _CHANNEL_INDEX_TO_COLOR as _CH_COLOR

# Reverse map: colour word → 1-based channel index
_COLOR_TO_CHANNEL = {v: k for k, v in _CH_COLOR.items()}


def _try_allele(label):
    """Return the allele number (float, rounded to 1 dp) if *label* encodes
    a numeric allele (integer or microvariant), otherwise return None.

    Accepted formats:
      "14"              — plain integer
      "9.3"             — plain microvariant
      "D3S1358:14"      — marker-prefixed integer
      "D3S1358:9.3"     — marker-prefixed microvariant
      "D3S1358:14*"     — virtual allele (strip the asterisk)
    Rejected (return None):
      "D3S1358:OL"      — off-ladder
      "ILS", ""         — special labels
    """
    if not label:
        return None
    # Strip marker prefix if present.
    allele_str = label.split(":", 1)[-1].rstrip("*")
    try:
        return round(float(allele_str), 1)
    except (ValueError, TypeError):
        return None


def apply_stutter_filter(peaksizes, peakheights, peakchannels, peakalleles,
                         panel, dye_names,):
    """Return a copy of *peakalleles* with stutter peaks relabelled.

    Parameters
    ----------
    peaksizes, peakheights, peakchannels : array-like
        Parallel peak data arrays from FileState.
    peakalleles : list[str]
        Allele labels as assigned by assign_alleles() and ILS stamping.
        Expected format: "MARKER_NAME:allele" or "MARKER_NAME:allele*".
    panel : dict
        The currently selected panel dict (marker_name → marker_info).
        marker_info must contain a "stutter" key with "minus"/"plus" floats;
        markers without stutter data are skipped.
    dye_names : list[str]
        state.Dye — dye names in channel order.

    Returns
    -------
    list[str]
        New allele label list; input is not mutated.
    """
    result = list(peakalleles)

    for marker, info in panel.items():
        stutter = info.get("stutter")
        m_thr = stutter.get("minus")
        p_thr = stutter.get("plus")
        # Skip markers with no stutter data at all.
        if m_thr is None and p_thr is None:
            continue
        ch_idx = _COLOR_TO_CHANNEL.get(info["dye"].lower())
        if ch_idx is None:
            continue
        # Indices of all peaks in this locus.
        locus_idx = [i for i, (ch, sz) in enumerate(zip(peakchannels,
                                                        peaksizes))
                     if int(ch) == ch_idx and
                     (info["min_size"] <= float(sz) <= info["max_size"])]
        if not locus_idx:
            continue

        # Build {allele_int: (peak_index, height)}.
        # When two peaks share the same integer allele number (rare but
        # possible with OL peaks nearby), keep the taller one as the parent
        # reference so the shorter one gets tested against it.
        allele_map = {}
        for i in locus_idx:
            a = _try_allele(result[i])
            if a is None:
                continue
            h = float(peakheights[i])
            if a not in allele_map or h > allele_map[a][1]:
                allele_map[a] = (i, h)

        # Apply stutter thresholds.
        # checks: (allele offset to parent, threshold)
        checks = []
        if m_thr is not None:
            checks.append((+1, m_thr))   # n-1 stutter: parent is one repeat up
        if p_thr is not None:
            checks.append((-1, p_thr))   # n+1 stutter: parent is one repeat down
        for i in locus_idx:
            a = _try_allele(result[i])
            if a is None:
                continue
            h = float(peakheights[i])
            for offset, thr in checks:
                parent = allele_map.get(round(a + offset, 1))
                if parent is not None and parent[0] != i and h / parent[1] < thr:
                    result[i] = ""
                    break

    return result
