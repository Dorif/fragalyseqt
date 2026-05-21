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

# Profile comparison orchestration for identity LR and kinship KI.
#
# allele_calls_from_state() extracts allele calls from a FileState object
# (duck-typed: needs .peakalleles and .peakheights parallel lists).
# Allele labels in peakalleles use the format "MARKER:ALLELE"; ILS, OL and
# empty entries are skipped. The two tallest distinct alleles per marker are
# taken as allele1/allele2.
#
# compare_identity(): per-locus LR = 1/RMP for matching loci, excludes
# mismatched loci with note "profile mismatch".
# compare_kinship(): per-locus KI from Fung & Hu Table 5.5.
# Both functions exclude loci absent from the frequency table.

from __future__ import annotations
from csv import writer
from io import StringIO
from math import log10
from .forensicstats import (AlleleCall, Relationship, RELATIONSHIPS, locus_ki,
                            locus_rmp, verbal_conclusion, ComparisonResult,
                            LocusResult,)
from .freqdb import FrequencyTable, get_allele_freq, _normalize_allele


_SEX_MARKERS = frozenset(['amel', 'amelogenin'])
_Y_MARKERS = frozenset(['yindel', 'y_indel'])


def _norm_marker(name: str) -> str:
    return name.strip().lower().replace(' ', '_').replace('-', '_')


def _absent_note(marker: str) -> str:
    n = _norm_marker(marker)
    if n in _SEX_MARKERS:
        return 'sex marker'
    if n in _Y_MARKERS:
        return 'Y-chromosome marker'
    return 'not in freq table'


def _marker_index(table: FrequencyTable) -> dict[str, str]:
    return {_norm_marker(k): k for k in table.markers}


def _resolve(index: dict[str, str], marker: str) -> str | None:
    if marker in index.values():
        return marker
    return index.get(_norm_marker(marker))


def allele_calls_from_state(state) -> list[AlleleCall]:
    by_marker: dict[str, list[tuple[str, float]]] = {}
    for allele_str, height in zip(state.peakalleles, state.peakheights):
        if not allele_str or allele_str in ('OL', 'ILS'):
            continue
        if ':' not in allele_str:
            continue
        marker, allele = allele_str.split(':', 1)
        by_marker.setdefault(marker, []).append((allele, float(height)))

    calls = []
    for marker, pairs in by_marker.items():
        pairs.sort(key=lambda x: -x[1])
        seen: list[str] = []
        for allele, _ in pairs:
            if allele not in seen:
                seen.append(allele)
            if len(seen) == 2:
                break
        if not seen:
            continue
        calls.append(AlleleCall(marker=marker, allele1=seen[0],
                     allele2=seen[1] if len(seen) > 1 else None,))
    return calls


def compare_identity(calls_q: list[AlleleCall], calls_r: list[AlleleCall],
                     table: FrequencyTable, theta: float = 0.0,
                     ) -> ComparisonResult:
    map_q = {c.marker: c for c in calls_q}
    map_r = {c.marker: c for c in calls_r}
    n_only_q = len(set(map_q) - set(map_r))
    n_only_r = len(set(map_r) - set(map_q))
    idx = _marker_index(table)

    loci: list[LocusResult] = []
    for marker in sorted(set(map_q) & set(map_r)):
        q = map_q[marker]
        r = map_r[marker]
        a1q = q.allele1
        a2q = q.allele2 or q.allele1
        a1r = r.allele1
        a2r = r.allele2 or r.allele1

        tkey = _resolve(idx, marker)
        if tkey is None:
            note = _absent_note(marker)
            alleles_q = frozenset([a1q, a2q])
            alleles_r = frozenset([a1r, a2r])
            if note == 'sex marker':
                if alleles_q != alleles_r:
                    note = 'sex marker — sex discordance'
            elif alleles_q != alleles_r:
                note += ' — profile mismatch'
            loci.append(LocusResult(marker, (a1q, q.allele2), (a1r, r.allele2),
                                    table.min_freq, None, 0.0, False, note))
            continue

        p1 = get_allele_freq(table, tkey, a1q)
        p2 = get_allele_freq(table, tkey, a2q) if a2q != a1q else None

        if {a1q, a2q} != {a1r, a2r}:
            rmp = locus_rmp(table, tkey, a1q, a2q, theta)
            loci.append(LocusResult(marker, (a1q, q.allele2), (a1r, r.allele2),
                                    p1, p2, 1.0 / rmp if rmp > 0 else 0.0,
                                    False, 'profile mismatch'))
            continue

        rmp = locus_rmp(table, tkey, a1q, a2q, theta)
        lr = 1.0 / rmp if rmp > 0 else float('inf')
        loci.append(LocusResult(marker, (a1q, q.allele2), (a1r, r.allele2),
                                p1, p2, lr, True, ''))

    return _build_result('identity', None, table, loci, n_only_q, n_only_r)


def compare_kinship(
        calls1: list[AlleleCall],
        calls2: list[AlleleCall],
        table: FrequencyTable,
        rel: Relationship,
        theta: float = 0.0,
) -> ComparisonResult:
    map1 = {c.marker: c for c in calls1}
    map2 = {c.marker: c for c in calls2}
    n_only_q = len(set(map1) - set(map2))
    n_only_r = len(set(map2) - set(map1))
    idx = _marker_index(table)

    loci: list[LocusResult] = []
    for marker in sorted(set(map1) & set(map2)):
        c1 = map1[marker]
        c2 = map2[marker]
        a1 = (c1.allele1, c1.allele2 or c1.allele1)
        a2 = (c2.allele1, c2.allele2 or c2.allele1)

        tkey = _resolve(idx, marker)
        if tkey is None:
            n_shared = len(frozenset(a1) & frozenset(a2))
            base_note = _absent_note(marker)
            if n_shared == 0:
                note = base_note + ' — no shared alleles'
            elif n_shared == 1:
                note = base_note + ' — 1 shared allele'
            else:
                note = base_note
            loci.append(LocusResult(marker, a1, a2, table.min_freq, None,
                                    0.0, False, note))
            continue

        p1 = get_allele_freq(table, tkey, c1.allele1)
        p2 = get_allele_freq(table, tkey, c1.allele2) if c1.allele2 else None

        ki = locus_ki(table, tkey, a1, a2, rel, theta)
        loci.append(LocusResult(marker, a1, a2, p1, p2, ki, True, ''))

    return _build_result('kinship', rel.name, table, loci, n_only_q, n_only_r)


def _build_result(mode: str, relationship: str | None, table: FrequencyTable,
                  loci: list[LocusResult],
                  n_only_q: int = 0, n_only_r: int = 0) -> ComparisonResult:
    included = [locus for locus in loci if locus.included]
    n_excl = len(loci) - len(included)

    if not included:
        combined = 0.0
        lg = float('-inf')
    else:
        combined = 1.0
        for locus in included:
            combined *= locus.locus_stat
        lg = log10(combined) if combined > 0 else float('-inf')

    return ComparisonResult(mode=mode, relationship=relationship,
                            freq_table=table.name, combined_stat=combined,
                            log10_stat=lg, n_loci=len(included),
                            n_excluded=n_excl,
                            verbal_scale=verbal_conclusion(lg), loci=loci,
                            n_only_q=n_only_q, n_only_r=n_only_r)


def search_profiles(db, calls: list[AlleleCall]) -> list[dict]:
    from .refprofile import list_profiles
    map_q = {
        _norm_marker(c.marker): frozenset([
            _normalize_allele(c.allele1),
            _normalize_allele(c.allele2 or c.allele1),
        ])
        for c in calls
    }
    results = []
    for p_info in list_profiles(db):
        alleles = db.get_reference_alleles(p_info['id'])
        map_p = {
            _norm_marker(a['marker']): frozenset([
                _normalize_allele(a['allele1']),
                _normalize_allele(a['allele2'] or a['allele1']),
            ])
            for a in alleles
        }
        common = set(map_q) & set(map_p)
        if not common:
            continue
        matched = sum(1 for m in common if map_q[m] == map_p[m])
        results.append({
            'id': p_info['id'],
            'name': p_info['name'],
            'role': p_info.get('role') or '',
            'matched': matched,
            'common': len(common),
        })
    results.sort(key=lambda x: (-x['matched'], -x['common']))
    return results


def export_comparison_csv(result: ComparisonResult) -> str:
    buf = StringIO()
    w = writer(buf)
    w.writerow(['mode', result.mode])
    if result.relationship:
        w.writerow(['relationship', result.relationship])
    w.writerow(['freq_table', result.freq_table])
    w.writerow(['combined_stat', result.combined_stat])
    w.writerow(['log10_stat', result.log10_stat])
    w.writerow(['verbal_scale', result.verbal_scale])
    w.writerow(['n_loci', result.n_loci])
    w.writerow(['n_excluded', result.n_excluded])
    w.writerow(['n_only_profile1', result.n_only_q])
    w.writerow(['n_only_profile2', result.n_only_r])
    w.writerow([])
    w.writerow(['marker', 'profile1', 'profile2', 'p_a1', 'p_a2', 'locus_stat',
                'included', 'note'])
    for locus in result.loci:
        q = '/'.join(str(a) for a in locus.alleles_q if a is not None)
        r = '/'.join(str(a) for a in locus.alleles_r if a is not None)
        w.writerow([locus.marker, q, r, f'{locus.freq_a1:.6f}',
                    f'{locus.freq_a2:.6f}' if locus.freq_a2 is not None else '',
                    f'{locus.locus_stat:.6f}' if locus.included else '',
                    'yes' if locus.included else 'no', locus.note,])
    return buf.getvalue()
