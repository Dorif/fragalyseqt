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

# Forensic STR statistics: RMP, identity LR, and kinship index (KI).
#
# Formulas follow Fung & Hu (2008) "Statistical DNA Forensics", Ch. 4–5.
# HWE case: Table 5.2. Subdivided-population (θ-corrected) case: Table 5.5.
#
# RELATIONSHIPS uses standard IBD notation (κ0, κ1, κ2).
# Fung's tables use (k0, 2k1, k2), so κ1 is halved when entering formulas.
# The helper _f(θ, p) = θ + (1−θ)p appears throughout both tables.

from __future__ import annotations
from collections import Counter
from dataclasses import dataclass, field

from .freqdb import FrequencyTable, get_allele_freq, _normalize_allele


@dataclass
class AlleleCall:
    marker: str
    allele1: str
    allele2: str | None


@dataclass
class Relationship:
    name: str
    k0: float
    k1: float
    k2: float


RELATIONSHIPS: dict[str, Relationship] = {
    'Identical twins': Relationship('Identical twins', 0, 0, 1),
    'Parent-Child': Relationship('Parent-Child', 0, 1, 0),
    'Full Siblings': Relationship('Full Siblings', 0.25, 0.5, 0.25),
    'Half Siblings': Relationship('Half Siblings', 0.5, 0.5, 0),
    'Grandparent-Grandchild': Relationship('Grandparent-Grandchild', 0.5, 0.5,
                                           0),
    'Uncle-Niece': Relationship('Uncle-Niece', 0.5, 0.5, 0),
    'First Cousins': Relationship('First Cousins', 0.75, 0.25, 0),
    'Unrelated': Relationship('Unrelated', 1, 0, 0),
}


@dataclass
class LocusResult:
    marker: str
    alleles_q: tuple
    alleles_r: tuple
    freq_a1: float
    freq_a2: float | None
    locus_stat: float
    included: bool
    note: str


@dataclass
class ComparisonResult:
    mode: str
    relationship: str | None
    freq_table: str
    combined_stat: float
    log10_stat: float
    n_loci: int
    n_excluded: int
    verbal_scale: str
    loci: list[LocusResult] = field(default_factory=list)
    n_only_q: int = 0
    n_only_r: int = 0


def verbal_conclusion(log10_stat: float) -> str:
    if log10_stat > 6:
        return 'Extremely strong support'
    if log10_stat > 4:
        return 'Very strong support'
    if log10_stat > 2:
        return 'Strong support'
    if log10_stat > 1:
        return 'Moderate support'
    if log10_stat > 0:
        return 'Limited support'
    if log10_stat >= 0:
        return 'Inconclusive'
    return 'Supports exclusion'


def _f(theta: float, p: float) -> float:
    return theta + (1.0 - theta) * p


def locus_rmp(table: FrequencyTable, marker: str,
              allele1: str, allele2: str,
              theta: float = 0.0) -> float:
    a1 = _normalize_allele(allele1)
    a2 = _normalize_allele(allele2)
    p1 = get_allele_freq(table, marker, a1)
    p2 = get_allele_freq(table, marker, a2)
    denom = (1.0 + theta) * (1.0 + 2.0 * theta)
    if a1 == a2:
        fi = _f(theta, p1)
        return (theta + fi) * (2.0 * theta + fi) / denom
    return 2.0 * _f(theta, p1) * _f(theta, p2) / denom


def combined_rmp(table: FrequencyTable,
                 calls: list[AlleleCall],
                 theta: float = 0.0) -> float:
    result = 1.0
    for c in calls:
        a2 = c.allele2 if c.allele2 is not None else c.allele1
        result *= locus_rmp(table, c.marker, c.allele1, a2, theta)
    return result


def identity_lr(table: FrequencyTable,
                calls: list[AlleleCall],
                theta: float = 0.0) -> float:
    return 1.0 / combined_rmp(table, calls, theta)


def locus_ki(table: FrequencyTable, marker: str, profile1: tuple[str, str],
             profile2: tuple[str, str], rel: Relationship,
             theta: float = 0.0) -> float:
    a1 = _normalize_allele(profile1[0])
    a2 = _normalize_allele(profile1[1])
    b1 = _normalize_allele(profile2[0])
    b2 = _normalize_allele(profile2[1])

    get = lambda a: get_allele_freq(table, marker, a)
    k0, k1, k2 = rel.k0, rel.k1, rel.k2
    t = theta

    shared = list((Counter([a1, a2]) & Counter([b1, b2])).elements())
    n = len(shared)

    if n == 0:
        return k0

    if n == 2:
        if a1 == a2:
            fi = _f(t, get(a1))
            return (k0 + k1 * (1 + 2*t) / (2*t + fi)
                    + k2 * (1 + t) * (1 + 2*t) / ((t + fi) * (2*t + fi)))
        fi = _f(t, get(a1))
        fj = _f(t, get(a2))
        return (k0 + k1 * (1 + 2*t) * (fi + fj) / (4 * fi * fj)
                + k2 * (1 + t) * (1 + 2*t) / (2 * fi * fj))

    fi = _f(t, get(shared[0]))
    if a1 == a2 or b1 == b2:
        return k0 + k1 * (1 + 2*t) / (2 * (t + fi))
    return k0 + k1 * (1 + 2*t) / (4 * fi)


def combined_ki(table: FrequencyTable,
                calls1: list[AlleleCall],
                calls2: list[AlleleCall],
                rel: Relationship,
                theta: float = 0.0) -> float:
    map1 = {c.marker: (c.allele1, c.allele2 or c.allele1) for c in calls1}
    map2 = {c.marker: (c.allele1, c.allele2 or c.allele1) for c in calls2}
    result = 1.0
    for marker in map1:
        if marker in map2:
            result *= locus_ki(table, marker, map1[marker], map2[marker], rel,
                               theta)
    return result
