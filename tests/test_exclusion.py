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

"""Identity comparison must treat a genetic contradiction as an exclusion.

One person cannot carry two different genotypes at the same marker. If loci
that contradict each other were simply dropped from the likelihood-ratio
product, a handful of matching loci would outvote any number of mismatches
and the dialog would report strong support for two obviously different
people.
"""

import pytest

from fragalyseqt.comparison import compare_identity, compare_kinship
from fragalyseqt.forensicstats import RELATIONSHIPS, AlleleCall
from fragalyseqt.freqdb import FrequencyTable

MARKERS = ['D3S1358', 'vWA', 'FGA', 'TH01', 'TPOX']


@pytest.fixture
def table():
    return FrequencyTable(
        name='test', panel='', population='test',
        markers={
            'D3S1358': {'12': 0.20, '14': 0.15, '15': 0.33, '17': 0.24},
            'vWA':     {'12': 0.20, '14': 0.20, '15': 0.04, '19': 0.15},
            'FGA':     {'12': 0.18, '14': 0.15, '15': 0.10, '17': 0.12},
            'TH01':    {'12': 0.22, '14': 0.18, '15': 0.11, '17': 0.14},
            'TPOX':    {'12': 0.25, '14': 0.19, '15': 0.13, '17': 0.16},
        }
    )


def _matching():
    return [AlleleCall(m, '12', '14') for m in MARKERS]


def _with_mismatches(n_match):
    return [AlleleCall(m, '12', '14') if i < n_match
            else AlleleCall(m, '15', '17')
            for i, m in enumerate(MARKERS)]


def test_full_match_supports_identity(table):
    result = compare_identity(_matching(), _matching(), table)
    assert result.n_excluded == 0
    assert result.log10_stat > 0
    assert result.verbal_scale != 'Supports exclusion'


@pytest.mark.parametrize('n_match', [4, 3, 2, 1, 0])
def test_any_mismatch_is_an_exclusion(n_match, table):
    # Regression: mismatching loci used to be dropped from the product, so
    # a single matching locus against nineteen mismatches still reported
    # 'Strong support'.
    result = compare_identity(_matching(), _with_mismatches(n_match), table)
    assert result.log10_stat == float('-inf')
    assert result.verbal_scale == 'Supports exclusion'


def test_one_mismatch_outweighs_many_matches(table):
    # The decisive case: four loci agree, one contradicts.
    result = compare_identity(_matching(), _with_mismatches(4), table)
    assert result.verbal_scale == 'Supports exclusion'
    assert result.n_loci == 4
    assert result.n_excluded == 1


def test_sex_discordance_is_an_exclusion(table):
    q = _matching() + [AlleleCall('AMEL', 'X', 'Y')]
    r = _matching() + [AlleleCall('AMEL', 'X', None)]
    result = compare_identity(q, r, table)
    assert result.verbal_scale == 'Supports exclusion'


def test_matching_sex_marker_is_not_an_exclusion(table):
    q = _matching() + [AlleleCall('AMEL', 'X', 'Y')]
    r = _matching() + [AlleleCall('AMEL', 'X', 'Y')]
    result = compare_identity(q, r, table)
    assert result.log10_stat > 0
    assert result.verbal_scale != 'Supports exclusion'


def test_marker_absent_from_table_is_not_an_exclusion(table):
    # A locus the frequency table cannot evaluate carries no contradiction
    # when the alleles agree, so it must not collapse the conclusion.
    q = _matching() + [AlleleCall('SE33', '19', '25.2')]
    r = _matching() + [AlleleCall('SE33', '19', '25.2')]
    result = compare_identity(q, r, table)
    assert result.log10_stat > 0
    assert result.verbal_scale != 'Supports exclusion'


def test_marker_absent_from_table_still_excludes_on_mismatch(table):
    # ...but differing alleles are a contradiction whether or not the
    # frequency table knows the marker.
    q = _matching() + [AlleleCall('SE33', '19', '25.2')]
    r = _matching() + [AlleleCall('SE33', '15', '16')]
    result = compare_identity(q, r, table)
    assert result.verbal_scale == 'Supports exclusion'


def test_kinship_is_unaffected_by_allele_differences(table):
    # Relatives legitimately differ at a locus: kinship must keep using the
    # full product and must never be collapsed by this rule.
    parent = [AlleleCall(m, '12', '14') for m in MARKERS]
    child = [AlleleCall(m, '12', '17') for m in MARKERS]
    result = compare_kinship(parent, child, table,
                             RELATIONSHIPS['Parent-Child'])
    assert result.log10_stat > 0
    assert result.verbal_scale != 'Supports exclusion'


def test_kinship_not_collapsed_by_unevaluable_locus(table):
    # A marker the frequency table does not know is annotated 'no shared
    # alleles' in kinship mode. That is a normal, informative outcome for
    # relatives and must never collapse the kinship statistic.
    parent = [AlleleCall(m, '12', '14') for m in MARKERS]
    child = [AlleleCall(m, '12', '17') for m in MARKERS]
    parent.append(AlleleCall('SE33', '19', '25.2'))
    child.append(AlleleCall('SE33', '15', '16'))
    result = compare_kinship(parent, child, table,
                             RELATIONSHIPS['Parent-Child'])
    assert result.log10_stat > 0
    assert result.verbal_scale != 'Supports exclusion'
