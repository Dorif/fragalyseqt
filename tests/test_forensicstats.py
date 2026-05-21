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

import pytest
from fragalyseqt.freqdb import FrequencyTable
from fragalyseqt.forensicstats import (
    AlleleCall, RELATIONSHIPS,
    verbal_conclusion,
    locus_rmp, combined_rmp, identity_lr,
    locus_ki, combined_ki,
)

HALF_SIB = RELATIONSHIPS['Half Siblings']
FIRST_COS = RELATIONSHIPS['First Cousins']
UNRELATED = RELATIONSHIPS['Unrelated']
IDENTICAL = RELATIONSHIPS['Identical twins']
PARENT_CHILD = RELATIONSHIPS['Parent-Child']


@pytest.fixture
def book_table():
    return FrequencyTable(
        name='book', panel='', population='book',
        markers={
            'D3S1358': {'15': 0.331, '17': 0.239},
            'vWA':     {'14': 0.200, '15': 0.035, '19': 0.150},
            'FGA':     {'22': 0.178, '23': 0.150},
        }
    )


# ---------------------------------------------------------------------------
# RELATIONSHIPS catalogue
# ---------------------------------------------------------------------------

def test_relationships_k_sum_to_one():
    for rel in RELATIONSHIPS.values():
        assert abs(rel.k0 + rel.k1 + rel.k2 - 1.0) < 1e-9, rel.name


# ---------------------------------------------------------------------------
# verbal_conclusion (SWGDAM 2016 scale)
# ---------------------------------------------------------------------------

def test_verbal_extremely_strong():
    assert verbal_conclusion(6.1) == 'Extremely strong support'

def test_verbal_very_strong():
    assert verbal_conclusion(5.0) == 'Very strong support'

def test_verbal_strong():
    assert verbal_conclusion(3.0) == 'Strong support'

def test_verbal_moderate():
    assert verbal_conclusion(1.5) == 'Moderate support'

def test_verbal_limited():
    assert verbal_conclusion(0.1) == 'Limited support'

def test_verbal_inconclusive():
    assert verbal_conclusion(0.0) == 'Inconclusive'

def test_verbal_exclusion():
    assert verbal_conclusion(-1.0) == 'Supports exclusion'


# ---------------------------------------------------------------------------
# locus_rmp
# ---------------------------------------------------------------------------

def test_locus_rmp_het_theta0(book_table):
    assert abs(locus_rmp(book_table, 'D3S1358', '15', '17') - 2*0.331*0.239) < 1e-9

def test_locus_rmp_homo_theta0(book_table):
    assert abs(locus_rmp(book_table, 'FGA', '22', '22') - 0.178**2) < 1e-9

def test_locus_rmp_theta0_equals_no_theta(book_table):
    a = locus_rmp(book_table, 'D3S1358', '15', '17', theta=0.0)
    b = locus_rmp(book_table, 'D3S1358', '15', '17')
    assert a == b

def test_locus_rmp_homo_theta_increases_probability(book_table):
    p0 = locus_rmp(book_table, 'FGA', '22', '22', theta=0.0)
    p1 = locus_rmp(book_table, 'FGA', '22', '22', theta=0.03)
    assert p1 > p0

def test_combined_rmp_is_product(book_table):
    calls = [AlleleCall('D3S1358', '15', '17'), AlleleCall('FGA', '22', '22')]
    expected = (locus_rmp(book_table, 'D3S1358', '15', '17') *
                locus_rmp(book_table, 'FGA', '22', '22'))
    assert abs(combined_rmp(book_table, calls) - expected) < 1e-12

def test_identity_lr_inverse_of_rmp(book_table):
    calls = [AlleleCall('D3S1358', '15', '17')]
    assert abs(identity_lr(book_table, calls) * combined_rmp(book_table, calls) - 1.0) < 1e-9


# ---------------------------------------------------------------------------
# locus_ki — book example (Fung & Hu 2008, Table 5.3, p.82), θ=0
# ---------------------------------------------------------------------------

def test_ki_d3s1358_half_sib(book_table):
    lr = locus_ki(book_table, 'D3S1358', ('15', '17'), ('15', '17'), HALF_SIB)
    assert abs(lr - 1.401) < 0.001

def test_ki_vwa_half_sib(book_table):
    lr = locus_ki(book_table, 'vWA', ('14', '15'), ('15', '19'), HALF_SIB)
    assert abs(lr - 4.071) < 0.001

def test_ki_fga_half_sib(book_table):
    lr = locus_ki(book_table, 'FGA', ('22', '22'), ('22', '23'), HALF_SIB)
    assert abs(lr - 1.904) < 0.001

def test_ki_d3s1358_first_cousin(book_table):
    lr = locus_ki(book_table, 'D3S1358', ('15', '17'), ('15', '17'), FIRST_COS)
    assert abs(lr - 1.2) < 0.001

def test_ki_vwa_first_cousin(book_table):
    lr = locus_ki(book_table, 'vWA', ('14', '15'), ('15', '19'), FIRST_COS)
    assert abs(lr - 2.536) < 0.001

def test_ki_fga_first_cousin(book_table):
    lr = locus_ki(book_table, 'FGA', ('22', '22'), ('22', '23'), FIRST_COS)
    assert abs(lr - 1.452) < 0.001

def test_ki_unrelated_no_shared_is_k0(book_table):
    lr = locus_ki(book_table, 'D3S1358', ('15', '17'), ('14', '13'), UNRELATED)
    assert abs(lr - 1.0) < 1e-9

def test_ki_identical_twins_no_shared_is_zero(book_table):
    assert locus_ki(book_table, 'D3S1358', ('15', '17'), ('14', '13'), IDENTICAL) == 0.0

def test_ki_unrelated_any_genotype_is_one(book_table):
    assert abs(locus_ki(book_table, 'D3S1358', ('15', '17'), ('15', '17'), UNRELATED) - 1.0) < 1e-9
    assert abs(locus_ki(book_table, 'FGA',     ('22', '22'), ('22', '22'), UNRELATED) - 1.0) < 1e-9

def test_ki_profile_order_symmetric(book_table):
    lr1 = locus_ki(book_table, 'vWA', ('14', '15'), ('15', '19'), HALF_SIB)
    lr2 = locus_ki(book_table, 'vWA', ('15', '19'), ('14', '15'), HALF_SIB)
    assert abs(lr1 - lr2) < 1e-9

def test_combined_ki_is_product(book_table):
    calls1 = [AlleleCall('D3S1358', '15', '17'), AlleleCall('FGA', '22', '22')]
    calls2 = [AlleleCall('D3S1358', '15', '17'), AlleleCall('FGA', '22', '23')]
    ki = combined_ki(book_table, calls1, calls2, HALF_SIB)
    assert abs(ki - 1.401 * 1.904) < 0.01
