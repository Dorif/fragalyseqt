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


# ---------------------------------------------------------------------------
# locus_ki — the k2 term of Fung & Hu Table 5.5 row 5 (both heterozygous,
# two shared alleles).  Every book example above uses a relationship with
# k2 = 0 (half sibs, first cousins, unrelated), so none of them exercises
# the k2 term; these tests cover that gap for Identical twins (k2 = 1) and
# Full Siblings (k2 = 0.25).
#
# Table 5.5 row 5 (p. 87):
#   LR = k0 + { k1(1+2θ)[2θ + (1−θ)(pi+pj)] + k2(1+θ)(1+2θ) }
#            / { 2[θ + (1−θ)pi][θ + (1−θ)pj] }
# ---------------------------------------------------------------------------

def test_ki_twins_het_equals_inverse_rmp(book_table):
    # For identical twins (k0=k1=0, k2=1) the kinship index must equal 1/RMP:
    # the two profiles are the same person's genotype.  Independent of how
    # Table 5.5 is transcribed.
    for marker, a1, a2 in (('D3S1358', '15', '17'), ('vWA', '14', '15')):
        ki = locus_ki(book_table, marker, (a1, a2), (a1, a2), IDENTICAL)
        rmp = locus_rmp(book_table, marker, a1, a2)
        assert abs(ki - 1.0 / rmp) < 1e-9, marker

def test_ki_twins_hom_equals_inverse_rmp(book_table):
    ki = locus_ki(book_table, 'FGA', ('22', '22'), ('22', '22'), IDENTICAL)
    rmp = locus_rmp(book_table, 'FGA', '22', '22')
    assert abs(ki - 1.0 / rmp) < 1e-9

def test_ki_twins_het_equals_inverse_rmp_with_theta(book_table):
    for theta in (0.01, 0.03):
        ki = locus_ki(book_table, 'D3S1358', ('15', '17'), ('15', '17'),
                      IDENTICAL, theta)
        rmp = locus_rmp(book_table, 'D3S1358', '15', '17', theta)
        assert abs(ki - 1.0 / rmp) < 1e-9, theta

def test_ki_twins_het_matches_table_5_5_row_5(book_table):
    # k0=0, k1=0, k2=1, θ=0  =>  LR = 1 / (2 p15 p17)
    p15, p17 = 0.331, 0.239
    ki = locus_ki(book_table, 'D3S1358', ('15', '17'), ('15', '17'), IDENTICAL)
    assert abs(ki - 1.0 / (2 * p15 * p17)) < 1e-9

def test_ki_full_sibs_het_matches_table_5_5_row_5(book_table):
    # Full sibs: kappa = (0.25, 0.5, 0.25); Fung's k1 = kappa1/2 = 0.25.
    # LR = k0 + k1(pi+pj)/(2 pi pj) + k2/(2 pi pj)
    full_sib = RELATIONSHIPS['Full Siblings']
    p15, p17 = 0.331, 0.239
    expect = (full_sib.k0
              + (full_sib.k1 / 2) * (p15 + p17) / (2 * p15 * p17)
              + full_sib.k2 / (2 * p15 * p17))
    ki = locus_ki(book_table, 'D3S1358', ('15', '17'), ('15', '17'), full_sib)
    assert abs(ki - expect) < 1e-9

def test_ki_full_sibs_het_matches_table_5_5_row_5_with_theta(book_table):
    full_sib = RELATIONSHIPS['Full Siblings']
    p15, p17 = 0.331, 0.239
    for t in (0.01, 0.03):
        fi = t + (1 - t) * p15
        fj = t + (1 - t) * p17
        expect = (full_sib.k0
                  + (full_sib.k1 / 2) * (1 + 2*t) * (2*t + (1-t)*(p15 + p17))
                  / (2 * fi * fj)
                  + full_sib.k2 * (1 + t) * (1 + 2*t) / (2 * fi * fj))
        ki = locus_ki(book_table, 'D3S1358', ('15', '17'), ('15', '17'),
                      full_sib, t)
        assert abs(ki - expect) < 1e-9, t

def test_ki_k2_term_dominates_for_twins_over_full_sibs(book_table):
    # Sanity ordering: identical twins > full sibs > half sibs > unrelated
    # at a shared heterozygous locus.
    args = ('D3S1358', ('15', '17'), ('15', '17'))
    twins = locus_ki(book_table, *args, IDENTICAL)
    sibs = locus_ki(book_table, *args, RELATIONSHIPS['Full Siblings'])
    half = locus_ki(book_table, *args, HALF_SIB)
    unrel = locus_ki(book_table, *args, UNRELATED)
    assert twins > sibs > half > unrel
