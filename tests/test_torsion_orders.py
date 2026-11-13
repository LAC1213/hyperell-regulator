"""Exact divisor orders on models with and without rational Weierstrass points."""

from unittest.mock import patch

import pytest
from sage.all import Infinity, PolynomialRing, QQ

from hyperell_regulator.components import (
    HyperellipticCurveWithReduction, divisor_class_order, odd_degree_model,
    _order_by_principality,
)
from hyperell_regulator.experiments import beilinson as b
from hyperell_regulator.regulator import torsion_order


x = PolynomialRing(QQ, 'x').gen()


def curve_388():
    # LMFDB 388.a.776.1: J(Q) = Z/21, with no rational Weierstrass point.
    return HyperellipticCurveWithReduction(x+2*x*x-x**4, 1+x+x**3)


@pytest.mark.parametrize('P,Q,expected', [
    ((0, 0), ('infinity', 0), 7),
    (('infinity', 0), (0, 0), 7),
    ((0, 0), ('infinity', -1), 21),
    (('infinity', 0), ('infinity', -1), 21),
    ((0, 0), (-1, 1), 7),
])
def test_even_model_without_weierstrass_point(P, Q, expected):
    C = curve_388()
    assert C.F.degree() == 6 and C.F.roots(QQ) == []
    assert odd_degree_model(C) is None
    assert divisor_class_order(C, P, Q) == expected


def test_even_and_odd_representations_agree():
    # 277.a.277.1 admits an independent Mumford arithmetic computation.
    C = HyperellipticCurveWithReduction(-x*x-x, x**3+x*x+x+1)
    for P, Q in [((0, 0), (0, -1)), (('infinity', 0), ('infinity', -1))]:
        assert _order_by_principality(C, P, Q, 30) == divisor_class_order(C, P, Q) == 15


def test_even_finite_points_with_no_rational_infinity():
    C = HyperellipticCurveWithReduction(2*x**6+1)
    assert not QQ(C.F.leading_coefficient()).is_square()
    # (y-1)/x^3 has divisor 3(0,1)-3(0,-1) on the original model.
    assert divisor_class_order(C, (0, 1), (0, -1)) == 3
    assert torsion_order(C.F, (0, 2), (0, -2)) == 3


def test_even_nontorsion_and_identity():
    C = HyperellipticCurveWithReduction(x**6+x+1)
    assert divisor_class_order(C, (0, 1), (0, -1)) is Infinity
    assert divisor_class_order(C, (0, 1), (0, 1), primes=3) == 1


def test_missing_good_primes_does_not_claim_nontorsion():
    assert divisor_class_order(curve_388(), (0, 0), ('infinity', 0), primes=3) is None


def test_already_odd_model_without_finite_weierstrass_point():
    C = HyperellipticCurveWithReduction(x**5+(x+1)**2)
    assert C.F.roots(QQ) == []
    assert odd_degree_model(C) is not None
    assert divisor_class_order(C, (0, -1), ('infinity', 0)) == 5


def test_odd_square_model_with_nonzero_original_infinity_coordinate():
    C = HyperellipticCurveWithReduction(x**5+(x+1)**2-x**6, 2*x**3)
    assert C.F.degree() == 5
    assert divisor_class_order(C, (0, -1), ('infinity', -1)) == 5


def test_explicit_infinities_bounded_search_and_known_multiple():
    F = curve_388().F
    assert torsion_order(F, (0, 1), ('infinity', 1), bound=6) is None
    assert torsion_order(F, (0, 1), ('infinity', 1), bound=1, multiple=21) == 7
    assert torsion_order(F, (0, 1), ('infinity', -1), multiple=21) == 21
    assert torsion_order(F, (0, 1), 'infinity', bound=7) == 7
    assert torsion_order(F, (0, 1), ('infinity', 1), multiple=3) is None
    with pytest.raises(ValueError, match='positive integer'):
        torsion_order(F, (0, 1), ('infinity', 1), multiple=0)


def test_nonmonic_square_model_infinity_coordinates():
    F = 9*curve_388().F
    assert torsion_order(F, (0, 3), ('infinity', 3), multiple=21) == 7
    assert torsion_order(F, (0, 3), ('infinity', -3), multiple=21) == 21


def test_invalid_or_nonrational_infinity_is_rejected():
    with pytest.raises(ValueError, match='point at infinity'):
        torsion_order(curve_388().F, (0, 1), ('infinity', 0))
    with pytest.raises(ValueError, match='rational points'):
        torsion_order(2*x**6+1, (0, 1), 'infinity')


def test_even_orders_recover_old_pair_records_without_rank_zero():
    C = curve_388()
    record = dict(pairs=[dict(P=P, Q=Q, same_component_everywhere=True,
                             order_of_P_minus_Q=None)
                        for P, Q in [((0, 0), ('infinity', 0)),
                                     ((-1, 1), ('infinity', 0))]])
    triples = b.find_triples(record, C, row=dict(mw_invs=[0, 21]))
    assert len(triples) == 1
    assert (triples[0]['order_P'], triples[0]['order_Q']) == (7, 7)
    assert triples[0]['torsion_certified_by'] == 'order in J(Q)'
    assert triples[0]['mw_rank'] == 1


def test_rank_zero_remains_a_fallback_for_unavailable_exact_orders():
    C = curve_388()
    record = dict(pairs=[dict(P=P, Q=Q, same_component_everywhere=True,
                             order_of_P_minus_Q=None)
                        for P, Q in [((0, 0), ('infinity', 0)),
                                     ((-1, 1), ('infinity', 0))]])
    with patch.object(b, 'divisor_class_order', return_value=None):
        triples = b.find_triples(record, C, row=dict(mw_invs=[21]))
        assert len(triples) == 1
        assert triples[0]['torsion_certified_by'] == 'LMFDB rank 0'
        assert triples[0]['torsion_exponent'] == 21
    with patch.object(b, 'divisor_class_order', return_value=Infinity):
        assert b.find_triples(record, C, row=dict(mw_invs=[21])) == []
