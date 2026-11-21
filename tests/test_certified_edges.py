"""Analytic identities and failure cases for certified edge integration."""

import pytest
from sage.all import ComplexBallField, QQ, QQbar, RealBallField

from hyperell_regulator.regulator.balls import edge_integrals
from hyperell_regulator.regulator.polynomials import ring


def test_endpoint_log_integrals_at_100_digits():
    # Beta(1/2,1/2) and its logarithmic derivative. These identities
    # exercise root cancellation, both singular tails and complex log seeds.
    x = ring().gen()
    result = edge_integrals(x * (1 - x), 0, 1,
                            centers=[0, 1], powers=[0], digits=100)
    C = ComplexBallField(result['bits'] + 80)
    period = result['periods'][0]
    left, right = (row[0] for row in result['logarithmic'])
    expected = -2 * C.pi() * C(2).log()
    assert period.overlaps(C.pi())
    assert left.overlaps(expected)
    assert right.overlaps(expected + C.gen(0) * C.pi()**2)
    for value in (period, left, right):
        assert RealBallField(result['bits'])(value.rad()) < QQ(1) / 10**100
    assert result['endpoint_bounds']['logarithmic'][0][0] > 0


def test_complex_branch_edge_and_regular_log_endpoint():
    x = ring().gen()
    F = -4*x**5 + 13*x**4 - 18*x**3 + 15*x**2 - 6*x + 1
    roots = F.roots(QQbar, multiplicities=False)
    endpoint = next(r for r in roots if r.imag() > 0 and r.real() < QQ(2)/5)
    lo = edge_integrals(F, 0, endpoint, centers=[0, 1], digits=20)
    hi = edge_integrals(F, 0, endpoint, centers=[0, 1], digits=35)
    for a, b in zip(lo['periods'] + sum(lo['logarithmic'], []),
                    hi['periods'] + sum(hi['logarithmic'], [])):
        assert a.overlaps(b)
        assert a.rad() < QQ(1)/10**20
        assert b.rad() < QQ(1)/10**35


def test_edge_rejects_interior_singularities():
    x = ring().gen()
    with pytest.raises(ValueError, match='another branch point'):
        edge_integrals(x*(x-1)*(x-2), 0, 2)
    with pytest.raises(ValueError, match='interior logarithmic center'):
        edge_integrals(x*(1-x), 0, 1, centers=[QQ(1)/2])
    with pytest.raises(ValueError, match='squarefree'):
        edge_integrals(x**2, 0, 1)
