"""Reduction onto the two branches of an übereven principal cluster."""

import pytest

from sage.all import PolynomialRing, ZZ

from hyperell_regulator import HyperellipticCurveWithReduction


@pytest.mark.parametrize('y_scale', [1, 7])
def test_constant_component_distinguishes_branches(y_scale):
    # Three twins of depth 1 at p=7, with no odd children or twins of depth
    # 1/2.  The component equation is constant, and the dependency returns a
    # scalar rather than a polynomial.  Scaling also exercises nu_s = 2.
    S = PolynomialRing(ZZ, 'x')
    x = S.gen()
    f = (x-1)*(x-64)*(x-2)*(x-128)*(x-3)*(x-192)
    C = HyperellipticCurveWithReduction(y_scale**2 * f)
    P = (0, y_scale * 3072)
    Q = (0, -y_scale * 3072)
    assert C.f(P[0]) == P[1]**2
    R = C.cluster_picture(7)
    assert C.is_semistable_at(7) is True
    assert R.is_ubereven()
    assert R.nu() == (0 if y_scale == 1 else 2)
    assert C._component_polynomial(R).degree() == 0

    result = C.same_component(P, Q, 7)
    assert result['answer'] is False
    rP, rQ = result['P'], result['Q']
    assert rP.is_smooth and rQ.is_smooth
    assert rP.component.sign == -rQ.component.sign
    X, Y = rP.point
    # The square model has y-coordinate 2*y; both reductions must satisfy
    # its equation and keep their opposite signs.
    assert X == 0 and Y == Y.parent()(2 * 3072)
    assert rQ.point == (X, -Y)
    assert Y**2 == C._component_polynomial(R)(X)
    assert C.same_component(P, P, 7)['answer'] is True


def test_constant_component_rejects_point_off_component():
    S = PolynomialRing(ZZ, 'x')
    x = S.gen()
    C = HyperellipticCurveWithReduction(
        (x-1)*(x-64)*(x-2)*(x-128)*(x-3)*(x-192))
    R = C.cluster_picture(7)
    K = R.roots()[0].parent()
    with pytest.raises(ArithmeticError, match='not on the component'):
        C._reduce_affine(R, K(0), K(1))


def test_ubereven_component_with_twins_of_half_depth():
    # Here G(X) is nonconstant and the splitting field has ramification 2.
    S = PolynomialRing(ZZ, 'x')
    x = S.gen()
    C = HyperellipticCurveWithReduction(
        ((x-1)**2-7)*((x-2)**2-7)*((x-3)**2-7))
    R = C.cluster_picture(7)
    assert R.roots()[0].parent().absolute_e() == 2
    assert C._component_polynomial(R).degree() == 6
    result = C.same_component((4, 6), (4, -6), 7)
    assert result['answer'] is False
    rP, rQ = result['P'], result['Q']
    assert rP.is_smooth and rQ.is_smooth
    assert rP.point == (4, 5)
    assert rP._y_datum == 2 and rQ._y_datum == -rP._y_datum
    # The dependency can represent this nonconstant component; its reduction
    # gives an independent check of our coordinates.
    pt = R.red(rP._raw)
    assert rP.point == (pt[0], pt[1])


def test_report_curve_with_constant_component():
    from hyperell_regulator.components import report_curve

    # The cached curve that triggers the Beilinson report's original failure.
    result = report_curve(dict(
        label='3447.a.837621.1', cond=3447,
        eqn='[[40,21,20,5,3],[0,0,1,1]]',
        rat_pts=[[1, -1, 0], [1, 0, 0]]), verbose=False)
    assert len(result['pairs']) == 1
    assert result['pairs'][0]['per_prime'] == {3: False, 383: True}
    assert result['pairs'][0]['verdict'] is False
