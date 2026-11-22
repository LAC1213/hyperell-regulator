"""Exact algebraic centers and logarithmic coincidences at branch points."""

from sage.all import PolynomialRing, QQ, QQbar

from hyperell_regulator.regulator import regulators_stokes
from hyperell_regulator.regulator.stokes import Cut


def test_algebraic_branch_center_is_one_cut_vertex():
    x = PolynomialRing(QQ, 'x').gen()
    F = x*(x*x-2)*(x*x-3)
    a = QQbar(2).sqrt()
    cut = Cut(F, [a], prec=25)
    assert len(cut.nodes) == len(cut.roots) == 5
    assert sum(z == cut.ar.C(a) for z in cut.nodes) == 1


def test_complex_conjugate_algebraic_centers_give_equal_real_moments():
    x = PolynomialRing(QQ, 'x').gen()
    I = QQbar(-1).sqrt()
    a = QQbar(2).sqrt()*I
    periods, cols = regulators_stokes(x**5-x, [a, a.conjugate()], prec=20)
    assert all(abs(u-v) < 1e-17 for u, v in zip(*cols))
    assert all(abs(v) < 1e6 for v in periods+cols[0])


def test_certified_moments_at_algebraic_centres_and_pole():
    from hyperell_regulator.regulator.balls import narrow

    x = PolynomialRing(QQ, 'x').gen()
    F = x*(x*x-2)*(x*x-3)
    pole = QQbar(2).sqrt()
    a = QQbar(-1).sqrt()/2
    per, cols = regulators_stokes(F, [pole, a, a.conjugate()], pole=pole,
                                  prec=15, certified=True)
    assert all(v == 0 for v in cols[0])
    assert any(not v.contains_zero() for v in cols[1])
    assert all(v.overlaps(w) for v, w in zip(cols[1], cols[2]))
    assert all(narrow(v, 15) for v in per+cols[1]+cols[2])
