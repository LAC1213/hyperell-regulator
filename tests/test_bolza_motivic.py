"""Graph traces agree with independently integrated logarithmic columns."""
from sage.all import PolynomialRing, QQ

from hyperell_regulator.experiments.bolza_motivic import regulator_matrix
from hyperell_regulator.regulator import regulators_stokes
from hyperell_regulator.regulator.balls import narrow


def test_graph_unit_regulators_match_logarithmic_representatives():
    M = regulator_matrix(20)
    x = PolynomialRing(QQ, 'x').gen()
    p, (r0, r1) = regulators_stokes(4*(x**5-x), [0, 1], prec=32, certified=True)
    B = M.base_ring()
    log2 = B(2).log()
    for i, w in enumerate((1, 2, 1)):
        expected = (w*(4*(r1[i]-r0[i])-2*log2*p[i]),
                    w*(4*r1[i]-2*log2*p[i]), w*p[i])
        for j, value in enumerate(expected):
            assert (M[i,j]-value).contains_zero()
            assert narrow(M[i,j], 20)
    assert narrow(M.det(), 20) and not M.det().contains_zero()


def test_graph_unit_determinant_matches_standard_lvalue():
    import pytest
    from hyperell_regulator.ball_lfunction import executable
    from hyperell_regulator.regulator.balls import CertificationError
    from hyperell_regulator.experiments.bolza_motivic import lvalue_comparison
    try:
        executable()
    except CertificationError:
        pytest.skip('optional native certified backend not built')
    data = lvalue_comparison(regulator_matrix(20), digits=20)
    out = data['comparison']
    assert out['rational'] == QQ(1)/12
    assert out['relative_agreement'].upper() < QQ(1)/10**20
    assert out['exact_relation_proved'] is False
