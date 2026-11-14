"""MST cycle bases, precision growth and integral changes of basis."""
import json
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pytest
from sage.all import ComplexField,RealField,ZZ,matrix
from hyperell_regulator.regulator import c_invariants
from hyperell_regulator.regulator._periods import period_paths
from hyperell_regulator.regulator.plane import period_matrix
from hyperell_regulator.regulator.polynomials import ring
from hyperell_regulator.regulator.precision import period_matrix_mp


def test_mst_period_edges_preserve_exact_branch_endpoints():
    C=ComplexField(240)
    roots=[C(-2),C(-1),C(0),C(1),C(3)]
    paths=period_paths(roots,2)
    assert len(paths)==4 and all(len(p)==2 for p in paths)
    assert all(z is roots[i] for p in paths for z,i in p)
    assert all(abs(p[0][0]-p[1][0]) in (1,2) for p in paths)


def test_even_degree_omits_a_leaf_edge():
    C=ComplexField(240);roots=[C(i) for i in (-3,-2,-1,0,1,2)]
    paths=period_paths(roots,2)
    assert len(paths)==4
    assert all(i!=0 for p in paths for z,i in p)  # farthest leaf -3 omitted


def test_mst_period_bases_are_integral_unimodular_changes():
    reference=json.loads((Path(__file__).parent/'data/period_basis_reference.json').read_text())
    R=RealField(180);C=ComplexField(180)
    real=lambda P:matrix(R,[[z.real() for z in rr]+[z.imag() for z in rr] for rr in P])
    for row in reference['cases']:
        old=[[C(R(z[0]),R(z[1])) for z in rr] for rr in row['periods']]
        new=period_matrix_mp(ring()(row['F']),30)
        change=real(new)*real(old).inverse()
        exact=matrix(ZZ,[[v.round() for v in rr] for rr in change.rows()])
        assert abs(exact.det())==1,row['F']
        assert max(abs(v) for v in (change-exact).list())<R('1e-27'),row['F']


@pytest.mark.parametrize('model', ['real', 'quintic', 'sextic', 'symmetric'])
def test_full_periods_and_minors_agree_across_precision(model):
    x = ring().gen()
    F = {
        'real': (x+2)*(x+1)*x*(x-1)*(x-3),
        'quintic': 4*x**5-7*x**4+4*x,
        'sextic': x**6+4*x**5+4*x**4+2*x**3+1,
        'symmetric': x**6-1,
    }[model]
    lower = period_matrix_mp(F, 20)
    lo = np.array(lower, complex)
    hi = period_matrix_mp(F, 35)
    assert np.max(abs(period_matrix(F, epsabs=1e-12, epsrel=1e-12)-lo)) < 1e-10
    assert np.max(abs(lo-np.array(hi, complex))) < 1e-14
    for row, refined in zip(lower, hi):
        for a, b in zip(row, refined):
            assert abs(a-b) < b.real().parent()('1e-20')
    for low, high in zip(c_invariants(F), c_invariants(F, prec=35)):
        assert abs(low-float(high)) < 1e-10
    # Full genus-two periods multiply a minor by four on fixed paths.
    from hyperell_regulator.regulator._periods import pivoted_minor
    R = hi[0][0].real().parent()
    for part, value in zip((lambda v:v.real(), lambda v:v.imag()),
                           c_invariants(F, prec=35)):
        half = [[part(v)/2 for v in row] for row in hi]
        assert abs(value-4*pivoted_minor(half, R('1e-38'))) < R('1e-32')


def test_c_invariants_integrates_period_matrix_once():
    from hyperell_regulator.regulator import plane
    from hyperell_regulator.regulator import precision

    F = ring().gen()**5 + 1
    expected = (float(plane.cplus(F)), float(plane.cminus(F)))
    with patch.object(precision, 'period_matrix_mp', wraps=precision.period_matrix_mp) as integrate:
        assert c_invariants(F) == expected
    integrate.assert_called_once_with(F, 15)


def test_float_period_tolerance_controls_internal_accuracy_and_retains_numpy_output():
    F = ring().gen()**5 + 1
    from hyperell_regulator.regulator import precision

    with patch.object(precision, 'period_matrix_mp', wraps=period_matrix_mp) as integrate:
        values = period_matrix(F, epsabs=0, epsrel=1e-20)
    integrate.assert_called_once_with(F, 20)
    assert values.shape == (4, 2) and values.dtype == np.complex128
    for tolerances in ((0, 0), (-1, 0), (np.nan, 1e-8), (np.inf, 1e-8)):
        with pytest.raises(ValueError, match='tolerances'):
            period_matrix(F, *tolerances)


def test_mst_period_minors_have_stable_signs_across_precision():
    x=ring().gen();F=x**6+2*x**5+3*x**4+4*x**3-x**2-2*x+1
    values=[c_invariants(F,prec=p) for p in (None,20,35,60,90)]
    for pair in values[:-1]:
        for low,high in zip(pair,values[-1]):
            assert abs(float(low)-float(high))<1e-10
