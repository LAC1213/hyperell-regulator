"""Semilunar topology, sheet continuation, and full cycle normalisation."""

from unittest.mock import patch

import mpmath as mp
import numpy as np
import pytest
from sage.all import ComplexField

from hyperell_regulator.regulator import c_invariants
from hyperell_regulator.regulator import _paths
from hyperell_regulator.regulator._periods import period_paths
from hyperell_regulator.regulator.cover import Cover, check_path_embedded
from hyperell_regulator.regulator.plane import period_matrix
from hyperell_regulator.regulator.polynomials import ring
from hyperell_regulator.regulator.precision import period_matrix_mp


def test_period_and_regulator_paths_share_validation_and_detours():
    roots = [-2, -1, 0, 1, 3]
    with patch.object(_paths, 'check_path_embedded', wraps=check_path_embedded) as check:
        paths = period_paths(roots, 2)
        assert check.called
        assert all(call.kwargs.get('include_ray') is False for call in check.call_args_list)
    assert paths[0][0][0] == 3
    assert len(paths[0]) > 2  # the closing chord would run through other roots
    nodes = [node for node, _ in paths[0]]
    nodes += [path[-1][0] for path in paths[1:]]
    check_path_embedded(nodes, crit=roots, include_ray=False)

    x = ring().gen()
    cov = Cover(x**5 - 5*x**3 + 4*x, x, 0)
    with patch.object(_paths, 'check_path_embedded', wraps=check_path_embedded) as check:
        nodes, kinds, ray = cov.cut_path()
        assert check.called
    assert 'none' in kinds
    check_path_embedded(nodes, ray, cov.critical_values())


def test_largest_modulus_start_and_exact_root_endpoints():
    C = ComplexField(240)
    roots = [C(-5), C(-1), C(0), C(1), C(2)]
    paths = period_paths(roots, 2)
    assert paths[0][0][0] == -5  # not the root with largest real part
    assert all(node is roots[label] for path in paths for node, label in path
               if label is not None)
    roots = [C(-2), C(-1), C(0), C(1), C(2)]
    assert period_paths(roots, 2)[0][0][0] == 2  # real-part tie-breaker


@pytest.mark.parametrize('nodes,crit', [
    ([0, 2], [0, 1, 2]),
    ([0, 1+1j, 1j, 1], [0, 1+1j, 1j, 1]),
    ([0, 1, 2, .5], [0, 1, 2, .5]),
])
def test_finite_path_validation_rejects_obstructions(nodes, crit):
    with pytest.raises(NotImplementedError):
        check_path_embedded(nodes, crit=crit, include_ray=False)


def test_real_branch_periods_match_independent_integrals_and_integral_basis():
    x = ring().gen()
    F = (x+2)*(x+1)*x*(x-1)*(x-3)
    periods = np.array(period_matrix_mp(F, 30), complex)
    roots = [-2, -1, 0, 1, 3]
    # Independent direct x-integrals along the ascending real chain.
    with mp.workdps(75):
        reference = []
        for a, b in zip(roots[:-1], roots[1:]):
            def f(t, k):
                return t**k / mp.sqrt(mp.fprod(t-r for r in roots))
            reference.append([complex(2*mp.quad(lambda t: f(t, k), [a, b]))
                              for k in range(2)])
    reference = np.array(reference)
    # After the detour, the remaining rows are the first three real segments.
    assert np.max(abs(periods[1:] - reference[:3])) < 1e-12
    # The detour must continue on the same sheet: the resulting cycles are
    # an integral unimodular change of the independently integrated basis.
    def real_matrix(P):
        return np.concatenate([P.real, P.imag], axis=1)
    change = real_matrix(periods) @ np.linalg.inv(real_matrix(reference))
    integral = np.rint(change)
    assert np.max(abs(change-integral)) < 1e-12
    assert abs(round(np.linalg.det(integral))) == 1


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


def test_277_period_pivot_sign_is_stable_under_quadrature_refinement():
    x = ring().gen()
    F = x**6+2*x**5+3*x**4+4*x**3-x**2-2*x+1
    values = [c_invariants(F, prec=digits) for digits in (None, 30, 50, 90)]
    # Conjugate rows have equal real pivots. A Newton-sized tie tolerance
    # previously made c_plus alternate signs at 30, 50, 70 and 90 digits.
    assert all(cp < 0 and cm > 0 for cp, cm in values)
    for pair in values[:-1]:
        for low, high in zip(pair, values[-1]):
            assert abs(float(low)-float(high)) < 1e-10
