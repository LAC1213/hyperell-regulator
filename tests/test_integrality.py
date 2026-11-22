"""Exact kernels, conservative unknowns, and LMFDB arithmetic regressions."""

import itertools
from types import SimpleNamespace
from unittest.mock import patch

import pytest
from sage.all import Infinity, NumberField, PolynomialRing, QQ, ZZ, vector

from hyperell_regulator import HyperellipticCurveWithReduction, TorsionPacketIntegrality
from hyperell_regulator.components import curve_from_row, curve_points
from hyperell_regulator.integrality import torsion_packets
from hyperell_regulator.number_field_integrality import NumberFieldTorsionPacketIntegrality


def curve_and_points():
    x = PolynomialRing(QQ, 'x').gen()
    C = HyperellipticCurveWithReduction(-x*x-x, x**3+x*x+x+1)
    return C, [('infinity', 0), ('infinity', -1), (0, 0), (0, -1)]


def synthetic_test(labels):
    """Two transverse component partitions, using actual rational points."""
    C, points = curve_and_points()

    class Model:
        def bad_primes(self):
            return sorted(labels)

        def is_semistable_at(self, p):
            return True

        def reduce_point(self, p, P):
            label = labels[p][points.index(P)]
            return SimpleNamespace(is_smooth=label is not None,
                                   component=SimpleNamespace(label=label),
                                   reason='contracted chain')

    with patch('hyperell_regulator.integrality.affine_model', return_value=(Model(), points, 0)), \
         patch('hyperell_regulator.integrality._smooth_affine_components',
               side_effect=NotImplementedError('no original affine patch in synthetic data')):
        return TorsionPacketIntegrality(C, points)


def test_balanced_sum_when_no_single_pair_passes():
    t = synthetic_test({5: ['A', 'A', 'B', 'B'], 7: ['C', 'D', 'C', 'D']})
    for i, j in itertools.combinations(range(4), 2):
        D = [0]*4
        D[i], D[j] = 1, -1
        result = t.check(D)
        assert result['criterion_satisfied'] is False
        assert result['q_integral'] is None  # not a nonintegrality assertion
    result = t.check([1, -1, -1, 1])
    assert result['q_integral'] is True
    assert all(set(d['multidegrees'].values()) == {0} for d in result['local'].values())
    kernel = t.kernel()
    assert kernel['complete'] and kernel['basis'].nrows() == 1
    assert list(kernel['divisors'][0]) == [1, -1, -1, 1]


def test_integer_kernel_in_alpha_coefficients_is_saturated():
    t = synthetic_test({5: ['A', 'A', 'B', 'B'], 7: ['C', 'D', 'C', 'D']})
    generators = [[-2, 2, 0, 0], [-3, 0, 3, 0], [-5, 0, 0, 5]]
    kernel = t.kernel(generators)
    row = kernel['basis'][0]
    assert row.base_ring() is ZZ
    assert list(row) == [15, 10, -6]
    assert t.check(t.combine(generators, row))['q_integral'] is True
    assert list(kernel['divisors'][0]) == [-30, 30, 30, -30]
    # Rationally scaled input generators must clear row denominators exactly.
    half = [[QQ(a)/2 for a in D] for D in generators]
    assert t.kernel(half)['basis'] == kernel['basis']
    assert t.kernel([])['basis'].nrows() == 0


def test_unknown_support_cancels_before_reduction_and_can_be_avoided():
    t = synthetic_test({5: ['A', 'A', None, None]})
    assert not t.complete
    assert t.check([1, -1, 0, 0])['q_integral'] is True
    result = t.check([0, 0, 1, -1])
    assert result['criterion_satisfied'] is None and result['q_integral'] is None
    assert set(result['local'][5]['unresolved']) == {2, 3}
    assert t.kernel()['basis'].nrows() == 1
    D = [0, 0, 1, -1]
    assert t.check(t.combine([D, D], [1, -1]))['q_integral'] is True
    relations = t.kernel([D, D])
    assert relations['basis'].nrows() == 1
    assert relations['divisors'].is_zero()


def test_higher_rank_integer_kernel_contains_mixed_parity_vectors():
    t = synthetic_test({5: ['A', 'A', 'B', 'B']})
    generators = [[-2, 0, 2, 0], [-1, 0, 1, 0], [-1, 0, 0, 1]]
    kernel = t.kernel(generators)
    # Clearing denominators in the usual rational echelon basis can miss
    # this vector, even though twice it is in the generated lattice.
    assert kernel['basis'].nrows() == 2
    assert vector(ZZ, [-1, 1, 1]) in kernel['basis'].row_module()


def test_bad_two_never_invents_component_data():
    t = synthetic_test({2: ['A']*4})
    assert t.check([1, -1, 0, 0])['criterion_satisfied'] is None
    assert t.kernel()['basis'].nrows() == 0
    assert t.check([0]*4)['q_integral'] is True


def test_number_field_chart_follows_sections_through_a_nonnodal_centre():
    x = PolynomialRing(QQ, 'x').gen()
    C = HyperellipticCurveWithReduction(2*x**5-3*x**4+x**3+x**2-x, x.parent()(1))
    K = NumberField(x, 'a')
    points = [('infinity', 0), (QQ(1)/2, -QQ(1)/2), (0, 0)]
    test = NumberFieldTorsionPacketIntegrality(C, points, field=K,
                                               torsion_multiples=[1, 2, 10])
    local = test.local_data[K.primes_above(2)[0]]
    assert local['components'][0] == local['components'][1]
    assert len(local['paths'][0]) == 2 and not local['paths'][0][0]['smooth']
    assert test.check([1, -1, 0])['q_integral'] is True
    assert test.check([1, 0, -1])['q_integral'] is None


def test_nontorsion_undetermined_and_invalid_inputs():
    C, points = curve_and_points()
    with patch('hyperell_regulator.integrality.divisor_class_order', return_value=Infinity):
        with pytest.raises(ValueError, match='torsion packet'):
            TorsionPacketIntegrality(C, points)
    with patch('hyperell_regulator.integrality.divisor_class_order', return_value=None):
        with pytest.raises(NotImplementedError, match='undetermined'):
            TorsionPacketIntegrality(C, points)
        packets, unknown = torsion_packets(C, points[:2])
        assert len(packets) == 2 and len(unknown) == 1
    for invalid in ([], [points[0], points[0]], [('infinity', 2)], [(0, 2)]):
        with pytest.raises(ValueError):
            TorsionPacketIntegrality(C, invalid)
    t = synthetic_test({5: ['A']*4})
    for invalid in ([1], [1]*4, {(1, 2): 0}):
        with pytest.raises(ValueError):
            t.check(invalid)
    with pytest.raises(ValueError):
        t.combine([[1, -1, 0, 0]], [])
    assert t.check({t.points[0]: 1, t.points[1]: -1})['q_integral'] is True


def test_249_nonminimal_discriminant_retains_chain_uncertainty():
    row = dict(label='249.a.6723.1', cond=249,
               eqn=[[2, 3, 1, 1, 0, -1], [1, 0, 0, 1]],
               rat_pts=[[-1, 0, 1], [0, -2, 1], [0, 1, 1],
                        [1, -1, 0], [1, 0, 0], [3, -14, 1]])
    C = curve_from_row(row)
    t = TorsionPacketIntegrality(C, curve_points(C, row))
    assert t.orders == [1, 28, 28, 7, 7, 2]
    assert t.local_data[3]['components'] == ['Gamma_R', None, None, 'Gamma_R', 'Gamma_R', None]
    assert t.check([1, 0, 0, -1, 0, 0])['q_integral'] is True
    assert t.check([0, 1, -1, 0, 0, 0])['criterion_satisfied'] is None


def test_256_integral_affine_points_certify_pairs_at_wild_two():
    C = curve_from_row(dict(eqn=[[0, -1, 1, 1, -3, 2], [1]]))
    points = [(0, -1), (0, 0), (QQ(1)/2, -QQ(1)/2), (1, -1), ('infinity', 0), (1, 0)]
    t = TorsionPacketIntegrality(C, points)
    assert t.orders == [1, 5, 10, 10, 10, 10]
    assert t.local_data[2]['components'] == ['affine_0', 'affine_0', None, 'affine_0', None, 'affine_0']
    assert t.check([1, 0, 0, -1, 0, 0])['q_integral'] is True
    assert t.check([1, 0, -1, 0, 0, 0])['criterion_satisfied'] is None
    assert t.kernel()['basis'].nrows() == 3
    assert not t.complete


def test_affine_patch_distinguishes_components_and_rejects_singular_points():
    from hyperell_regulator.integrality import _smooth_affine_components
    x = PolynomialRing(QQ, 'x').gen()
    # At 2 this fibre is y(y+1)=0: equal x need not mean equal component.
    C = HyperellipticCurveWithReduction(2*x**5, 1)
    labels, reasons, _ = _smooth_affine_components(C, [(0, 0), (0, -1)], 2)
    assert labels[0] != labels[1] and not reasons
    C = HyperellipticCurveWithReduction(x**5-x)
    labels, reasons, _ = _smooth_affine_components(C, [(0, 0), (1, 0)], 2)
    assert labels[0] is not None and labels[1] is None
    assert 'singular' in reasons[1]
