"""Explicit bad factors, conductor exponents, and cache invalidation."""

import json
from unittest.mock import patch

import pytest
from sage.all import PolynomialRing, QQ, RealBallField, RealField

from hyperell_regulator import HyperellipticCurveWithReduction
from hyperell_regulator.experiments import beilinson as b
from hyperell_regulator.experiments.beilinson_local_data import (
    LOCAL_DATA, lfunction_options, local_data_metadata,
)
from hyperell_regulator.lfunction import StandardLFunction


def curve(label):
    R = PolynomialRing(QQ, 'x')
    f, h = LOCAL_DATA[label]['model']
    return HyperellipticCurveWithReduction(R(f), R(h))


@pytest.mark.parametrize('label,conductor', [
    ('324.a.648.1', 26244), ('576.a.576.1', 36864), ('784.c.614656.1', 9604),
    ('688.a.2752.1', 29584), ('688.a.704512.1', 29584), ('688.a.704512.2', 29584),
    ('1296.a.20736.1', 26244), ('2704.a.43264.1', 114244),
    ('3969.c.35721.1', 15752961), ('4096.b.65536.1', 4096),
    ('8649.c.700569.1', 74805201),
    ('3125.a.3125.1', 15625),
])
def test_manual_factors_and_full_standard_conductors(label, conductor):
    C = curve(label)
    h1_conductor = 576 if label == '576.a.576.1' else None
    L = StandardLFunction(C, h1_conductor=h1_conductor, **lfunction_options(label, C))
    assert L.conductor() == conductor
    for p, a in LOCAL_DATA[label]['standard_factors'].items():
        assert L.euler_factor(int(p)) == L._RT(a)


def test_manual_data_are_guarded_by_the_curve_equation():
    with pytest.raises(ValueError, match='require its recorded model'):
        lfunction_options('324.a.648.1', curve('4096.b.65536.1'))
    assert lfunction_options('unlisted', None) == {}
    assert local_data_metadata('unlisted') is None


@pytest.mark.parametrize('overrides', [{4: 2}, {2: -1}, {2: QQ(3)/2}, {3: 1}])
def test_conductor_override_validation(overrides):
    with pytest.raises((ValueError, TypeError)):
        StandardLFunction(curve('4096.b.65536.1'), conductor_exponents=overrides)


def test_wild_conductor_is_independent_of_euler_polynomial_degree():
    label = '324.a.648.1'
    C = curve(label)
    L = StandardLFunction(C, **lfunction_options(label, C))
    assert 5-L.euler_factor(3).degree() == 4
    assert L.conductor_exponents[3] == 8
    assert L.check_functional_equation(digits=12) < 1e-18


@pytest.mark.parametrize('certified', [False, True])
def test_both_lvalue_paths_use_manual_data_and_invalidate_stale_caches(tmp_path, certified):
    label = '4096.b.65536.1'
    C = curve(label)
    options = lfunction_options(label, C)
    metadata = local_data_metadata(label)
    path = tmp_path/'values.json'
    constructions = []
    evaluations = []

    class LF:
        def __init__(self, C, **kwargs):
            constructions.append(kwargs)
            self._certified_metadata = {'coefficient_cutoff': 100}

        def conductor(self):
            return 4096

        def _coeffs_needed(self, digits):
            return 100

        def lfun(self, *args):
            return None

        def check_functional_equation(self, *args):
            return 1e-30

        def value(self, *args, **kwargs):
            evaluations.append(kwargs)
            return RealBallField(200).pi() if kwargs.get('certified') else RealField(200).pi()

    with patch.object(b, 'StandardLFunction', LF):
        first = b.standard_lvalue_cached(label, C, digits=20, path=path, certified=certified)
        assert first.get('L_status') is None
        assert not first['L_cached'] and first['L_local_data'] == metadata
        second = b.standard_lvalue_cached(label, C, digits=20, path=path, certified=certified)
        assert second['L_cached'] and len(evaluations) == 1
        saved = json.loads(path.read_text())
        key = label + (':certified' if certified else '')
        # Even with the same conductor, old Euler polynomials cannot be reused.
        del saved[key]['L_local_data']
        path.write_text(json.dumps(saved))
        third = b.standard_lvalue_cached(label, C, digits=20, path=path, certified=certified)
        assert not third['L_cached'] and len(evaluations) == 2
    assert all(kwargs['bad_factors'] == options['bad_factors'] and
               kwargs['conductor_exponents'] == options['conductor_exponents']
               for kwargs in constructions)


def test_report_cache_revalidates_local_data():
    label = '4096.b.65536.1'
    settings = dict(with_lvalue=True, regulator_digits=20, lvalue_digits=20)
    record = dict(label=label, status='ok', settings=settings,
                  L_local_data=local_data_metadata(label), L_second_derivative='1')
    assert b._stored_result_satisfies(record, settings)
    stale = dict(record, L_local_data=None)
    assert not b._stored_result_satisfies(stale, settings)
