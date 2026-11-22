"""Validated public APIs, topology, precision checks and experiment storage."""
import json
from unittest.mock import patch

import pytest
from sage.all import PolynomialRing, QQ, RealBallField, Matrix, ZZ

from hyperell_regulator.certified import encode_ball, decode_ball
from hyperell_regulator.experiments import beilinson as b
from hyperell_regulator.regulator import regulators_stokes, c_invariants, IntegrationMonitor, IntegrationLimitError
from hyperell_regulator.regulator.balls import narrow, CertificationError

R = PolynomialRing(QQ, 'x')
x = R.gen()


def test_parallel_slits_agree_with_independent_mst_and_bound_all_entries():
    F = x**5+(x+1)**2
    details = {}
    per, cols = regulators_stokes(F, [0, -1], prec=10, certified=True,
                                  certified_cut='parallel', details=details)
    expected_per, expected_cols = regulators_stokes(F, [0, -1], prec=20)
    assert details['certified'] and details['cut_kind'] == 'parallel rays'
    for actual, expected in zip(per+sum(cols, []), expected_per+sum(expected_cols, [])):
        assert narrow(actual, 10)
        assert abs(actual-actual.parent()(expected)) < QQ(1)/10**16


def test_all_real_roots_weierstrass_centres_and_finite_pole():
    F = x*(x*x-1)*(x*x-4)
    per, cols = regulators_stokes(F, [0, 1, 0], pole=0, prec=8, certified=True)
    assert all(v.is_exact() and v == 0 for v in cols[0]+cols[2])
    assert all(narrow(v, 8) for v in per+sum(cols, []))
    # Reflection gives the middle unweighted moment exactly zero.
    assert per[1].contains_zero()
    old_per, old = regulators_stokes(F, [1], pole=0, prec=16)
    for actual, expected in zip(cols[1], old[0]):
        assert abs(actual-actual.parent()(expected)) < QQ(1)/10**10


@pytest.mark.parametrize('F', [x**5+(x+1)**2, x**6+2*x**3+4*x*x+4*x+1])
def test_signed_period_minors_match_existing_basis(F):
    actual = c_invariants(F, prec=12, certified=True)
    expected = c_invariants(F, prec=20)
    for a, e in zip(actual, expected):
        assert narrow(a, 12) and abs(a-a.parent()(e)) < QQ(1)/10**17


def test_certified_integrator_obeys_work_budget():
    with pytest.raises(IntegrationLimitError):
        with IntegrationMonitor(max_edge_steps=1):
            regulators_stokes(x**5+x+1, [], prec=8, certified=True)


def test_floating_inputs_are_coerced_exactly_without_rational_reconstruction():
    from sage.all import RealField
    from hyperell_regulator.regulator.balls import rational, _model
    n, d = (0.1).as_integer_ratio()
    assert rational(0.1) == QQ(n)/d != QQ(1)/10
    assert rational(RealField(53)(0.1)) == QQ(n)/d
    RF = PolynomialRing(RealField(53), 'x')
    f = RF.gen()**5+RF(0.1)*RF.gen()+1
    assert _model(f)[1] == QQ(n)/d
    from hyperell_regulator.regulator import balls
    BF = RealBallField(100)
    with patch.object(balls, 'compute', return_value=([BF(1)], [[BF(2)], [BF(3)]])) as compute:
        _, cols = regulators_stokes(x**5+x+1, [0.1], pole=QQ(1)/10,
                                     prec=10, certified=True)
    assert len(compute.call_args.args[1]) == 2
    assert cols[0][0] == -1


def test_lossless_interval_storage_and_safe_comparison():
    F = RealBallField(450)
    a = F.pi().add_error(F(10)**-110)
    record = json.loads(json.dumps(encode_ball(a)))
    assert a in decode_ball(record)
    out = b.beilinson_comparison(F(1), F(1), F(1), F(1)/F.pi(),
                                 digits=100, certified=True)
    assert 1 in out['ratio']
    assert out['exact_relation_proved'] is False
    assert out['relative_agreement'].upper() < QQ(1)/10**100
    assert b.beilinson_comparison(F(0).add_error(1), F(1), F(1), F(1),
                                 certified=True)['comparison_error']
    with pytest.raises(CertificationError, match='requires balls'):
        b.beilinson_comparison(1., F(1), F(1), F(1), certified=True)


def test_experiment_cli_propagates_certified():
    with patch.object(b, 'run') as run:
        b.main(['250', '--certified', '--digits', '100', '--reg-digits', '100'])
    assert run.call_args.kwargs['certified'] is True
    assert run.call_args.kwargs['out'].endswith('_certified.json')
    assert run.call_args.kwargs['digits'] == run.call_args.kwargs['reg_digits'] == 100


def test_report_cache_rejects_numerical_records_and_wide_balls():
    F = RealBallField(150)
    settings = dict(certified=True, regulator_digits=15, lvalue_digits=15, with_lvalue=False)
    old = dict(status='ok', certified=True, regulator_cut='tree', det=encode_ball(F.pi()), settings=settings)
    assert b._stored_result_satisfies(old, settings)
    assert not b._stored_result_satisfies(dict(old,regulator_cut='parallel_rays'),settings)
    assert not b._stored_result_satisfies(old, dict(settings, certified=False))
    assert not b._stored_result_satisfies(dict(old, det=encode_ball(F(1).add_error(1))), settings)


def test_certified_cache_is_separate_and_revalidates_radius(tmp_path):
    F = RealBallField(200)
    path = tmp_path/'values.json'
    path.write_text(json.dumps({'example': {'L_second_derivative': '999', 'L_digits': 1000}}))
    class LF:
        def __init__(self, *args, **kwargs):
            self._certified_metadata = {'coefficient_cutoff': 100}
        def conductor(self):
            return 17
        def value(self, *args, **kwargs):
            assert kwargs['certified']
            return F.pi()
    with patch.object(b, 'StandardLFunction', LF):
        first = b.standard_lvalue_cached('example', None, digits=15, path=path, certified=True)
        assert first['L_certified'] and not first['L_cached']
        second = b.standard_lvalue_cached('example', None, digits=15, path=path, certified=True)
        assert second['L_cached']
        records = json.loads(path.read_text())
        assert records['example']['L_second_derivative'] == '999'
        records['example:certified']['L_second_derivative'] = encode_ball(F(1).add_error(1))
        path.write_text(json.dumps(records))
        assert not b.standard_lvalue_cached('example', None, digits=15, path=path,
                                           certified=True)['L_cached']


def test_newton_lift_matches_exact_characteristic_polynomial():
    from hyperell_regulator.frobenius import _lift
    M = Matrix(ZZ, [[12, 3, -4, 8], [1, 6, 7, 9], [4, 8, 14, 0], [9, 5, 3, -1]])
    p, N = 101, 3
    coefficients = M.charpoly().list()[2:]
    q = p**N
    a = [int(c)%q for c in coefficients]
    a = [c if 2*c < q else c-q for c in a]
    expected = R([p*p*a[2], p*a[1]]+a)
    assert _lift(M.rows(), p, 2, N) == expected


def test_certified_lvalue_fails_closed_when_backend_missing(monkeypatch):
    from hyperell_regulator import HyperellipticCurveWithReduction, StandardLFunction
    L = StandardLFunction(HyperellipticCurveWithReduction(x*x+x, x**3+1))
    monkeypatch.setenv('HYPERELL_CERTIFIED_LVALUE', '/missing/certified-lvalue')
    with pytest.raises(CertificationError, match='backend unavailable'):
        L.value(digits=10, certified=True)
    with pytest.raises(NotImplementedError, match='binary64'):
        L.value(QQ(1)/3, derivative=0, certified=True)
    assert L.value(derivative=0, certified=True).is_zero()
