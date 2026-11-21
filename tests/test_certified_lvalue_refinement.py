"""Refinement checks the final absolute radius, including normalization."""
import json
from pathlib import Path

import pytest
from sage.all import PolynomialRing, QQ, RealBallField

from hyperell_regulator import HyperellipticCurveWithReduction, StandardLFunction
from hyperell_regulator import ball_lfunction as ball
from hyperell_regulator.regulator.balls import CertificationError, narrow


def backend(monkeypatch, accuracies):
    evaluations = []
    plans = []
    def run(args, **kwargs):
        if args[3] == 'info':
            plans.append(int(args[1]))
            return {'coefficient_cutoff': [str(int(args[1]))]}
        accuracy = accuracies[min(len(evaluations), len(accuracies)-1)]
        evaluations.append(list(args))
        center = 2**accuracy
        return dict(value=[str(center-1), str(center+1), str(-accuracy)],
                    working_bits=['1200'], root_number=['1'], flint_version=['test'])
    monkeypatch.setattr(ball, 'executable', lambda: str(Path(__file__)))
    monkeypatch.setattr(ball, '_run', run)
    return evaluations, plans


class FakeL:
    def __init__(self):
        self._cache = {}
        self._RT = PolynomialRing(QQ, 'T')
    def conductor(self):
        return 1
    def precompute_h1(self, primes, **kwargs):
        pass
    def euler_factor(self, p):
        return self._RT(1)


def test_measured_shortfall_increases_coefficient_budget(monkeypatch):
    evaluations, plans = backend(monkeypatch, [280, 350])
    L = FakeL()
    result = ball.value(L, s=2, digits=100)
    assert narrow(result, 100)
    assert len(evaluations) == 2 and plans[1]-plans[0] > 40
    assert L._certified_metadata['attempts'][0]['achieved_absolute_bits'] < 333


def test_refinement_can_continue_beyond_three_attempts(monkeypatch):
    evaluations, _ = backend(monkeypatch, [320, 320, 320, 350])
    L = FakeL()
    assert narrow(ball.value(L, s=2, digits=100), 100)
    assert len(evaluations) == 4
    assert len(L._certified_metadata['attempts']) == 4


def test_normalization_is_included_in_radius_check(monkeypatch):
    x = PolynomialRing(QQ, 'x').gen()
    L = StandardLFunction(HyperellipticCurveWithReduction(x*x+x, x**3+1))
    monkeypatch.setattr(L, 'precompute_h1', lambda primes, **kwargs: None)
    monkeypatch.setattr(L, 'euler_factor', lambda p: L._RT(1))
    evaluations, _ = backend(monkeypatch, [334, 350])
    # The first L(2) radius is below 10^-100, but multiplying by
    # sqrt(62001)/(4*pi^2) makes the standard L''(1) radius too wide.
    result = ball.value(L, s=1, derivative=2, digits=100)
    assert len(evaluations) == 2 and narrow(result, 100)


def test_failed_explicit_cutoff_records_accuracy_and_still_fails_closed(monkeypatch):
    evaluations, _ = backend(monkeypatch, [320])
    L = FakeL()
    with pytest.raises(CertificationError, match='last radius.*coefficient cutoff'):
        ball.value(L, s=2, digits=100, num_coeffs=200)
    assert len(evaluations) == 1
    assert not L._certified_metadata['certified']
    assert L._certified_metadata['attempts'][0]['coefficient_cutoff'] == 200


def test_experiment_preserves_failed_enclosure_diagnostics(monkeypatch, tmp_path):
    from hyperell_regulator.experiments import beilinson as b
    attempts = [dict(radius_display='1.5e-97', coefficient_cutoff=3183493)]
    class FailedL:
        def __init__(self, *args, **kwargs):
            self._certified_metadata = dict(certified=False, attempts=attempts)
        def value(self, *args, **kwargs):
            raise CertificationError('last radius 1.5e-97')
    monkeypatch.setattr(b, 'StandardLFunction', FailedL)
    path = tmp_path/'values.json'
    result = b.standard_lvalue_cached('997.a.997.1', None, digits=100,
                                     path=path, certified=True)
    assert not result['L_certified'] and result['L_second_derivative'] is None
    assert result['L_metadata']['attempts'] == attempts
    assert json.loads(path.read_text())['997.a.997.1:certified']['L_metadata']['attempts'] == attempts


def test_later_backend_failure_retains_completed_attempt(monkeypatch):
    evaluations, _ = backend(monkeypatch, [320])
    first_run = ball._run
    def run(args, **kwargs):
        if args[3] != 'info' and evaluations:
            raise CertificationError('native timeout')
        return first_run(args, **kwargs)
    monkeypatch.setattr(ball, '_run', run)
    L = FakeL()
    with pytest.raises(CertificationError, match='native timeout'):
        ball.value(L, s=2, digits=100)
    assert not L._certified_metadata['certified']
    assert len(L._certified_metadata['attempts']) == 1

