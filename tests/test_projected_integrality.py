"""Replay the recorded formulas through the public boundary and kernel API."""

import json
from pathlib import Path

import pytest
from sage.all import vector

from hyperell_regulator import ProjectedTorsionPacketIntegrality
from hyperell_regulator.components import curve_from_row
from hyperell_regulator import projected_integrality
from reports.solve_relation_integral_kernels import recorded_packet


REPORTS = Path(__file__).resolve().parents[1]/'reports'
AUDIT = json.loads((REPORTS/'relation_integrality_below_1000.json').read_text())
DECISIONS = {r['row']: r['unramified']
             for name in ('semistable_relation_boundaries', 'dyadic_relation_boundaries')
             for r in json.loads((REPORTS/(name+'.json')).read_text())}


@pytest.mark.parametrize('row', AUDIT, ids=[f"{r['index']}-{r['label']}" for r in AUDIT])
def test_recorded_formula(row):
    # The transported pair is represented by its exact source classes.
    if 'source_index' in row:
        source = AUDIT[row['source_index']-1]
    else:
        source = row
    test, generators = recorded_packet(source)
    expected = DECISIONS.get(row['index'], len(row['integral_basis']) == row['determinant_rank'])
    answers = [test.check(D)['q_integral'] for D in generators]
    assert None not in answers
    assert all(answers) == expected
    solution = test.kernel(generators)
    assert solution['complete']
    assert test.kernel()['complete']
    assert solution['divisors'].rank() == row['determinant_rank']-int(not expected)
    assert solution['abscissa_basis'].nrows() == solution['divisors'].rank()
    assert all(test.check(D)['q_integral'] is True for D in solution['divisors'])


def test_abscissa_constructor_and_cancellation():
    row = AUDIT[2]
    C = curve_from_row(dict(label=row['label'], eqn=row['equation']))
    test = ProjectedTorsionPacketIntegrality.from_abscissae(C, [3, -1, None, 0])
    good = test.abscissa_divisor([3, 0, 1, -4])
    bad = test.abscissa_divisor([1, -1, 0, 0])
    assert test.check(good)['q_integral'] is True
    assert test.check(bad)['q_integral'] is False
    assert test.kernel([good/3, bad])['basis'].rows() == [vector([1, 0])]
    assert test.kernel([])['basis'].nrows() == 0


def test_intermediate_dyadic_chart_distinguishes_the_two_groups():
    test, _ = recorded_packet(AUDIT[42])
    original = next(P for P in test.points if P[0] == 0)
    intermediate = next(P for P in test.points if P[0] == 1)
    result = test.check({intermediate: 1, original: -1})
    at_two = next(data for p, data in result['local'].items() if p.smallest_integer() == 2)
    assert at_two['method'] == 'intermediate_chart'
    assert at_two['q_integral'] is False
    assert test.kernel()['complete']


def test_missing_local_model_gives_kernel_bounds(monkeypatch):
    solve = projected_integrality.semistable_constraints

    def missing(local, points, groups):
        if local.k.characteristic() == 3:
            raise NotImplementedError('test an unavailable local graph')
        return solve(local, points, groups)

    monkeypatch.setattr(projected_integrality, 'semistable_constraints', missing)
    test, (mixed, balanced) = recorded_packet(AUDIT[2])
    assert test.check(mixed)['q_integral'] is None
    assert test.check(balanced)['q_integral'] is True
    assert not test.kernel([mixed, balanced])['complete']
    assert test.kernel([balanced])['complete']
