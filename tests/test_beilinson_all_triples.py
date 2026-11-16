"""Exhaustive torsion triples without component checks, up to permutations."""
import itertools
import json
from types import SimpleNamespace
from unittest.mock import patch

from hyperell_regulator.experiments import beilinson as b
from hyperell_regulator.regulator.polynomials import ring


def test_cli_enables_exhaustive_mode():
    with patch.object(b, 'run', return_value=[]) as run:
        b.main(['--ignore-component-criterion'])
        assert run.call_args.kwargs['ignore_component_criterion'] is True
        b.main([])
        assert run.call_args.kwargs['ignore_component_criterion'] is False


def test_failed_components_recover_orders_and_keep_zero_filters():
    x = ring().gen()
    C = SimpleNamespace(F=x**6+x+1, h=ring()(0), genus=lambda: 2)
    points = [(0, 1), (0, -1), (-1, 1), (-1, -1),
              ('infinity', 1), ('infinity', -1)]
    pairs = [dict(P=P, Q=Q, same_component_everywhere=False,
                  order_of_P_minus_Q=None)
             for P, Q in itertools.combinations(points, 2)]
    # Missing orders on odd-degree models must also be recovered in this mode.
    with patch.object(b, 'odd_degree_model', return_value=object()), \
            patch.object(b, 'divisor_class_order', return_value=5) as order:
        assert b.find_triples(dict(pairs=pairs), C, row={'mw_invs': [5]}) == []
        order.assert_not_called()
        triples = b.find_triples(dict(pairs=pairs), C, row={'mw_invs': [5]},
                                 allow_undecided=False, ignore_component_criterion=True)
    assert order.call_count == 12  # Conjugate pairs never need a torsion test.
    assert len(triples) == 24
    assert len({frozenset(t[k] for k in ('P', 'Q', 'R')) for t in triples}) == 8
    for t in triples:
        assert len({t[k][0] for k in ('P', 'Q', 'R')}) == 3
        assert t['order_P'] == t['order_Q'] == 5
        assert t['verdict_P'] is t['verdict_Q'] is False


def test_exhaustive_mode_still_rejects_nontorsion_and_uncertified_pairs():
    x = ring().gen()
    C = SimpleNamespace(F=x**6+x+1, h=ring()(0), genus=lambda: 2)
    pairs = [dict(P=P, Q=Q, same_component_everywhere=False,
                  order_of_P_minus_Q=order)
             for (P, Q), order in zip(itertools.combinations([(0, 1), (1, 1), (2, 1)], 2),
                                      (5, 'Infinity', None))]
    with patch.object(b, 'odd_degree_model', return_value=None), \
            patch.object(b, 'divisor_class_order', return_value=None):
        assert b.find_triples(dict(pairs=pairs), C, row={'mw_invs': [0, 5]},
                              ignore_component_criterion=True) == []


def test_records_skip_component_checks_and_keep_cached_metadata(tmp_path):
    row = dict(label='test', cond=1, is_simple_base=True, is_simple_geom=True)
    P, Q, R = (0, 1), (1, 1), ('infinity', 1)
    known = dict(label='test', cond=1, pairs=[dict(
        P=P, Q=Q, same_component_everywhere=False, order_of_P_minus_Q=5)])
    cache, records = tmp_path/'curves.json', tmp_path/'components.json'
    cache.write_text(json.dumps([row]))
    records.write_text(json.dumps([known]))
    before = records.read_text()
    with patch.object(b, 'curve_from_row'), \
            patch('hyperell_regulator.components.curve_points', return_value=[P, Q, R, P, None]), \
            patch('hyperell_regulator.components.report_curve') as report, \
            patch('hyperell_regulator.components.same_component_all_primes') as check:
        prepared = b.ensure_records(1000, cache=cache, results=records,
                                    ignore_component_criterion=True)
        assert records.read_text() == before
        records.unlink()
        fresh = b.ensure_records(1000, cache=cache, results=records,
                                 ignore_component_criterion=True)
    report.assert_not_called()
    check.assert_not_called()
    assert len(prepared[0]['pairs']) == len(fresh[0]['pairs']) == 3
    pair = prepared[0]['pairs'][0]
    assert pair['same_component_everywhere'] is False
    assert pair['order_of_P_minus_Q'] == 5
    assert all(p['same_component_everywhere'] is None for p in fresh[0]['pairs'])
    assert not records.exists()


def test_all_unordered_triples_are_evaluated_and_failures_do_not_stop_the_run(tmp_path, capsys):
    row = dict(label='test', cond=1, is_simple_geom=True, mw_invs=[5])
    cache = tmp_path/'curves.json'
    cache.write_text(json.dumps([row]))
    points = [(0, 1), (1, 1), (2, 1), (3, 1)]
    triples = []
    for group in itertools.combinations(points, 3):
        for P, Q, R in itertools.permutations(group):
            triples.append(dict(P=P, Q=Q, R=R, order_P=5, order_Q=5))
    # The first ranked orientation of each unordered set is the representative.
    ranked = list(reversed(triples))
    expected = {}
    for t in ranked:
        expected.setdefault(frozenset(t[k] for k in ('P', 'Q', 'R')), t)
    calls = []

    def calculate(record, C, t, n_triples, **kwargs):
        calls.append(t)
        assert n_triples == 4
        if len(calls) == 2:
            raise ValueError('bad model')
        return dict(record, triple=t, status='numerical_failure' if len(calls) == 3 else 'ok')

    with patch.object(b, 'curve_from_row'), \
            patch.object(b, 'find_triples', return_value=triples) as find, \
            patch.object(b, 'rank_triples', return_value=ranked) as rank, \
            patch.object(b, '_analyse_triple', side_effect=calculate):
        result = b.analyse(dict(label='test', cond=1), cache=cache, verbose=False,
                           tries=1, ignore_component_criterion=True)
    assert find.call_args.kwargs['ignore_component_criterion'] is True
    assert rank.call_args.kwargs['usable_only'] is False
    assert calls == list(expected.values())
    assert result['status'] == 'ok'
    assert result['n_triples'] == result['n_triples_attempted'] == 4
    assert [r['status'] for r in result['triple_results']] == [
        'ok', 'preparation_error', 'numerical_failure', 'ok']
    assert result['triple_results'][1]['error'] == 'ValueError: bad model'
    # Exercise printing after a JSON round trip, as when reusing a saved report.
    b._report_row(json.loads(json.dumps(result)), stored=True)
    assert capsys.readouterr().out.count('  P ') == 4


def test_mode_is_forwarded_and_saved_and_invalidates_other_cached_policies(tmp_path):
    record = dict(label='test', cond=1)
    triple = dict(P=(0, 1), Q=(1, 1), R=(2, 1))
    outcome = dict(record, status='ok', triple=triple, det=1)
    result = dict(outcome, triple_results=[outcome], n_triples=1, n_triples_attempted=1)
    out = tmp_path/'results.json'
    with patch.object(b, 'ensure_records', return_value=[record]) as prepare, \
            patch.object(b, 'analyse', return_value=result) as analyse:
        b.run(out=out, with_lvalue=False, ignore_component_criterion=True)
        assert prepare.call_args.kwargs['ignore_component_criterion'] is True
        assert analyse.call_args.kwargs['ignore_component_criterion'] is True
        b.run(out=out, with_lvalue=False, ignore_component_criterion=True)
        assert analyse.call_count == 1
        b.run(out=out, with_lvalue=False)
        assert analyse.call_count == 2
    saved = json.loads(out.read_text())[0]
    settings = dict(saved['settings'], ignore_component_criterion=True)
    good = dict(saved, settings=settings)
    assert b._stored_result_satisfies(good, settings)
    failed = dict(good, triple_results=[outcome, dict(outcome, status='numerical_failure')])
    assert not b._stored_result_satisfies(failed, settings)
    assert not b._stored_result_satisfies(saved, settings)
    lvalue_settings = dict(settings, with_lvalue=True)
    failed_lvalue = dict(good, settings=lvalue_settings, triple_results=[
        outcome, dict(outcome, L_status='failed', L_second_derivative=None)])
    assert not b._stored_result_satisfies(failed_lvalue, lvalue_settings)
