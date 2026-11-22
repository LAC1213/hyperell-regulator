"""Load exact LMFDB labels and run one divisor-basis experiment at a time."""

import io
import json
from urllib.parse import parse_qs, urlparse
from unittest.mock import patch

import pytest
from sage.all import PolynomialRing, QQ

from hyperell_regulator import HyperellipticCurveWithReduction, load_curve, beilinson_experiment
from hyperell_regulator import components
from hyperell_regulator.experiments import beilinson as b


ROW = dict(label='394.a.394.1', cond=394, end_alg='Q', is_simple_base=True,
           is_simple_geom=True, eqn=[[-9, 17, 0, -12, 1, 2], [0, 1, 0, 1]],
           rat_pts=[[1, -1, 0], [1, -1, 1], [1, 0, 0]], mw_invs=[10])


def api_response(rows):
    return io.BytesIO(json.dumps(dict(data=rows)).encode())


def test_cached_label_returns_curve_with_reduction_and_metadata(tmp_path):
    cache = tmp_path/'curves.json'
    cache.write_text(json.dumps([ROW]))
    with patch('urllib.request.urlopen') as network:
        C = load_curve(ROW['label'], cache=cache, fetch=False, prec=180)
    network.assert_not_called()
    assert C.genus() == 2 and C.prec == 180
    assert C.bad_primes() == [2, 197]
    assert C.lmfdb_label == ROW['label'] and C.lmfdb_record == ROW


def test_fetch_queries_only_requested_label_and_preserves_other_cached_curves(tmp_path):
    cache = tmp_path/'curves.json'
    other = dict(ROW, label='other', cond=1)
    cache.write_text(json.dumps([other]))
    curve_row = {key: value for key, value in ROW.items() if key not in ('rat_pts', 'mw_invs')}
    points = dict(label=ROW['label'], rat_pts=ROW['rat_pts'], mw_invs=[10],
                  num_rat_pts=3, rat_pts_v=True)
    with patch('urllib.request.urlopen', side_effect=[api_response([curve_row]), api_response([points])]) as network:
        C = load_curve(ROW['label'], cache=cache)
    assert network.call_count == 2
    for call, table in zip(network.call_args_list, ('g2c_curves', 'g2c_ratpts')):
        request = call.args[0]
        assert urlparse(request.full_url).path == '/api/'+table+'/'
        query = parse_qs(urlparse(request.full_url).query)
        assert query['label'] == ['s'+ROW['label']]
        assert '_offset' not in query
    saved = {row['label']: row for row in json.loads(cache.read_text())}
    assert saved['other'] == other
    assert C.lmfdb_record == saved[ROW['label']]
    assert saved[ROW['label']]['rat_pts_verified'] is True
    assert saved[ROW['label']]['mw_invs'] == [10]
    with patch('urllib.request.urlopen') as network:
        assert load_curve(ROW['label'], cache=cache).lmfdb_record == C.lmfdb_record
    network.assert_not_called()


@pytest.mark.parametrize('label', ['394.a', '394.a.394', '0.a.1.1', '../394.a.394.1', None])
def test_invalid_labels_are_rejected_before_network_or_cache_writes(tmp_path, label):
    with patch('urllib.request.urlopen') as network, pytest.raises(ValueError, match='label'):
        load_curve(label, cache=tmp_path/'absent.json')
    network.assert_not_called()
    assert not list(tmp_path.iterdir())


def test_unknown_label_and_captcha_do_not_create_cache_entries(tmp_path):
    cache = tmp_path/'curves.json'
    with patch('urllib.request.urlopen', return_value=api_response([])) as network:
        with pytest.raises(ValueError, match='not found'):
            load_curve(ROW['label'], cache=cache)
    assert network.call_count == 1 and not cache.exists()
    with patch('urllib.request.urlopen', return_value=io.BytesIO(b'<html>challenge</html>')):
        with pytest.raises(RuntimeError, match='HTML'):
            load_curve(ROW['label'], cache=cache)
    assert not cache.exists()
    with patch('urllib.request.urlopen') as network:
        with pytest.raises(ValueError, match='without fetching'):
            load_curve(ROW['label'], cache=cache, fetch=False)
    network.assert_not_called()


def test_failed_refresh_preserves_existing_cache(tmp_path):
    cache = tmp_path/'curves.json'
    cache.write_text(json.dumps([ROW]))
    before = cache.read_bytes()
    with patch('urllib.request.urlopen', side_effect=[api_response([ROW]), OSError('offline')]):
        with pytest.raises(OSError, match='offline'):
            load_curve(ROW['label'], cache=cache, refresh=True)
    assert cache.read_bytes() == before
    with patch('urllib.request.urlopen', return_value=api_response([dict(ROW, label='395.a.395.1')])):
        with pytest.raises(ValueError, match='unexpected curve'):
            load_curve(ROW['label'], cache=cache, refresh=True)
    assert cache.read_bytes() == before


def test_single_curve_runs_real_integral_kernel_and_saves_one_result(tmp_path, monkeypatch):
    monkeypatch.setenv('HYPERELL_REGULATOR_DATA', str(tmp_path))
    cache, out = tmp_path/'curves.json', tmp_path/'experiment.json'
    cache.write_text(json.dumps([ROW]))
    C = load_curve(ROW['label'], cache=cache, fetch=False)
    with patch('urllib.request.urlopen') as network:
        result = beilinson_experiment(C, packet='rational', reg_digits=20,
                                      with_lvalue=False, out=out)
    network.assert_not_called()
    assert result['label'] == ROW['label']
    assert result['status'] == 'insufficient_integral_rank'
    assert result['integral_rank'] == 1
    assert len(result['integral_kernels'][0]['basis']) == 1
    assert result['settings']['packet'] == 'rational'
    assert result['settings']['regulator_digits'] == 20
    assert result['edges_attempted'] == 0 and result['seconds'] >= 0
    assert json.loads(out.read_text()) == [result]


def test_custom_polynomial_curves_need_no_lmfdb_cache_and_have_distinct_identities(tmp_path):
    x = PolynomialRing(QQ, 'x').gen()
    curves = [HyperellipticCurveWithReduction(x**5+(x+1)**2),
              HyperellipticCurveWithReduction(x**5+x+1)]
    seen = []

    def experiment(row, curve, metadata, selection, **options):
        assert curve is curves[len(seen)]
        assert metadata['eqn'] == [list(map(str, curve.f)), list(map(str, curve.h))]
        seen.append(row['label'])
        return dict(label=row['label'], cond=row['cond'], status='insufficient_integral_rank')

    with patch('hyperell_regulator.experiments._beilinson_divisors.analyse_divisors', side_effect=experiment), \
            patch('urllib.request.urlopen') as network:
        for C in curves:
            b.run_curve(C, cache=tmp_path/'does-not-exist.json', with_lvalue=False)
    network.assert_not_called()
    assert len(set(seen)) == 2 and all(label.startswith('custom.') for label in seen)


def test_explicit_curve_is_not_filtered_by_survey_endomorphism_policy(tmp_path):
    C = components.curve_from_row(dict(ROW, label='custom', end_alg='RM'))
    with patch('hyperell_regulator.experiments._beilinson_divisors.analyse_divisors',
               return_value=dict(status='insufficient_integral_rank')) as experiment:
        result = b.run_curve(C, with_lvalue=False)
    assert result['status'] == 'insufficient_integral_rank'
    experiment.assert_called_once()


def test_single_label_above_survey_bound_uses_exact_curve_and_options(tmp_path):
    row = dict(ROW, label='1394.a.1394.1', cond=1394)
    cache = tmp_path/'curves.json'
    cache.write_text(json.dumps([ROW, row]))
    with patch.object(b, 'analyse', return_value=dict(status='ok', det='2')) as analyse:
        result = b.run_curve(row['label'], cache=cache, fetch=False,
                             certified=True, digits=20, with_lvalue=False)
    assert result['label'] == row['label']
    assert analyse.call_args.args[0].lmfdb_label == row['label']
    assert analyse.call_args.kwargs['certified']
    assert analyse.call_args.kwargs['reg_digits'] == 20
    assert analyse.call_args.kwargs['end_algs'] is None


def test_failure_is_recorded_and_single_curve_success_cache_is_reused(tmp_path):
    C = components.curve_from_row(ROW)
    out = tmp_path/'experiment.json'
    with patch.object(b, 'analyse', side_effect=RuntimeError('computation failed')):
        failed = b.run_curve(C, with_lvalue=False, out=out)
    assert failed['status'] == 'preparation_error' and 'computation failed' in failed['error']
    with patch.object(b, 'analyse', return_value=dict(status='ok', det='2')) as analyse:
        successful = b.run_curve(C, with_lvalue=False, out=out)
        assert analyse.call_count == 1
        assert b.run_curve(C, with_lvalue=False, out=out) == successful
        assert analyse.call_count == 1
        b.run_curve(C, with_lvalue=False, out=out, ignore_integrality=True)
        assert analyse.call_count == 2
    assert len(json.loads(out.read_text())) == 1


def test_cli_defaults_and_single_label_dispatch():
    with patch.object(b, 'run', return_value=[]) as run:
        b.main([])
        assert run.call_args.args == (1000,)
        assert run.call_args.kwargs['packet'] == 'geometric'
        assert run.call_args.kwargs['reg_digits'] is None
        assert run.call_args.kwargs['ignore_integrality'] is False
        b.main(['--ignore-integrality'])
        assert run.call_args.kwargs['ignore_integrality'] is True
    with patch.object(b, 'run') as survey, patch.object(b, 'run_curve', return_value={}) as one:
        b.main(['--label', ROW['label'], '--certified', '--digits', '20', '--fetch'])
        survey.assert_not_called()
        assert one.call_args.args == (ROW['label'],)
        assert one.call_args.kwargs['refresh'] and one.call_args.kwargs['certified']
        assert one.call_args.kwargs['out'].endswith('beilinson_394.a.394.1_certified.json')


@pytest.mark.parametrize('args', [['500', '--label', ROW['label']], ['--packet', 'legacy'],
                                  ['--exclude-undecided'], ['--digits', '0'], ['--edge-seconds', 'nan'],
                                  ['--curve-seconds', '-1']])
def test_cli_rejects_conflicting_and_removed_legacy_options(args):
    with pytest.raises(SystemExit) as error:
        b.main(args)
    assert error.value.code == 2


def test_invalid_options_are_rejected_before_loading_an_external_curve():
    with patch.object(b, 'load_curve') as loader:
        with pytest.raises(ValueError, match='positive'):
            b.run_curve(ROW['label'], digits=0)
    loader.assert_not_called()


def test_genus_two_requirement_is_checked_before_experiment():
    x = PolynomialRing(QQ, 'x').gen()
    with pytest.raises(ValueError, match='genus-two'):
        b.run_curve(HyperellipticCurveWithReduction(x**7-x+1), with_lvalue=False)


def test_survey_bound_is_strict_and_needs_no_component_records(tmp_path):
    cache = tmp_path/'curves.json'
    rows = [dict(ROW, label=str(n), cond=n) for n in (999, 1000, 1001)]
    cache.write_text(json.dumps(rows))
    assert [row['cond'] for row in b.survey_curves(1000, cache=cache)] == [999]


def test_survey_records_timeouts_and_continues(tmp_path):
    cache, out = tmp_path/'curves.json', tmp_path/'experiment.json'
    rows = [dict(ROW, label=f'{n}.a.1.1', cond=n) for n in (394, 395)]
    cache.write_text(json.dumps(rows))

    result = b.run(400, cache=cache, out=out, with_lvalue=False, curve_seconds=0.01)
    assert [r['label'] for r in result] == [r['label'] for r in rows]
    assert [r['status'] for r in result] == ['curve_timeout', 'curve_timeout']
    assert json.loads(out.read_text()) == result
