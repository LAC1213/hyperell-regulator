"""Budgeted geometric searches retain exact information about support."""

from subprocess import TimeoutExpired
from types import SimpleNamespace

from hyperell_regulator.components import curve_from_row
from hyperell_regulator.experiments import _geometric_kernel as worker
from hyperell_regulator.geometric_torsion import GeometricTorsionPacket


def test_geometric_kernel_timeout_retries_on_rational_abscissae(monkeypatch):
    packet = SimpleNamespace(size=22, metadata={'cached': True},
                             factors=[(SimpleNamespace(degree=lambda: 6), 1)])
    monkeypatch.setattr(worker, 'geometric_torsion_packet', lambda *args, **kwargs: packet)
    calls = []
    fallback = dict(points=[1, 2, 3], basis=[], complete=True, places={})

    def search(row, packet, rational_only, seconds, progress):
        calls.append((rational_only, seconds))
        if not rational_only:
            raise TimeoutExpired('sage', seconds)
        return fallback

    monkeypatch.setattr(worker, '_kernel_in_worker', search)
    messages = []
    kernel, info = worker.geometric_kernel(None, {}, seconds=3, progress=messages.append)
    assert kernel == fallback
    assert calls == [(False, 3), (True, 3)]
    assert info['geometric_points'] == 22 and not info['packet_support_complete']
    assert 'exceeded 3s' in info['geometric_issue']
    assert any('rational-abscissa support' in message for message in messages)


def test_poonen_timeout_is_reported_without_an_empty_kernel(monkeypatch):
    def timeout(*args, **kwargs):
        raise TimeoutExpired('gp', kwargs['timeout'])

    monkeypatch.setattr(worker, 'geometric_torsion_packet', timeout)
    kernel, info = worker.geometric_kernel(None, {}, seconds=3, progress=lambda message: None)
    assert kernel is None and info['status'] == 'packet_timeout'


def test_geometric_kernel_worker_returns_and_caches_exact_point_encodings(tmp_path, monkeypatch):
    monkeypatch.setenv('HYPERELL_REGULATOR_DATA', str(tmp_path))
    row = dict(label='test', eqn=[[0, -1, 0, 0, 0, 1], []])
    data = [6, [[[0, 1], 1], [[1, 1], 1], [[-1, 1], 1], [[1, 0, 1], 1], [[], 1]]]
    packet = GeometricTorsionPacket(curve_from_row(row), data, metadata={'data': data})
    kernel = worker._kernel_in_worker(row, packet, True, 30, lambda message: None)
    assert len(kernel['points']) == 4
    assert all(P['field_polynomial'] is None for P in kernel['points'])
    assert ['infinity', '0'] in [P['coordinates'] for P in kernel['points']]
    assert list((tmp_path/'results'/'integral_bases').glob('*.json'))

    def unused(*args, **kwargs):
        raise AssertionError('the cached kernel should bypass the Sage worker')

    monkeypatch.setattr(worker._subprocess, 'run', unused)
    assert worker._kernel_in_worker(row, packet, True, 1, lambda message: None) == kernel


def test_certified_survey_uses_only_integral_traces_on_277(tmp_path, monkeypatch):
    import json

    from hyperell_regulator.experiments import beilinson

    monkeypatch.setenv('HYPERELL_REGULATOR_DATA', str(tmp_path))
    row = dict(label='277.a.277.2', cond=277, end_alg='Q', mw_invs=[5],
               eqn=[[-6, 11, -19, 14, -9, 1], [1]], rat_pts=[[1, 0, 0]])
    data = [10, [[[-23, 44, -76, 56, -36, 4], 1], [[], 1], [[1, -1, 1], 15]]]
    packet = GeometricTorsionPacket(curve_from_row(row), data, metadata={'complete': True, 'data': data})
    monkeypatch.setattr(worker, 'geometric_torsion_packet', lambda *args, **kwargs: packet)
    cache = tmp_path/'curves.json'
    cache.write_text(json.dumps([row]))
    result = beilinson.analyse(row, cache=cache, packet='geometric', certified=True, digits=20,
                               with_lvalue=False, with_periods=False, verbose=False)
    assert result['packet_support_complete'] and result['integral_rank_complete']
    assert result['integral_rank'] == 1 and result['geometric_points'] == 10
    kernel, = result['integral_kernels']
    assert kernel['places']['277']['method'] == 'simple_node_galois_pairs'
    assert all(not any(map(int, D[:5])) for D in kernel['basis'])
    assert result['status'] == 'insufficient_integral_rank'

    monkeypatch.setattr(beilinson, 'standard_lvalue_cached', lambda *args, **kwargs: dict(
        L_second_derivative='6.547289332454137727784192', L_digits=20))
    result = beilinson.analyse(row, cache=cache, packet='geometric', digits=20, reg_digits=20,
                               ignore_integrality=True, verbose=False)
    assert result['status'] == 'ok' and result['packet_used'] == 'geometric'
    assert result['integral_rank'] == 1
    kernel, = result['integral_kernels']
    assert len(kernel['divisor_basis']) == 2
    for D in kernel['divisor_basis']:
        assert len(set(D[:5])) == 1
        assert len(set(D[6:])) == 1
    assert result['pair_searches'] and not result['linear_comparison'].get('terms')


def test_geometric_survey_finds_full_integral_rank_on_353(tmp_path, monkeypatch):
    import json

    from hyperell_regulator.experiments import beilinson

    monkeypatch.setenv('HYPERELL_REGULATOR_DATA', str(tmp_path))
    row = dict(label='353.a.353.1', cond=353, end_alg='Q', mw_invs=[11],
               eqn=[[0, 0, 1], [1, 1, 0, 1]],
               rat_pts=[[0, -1, 1], [0, 0, 1], [1, -1, 0], [1, 0, 0]])
    data = [10, [[[1, 2, 5, 2, 2, 0, 1], 1], [[], 11], [[0, 1], 11]]]
    packet = GeometricTorsionPacket(curve_from_row(row), data, metadata={'complete': True, 'data': data})
    monkeypatch.setattr(worker, 'geometric_torsion_packet', lambda *args, **kwargs: packet)
    cache = tmp_path/'curves.json'
    cache.write_text(json.dumps([row]))
    result = beilinson.analyse(row, cache=cache, packet='geometric', with_lvalue=False,
                               with_periods=False, verbose=False)
    assert result['packet_support_complete'] and result['integral_rank_complete']
    assert result['status'] == 'insufficient_integral_rank'
    assert result['integral_rank'] == 1 and result['geometric_points'] == 10
    kernel = result['integral_kernels'][-1]
    assert kernel['places']['353']['collision_block_size'] == 3
    divisor, = kernel['basis']
    assert divisor[:6] == ['0']*6
    assert divisor[6:] in (['1', '1', '-1', '-1'], ['-1', '-1', '1', '1'])


def test_geometric_survey_handles_mixed_branch_collision_on_427(tmp_path, monkeypatch):
    import json

    from hyperell_regulator.experiments import beilinson

    monkeypatch.setenv('HYPERELL_REGULATOR_DATA', str(tmp_path))
    row = dict(label='427.a.2989.1', cond=427, end_alg='Q', mw_invs=[14],
               eqn=[[-4, 4, 4, -5, -1, 1], [1, 0, 0, 1]],
               rat_pts=[[-3, 13, 1], [1, -1, 0], [1, -1, 1], [1, 0, 0]])
    data = [8, [[[-1, 1], 1], [[3, 1], 1], [[5, -2, -5, 2, 1], 1], [[], 7]]]
    packet = GeometricTorsionPacket(curve_from_row(row), data, metadata={'complete': True, 'data': data})
    monkeypatch.setattr(worker, 'geometric_torsion_packet', lambda *args, **kwargs: packet)
    cache = tmp_path/'curves.json'
    cache.write_text(json.dumps([row]))
    result = beilinson.analyse(row, cache=cache, packet='geometric', with_lvalue=False,
                               with_periods=False, verbose=False)
    assert result['packet_support_complete'] and result['integral_rank_complete']
    assert result['status'] == 'insufficient_integral_rank'
    assert result['integral_rank'] == 1 and result['geometric_points'] == 8
    kernel = result['integral_kernels'][-1]
    assert kernel['places']['7']['branch_factors'] == ['x + 3', 'x^4 + 2*x^3 - 5*x^2 - 2*x + 5']
    assert kernel['places']['7']['node_thickness'] == 2
    assert kernel['places']['61']['collision_block_size'] == 4
    divisor, = kernel['basis']
    assert divisor[1:6] == ['0']*5
    assert [int(divisor[i]) for i in (0, 6, 7)] in ([2, -1, -1], [-2, 1, 1])


def test_geometric_survey_proves_zero_rank_on_604_at_one_prime(tmp_path, monkeypatch):
    import json

    from hyperell_regulator.experiments import beilinson

    monkeypatch.setenv('HYPERELL_REGULATOR_DATA', str(tmp_path))
    row = dict(label='604.a.9664.1', cond=604, end_alg='Q', mw_invs=[],
               eqn=[[-21, -53, -4, 48, 9, 4], [1, 1, 1]], rat_pts=[[1, 0, 0]])
    data = [6, [[[-83, -210, -13, 194, 37, 16], 1], [[], 1]]]
    packet = GeometricTorsionPacket(curve_from_row(row), data, metadata={'complete': True, 'data': data})
    monkeypatch.setattr(worker, 'geometric_torsion_packet', lambda *args, **kwargs: packet)
    cache = tmp_path/'curves.json'
    cache.write_text(json.dumps([row]))
    result = beilinson.analyse(row, cache=cache, packet='geometric', with_lvalue=False,
                               with_periods=False, verbose=False)
    assert result['packet_support_complete'] and result['integral_rank_complete']
    assert result['status'] == 'insufficient_integral_rank'
    assert result['integral_rank'] == 0 and result['geometric_points'] == 6
    kernel, = result['integral_kernels']
    assert kernel['places']['151']['collision_block_size'] == 5
    assert kernel['places']['151']['exact']
    assert kernel['places']['2']['method'] == 'skipped_after_zero_kernel'
    assert not kernel['basis']
    assert kernel['divisor_integrality'] and all(v is False for v in kernel['divisor_integrality'])


def test_geometric_kernel_restricts_dyadic_arithmetic_to_surviving_support(tmp_path, monkeypatch):
    monkeypatch.setenv('HYPERELL_REGULATOR_DATA', str(tmp_path))
    row = dict(label='708.a.19116.1', cond=708,
               eqn=[[-1, 4, 4, 0, 0, -1], [1, 0, 0, 1]])
    data = [8, [[[-3, 1], 1], [[1, 1], 1], [[1, -6, -1, -2, 1], 1], [[], 5]]]
    C = curve_from_row(row)
    packet = GeometricTorsionPacket(C, data, metadata={'complete': True, 'data': data})
    monkeypatch.setattr(worker, 'geometric_torsion_packet', lambda *args, **kwargs: packet)

    kernel, info = worker.geometric_kernel(C, row, seconds=15, progress=lambda message: None)
    assert info['packet_support_complete'] and kernel['complete']
    assert len(kernel['points']) == 8
    assert len(kernel['basis']) == kernel['rank_upper_bound'] == 1
    assert kernel['places']['Fractional ideal (2)']['support_indices'] == [0, 6, 7]
