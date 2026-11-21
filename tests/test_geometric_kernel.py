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
    kernel, info = worker.geometric_kernel(None, {}, seconds=3, certified=False,
                                          progress=messages.append)
    assert kernel == fallback
    assert calls == [(False, 3), (True, 3)]
    assert info['geometric_points'] == 22 and not info['packet_support_complete']
    assert 'exceeded 3s' in info['geometric_issue']
    assert any('rational-abscissa support' in message for message in messages)


def test_poonen_timeout_is_reported_without_an_empty_kernel(monkeypatch):
    def timeout(*args, **kwargs):
        raise TimeoutExpired('gp', kwargs['timeout'])

    monkeypatch.setattr(worker, 'geometric_torsion_packet', timeout)
    kernel, info = worker.geometric_kernel(None, {}, seconds=3, certified=False,
                                          progress=lambda message: None)
    assert kernel is None and info['status'] == 'packet_timeout'


def test_geometric_kernel_worker_returns_exact_point_encodings():
    row = dict(label='test', eqn=[[0, -1, 0, 0, 0, 1], []])
    data = [6, [[[0, 1], 1], [[1, 1], 1], [[-1, 1], 1], [[1, 0, 1], 1], [[], 1]]]
    packet = GeometricTorsionPacket(curve_from_row(row), data, metadata={'data': data})
    kernel = worker._kernel_in_worker(row, packet, True, 30, lambda message: None)
    assert len(kernel['points']) == 4
    assert all(P['field_polynomial'] is None for P in kernel['points'])
    assert ['infinity', '0'] in [P['coordinates'] for P in kernel['points']]
