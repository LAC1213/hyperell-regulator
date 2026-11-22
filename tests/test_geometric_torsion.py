"""Geometric packets, exact algebraic orders, and failure-safe GP wrapping."""

import hashlib
import json
from pathlib import Path
import shutil
import subprocess
from unittest.mock import patch

import pytest
from sage.all import PolynomialRing, QQ, QQbar

from hyperell_regulator.components import curve_from_row
from hyperell_regulator.geometric_torsion import (
    GeometricTorsionPacket, _validate_result, geometric_torsion_packet,
    install_poonen_backend, RESULT_MARKER,
)


DATA_249 = [10, [[[-3, 1], 1], [[1, 1], 1], [[-3, -2, -1, -2, 1], 1],
                 [[], 7], [[0, 1], 14]]]
DATA_256 = [14, [[[-1, 2], 1], [[1, -2, 2], 1], [[-1, 0, 2], 1],
                 [[], 1], [[-1, 1], 5], [[0, 1], 5], [[1, -4, 2], 5]]]


def curve_249():
    return curve_from_row(dict(eqn=[[2, 3, 1, 1, 0, -1], [1, 0, 0, 1]]))


def curve_256():
    return curve_from_row(dict(eqn=[[0, -1, 1, 1, -3, 2], [1]]))


def test_249_complete_packet_has_four_new_order_two_points():
    T = GeometricTorsionPacket(curve_249(), DATA_249, basepoint=(-1, 0))
    assert not T.complete  # bare polynomial data are not a completeness certificate
    points = T.points()
    assert T.size == len(points) == 10
    new = [P for P in points if P[0] != 'infinity' and P[0] not in QQ]
    assert len(new) == 4
    x = PolynomialRing(QQ, 'x').gen()
    g = x**4-2*x**3-x*x-2*x-3
    for P in new:
        assert g(P[0]) == 0
        assert P[1] == -(1+P[0]**3)/2
        assert T.order(P) == 2
    assert T.order((0, 1)) == 28
    assert T.order(('infinity', 0)) == 7
    assert T.order((-1, 0)) == 1


def test_256_packet_and_change_to_nonweierstrass_algebraic_base():
    T = GeometricTorsionPacket(curve_256(), DATA_256, basepoint=(0, -1))
    points = T.points()
    assert T.size == len(points) == 14
    new = [P for P in points if P[0] != 'infinity' and P[0] not in QQ]
    assert len(new) == 8
    assert all(T.order(P) == 10 for P in new)
    T2 = GeometricTorsionPacket(curve_256(), DATA_256, basepoint=new[0])
    assert set(T2.points()) == set(points)
    assert T2.order(new[0]) == 1
    assert T2.order((0, -1)) == 10


def test_nonhyperelliptic_base_is_not_silently_given_the_wrong_packet():
    C = curve_from_row(dict(eqn=[[1, -1, 0, 0, 0, 1], []]))
    with pytest.raises(NotImplementedError, match='outside the hyperelliptic torsion packet'):
        GeometricTorsionPacket(C, [6, [[[1, -1, 0, 0, 0, 1], 1], [[], 1]]], basepoint=(0, 1))
    T = GeometricTorsionPacket(curve_249(), DATA_249)
    with pytest.raises(ValueError, match='not on'):
        T.contains((0, 0))
    with pytest.raises(ValueError, match='which point'):
        T.contains('infinity')


@pytest.mark.parametrize('mutation', ['count', 'missing', 'duplicate', 'order'])
def test_malformed_backend_result_is_rejected(mutation):
    data = json.loads(json.dumps(DATA_249))
    if mutation == 'count':
        data[0] += 1
    elif mutation == 'missing':
        data[1].pop(0)
    elif mutation == 'duplicate':
        data[1].append(data[1][0])
    else:
        data[1][0][1] = 0
    with pytest.raises(ValueError):
        _validate_result(curve_249().F.change_ring(QQ), data)


def fake_backend(tmp_path):
    path = tmp_path/'torsion.gp'
    path.write_text('fixture backend')
    return path, hashlib.sha256(path.read_bytes()).hexdigest()


def test_backend_requires_success_and_checks_cache_and_source(tmp_path, monkeypatch):
    monkeypatch.setenv('HYPERELL_REGULATOR_DATA', str(tmp_path))
    backend, digest = fake_backend(tmp_path)
    result = subprocess.CompletedProcess([], 0,
        'SUCCESS!\n'+RESULT_MARKER+json.dumps(DATA_249)+'\n', '')
    with patch('hyperell_regulator.geometric_torsion.SOURCE_SHA256', digest), \
         patch('hyperell_regulator.geometric_torsion._subprocess.run', return_value=result) as run:
        first = geometric_torsion_packet(curve_249(), backend=backend, gp='gp')
        second = geometric_torsion_packet(curve_249(), backend=backend, gp='gp')
        assert first.complete and second.metadata['cached']
        assert first.points() == second.points() and run.call_count == 1
        result.stdout = RESULT_MARKER+json.dumps(DATA_249)+'\n'
        with pytest.raises(ArithmeticError, match='did not certify'):
            geometric_torsion_packet(curve_249(), backend=backend, gp='gp', cache=False)
    with pytest.raises(ValueError, match='pinned'):
        geometric_torsion_packet(curve_249(), backend=backend, gp='gp')
    with pytest.raises(ValueError, match='checksum'):
        install_poonen_backend(tmp_path/'destination.gp', source_path=backend)
    assert not (tmp_path/'destination.gp').exists()


def test_timeout_propagates_and_does_not_create_a_success_cache(tmp_path, monkeypatch):
    monkeypatch.setenv('HYPERELL_REGULATOR_DATA', str(tmp_path))
    backend, digest = fake_backend(tmp_path)
    with patch('hyperell_regulator.geometric_torsion.SOURCE_SHA256', digest), \
         patch('hyperell_regulator.geometric_torsion._subprocess.run',
               side_effect=subprocess.TimeoutExpired('gp', 1)):
        with pytest.raises(subprocess.TimeoutExpired):
            geometric_torsion_packet(curve_249(), backend=backend, gp='gp', timeout=1)
    assert not list(tmp_path.rglob('*.json'))


def test_live_upstream_algorithm_when_backend_is_installed():
    backend = Path('results/poonen/torsion.gp')
    if not backend.is_file() or not shutil.which('gp'):
        pytest.skip('optional upstream GP backend is not installed')
    # A separate curve from the saved LMFDB examples exercises the complete
    # source, mod-p^2 sieve, lifting, and exact verification in a fresh GP.
    C = curve_from_row(dict(eqn=[[1, 0, 0, 0, 0, 0, 1], []]))
    T = geometric_torsion_packet(C, backend=backend, cache=False, timeout=120)
    assert T.complete and T.size == 10 and len(T.points()) == 10
