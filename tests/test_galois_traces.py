"""Integral bases over QQ use full Galois traces, not individual embeddings."""

import json

from sage.all import PolynomialRing, QQ, QQbar, vector

from hyperell_regulator._galois import galois_trace_basis
from hyperell_regulator.components import curve_from_row
from hyperell_regulator.experiments import _geometric_kernel as worker, beilinson
from hyperell_regulator.experiments._kernel_cache import KERNEL_CACHE_VERSION
from hyperell_regulator.geometric_torsion import GeometricTorsionPacket
from hyperell_regulator.motivic import MotivicPacket


def test_trace_basis_uses_all_conjugates_including_nonreal_embeddings():
    x = PolynomialRing(QQ, 'x').gen()
    quadratic = (x*x-2).roots(QQbar, multiplicities=False)
    cubic = (x**3-2).roots(QQbar, multiplicities=False)
    xs = [*quadratic, *cubic, 'infinity']
    traces = galois_trace_basis(xs)
    assert traces.nrows() == 2
    for row in traces:
        assert sum(row) == 0
        assert row[0] == row[1]
        assert row[2] == row[3] == row[4]
    assert vector(QQ, [1, -1, 0, 0, 0, 0]) not in traces.row_space()
    assert vector(QQ, [1, 1, 0, 0, 0, -2]) in traces.row_space()
    assert vector(QQ, [0, 0, 1, 1, 1, -3]) in traces.row_space()


def test_incomplete_orbits_contribute_no_trace_direction():
    root = QQbar(2).sqrt()
    traces = galois_trace_basis([root, 0, 'infinity'])
    assert traces.rows() == [vector(QQ, [0, -1, 1])]
    assert galois_trace_basis([root, 'infinity']).dimensions() == (0, 2)
    assert galois_trace_basis([]).dimensions() == (0, 0)


def test_rational_support_keeps_the_full_degree_zero_space():
    traces = galois_trace_basis([0, 1, -1, 'infinity'])
    assert traces.nrows() == 3
    assert all(sum(row) == 0 for row in traces)


def test_394_trace_kernel_excludes_the_locally_integral_quadratic_difference(tmp_path, monkeypatch):
    monkeypatch.setenv('HYPERELL_REGULATOR_DATA', str(tmp_path))
    row = dict(label='394.a.394.1', cond=394, end_alg='Q', mw_invs=[10],
               eqn=[[-9, 17, 0, -12, 1, 2], [0, 1, 0, 1]], rat_pts_verified=True,
               rat_pts=[[1, -1, 0], [1, -1, 1], [1, 0, 0]])
    data = [8, [[[-1, 1], 1], [[-4, 4, 1], 1], [[-9, -1, 5, 1], 1], [[], 5]]]
    packet = GeometricTorsionPacket(curve_from_row(row), data,
                                    metadata={'complete': True, 'data': data})
    motivic = MotivicPacket.from_geometric_packet(packet)
    algebraic = motivic.divisor([0, 1, -1, 0, 0, 0, 0, 0])
    assert motivic.integrality.check(algebraic.coefficients)['q_integral'] is True
    basis, complete = motivic.integral_basis()
    assert complete and len(basis) == 1
    D, = basis
    assert D.coefficients[1:6] == (0, 0, 0, 0, 0)
    assert D.coefficients[0] == -2*D.coefficients[6] == -2*D.coefficients[7] != 0

    monkeypatch.setattr(worker, 'geometric_torsion_packet', lambda *args, **kwargs: packet)

    def no_regulator(*args, **kwargs):
        raise AssertionError('an anti-invariant divisor must not trigger a regulator pair')

    from hyperell_regulator.experiments import _beilinson_divisors
    monkeypatch.setattr(_beilinson_divisors, '_divisor_evaluator', no_regulator)
    cache = tmp_path/'curves.json'
    cache.write_text(json.dumps([row]))
    result = beilinson.analyse(row, cache=cache, packet='geometric', with_lvalue=False,
                               with_periods=False, verbose=False)
    assert result['status'] == 'insufficient_integral_rank'
    assert result['integral_rank'] == 1 and result['integral_rank_complete']
    kernel = result['integral_kernels'][-1]
    assert kernel['divisor_descent'] == 'galois_trace' and kernel['base_field'] == 'Q'
    assert len(kernel['basis']) == 1
    assert len(kernel['divisor_basis']) == 3
    for D in kernel['divisor_basis']:
        assert D[1] == D[2]
        assert D[3] == D[4] == D[5]


def test_saved_number_field_basis_is_not_reused_as_a_trace_basis():
    settings = dict(packet='geometric', regulator_digits=20, lvalue_digits=20,
                    with_lvalue=False, tries=1, integral_basis_cache_version=KERNEL_CACHE_VERSION)
    old = dict(status='ok', settings=dict(settings, integral_basis_cache_version=2))
    assert not beilinson._stored_result_satisfies(old, settings)
    old['settings'] = dict(settings)
    assert beilinson._stored_result_satisfies(old, settings)
