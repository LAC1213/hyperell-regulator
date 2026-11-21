"""Regulator and graph calculations from a common torsion-packet divisor."""

from sage.all import PolynomialRing, QQ, RealField

from hyperell_regulator import HyperellipticCurveWithReduction, MotivicPacket, stable_model
from hyperell_regulator.components import curve_from_row
from hyperell_regulator.experiments import beilinson


def test_weierstrass_divisor_and_tame_stable_model():
    x = PolynomialRing(QQ, 'x').gen()
    C = HyperellipticCurveWithReduction(x*(x-1)*(x-2)*(x-3)*(x-4))
    packet = MotivicPacket(C)
    assert len(packet.points) == 6
    left = packet.divisor({(0, 0): 1, (1, 0): -1})
    right = packet.divisor({(2, 0): 1, (3, 0): -1})
    det = left.determinant(right, digits=20)
    assert abs(det) > QQ(1)/10
    twice = packet.divisor({(0, 0): 2, (1, 0): -2})
    assert abs(twice.determinant(right, digits=20) - 2*det) < QQ(10)**-19
    model = stable_model(C, 3, points=packet.points)
    assert model.prime.ramification_index() == 2
    assert model.stable_graph.genera == (0,)
    assert model.stable_graph.cycles.nrows() == 2
    assert len(model.node_coordinates) == len(model.graph.edges)
    assert not model.boundary([1, 0, 0, -1, 0, 0]).projected
    assert model.boundary(left.coefficients).projected


def test_dyadic_smooth_chart_over_finite_extension():
    C = curve_from_row(dict(label='688.a.2752.1', eqn=[[0,-1,0,4,-5,2], [1]]))
    model = stable_model(C, 2, points=[(0, 0), ('infinity', 0)])
    assert model.prime.ramification_index() == 5
    assert model.stable_graph.genera == (2,)
    assert not model.stable_graph.edges
    assert not model.boundary([1, -1]).projected


def test_dyadic_compact_type_graph_resolves_the_integral_kernel():
    C = curve_from_row(dict(label='464.a.29696.2', eqn=[[0,1,16,72,33,4], [0,1]]))
    packet = MotivicPacket(C, [(-4, 2), (0, 0), ('infinity', 0)], torsion_multiples=2)
    model = packet.stable_model(2)
    assert model.stable_graph.genera == (1, 1)
    assert model.stable_graph.edges == ((0, 1, QQ(5)/6),)
    assert not packet.divisor([1, -1, 0]).boundary(2).projected
    assert packet.divisor([0, -1, 1]).boundary(2).projected
    basis, complete = packet.integral_basis()
    assert complete and len(basis) == 1


def test_geometric_experiment_uses_integer_divisors_and_skips_unneeded_search(tmp_path, monkeypatch):
    import json
    from hyperell_regulator.experiments import _geometric_kernel

    row = dict(label='249.a.6723.1', cond=249, end_alg='Q', mw_invs=[28],
               eqn=[[2,3,1,1,0,-1], [1,0,0,1]], rat_pts_verified=True,
               rat_pts=[[-1,0,1],[0,-2,1],[0,1,1],[1,-1,0],[1,0,0],[3,-14,1]])
    cache = tmp_path/'curves.json'
    cache.write_text(json.dumps([row]))

    def unused(*args, **kwargs):
        raise AssertionError('Poonen should be skipped when the rational kernel has rank two')

    monkeypatch.setattr(_geometric_kernel, 'geometric_kernel', unused)
    result = beilinson.analyse(row, cache=cache, packet='geometric', digits=20, reg_digits=20,
                               with_periods=False, with_lvalue=False, verbose=False)
    assert result['status'] == 'ok'
    assert result['packet_selection'] == 'geometric' and result['packet_used'] == 'rational'
    assert result['integral'] and result['basis_complete']
    basis = result['integral_kernels'][0]['basis']
    assert len(basis) == 2
    assert all(QQ(a).is_integer() for d in basis for a in d)
    assert list(result['divisors']) == basis
    assert any(sum(bool(QQ(a)) for a in d) > 2 for d in result['divisors'])
    assert abs(abs(RealField(100)(result['det']))-RealField(100)('66.9754630253291197696')) < QQ(10)**-18


def test_geometric_experiment_keeps_rank_one_rational_kernel_before_fallback(tmp_path, monkeypatch):
    import json
    from hyperell_regulator.experiments import _geometric_kernel

    row = dict(label='464.a.29696.2', cond=464, end_alg='Q', mw_invs=[2, 2],
               eqn=[[0,1,16,72,33,4], [0,1]], rat_pts_verified=True,
               rat_pts=[[-4,2,1],[0,0,1],[1,0,0]])
    cache = tmp_path/'curves.json'
    cache.write_text(json.dumps([row]))
    monkeypatch.setattr(_geometric_kernel, 'geometric_kernel',
                        lambda *args, **kwargs: (None, {'status': 'packet_timeout'}))
    result = beilinson.analyse(row, cache=cache, packet='geometric', with_periods=False,
                               with_lvalue=False, verbose=False)
    assert result['packet_used'] == 'geometric' and result['status'] == 'packet_timeout'
    assert len(result['integral_kernels'][0]['basis']) == 1
