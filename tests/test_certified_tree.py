"""Exact tree topology and independent checks of validated Stokes integrals."""
from collections import Counter
import json
from pathlib import Path
from unittest.mock import patch
from types import SimpleNamespace

import pytest
from sage.all import QQ, QQbar

from hyperell_regulator.regulator import regulators_stokes, IntegrationMonitor
from hyperell_regulator.regulator.balls import CertificationError, fields, narrow
from hyperell_regulator.regulator.polynomials import ring
from hyperell_regulator.regulator.tree_balls import TreeCut, _cyclic_order, _validate_embedding, _ray_tail

x = ring().gen()


@pytest.mark.parametrize('F,centers', [
    (x**5-x, [-1, 0]),
    (x**5+(x+1)**2, [-1, 0]),
    (x**6+4*x**5+4*x**4+2*x**3+1, [-1, 0]),
    (x*(x-1)*(x+1)*(x-2)*(x+2), [QQ(1)/2]),
])
def test_turns_and_sheet_gluing_match_exact_tree_components(F, centers):
    cut = TreeCut(F, centers, 12)
    assert len(cut.edges) == len(cut.nodes)
    assert Counter(cut.boundary_word) == Counter((i,s) for i in range(len(cut.edges)) for s in (-1,1))
    assert all(b == c for (_,b),(c,_) in zip(cut.banks,cut.banks[1:]+cut.banks[:1]))
    assert cut.edges[-1][1] is None
    parents = dict(cut.edges)
    for i, (a, _) in enumerate(cut.edges):
        descendants = []
        for vertex in range(len(cut.nodes)):
            ancestor = vertex
            while ancestor is not None and ancestor != a:
                ancestor = parents[ancestor]
            if ancestor == a:
                descendants.append(vertex)
        differences = [u-v for u,v in zip(cut.turns[i,1],cut.turns[i,-1])]
        assert differences == [-1 if vertex in descendants else 0 for vertex in range(len(cut.nodes))]
        assert sum(differences[:len(cut.roots)]) % 2 == sum(v < len(cut.roots) for v in descendants) % 2


def test_embedding_rejects_crossings_and_interior_vertices():
    C, _ = fields(10)
    ii = QQbar.gen()
    nodes = [-1-ii,1+ii,-1+ii,1-ii]
    with pytest.raises(CertificationError, match='disjoint'):
        _validate_embedding(nodes, [(0,1),(2,3)], C)
    with pytest.raises(CertificationError, match='another cut vertex'):
        _validate_embedding(list(map(QQbar,[0,1,2])),[(0,2)], C)
    with pytest.raises(CertificationError, match='overlap'):
        _cyclic_order({0:QQbar(1),1:QQbar(2)})


def test_clockwise_order_resolves_directions_collapsed_in_float64():
    eps = QQ(1)/10**100
    directions = {0:QQbar(1), 1:QQbar(1)+QQbar.gen()*eps,
                  2:QQbar(-1), 3:-QQbar.gen()}
    assert _cyclic_order(directions) == [0,1,2,3]


def test_tree_proposal_handles_distinct_vertices_collapsed_in_float64():
    from hyperell_regulator.regulator import tree_balls
    with patch.object(tree_balls,'_precision_tree',wraps=tree_balls._precision_tree) as fallback:
        cut = TreeCut(x**5+x+1,[QQ(1),QQ(1)+QQ(1)/10**70],100)
    fallback.assert_called_once()
    assert len(cut.nodes) == 7 and len(cut.edges) == 7


@pytest.mark.parametrize('reference', json.loads(
    (Path(__file__).parent/'data/semilunar_reference.json').read_text())['cases'])
def test_default_tree_never_uses_sampled_topology(reference):
    F = ring()(reference['F'])
    details = {}
    with patch('hyperell_regulator.regulator.stokes.Cut.topology', side_effect=AssertionError('sampled topology')):
        with IntegrationMonitor(max_edge_seconds=60) as mon:
            per, cols = regulators_stokes(F,reference['centers'],pole=reference['pole'],
                                          prec=10,certified=True,details=details)
    assert details['certified'] and details['cut_kind'] == 'tree'
    assert len(mon.records) == len(details['cut'].edges)
    assert sum('ray' in record['edge'] for record in mon.records) == 1
    # Frozen values were computed by the independent semilunar backend.
    expected = [reference['periods']]+reference['regulators']
    for actual, previous in zip([per]+cols, expected):
        for a,b in zip(actual,previous):
            assert narrow(a,10) and abs(a-a.parent()(b)) < QQ(1)/10**18


def test_uncertain_tree_topology_fails_without_parallel_fallback():
    with patch('hyperell_regulator.regulator.tree_balls._validate_embedding',
               side_effect=CertificationError('uncertain embedding')):
        with pytest.raises(CertificationError,match='uncertain embedding'):
            regulators_stokes(x**5+x+1,[0],prec=10,certified=True)
    with pytest.raises(ValueError,match='certified_cut'):
        regulators_stokes(x**5+x+1,[0],certified=True,certified_cut='unknown')


@pytest.mark.parametrize('degree', [5,6])
def test_ray_endpoint_bounds_against_analytic_beta_integrals_at_100_digits(degree):
    # Repeated roots outside the ray are intentional: this isolates the ray
    # kernel against elementary beta integrals, not a smooth curve regulator.
    cut = SimpleNamespace(F=(x+1)**degree, roots=[QQbar(-1)]*degree,
                          nodes=[QQbar(0)], ray_vertex=0, span=QQ(1),
                          references=[QQbar(1)], g=2, centers=[QQ(-1)])
    with IntegrationMonitor(max_edge_steps=100000,max_edge_seconds=60) as mon:
        with mon.edge('analytic ray'):
            result = _ray_tail(cut,100)
    C, _ = fields(120)
    if degree == 5:
        root2 = C(2).sqrt()
        log0 = root2*(C(2).log()/6+C(1)/9)
        expected = [root2/6,5*root2/6,log0,root2*(C(2).log()+2)-log0]
    else:
        expected = [C(1)/8,C(3)/8,C(2).log()/8+C(1)/16,
                    3*C(2).log()/8+C(7)/16]
    for actual, value in zip(result['periods']+sum(result['logarithmic'],[]),expected):
        assert actual.overlaps(value) and narrow(actual,100)
