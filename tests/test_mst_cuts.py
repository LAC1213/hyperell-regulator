"""Minimum length, tree boundary order, gluing parity and numerical agreement."""
from collections import Counter
import json
from pathlib import Path
import numpy as np
import pytest
from scipy.sparse.csgraph import minimum_spanning_tree as scipy_mst
from hyperell_regulator.regulator import regulators_stokes
from hyperell_regulator.regulator._tree import minimum_spanning_tree, orient_to_root, boundary_word
from hyperell_regulator.regulator.polynomials import ring
from hyperell_regulator.regulator.stokes import Cut


def test_mst_agrees_with_independent_graph_algorithm():
    points=np.array([0, 2+.2j, -.7+1j, -.2-1.1j, 3+2j, 4-.4j])
    edges=minimum_spanning_tree(points)
    distance=abs(points[:,None]-points[None,:])
    assert len(edges)==len(points)-1
    assert abs(sum(distance[a,b] for a,b in edges)-scipy_mst(distance).sum())<1e-12
    assert edges==minimum_spanning_tree(points)


def test_four_way_junction_visits_each_bank_once():
    points=[0j,1+0j,1j,-1+0j,-1j]
    edges=orient_to_root(5,minimum_spanning_tree(points),1)
    _,_,word,banks=boundary_word(points,edges,1)
    assert len(word)==2*len(edges)
    assert Counter(word)==Counter((i,s) for i in range(len(edges)) for s in (1,-1))
    assert all(a[1]==b[0] for a,b in zip(banks,banks[1:]+banks[:1]))
    assert sum(abs(points[a]-points[b]) for a,b in edges[:-1])==4


def subtree(cut, edge):
    children={a:[] for a in range(len(cut.nodes))}
    for a,b in cut.edges[:-1]:children[b].append(a)
    descendants=[cut.edges[edge][0]]
    for a in descendants:descendants.extend(children[a])
    return descendants


@pytest.mark.parametrize('expression', ['x^5-x', 'x^6+4*x^5+4*x^4+2*x^3+1'])
def test_log_jumps_and_sheet_gluing_follow_tree_subtrees(expression):
    cut=Cut(ring()(expression),[-1,0])
    assert len(cut.edges)==len(cut.nodes)
    assert cut.edges[-1][1] is None and cut.direction==1
    for i in range(len(cut.edges)):
        descendants=subtree(cut,i)
        branch_count=sum(cut.labels[v] is not None for v in descendants)
        permutation=[cut.lower[i].index(cut.upper[i][s]) for s in (0,1)]
        assert permutation==([1,0] if branch_count%2 else [0,1])
        for k,a in enumerate(cut.centers):
            expected=-2*cut.ar.R.pi() if any(cut.nodes[v]==a for v in descendants) else 0
            assert abs(cut.jumps[k][i].imag()-expected)<cut.ar.R('1e-30')


def test_mst_integrals_match_frozen_independent_backend():
    reference=json.loads((Path(__file__).parent/'data/semilunar_reference.json').read_text())
    for row in reference['cases']:
        per,cols=regulators_stokes(ring()(row['F']),row['centers'],pole=row['pole'],prec=25)
        R=per[0].parent()
        for actual,expected in zip([per]+cols,[row['periods']]+row['regulators']):
            assert max(abs(a-R(b)) for a,b in zip(actual,expected))<R('1e-23')


def test_real_vertices_need_no_detour_or_extra_edges():
    x=ring().gen();F=x*(x-1)*(x+1)*(x-2)*(x+2)
    cut=Cut(F,[0])
    assert len(cut.nodes)==5 and len(cut.edges)==5
    assert all(z.imag()==0 for z in cut.nodes)
    assert all(abs(cut.nodes[a]-cut.nodes[b])==1 for a,b in cut.edges[:-1])


def test_removed_backend_selector_is_not_accepted():
    with pytest.raises(TypeError):
        regulators_stokes(ring().gen()**5+1,[0],cut='semilunar')
