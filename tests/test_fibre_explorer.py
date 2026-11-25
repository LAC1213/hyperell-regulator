"""Regular fibre graphs and the exact local charge data shown by the viewer."""

import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from sage.all import QQ, matrix, vector

from hyperell_regulator.experiments.fibre_explorer import compute_prime, regular_graph
from hyperell_regulator.graph_boundary import FibreGraph
from hyperell_regulator.stable_model import StableModel


ROOT = Path(__file__).resolve().parents[1]


def test_self_node_becomes_a_strict_normal_crossing_cycle():
    source = StableModel(None, None, FibreGraph([1], [(0, 0, 1)]), ((0,),))
    graph, scale = regular_graph(source)
    assert scale == 2
    assert graph.genera == (1, 0)
    assert graph.edges == ((0, 1, 1), (1, 0, 1))
    assert graph.cycles.nrows() == 1
    assert graph.boundary([-1, 1]).potential == vector(QQ, [0, QQ(1)/2])


def test_metric_resolution_and_wild_ramification_preserve_marked_vertices():
    source = StableModel(None, None, FibreGraph([1, 1], [(0, 1, QQ(3)/2)]), ((0,), (1,)))
    graph, scale = regular_graph(source)
    assert scale == 2 and graph.genera == (1, 1, 0, 0)
    assert graph.boundary([-2, 2, 0, 0]).potential[1] == 6
    prime = SimpleNamespace(ramification_index=lambda: 2)
    chart_prime = SimpleNamespace(ramification_index=lambda: 8)
    source = StableModel(None, prime, FibreGraph([1, 1], [(0, 1, QQ(3)/8)]), (),
                         charts=({'centre_prime': chart_prime, 'rho_power': (3, 8)},),
                         length_normalization='v(2)=1')
    graph, scale = regular_graph(source)
    assert scale == 96  # Include each extension once, then clear edge denominators.
    assert len(graph.edges) == 36 and all(a != b and n == 1 for a, b, n in graph.edges)


def test_paper_example_distinguishes_the_two_tangent_branches():
    catalog = json.loads((ROOT/'docs/_extra/explorer/catalog.json').read_text())
    entry = next(c for c in catalog['curves'] if c['label'] == '249.a.6723.1')
    result = compute_prime(dict(entry=entry, prime=3))
    data, = result['places']
    assert data['status'] == 'available'
    marks = {m['point_id']: m['components'] for m in data['specializations']}
    assert marks['P2'] == [3] and marks['P3'] == [1] and marks['P6'] == [2]
    graph = FibreGraph(data['genera'], [(a, b, 1) for a, b in data['edges']])
    charge = [QQ(0)]*len(data['genera'])
    for point in entry['support']:
        charge[marks[point['id']][0]] -= QQ(point['coefficients'].get('D2', 0))
    boundary = graph.boundary(charge)
    assert boundary.potential == vector(QQ, [0, -1, 2, -1])
    assert boundary.cycle_form == matrix(QQ, [[0]])
    assert not boundary.projected


@pytest.mark.parametrize('label', ['688.a.704512.1', '688.a.704512.2'])
def test_verified_square_completion_has_good_dyadic_reduction(label):
    catalog = json.loads((ROOT/'docs/_extra/explorer/catalog.json').read_text())
    entry = next(c for c in catalog['curves'] if c['label'] == label)
    result = compute_prime(dict(entry=entry, prime=2))
    assert result['status'] == 'available'
    for place in result['places']:
        assert place['genera'] == [2] and not place['edges']
        assert all(mark['components'] == [0] for mark in place['specializations'])


def test_exported_models_have_unit_edges_and_honest_specialization_status():
    catalog = json.loads((ROOT/'docs/_extra/explorer/catalog.json').read_text())
    for curve in catalog['curves']:
        record = json.loads((ROOT/'docs/_extra/explorer/reduction'/(curve['label']+'.json')).read_text())
        assert record['label'] == curve['label'] and record['schema_version'] == 1
        assert sorted(f['prime'] for f in record['fibres']) == sorted(record['bad_primes'])
        for fibre in record['fibres']:
            for place in fibre['places']:
                if place['status'] != 'available':
                    assert place['reason']
                    continue
                assert all(a != b for a, b in place['edges'])
                graph = FibreGraph(place['genera'], [(a, b, 1) for a, b in place['edges']])
                assert sum(graph.genera)+graph.cycles.nrows() == 2
                assert matrix(QQ, place['laplacian']) == graph.laplacian
                assert matrix(QQ, len(place['cycles']), len(place['edges']), place['cycles']) == graph.cycles
                assert [p['point_id'] for p in place['specializations']] == [p['id'] for p in curve['support']]
                for mark in place['specializations']:
                    assert all(0 <= i < len(graph.genera) for i in mark['components'])
                    assert (mark['status'] == 'known') == (len(mark['components']) == 1)
                marks = {mark['point_id']: mark for mark in place['specializations']}
                for divisor in curve['divisors']:
                    support = [p for p in curve['support'] if QQ(p['coefficients'].get(divisor['id'], 0))]
                    if any(marks[p['id']]['status'] != 'known' for p in support):
                        continue
                    charge = [QQ(0)]*len(graph.genera)
                    for p in support:
                        charge[marks[p['id']]['components'][0]] -= QQ(p['coefficients'][divisor['id']])
                    # These generators come from the saved certified integral bases.
                    assert not graph.boundary(charge).projected, (curve['label'], fibre['prime'], divisor['id'])
