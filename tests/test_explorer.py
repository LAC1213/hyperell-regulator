"""Numerical/export checks independent of the browser's plot rendering."""

import json
import hashlib
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
from sage.all import PolynomialRing, QQ

from hyperell_regulator.regulator.explorer import parse_request
from hyperell_regulator.regulator.surfaces import _bank_radius, _ray_endpoint


def test_bank_probe_stays_between_long_nearby_cut_segments():
    # At the midpoint the other bank is only one unit away, although every
    # branch value is roughly half a billion units away (the 603 geometry).
    cut = np.array([0j, 1e9+0j, 2j, -1e9+2j])
    midpoint = 5e8+0j
    radius = _bank_radius(midpoint, cut[:-1], list(zip(cut[:-1],cut[1:])), 0)
    assert 0 < radius < .5


def test_diagonal_ray_is_attached_exactly_to_the_mesh_boundary():
    start = np.complex128(-0.3631747176028129-0.5361912224070701j)
    direction = complex(0.5850690806785093, -0.8109834590384712)
    bounds = [-1.8954284433231354, 1.3755735373634823,
              -2.7982054841030943, 0.9197943833330182]
    end = _ray_endpoint(start, direction, bounds)
    # A gap of even one ulp lets a mesh cycle go around all five branch
    # values, reversing a hyperelliptic sheet (conductor 909).
    assert end.imag == bounds[2]
    assert bounds[0] <= end.real <= bounds[1]
    assert abs(((end-start)/direction).imag) < 1e-14


@pytest.mark.parametrize('expression', ["__import__('os').getcwd()", "x.__class__", "[x]", "x**10000"])
def test_local_export_parser_accepts_only_rational_arithmetic(expression):
    with pytest.raises(ValueError):
        parse_request(dict(F=expression, centers=["0"], prec=15))


def test_local_export_rejects_singular_curve_and_wrong_precision():
    with pytest.raises(ValueError, match='squarefree'):
        parse_request(dict(F="x^5", centers=["0"], prec=15))
    with pytest.raises(ValueError, match='prec'):
        parse_request(dict(F="x^5+x+1", centers=["0"], prec=3))


def _complex(z):
    return complex(float(z[0]), float(z[1]))


def _assemble_export(data, logarithmic, weight=0):
    """Stokes via complete sheet boundaries reconstructed from exported banks.

    This checks the serialization's sheet permutations, detour aggregation,
    bank jump and infinite edge against the independently exported pairings.
    """
    edges, g, n = data['edges'], data['genus'], data['degree']
    r = len(edges)
    u = np.array([[[_complex(lift['ordinary'][a]) for lift in e['lifts']]
                   for e in edges] for a in range(g)])
    v = np.array([[[_complex(lift['logs'][weight]['minus'][a]) for lift in e['lifts']]
                   for e in edges] for a in range(g)]) if logarithmic else u
    jumps = [_complex(e['jumps'][weight]) if logarithmic else 0j for e in edges]
    boundaries = []
    word=data.get('boundary_word',[(j,1) for j in range(r)]+[(j,-1) for j in reversed(range(r))])
    for sheet in range(n):
        walk=[]
        for j,side in word:
            s=sheet if side>0 else next(i for i,lift in enumerate(edges[j]['lifts']) if lift['minus_sheet']==sheet+1)
            walk.append((j,s,side))
        boundaries.append(walk)
    H = np.zeros((g, g), complex)
    for a in range(g):
        primitives = []
        starts, ends = {}, {}
        for walk in boundaries:
            values = [0j]
            for j, s, sign in walk:
                (starts if sign == 1 else ends)[j, s] = values[-1]
                values.append(values[-1]+(v[a, j, s]+jumps[j]*u[a, j, s] if sign == 1 else -v[a, j, s]))
            primitives.append(values)
        assert max(abs(p[-1]) for p in primitives) < 1e-10*max(1, abs(v).max())
        for j in range(r):
            for s in range(n):
                coefficient = starts[j, s]-ends[j, s]+v[a, j, s]
                H[a] += coefficient*np.conj(u[:, j, s])
    matrix = (H.imag+H.imag.T)/2
    if logarithmic:
        matrix += sum(jumps[j].imag/2*np.real(np.einsum('am,bm->ab', u[:,j], np.conj(u[:,j]))) for j in range(r))
    return np.array([-matrix[a, b]/4 for a in range(g) for b in range(a, g)])


def test_catalog_has_one_two_sheet_plot_per_curve_and_reassembles_all_weights():
    root=Path(__file__).resolve().parents[1];site=root/'docs/_extra/explorer'
    catalog=json.loads((site/'catalog.json').read_text())
    inputs=json.loads((root/'docs/explorer_inputs.json').read_text())
    assert catalog['schema_version']==2
    assert [c['label'] for c in catalog['curves']]==[c['label'] for c in inputs['curves']]
    assert len(catalog['curves'])==21
    assert len({c['file'] for c in catalog['curves']})==21
    assert catalog['source']['report_sha256']==hashlib.sha256((root/'reports/beilinson_100_digits.json').read_bytes()).hexdigest()
    for curve in catalog['curves']:
        assert 'functions' not in curve and 'exact_functions' not in curve
        assert float(curve['validation']['determinant_relative_error'])<1e-10
        data=json.loads((site/curve['file']).read_text())
        assert data['schema_version']==2 and data['degree']==2
        assert data['input']['cut']=='mst'
        assert len(data['boundary_word'])==2*len(data['edges'])
        assert len(data['edges'])==len(data['nodes'])
        assert len(data['weights'])>=2
        expected_periods=[float(p['period']) for p in data['weights'][0]['pairings']]
        np.testing.assert_allclose(_assemble_export(data,False),expected_periods,atol=1e-10,rtol=1e-10)
        for k,weight in enumerate(data['weights']):
            expected=[float(p['regulator']) for p in weight['pairings']]
            np.testing.assert_allclose(_assemble_export(data,True,k),expected,atol=1e-10,rtol=1e-10,err_msg=curve['label'])
            for edge in data['edges']:
                for lift in edge['lifts']:
                    for a in range(data['genus']):
                        delta=_complex(lift['logs'][k]['plus'][a])-_complex(lift['logs'][k]['minus'][a])
                        expected_delta=_complex(edge['jumps'][k])*_complex(lift['ordinary'][a])
                        assert abs(delta-expected_delta)<1e-10*max(1,abs(expected_delta))
        mesh=json.loads((site/curve['mesh']).read_text())
        assert mesh['degree']==2 and len(mesh['sheets'])==2
        assert mesh['integral_data_sha256']==hashlib.sha256((site/curve['file']).read_bytes()).hexdigest()
        for key in ('max_relative_sheet_cycle_defect','max_relative_sheet_boundary_defect','max_relative_curve_residual'):
            assert mesh['validation'][key]<1e-6,(curve['label'],key)
        matrix=np.array(curve['normalized_matrix'],float)
        assert abs(np.linalg.det(matrix)/float(curve['determinant'])-1)<1e-10
