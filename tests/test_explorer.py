"""Numerical/export checks independent of the browser's plot rendering."""

import json
import hashlib
from pathlib import Path

import numpy as np
import pytest
from sage.all import PolynomialRing, QQ

from hyperell_regulator.regulator._mp_continuation import Arithmetic, LogFunctionSamples, quadrature
from hyperell_regulator.regulator.explorer import parse_request


def test_endpoint_logarithm_is_continued_across_principal_cut():
    ar = Arithmetic(25)
    angular_rate = ar.C.gen()*5*ar.R.pi()
    def evaluate(t):
        return [ar.C(1)], (angular_rate*t).exp().log()
    result = quadrature(ar, LogFunctionSamples(ar, evaluate, ar.C(0)), 2)
    assert abs(result[0]-1) < ar.R('1e-25')
    assert abs(result[1]-angular_rate/2) < ar.R('1e-25')


@pytest.mark.parametrize('expression', ["__import__('os').getcwd()", "x.__class__", "[x]", "x**10000"])
def test_local_export_parser_accepts_only_rational_arithmetic(expression):
    with pytest.raises(ValueError):
        parse_request(dict(F=expression, A="x", B="0", prec=15))


def test_local_export_rejects_singular_curve_and_wrong_precision():
    with pytest.raises(ValueError, match='squarefree'):
        parse_request(dict(F="x^5", A="x", B="0", prec=15))
    with pytest.raises(ValueError, match='prec'):
        parse_request(dict(F="x^5+x+1", A="x", B="0", prec=3))


def _complex(z):
    return complex(float(z[0]), float(z[1]))


def _assemble_export(data, logarithmic):
    """Stokes via complete sheet boundaries reconstructed from exported banks.

    This checks the serialization's sheet permutations, detour aggregation,
    bank jump and infinite edge against the independently exported pairings.
    """
    edges, g, n = data['edges'], data['genus'], data['degree']
    r = len(edges)
    u = np.array([[[_complex(lift['ordinary'][a]) for lift in e['lifts']]
                   for e in edges] for a in range(g)])
    v = np.array([[[_complex(lift['log_minus'][a]) for lift in e['lifts']]
                   for e in edges] for a in range(g)]) if logarithmic else u
    jump = _complex(data['jump']) if logarithmic else 0j
    boundaries = []
    for sheet in range(n):
        forward = [(j, sheet, 1) for j in range(r)]
        backward = [(j, next(i for i, lift in enumerate(edges[j]['lifts'])
                             if lift['minus_sheet'] == sheet+1), -1) for j in reversed(range(r))]
        boundaries.append(forward+backward)
    H = np.zeros((g, g), complex)
    for a in range(g):
        primitives = []
        starts, ends = {}, {}
        for walk in boundaries:
            values = [0j]
            for j, s, sign in walk:
                (starts if sign == 1 else ends)[j, s] = values[-1]
                values.append(values[-1]+(v[a, j, s]+jump*u[a, j, s] if sign == 1 else -v[a, j, s]))
            primitives.append(values)
        assert max(abs(p[-1]) for p in primitives) < 1e-10*max(1, abs(v).max())
        for j in range(r):
            for s in range(n):
                coefficient = starts[j, s]-ends[j, s]+v[a, j, s]
                H[a] += coefficient*np.conj(u[:, j, s])
    matrix = (H.imag+H.imag.T)/2
    if logarithmic:
        matrix += jump.imag/2*np.real(np.einsum('aim,bim->ab', u, np.conj(u)))
    return np.array([-matrix[a, b]/4 for a in range(g) for b in range(a, g)])


def test_complete_small_conductor_catalog_reassembles_all_pairings():
    root = Path(__file__).resolve().parents[1]
    site = root/'docs/_extra/explorer'
    catalog = json.loads((site/'catalog.json').read_text())
    inputs = json.loads((root/'reports/small_conductor_explorer_inputs.json').read_text())
    assert [c['label'] for c in catalog['curves']] == [c['label'] for c in inputs['curves']]
    assert len(catalog['curves']) == 16 and catalog['display_digits'] == 10
    assert catalog['source_sha256']
    assert all(len(checksum) == 64 for checksum in catalog['source_sha256'].values())
    for curve in catalog['curves']:
        assert curve['validation']['determinant_rounds_as_printed']
        assert len(curve['functions']) == 2
        assert [f['name'] for f in curve['functions']] == ['P', 'R']
        for f in curve['functions']:
            data = json.loads((site/f['file']).read_text())
            assert data['degree'] == f['order']
            assert data['collection']['pole'] == 'Q'
            assert data['collection']['function'] == f['name']
            assert data['source_sha256'] == catalog['source_sha256']
            assert not data['warnings'], (curve['label'], f['name'], data['warnings'])
            for edge in data['edges']:
                assert all(lift['xyz'] and None not in lift['xyz'] for lift in edge['lifts'])
                for lift in edge['lifts']:
                    for a in range(data['genus']):
                        difference = _complex(lift['log_plus'][a]) - _complex(lift['log_minus'][a])
                        expected = _complex(data['jump'])*_complex(lift['ordinary'][a])
                        assert abs(difference-expected) < 1e-10*max(1,abs(expected))
            for log, key in ((False, 'period'), (True, 'regulator')):
                actual = _assemble_export(data, log)
                expected = [float(p[key]) for p in data['pairings']]
                np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-10,
                                           err_msg=curve['label']+' '+f['name'])


def test_all_surface_meshes_match_the_unchanged_integral_data():
    site = Path(__file__).resolve().parents[1]/'docs/_extra/explorer'
    catalog = json.loads((site/'catalog.json').read_text())
    for curve in catalog['curves']:
        for function in curve['functions']:
            integral_bytes = (site/function['file']).read_bytes()
            mesh = json.loads((site/function['mesh']).read_text())
            assert hashlib.sha256(integral_bytes).hexdigest() == mesh['integral_data_sha256']
            assert mesh['degree'] == function['order']
            assert len(mesh['sheets']) == function['order']
            triangles = np.array(mesh['triangles'])
            assert len(triangles) > 100 and triangles.shape[1] == 3
            for index, sheet in enumerate(mesh['sheets']):
                assert sheet['sheet'] == index+1
                vertices = np.array(sheet['vertices'])
                assert vertices.ndim == 2 and vertices.shape[1] == 3
                assert np.isfinite(vertices).all()
                assert triangles.min() >= 0 and triangles.max() < len(vertices)
            checks = mesh['validation']
            assert checks['checked_continuation_cycles'] >= 80
            assert checks['triangles_per_sheet'] == len(triangles)
            for key in ('max_relative_sheet_cycle_defect', 'max_relative_sheet_boundary_defect',
                        'max_relative_curve_residual', 'max_relative_function_residual'):
                assert checks[key] < 1e-6, (curve['label'], function['name'], key)
