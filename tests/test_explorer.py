"""Numerical/export checks independent of the browser's plot rendering."""

import json
import hashlib
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest
from sage.all import PolynomialRing, QQ

from hyperell_regulator.regulator._mp_continuation import Arithmetic, LogFunctionSamples, quadrature
from hyperell_regulator.regulator.explorer import parse_request
from hyperell_regulator.regulator.surfaces import _bank_radius, _label_seed, _ray_endpoint


def test_bank_probe_stays_between_long_nearby_cut_segments():
    # At the midpoint the other bank is only one unit away, although every
    # branch value is roughly half a billion units away (the 603 geometry).
    cut = np.array([0j, 1e9+0j, 2j, -1e9+2j])
    midpoint = 5e8+0j
    radius = _bank_radius(midpoint, cut[:-1], cut, 0)
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


def test_mesh_seed_matches_the_integrators_racetrack():
    from hyperell_regulator.regulator.stokes import racetrack
    nodes = np.array([0j, 1+2j, -1+3j])
    assert _label_seed(nodes) == racetrack(nodes, 1j)[0][0]


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
    inputs = json.loads((root/'docs/explorer_inputs.json').read_text())
    assert [c['label'] for c in catalog['curves']] == [c['label'] for c in inputs['curves']]
    report = json.loads((root/'reports/beilinson_100_digits.json').read_text())
    successful = {r['label']: r for r in report if r['status'] == 'ok' and int(r['cond']) != 1116}
    assert any(int(r['cond']) == 1116 for r in report)
    assert all(int(c['cond']) != 1116 for c in catalog['curves'])
    assert inputs['source']['excluded_conductors'] == [1116]
    assert [c['label'] for c in catalog['curves']] == list(successful)
    assert catalog['display_digits'] == 10
    assert catalog['source']['report_sha256'] == hashlib.sha256((root/'reports/beilinson_100_digits.json').read_bytes()).hexdigest()
    assert catalog['source_sha256']
    assert catalog['conventions']['determinant_weights'] == [1, 2, 1]
    assert all(len(checksum) == 64 for checksum in catalog['source_sha256'].values())
    for curve in catalog['curves']:
        reference = successful[curve['label']]
        assert curve['cond'] < 1300
        assert curve['numerical']['det'] == reference['det']
        assert curve['numerical']['L_second_derivative'] == reference.get('L_second_derivative')
        assert curve['regulator_digits'] == 100
        if curve.get('visualization_status') == 'not_exported':
            assert curve['label'] in catalog['reference_only']
            assert curve['functions'] == []
            continue
        assert float(curve['validation']['determinant_relative_error']) < 1e-10
        assert len(curve['functions']) == 2
        assert [f['name'] for f in curve['functions']] == ['P', 'R']
        if curve['comparison_verified']:
            assert curve['q'] == reference['comparison']['rational']
            if curve['integrality_verified']:
                assert Fraction(curve['q']).denominator == 1
        assert curve['integrality_verified'] == all(reference['triple'].get(k) is True
                                                   for k in ('verdict_P', 'verdict_Q'))
        matrix = np.array(curve['normalized_matrix'], float)
        assert abs(np.linalg.det(matrix)/float(curve['determinant'])-1) < 1e-10
        for f in curve['functions']:
            data = json.loads((site/f['file']).read_text())
            assert data['degree'] == f['order']
            assert data['collection']['pole'] == 'Q'
            assert data['collection']['function'] == f['name']
            assert data['normalization_source_sha256'] == catalog['normalization_source_sha256']
            for module, checksum in data['source_sha256'].items():
                if module.startswith('hyperell_regulator/regulator/') and not module.endswith(('/api.py', '/plane.py', '/surfaces.py')):
                    assert checksum == catalog['source_sha256'][module]
            for w, pair in zip((1, 2, 1), data['pairings']):
                assert abs(float(pair['determinant_period'])-w*float(pair['period'])) < 1e-12
            column = curve['functions'].index(f)
            np.testing.assert_allclose(matrix[:, column], [float(p['normalized_regulator']) for p in data['pairings']], rtol=1e-12)
            np.testing.assert_allclose(matrix[:, 2], [float(p['determinant_period']) for p in data['pairings']], rtol=1e-10)
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
    missing = []
    for curve in catalog['curves']:
        for function in curve['functions']:
            if 'mesh' not in function:
                assert function.get('mesh_error')
                missing.append((curve['label'], function['name']))
                assert (site/function['file']).is_file()
                continue
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
    # These nearly coincident cuts remain usable through their exact
    # functions, paths and integrals. At 997 the real-axis detour has only
    # a roughly 1e-13 imaginary separation, disconnecting the float mesh.
    assert set(missing) <= {('997.a.997.1', 'P'), ('1147.a.35557.2', 'P')}
