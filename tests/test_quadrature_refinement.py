"""Sheet-local convergence must preserve continuation data and sheet labels."""

from unittest.mock import patch
import numpy as np

from hyperell_regulator.regulator import _legs, IntegrationMonitor


def test_quadrature_refines_only_unresolved_sheets():
    # The last two sheets have the same difficult analytic factor, once in
    # I and once in L.  The first sheet must stop after the first comparison.
    width, center = 0.1, 0.43
    exact = width * (np.arctan((1-center)/width) - np.arctan(-center/width))

    class Cover:
        g = 2

    for logarithmic in (False, True):
        calls, cache_ids = [], {}
        charts = [object() for _ in range(3)]

        def integrate(cov, coord, rate, x, y, lam, nodes, wts, **kwargs):
            ids = kwargs["sheet_ids"]
            calls.append(ids.copy())
            t = -np.expm1(-nodes) if logarithmic else nodes
            weights = wts * np.exp(-nodes) if logarithmic else wts
            smooth = np.sum(weights)
            hard = np.sum(weights / (1 + ((t-center)/width)**2))
            I, L = np.empty((2, len(ids)), complex), np.empty((2, len(ids)), complex)
            for local, original in enumerate(ids):
                assert x[local] == original
                assert y[local] == original + 5
                assert lam[local] == original + 10
                assert kwargs["chosen"][local] is charts[original]
                cache = kwargs["samples"][local]
                if original in cache_ids:
                    assert id(cache) == cache_ids[original]
                    assert cache["sheet"] == original
                else:
                    cache_ids[original] = id(cache)
                    cache["sheet"] = original
                scale = np.array([1, original + 2])
                I[:, local] = scale * (hard if original == 2 else smooth)
                L[:, local] = scale * (hard if original == 1 else (original+2)*smooth)
            return x+20, y+20, lam+20, I, L

        with patch.object(_legs.walks, "choose", return_value=charts):
            with patch.object(_legs.walks, "leg", integrate):
                out = _legs._walk_quadrature(
                    Cover(), lambda t: t, lambda t: 1, np.arange(3),
                    np.arange(3)+5, np.arange(3)+10, charts, level=3,
                    logarithmic=logarithmic)
        assert calls == [[0, 1, 2], [0, 1, 2], [1, 2], [1, 2]]
        for value, offset in zip(out[:3], (20, 25, 30)):
            np.testing.assert_array_equal(value, np.arange(3)+offset)
        scale = np.array([[1, 1, 1], [2, 3, 4]])
        np.testing.assert_allclose(out[3], scale*np.array([1, 1, exact]), rtol=1e-12, atol=1e-13)
        np.testing.assert_allclose(out[4], scale*np.array([2, exact, 4]), rtol=1e-12, atol=1e-13)


def test_sheetwise_quadrature_matches_small_curve_plane_reference():
    from hyperell_regulator.regulator import periods_and_regulator_stokes
    from hyperell_regulator.regulator.polynomials import ring

    x = ring().gen()
    # Stored independent plane integrals; their displayed precision is 1e-8.
    with IntegrationMonitor(max_edge_seconds=15) as monitor:
        periods, regulator = periods_and_regulator_stokes(x**5+x*x+x, (x, 0))
    np.testing.assert_allclose(periods, [9.22548668, -1.22108503, 8.98356248], atol=2e-8, rtol=0)
    np.testing.assert_allclose(regulator, [-5.80825813, 0.45209334, 5.53626738], atol=2e-8, rtol=0)
    assert monitor.records and all(record["status"] == "ok" for record in monitor.records)


def test_nearby_complex_singularity_can_require_more_than_four_levels():
    """An analytic integrand with known integral and nearby complex poles."""
    width, center = 0.012, 0.43
    exact = width*(np.arctan((1-center)/width)+np.arctan(center/width))

    class Cover:
        g = 1

    counts = []
    def integrate(cov, coord, rate, x, y, lam, nodes, wts, **kwargs):
        counts.append(len(nodes))
        value = np.sum(wts/(1+((nodes-center)/width)**2))
        return x, y, lam, np.array([[value]]), np.array([[2*value]])

    with patch.object(_legs.walks, 'choose', return_value=[object()]):
        with patch.object(_legs.walks, 'leg', integrate):
            out = _legs._walk_quadrature(Cover(), lambda t:t, lambda t:1,
                                         [0], [1], [0], [], level=5)
            assert abs(out[3][0,0]-exact) < 1e-12
            assert counts[-1] > len(_legs._tanh_sinh(8)[0])
            try:
                _legs._walk_quadrature(Cover(), lambda t:t, lambda t:1,
                                       [0], [1], [0], [], level=5, max_level=8)
            except _legs.QuadratureFailure as exc:
                assert 'last scaled period/log errors' in str(exc)
            else:
                raise AssertionError('expected quadrature resource limit')


def test_quadrature_failure_does_not_retry_endpoint_geometry():
    class Cover:
        def endpoint_ramification(self, z, kind):
            return 2

    error = _legs.QuadratureFailure('unresolved integral')
    with patch.object(_legs, 'leg_in_w', side_effect=error) as integrate:
        with patch.object(_legs, 'certified_reach') as retry:
            try:
                _legs.side(Cover(), 1, [1], [1], [0], 0, 'simple', False)
            except _legs.QuadratureFailure as exc:
                assert exc is error
            else:
                raise AssertionError('quadrature failure was swallowed')
            assert integrate.call_count == 1
            retry.assert_not_called()


def test_degree_three_ramification_refines_both_endpoint_lifts():
    from hyperell_regulator.regulator import function_with_divisor
    from hyperell_regulator.regulator.cover import Cover
    from hyperell_regulator.regulator.polynomials import ring

    x = ring().gen()
    F = x**6 - 2*x**5 - x**4 + 4*x**3 + 3*x**2 + 2*x + 1
    cov = Cover(F, *function_with_divisor(F, (0, -1), 'infinity', 3))
    z = min(cov.critical_values(), key=lambda w: abs(w-1.55133))
    xs, ys = cov.fibre(z+0.001j)
    targets = [_legs._target(cov, a, b, z) for a, b in zip(xs, ys)]
    # Two lifts must end at the same critical point, accurately, and the
    # third lift must retain its distinct, unramified endpoint.
    i, j = min(((i,j) for i in range(3) for j in range(i)),
               key=lambda ij: abs(targets[ij[0]][0]-targets[ij[1]][0]))
    assert abs(targets[i][0]-targets[j][0]) < 1e-12
    k = next(k for k in range(3) if k not in (i,j))
    assert abs(targets[k][0]-targets[i][0]) > 0.1
    for a,b in targets:
        assert abs(complex(cov.phi(np.array([a]),np.array([b]))[0])-z) < 1e-12
