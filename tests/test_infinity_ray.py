"""Odd-degree rays must end in a genuine infinity uniformizer."""
from unittest.mock import patch
import numpy as np
from hyperell_regulator.regulator import IntegrationMonitor, function_with_divisor
from hyperell_regulator.regulator import _legs, stokes
from hyperell_regulator.regulator.cover import Cover
from hyperell_regulator.regulator.polynomials import ring


def test_709_ray_is_independent_of_infinity_handover():
    x = ring().gen()
    F = 4*x**5-7*x**4+4*x
    cov = Cover(F, *function_with_divisor(F, (1,-1), 'infinity', 8))
    z, kinds, d = cov.cut_path()
    _, _, lam, _ = stokes.cut_data(cov, z, d)
    hold = stokes.hold_radii(z)[-1]
    radius, bound = _legs._odd_infinity_region(cov)
    radius2 = 2*radius
    bound2 = (np.polyval(abs(cov.a),radius2)/abs(cov.ad[0])
              + np.polyval(abs(cov.b),radius2)/abs(cov.bd[0])
              * np.sqrt(np.polyval(abs(cov.F),radius2)))
    legs = []
    for region in ((radius,bound),(radius2,bound2)):
        with patch.object(_legs, '_odd_infinity_region', return_value=region):
            with IntegrationMonitor(max_edge_seconds=15) as monitor:
                _, leg = _legs.ray_leg(cov,z[-1],d,lam[-1],kind_in=kinds[-1],
                                       span=stokes.cut_scale(z)[0],hold=hold)
        legs.append(leg)
        assert max(abs(leg.dI.sum(axis=1))) < 1e-8
        assert max(abs(leg.dL.sum(axis=1))) < 1e-8
        # The previous infinity chart search exhausted 3*4096 attempts
        # on the first sheet. All sheets together should now take <2000.
        searches = [p for p in monitor.records[0]['phases']
                    if p['path'].endswith('odd_infinity_ray/walk_quadrature/choose_charts')]
        assert searches
        assert sum(p['counters']['walk_attempts'] for p in searches) < 2000
    np.testing.assert_allclose(legs[0].dI,legs[1].dI,atol=1e-9,rtol=0)
    np.testing.assert_allclose(legs[0].dL,legs[1].dL,atol=1e-9,rtol=0)
