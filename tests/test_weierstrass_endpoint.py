"""Logarithmic endpoint charts must preserve the order of vanishing of phi."""

import numpy as np

from hyperell_regulator.regulator import IntegrationMonitor, function_with_divisor, stokes
from hyperell_regulator.regulator._legs import _weierstrass_chart
from hyperell_regulator.regulator.cover import Cover
from hyperell_regulator.regulator.polynomials import ring


def _first_cover(shift=0):
    from sage.all import QQ

    x = ring().gen()
    F = x**6 + 4*x**5 + 4*x**4 + 2*x**3 + 1
    A, B = function_with_divisor(F, (-1, 0), "infinity", 7)
    argument = x - QQ(shift)
    return Cover(F(argument), A(argument), B(argument))


def test_weierstrass_logarithmic_valuation_survives_rounding():
    from sage.all import QQ

    # The translated case has endpoint 1/3, which is not a binary float.
    for shift in (0, QQ(4)/3):
        cov = _first_cover(shift)
        endpoint = complex(shift - 1)
        center, local = _weierstrass_chart(cov, endpoint + 2e-15 + 3e-16j)
        assert center == endpoint
        assert local.F[-1] == 0
        assert np.all(local.a[-4:] == 0)
        assert np.all(local.b[-3:] == 0)
        G = local.F[:-1]
        for y in (1e-6, 1e-10, 1e-20):
            u = y*y / G[-1]
            for _ in range(4):
                u = y*y / np.polyval(G, u)
            dy = -y
            dx = 2*y*dy / np.polyval(local.Fp, u)
            rate = local.dlogphi(u, y, dx, dy)
            # Along y=y0 exp(-s), a zero of order seven gives
            # d(log phi)/ds = -7 + O(y), with no rounded constant term.
            assert abs(rate + 7) < 20*y + 1e-13, (shift, y, rate)


def test_first_edge_logarithmic_trace_closes():
    cov = _first_cover()
    nodes, kinds, direction = cov.cut_path()
    _, _, lam, _ = stokes.cut_data(cov, nodes, direction)
    holds = stokes.hold_radii(nodes)
    with IntegrationMonitor(max_edge_seconds=15) as monitor:
        _, leg = stokes.edge_leg(cov, nodes[0], nodes[1], lam[0], kinds[0], kinds[1],
                                hold_a=holds[0], hold_b=holds[1])
    # Previously the period trace was ~2e-11 but the logarithmic trace
    # was 2.3e-4, from rounded high-order zeros in the y-chart at P.
    assert np.max(np.abs(leg.dI.sum(axis=1))) < 1e-10
    assert np.max(np.abs(leg.dL.sum(axis=1))) < 1e-9
    assert monitor.records[0]["status"] == "ok"
    print("P7 edge 0:", monitor.records[0]["seconds"], "seconds;",
          monitor.records[0]["counters"], "traceI", np.max(np.abs(leg.dI.sum(axis=1))),
          "traceL", np.max(np.abs(leg.dL.sum(axis=1))))
