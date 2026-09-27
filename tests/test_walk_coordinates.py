"""Regressions for endpoint charts and resumable continuation."""
import numpy as np

from hyperell_regulator.regulator import function_with_divisor, IntegrationMonitor
from hyperell_regulator.regulator.cover import Cover
from hyperell_regulator.regulator.polynomials import ring
from hyperell_regulator.regulator.walks import Chart, charts_over, race, walker

x = ring().gen()
F249 = x**6 + 4*x**5 + 4*x**4 + 2*x**3 + 1


def test_inverse_differential_cancels_weierstrass_singularity():
    cov = Cover(F249, *function_with_divisor(F249, (0, -1), "infinity", 14))
    for sigma in (False, True):
        chart = Chart(cov, "inv", 0, sigma=sigma)
        zeta = -0.5 if sigma else -2
        num, den = chart.integrand(-1, 0, zeta)
        dzdw = -zeta*zeta if sigma else 1
        expected = np.array([1, -1]) / cov.y_dphi_dx(-1, 0)
        assert np.all(np.isfinite(num / den))
        assert np.max(abs(num / den * dzdw - expected)) < 1e-14


def test_inverse_differential_stays_scaled_at_infinity():
    cov = Cover(x**6+1, x**3, 1)
    chart = Chart(cov, "inv", 0, sigma=True)
    # phi = (1+eta)/u^3, eta^2=1+u^6. A rational derivative
    # would overflow at this u, but d(1/phi)/du is 3u^2/2.
    u = 1e-100
    eta = 1.0
    sigma = u**3 / 2
    num, den = chart.integrand(u, eta, sigma)
    assert abs(den / (1.5*u*u) - 1) < 1e-14
    assert abs(num[1] + 1) < 1e-14


def test_chart_race_resumes_the_same_walk():
    cov = Cover(x**5+(x+1)**2, x+1, 1)
    xx, yy = cov.fibre(1+1j)
    chart = Chart(cov, "x", 0)
    coord = lambda t: 1+1j + t*(0.8+0.3j)
    with IntegrationMonitor() as monitor:
        with monitor.edge("test", degree=cov.N):
            xe, ye, steps, _ = race(cov, coord, xx[0], yy[0], [chart],
                                    budgets=(1, 2, 8, 64, 256))
    counts = monitor.records[0]["counters"]
    assert counts["chart_trials"] > 1
    assert counts["walk_attempts"] == steps
    assert abs(cov.phi_stable(xe, ye) - coord(1)) < 1e-12
    direct = walker(cov, coord, xx[0], yy[0], chart)
    direct.to(1)
    assert abs(direct.u-xe) < 1e-12


def test_infinity_pole_uses_endpoint_chart_only():
    cov = Cover(F249, *function_with_divisor(F249, (0, -1), "infinity", 14))
    charts = charts_over(cov, None, sigma=True)
    assert len(charts) == 1 and charts[0].kind == "inv"


def test_inverse_chart_can_change_center_when_path_crosses_its_pole():
    cov = Cover(x**5+x*x+1, x, 0)
    xx, yy = cov.fibre(-0.1)
    coord = lambda t: -0.1+0.2*t
    chart = Chart(cov, "inv", 0)
    xe, ye, _, chosen = race(cov, coord, xx[0], yy[0], [chart],
                             budgets=(8, 32, 128))
    assert chosen.kind == "inv" and chosen.c != 0
    assert abs(xe-0.1) < 1e-12
    assert abs(ye*ye-complex(cov.F[0]*xe**5+xe*xe+1)) < 1e-12
