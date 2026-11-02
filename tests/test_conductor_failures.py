"""Regressions for the conductor-997 and conductor-1116 survey failures."""

from sage.all import PolynomialRing, QQ

from hyperell_regulator.regulator import IntegrationMonitor, function_with_divisor
from hyperell_regulator.regulator.cover import Cover
from hyperell_regulator.regulator._mp_continuation import Arithmetic, FunctionSamples, quadrature
from hyperell_regulator.regulator._mp_stokes import Integrator
from hyperell_regulator.regulator.stokes import cut_data, cut_scale, hold_radii


def test_997_endpoint_approach_matches_analytic_square_root_continuation():
    x = PolynomialRing(QQ, 'x').gen()
    F = 4*x**5-32*x**4+64*x**3+x*x-4*x
    ar = Arithmetic(40)
    engine = Integrator(Cover(F, x-4, 0), ar)
    nodes, kinds, direction = engine.cov.cut_path()
    z = min(engine.critical_values, key=lambda v: abs(v-ar.C(nodes[-1])))
    a = z+ar.R(cut_scale(nodes)[0])*ar.C(direction)
    fraction = ar.R(hold_radii(nodes)[-1])/abs(a-z)
    b = z+fraction*(a-z)
    points = engine.fibre(a)
    lam = a.log()
    with IntegrationMonitor(max_edge_seconds=30, max_edge_steps=100000) as monitor:
        with monitor.edge('997 finite approach', degree=2):
            _, _, ordinary, logarithmic = engine.finite(a, b, points, lam, toward=(z, 2))

    # Independent sheet continuation: on this straight x-path, each root
    # factor has argument change strictly between -pi and pi. Their square
    # root product gives the analytic y branch without the adaptive walker
    # or any of its cached quadrature states.
    x0, y0 = points[0]
    end_scale = fraction.sqrt()
    def evaluate(t):
        scale = 1+(end_scale-1)*t
        w = z+(a-z)*scale**2
        xx = w+4
        yy = y0
        for root in engine.global_chart.roots:
            yy *= ((xx-root)/(x0-root)).sqrt()
        dx = 2*(a-z)*scale*(end_scale-1)
        values = [xx**j*dx/yy for j in range(2)]
        ll = lam+(w/a).log()
        return values+[ll*v for v in values]
    reference = quadrature(ar, FunctionSamples(evaluate), 4)
    for j in range(2):
        assert abs(ordinary[j][0]-reference[j]) < ar.R('1e-38')
        assert abs(logarithmic[j][0]-reference[j+2]) < ar.R('1e-38')
        assert abs(sum(ordinary[j])) < ar.R('1e-45')
    assert monitor.records[0]['steps'] < 100000


def test_1116_reference_labels_continue_the_point_near_infinity():
    x = PolynomialRing(QQ, 'x').gen()
    F = x**6+4*x**4+10*x**3+4*x*x-4*x+1
    cov = Cover(F, *function_with_divisor(F, (-1, -2), 'infinity', 13))
    nodes, kinds, direction = cov.cut_path()
    # Previously failed at sigma_27 even after four density refinements.
    upper, lower, _, jump = cut_data(cov, nodes, direction, tries=1)
    assert len(upper) == len(nodes) == len(lower)
    for permutation in upper+lower:
        assert sorted(permutation) == list(range(13))
    assert abs(abs(jump)-2*3.141592653589793) < 1e-12


def test_1116_even_infinity_ray_matches_logarithmic_reference():
    x = PolynomialRing(QQ, 'x').gen()
    F = x**6+4*x**4+10*x**3+4*x*x-4*x+1
    A, B = function_with_divisor(F, (-1, -2), 'infinity', 13)
    # The raw logarithmic integral and its float64 starting point depend on
    # phi's scalar. Preserve the scalar of this independent reference even
    # when the Riemann--Roch basis is subsequently normalized.
    leading = A.numerator().leading_coefficient()/A.denominator().leading_coefficient()
    scale = (QQ(41875329343225)/373360370073984)/leading
    cov = Cover(F, scale*A, scale*B)
    ar = Arithmetic(25)
    engine = Integrator(cov, ar)
    nodes, _, direction = cov.cut_path()
    z = min(engine.critical_values, key=lambda v: abs(v-ar.C(nodes[-1])))
    span, direction = ar.R(cut_scale(nodes)[0]), ar.C(direction)
    base = z+span*direction
    with IntegrationMonitor(max_edge_seconds=30, max_edge_steps=100000) as monitor:
        with monitor.edge('1116 infinity ray', degree=13):
            ordinary, logarithmic = engine.ray(z, direction, span, engine.fibre(base), base.log())
    # Independently computed with the previous, logarithmic parametrization
    # at 30 digits. Match by the integral rather than numerical fibre order.
    reference = [
        ar.C('0.198721357442399523171401443938713583', '-0.124330010029909223861306466082166028'),
        ar.C('-0.178308445582085952175036376582492277', '-0.796227881805596982976396880451753178'),
        ar.C('1.734392421735152593107331975021509029', '-1.842808791183529568017056624023410464'),
        ar.C('-5.558263179834582652864344113139157837', '-11.483185120815845226575013159914073235')]
    sheet = min(range(13), key=lambda j: abs(ordinary[0][j]-reference[0]))
    for row, expected in zip(ordinary+logarithmic, reference):
        assert abs(row[sheet]-expected) < ar.R('1e-24')
        # The trace of a holomorphic differential to the base sphere vanishes.
        assert abs(sum(row)) < ar.R('1e-24')
