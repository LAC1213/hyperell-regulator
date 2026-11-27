"""Verified auxiliary models for specific wild genus 2 equations at 2.

The chart parameters refer to equations, independently of database labels.
Every chart used for a constraint is checked by exact arithmetic.
"""

from functools import lru_cache
from math import factorial

from sage.all import Infinity, NumberField, PolynomialRing, QQ, matrix

from ._projected_reduction import LocalConstraints, graph_constraints, node_thickness
from .graph_boundary import FibreGraph
from .stable_model import StableModel


# Reciprocal coordinate, depth of the elliptic chart centre relative to
# the singular abscissa, and disc radius, in the normalization v(2)=1.
_ELLIPTIC_CHARTS = {
    ((0,-1,0,0,0,1), ()): (False, QQ(1)/2, QQ(1)),
    ((0,-1,1,1,-3,2), (1,)): (True, QQ(1)/2, QQ(3)/2),
    ((0,0,-2,-4,3,8), (1,1)): (True, QQ(3)/2, QQ(5)/3),
    ((0,1,16,72,33,4), (0,1)): (True, QQ(1), QQ(5)/3),
    ((0,0,0,-1,-2,-2,-1), (1,1)): (True, QQ(5)/8, QQ(2)/3),
    ((0,1,0,-2,-1,1), (1,0,1)): (False, QQ(1)/2, QQ(2)/3),
    ((0,-1,0,-1), (1,1,1,1)): (False, Infinity, QQ(1)),
    ((0,0,-2,-1,4,4), (1,)): (True, QQ(1), QQ(2)),
    ((0,1,0,-3,-1,2), (1,)): (True, QQ(1)/2, QQ(3)/2),
    ((-1,-1,-1,-1), (1,1,1,1)): (False, QQ(1)/2, QQ(3)/2),
    ((0,1,-1,-4,0,3), (1,0,1)): (False, QQ(3)/4, QQ(5)/6),
}
_INTERMEDIATE_CHART = ((0,1,-1,-4,0,3), (1,0,1))
_SPLIT_ANNULUS = ((0,0,-1,2,0,-2,1), (1,1))
_TWO_BRANCH_LOOP = ((0,1,0,-1,-1), (0,1,0,1))
_BOLZA = ((0,-1,0,0,0,1), ())


def _valuation(a, prime):
    return a.valuation(prime)/prime.ramification_index() if a else Infinity


@lru_cache(maxsize=32)
def elliptic_chart(equation, reciprocal, depth, radius):
    """Verify a square completion whose residual Artin–Schreier curve has genus one."""
    R = PolynomialRing(QQ, 't')
    t = R.gen()
    f, h = [R(list(coeffs)) for coeffs in equation]
    F = h*h+4*f
    if reciprocal:
        F = R([F[6-i] for i in range(7)])
    A = [F.derivative(i)/factorial(i) for i in range(7)]
    V = 4*F*A[2]-A[1]**2
    singular = 0 if reciprocal else 1
    factors = ([t-singular] if depth == Infinity else
               [g for g, _ in (64*F**3*A[4]-V**2).factor()])
    for factor in factors:
        K = NumberField(factor.monic(), 'c')
        c = K.gen()
        a = [g(c) for g in A]
        if not a[0]:
            continue
        b = a[1]/(2*a[0])
        d = (a[2]/a[0]-b*b)/2
        G = (a[3]-2*a[0]*b*d, a[4]-a[0]*d*d, a[5], a[6])
        for prime in K.primes_above(2):
            if _valuation(c-singular, prime) != depth:
                continue
            threshold = _valuation(a[0], prime)+2
            w = [_valuation(g, prime)+i*radius for i, g in enumerate(G, 3)]
            if (min(w) != threshold or _valuation(b, prime)+radius <= 0
                    or _valuation(d, prime)+2*radius <= 0):
                continue
            support = {i for i, v in enumerate(w, 3) if v == threshold}
            if 5 in support or not {3, 6} & support:
                continue
            if {3, 6} <= support and _valuation(G[0]**2/(4*a[0])-G[3], prime) != _valuation(G[3], prime):
                continue
            return c, prime, (a[0], b, d), G
    raise NotImplementedError('no verified elliptic chart at the specified radius')


def _depth(x, centre, prime, reciprocal):
    if reciprocal:
        x = 0 if x == 'infinity' else (1/x if x else 'infinity')
    if x == 'infinity':
        return -Infinity
    K = centre.parent()
    if x in QQ:
        return _valuation(K(QQ(x))-centre, prime)

    # If the Newton polygon of the shifted minimal polynomial has one
    # slope, the distance is independent of the local embedding of x.
    R = PolynomialRing(K, 'T')
    g = R(x.minpoly().list())
    shifted = g(R.gen()+centre)
    if not shifted[0]:
        return None
    d = shifted.degree()
    v = _valuation(shifted[0], prime)/d
    if all(_valuation(shifted[i], prime) >= (d-i)*v for i in range(1, d)):
        return v
    return None


def _mean_clipped_depth(polynomial, centre, prime, radius, floor=QQ(0)):
    """Average the clipped root distances using the exact lower Newton polygon."""
    R = PolynomialRing(centre.parent(), 'z')
    shifted = R(polynomial.list())(R.gen()+centre)
    vertices = []
    for i, coefficient in enumerate(shifted):
        if not coefficient:
            continue
        point = (i, _valuation(coefficient, prime))
        while len(vertices) > 1:
            a, b = vertices[-2:]
            if (b[1]-a[1])*(point[0]-b[0]) < (point[1]-b[1])*(b[0]-a[0]):
                break
            vertices.pop()
        vertices.append(point)
    total = vertices[0][0]*radius  # Roots equal to the centre have infinite depth.
    for (i, a), (j, b) in zip(vertices, vertices[1:]):
        depth = (a-b)/(j-i)
        total += (j-i)*max(floor, min(radius, depth))
    return total/polynomial.degree()


def _bolza_constraints(local, points, groups, chart, radius):
    """Use the two disjoint elliptic discs exchanged by x -> i*x.

    The original characteristic-two normalization is rational, so the
    outside of the discs retracts to the bridge midpoint, not an endpoint.
    A full rational orbit is evaluated by its Newton-polygon averages.
    """
    centre, prime, abc, remainder = chart
    K = centre.parent()
    R = PolynomialRing(K, 'i')
    relative = K.extension(R.gen()**2+1, 'i')
    field = relative.absolute_field('b')
    _, embed = field.structure()
    lift = lambda x: embed(relative(x))
    c, other = lift(centre), lift(centre)*lift(relative.gen())
    # F(i*x)=i*F(x) transports the verified elliptic square completion
    # to the second disc, with the same radius and residual genus.
    q = field.ideal([lift(g) for g in prime.gens()]).prime_factors()[0]
    join = _valuation(c-other, q)
    if join != QQ(1)/2 or not join < radius:
        raise ArithmeticError('the two Bolza elliptic discs did not separate')
    half_span = radius-join
    span = 2*half_span
    transformed = (c, q, tuple(lift(a) for a in abc), tuple(lift(a) for a in remainder))
    weights = []
    for x, _ in points:
        a, b = _depth(x, c, q, False), _depth(x, other, q, False)
        weights.append(None if a is None or b is None else
                       half_span+max(join, min(radius, a))-max(join, min(radius, b)))
    model = _compact_type_model(local, groups, weights, span, transformed, False)
    data = graph_constraints(model, 'bolza_two_elliptic_discs')
    data.evidence.update(centre=c, other_centre=other, radius=radius,
                         join_depth=join, weights=weights)
    orbits, averages = {}, [half_span]*len(groups)
    for i, group in enumerate(groups):
        x = points[group[0]][0]
        key = None if x == 'infinity' else tuple(x.minpoly())
        orbits.setdefault(key, []).append(i)
    for key, indices in orbits.items():
        if key is None:
            continue
        polynomial = PolynomialRing(QQ, 'z')(key)
        if len(indices) != polynomial.degree():
            return data  # A partial orbit cannot be replaced by its average.
        mean = (half_span+_mean_clipped_depth(polynomial, c, q, radius, join)
                -_mean_clipped_depth(polynomial, other, q, radius, join))
        for i in indices:
            averages[i] = mean
    endpoint = data.evidence['model'].graph.vertex_boundaries[-1]
    columns = [a/span*endpoint for a in averages]
    data.evidence.update(galois_trace_matrix=matrix(QQ, columns).transpose(),
                         galois_trace_orbits=tuple(tuple(g) for g in orbits.values()),
                         galois_trace_depths=averages)
    return data


def _weight_constraints(method, weights, groups, evidence):
    cols = [weights[group[0]] for group in groups]
    unresolved = tuple(i for i, value in enumerate(cols) if value is None)
    return LocalConstraints(method, matrix(QQ, [[0 if x is None else x for x in cols]]),
                            True, unresolved, evidence)


def _compact_type_model(local, groups, weights, radius, chart, reciprocal):
    centre, prime, (a, b, d), G = chart
    # In the double cover the bridge has half the length of its image.
    # Work with v(2)=1, so adjoining the chart parameters merely rescales
    # the graph uniformly when using a uniformizer of the extension.
    depths = sorted({QQ(0), radius, *(r for r in weights if r is not None)})
    genera = [1]+[0]*(len(depths)-2)+[1]
    edges = [(i, i+1, (s-r)/2) for i, (r, s) in enumerate(zip(depths, depths[1:]))]
    marks = tuple(None if weights[g[0]] is None else (depths.index(weights[g[0]]),)
                  for g in groups)
    parameters = dict(centre=centre, centre_prime=prime, reciprocal=reciprocal,
                      radius=radius, rho_power=(radius.denominator(), 2**radius.numerator()),
                      square_root=a, square_completion=(b, d), remainder=G)
    return StableModel(local.K, local.prime, FibreGraph(genera, edges), marks,
                       charts=(parameters,), length_normalization='v(2)=1')


def _intermediate(local, points, groups, depths, chart):
    centre, prime, (a, b, d), G = chart
    radius = QQ(3)/4
    w = [_valuation(g, prime)+i*radius for i, g in enumerate(G, 3)]
    if not (w[0] < min(w[1:]) and
            1+_valuation(a, prime)/2+min(0, _valuation(b, prime)+radius,
                                         _valuation(d, prime)+2*radius) > w[0]/2):
        raise NotImplementedError('the intermediate cusp chart was not certified')
    # At x=1+z the original Artin–Schreier equation becomes
    # Y^2+Y=z+1 after a shift by z^-2. The two branches over z=0
    # are interchanged. Thus the original component is on the loop,
    # opposite the elliptic bridge's attachment.
    z = local.R.gen()
    f, h = map(local.reduce_poly, (local.f, local.h))
    if h(z+1) != z**2 or f(z+1) != z**4*(z+1)+z**2+1:
        raise NotImplementedError('the two branches on the original component were not certified')
    # All unit points on the intermediate chart have the same
    # projected signature. It differs from the original component's
    # signature because their t-c potential is positive there.
    weights = []
    for point, depth in zip(points, depths):
        label, _, path = local.point(point)
        weights.append(1 if depth == radius else
                       0 if label is not None and len(path) == 1 else None)
    return _weight_constraints('intermediate_chart', weights, groups,
                               dict(radius=radius, cubic_valuation=w[0],
                                    centre=centre.minpoly(), weights=weights))


@lru_cache(maxsize=1)
def _verify_split_annulus():
    R = PolynomialRing(QQ, 's')
    K = NumberField(R.gen()**2-2, 's')
    s = K.gen()
    R = PolynomialRing(K, 'u')
    u = R.gen()
    f, h = [R(list(coeffs)) for coeffs in _SPLIT_ANNULUS]
    F = h*h+4*f
    F = R([F[6-i] for i in range(7)])
    H = 1+u**2+s*(u**3+u)
    if F(s*u)/4-H**2 != -4*(2*u**4+u**2+s*u):
        raise ArithmeticError('the split annulus identity failed')


def _split_annulus(local, points, groups):
    _verify_split_annulus()
    e = local.prime.ramification_index()
    n = QQ(node_thickness(local, local.f, local.h, local.k(1), local.k(0)))/e
    if n != 6:
        raise NotImplementedError('unexpected length of the finite arc')
    weights = []
    for P in points:
        x = P[0]
        if x != 'infinity' and local.valuation(x) >= 0:
            label, _, path = local.point(P)
            if label is None:
                weights.append(None)
            else:
                depth = QQ(len(path)-1)/e
                weights.append(depth*(n-depth))
            continue
        t = local.K(0) if x == 'infinity' else 1/x
        depth = min(QQ(1)/2, _valuation(t, local.prime))
        # On the inner chart u=t/s the only omitted point is u=1.
        if _valuation(t, local.prime) == QQ(1)/2 and _valuation(t*t-2, local.prime) > 1:
            weights.append(None)
        else:
            weights.append(-n*depth-depth**2)
    return _weight_constraints('split_annulus', weights, groups,
                               dict(finite_arc=n, split_arc=QQ(1)/2, weights=weights))


def _two_branch_loop(local, points, groups):
    x = local.R.gen()
    z = x+1
    f, h = map(local.reduce_poly, (local.f, local.h))
    # Shifting the reduced Artin–Schreier ordinate by z^-2+z^-1
    # gives Y^2+Y=1/x, split at x=1 over the geometric residue field.
    if h != x*z**2 or f != x*z**4+x*x+x*x*z**3:
        raise NotImplementedError('the two branches of the loop were not certified')
    reference, _, sections = local.point((local.K(1), local.K(-1)))
    if reference is None or len(sections) < 2:
        raise NotImplementedError('the reference section on the loop was not resolved')
    # The involution fixes the original rational component and swaps
    # its two branches, so it is opposite the sole elliptic component.
    # A section over the singular point retracts away from the antipode.
    # All points sharing the reference component have the same weight.
    weights = []
    for P in points:
        component, _, path = local.point(P)
        weights.append(0 if component == 'main:0' and len(path) == 1 else
                       1 if component == reference else None)
    return _weight_constraints('two_branch_loop', weights, groups,
                               dict(reference_component=reference, weights=weights))


def dyadic_constraints(local, points, groups, stable_type):
    try:
        key = tuple(tuple(QQ(x) for x in g) for g in (local.f, local.h))
    except TypeError:
        raise NotImplementedError('no verified wild model for this equation') from None
    if key == _SPLIT_ANNULUS and stable_type == 'VI':
        return _split_annulus(local, points, groups)
    if key == _TWO_BRANCH_LOOP and stable_type == 'II':
        return _two_branch_loop(local, points, groups)
    if key not in _ELLIPTIC_CHARTS or stable_type != ('VI' if key == _INTERMEDIATE_CHART else 'V'):
        raise NotImplementedError('no verified wild model for this equation')
    reciprocal, depth, radius = _ELLIPTIC_CHARTS[key]
    chart = elliptic_chart(key, reciprocal, depth, radius)
    if key == _BOLZA:
        return _bolza_constraints(local, points, groups, chart, radius)
    centre, prime, _, _ = chart
    depths = [_depth(x, centre, prime, reciprocal) for x, _ in points]
    if key == _INTERMEDIATE_CHART:
        return _intermediate(local, points, groups, depths, chart)
    weights = [None if d is None else max(0, min(radius, d)) for d in depths]
    model = _compact_type_model(local, groups, weights, radius, chart, reciprocal)
    data = graph_constraints(model, 'compact_type')
    data.evidence.update(centre=centre.minpoly(), radius=radius, weights=weights)
    return data
