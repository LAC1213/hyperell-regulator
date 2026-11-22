r"""Validated integration on algebraic curves using Arb balls.

Plane integrals default to a certified embedded tree with one infinity ray.
An independent implementation using parallel slits is also available. Both
determine sheet gluing and logarithmic jumps without sampled winding tests.
"""

import math
from contextlib import nullcontext

from sage.all import QQ, QQbar, PolynomialRing, ComplexBallField, RealBallField, NaN

from .polynomials import polynomial
from .diagnostics import _active, checkpoint


class CertificationError(ArithmeticError):
    """The requested rigorous enclosure could not be established."""


def digits_value(digits):
    digits = 15 if digits is None else digits
    if int(digits) != digits or digits < 5:
        raise ValueError('precision must be an integer of at least 5 decimal digits')
    return int(digits)


def rational(value):
    """Coerce point inputs exactly, avoiding QQ's float reconstruction."""
    if hasattr(value, 'exact_rational'):
        return QQ(value.exact_rational())
    if hasattr(value, 'as_integer_ratio'):
        numerator, denominator = value.as_integer_ratio()
        return QQ(numerator)/denominator
    return QQ(value)


def algebraic_center(value):
    """Keep exact algebraic centres and preserve exact rational coercion."""
    try:
        return rational(value)
    except (TypeError, ValueError):
        return QQbar(value)


def _polynomial(value):
    if hasattr(value, 'list'):
        return PolynomialRing(QQ, 'x')([rational(c) for c in value.list()])
    if isinstance(value, (list, tuple)):
        return PolynomialRing(QQ, 'x')([rational(c) for c in reversed(value)])
    return polynomial(value)


def fields(digits, guard=40):
    bits = math.ceil((digits_value(digits) + guard) * math.log2(10))
    return ComplexBallField(bits), RealBallField(bits)


def narrow(value, digits):
    """Check the actual absolute radius, never just a requested tolerance."""
    if hasattr(value, 'imag'):
        finite = value.real().is_finite() and value.imag().is_finite()
    else:
        finite = value.is_finite()
    R = RealBallField(max(64, value.parent().precision()))
    return finite and R(value.rad()) < R(10)**(-digits)


def principal(z, analytic, operation):
    """Report failure of holomorphy across the principal negative-axis cut."""
    if analytic and z.imag().contains_zero() and not z.real() > 0:
        return z.parent()(NaN)
    return getattr(z, operation)()


def _model(F):
    F = _polynomial(F)
    if F is None or F.degree() < 5 or F.gcd(F.derivative()).degree() != 0:
        raise ValueError('F must be an exact squarefree polynomial of degree at least five')
    return F


def _parallel_geometry(roots, centers):
    nodes = list(roots)
    for c in centers:
        if QQbar(c) not in nodes:
            nodes.append(QQbar(c))
    # There are at most n(n-1)/2 forbidden slopes. AA comparisons are exact.
    for k in range(1, len(nodes)**2 + 2):
        heights = [z.imag() - k*z.real() for z in nodes]
        if all(heights[i] != heights[j] for i in range(len(nodes)) for j in range(i)):
            return nodes, QQbar(1) + QQbar.gen()*k, sorted(range(len(nodes)), key=lambda j: heights[j])
    raise CertificationError('could not separate parallel slit vertices')


def _ray_integrals(F, roots, vertex, direction, centers, digits):
    r"""Enclose x^j dx/y and complex-log weighted integrals on one slit.

    x = z + D*t^2/(1-t)^2 regularizes both infinity and a branch endpoint.
    The reference square root is the upper-bank value of
    sqrt(lc*(-D)^d) * product(sqrt(-(x-r)/D)). Logarithms use lower banks.
    """
    C, R = fields(digits)
    g, d = (F.degree()-1)//2, F.degree()
    z, D = C(vertex), C(direction)
    K = C(QQbar(F.leading_coefficient()*(-direction)**d).sqrt())
    branch = vertex in roots
    factors = []
    for root in roots:
        if root == vertex:
            continue
        a = -(vertex-root)/direction
        if a.imag() == 0:
            raise CertificationError('parallel slit separation failed')
        phase = C.gen(0) if a.imag() > 0 else -C.gen(0)
        factors.append((C(a), phase))
    log_constant = C(-direction).log()
    logs = []
    for center in centers:
        own = vertex == center
        a = -(vertex-center)/direction
        phase = 1 if own or a.imag() > 0 else -1
        logs.append((own, C(a), log_constant + phase*C.gen(0)*C.pi()))

    def differential(t, analytic, j):
        v = 1-t
        X = z*v*v + D*t*t
        y = K * (-C.gen(0) if branch else 1)
        for a, phase in factors:
            y *= phase*principal(t*t-a*v*v, analytic, 'sqrt')
        return 2*D*(1 if branch else t)*v**(d-3-2*j)*X**j/y

    def log_value(t, analytic, record, omit=None):
        own, a, constant = record
        value = constant
        if own:
            if omit != 'left':
                value += 2*principal(t, analytic, 'log')
        else:
            value += principal(t*t-a*(1-t)**2, analytic, 'log')
        if omit != 'right':
            value -= 2*principal(1-t, analytic, 'log')
        return value

    eps = QQ(1)/10**(digits+10)
    tol = R(10)**(-digits-5)
    left = C(R(eps/2).add_error(R(eps/2)))
    right = 1-left

    def integrate(j, record=None):
        error = R(0)
        for t, side in ((left, 'left'), (right, 'right')):
            amplitude = abs(differential(t, False, j))
            if record is None:
                bound = R(1)
            else:
                singular = side == 'right' or record[0]
                bound = abs(log_value(t, False, record, side if singular else None))
                if singular:
                    bound += 2*(1-R(eps).log())
            error += R(eps)*amplitude*bound
        def evaluate(t, analytic):
            checkpoint('quadrature_nodes', phase='certified ray integration')
            value = differential(t, analytic, j)
            return value if record is None else value*log_value(t, analytic, record)
        result = C.integral(evaluate, eps, 1-eps, abs_tol=tol, rel_tol=tol).add_error(error)
        if not narrow(result, digits):
            raise CertificationError('parallel-ray integral did not reach the requested radius')
        return result

    periods = [integrate(j) for j in range(g)]
    weighted = [[integrate(j, record) for j in range(g)] for record in logs]
    return periods, weighted


def _assemble(C, R, g, nodes, roots, centers, order, edge_periods, weighted, log_index=None):
    """The Stokes bilinear identity with exact parallel-slit bank data."""
    n = len(nodes)
    word = [(j, side) for j in order for side in (-1, 1)]
    upper = [[0, 1] for _ in nodes]
    lower = [[1, 0] if z in roots else [0, 1] for z in nodes]
    walks = [[(j, (upper if side > 0 else lower)[j][sheet], side)
              for j, side in word] for sheet in (0, 1)]
    u = [[[row[a], -row[a]] for row in edge_periods] for a in range(g)]
    if log_index is None:
        ell = u
        jumps = [C(0)]*n
    else:
        ell = [[[row[log_index][a], -row[log_index][a]] for row in weighted] for a in range(g)]
        jumps = [-2*C.gen(0)*C.pi() if z == centers[log_index] else C(0) for z in nodes]
    return _assemble_walks(C, R, g, walks, u, ell, jumps)


def _assemble_walks(C, R, g, walks, u, ell, jumps):
    """Stokes assembly for a certified boundary word and bank data."""
    n = len(jumps)
    forward, backward = {}, {}
    for k, walk in enumerate(walks):
        for j, (i, sheet, side) in enumerate(walk):
            (forward if side > 0 else backward)[i, sheet] = k, j
    H = [[C(0) for _ in range(g)] for _ in range(g)]
    for a in range(g):
        primitives = []
        for walk in walks:
            row = [C(0)]
            for i, sheet, side in walk:
                value = ell[a][i][sheet]
                row.append(row[-1] + (value+jumps[i]*u[a][i][sheet] if side > 0 else -value))
            if not row[-1].real().contains_zero() or not row[-1].imag().contains_zero():
                raise CertificationError('certified slit boundary does not close')
            primitives.append(row)
        for i in range(n):
            for sheet in (0, 1):
                k, j = forward[i, sheet]
                l, jj = backward[i, sheet]
                constant = primitives[k][j]-primitives[l][jj]+ell[a][i][sheet]
                for b in range(g):
                    H[a][b] += constant*u[b][i][sheet].conjugate()
    result = []
    for a in range(g):
        for b in range(a, g):
            moment = sum((jumps[i].imag()/2*(u[a][i][s]*u[b][i][s].conjugate()).real()
                          for i in range(n) for s in (0, 1)), R(0))
            result.append(-((H[a][b].imag()+H[b][a].imag())/2+moment)/4)
    return result


def compute_parallel(F, centers=(), prec=15, *, details=None):
    """Independent certified parallel-slit regulator backend."""
    digits = digits_value(prec)
    F = _model(F)
    centers = list(dict.fromkeys(algebraic_center(c) for c in centers))
    roots = F.roots(QQbar, multiplicities=False)
    nodes, direction, order = _parallel_geometry(roots, centers)
    for guard in (10, 25, 50):
        target = digits+guard
        C, R = fields(target)
        edge_periods, weighted = [], []
        for i, vertex in enumerate(nodes):
            monitor = _active.get()
            with monitor.edge('certified ray %d' % i) if monitor else nullcontext():
                u, ell = _ray_integrals(F, roots, vertex, direction, centers, target)
            edge_periods.append(u)
            weighted.append(ell)
        g = (F.degree()-1)//2
        per = _assemble(C, R, g, nodes, roots, centers, order, edge_periods, weighted)
        reg = [_assemble(C, R, g, nodes, roots, centers, order, edge_periods, weighted, k)
               for k in range(len(centers))]
        if all(narrow(v, digits) for v in per + [v for row in reg for v in row]):
            if details is not None:
                details.update(certified=True, cut_kind='parallel rays',
                               exact_vertices=nodes, exact_direction=direction,
                               transverse_order=order, working_bits=C.precision())
            return per, reg
    raise CertificationError('assembled regulator did not reach the requested radius')


def compute(F, centers=(), prec=15, *, details=None, cut='tree'):
    """Certified plane integrals; the embedded tree is the default cut."""
    if cut == 'parallel':
        return compute_parallel(F, centers, prec, details=details)
    if cut != 'tree':
        raise ValueError("certified_cut must be 'tree' or 'parallel'")
    from .tree_balls import compute_tree
    return compute_tree(F, centers, prec, details=details)


def edge_integrals(F, a, b, *, centers=(), powers=(0, 1), digits=30):
    r"""Enclose integrals of ``x^j dx/y`` and ``log(x-c) x^j dx/y``.

    ``a`` and ``b`` are exact algebraic endpoints, each either a branch point
    or a regular point. The sheet and logarithm branches are fixed at the
    midpoint and continued along the segment. Returned logarithms are complex
    logarithms, as used in the Stokes code, not ``log(abs(x-c))``.

    A failed analyticity check or an insufficient final absolute radius raises
    ``ArithmeticError``. In particular, the integrator's requested tolerance
    is never interpreted as its achieved error bound.
    """
    if int(digits) != digits or digits < 5:
        raise ValueError("digits must be an integer of at least 5")
    F = _polynomial(F)
    if F is None or F.degree() < 2 or F.gcd(F.derivative()).degree() != 0:
        raise ValueError("F must be an exact squarefree polynomial of degree >= 2")
    powers = tuple(powers)
    if not powers or any(int(j) != j or j < 0 for j in powers):
        raise ValueError("powers must be nonnegative integers")
    bits = math.ceil((int(digits) + 40) * math.log2(10))
    C, R = ComplexBallField(bits), RealBallField(bits)
    aa, bb = QQbar(a), QQbar(b)
    if aa == bb:
        raise ValueError("edge endpoints must differ")
    roots = F.roots(QQbar, multiplicities=False)
    others = [r for r in roots if r != aa and r != bb]
    for r in others:
        position = (r - aa) / (bb - aa)
        if position.imag() == 0 and 0 <= position.real() <= 1:
            raise ValueError("edge passes through another branch point")
    centers = tuple(algebraic_center(c) for c in centers)
    for c in centers:
        position = (QQbar(c) - aa) / (bb - aa)
        if position.imag() == 0 and 0 < position.real() < 1:
            raise ValueError("split the edge at an interior logarithmic center")

    a, b = C(aa), C(bb)
    delta, mid = b - a, (a + b) / 2
    branch_a, branch_b = aa in roots, bb in roots
    others = [C(r) for r in others]
    qmid = C(F.leading_coefficient())
    if branch_a:
        qmid *= delta
    if branch_b:
        qmid *= -delta
    for r in others:
        qmid *= mid - r
    # Avoid the principal sqrt discontinuity at a negative real seed.
    if qmid.real() < 0:
        if qmid.imag() < 0:
            side = -1
        elif qmid.imag() > 0:
            side = 1
        else:
            exact_q = 2**(int(branch_a)+int(branch_b))*F((aa+bb)/2)
            side = -1 if exact_q.imag() < 0 else 1
        ymid = side*C.gen(0)*(-qmid).sqrt()
    else:
        ymid = qmid.sqrt()
    root_factors = [(r, 1/(mid-r)) for r in others]
    half_pi, log_two = C.pi()/2, C(2).log()
    seed_logs = {}
    log_positions = {}
    centre_balls = {c: C(c) for c in centers}
    for center in centers:
        z = mid - centre_balls[center]
        exact_z = (aa + bb) / 2 - center
        if exact_z.real() < 0:
            # Select a single logarithm branch even if root enclosures
            # straddle the negative real axis at the midpoint.
            side = -1 if exact_z.imag() < 0 else 1
            seed_logs[center] = (-z).log() + side * C.gen(0) * C.pi()
        else:
            seed_logs[center] = z.log()
        log_positions[center] = 0 if QQbar(center) == aa else 1 if QQbar(center) == bb else 2

    def parameters(t):
        angle = half_pi*t
        return angle.sin(), angle.cos()

    def principal(z, analytic, operation):
        if analytic and z.imag().contains_zero() and not z.real() > 0:
            return C(NaN)
        return getattr(z, operation)()

    def differential(t, analytic, j, sc=None):
        s, c = parameters(t) if sc is None else sc
        x = a + delta * s**2
        y = ymid
        for r, inverse in root_factors:
            y *= principal((x-r)*inverse, analytic, "sqrt")
        rate = 2*half_pi*delta/y
        if not branch_a:
            rate *= s
        if not branch_b:
            rate *= c
        return rate * x**j

    def logarithm(t, analytic, center, sc=None):
        s, c = parameters(t) if sc is None else sc
        lam = seed_logs[center]
        if log_positions[center] == 0:
            return lam + log_two + 2 * principal(s, analytic, "log")
        if log_positions[center] == 1:
            return lam + log_two + 2 * principal(c, analytic, "log")
        x = a + delta * s**2
        c = centre_balls[center]
        return lam + principal((x - c) / (mid - c), analytic, "log")

    eps = QQ(1) / 10 ** (int(digits) + 10)
    tolerance = R(10) ** (-int(digits) - 5)

    def tail(j, at_right, center):
        # This ball contains [0, eps], or [1-eps, 1]. The differential
        # extends analytically to both endpoints after root cancellation.
        t = C(R(eps / 2).add_error(R(eps / 2)))
        if at_right:
            t = 1 - t
        amplitude = abs(differential(t, False, j))
        if center is None:
            return R(eps) * amplitude
        endpoint = bb if at_right else aa
        if QQbar(center) == endpoint:
            # sin(pi*u/2) >= u for 0 <= u <= 1. Thus
            # integral_0^eps |log(sin(pi*u/2))| du
            # <= eps*(1-log(eps)); use the same bound for cos at t=1.
            constant = abs(seed_logs[center] + C(2).log())
            return R(eps) * amplitude * (constant + 2 * (1 - R(eps).log()))
        return R(eps) * amplitude * abs(logarithm(t, False, center))

    def integrate(j, center=None):
        def evaluate(t, analytic):
            checkpoint('quadrature_nodes', phase='certified period integration')
            sc = parameters(t)
            value = differential(t, analytic, j, sc)
            return value if center is None else value * logarithm(t, analytic, center, sc)

        error = tail(j, False, center) + tail(j, True, center)
        value = C.integral(evaluate, eps, 1 - eps,
                           abs_tol=tolerance, rel_tol=tolerance)
        value = value.add_error(error)
        if (not value.real().is_finite() or not value.imag().is_finite()
                or not R(value.rad()) < R(10) ** (-int(digits))):
            raise ArithmeticError("edge enclosure did not reach the requested absolute radius")
        return value, error

    periods = [integrate(j) for j in powers]
    logs = [[integrate(j, center) for j in powers] for center in centers]
    return dict(periods=[v for v, _ in periods],
                logarithmic=[[v for v, _ in row] for row in logs],
                reference_y=ymid/C(2).sqrt()**(int(branch_a)+int(branch_b)),
                reference_logs=[seed_logs[c] for c in centers],
                endpoint_bounds=dict(periods=[e for _, e in periods],
                                     logarithmic=[[e for _, e in row] for row in logs]),
                bits=bits, digits=int(digits),
                scope="finite edge only; locally chosen sheet and logarithm branches")


def _period_edges(roots, g, digits):
    """Propose the existing tree basis, then certify its embedded geometry.

    Minimality is not needed for the integral-cycle basis theorem. Root
    approximations choose the edges; all geometric claims are checked on
    enclosures of the exact roots before they are used.
    """
    from sage.all import ComplexField
    from ._periods import period_paths
    P = ComplexField(math.ceil((digits+40)*math.log2(10)))
    approx = [P(z) for z in roots]
    try:
        paths = period_paths(approx, g, root_tol=P(10).real()**(-digits-20))
        edges = [(a[1], b[1]) for a, b in paths]
    except ValueError:
        # Arbitrarily close roots need not be distinguishable in float64.
        edges = _precision_tree(approx, g)
    parent = list(range(len(roots)))
    def find(i):
        while parent[i] != i:
            i = parent[i]
        return i
    for a, b in edges:
        if a == b or find(a) == find(b):
            raise CertificationError('cycle paths are not a tree')
        parent[find(a)] = find(b)
    vertices = {i for edge in edges for i in edge}
    if len(edges) != 2*g or len(vertices) != 2*g+1 or len({find(i) for i in vertices}) != 1:
        raise CertificationError('cycle tree does not have the required rank')
    C, _ = fields(digits)
    z = [C(r) for r in roots]
    def separated(a, b, c, d):
        for part in (lambda w: w.real(), lambda w: w.imag()):
            ab = [part(z[a]), part(z[b])]
            cd = [part(z[c]), part(z[d])]
            if max(v.upper() for v in ab) < min(v.lower() for v in cd):
                return True
            if max(v.upper() for v in cd) < min(v.lower() for v in ab):
                return True
        def cross(i, j, k):
            return ((z[j]-z[i]).conjugate()*(z[k]-z[i])).imag()
        x, y = cross(a, b, c), cross(a, b, d)
        if (x > 0 and y > 0) or (x < 0 and y < 0):
            return True
        x, y = cross(c, d, a), cross(c, d, b)
        return (x > 0 and y > 0) or (x < 0 and y < 0)
    for i, (a, b) in enumerate(edges):
        for c, d in edges[:i]:
            if len({a, b, c, d}) == 4 and not separated(a, b, c, d):
                raise CertificationError('could not certify that cycle-tree edges are disjoint')
    return edges


def _precision_tree(roots, g):
    """High-precision proposal when float64 cannot distinguish vertices."""
    parent = list(range(len(roots)))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i
    pairs = sorted((abs(roots[i]-roots[j]), i, j)
                   for i in range(len(roots)) for j in range(i))
    edges = []
    for _, i, j in pairs:
        a, b = find(i), find(j)
        if a != b:
            parent[a] = b
            edges.append((i, j))
    if len(roots) == 2*g+2:
        degrees = [0]*len(roots)
        for a, b in edges:
            degrees[a] += 1
            degrees[b] += 1
        leaves = [i for i, d in enumerate(degrees) if d == 1]
        omitted = max(leaves, key=lambda i: (abs(roots[i]), roots[i].real(), roots[i].imag()))
        edges = [(a, b) for a, b in edges if omitted not in (a, b)]
    return edges


def period_matrix(F, prec=15):
    """Full cycle periods on the exact tree basis as complex balls."""
    digits = digits_value(prec)
    F = _model(F)
    g = (F.degree()-1)//2
    roots = F.roots(QQbar, multiplicities=False)
    edges = _period_edges(roots, g, digits)
    for guard in (8, 20, 40):
        result = []
        for a, b in edges:
            monitor = _active.get()
            with monitor.edge('certified period %d-%d' % (a, b)) if monitor else nullcontext():
                result.append([2*v for v in edge_integrals(F, roots[a], roots[b],
                                              powers=range(g), digits=digits+guard)['periods']])
        if all(narrow(z, digits) for row in result for z in row):
            return result
    raise CertificationError('cycle periods did not reach the requested radius')


def c_invariants(F, prec=15):
    """Signed period minors; every division uses a certified nonzero pivot."""
    digits = digits_value(prec)
    for guard in (8, 20, 40):
        periods = period_matrix(F, digits+guard)
        R = periods[0][0].real().parent()
        def minor(part):
            rows = [[part(z) for z in row] for row in periods]
            det = R(1)
            m = len(rows[0])
            for k in range(m):
                candidates = [i for i in range(k, len(rows)) if not rows[i][k].contains_zero()]
                if not candidates:
                    raise CertificationError('period minor has no certified nonzero pivot')
                largest = max(abs(rows[i][k].mid()) for i in candidates)
                tol = largest.parent()(10)**(-digits-3)*largest
                p = next(i for i in candidates if largest-abs(rows[i][k].mid()) <= tol)
                if p != k:
                    rows[k], rows[p] = rows[p], rows[k]
                    det = -det
                det *= rows[k][k]
                for i in range(k+1, len(rows)):
                    factor = rows[i][k]/rows[k][k]
                    for j in range(k+1, m):
                        rows[i][j] -= factor*rows[k][j]
            return det
        values = (minor(lambda z: z.real()), minor(lambda z: z.imag()))
        if all(narrow(v, digits) for v in values):
            return values
    raise CertificationError('period minors did not reach the requested radius')

