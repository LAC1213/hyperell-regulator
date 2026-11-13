"""MST cycle paths and regularised square-root differentials."""

from ._tree import minimum_spanning_tree


def _part(z, name):
    value = getattr(z, name)
    return value() if callable(value) else value


def period_paths(roots, g, root_tol=0):
    """An integral cycle basis from doubled branch-to-branch MST edges.

    With 2g+1 finite branch values, all 2g edges give a basis. With 2g+2,
    the lifted edges have one relation. A leaf edge has coefficient +/-1
    in that relation, so omitting it gives an integral basis. Omit the
    leaf of largest modulus, with real/imaginary parts breaking ties.
    Exact input roots are retained; float64 chooses only the tree topology.
    """
    if len(roots) not in (2*g+1, 2*g+2) or g < 1:
        raise ValueError('branch values do not match the genus')
    # Input root solvers need not return the same ordering. Give tied MST
    # edges a geometric order, independent of the solver's list order.
    order = sorted(range(len(roots)), key=lambda i: (float(_part(roots[i], 'real')),
                                                    float(_part(roots[i], 'imag'))))
    edges = [(order[a], order[b]) for a,b in minimum_spanning_tree([roots[i] for i in order])]
    if len(roots) == 2*g+2:
        degree = [0]*len(roots)
        for a, b in edges:
            degree[a] += 1
            degree[b] += 1
        candidates = [i for i, d in enumerate(degree) if d == 1]
        for key in (abs, lambda z: _part(z, 'real'), lambda z: _part(z, 'imag')):
            largest = max(key(roots[i]) for i in candidates)
            tolerance = root_tol*max(abs(roots[i]) for i in candidates)
            candidates = [i for i in candidates if largest-key(roots[i]) <= tolerance]
        omitted = candidates[0]
        edges = [(a,b) for a,b in edges if omitted not in (a,b)]
    return [[(roots[a],a), (roots[b],b)] for a,b in edges]


def pivoted_minor(rows, tol):
    """Signed maximal minor; working-precision pivot ties use the first row.

    Conjugate cycles often give equal pivot magnitudes. Resolving these ties
    consistently prevents rounding from changing signs across backends.
    """
    a = [list(row) for row in rows]
    n, m = len(a), len(a[0])
    det = a[0][0] * 0 + 1
    for k in range(m):
        largest = max(abs(a[i][k]) for i in range(k, n))
        if not largest:
            return det * 0
        p = next(i for i in range(k, n) if largest - abs(a[i][k]) <= tol * largest)
        if p != k:
            a[k], a[p] = a[p], a[k]
            det = -det
        det *= a[k][k]
        for i in range(k + 1, n):
            factor = a[i][k] / a[k][k]
            for j in range(k + 1, m):
                a[i][j] -= factor * a[k][j]
    return det


def segment_differential(roots, leading, a, b, ia, ib, g, *,
                         sqrt, sin, cos, pi, sign_tol):
    """Regularise endpoints and continue the square root along an edge.

    Set ``x=a+(b-a)*sin(pi*t/2)^2``. Endpoint root factors cancel
    analytically. All other factors are continued relative to the midpoint;
    their principal square roots are continuous on a validated segment.
    """
    delta, mid = b - a, (a + b) / 2
    others = [r for i, r in enumerate(roots) if i not in (ia, ib)]
    qmid = leading
    if ia is not None:
        qmid *= delta
    if ib is not None:
        qmid *= -delta
    for r in others:
        qmid *= mid - r
    ymid = sqrt(qmid)
    if abs(_part(ymid, 'real')) <= sign_tol * abs(ymid) and _part(ymid, 'imag') < 0:
        ymid = -ymid

    def scaled_y(x):
        value = ymid
        for r in others:
            value *= sqrt((x - r) / (mid - r))
        return value

    def evaluate(t):
        s, c = sin(pi * t / 2), cos(pi * t / 2)
        x = a + delta * s**2
        rate = pi * delta / scaled_y(x)
        if ia is None:
            rate *= s
        if ib is None:
            rate *= c
        return [rate * x**i for i in range(g)]

    return evaluate, scaled_y(b) if ib is None else None
