"""Period paths and regularised differentials for both precision backends."""

from ._paths import semilunar_path


def _part(z, name):
    value = getattr(z, name)
    return value() if callable(value) else value


def period_paths(roots, g, root_tol=0):
    """First ``2*g`` branch-to-branch paths, starting at the largest modulus.

    Real and imaginary parts break modulus ties. Root endpoints retain their
    input precision; only ordinary detour corners come from float64 geometry.
    All roots, including the unused sextic endpoint, enter path validation.
    """
    candidates = list(range(len(roots)))
    for key in (abs, lambda z: _part(z, 'real'), lambda z: _part(z, 'imag')):
        largest = max(key(roots[i]) for i in candidates)
        tolerance = root_tol * max(abs(roots[i]) for i in candidates)
        candidates = [i for i in candidates if largest - key(roots[i]) <= tolerance]
    start = candidates[0]
    nodes, labels, _ = semilunar_path(roots, start, include_ray=False)
    paths, current = [], []
    for node, label in zip(nodes, labels):
        current.append((roots[label] if label is not None else node, label))
        if label is not None and len(current) > 1:
            paths.append(current)
            if len(paths) == 2 * g:
                break
            current = [current[-1]]
    return paths


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
                         sqrt, sin, cos, pi, sign_tol, previous_y=None):
    """Regularise endpoints and continue the square root through corners.

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

    if previous_y is not None:
        ya = scaled_y(a)
        if abs(ya - previous_y) > abs(ya + previous_y):
            ymid = -ymid

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
