"""Galois traces of divisor directions on the rational standard motive."""

from sage.all import QQ, QQbar, matrix, vector


def galois_trace_basis(abscissae):
    """Span degree-zero normalized traces in the supplied abscissa order.

    The standard projection identifies the two ordinates over each x.
    Its Galois orbits are therefore the roots of the minimal polynomials
    of x over QQ, with infinity a rational singleton. A normalized orbit
    sum is the trace of a point divided by its orbit size. Differences
    of these sums span all degree-zero traces over QQ.

    An incomplete orbit cannot support a nonzero invariant coefficient,
    so it contributes no generator. Arithmetic is exact, including for
    nonreal conjugates. Callers clear denominators after lifting back to
    point divisors.
    """
    xs = tuple('infinity' if x == 'infinity' else QQbar(x) for x in abscissae)
    if len(set(xs)) != len(xs):
        raise ValueError('Galois traces require distinct abscissae')
    orbits = {}
    for i, x in enumerate(xs):
        key = None if x == 'infinity' else tuple(x.minpoly())
        orbits.setdefault(key, []).append(i)
    averages = [vector(QQ, len(xs), {i: QQ(1)/len(indices) for i in indices})
                for key, indices in orbits.items()
                if key is None or len(indices) == len(key)-1]
    rows = [v-averages[0] for v in averages[1:]]
    return matrix(QQ, len(rows), len(xs), rows)
