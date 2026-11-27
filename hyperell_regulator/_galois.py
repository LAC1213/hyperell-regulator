"""Galois traces of divisor directions on the rational standard motive."""

from sage.all import QQ, QQbar, matrix, vector


def galois_traces(abscissae, divisors):
    """Return normalized full Galois traces in the supplied abscissa coordinates.

    Rows have rational coefficients after identifying symmetric ordinates.
    Averaging each complete orbit is the exact Reynolds projection over QQ.
    A trace requiring a missing conjugate cannot be represented here.
    """
    xs = tuple('infinity' if x == 'infinity' else QQbar(x) for x in abscissae)
    if len(set(xs)) != len(xs):
        raise ValueError('Galois traces require distinct abscissae')
    rows = [tuple(map(QQ, row)) for row in divisors]
    if any(len(row) != len(xs) for row in rows):
        raise ValueError('one coefficient per abscissa is required')
    orbits = {}
    for i, x in enumerate(xs):
        key = None if x == 'infinity' else tuple(x.minpoly())
        orbits.setdefault(key, []).append(i)
    output = matrix(QQ, len(rows), len(xs))
    for key, indices in orbits.items():
        for j, row in enumerate(rows):
            total = sum(row[i] for i in indices)
            if key is not None and len(indices) != len(key)-1 and total:
                raise ValueError('the trace requires conjugates outside the supplied support')
            for i in indices:
                output[j, i] = total/len(indices)
    return output


def galois_trace_basis(abscissae):
    """Span degree-zero normalized traces in the supplied abscissa order.

    The standard projection identifies the two ordinates over each x, so
    the orbits are the roots of the minimal polynomials of x over QQ, with
    infinity a singleton. Differences of orbit sums divided by orbit size
    span the degree-zero traces. Incomplete orbits contribute nothing.
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
