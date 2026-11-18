"""Exact orders of divisor classes by principality tests."""

def _field(F):
    r"""`(\QQ(x), \QQ(x)[y]/(y^2 - F))` for the curve `y^2 = F(x)`."""
    from sage.all import FunctionField, PolynomialRing, QQ

    from .polynomials import polynomial

    f = polynomial(F)
    if f is None:
        raise ValueError("F must be a polynomial over QQ")
    K = FunctionField(QQ, "x")
    R = PolynomialRing(K, "Y")
    Y = R.gen()
    return K, K.extension(Y ** 2 - K(f), "y")


def _places(K, L, point):
    r"""Places of ``L`` over a finite point or infinity.

    ``("infinity", eta)`` selects `\lim y/x^{g+1}=eta`;
    bare ``"infinity"`` returns all infinite places.
    """
    from sage.all import QQ

    x = K.gen()
    y = L.gen()
    if point == "infinity":
        return list(L.places_infinite())
    if point[0] == "infinity":
        eta = QQ(point[1])
        f = -L.polynomial()[0]
        degree = f.numerator().degree()
        value = y/x**((degree+1)//2)
        places = [pl for pl in L.places_infinite()
                  if pl.degree() == 1 and (value-eta).valuation(pl) > 0]
        if not places:
            raise ValueError("%r is not a rational point at infinity" % (point,))
        return places
    x0, y0 = QQ(point[0]), QQ(point[1])
    below = K.maximal_order().ideal(x - x0).place()
    for pl in L.places_above(below):
        if (y - y0).valuation(pl) > 0:
            return [pl]
    raise ValueError("%r is not a point of the curve" % (point,))


def _place(K, L, point):
    r"""The unique place at ``point``; raise if there is more than one."""
    pls = _places(K, L, point)
    if len(pls) != 1 or pls[0].degree() != 1:
        raise ValueError("point must specify a single rational place of this model")
    return pls[0]


def _order_from_multiple(divisor, multiple):
    """The exact divisor-class order if it divides ``multiple``; otherwise ``None``."""
    from sage.all import ZZ

    order = ZZ(multiple)
    if order <= 0:
        raise ValueError("multiple must be a positive integer")
    if not (order*divisor).basis_function_space():
        return None
    for prime in order.prime_divisors():
        while order % prime == 0:
            smaller = order//prime
            if not (smaller*divisor).basis_function_space():
                break
            order = smaller
    return order


def torsion_order(F, P, Q, bound=24, *, multiple=None):
    r"""The least `N\le` ``bound`` with `N(P-Q)` principal, or ``None``.

    With ``multiple``, find the exact order among its divisors, ignoring ``bound``.
    ``("infinity", eta)`` selects `\lim y/x^{g+1}=eta`; bare ``"infinity"``
    returns the smallest order over rational infinite places.

    EXAMPLES::

        sage: from hyperell_regulator.regulator import torsion_order
        sage: from hyperell_regulator.regulator.polynomials import ring
        sage: x = ring().gen()
        sage: torsion_order(x^5 + x^2 + 2*x + 1, (0, -1), "infinity")
        5
        sage: F = (1 + x + x^3)^2 + 4*(x + 2*x^2 - x^4)
        sage: torsion_order(F, (0, 1), ("infinity", -1), bound=1, multiple=21)
        21
    """
    K, L = _field(F)
    ps = [p for p in _places(K, L, P) if p.degree() == 1]
    qs = [q for q in _places(K, L, Q) if q.degree() == 1]
    if not ps or not qs:
        raise ValueError("P and Q must be rational points")
    differences = [q.divisor()-p.divisor() for p in ps for q in qs]
    if multiple is not None:
        orders = [_order_from_multiple(D, multiple) for D in differences]
        return min((n for n in orders if n is not None), default=None)
    for N in range(1, bound + 1):
        if any((N*D).basis_function_space() for D in differences):
            return N
    return None
