r"""
Producing a function with a prescribed divisor `N P - N Q`.

`\operatorname{div}(\varphi) = NP - NQ` says exactly that
`\varphi \in L(NQ - NP)`, a space of degree zero: it is nonzero precisely when
`N(P-Q)` is principal, and then one dimensional, spanned by the function
wanted.  Sage's function fields compute that space, and the result is read
back as the pair `(A, B)` with `\varphi = A(x) + B(x)\,y` that the rest of the
package takes -- two rational functions of `x`, as they come out of the
function field.
"""

__all__ = ["function_with_divisor", "torsion_order"]


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
    r"""
    The places of `L` at a finite point or at infinity.

    ``("infinity", eta)`` selects the branch with
    `y/x^{g+1} \to eta`. Bare ``"infinity"`` returns every infinite place;
    :func:`function_with_divisor` tries them in turn.
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
    r"""The place at ``point``, which must be a single one."""
    pls = _places(K, L, point)
    if len(pls) != 1 or pls[0].degree() != 1:
        raise ValueError("point must specify a single rational place of this model")
    return pls[0]


def _in_x(c):
    r"""
    An element of the rational function field as an element of `\QQ(x)`.

    The two live in different parents -- one is a function field, the other
    the fraction field of a polynomial ring -- and it is the second that the
    rest of the package works in.
    """
    from .polynomials import ring

    R = ring()
    K = R.fraction_field()
    return K(R(c.numerator().list())) / K(R(c.denominator().list()))


def function_with_divisor(F, P, Q, N):
    r"""
    `\varphi = (A, B)` with `\operatorname{div}(\varphi) = N P - N Q` on
    `y^2 = F(x)`, as a pair of elements of `\QQ(x)`.

    `F` is a polynomial over `\QQ`; `P` and `Q` are finite pairs `(x_0, y_0)`,
    ``("infinity", eta)`` with `eta = \lim y/x^{g+1}`, or bare ``"infinity"``
    to try every infinite place. Raises :exc:`ValueError` when `N(P-Q)` is not
    principal, which is the same as `L(NQ - NP)` being zero.

    EXAMPLES::

        sage: from hyperell_regulator.regulator import function_with_divisor
        sage: from hyperell_regulator.regulator.polynomials import ring
        sage: x = ring().gen()
        sage: function_with_divisor(x^5 + x^2 + 2*x + 1, (0, -1), "infinity", 5)
        (x + 1, 1)
        sage: function_with_divisor(x^5 + x^2 + x, (0, 0), "infinity", 2)
        (x, 0)
    """
    K, L = _field(F)
    last = None
    for q in _places(K, L, Q):
        for p in _places(K, L, P):
            D = N * q.divisor() - N * p.divisor()
            basis = D.basis_function_space()
            if len(basis) == 1:
                A, B = basis[0].list()
                return _in_x(A), _in_x(B)
            last = len(basis)
    raise ValueError(
        "L(%d Q - %d P) has dimension %s, so %d(P-Q) is not principal"
        % (N, N, last, N))


def _order_from_multiple(divisor, multiple):
    """Exact order of a degree-zero divisor, if it divides ``multiple``."""
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
    r"""
    The least `N \le` ``bound`` with `N(P-Q)` principal, or ``None``.

    Both odd and even degree models are supported. Use
    ``("infinity", eta)`` to distinguish the infinities in even degree,
    where `eta = \lim y/x^{g+1}`. Bare ``"infinity"`` tries all rational
    infinite places and returns the smallest order among the choices.

    If ``multiple`` is supplied, test that multiple and remove its prime
    factors to find the exact order, without the ``bound`` cutoff. Returns
    ``None`` if that multiple does not annihilate the divisor class.

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
