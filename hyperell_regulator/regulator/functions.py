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
    The places of `L` at ``point``, either ``"infinity"`` or `(x_0, y_0)`.

    In even degree infinity is two places, and which of them a function with
    divisor `N P - N Q` wants is settled by trying: at most one of them makes
    `N(P-Q)` principal, so :func:`function_with_divisor` takes whichever does
    rather than being told which.
    """
    from sage.all import QQ

    x = K.gen()
    y = L.gen()
    if point == "infinity":
        return list(L.places_infinite())
    x0, y0 = QQ(point[0]), QQ(point[1])
    below = K.maximal_order().ideal(x - x0).place()
    for pl in L.places_above(below):
        if (y - y0).valuation(pl) > 0:
            return [pl]
    raise ValueError("%r is not a point of the curve" % (point,))


def _place(K, L, point):
    r"""The place at ``point``, which must be a single one."""
    pls = _places(K, L, point)
    if len(pls) != 1:
        raise ValueError("infinity is not a single place of this model")
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

    `F` is a polynomial over `\QQ`; `P` and `Q` are ``"infinity"`` or a pair
    `(x_0, y_0)` of rationals.  Raises :exc:`ValueError` when `N(P-Q)` is not
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


def torsion_order(F, P, Q, bound=24):
    r"""
    The least `N \le` ``bound`` with `N(P-Q)` principal, or ``None``.

    EXAMPLES::

        sage: from hyperell_regulator.regulator import torsion_order
        sage: from hyperell_regulator.regulator.polynomials import ring
        sage: x = ring().gen()
        sage: torsion_order(x^5 + x^2 + 2*x + 1, (0, -1), "infinity")
        5
    """
    K, L = _field(F)
    p, q = _place(K, L, P), _place(K, L, Q)
    for N in range(1, bound + 1):
        if len((N * q.divisor() - N * p.divisor()).basis_function_space()) == 1:
            return N
    return None
