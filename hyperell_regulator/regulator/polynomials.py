r"""
Polynomials in and out of the numeric layer.

The mathematics here is over `\QQ`: the curve is `y^2 = F(x)` and a function
on it is `\varphi = A(x) + B(x)\,y` with `A, B \in \QQ(x)`, so that is what
this package's entry points take -- Sage polynomials and rational functions,
not coefficient lists.

Underneath, the quadrature is numpy's and wants float64 coefficient arrays in
descending degree, which is what :func:`array` and :func:`arrays` produce.
The exact objects are kept beside them, because there are places where
rounding is fatal: the cover polynomial's coefficients are products of `F`,
`A` and `B`'s and cancel over several orders of magnitude, and the resultant
that gives the critical values has multiple roots that a numerical root finder
spreads over `\varepsilon^{1/m}`.

Coefficient lists in descending degree are still accepted everywhere a
polynomial is, together with a ``(numerator, denominator)`` pair of them where
a rational function is; that is how the package used to be addressed, and it
is also the only way to hand over the *inexact* data that
:meth:`~.cover.Cover.shifted` produces.

EXAMPLES::

    sage: from hyperell_regulator.regulator.polynomials import polynomial, rational
    sage: R = polynomial([1, 0, 0, 1, 2, 1]).parent(); R
    Univariate Polynomial Ring in x over Rational Field
    sage: polynomial([1, 0, 0, 1, 2, 1])
    x^5 + x^2 + 2*x + 1
    sage: rational(([1, 0], [1, -1]))
    x/(x - 1)
"""

import numpy as np

__all__ = ["ring", "polynomial", "rational", "descending", "array", "arrays"]


def ring(name="x"):
    r"""The polynomial ring `\QQ[x]` this package works in."""
    from sage.all import PolynomialRing, QQ

    return PolynomialRing(QQ, name)


def rational(v, name="x"):
    r"""
    ``v`` as an element of `\QQ(x)`, or ``None`` when it is not rational.

    ``v`` is a Sage polynomial or rational function, a rational number, a
    coefficient list in descending degree, or a ``(numerator, denominator)``
    pair of such lists.

    EXAMPLES::

        sage: from hyperell_regulator.regulator.polynomials import rational
        sage: rational([1, 1]), rational(0), rational(([1, 0], [1, -1]))
        (x + 1, 0, x/(x - 1))
        sage: rational([1j, 0]) is None
        True
    """
    from sage.all import QQ

    R = ring(name)
    K = R.fraction_field()
    if hasattr(v, "parent"):                       # already a Sage object
        try:
            return K(v)
        except (TypeError, ValueError, ZeroDivisionError, ArithmeticError):
            return None
    if isinstance(v, tuple) and len(v) == 2:
        n, d = rational(v[0], name), rational(v[1], name)
        return None if n is None or d is None or d.is_zero() else n / d
    try:
        c = [QQ(t) for t in
             np.atleast_1d(np.asarray(v, dtype=object)).ravel().tolist()]
    except (TypeError, ValueError, ArithmeticError):
        return None
    return K(R(list(reversed(c))))


def polynomial(v, name="x"):
    r"""
    ``v`` as an element of `\QQ[x]`, or ``None`` when it is not rational.

    Raises :exc:`ValueError` for a rational function with a genuine
    denominator, since the places that ask for this -- `F`, above all -- want
    a polynomial.

    EXAMPLES::

        sage: from hyperell_regulator.regulator.polynomials import polynomial
        sage: polynomial([1, 2, 1])
        x^2 + 2*x + 1
    """
    r = rational(v, name)
    if r is None:
        return None
    if r.denominator() != 1:                       # the fraction field is
        raise ValueError("%s is not a polynomial" % r)   # normalised: monic
    return r.numerator()


def descending(p):
    r"""
    The coefficients of ``p``, highest degree first, exactly.

    EXAMPLES::

        sage: from hyperell_regulator.regulator.polynomials import descending, ring
        sage: x = ring().gen()
        sage: descending(x^3 + 2)
        [1, 0, 0, 2]
    """
    return list(reversed(p.list())) or [p.base_ring().zero()]


def array(v):
    r"""
    A numpy coefficient array, highest degree first and with its leading
    zeros dropped.

    EXAMPLES::

        sage: from hyperell_regulator.regulator.polynomials import array, ring
        sage: array(ring().gen()^2 - 1)
        array([ 1.+0.j,  0.+0.j, -1.+0.j])
    """
    if hasattr(v, "parent") and hasattr(v, "list"):
        v = [complex(c) for c in descending(v)]
    a = np.atleast_1d(np.asarray(v)).ravel().astype(complex)
    nz = np.nonzero(a)[0]
    return a[nz[0]:] if len(nz) else np.zeros(1, complex)


def arrays(v):
    r"""
    ``(numerator, denominator)`` coefficient arrays of the rational function
    ``v``.

    EXAMPLES::

        sage: from hyperell_regulator.regulator.polynomials import arrays
        sage: [[float(t) for t in c.real] for c in arrays(([1, 0], [1, -1]))]
        [[1.0, 0.0], [1.0, -1.0]]
    """
    r = rational(v)
    if r is not None:
        return array(r.numerator()), array(r.denominator())
    if isinstance(v, tuple) and len(v) == 2:       # inexact: a pair of lists
        return array(v[0]), array(v[1])
    return array(v), np.ones(1, complex)


def evaluate(coefficients, value):
    """Horner evaluation without allocating numpy arrays for scalar points."""
    if np.ndim(value):
        return np.polyval(coefficients, value)
    z = complex(value)
    out = 0j
    for coefficient in coefficients:
        out = out * z + complex(coefficient)
    # Retain numpy's non-finite arithmetic at poles (rather than raising
    # Python's ZeroDivisionError before the chart can handle the pole).
    return np.complex128(out)
