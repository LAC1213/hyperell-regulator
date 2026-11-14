"""Exact curve polynomials; coefficient lists use descending degree."""
from sage.all import PolynomialRing, QQ


def ring(name="x"):
    return PolynomialRing(QQ, name)


def polynomial(value, name="x"):
    """Coerce a Sage polynomial or descending coefficient list to QQ[x]."""
    R = ring(name)
    try:
        if hasattr(value, 'parent') or isinstance(value, str):
            return R(value)
        if isinstance(value, (list, tuple)):
            return R(list(reversed(value)))
        return R(value)
    except (TypeError, ValueError, ArithmeticError):
        return None
