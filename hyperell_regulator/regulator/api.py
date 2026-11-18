r"""Plane periods and integrals weighted by `\log|x-a|` on `y^2=F(x)`.

Entries have weights `\Re(x^i\bar x^j)/|F(x)|`, with `i\le j`.
In genus two the order is `(1,\Re(x),|x|^2)`; double the middle entry
for the Beilinson determinant. A finite pole `b` gives `\log|(x-a)/(x-b)|`.
"""
import numpy as np
from sage.all import QQ
from . import plane
from .stokes import compute

__all__ = ['c_invariants', 'periods_plane', 'periods_stokes',
           'periods_and_regulator_stokes', 'regulators_stokes',
           'regulator_plane', 'regulator_stokes']


def _output(values, prec):
    return np.asarray(values, float) if prec is None else values


def regulators_stokes(F, centers, *, pole=None, prec=None, details=None):
    """Periods and logarithmic vectors on a common cut.

    Centers and finite poles are rational; ``pole=None`` means infinity.
    ``prec`` is the decimal accuracy target; ``None`` returns float64 arrays.
    """
    centers = [QQ(a) for a in centers]
    poles = [] if pole is None else [QQ(pole)]
    unique = list(dict.fromkeys(centers+poles))
    per, reg = compute(F, unique, 15 if prec is None else prec, details=details)
    lookup = dict(zip(unique, reg))
    zero = [v*0 for v in per]
    base = zero if pole is None else lookup[QQ(pole)]
    return _output(per, prec), [_output([a-b for a, b in zip(lookup[c], base)], prec) for c in centers]


def periods_and_regulator_stokes(F, a=0, *, pole=None, prec=None):
    """Periods and the logarithmic vector for ``a`` and an optional finite pole."""
    per, reg = regulators_stokes(F, [a], pole=pole, prec=prec)
    return per, reg[0]


def periods_stokes(F, *, prec=None):
    """The unweighted plane period vector by Stokes integration."""
    return regulators_stokes(F, [], prec=prec)[0]


def regulator_stokes(F, a=0, *, pole=None, prec=None):
    """The logarithmic plane vector by Stokes integration."""
    return periods_and_regulator_stokes(F, a, pole=pole, prec=prec)[1]


def periods_plane(F, **kwds):
    """The unweighted plane vector by float64 quadrature."""
    return np.asarray(plane.regulator_vector_2(F, **kwds)[0], float)


def regulator_plane(F, a=0, *, pole=None, **kwds):
    """The logarithmic plane vector by float64 quadrature."""
    reg = np.asarray(plane.regulator_vector(F, a=float(a), **kwds)[0], float)
    if pole is not None:
        reg -= np.asarray(plane.regulator_vector(F, a=float(pole), **kwds)[0], float)
    return reg


def c_invariants(F, prec=None):
    """The signed real and imaginary minors of full cycle periods."""
    if prec is None:
        return plane._c_invariants(F)
    from .precision import c_invariants_mp
    return c_invariants_mp(F, prec)
