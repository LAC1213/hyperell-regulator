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


def _certified_options(certified, cut):
    if not certified:
        return {}
    return dict(certified=True, **({"certified_cut": cut} if cut != "tree" else {}))


def _output(values, prec):
    return np.asarray(values, float) if prec is None else values


def regulators_stokes(F, centers, *, pole=None, prec=None, details=None, certified=False,
                      certified_cut="tree"):
    """Periods and logarithmic vectors on a common cut.

    Centers and finite poles are rational; ``pole=None`` means infinity.
    ``prec`` is the decimal accuracy target; ``None`` returns float64 arrays.
    With ``certified=True`` return real balls with checked absolute radii
    below ``10^-prec`` (15 digits when ``prec`` is omitted). The default
    certified cut is a tree with one infinity ray. ``certified_cut="parallel"``
    selects the independent parallel-ray implementation.
    """
    if certified:
        from . import balls
        coerce = balls.rational
    else:
        coerce = QQ
    centers = [coerce(a) for a in centers]
    pole = None if pole is None else coerce(pole)
    poles = [] if pole is None else [pole]
    unique = list(dict.fromkeys(centers+poles))
    if certified:
        digits = balls.digits_value(prec)
        per, reg = balls.compute(F, unique, digits+(5 if pole is not None else 0), details=details,
                                **({"cut": certified_cut} if certified_cut != "tree" else {}))
    else:
        per, reg = compute(F, unique, 15 if prec is None else prec, details=details)
    lookup = dict(zip(unique, reg))
    zero = [v*0 for v in per]
    base = zero if pole is None else lookup[QQ(pole)]
    columns = [([v.parent()(0) for v in per] if certified and c == pole else
                [a-b for a, b in zip(lookup[c], base)]) for c in centers]
    if certified:
        if not all(balls.narrow(v, digits) for v in per+[v for col in columns for v in col]):
            raise balls.CertificationError('finite-pole subtraction did not reach the requested radius')
        return per, columns
    return _output(per, prec), [_output(col, prec) for col in columns]


def periods_and_regulator_stokes(F, a=0, *, pole=None, prec=None, certified=False, certified_cut="tree"):
    """Periods and the logarithmic vector for ``a`` and an optional finite pole."""
    options = _certified_options(certified, certified_cut)
    per, reg = regulators_stokes(F, [a], pole=pole, prec=prec, **options)
    return per, reg[0]


def periods_stokes(F, *, prec=None, certified=False, certified_cut="tree"):
    """The unweighted plane period vector by Stokes integration."""
    options = _certified_options(certified, certified_cut)
    return regulators_stokes(F, [], prec=prec, **options)[0]


def regulator_stokes(F, a=0, *, pole=None, prec=None, certified=False, certified_cut="tree"):
    """The logarithmic plane vector by Stokes integration."""
    options = _certified_options(certified, certified_cut)
    return periods_and_regulator_stokes(F, a, pole=pole, prec=prec, **options)[1]


def periods_plane(F, *, certified=False, prec=None, certified_cut="tree", **kwds):
    """The unweighted plane vector by float64 quadrature."""
    if certified:
        if kwds:
            raise ValueError('float64 plane-quadrature options do not apply to certified integration')
        return periods_stokes(F, prec=prec, **_certified_options(True, certified_cut))
    return np.asarray(plane.regulator_vector_2(F, **kwds)[0], float)


def regulator_plane(F, a=0, *, pole=None, certified=False, prec=None, certified_cut="tree", **kwds):
    """The logarithmic plane vector by float64 quadrature."""
    if certified:
        if kwds:
            raise ValueError('float64 plane-quadrature options do not apply to certified integration')
        return regulator_stokes(F, a, pole=pole, prec=prec, **_certified_options(True, certified_cut))
    reg = np.asarray(plane.regulator_vector(F, a=float(a), **kwds)[0], float)
    if pole is not None:
        reg -= np.asarray(plane.regulator_vector(F, a=float(pole), **kwds)[0], float)
    return reg


def c_invariants(F, prec=None, *, certified=False):
    """The signed real and imaginary minors of full cycle periods."""
    if certified:
        from .balls import c_invariants as certified_invariants
        return certified_invariants(F, prec)
    if prec is None:
        return plane._c_invariants(F)
    from .precision import c_invariants_mp
    return c_invariants_mp(F, prec)
