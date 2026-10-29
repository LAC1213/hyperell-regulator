"""Arbitrary-precision periods and regulators using the shared Stokes algorithm.

``prec`` is a decimal accuracy target. Exact rational coefficients are read
without a float64 round trip. The discrete cut topology, chart preflight, and initial sheet
labels are shared with float64; Newton continuation, endpoint coordinates,
quadrature, error checks and assembly use guarded arbitrary precision.

Both vectors are computed from one integration of the cut edges. There is
no float64 fallback and a failed convergence test raises an exception.
"""
from sage.all import PolynomialRing
from ._mp_continuation import Arithmetic, FunctionSamples, quadrature
from ._mp_stokes import compute
from .polynomials import polynomial

__all__ = ['periods_and_regulator_stokes_mp','periods_stokes_mp',
           'regulator_stokes_mp','period_matrix_mp','c_invariants_mp']


def periods_and_regulator_stokes_mp(F,phi,prec,check_tol=1e-6):
    """Return both vectors in a Sage real field, integrating each edge once."""
    return compute(F,phi,prec,check_tol)


def periods_stokes_mp(F,phi,prec,check_tol=1e-6):
    """Unweighted period vector, to ``prec`` decimal digits."""
    return periods_and_regulator_stokes_mp(F,phi,prec,check_tol)[0]


def regulator_stokes_mp(F,phi,prec,check_tol=1e-6):
    """Logarithmically weighted regulator vector, to ``prec`` decimal digits."""
    return periods_and_regulator_stokes_mp(F,phi,prec,check_tol)[1]


def period_matrix_mp(F,prec):
    """Full cycle periods on the shared, validated semilunar paths.

    The first root has largest absolute value, with real and imaginary parts
    breaking ties at working precision. Endpoint singularities are removed
    analytically, and the square root is continued through detour corners.
    """
    from ._periods import period_paths, segment_differential

    ar=Arithmetic(prec);C,R=ar.C,ar.R
    exact=polynomial(F)
    if exact is None:
        raise ValueError('arbitrary precision requires exact polynomial coefficients')
    p=PolynomialRing(C,'x')(exact)
    g=(p.degree()-1)//2
    roots=p.roots(C,multiplicities=False)
    if len(roots)!=p.degree():raise ValueError('the curve polynomial must be squarefree')
    output=[]
    for path in period_paths(roots,g,root_tol=ar.newton_tol):
        row=[C(0) for _ in range(g)]
        previous_y=None
        for (a,ia),(b,ib) in zip(path[:-1],path[1:]):
            evaluate,previous_y=segment_differential(
                roots,p.leading_coefficient(),C(a),C(b),ia,ib,g,
                sqrt=lambda v:v.sqrt(),sin=lambda v:v.sin(),
                cos=lambda v:v.cos(),pi=R.pi(),sign_tol=ar.newton_tol,
                previous_y=previous_y)
            part=quadrature(ar,FunctionSamples(evaluate),g)
            row=[u+v for u,v in zip(row,part)]
        output.append([2*v for v in row])
    return output


def c_invariants_mp(F,prec):
    """Signed pivoted minors of full cycle periods at the requested precision."""
    from ._periods import pivoted_minor

    ar=Arithmetic(prec)
    P=period_matrix_mp(F,prec)
    # Equal period pivots can differ by the quadrature error, which is much
    # larger than Newton's tolerance. Use the integration accuracy for ties.
    return tuple(ar.output(pivoted_minor([[part(v) for v in row] for row in P],ar.tol))
                 for part in (lambda z:z.real(),lambda z:z.imag()))
