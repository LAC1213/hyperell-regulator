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
    """Branch-segment periods with endpoint singularities removed analytically.

    Set x=a+(b-a)*sin(pi*t/2)^2. Factoring out the two endpoint
    roots cancels the square-root singularities exactly, and continuing
    each remaining linear factor fixes the square-root branch on the segment.
    """
    ar=Arithmetic(prec);C,R=ar.C,ar.R
    exact=polynomial(F)
    if exact is None:
        raise ValueError('arbitrary precision requires exact polynomial coefficients')
    p=PolynomialRing(C,'x')(exact)
    g=(p.degree()-1)//2
    roots=sorted(p.roots(C,multiplicities=False),key=lambda z:(float(z.real()),float(z.imag())))
    if len(roots)!=p.degree():raise ValueError('the curve polynomial must be squarefree')
    output=[]
    for j in range(2*g):
        a,b=roots[j],roots[j+1];delta=b-a;mid=(a+b)/2
        others=[v for k,v in enumerate(roots) if k not in (j,j+1)]
        qmid=-delta**2*p.leading_coefficient()
        for v in others:qmid*=mid-v
        ymid=qmid.sqrt()
        # Conjugate endpoints make qmid real. Tiny rounding errors on the
        # negative real axis must not reverse this period's orientation.
        if abs(ymid.real())<=ar.newton_tol*abs(ymid) and ymid.imag()<0:
            ymid=-ymid
        def evaluate(t):
            x=a+delta*(R.pi()*t/2).sin()**2
            y=ymid
            for v in others:y*=((x-v)/(mid-v)).sqrt()
            return [R.pi()*delta*x**i/y for i in range(g)]
        row=quadrature(ar,FunctionSamples(evaluate),g)
        output.append(row)
    return output


def _pivot_det(rows):
    a=[list(row) for row in rows]
    n,m=len(a),len(a[0]);det=a[0][0].parent()(1)
    for k in range(m):
        p=max(range(k,n),key=lambda i:abs(a[i][k]))
        if p!=k:a[k],a[p]=a[p],a[k];det=-det
        if not a[k][k]:return det.parent()(0)
        det*=a[k][k]
        for i in range(k+1,n):
            factor=a[i][k]/a[k][k]
            for j in range(k+1,m):a[i][j]-=factor*a[k][j]
    return det


def c_invariants_mp(F,prec):
    """Pivoted real and imaginary period minors at the requested precision."""
    ar=Arithmetic(prec)
    P=period_matrix_mp(F,prec)
    return tuple(ar.output(_pivot_det([[part(v) for v in row] for row in P]))
                 for part in (lambda z:z.real(),lambda z:z.imag()))
