"""Guarded full cycle periods shared by both numerical output modes."""
from .quadrature import Arithmetic, quadrature
from .polynomials import polynomial
from ._periods import period_paths, segment_differential, pivoted_minor


def period_matrix_mp(F, prec):
    """Full periods, twice the lifted branch-to-branch integrals."""
    ar = Arithmetic(prec)
    exact = polynomial(F)
    if exact is None or exact.gcd(exact.derivative()).degree() > 0:
        raise ValueError('F must be an exact squarefree polynomial')
    g = (exact.degree()-1)//2
    roots = exact.roots(ar.C, multiplicities=False)
    output = []
    for (a, ia), (b, ib) in period_paths(roots, g, root_tol=ar.root_tol):
        evaluate, _ = segment_differential(
            roots, ar.C(exact.leading_coefficient()), a, b, ia, ib, g,
            sqrt=lambda z: z.sqrt(), sin=lambda z: z.sin(),
            cos=lambda z: z.cos(), pi=ar.R.pi(), sign_tol=ar.root_tol)
        output.append([2*v for v in quadrature(ar, evaluate, g)])
    return output


def c_invariants_mp(F, prec):
    ar = Arithmetic(prec)
    periods = period_matrix_mp(F, prec)
    return tuple(ar.output(pivoted_minor([[part(v) for v in row] for row in periods], ar.tol))
                 for part in (lambda z: z.real(), lambda z: z.imag()))
