"""Full cycle periods with explicit branch-point endpoint substitutions."""
from sage.all import PolynomialRing
from .quadrature import Arithmetic, quadrature
from .polynomials import polynomial
from ._periods import period_paths, segment_differential, pivoted_minor


def period_matrix_mp(F, prec):
    """Full periods, twice the lifted branch-to-branch integrals."""
    ar = Arithmetic(prec)
    exact = polynomial(F)
    if exact is None or exact.gcd(exact.derivative()).degree() > 0:
        raise ValueError('F must be an exact squarefree polynomial')
    p = PolynomialRing(ar.C, 'x')(exact)
    g = (p.degree()-1)//2
    roots = p.roots(ar.C, multiplicities=False)
    output = []
    for path in period_paths(roots, g, root_tol=ar.root_tol):
        row = [ar.C(0) for _ in range(g)]
        previous_y = None
        for (a, ia), (b, ib) in zip(path[:-1], path[1:]):
            evaluate, previous_y = segment_differential(
                roots, p.leading_coefficient(), ar.C(a), ar.C(b), ia, ib, g,
                sqrt=lambda z: z.sqrt(), sin=lambda z: z.sin(),
                cos=lambda z: z.cos(), pi=ar.R.pi(),
                sign_tol=ar.root_tol, previous_y=previous_y)
            part = quadrature(ar, evaluate, g)
            row = [u+v for u, v in zip(row, part)]
        output.append([2*v for v in row])
    return output


def c_invariants_mp(F, prec):
    ar = Arithmetic(prec)
    periods = period_matrix_mp(F, prec)
    return tuple(ar.output(pivoted_minor([[part(v) for v in row] for row in periods], ar.tol))
                 for part in (lambda z: z.real(), lambda z: z.imag()))
