"""Verify two explicit integral motivic classes for the Bolza curve.

Over K=Q(sqrt(2)), put e=1+sqrt(2) and let t1,t2 be the involutions
    t1(x,y)=((x+1)/(x-1), 2*sqrt(2)*y/(x-1)^3),
    t2(x,y)=((1-x)/(x+1), 2*sqrt(2)*y/(x+1)^3).
The classes are Cor_{K/Q}(Graph(tj),e) in CH^2(C x C,1).
The unit e makes these classes integral at every finite prime.
Their nonzero, independent primitive regulators are proved by the exact
correspondence matrices below. See reports/bolza_motivic_class.rst.

Run: sage -python -m hyperell_regulator.experiments.bolza_motivic
"""
from sage.all import (QuadraticField, PolynomialRing,
                      QQ, identity_matrix, matrix, vector)

from hyperell_regulator.regulator import periods_stokes
from hyperell_regulator.regulator.balls import CertificationError, digits_value, narrow


def _class_data():
    """Exact automorphisms, unit and pullback matrices, with their checks."""
    K = QuadraticField(2, 's')
    s = K.gen()
    epsilon = 1+s
    assert epsilon.is_integral() and epsilon.norm() == -1
    assert epsilon*(s-1) == 1

    A = PolynomialRing(K, 'x')
    Q = A.fraction_field()
    x = Q(A.gen())
    f = x**5-x
    maps = [((x+1)/(x-1), 2*s/(x-1)**3),
            ((1-x)/(x+1), 2*s/(x+1)**3)]
    matrices = []
    for u, scale in maps:
        # Exact automorphism identity, including the y-coordinate scale.
        assert u**5-u == scale**2*f
        a, b = A(u.derivative()/scale), A(u*u.derivative()/scale)
        M = matrix(K, [[a[0], a[1]], [b[0], b[1]]])
        assert M**2 == identity_matrix(K, 2) and M.trace() == 0
        assert matrix(K, 2, [c.galois_conjugate() for c in M.list()]) == -M
        matrices.append(M)

    assert matrix(K, [M.list() for M in matrices]).rank() == 2
    return epsilon, matrices


def regulator_matrix(digits=20):
    """Certified 3x3 matrix: Xi_1, Xi_2, and the period column.

    Uses Y^2=4(x^5-x) and weights (1,2 Re(x),|x|^2), exactly as the
    Beilinson experiment. The graph pairings are 2*log(epsilon)*G*T^t,
    where G is the unweighted period Gram matrix and T a pullback matrix.
    The factor two is the sum over the two real embeddings of K.
    """
    digits = digits_value(digits)
    _, matrices = _class_data()
    x = PolynomialRing(QQ, 'x').gen()
    periods = periods_stokes(4*(x**5-x), prec=digits+12, certified=True)
    B = periods[0].parent()
    G = matrix(B, [[periods[0], periods[1]], [periods[1], periods[2]]])
    root = B(2).sqrt()
    log_unit = (1+root).log()

    def embed(a):
        # Exact rational coefficients; choose the real embedding sqrt(2)>0.
        p = a.polynomial()
        return B(p[0])+B(p[1])*root

    columns = []
    for exact in matrices:
        T = matrix(B, 2, [embed(a) for a in exact.list()])
        H = 2*log_unit*G*T.transpose()
        columns.append(vector(B, [H[0,0], H[0,1]+H[1,0], H[1,1]]))
    columns.append(vector(B, [periods[0], 2*periods[1], periods[2]]))
    M = matrix(B, columns).transpose()
    if not all(narrow(a, digits) for a in M.list()+[M.det()]):
        raise CertificationError('graph regulator matrix did not attain the requested radius')
    if M.det().contains_zero():
        raise CertificationError('graph regulator determinant enclosure contains zero')
    return M


def lvalue_comparison(M, digits=20):
    """Certified comparison with pi*c+*c-*L''(std,1), in the same model."""
    from hyperell_regulator import HyperellipticCurveWithReduction
    from hyperell_regulator.regulator import c_invariants
    from .bolza_triple import BolzaStandardL
    from .beilinson import beilinson_comparison

    digits = digits_value(digits)
    work = digits+12
    x = PolynomialRing(QQ, 'x').gen()
    C = HyperellipticCurveWithReduction(x**5-x)
    T = PolynomialRing(QQ, 'T').gen()
    L = BolzaStandardL(C, bad_factors={2: 1-4*T**2}, h1_conductor=4096)
    ell = L.value(1, 2, digits=work, certified=True)
    cp, cm = c_invariants(C.F, prec=work, certified=True)
    comparison = beilinson_comparison(M.det(), cp, cm, ell,
                                     digits=digits, certified=True)
    if not narrow(comparison['ratio'], digits):
        raise CertificationError('L-value comparison did not attain the requested radius')
    return dict(L_second_derivative=ell, c_plus=cp, c_minus=cm,
                L_metadata=L._certified_metadata, comparison=comparison)


def main(argv=None):
    import argparse
    import hashlib
    import json
    from pathlib import Path
    from sage.all import RealField
    from hyperell_regulator.certified import encode_result

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--digits', type=int, default=20)
    parser.add_argument('--with-lvalue', action='store_true',
                        help='also compare the determinant with the standard L-value')
    parser.add_argument('--json', type=Path, help='save exact rational ball endpoints')
    args = parser.parse_args(argv)
    _, matrices = _class_data()
    for j, M in enumerate(matrices, 1):
        print('Class %d: Cor_{K/Q}(Graph(t%d), 1+sqrt(2))' % (j, j))
        print('  Pullback on (dx/y, x dx/y):\n', M)
    print('Both classes are integral: their function is an algebraic unit.')
    M = regulator_matrix(args.digits)
    B = M.base_ring()
    log_unit = (1+B(2).sqrt()).log()
    print('log(1+sqrt(2)):', log_unit)
    print('Columns: Xi_1, Xi_2, period; rows: 1, 2 Re(x), |x|^2.')
    display = RealField(B.precision())
    for row in M.rows():
        print('[ '+'  '.join(format(display(a.mid()), '.%df' % args.digits)
                              for a in row)+' ]')
    print('Determinant:', M.det())
    extra = lvalue_comparison(M, args.digits) if args.with_lvalue else {}
    if extra:
        comparison = extra['comparison']
        print("L''(std,1):", extra['L_second_derivative'])
        print("pi*c+*c-*L'' / det:", comparison['ratio'])
        print('Rational candidate:', comparison.get('rational'))
        print('Relative agreement bound:', comparison.get('relative_agreement'))
    if args.json:
        record = dict(curve='y^2=x^5-x', integration_model='Y^2=4(x^5-x)',
                      differentials=['dx/Y', 'x dx/Y'], weights=['1', '2 Re(x)', '|x|^2'],
                      columns=['Xi_1', 'Xi_2', 'period'], digits=args.digits, certified=True,
                      log_unit=log_unit, matrix=[list(row) for row in M.rows()],
                      regulator_matrix=[list(row[:2]) for row in M.rows()], det=M.det(),
                      normalization='(i/4) integral over the graph; both embeddings included',
                      source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        record.update(extra)
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(encode_result(record), indent=2, default=str)+'\n')
        print('Saved', args.json)


if __name__ == '__main__':
    main()
