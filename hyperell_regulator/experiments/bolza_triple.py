"""Compare L''(std,1) with a formal logarithmic determinant for a triple.

Run with: sage -python -m hyperell_regulator.experiments.bolza_triple
Uses the Beilinson experiment's completed-square model and signed periods.
L-value bounds are conditional on the backend's stated analytic assumptions.
The triple ((0,0), infinity, (2,sqrt(30))) fails the torsion condition for
the two-column construction: [(2,sqrt(30))-infinity] is nontorsion.
Thus D below is a formal integral determinant, not a certified motivic class.
A single Bloch cycle from the pair ((0,0),infinity) does exist; a third
point used only to normalize that cycle need not have torsion differences.
"""
import os
from fractions import Fraction

os.environ.setdefault('HYPELLFROB_THREADS', '2')

from sage.all import PolynomialRing, QQ, ZZ, RealField, matrix
from hyperell_regulator import HyperellipticCurveWithReduction, StandardLFunction
from hyperell_regulator.regulator import regulators_stokes, c_invariants
from hyperell_regulator.regulator.balls import narrow


class BolzaStandardL(StandardLFunction):
    # CM decomposition: L(std,s) = L(f_8,s) L(chi_8,s-1)^2 L(chi_-8,s-1),
    # where f_8 is the weight-3 newform 8.3.d.a. Thus N_std = 8^4.
    def conductor(self):
        return ZZ(4096)


def main():
    print("This triple fails the torsion condition; D is a formal logarithmic determinant.")
    digits, work = 20, 32
    x = PolynomialRing(QQ, 'x').gen()
    C = HyperellipticCurveWithReduction(x**5-x)
    T = PolynomialRing(QQ, 'T').gen()
    # bad_factors expects the exterior-square polynomial, BEFORE removing
    # (1-2T). The standard factor at 2 is consequently 1+2T.
    L = BolzaStandardL(C, bad_factors={2: 1-4*T**2}, h1_conductor=4096)
    ell = L.value(1, 2, digits=work, certified=True)

    # Only the rational x-coordinates enter these integrals. For common pole
    # R=(2,sqrt(30)), the logarithms are log|x/(x-2)| and -log|x-2|.
    per, (r0, r2) = regulators_stokes(C.F, [0, 2], prec=work, certified=True)
    B = per[0].parent()
    M = matrix(B, [[w*(r0[i]-r2[i]), -w*r2[i], w*per[i]]
                   for i, w in enumerate((1, 2, 1))])
    D = M.det()
    cp, cm = c_invariants(C.F, prec=work, certified=True)
    assert not D.contains_zero()
    q = B.pi()*B(cp)*B(cm)*B(ell)/D
    assert all(narrow(v, digits) for v in (D, cp, cm, ell, q))

    display = RealField(140)
    for name, value in (("D", D), ("L''(std,1)", ell),
                        ("pi*c+*c-*L''/D", q)):
        print(f"{name} = {display(value.mid()):.20f}")
    print("Ratio enclosure:", q)
    candidate = Fraction(str(display(q.mid()))).limit_denominator(10**6)
    rational = QQ(candidate.numerator)/candidate.denominator
    if abs(q-B(rational)).upper() < QQ(1)/10**digits:
        print("Rational candidate agreeing to 20 digits:", rational)
    else:
        print("No rational candidate with denominator <= 10^6 agrees to 20 digits.")


if __name__ == '__main__':
    main()
