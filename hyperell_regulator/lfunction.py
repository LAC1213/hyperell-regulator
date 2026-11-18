r"""Curve and standard `L`-functions of genus two curves over `\QQ`.

The standard motive is `\ker(\wedge^2 h^1(C) \to \QQ(-1))`, so
`L(\wedge^2 h^1(C),s)=L(M_{\mathrm{std}}(C),s)\zeta(s-1)`.

EXAMPLES::

    sage: from hyperell_regulator.lfunction import StandardLFunction
    sage: from hyperell_regulator.reduction import HyperellipticCurveWithReduction
    sage: R.<x> = ZZ[]
    sage: C = HyperellipticCurveWithReduction(x^5 + x^3 + 1)
    sage: L = StandardLFunction(C, bad_factors={2: None})  # doctest: +ELLIPSIS
    ...Warning: unknown valuation of conductor at 2.
    sage: L.euler_factor(53).degree()
    3
"""

from contextlib import contextmanager
from struct import calcsize

from cypari2 import PariError

from sage.all import (
    prime_range,
    GF,
    HyperellipticCurve,
    Matrix,
    PolynomialRing,
    PowerSeriesRing,
    QQ,
    RR,
    SageObject,
    ZZ,
    gp,
    pari,
)
from sage.modules.with_basis.subquotient import SubmoduleWithBasis

from hyperell_regulator.reduction import potential_toric_rank


def _bits(digits):
    r"""Bit precision for ``digits`` decimal digits, with guard bits."""
    return int(3.33 * int(digits)) + 64


def _pari_lfun_call(function, *args, coefficient_count=0, digits=10, **kwargs):
    """Evaluate with a PARI stack budget proportional to the coefficient count.

    On stack exhaustion, double the ceiling and retry.
    """
    if coefficient_count:
        word_bytes = calcsize('l')
        precision_words = (_bits(digits) + 8*word_bytes - 1) // (8*word_bytes)
        budget = int(coefficient_count) * word_bytes * (64 + 8*precision_words)
        ceiling = 1 << max(26, (budget-1).bit_length())
        if ceiling > pari.stacksizemax():
            pari.allocatemem(pari.stacksize(), ceiling, silent=True)
    while True:
        try:
            return function(*args, **kwargs)
        except PariError as exc:
            if str(pari.errname(exc.errdata())) != 'e_STACK':
                raise
            size = pari.stacksize()
            ceiling = max(size, pari.stacksizemax())
            pari.allocatemem(size, 2*ceiling, silent=True)


def _real_field(digits):
    r"""A real field with guard bits for ``digits`` decimal digits.

    EXAMPLES::

        sage: from hyperell_regulator.lfunction import _real_field
        sage: _real_field(30).precision() >= 100
        True
    """
    from sage.all import RealField

    return RealField(int(RR(digits) * RR(10).log2()) + 12)


@contextmanager
def _pari_precision(digits):
    """Set PARI precision for the context and restore it on exit."""
    previous = pari.get_real_precision()
    pari.set_real_precision(digits)
    try:
        yield
    finally:
        pari.set_real_precision(previous)


def _coefficient_count(conductor, gamma, weight, digits, derivative):
    """PARI's coefficient count at the requested decimal precision."""
    previous = gp.eval('default(realprecision)')
    gp.eval('default(realprecision, %d)' % int(digits))
    try:
        L = pari.lfuncreate([[QQ(0)], 0, gamma, weight, conductor, 1])
        cost = gp.eval('lfuncost(%s, [1,0,0], %d)' % (L, int(derivative)))
        return int(ZZ(str(cost).strip('[]').split(',')[0].strip()))
    finally:
        gp.eval('default(realprecision, %s)' % previous)


def _dirichlet_coefficients(L, num):
    """Dirichlet coefficients from Euler factors, with a zero sentinel at index 0."""
    a = [ZZ(0)] * (num + 1)
    a[1] = ZZ(1)
    L.precompute_h1(prime_range(num + 1))

    # Expand 1/P_p once per prime, up to the greatest power of p <= num.
    inv = {}

    def local(p, e):
        if p not in inv:
            depth = ZZ(num).exact_log(p)
            RS = PowerSeriesRing(QQ, 'x', default_prec=depth + 1)
            x = RS.gen()
            inv[p] = (1 / RS(L.euler_factor(p)(x))).padded_list(depth + 1)
        return inv[p][e]

    for n in range(2, num + 1):
        coeff = QQ(1)
        for p, e in ZZ(n).factor():
            coeff *= local(p, e)
        a[n] = coeff
    return a


def _precompute_h1(L, primes, nthreads):
    """Cache `h^1` factors at good odd primes in one batch."""
    from .frobenius import frobenius_polynomials

    bad = set(L.C.bad_primes())
    todo = [int(p) for p in primes
            if int(p) != 2 and ZZ(p).is_prime() and ZZ(p) not in bad
            and ZZ(p) not in L._h1cache]
    if not todo:
        return
    for p, cp in frobenius_polynomials(L.C.F, todo, nthreads).items():
        L._h1cache[ZZ(p)] = L._RT(list(reversed(cp.list())))


def wedge2_euler_factor(P):
    r"""The Euler polynomial `\prod_{i<j}(1-\alpha_i\alpha_j T)`
    from ``P`` equal to `\prod_i(1-\alpha_i T)` of degree four.

    EXAMPLES::

        sage: from hyperell_regulator.lfunction import wedge2_euler_factor
        sage: S = PolynomialRing(QQ, ['a1','a2','a3','a4'])
        sage: a1, a2, a3, a4 = S.gens()
        sage: RT = PolynomialRing(S, 'T'); T = RT.gen()
        sage: P = prod(1 - a*T for a in (a1, a2, a3, a4))
        sage: W = prod(1 - a*b*T for a, b in Combinations([a1,a2,a3,a4], 2))
        sage: wedge2_euler_factor(P) == W
        True
    """
    RT = P.parent()
    T = RT.gen()
    c = [P[i] for i in range(5)]
    s1, s2, s3, s4 = -c[1], c[2], -c[3], c[4]
    d = [1,
         -s2,
         s1 * s3 - s4,
         -(s1 ** 2 * s4 - 2 * s2 * s4 + s3 ** 2),
         s1 * s3 * s4 - s4 ** 2,
         -s2 * s4 ** 2,
         s4 ** 3]
    return sum(d[i] * T ** i for i in range(7))


class StandardLFunction(SageObject):
    r"""The degree five `L`-function of `M_{\mathrm{std}}(C)`.

    Its local factor is `P_p^{\wedge^2}(T)/(1-pT)`. Bad factors use
    monodromy and inertia data; unsupported cases require ``bad_factors``.
    ``h1_conductor`` supplies the exact Jacobian conductor when needed at 2.
    The functional equation is `s \leftrightarrow 3-s`, with sign `+1`.

    EXAMPLES::

        sage: from hyperell_regulator.lfunction import StandardLFunction
        sage: from hyperell_regulator.reduction import HyperellipticCurveWithReduction
        sage: R.<x> = ZZ[]
        sage: C = HyperellipticCurveWithReduction(x^5 + x^3 + 1)
        sage: L = StandardLFunction(C, bad_factors={2: None})   # doctest: +ELLIPSIS
        ...Warning: unknown valuation of conductor at 2.
        sage: L.euler_factor(53).degree()
        3
    """

    def __init__(self, C, bad_factors=None, prec=None, h1_conductor=None):
        self.C = C
        if C.genus() != 2:
            raise ValueError("the standard L-function here is for genus 2; "
                             "this curve has genus %s" % C.genus())
        self.bad_factors = dict(bad_factors or {})
        self.h1_conductor = None if h1_conductor is None else ZZ(h1_conductor)
        if self.h1_conductor is not None and self.h1_conductor < 1:
            raise ValueError("h1_conductor must be positive")
        RQ = PolynomialRing(QQ, 'x')
        self._pari_L = pari.lfungenus2([pari(str(RQ(C.f))), pari(str(RQ(C.h)))])
        self._RT = PolynomialRing(QQ, 'T')
        self._cache = {}
        self._h1cache = {}

    def _repr_(self):
        return "Standard L-function of %s" % self.C

    # -- Euler factors -----------------------------------------------------

    def h1_euler_factor(self, p):
        r"""The `h^1` Euler factor at ``p``; bad primes and 2 use PARI."""
        p = ZZ(p)
        if p in self._h1cache:
            return self._h1cache[p]
        if p != 2 and p not in self.C.bad_primes():
            H = HyperellipticCurve(self.C.F.change_ring(GF(p)))
            cp = list(H.frobenius_polynomial())        # ascending, leading coeff 1
            out = self._RT(list(reversed(cp)))          # -> Euler-factor convention
            self._h1cache[p] = out
            return out
        e = self._pari_L.lfuneuler(p)
        num, den = pari(e).numerator(), pari(e).denominator()
        # lfuneuler returns 1/P_p(x); take the denominator
        poly = den if str(num).strip() in ("1", "1.0") else den / num
        return self._RT([QQ(c) for c in pari(poly).Vecrev()])

    def precompute_h1(self, primes, nthreads=None):
        r"""Cache `h^1` factors at good odd primes in one batch.

        EXAMPLES::

            sage: from sage.all import ZZ, PolynomialRing
            sage: from hyperell_regulator import HyperellipticCurveWithReduction
            sage: from hyperell_regulator.lfunction import StandardLFunction
            sage: R = PolynomialRing(ZZ, 'x'); x = R.gen()
            sage: L = StandardLFunction(HyperellipticCurveWithReduction(x**5 + x**3 + 1))
            ...Warning: unknown valuation of conductor at 2.
            sage: L.precompute_h1([11, 13, 17])
            sage: sorted(L._h1cache)
            [11, 13, 17]
        """
        _precompute_h1(self, primes, nthreads)

    def toric_factor(self, p):
        r"""`\det(1-\mathrm{Frob}_p T \mid H_1(\Gamma))` for the dual graph `\Gamma`."""
        R = self.C.cluster_picture(p)
        H1, M, frob = R.homology_of_special_fibre()
        basis = list(H1.basis())
        t = len(basis)
        if t == 0:
            return self._RT(1)
        rows = []
        for b in basis:
            if isinstance(H1, SubmoduleWithBasis):
                # The cluster package defines its map on the generators
                # passed to submodule(), but Sage echelonizes that basis.
                # Act on the actual lifted cycles, then retract to H1.
                # Merely retracting frob(b) can give an invertible but
                # incorrect matrix, depending on the cluster ordering.
                image = H1.ambient().sum_of_terms(
                    (s.frobenius(), n*s.epsilon(
                        lambda c: c.frobenius(), R.field_frobenius()))
                    for s, n in H1.lift(b))
                image = H1.retract(image)
            else:
                image = frob(b)
            rows.append([QQ(image.coefficient(c.support()[0])) for c in basis])
        Mfrob = Matrix(QQ, t, t, rows)
        if Mfrob.det() not in (-1, 1):
            raise ArithmeticError("Frobenius on graph homology is not invertible at %s" % p)
        return self._RT(Mfrob.charpoly().reverse())

    def _wedge2_totally_toric(self, p, U):
        r"""The Euler factor of `\wedge^2 U \oplus \mathrm{Sym}^2 U(-1)`."""
        if U.degree() != 2 or U[0] != 1 or U[2] not in (-1, 1):
            raise ArithmeticError("invalid rank-two toric Euler factor at %s: %s" % (p, U))
        T = self._RT.gen()
        trace, det = -U[1], U[2]
        X = p*T
        symmetric = 1-(trace**2-det)*X+det*(trace**2-det)*X**2-det**3*X**3
        return (1-det*T)*symmetric

    def _wedge2_at_two(self, P):
        r"""The exterior-square factor at 2 when `h^1` has nonzero inertia invariants.

        For a ramified symplectic complement of dimension two, both determinant
        lines have eigenvalue 2, giving `(1-2T)^2`.
        """
        T = self._RT.gen()
        if P.degree() == 4:
            return wedge2_euler_factor(P)
        if P.degree() == 3:
            for gamma in (1, -1):
                A, remainder = P.quo_rem(1-gamma*T)
                if not remainder and A.degree() == 2 and A[2] == 2 and abs(A[1]) <= 2:
                    return A(gamma*T)*(1-2*T)**2
        if P in ((1-T)**2, (1+T)**2, 1-T**2, 1+T**2,
                 1+T+T**2, 1-T+T**2):
            return self._wedge2_totally_toric(2, P)
        if ((P.degree() == 2 and P[2] == 2 and abs(P[1]) <= 2)
                or (P.degree() == 1 and P[1] in (-1, 1))):
            self._cache["two_adic_ramified_complement"] = True
            return (1-2*T)**2
        raise NotImplementedError(
            "the h^1 factor at 2 is %s; the exterior-square factor requires "
            "additional inertia data in this case; supply bad_factors[2]" % P)

    def _h1_conductor_exponent(self, p):
        """The local Jacobian conductor exponent, using the supplied conductor if needed."""
        if self.h1_conductor is not None:
            return self.h1_conductor.valuation(p)
        exponent = self.C.local_data(p)["conductor_exponent"]
        if exponent is None:
            raise NotImplementedError(
                "the wild conductor at %s requires the exact Jacobian "
                "conductor; pass h1_conductor to StandardLFunction" % p)
        return ZZ(exponent)

    def wedge2_euler_factor(self, p):
        r"""
        The `\wedge^2 h^1` Euler factor at ``p``.
        """
        if p in self.bad_factors and self.bad_factors[p] is not None:
            return self._RT(self.bad_factors[p])
        T = self._RT.gen()
        p = ZZ(p)
        P = self.h1_euler_factor(p)
        if p not in self.C.bad_primes():
            return wedge2_euler_factor(P)

        if p == 2:
            return self._wedge2_at_two(P)
        if self.C.is_semistable_at(p) is not True:
            return self._wedge2_potentially_good(p, P)

        U = self.toric_factor(p)
        t = U.degree()
        if t == 0:
            W = wedge2_euler_factor(P)
        elif t == 1:
            gamma = -U[1]                       # U(T) = 1 - gamma T
            A = P // U                          # abelian factor, degree 2
            if A * U != P:
                raise ArithmeticError(
                    "at %s the toric factor does not divide the h^1 factor "
                    "(%s does not divide %s)" % (p, U, P))
            b2 = A[2]
            W = A(gamma * T) * (1 - b2 * T) * (1 - gamma ** 2 * p * T)
        elif t == 2:
            W = self._wedge2_totally_toric(p, U)
        else:
            raise NotImplementedError("toric rank %s at %s" % (t, p))
        return W

    def _wedge2_potentially_good(self, p, P):
        r"""The exterior-square factor for tame potentially good reduction.

        Implemented for conductor exponent 0 or 2.
        """
        T = self._RT.gen()
        if self.C.is_tame(p) is not True:
            raise NotImplementedError(
                "C is neither semistable nor known to be tame at %s "
                "(is_tame = %s); supply bad_factors[%s]"
                % (p, self.C.is_tame(p), p))
        stable_type = self.C.local_data(p)["stable_type"]
        t = potential_toric_rank(stable_type)
        if t != 0:
            raise NotImplementedError(
                "C is not semistable at %s and its potential stable reduction "
                "has toric rank %s (Liu type (%s)); only potentially good "
                "reduction is implemented in the non-semistable case, since "
                "toric rank > 0 needs (ker N)^Phi.  Supply bad_factors[%s]"
                % (p, t, stable_type, p))

        fp = ZZ(self.C.conductor()).valuation(p)
        if fp != 4 - P.degree():
            raise ArithmeticError(
                "at %s the conductor exponent %s and deg P_p = %s disagree; "
                "for tame reduction f_p = 4 - dim V^I" % (p, fp, P.degree()))
        if fp == 0:
            return wedge2_euler_factor(P)
        if fp == 2:
            d = P[2]                        # det(Frob | V^I), P_p of degree 2
            if d != p:
                raise ArithmeticError(
                    "at %s the potentially good wedge^2 factor needs "
                    "det(Frob | V^I) = %s, but P_p = %s has leading "
                    "coefficient %s" % (p, p, P, d))
            return (1 - p * T) ** 2
        raise NotImplementedError(
            "potentially good reduction at %s with f_%s = %s: the wedge^2 "
            "invariants are not determined by P_p (which is %s) -- they need "
            "the eigenvalues of inertia on V_ell.  Supply bad_factors[%s]"
            % (p, p, fp, P, p))

    def euler_factor(self, p):
        r"""
        The standard Euler factor `P_p^{\wedge^2}(T)/(1 - pT)`.
        """
        p = ZZ(p)
        if p in self._cache:
            return self._cache[p]
        T = self._RT.gen()
        W = self.wedge2_euler_factor(p)
        q, r = W.quo_rem(1 - p * T)
        if r != 0:
            raise ArithmeticError(
                "(1 - %s T) does not divide the wedge^2 factor %s at %s; the "
                "Q(-1) summand of wedge^2 h^1 must contribute it"
                % (p, W, p))
        self._cache[p] = q
        return q

    def conductor(self):
        r"""The standard conductor. Tame exponents are `5-\deg P_p^{\mathrm{std}}`.

        For the supported ramified complements at 2, add twice the `h^1` Swan conductor.
        """
        N = ZZ(1)
        for p in self.C.bad_primes():
            exponent = 5-self.euler_factor(p).degree()
            if p == 2 and self._cache.get("two_adic_ramified_complement"):
                swan = self._h1_conductor_exponent(p)-(4-self.h1_euler_factor(p).degree())
                if swan < 0:
                    raise ArithmeticError("the supplied H1 conductor at 2 is smaller than its tame part")
                exponent += 2*swan
            N *= ZZ(p)**exponent
        return N

    # -- the L-function ----------------------------------------------------

    def dirichlet_coefficients(self, num):
        r"""The Dirichlet coefficients indexed from 1 through ``num``."""
        return _dirichlet_coefficients(self, num)[1:]

    def _coeffs_needed(self, digits, der=2):
        r"""PARI's coefficient count for derivative ``der`` at ``digits`` decimal digits."""
        return _coefficient_count(self.conductor(), [-1, -1, 0, 0, 1], 3, digits, der)

    def lfun(self, num_coeffs=None, digits=10):
        r"""The PARI object with gamma parameters `[-1,-1,0,0,1]` and weight 3.

        ``digits`` is the decimal working precision.
        """
        key = ("lfun", num_coeffs, digits)
        if key in self._cache:
            return self._cache[key]
        N = self.conductor()
        if num_coeffs is None:
            num_coeffs = max(200, self._coeffs_needed(digits))
        an = self.dirichlet_coefficients(num_coeffs)
        with _pari_precision(digits):
            L = _pari_lfun_call(
                pari.lfuncreate, [[QQ(c) for c in an], 0,
                                 [-1, -1, 0, 0, 1], 3, N, 1],
                coefficient_count=len(an), digits=digits)
        self._cache[key] = L
        return L

    def check_functional_equation(self, num_coeffs=None, digits=10):
        r"""The functional equation residual ``2^lfuncheckfeq``."""
        L = self.lfun(num_coeffs, digits)
        K = _real_field(digits)
        with _pari_precision(digits):
            return K(2) ** K(_pari_lfun_call(
                pari.lfuncheckfeq, L, precision=_bits(digits)))

    def vanishing_order(self, s=1, max_order=4, num_coeffs=None, digits=10,
                        tol=None):
        r"""The numerical order of vanishing at ``s``.

        Derivatives below ``tol`` are treated as zero; the default is `10^{-\mathrm{digits}/2}`.
        """
        if tol is None:
            tol = _real_field(digits)(10) ** (-_real_field(digits)(digits) / 2)
        for k in range(max_order + 1):
            if abs(self.value(s, k, num_coeffs, digits)) > tol:
                return k
        return None

    def value(self, s=1, derivative=2, num_coeffs=None, digits=10):
        r"""`L^{(k)}(M_{\mathrm{std}}(C),s)`; by default `L^{\prime\prime}(M_{\mathrm{std}}(C),1)`.

        EXAMPLES::

            sage: from hyperell_regulator.lfunction import StandardLFunction
            sage: from hyperell_regulator.reduction import HyperellipticCurveWithReduction
            sage: R.<x> = ZZ[]
            sage: C = HyperellipticCurveWithReduction(x^5 + x^3 + 1)
        """
        L = self.lfun(num_coeffs, digits)
        K = _real_field(digits)
        with _pari_precision(digits):
            return K(_pari_lfun_call(
                pari.lfun, L, s, derivative, precision=_bits(digits)))


def standard_lvalue(C, bad_factors=None, derivative=2, num_coeffs=None,
                    digits=10, verbose=True, h1_conductor=None):
    r"""`L^{\prime\prime}(M_{\mathrm{std}}(C),1)` and its numerical data.

    Return the value, conductor, Euler factors and functional equation residual.
    ``h1_conductor`` supplies the exact Jacobian conductor when needed at 2.
    """
    L = StandardLFunction(C, bad_factors=bad_factors, h1_conductor=h1_conductor)
    out = {"curve": str(C), "bad_primes": C.bad_primes(),
           "euler_factors": {}, "conductor": None, "value": None,
           "derivative": derivative, "feq_check": None}
    for p in C.bad_primes():
        out["euler_factors"][p] = str(L.euler_factor(p))
    out["conductor"] = L.conductor()
    out["digits"] = digits
    # Build and cache the coefficients once for all evaluations below.
    L.lfun(num_coeffs, digits)
    out["feq_check"] = L.check_functional_equation(num_coeffs, digits)
    out["value"] = L.value(1, derivative, num_coeffs, digits)
    out["L_at_1"] = L.value(1, 0, num_coeffs, digits)
    out["Lprime_at_1"] = L.value(1, 1, num_coeffs, digits)
    out["vanishing_order"] = L.vanishing_order(1, num_coeffs=num_coeffs,
                                               digits=digits)
    if verbose:
        print(C)
        print("  bad primes %s, standard conductor %s"
              % (C.bad_primes(), out["conductor"]))
        for p, e in out["euler_factors"].items():
            print("    P_%s^std(T) = %s" % (p, e))
        print("  functional equation check (small is good): %.3e"
              % out["feq_check"])
        print("  L(std,1) = %.3e   L'(std,1) = %.3e   order of vanishing %s"
              % (out["L_at_1"], out["Lprime_at_1"], out["vanishing_order"]))
        print("  L^(%s)(std, 1) = %s   (%s digits)"
              % (derivative, out["value"], digits))
    return out


# ---------------------------------------------------------------------------
# L(C, s) itself
# ---------------------------------------------------------------------------

class CurveLFunction(SageObject):
    r"""The degree four `L`-function `L(C,s)=L(h^1(C),s)`.

    Its gamma factor is `\Gamma_{\mathbb C}(s)^2`;
    the functional equation is `s \leftrightarrow 2-s`.

    EXAMPLES::

        sage: from sage.all import ZZ, PolynomialRing
        sage: from hyperell_regulator import HyperellipticCurveWithReduction
        sage: from hyperell_regulator.lfunction import CurveLFunction
        sage: R = PolynomialRing(ZZ, 'x'); x = R.gen()
        sage: L = CurveLFunction(HyperellipticCurveWithReduction(x**5 + x**3 + 1))  # doctest: +ELLIPSIS
        ...Warning: unknown valuation of conductor at 2.
        sage: L.euler_factor(11).degree()
        4
        sage: L.conductor()
        3233
    """

    def __init__(self, C):
        self.C = C
        self._RT = PolynomialRing(QQ, 'T')
        self._cache = {}
        self._h1cache = {}
        RQ = PolynomialRing(QQ, 'x')
        self._pari_L = pari.lfungenus2([pari(str(RQ(C.f))), pari(str(RQ(C.h)))])
        self._eps = None

    def _repr_(self):
        return "L(C, s) for %s" % (self.C,)

    def precompute_h1(self, primes, nthreads=None):
        r"""Cache the Euler factors at good odd primes in one batch."""
        _precompute_h1(self, primes, nthreads)

    def euler_factor(self, p):
        r"""`P_p(T)` with `L(C,s)=\prod_p P_p(p^{-s})^{-1}`.

        Bad primes and 2 use PARI.

        EXAMPLES::

            sage: from sage.all import ZZ, PolynomialRing
            sage: from hyperell_regulator import HyperellipticCurveWithReduction
            sage: from hyperell_regulator.lfunction import CurveLFunction
            sage: R = PolynomialRing(ZZ, 'x'); x = R.gen()
            sage: L = CurveLFunction(HyperellipticCurveWithReduction(x**5 + x**3 + 1))
            ...Warning: unknown valuation of conductor at 2.
            sage: L.euler_factor(11)[0]
            1
        """
        p = ZZ(p)
        if p in self._h1cache:
            return self._h1cache[p]
        if p != 2 and p not in self.C.bad_primes():
            H = HyperellipticCurve(self.C.F.change_ring(GF(p)))
            out = self._RT(list(reversed(list(H.frobenius_polynomial()))))
            self._h1cache[p] = out
            return out
        e = self._pari_L.lfuneuler(p)
        num, den = pari(e).numerator(), pari(e).denominator()
        poly = den if str(num).strip() in ("1", "1.0") else den / num
        out = self._RT([QQ(c) for c in pari(poly).Vecrev()])
        self._h1cache[p] = out
        return out

    def conductor(self):
        r"""The conductor of `C`."""
        return ZZ(self.C.conductor())

    def dirichlet_coefficients(self, num):
        r"""The first ``num`` Dirichlet coefficients, with a zero sentinel at index 0."""
        return _dirichlet_coefficients(self, num)

    def root_number(self, num_coeffs=200, digits=10):
        r"""The sign `\pm1` minimizing the functional equation residual."""
        if self._eps is not None:
            return self._eps
        an = self.dirichlet_coefficients(num_coeffs)[1:]
        N = self.conductor()
        best, score = 1, None
        with _pari_precision(digits):
            for eps in (1, -1):
                L = _pari_lfun_call(
                    pari.lfuncreate, [[QQ(c) for c in an], 0, [0, 0, 1, 1],
                                     2, N, eps], coefficient_count=len(an), digits=digits)
                K = _real_field(digits)
                v = K(2) ** K(_pari_lfun_call(
                    pari.lfuncheckfeq, L, precision=_bits(digits)))
                if score is None or v < score:
                    best, score = eps, v
        self._eps = best
        return best

    def _coeffs_needed(self, digits, der=0):
        r"""PARI's coefficient count at the requested decimal precision."""
        return _coefficient_count(self.conductor(), [0, 0, 1, 1], 2, digits, der)

    def lfun(self, num_coeffs=None, digits=10, der=2):
        r"""The PARI object with gamma parameters `[0,0,1,1]` and weight 2."""
        if num_coeffs is None:
            num_coeffs = self._coeffs_needed(digits, der)
        key = ("lfun", int(num_coeffs), int(digits))
        if key in self._cache:
            return self._cache[key]
        an = self.dirichlet_coefficients(num_coeffs)[1:]
        eps = self.root_number(min(num_coeffs, 400), digits)
        with _pari_precision(digits):
            L = _pari_lfun_call(
                pari.lfuncreate, [[QQ(c) for c in an], 0, [0, 0, 1, 1],
                                 2, self.conductor(), eps],
                coefficient_count=len(an), digits=digits)
        self._cache[key] = L
        return L

    def check_functional_equation(self, num_coeffs=None, digits=10):
        r"""The functional equation residual ``2^lfuncheckfeq``."""
        L = self.lfun(num_coeffs, digits)
        K = _real_field(digits)
        with _pari_precision(digits):
            return K(2) ** K(_pari_lfun_call(
                pari.lfuncheckfeq, L, precision=_bits(digits)))

    def value(self, s=1, derivative=0, num_coeffs=None, digits=10):
        r"""
        `L^{(k)}(C, s)` at ``digits`` decimal digits.

        EXAMPLES::

            sage: from sage.all import ZZ, PolynomialRing
            sage: from hyperell_regulator import HyperellipticCurveWithReduction
            sage: from hyperell_regulator.lfunction import CurveLFunction
            sage: R = PolynomialRing(ZZ, 'x'); x = R.gen()
            sage: L = CurveLFunction(HyperellipticCurveWithReduction(x**5 + x**3 + 1))  # doctest: +ELLIPSIS
            ...Warning: unknown valuation of conductor at 2.
            sage: bool(L.value(2, digits=10) > 0)            # long time
            True
        """
        L = self.lfun(num_coeffs, digits)
        K = _real_field(digits)
        with _pari_precision(digits):
            return K(_pari_lfun_call(
                pari.lfun, L, s, derivative, precision=_bits(digits)))


def curve_lvalue(C, s=1, derivative=0, num_coeffs=None, digits=10,
                 verbose=False):
    r"""`L^{(k)}(C,s)` to ``digits`` decimal digits.

    Return the value, conductor, root number and functional equation residual.

    EXAMPLES::

        sage: from sage.all import ZZ, PolynomialRing
        sage: from hyperell_regulator import HyperellipticCurveWithReduction
        sage: from hyperell_regulator.lfunction import curve_lvalue
        sage: R = PolynomialRing(ZZ, 'x'); x = R.gen()
        sage: d = curve_lvalue(HyperellipticCurveWithReduction(x**5 + x**3 + 1),  # doctest: +ELLIPSIS
        ....:                  s=2, digits=12)                      # long time
        ...Warning: unknown valuation of conductor at 2.
        ...
        sage: d["conductor"]                                         # long time
        3233
    """
    L = CurveLFunction(C)
    out = {"curve": str(C), "conductor": L.conductor(),
           "root_number": L.root_number(),
           "s": s, "derivative": derivative, "digits": digits}
    out["feq_check"] = L.check_functional_equation(num_coeffs, digits)
    out["value"] = L.value(s, derivative, num_coeffs, digits)
    if verbose:
        print(C)
        print("  conductor %s, root number %+d"
              % (out["conductor"], out["root_number"]))
        print("  L^(%d)(C, %s) = %s" % (derivative, s, out["value"]))
        print("  functional equation check %s" % out["feq_check"])
    return out
