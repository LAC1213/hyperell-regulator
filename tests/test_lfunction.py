r"""Frobenius polynomials, Euler factors and numerical `L`-values."""

import os
import subprocess
from unittest.mock import Mock, patch

import pytest

from sage.all import (GF, HyperellipticCurve, PolynomialRing, QQ, RR, ZZ,
                      pari, prime_range)

from hyperell_regulator import HyperellipticCurveWithReduction
from hyperell_regulator.frobenius import (
    euler_binary,
    frobenius_polynomial,
    frobenius_polynomials,
    odd_model,
)
from hyperell_regulator.lfunction import (CurveLFunction, StandardLFunction,
                                          _pari_lfun_call, _pari_precision, _real_field)

R = PolynomialRing(ZZ, "x")
X = R.gen()
F = X ** 5 + X ** 3 + 1                      # conductor 3233, bad 53 and 61
PRIMES = [p for p in prime_range(5, 300) if p not in (53, 61)]


def _sage_frobenius(f, p):
    H = HyperellipticCurve(PolynomialRing(GF(p), "x")(f))
    return R(H.frobenius_polynomial())


def test_batched_frobenius_matches_sage():
    got = frobenius_polynomials(F, PRIMES, nthreads=4)
    assert set(got) == set(PRIMES)
    for p in PRIMES:
        assert got[p] == _sage_frobenius(F, p), p


def test_frobenius_is_a_weil_polynomial():
    r"""Check the degree, constant coefficient and functional equation."""
    for p, cp in frobenius_polynomials(F, PRIMES[:8]).items():
        c = cp.list()
        assert cp.degree() == 4 and c[4] == 1 and c[0] == p ** 2
        assert c[3] * p == c[1]


def test_falls_back_when_the_program_is_missing():
    r"""Sage fallback agrees with native Frobenius polynomials."""
    old = os.environ.get("HYPELLFROB_EULER")
    os.environ["HYPELLFROB_EULER"] = "/nonexistent/euler"
    try:
        assert euler_binary() is None
        got = frobenius_polynomials(F, PRIMES[:6])
        for p in PRIMES[:6]:
            assert got[p] == _sage_frobenius(F, p), p
    finally:
        if old is None:
            del os.environ["HYPELLFROB_EULER"]
        else:
            os.environ["HYPELLFROB_EULER"] = old


def test_odd_model():
    r"""The odd model is monic and integral when a rational Weierstrass point exists."""
    assert odd_model(F) == F                       # already monic and odd
    for f in ((X - 1) * (X ** 5 + X + 3), (X + 2) * (X ** 5 - X + 1)):
        g = odd_model(f)
        assert g.degree() == 5 and g.is_monic()
        assert all(c in ZZ for c in g.list())
    assert odd_model(X ** 6 + X + 1) is None
    assert odd_model(4 * (X ** 5 - X + 1) + (X ** 3 + X + 1) ** 2) is None


def test_even_degree_with_a_weierstrass_point_takes_the_fast_path():
    r"""The odd model gives the Frobenius polynomial of the original curve."""
    import time

    f = (X - 1) * (X ** 5 + X + 3)
    assert f.degree() % 2 == 0 and odd_model(f) is not None
    ps = [p for p in prime_range(7, 150) if f.discriminant() % p]
    t0 = time.time()
    got = frobenius_polynomials(f, ps, nthreads=4)
    assert time.time() - t0 < 30
    for p in ps:
        assert got[p] == _sage_frobenius(f, p), p


def test_even_degree_without_a_weierstrass_point():
    r"""Even degree Frobenius polynomials agree with Sage."""
    f = 4 * (X ** 5 - X + 1) + (X ** 3 + X + 1) ** 2
    assert odd_model(f) is None
    ps = [p for p in prime_range(7, 90) if f.discriminant() % p]
    got = frobenius_polynomials(f, ps, nthreads=4)
    for p in ps:
        assert got[p] == _sage_frobenius(f, p), p


def test_single_prime_convenience():
    assert frobenius_polynomial(F, 37) == _sage_frobenius(F, 37)


def test_curve_lfunction_matches_pari():
    r"""Curve Dirichlet coefficients and values agree with PARI."""
    C = HyperellipticCurveWithReduction(F)
    L = CurveLFunction(C)
    assert L.conductor() == 3233
    assert L.root_number() in (1, -1)
    ours = [ZZ(c) for c in L.dirichlet_coefficients(40)[1:]]
    pari.set_real_precision(20)
    theirs = [ZZ(c) for c in pari.lfunan(L._pari_L, 40)]
    assert ours == theirs
    params = list(pari.lfunparams(L._pari_L))
    assert ZZ(params[0]) == L.conductor() and ZZ(params[1]) == 2


def test_curve_lvalue_matches_pari():
    r"""Compare curve values with PARI at equal precision."""
    C = HyperellipticCurveWithReduction(F)
    L = CurveLFunction(C)
    # the Dirichlet coefficients are the part that is ours, and they agree
    # with PARI's exactly
    ours = L.dirichlet_coefficients(40)
    theirs = pari.lfunan(L._pari_L, 40)
    assert [QQ(ours[n]) for n in range(1, 41)] == [QQ(v) for v in theirs]
    # the values agree only as well as the underlying object allows: PARI's
    # own lfungenus2 object for this curve satisfies its functional equation
    # to 2^-7, and ours -- built from the same coefficients and the same
    # conductor -- to the same, so a few digits is what either is worth
    bits = 128
    K = _real_field(30)
    feq = K(2) ** K(pari.lfuncheckfeq(L._pari_L, precision=bits))
    assert feq > 1e-4, feq            # PARI's own object, not ours, says so
    v_ours = L.value(2, 0, digits=int(bits * 0.301))
    v_pari = pari.lfun(L._pari_L, 2, precision=bits)
    assert abs(RR(v_ours) - RR(v_pari)) < 1e-3


def test_arbitrary_precision_really_is():
    r"""Values and functional equation residuals improve with precision."""
    # LMFDB 249.a.249.1, whose bad primes 3 and 83 are both odd, so the
    # cluster pictures give its standard L-function outright
    R = PolynomialRing(ZZ, "x")
    x = R.gen()
    L = StandardLFunction(
        HyperellipticCurveWithReduction((x ** 3 + 1) ** 2 + 4 * (x ** 2 + x)))
    lo, hi = 15, 30
    nlo, nhi = L._coeffs_needed(lo), L._coeffs_needed(hi)
    vlo, vhi = L.value(1, 2, nlo, lo), L.value(1, 2, nhi, hi)
    assert vlo.parent().precision() < vhi.parent().precision()
    assert vhi.parent().precision() >= 100
    flo = L.check_functional_equation(nlo, lo)
    fhi = L.check_functional_equation(nhi, hi)
    assert flo < 1e-20, flo              # not the 3.6e-15 of the fixed value
    assert fhi < flo, (flo, fhi)         # and it improves with precision
    assert abs(RR(vlo) - RR(vhi)) < 1e-14


def test_real_field_is_wide_enough():
    for d in (10, 30, 60):
        assert _real_field(d).precision() >= int(d * 3.32)


def test_large_coefficient_budget_can_exceed_one_gib_without_shrinking_stack():
    """The coefficient budget may require a PARI stack larger than 1 GiB."""
    engine = Mock(wraps=pari)
    engine.allocatemem = Mock()
    engine.stacksize.return_value = 8_000_000
    engine.stacksizemax.return_value = 2**30
    with patch('hyperell_regulator.lfunction.pari', engine):
        assert _pari_lfun_call(lambda: 42, coefficient_count=2_329_751, digits=100) == 42
    size, ceiling = engine.allocatemem.call_args.args
    assert size == 8_000_000 and ceiling > 2**30
    engine.allocatemem.reset_mock()
    engine.stacksize.return_value = 2**30
    engine.stacksizemax.return_value = 2**34
    with patch('hyperell_regulator.lfunction.pari', engine):
        _pari_lfun_call(lambda: 42, coefficient_count=2_329_751, digits=100)
    engine.allocatemem.assert_not_called()


def test_pari_stack_exhaustion_retries_the_operation():
    """PARI stack exhaustion triggers repeated enlargement."""
    previous = pari.stacksize(), pari.stacksizemax()
    attempts = 0

    def build_vector():
        nonlocal attempts
        attempts += 1
        return pari('vector(100000,i,i)')

    try:
        pari.allocatemem(2**18, 2**19, silent=True)
        vector = _pari_lfun_call(build_vector)
        assert len(vector) == 100_000 and vector[99_999] == 100_000
        assert attempts > 1
    finally:
        pari.allocatemem(*previous, silent=True)


def test_other_pari_errors_propagate_without_retry():
    from cypari2 import PariError

    attempts = 0
    previous = pari.stacksize(), pari.stacksizemax(), pari.get_real_precision()

    def fail():
        nonlocal attempts
        attempts += 1
        return pari('1/0')

    with pytest.raises(PariError, match='impossible inverse'):
        with _pari_precision(100):
            _pari_lfun_call(fail)
    assert attempts == 1
    assert (pari.stacksize(), pari.stacksizemax(), pari.get_real_precision()) == previous


@pytest.mark.parametrize('cls', [CurveLFunction, StandardLFunction])
def test_coefficient_expansion_matches_divisor_sums(cls):
    r"""The coefficients of `\zeta(s)\zeta(s-1)` are `\sigma_1(n)`."""
    L = cls(HyperellipticCurveWithReduction(F))
    T = PolynomialRing(QQ, 'T').gen()
    with patch.object(cls, 'precompute_h1'), patch.object(
            cls, 'euler_factor', side_effect=lambda p: (1-T)*(1-p*T)):
        coefficients = L.dirichlet_coefficients(40)
    if cls is CurveLFunction:
        assert coefficients[0] == 0
        coefficients = coefficients[1:]
    assert coefficients == [sum(ZZ(n).divisors()) for n in range(1, 41)]


def test_pari_precision_restored_after_nested_failure():
    previous = pari.get_real_precision()
    with _pari_precision(30):
        outer = pari.get_real_precision()
        with pytest.raises(RuntimeError, match='failed'):
            with _pari_precision(60):
                assert pari.get_real_precision() > outer
                raise RuntimeError('failed')
        assert pari.get_real_precision() == outer
    assert pari.get_real_precision() == previous


@pytest.mark.parametrize('shift', [0, 1])
def test_conductor_388_factor_at_two_is_detected_from_reduction(shift):
    x = PolynomialRing(QQ, 'x').gen()
    f, h = x+2*x**2-x**4, x**3+x+1
    L = StandardLFunction(HyperellipticCurveWithReduction(f(x+shift), h(x+shift)))
    T = L._RT.gen()
    assert L.h1_euler_factor(2) == 1+T+T**2
    assert L.wedge2_euler_factor(2) == (1-T)*(1-8*T**3)
    assert L.euler_factor(2) == (1-T)*(1+2*T+4*T**2)
    assert L.conductor() == 4*97**2
    # An explicit wedge-square override still takes precedence.
    L.bad_factors[2] = (1-2*T)*(1+T)
    assert L.wedge2_euler_factor(2) == (1-2*T)*(1+T)


def test_unsupported_reduction_at_two_still_requires_override():
    L = StandardLFunction(HyperellipticCurveWithReduction(F))
    with pytest.raises(NotImplementedError, match='additional inertia data'):
        L.wedge2_euler_factor(2)


@pytest.mark.parametrize('f,h,conductor,local_factor,standard_conductor', [
    (-X**5, 1+X, 394, [1, 0, -2, -4], 4*197**2),
    (3-4*X-4*X**2+3*X**3+X**4+X**5, X+X**2, 826,
     [1, -1, 0, -4], 4*413**2),
    (-1-X+5*X**2-4*X**3-X**4+X**5, 1+X**3, 644,
     [1, 3, 0, -4], 4*161**2),
    (-2*X**2-4*X**3+3*X**4+8*X**5, 1+X, 464,
     [1, -2], 256*29**2),
    (-2*X**2-X**3+4*X**4+4*X**5, 1, 704,
     [1, -2], 4096*11**2),
    (X-X**3-X**4, X+X**3, 856, [1, -2], 16*107**2),
])
def test_two_adic_factors_and_wild_conductors_satisfy_functional_equation(
        f, h, conductor, local_factor, standard_conductor):
    """The factors at 2 and wild conductors satisfy the functional equation."""
    L = StandardLFunction(HyperellipticCurveWithReduction(f, h), h1_conductor=conductor)
    assert L.euler_factor(2) == L._RT(local_factor)
    assert L.conductor() == standard_conductor
    assert L.check_functional_equation(digits=10) < 1e-15


def test_wild_conductor_is_not_silently_assumed_tame_when_pari_omits_it():
    C = HyperellipticCurveWithReduction(-2*X**2-4*X**3+3*X**4+8*X**5, 1+X)
    L = StandardLFunction(C)
    T = L._RT.gen()
    assert L.euler_factor(2) == 1-2*T
    with pytest.raises(NotImplementedError, match='pass h1_conductor'):
        L.conductor()
    inconsistent = StandardLFunction(C, h1_conductor=29)
    with pytest.raises(ArithmeticError, match='smaller than its tame part'):
        inconsistent.conductor()


def test_1116_toric_homology_uses_submodule_coordinates():
    C = HyperellipticCurveWithReduction(-X+X**2+2*X**3+X**4, 1+X**3)
    L = StandardLFunction(C)
    T = L._RT.gen()
    # The totally toric h^1 factor is independently supplied by PARI.
    assert L.toric_factor(3) == L.h1_euler_factor(3) == 1+T+T**2
    assert L.wedge2_euler_factor(3) == (1-T)*(1-27*T**3)
    assert L.euler_factor(3) == (1-T)*(1+3*T+9*T**2)
    assert L.conductor() == 4*9*31**2


@pytest.mark.parametrize('support_order', [(0, 1, 2), (2, 1, 0)])
def test_toric_frobenius_uses_echelonized_cycle_basis(support_order):
    """Order-three toric Frobenius is independent of the cycle basis."""
    C = HyperellipticCurveWithReduction(-X+X**2+2*X**3+X**4, 1+X**3)
    L = StandardLFunction(C)
    R = C.cluster_picture(3)
    H1, _, _ = R.homology_of_special_fibre()
    ambient = H1.ambient()
    edges = list(ambient.basis())
    generators = [edges[1]-edges[0], edges[2]-edges[0]]
    cycles = ambient.submodule(
        generators, support_order=[edges[i].support()[0] for i in support_order])
    # Reproduce the upstream map: its generator indices do not necessarily
    # match the echelonized basis returned by Sage's submodule constructor.
    def on_basis(i):
        return ambient.sum_of_terms(
            (s.frobenius(), n*s.epsilon(
                lambda c: c.frobenius(), R.field_frobenius()))
            for s, n in generators[i])
    frob = cycles.module_morphism(on_basis=on_basis, codomain=cycles)
    with patch.object(type(R), 'homology_of_special_fibre',
                      return_value=(cycles, None, frob)):
        assert L.toric_factor(3) == L.h1_euler_factor(3)


@pytest.mark.parametrize('entries', [
    [1, 0, 0, 1], [-1, 0, 0, -1], [1, 0, 0, -1],
    [0, -1, 1, -1], [0, -1, 1, 0], [0, -1, 1, 1],
])
def test_toric_factor_matches_exterior_square_monodromy_kernel(entries):
    """Compare toric factors with the exterior-square monodromy kernel."""
    from itertools import combinations
    from sage.all import Matrix, block_diagonal_matrix, block_matrix, identity_matrix, zero_matrix

    M = Matrix(QQ, 2, 2, entries)
    p = 7
    frobenius = block_diagonal_matrix(M, p*M)
    N = block_matrix([[zero_matrix(QQ, 2), identity_matrix(QQ, 2)],
                      [zero_matrix(QQ, 2), zero_matrix(QQ, 2)]])
    pairs = list(combinations(range(4), 2))
    Fw = Matrix(QQ, [[frobenius[a, i]*frobenius[b, j]-frobenius[b, i]*frobenius[a, j]
                     for i, j in pairs] for a, b in pairs])
    Nw = Matrix(QQ, [[N[a, i]*(b == j)-N[b, i]*(a == j)
                     +(a == i)*N[b, j]-(b == i)*N[a, j]
                     for i, j in pairs] for a, b in pairs])
    invariant = Nw.right_kernel()
    restricted = Matrix(QQ, [invariant.coordinate_vector(Fw*v) for v in invariant.basis()])
    L = StandardLFunction(HyperellipticCurveWithReduction(F))
    U = L._RT(M.charpoly().reverse())
    expected = L._RT(restricted.charpoly().reverse())
    assert L._wedge2_totally_toric(p, U) == expected


@pytest.mark.parametrize("workers", [2, 4, 8, None])
@pytest.mark.parametrize("f", [F, 4 * (X**5 - X + 1) + (X**3 + X + 1)**2])
def test_native_parallel_matrices_match_serial_and_sage(workers, f):
    """Native matrices agree across worker counts and with Sage."""
    binary = euler_binary()
    if binary is None:
        pytest.skip("Frobenius accelerator not built")
    from hyperell_regulator.frobenius import _lift
    ps = [int(p) for p in prime_range(41, 60) if f.discriminant() % p]
    batch = " ".join(map(str, f.list())) + "\n"
    batch += "".join("%d 2\n" % p for p in ps)
    def run(n):
        command = [binary] if n is None else [binary, str(n)]
        result = subprocess.run(command, input=batch, text=True,
                                capture_output=True, check=True, timeout=30)
        rows = [list(map(int, line.split())) for line in result.stdout.splitlines()]
        assert len(rows) == len(ps)
        for p, row in zip(ps, rows):
            assert row[:3] == [p, 1, 4]
            assert len(row) == 19
            matrix = [row[3 + 4*i:7 + 4*i] for i in range(4)]
            assert _lift(matrix, p, 2, 2) == _sage_frobenius(f, p)
        return rows
    serial = run(1)
    for _ in range(3):
        assert run(workers) == serial
