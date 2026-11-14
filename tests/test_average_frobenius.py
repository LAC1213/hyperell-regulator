"""Simultaneous Frobenius against independent Sage matrices and Euler factors.

Native calls never use fallback. Tests require matrices for supported jobs,
including even degree, and exercise the Sage fallback for unsupported jobs.
"""
import random
import os
import subprocess

import pytest
from sage.all import GF, HyperellipticCurve, PolynomialRing, ZZ, matrix, prime_range
from sage.schemes.hyperelliptic_curves.hypellfrob import hypellfrob

from hyperell_regulator.frobenius import (
    _bound, _lift, euler_binary, frobenius_polynomials, odd_model,
)

R = PolynomialRing(ZZ, "x")
x = R.gen()


@pytest.fixture
def binary():
    result = euler_binary()
    if result is None:
        pytest.skip("native Frobenius accelerator not built")
    return result


def native(binary, q, jobs, backend="rforest-only", workers=1):
    batch = " ".join(map(str, q.list())) + "\n"
    batch += "".join("%d %d\n" % job for job in jobs)
    result = subprocess.run([binary, str(workers), backend], input=batch,
                            text=True, capture_output=True, check=True,
                            timeout=90)
    rows = [list(map(int, line.split())) for line in result.stdout.splitlines()]
    assert len(rows) == len(jobs)
    assert [row[0] for row in rows] == [p for p, _ in jobs]
    return rows


def sage_matrices(q, jobs):
    """Independent full matrices for odd degree, in the same de Rham basis."""
    rows = []
    for p, n in jobs:
        m = matrix(ZZ, hypellfrob(p, n, q))
        rows.append([p, 1, m.nrows()] + [int(v) % p**n for v in m.list()])
    return rows


def check_characteristic_polynomials(q, jobs, rows):
    """Compare every coefficient modulo p^N, without imposing reciprocity."""
    genus = (q.degree()-1)//2
    for (p, n), row in zip(jobs, rows):
        assert row[:3] == [p, 1, 2*genus]
        m = matrix(ZZ, 2*genus, 2*genus, row[3:])
        expected = HyperellipticCurve(q.change_ring(GF(p))).frobenius_polynomial()
        assert all((a-b) % p**n == 0 for a, b in zip(m.charpoly(), expected))


def sage_polynomials(q, primes):
    return {p: R(HyperellipticCurve(q.change_ring(GF(p))).frobenius_polynomial())
            for p in primes}


def test_native_default_uses_rforest_and_declines_unsupported_jobs(binary):
    batch = "1 6 9 2 -2 0 1\n41 2\n43 2\n31 2\n"
    env = dict(os.environ, HYPELLFROB_PROFILE="1")
    default = subprocess.run([binary, "4"], input=batch, text=True,
                             capture_output=True, check=True, timeout=30, env=env)
    explicit = subprocess.run([binary, "4", "rforest-only"], input=batch, text=True,
                              capture_output=True, check=True, timeout=30, env=env)
    assert "rforest a=" in default.stderr
    assert default.stdout == explicit.stdout
    assert [line.split()[1] for line in default.stdout.splitlines()] == ["1", "1", "0"]


def test_native_rejects_removed_backend(binary):
    result = subprocess.run([binary, "4", "legacy"], input="1 0 0 1 0 1\n101 2\n",
                            text=True, capture_output=True)
    assert result.returncode == 2
    assert "backend must be rforest" in result.stderr


def test_python_default_monicizes_even_models(binary, monkeypatch):
    from hyperell_regulator import frobenius
    monkeypatch.setenv("HYPELLFROB_EULER", binary)
    monkeypatch.delenv("HYPELLFROB_BACKEND", raising=False)
    real_run = subprocess.run
    calls = []

    def capture(command, **kwargs):
        if command and command[0] == binary:
            calls.append((command, kwargs["input"].splitlines()[0]))
        return real_run(command, **kwargs)

    monkeypatch.setattr(frobenius.subprocess, "run", capture)
    q = 2*(x**6-2*x**4+2*x**3+9*x**2+6*x+1)
    ps = [41, 43]
    expected = sage_polynomials(q, ps)
    assert frobenius_polynomials(q, ps, nthreads=4) == expected
    assert calls[-1][0] == [binary, "4", "rforest"]
    assert calls[-1][1].split()[-1] == "1"
    monkeypatch.setenv("HYPELLFROB_BACKEND", "rforest")
    assert frobenius_polynomials(q, ps, nthreads=4) == expected


@pytest.mark.parametrize("backend", ["legacy", "unknown"])
def test_python_rejects_removed_or_unknown_backend(monkeypatch, backend):
    monkeypatch.setenv("HYPELLFROB_EULER", "/nonexistent/euler")
    with pytest.raises(ValueError, match="backend"):
        frobenius_polynomials(x**5+x+1, [101], backend=backend)
    monkeypatch.setenv("HYPELLFROB_BACKEND", backend)
    with pytest.raises(ValueError, match="backend"):
        frobenius_polynomials(x**5+x+1, [101])


@pytest.mark.parametrize("q", [
    x**3 - 2*x + 3,
    x**5 + x**3 + 1,
    x**5 + x**2 + x,                 # c0=0, so no vertical reductions
    x**5 - 7*x**4 + 3*x**2 - 5*x - 2,
    x**7 - x + 1,
])
@pytest.mark.parametrize("precision", [1, 2, 3, 4])
@pytest.mark.parametrize("backend", ["rforest-only"])
def test_full_matrices(binary, q, precision, backend):
    jobs = [(int(p), precision) for p in prime_range(151, 200)
            if q.discriminant() % p]
    expected = sage_matrices(q, jobs)
    actual = native(binary, q, jobs, backend)
    assert all(row[1] == 1 for row in actual), actual
    assert actual == expected


@pytest.mark.parametrize("backend", ["rforest-only"])
def test_unsorted_duplicate_jobs_mixed_precision_and_threads(binary, backend):
    q = x**5 + x**3 + 1
    jobs = [(211, 2), (101, 3), (211, 2), (151, 1), (211, 3), (103, 3)]
    expected = sage_matrices(q, jobs)
    for workers in (1, 4, 16):
        assert native(binary, q, jobs, backend, workers=workers) == expected


@pytest.mark.parametrize("backend", ["rforest-only"])
def test_dense_batch_and_deterministic_random_curves(binary, backend):
    rng = random.Random(12108239)
    for genus in (1, 2, 3):
        for _ in range(2):
            q = x**(2*genus+1) + R([rng.randrange(-5, 6)
                                     for _ in range(2*genus+1)])
            if not q.is_squarefree():
                continue
            jobs = [(int(p), 2) for p in prime_range(101, 500)
                    if q.discriminant() % p]
            actual = native(binary, q, jobs, backend)
            assert all(row[1] == 1 for row in actual)
            assert actual == sage_matrices(q, jobs)


@pytest.mark.parametrize("q,jobs", [
    (x**5 + x**3 + 1, [(3, 2), (31, 2), (53, 2), (61, 2), (101, 2)]),
    (x**5 + x + 101, [(101, 2), (103, 2)]),  # c0 divisible by a good prime
    (x**6 + x + 1, [(31, 2), (101, 2)]),
    (2*x**5 + x + 1, [(101, 2)]),
    ((x+1)**2*(x**3+2), [(101, 2)]),
])
@pytest.mark.parametrize("backend", ["rforest-only"])
def test_declines_without_native_fallback(binary, q, jobs, backend):
    strict = native(binary, q, jobs, backend)
    # Each case deliberately includes at least one unsupported job.
    assert any(row[1] == 0 for row in strict)
    if q.is_monic() and q.is_squarefree():
        assert native(binary, q, jobs, "rforest") == strict


@pytest.mark.parametrize("backend", ["rforest-only"])
def test_strict_input_validation(binary, backend):
    q = x**5 + x**3 + 1
    jobs = [(101, 0), (101, -1), (1, 2), (2, 2), (121, 2)]
    assert all(row[1] == 0 for row in native(binary, q, jobs, backend))
    assert native(binary, q, [], backend) == []


@pytest.mark.parametrize("q", [x**3-x+1, x**5+x**3+1, x**7-x+1])
@pytest.mark.parametrize("backend", ["rforest-only"])
def test_euler_factors_against_sage(binary, q, backend):
    genus = (q.degree()-1)//2
    ps = [int(p) for p in prime_range(151, 230) if q.discriminant() % p]
    jobs = [(p, _bound(p, genus)) for p in ps]
    for (p, n), row in zip(jobs, native(binary, q, jobs, backend)):
        assert row[:3] == [p, 1, 2*genus]
        dim = row[2]
        matrix = [row[3+i*dim:3+(i+1)*dim] for i in range(dim)]
        expected = HyperellipticCurve(PolynomialRing(GF(p), "x")(q)).frobenius_polynomial()
        assert _lift(matrix, p, genus, n) == expected


@pytest.mark.parametrize("backend", ["rforest"])
def test_python_api_models_and_environment(binary, monkeypatch, backend):
    monkeypatch.setenv("HYPELLFROB_EULER", binary)
    ps = [101, 103, 107, 109]
    for q in (x**5+x**3+1, (x-1)*(x**5+x+3), x**6+x+1):
        good = [p for p in ps if q.discriminant() % p and
                (odd_model(q) is None or odd_model(q).discriminant() % p)]
        expected = sage_polynomials(q, good)
        assert frobenius_polynomials(q, good, backend=backend) == expected
        monkeypatch.setenv("HYPELLFROB_BACKEND", backend)
        assert frobenius_polynomials(q, good) == expected
    with pytest.raises(ValueError, match="backend"):
        frobenius_polynomials(x**5+x+1, ps, backend="unknown")


@pytest.mark.parametrize("backend", ["rforest"])
def test_lfunction_coefficients_use_new_backend(binary, monkeypatch, backend):
    from hyperell_regulator import HyperellipticCurveWithReduction
    from hyperell_regulator.lfunction import CurveLFunction, StandardLFunction
    monkeypatch.setenv("HYPELLFROB_EULER", binary)
    # Conductor 249, with only odd bad primes; the standard L-function has
    # no automatic bad factor at 2 for the other test curve.
    q = (x**3+1)**2 + 4*(x**2+x)
    for cls in (CurveLFunction, StandardLFunction):
        monkeypatch.setenv("HYPELLFROB_EULER", "/nonexistent/euler")
        expected = cls(HyperellipticCurveWithReduction(q)).dirichlet_coefficients(250)
        monkeypatch.setenv("HYPELLFROB_EULER", binary)
        monkeypatch.setenv("HYPELLFROB_BACKEND", backend)
        actual = cls(HyperellipticCurveWithReduction(q)).dirichlet_coefficients(250)
        assert actual == expected


def test_conductor_388_uses_even_degree_forest(binary, monkeypatch):
    """The even forest handles this curve without a rational branch point."""
    from sage.all import QQ, pari
    from hyperell_regulator import HyperellipticCurveWithReduction
    from hyperell_regulator.lfunction import CurveLFunction, StandardLFunction

    f, h = -x**4+2*x**2+x, x**3+x+1
    q = 4*f+h**2
    assert q.roots(QQ) == []
    assert odd_model(q) is None
    assert q.change_ring(GF(11)).is_irreducible()  # no root over F_(11^2)
    jobs = [(p, 2) for p in (41, 43, 47, 71, 101)]
    strict = native(binary, q, jobs, "rforest-only", workers=4)
    assert all(row[:3] == [p, 1, 4] for (p, _), row in zip(jobs, strict))
    actual = native(binary, q, jobs, "rforest", workers=4)
    assert actual == strict
    check_characteristic_polynomials(q, jobs, actual)
    reference = pari.lfungenus2([pari(str(f)), pari(str(h))])
    for (p, n), row in zip(jobs, actual):
        assert row[:3] == [p, 1, 4]
        matrix = [row[3+4*i:7+4*i] for i in range(4)]
        factor = R(reference.lfuneuler(p).denominator())
        assert _lift(matrix, p, 2, n) == factor.reverse()

    monkeypatch.setenv("HYPELLFROB_EULER", binary)
    for cls in (CurveLFunction, StandardLFunction):
        monkeypatch.setenv("HYPELLFROB_EULER", "/nonexistent/euler")
        expected = cls(HyperellipticCurveWithReduction(f, h)).dirichlet_coefficients(128)
        monkeypatch.setenv("HYPELLFROB_EULER", binary)
        monkeypatch.setenv("HYPELLFROB_BACKEND", "rforest")
        actual = cls(HyperellipticCurveWithReduction(f, h)).dirichlet_coefficients(128)
        assert actual == expected


@pytest.mark.parametrize("q", [
    x**4-x+1,
    x**6-2*x**4+2*x**3+9*x**2+6*x+1,
    x**6+x+1,
    x**6+x**2+x,
    x**6-3*x**5+2*x**4-5*x**2+7*x-2,
    x**8-x+1,
])
@pytest.mark.parametrize("precision", [1, 2, 3, 4])
def test_even_degree_characteristic_polynomials_at_full_precision(binary, q, precision):
    # One good prime beyond the common p-admissibility bound. Compare all
    # characteristic coefficients, including higher p-adic digits.
    bound = q.degree()*(4*precision-1)
    p = next(p for p in prime_range(bound+1, bound+100) if q.discriminant()%p)
    jobs = [(int(p), precision)]
    actual = native(binary, q, jobs, "rforest-only", workers=4)
    assert actual[0][:3] == [p, 1, q.degree()-2]
    check_characteristic_polynomials(q, jobs, actual)


@pytest.mark.parametrize("q", [x**4-x+1, x**6+x+1, x**6+x**2+x, x**8-x+1])
def test_even_degree_dense_euler_factors(binary, q):
    genus = (q.degree()-2)//2
    jobs = [(int(p), _bound(p, genus)) for p in prime_range(151, 210)
            if q.discriminant()%p]
    actual = native(binary, q, jobs, "rforest-only", workers=4)
    for (p, n), row in zip(jobs, actual):
        assert row[:3] == [p, 1, 2*genus]
        matrix = [row[3+i*2*genus:3+(i+1)*2*genus] for i in range(2*genus)]
        expected = HyperellipticCurve(q.change_ring(GF(p))).frobenius_polynomial()
        assert _lift(matrix, p, genus, n) == expected


def test_even_degree_mixed_jobs_and_threads(binary):
    q = x**6-2*x**4+2*x**3+9*x**2+6*x+1
    jobs = [(101, 2), (83, 3), (101, 2), (71, 1), (101, 3), (89, 3)]
    expected = native(binary, q, jobs, workers=1)
    check_characteristic_polynomials(q, jobs, expected)
    assert all(row[1] for row in expected)
    for workers in (1, 4, 16):
        assert native(binary, q, jobs, "rforest-only", workers=workers) == expected


@pytest.mark.parametrize("q", [
    2*(x**6+x+1),
    4*(x**6+x+1),
    -3*(x**6+x+1),
    -2*x**6+x**3+x+1,
    5*x**6+x**5+x+1,              # p=5 remains good after degree drops
    3*x**4+x+1,
])
def test_even_degree_nonmonic_python_models(binary, monkeypatch, q):
    monkeypatch.setenv("HYPELLFROB_EULER", binary)
    ps = [p for p in (3, 5, 7, 43, 47, 59, 71, 83) if q.discriminant()%p]
    actual = frobenius_polynomials(q, ps, nthreads=4, backend="rforest")
    for p in ps:
        expected = HyperellipticCurve(q.change_ring(GF(p))).frobenius_polynomial()
        assert actual[p] == expected


@pytest.mark.parametrize("kappa", [0, 2, 20])
@pytest.mark.parametrize("evaluator", ["differences", "powers"])
def test_even_degree_prefix_and_subtrees(binary, monkeypatch, kappa, evaluator):
    q = x**6-2*x**4+2*x**3+9*x**2+6*x+1
    jobs = [(int(p), 2) for p in prime_range(41, 170) if q.discriminant()%p]
    jobs += [(101, 2), (163, 2)]
    random.Random(388).shuffle(jobs)
    expected = native(binary, q, jobs, "rforest-only", workers=4)
    assert all(row[1] for row in expected)
    monkeypatch.setenv("HYPELLFROB_RFOREST_EVAL", evaluator)
    batch = " ".join(map(str, q.list()))+"\n"+"".join("%d %d\n"%job for job in jobs)
    actual = subprocess.run([binary, "4", "rforest-only", str(kappa)],
                            input=batch, text=True, capture_output=True,
                            check=True, timeout=90)
    assert [list(map(int, line.split())) for line in actual.stdout.splitlines()] == expected


@pytest.mark.parametrize("kappa", [0, 1, 3, 20])
def test_rforest_subtree_boundaries(binary, kappa):
    q = x**5+x**3+1
    jobs = [(int(p), 2) for p in prime_range(37, 230) if q.discriminant() % p]
    jobs += [(101, 2), (211, 2)]
    random.Random(10).shuffle(jobs)
    expected = sage_matrices(q, jobs)
    batch = " ".join(map(str, q.list())) + "\n" + "".join("%d %d\n" % job for job in jobs)
    result = subprocess.run([binary, "4", "rforest-only", str(kappa)], input=batch,
                            capture_output=True, text=True, check=True, timeout=90)
    assert [list(map(int, line.split())) for line in result.stdout.splitlines()] == expected


@pytest.mark.parametrize("evaluator", ["differences", "powers"])
@pytest.mark.parametrize("fft", ["0", "1"])
@pytest.mark.parametrize("kappa", [1, 20])
def test_rforest_evaluators_and_fft_controls(binary, monkeypatch, evaluator, fft, kappa):
    monkeypatch.setenv("HYPELLFROB_RFOREST_EVAL", evaluator)
    monkeypatch.setenv("HYPELLFROB_RFOREST_FFT", fft)
    q = x**5+x**3+1
    jobs = [(37, 2), (101, 2), (101, 2), (163, 3), (211, 2)]
    expected = sage_matrices(q, jobs)
    batch = " ".join(map(str, q.list())) + "\n" + "".join("%d %d\n" % job for job in jobs)
    result = subprocess.run([binary, "4", "rforest-only", str(kappa)], input=batch,
                            capture_output=True, text=True, check=True, timeout=90)
    assert [list(map(int, line.split())) for line in result.stdout.splitlines()] == expected
