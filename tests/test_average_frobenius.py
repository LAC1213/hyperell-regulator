"""Harvey 2013 against full 2007 matrices and independent Sage Euler factors.

The strict native mode is essential: fallback must not hide an arithmetic bug
or an implementation that simply declines all the new backend's work.
"""
import random
import subprocess

import pytest
from sage.all import GF, HyperellipticCurve, PolynomialRing, ZZ, prime_range

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


def native(binary, q, jobs, backend="remainder-only", workers=1):
    batch = " ".join(map(str, q.list())) + "\n"
    batch += "".join("%d %d\n" % job for job in jobs)
    result = subprocess.run([binary, str(workers), backend], input=batch,
                            text=True, capture_output=True, check=True,
                            timeout=90)
    rows = [list(map(int, line.split())) for line in result.stdout.splitlines()]
    assert len(rows) == len(jobs)
    assert [row[0] for row in rows] == [p for p, _ in jobs]
    return rows


@pytest.mark.parametrize("q", [
    x**3 - 2*x + 3,
    x**5 + x**3 + 1,
    x**5 + x**2 + x,                 # c0=0, so no vertical reductions
    x**5 - 7*x**4 + 3*x**2 - 5*x - 2,
    x**7 - x + 1,
])
@pytest.mark.parametrize("precision", [1, 2, 3, 4])
def test_full_matrices(binary, q, precision):
    jobs = [(int(p), precision) for p in prime_range(151, 200)
            if q.discriminant() % p]
    expected = native(binary, q, jobs, "legacy")
    actual = native(binary, q, jobs)
    assert all(row[1] == 1 for row in actual), actual
    assert actual == expected


def test_unsorted_duplicate_jobs_mixed_precision_and_threads(binary):
    q = x**5 + x**3 + 1
    jobs = [(211, 2), (101, 3), (211, 2), (151, 1), (211, 3), (103, 3)]
    expected = native(binary, q, jobs, "legacy")
    for workers in (1, 4, 16):
        assert native(binary, q, jobs, workers=workers) == expected


def test_dense_batch_and_deterministic_random_curves(binary):
    rng = random.Random(12108239)
    for genus in (1, 2, 3):
        for _ in range(2):
            q = x**(2*genus+1) + R([rng.randrange(-5, 6)
                                     for _ in range(2*genus+1)])
            if not q.is_squarefree():
                continue
            jobs = [(int(p), 2) for p in prime_range(101, 500)
                    if q.discriminant() % p]
            actual = native(binary, q, jobs)
            assert all(row[1] == 1 for row in actual)
            assert actual == native(binary, q, jobs, "legacy")


@pytest.mark.parametrize("q,jobs", [
    (x**5 + x**3 + 1, [(3, 2), (31, 2), (53, 2), (61, 2), (101, 2)]),
    (x**5 + x + 101, [(101, 2), (103, 2)]),  # c0 divisible by a good prime
    (x**6 + x + 1, [(101, 2)]),
    (2*x**5 + x + 1, [(101, 2)]),
    ((x+1)**2*(x**3+2), [(101, 2)]),
])
def test_declines_and_native_fallback(binary, q, jobs):
    strict = native(binary, q, jobs)
    # Each case deliberately includes at least one unsupported job.
    assert any(row[1] == 0 for row in strict)
    if q.is_monic() and q.is_squarefree():
        assert native(binary, q, jobs, "remainder") == native(binary, q, jobs, "legacy")


def test_strict_input_validation(binary):
    q = x**5 + x**3 + 1
    jobs = [(101, 0), (101, -1), (1, 2), (2, 2), (121, 2)]
    assert all(row[1] == 0 for row in native(binary, q, jobs))
    assert native(binary, q, []) == []


@pytest.mark.parametrize("q", [x**3-x+1, x**5+x**3+1, x**7-x+1])
def test_euler_factors_against_sage(binary, q):
    genus = (q.degree()-1)//2
    ps = [int(p) for p in prime_range(151, 230) if q.discriminant() % p]
    jobs = [(p, _bound(p, genus)) for p in ps]
    for (p, n), row in zip(jobs, native(binary, q, jobs)):
        assert row[:3] == [p, 1, 2*genus]
        dim = row[2]
        matrix = [row[3+i*dim:3+(i+1)*dim] for i in range(dim)]
        expected = HyperellipticCurve(PolynomialRing(GF(p), "x")(q)).frobenius_polynomial()
        assert _lift(matrix, p, genus, n) == expected


def test_python_api_models_and_environment(binary, monkeypatch):
    monkeypatch.setenv("HYPELLFROB_EULER", binary)
    ps = [101, 103, 107, 109]
    for q in (x**5+x**3+1, (x-1)*(x**5+x+3), x**6+x+1):
        good = [p for p in ps if q.discriminant() % p and
                (odd_model(q) is None or odd_model(q).discriminant() % p)]
        expected = frobenius_polynomials(q, good, backend="legacy")
        assert frobenius_polynomials(q, good, backend="remainder") == expected
        monkeypatch.setenv("HYPELLFROB_BACKEND", "remainder")
        assert frobenius_polynomials(q, good) == expected
    with pytest.raises(ValueError, match="backend"):
        frobenius_polynomials(x**5+x+1, ps, backend="unknown")


def test_lfunction_coefficients_use_new_backend(binary, monkeypatch):
    from hyperell_regulator import HyperellipticCurveWithReduction
    from hyperell_regulator.lfunction import CurveLFunction, StandardLFunction
    monkeypatch.setenv("HYPELLFROB_EULER", binary)
    # Conductor 249, with only odd bad primes; the standard L-function has
    # no automatic bad factor at 2 for the other test curve.
    q = (x**3+1)**2 + 4*(x**2+x)
    for cls in (CurveLFunction, StandardLFunction):
        monkeypatch.setenv("HYPELLFROB_BACKEND", "legacy")
        expected = cls(HyperellipticCurveWithReduction(q)).dirichlet_coefficients(250)
        monkeypatch.setenv("HYPELLFROB_BACKEND", "remainder")
        actual = cls(HyperellipticCurveWithReduction(q)).dirichlet_coefficients(250)
        assert actual == expected
