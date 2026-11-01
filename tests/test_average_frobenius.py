"""Harvey 2013 against full 2007 matrices and independent Sage Euler factors.

The strict native mode is essential: fallback must not hide an arithmetic bug
or an implementation that simply declines all the new backend's work.
"""
import random
import os
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


def test_native_default_uses_rforest_and_retains_fallback(binary):
    batch = "1 6 9 2 -2 0 1\n41 2\n43 2\n31 2\n"
    env = dict(os.environ, HYPELLFROB_PROFILE="1")
    default = subprocess.run([binary, "4"], input=batch, text=True,
                             capture_output=True, check=True, timeout=30, env=env)
    legacy = subprocess.run([binary, "4", "legacy"], input=batch, text=True,
                            capture_output=True, check=True, timeout=30, env=env)
    assert "rforest a=" in default.stderr
    assert "rforest a=" not in legacy.stderr
    assert default.stdout == legacy.stdout
    assert all(line.split()[1] == "1" for line in default.stdout.splitlines())


def test_python_default_and_legacy_overrides(binary, monkeypatch):
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
    actual = frobenius_polynomials(q, ps, nthreads=4)
    assert calls[-1][0] == [binary, "4", "rforest"]
    assert calls[-1][1].split()[-1] == "1"  # default monicizes even input
    expected = frobenius_polynomials(q, ps, nthreads=4, backend="legacy")
    assert calls[-1][0] == [binary, "4", "legacy"]
    assert actual == expected
    monkeypatch.setenv("HYPELLFROB_BACKEND", "legacy")
    assert frobenius_polynomials(q, ps, nthreads=4) == expected
    assert calls[-1][0] == [binary, "4", "legacy"]
    assert frobenius_polynomials(q, ps, nthreads=4, backend="rforest") == expected
    assert calls[-1][0] == [binary, "4", "rforest"]


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
    expected = native(binary, q, jobs, "legacy")
    actual = native(binary, q, jobs, backend)
    assert all(row[1] == 1 for row in actual), actual
    assert actual == expected


@pytest.mark.parametrize("backend", ["rforest-only"])
def test_unsorted_duplicate_jobs_mixed_precision_and_threads(binary, backend):
    q = x**5 + x**3 + 1
    jobs = [(211, 2), (101, 3), (211, 2), (151, 1), (211, 3), (103, 3)]
    expected = native(binary, q, jobs, "legacy")
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
            assert actual == native(binary, q, jobs, "legacy")


@pytest.mark.parametrize("q,jobs", [
    (x**5 + x**3 + 1, [(3, 2), (31, 2), (53, 2), (61, 2), (101, 2)]),
    (x**5 + x + 101, [(101, 2), (103, 2)]),  # c0 divisible by a good prime
    (x**6 + x + 1, [(31, 2), (101, 2)]),
    (2*x**5 + x + 1, [(101, 2)]),
    ((x+1)**2*(x**3+2), [(101, 2)]),
])
@pytest.mark.parametrize("backend", ["rforest-only"])
def test_declines_and_native_fallback(binary, q, jobs, backend):
    strict = native(binary, q, jobs, backend)
    # Each case deliberately includes at least one unsupported job.
    assert any(row[1] == 0 for row in strict)
    if q.is_monic() and q.is_squarefree():
        assert native(binary, q, jobs, backend.removesuffix("-only")) == native(binary, q, jobs, "legacy")


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
        expected = frobenius_polynomials(q, good, backend="legacy")
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
        monkeypatch.setenv("HYPELLFROB_BACKEND", "legacy")
        expected = cls(HyperellipticCurveWithReduction(q)).dirichlet_coefficients(250)
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
    assert actual == native(binary, q, jobs, "legacy", workers=4)
    reference = pari.lfungenus2([pari(str(f)), pari(str(h))])
    for (p, n), row in zip(jobs, actual):
        assert row[:3] == [p, 1, 4]
        matrix = [row[3+4*i:7+4*i] for i in range(4)]
        factor = R(reference.lfuneuler(p).denominator())
        assert _lift(matrix, p, 2, n) == factor.reverse()

    monkeypatch.setenv("HYPELLFROB_EULER", binary)
    for cls in (CurveLFunction, StandardLFunction):
        monkeypatch.setenv("HYPELLFROB_BACKEND", "legacy")
        expected = cls(HyperellipticCurveWithReduction(f, h)).dirichlet_coefficients(128)
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
def test_even_degree_full_matrices(binary, q, precision):
    # One good prime beyond the common p-admissibility bound. Full entries,
    # including higher p-adic digits, agree in the reference's residue basis.
    bound = q.degree()*(4*precision-1)
    p = next(p for p in prime_range(bound+1, bound+100) if q.discriminant()%p)
    jobs = [(int(p), precision)]
    actual = native(binary, q, jobs, "rforest-only", workers=4)
    assert actual[0][:3] == [p, 1, q.degree()-2]
    assert actual == native(binary, q, jobs, "legacy", workers=4)


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
    expected = native(binary, q, jobs, "legacy", workers=4)
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
    expected = native(binary, q, jobs, "legacy", workers=4)
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
    expected = native(binary, q, jobs, "legacy", workers=4)
    batch = " ".join(map(str, q.list())) + "\n" + "".join("%d %d\n" % job for job in jobs)
    result = subprocess.run([binary, "4", "rforest-only", str(kappa)], input=batch,
                            capture_output=True, text=True, check=True, timeout=90)
    assert [list(map(int, line.split())) for line in result.stdout.splitlines()] == expected
