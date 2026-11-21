r"""Frobenius polynomials by Harvey's simultaneous algorithm.

The rforest implementation handles odd and even degree models.
Unsupported primes are computed by Sage.
"""

import os
import subprocess

__all__ = ["euler_binary", "odd_model", "frobenius_polynomials",
           "frobenius_polynomial", "FrobeniusBackendError"]


class FrobeniusBackendError(ArithmeticError):
    """A bulk native computation failed; per-prime fallback would be too slow."""


def euler_binary():
    r"""The Frobenius executable, or ``None`` if unavailable.

    EXAMPLES::

        sage: from hyperell_regulator.frobenius import euler_binary
        sage: euler_binary() is None or euler_binary().endswith("euler")
        True
    """
    env = os.environ.get("HYPELLFROB_EULER")
    if env:
        return env if os.path.isfile(env) and os.access(env, os.X_OK) else None
    here = os.path.dirname(os.path.abspath(__file__))
    # installed inside the package by setup.py; then the source tree, so a
    # checkout that has run ``make`` in hyperell-frobenius works uninstalled
    for rel in (("bin", "euler"),
                ("..", "hyperell-frobenius", "build", "euler")):
        guess = os.path.normpath(os.path.join(here, *rel))
        if os.path.isfile(guess) and os.access(guess, os.X_OK):
            return guess
    return None


def odd_model(f):
    r"""A monic integral odd degree model, or ``None`` if none exists over `\QQ`.

    EXAMPLES::

        sage: from hyperell_regulator.frobenius import odd_model
        sage: from sage.all import ZZ, PolynomialRing
        sage: R = PolynomialRing(ZZ, 'x'); x = R.gen()
        sage: odd_model(x**5 + x**3 + 1)             # already odd and monic
        x^5 + x^3 + 1
        sage: f = (x - 1) * (x**5 + x + 3)           # even, with a rational root
        sage: g = odd_model(f); g.degree(), g.is_monic()
        (5, True)

        sage: odd_model(x**6 + x + 1) is None
        True
    """
    from sage.all import QQ, PolynomialRing, ZZ, lcm

    R = PolynomialRing(QQ, "x")
    f = R(f)
    d = f.degree()
    if d % 2 == 0:
        roots = f.roots(QQ, multiplicities=False)
        if not roots:
            return None
        r = roots[0]
        u = R.gen()
        f = sum(f[i] * u ** (d - i) * (r * u + 1) ** i for i in range(d + 1))
        d = f.degree()
        if d % 2 == 0:                       # r was a multiple root
            return None
    c = f.leading_coefficient()
    q = [f[i] * c ** (d - 1 - i) for i in range(d + 1)]
    m = lcm([QQ(v).denominator() for v in q] + [1])
    q = [q[i] * m ** (2 * (d - i)) for i in range(d + 1)]
    return PolynomialRing(ZZ, "x")(q)


def _bound(p, g):
    r"""Kedlaya's bound `p^N \ge 2\binom{2g}{g}p^{g/2}`."""
    from sage.all import RR, ZZ, binomial

    M = 2 * binomial(2 * g, g) * RR(p).sqrt() ** g
    B = ZZ(M.ceil()).exact_log(p)
    if p ** B < M:
        B += 1
    return int(B)


def _lift(M, p, g, N):
    r"""Lift the Frobenius characteristic polynomial from a matrix modulo `p^N`."""
    from sage.all import ZZ, matrix

    if g == 2:
        # Newton identities avoid constructing a symbolic characteristic
        # polynomial for every one of tens of thousands of 4x4 matrices.
        trace = sum(ZZ(M[i][i]) for i in range(4))
        trace2 = sum(ZZ(M[i][j])*ZZ(M[j][i]) for i in range(4) for j in range(4))
        a = [(trace**2-trace2)//2, -trace, ZZ(1)]
    else:
        a = list(matrix(ZZ, M).charpoly().list()[g:2 * g + 1])
    ppow = p ** N
    a = [x % ppow for x in a]
    a = [x if 2 * x < ppow else x - ppow for x in a]
    a = [a[g - i] * p ** (g - i) for i in range(g)] + a
    return ZZ["x"](a)


def _run(fcoeffs, jobs, nthreads, backend=None, progress=None):
    """Frobenius matrices for the jobs ``(p, N)``, keyed by prime."""
    from ._subprocess import run as _native_run
    backend = backend or os.environ.get("HYPELLFROB_BACKEND", "rforest")
    if backend != "rforest":
        raise ValueError("backend must be 'rforest'; the per-prime backends have been removed")
    binary = euler_binary()
    if binary is None or not jobs:
        if jobs and len(jobs) > 1000:
            raise FrobeniusBackendError('native Frobenius executable unavailable for %d primes; '
                'refusing bulk Sage fallback; build/install the Frobenius accelerator' % len(jobs))
        return {}
    inp = " ".join(str(int(c)) for c in fcoeffs) + "\n"
    inp += "".join("%d %d\n" % (p, N) for p, N in jobs)
    if nthreads is None:
        configured = os.environ.get("HYPELLFROB_THREADS")
        workers = max(1, int(configured)) if configured else (os.cpu_count() or 1)
    else:
        workers = max(1, int(nthreads))
    args = [binary, str(workers), backend]
    try:
        out = _native_run(args, input=inp, check=True,
                          timeout=float(os.environ.get('HYPELLFROB_TIMEOUT', '900')),
                          progress=progress, label='Euler factors').stdout
    except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired) as exc:
        if len(jobs) > 1000:
            details = str(exc)
            if isinstance(exc, subprocess.CalledProcessError):
                details += '; '+(exc.stderr or '').strip()[-1000:]
            raise FrobeniusBackendError('native Frobenius computation failed for %d primes: %s; '
                'refusing bulk Sage fallback' % (len(jobs), details)) from exc
        return {}
    got = {}
    for line in out.splitlines():
        parts = line.split()
        if len(parts) < 2 or parts[1] != "1":
            continue
        p, dim = int(parts[0]), int(parts[2])
        vals = [int(v) for v in parts[3:3 + dim * dim]]
        got[p] = [vals[r * dim:(r + 1) * dim] for r in range(dim)]
    return got


def frobenius_polynomials(f, primes, nthreads=None, *, backend=None, progress=None):
    r"""Frobenius characteristic polynomials of `y^2=f(x)`, keyed by prime.

    ``f`` has integral coefficients and ``primes`` must have good reduction.
    The native backend is ``rforest``; unsupported jobs use Sage.
    ``nthreads`` overrides ``HYPELLFROB_THREADS``, then the logical CPU count.

    EXAMPLES::

        sage: from hyperell_regulator.frobenius import frobenius_polynomials
        sage: from sage.all import ZZ, PolynomialRing
        sage: R = PolynomialRing(ZZ, 'x'); x = R.gen()
        sage: P = frobenius_polynomials(x**5 + x**3 + 1, [11, 13, 17])
        sage: [P[p].degree() for p in (11, 13, 17)]
        [4, 4, 4]
        sage: all(P[p][0] == p**2 for p in (11, 13, 17))     # a_0 = q^g
        True
    """
    from sage.all import GF, HyperellipticCurve, PolynomialRing, ZZ

    R = PolynomialRing(ZZ, "x")
    fz = R(f)
    g = (fz.degree() - 1) // 2
    primes = [int(p) for p in primes]
    jobs = [(p, _bound(p, g)) for p in primes]
    # A rational Weierstrass point gives a smaller odd degree model,
    # isomorphic over QQ. Otherwise use simultaneous even degree reductions.
    selected_backend = backend or os.environ.get("HYPELLFROB_BACKEND", "rforest")
    model = odd_model(fz)
    twist = ZZ(1)
    if model is None and not fz.is_monic():
        # c^(d-1)*f(X/c) is monic and integral. In even degree it is
        # isomorphic over Q to the c-twist of the input curve. The common
        # rational transformation preserves sharing between primes; undo
        # the twist on characteristic polynomials via (c/p).
        d = fz.degree()
        twist = fz.leading_coefficient()
        model = R([fz[i]*twist**(d-1-i) for i in range(d)]+[1])
    options = {'progress': progress} if progress is not None else {}
    got = _run((model if model is not None else fz).list(), jobs, nthreads, selected_backend, **options)

    missing = [p for p in primes if p not in got or not twist%p]
    if len(primes) > 1000 and len(missing) > 64:
        raise FrobeniusBackendError('native Frobenius computation declined %d of %d primes; '
            'refusing bulk Sage fallback' % (len(missing), len(primes)))

    out = {}
    for p, N in jobs:
        if p in got and twist%p:
            cp = _lift(got[p], p, g, N)
            out[p] = cp(-R.gen()) if twist != 1 and twist.kronecker(p) == -1 else cp
        else:
            H = HyperellipticCurve(PolynomialRing(GF(p), "x")(fz))
            out[p] = R(H.frobenius_polynomial())
    return out


def frobenius_polynomial(f, p, *, backend=None):
    r"""The Frobenius characteristic polynomial at one good prime.

    EXAMPLES::

        sage: from hyperell_regulator.frobenius import frobenius_polynomial
        sage: from sage.all import ZZ, PolynomialRing
        sage: R = PolynomialRing(ZZ, 'x'); x = R.gen()
        sage: frobenius_polynomial(x**5 + x**3 + 1, 37).degree()
        4
    """
    return frobenius_polynomials(f, [p], backend=backend)[int(p)]


# Fixed full-matrix references for odd and even degree, verified against
# the former independent implementations and Sage characteristic polynomials.
# Kept here so build isolation can check arithmetic without importing Sage.
_CHECK_CASES = (
    (
        (
            '1 0 0 1 0 1\n'
            '101 3\n'
            '103 3\n'
            '107 3\n'
            '109 3\n'
            '31 3\n'
        ),
        (
            '101 1 4 1011616 996870 945785 288208 460661 90799 98296 265518 397738 476518 632995 265953 658419 476922 143336 325199\n'
            '103 1 4 689070 435175 295697 185257 399846 888272 442156 365773 1000336 732536 812626 320027 407777 109592 101499 888207\n'
            '107 1 4 1037579 558968 609797 378854 503863 1165979 895436 1133730 193028 141775 1186127 1210054 474010 451540 119495 285437\n'
            '109 1 4 1170878 640702 152452 271988 507286 748612 1080818 597880 992663 888786 543370 1020913 877559 412129 531249 127182\n'
            '31 0\n'
        ),
    ),
    (
        (
            '1 6 9 2 -2 0 1\n'
            '41 2\n'
            '43 2\n'
            '31 2\n'
        ),
        (
            '41 1 4 205 205 725 396 164 328 1490 402 779 861 5 40 1025 82 1659 1144\n'
            '43 1 4 645 1763 1235 43 1247 129 1247 565 1075 0 1204 1806 0 1548 645 1720\n'
            '31 0\n'
        ),
    ),
)


def _check_binary(binary):
    """Check native matrices in both degree parities and with several worker counts."""
    workers = tuple(dict.fromkeys((1, 4, os.cpu_count() or 1)))
    for batch, expected in _CHECK_CASES:
        for count in workers:
            result = subprocess.run([binary, str(count), "rforest"], input=batch, text=True,
                                    capture_output=True, check=True, timeout=30)
            if result.stdout != expected:
                raise RuntimeError("Frobenius matrices differ from verified reference: " +
                                   result.stdout)
    return workers


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Check the installed Frobenius accelerator.")
    parser.add_argument("--check", action="store_true", required=True)
    parser.add_argument("--binary", help="check this executable instead of the installed one")
    args = parser.parse_args()
    binary = args.binary or euler_binary()
    if binary is None:
        parser.exit(1, "Frobenius accelerator not found; computations use Sage.\n")
    try:
        workers = _check_binary(binary)
    except (OSError, subprocess.SubprocessError, RuntimeError) as exc:
        parser.exit(1, "Accelerator failed: %s\n%s\n" %
                    (exc, getattr(exc, "stderr", "")))
    print("Rforest Frobenius checks passed with worker counts %s: %s" %
          (", ".join(map(str, workers)), binary))
