# rforest test on conductor 388, 1 November 2026

This is the historical check of commit `434e654`, before even-degree
support was added. The [even-degree extension and new results](frobenius_rforest_even_20261101.md)
supersede its backend limitation; the measurements below are preserved.

The new backend does **not accelerate `388.a.776.1`**. It supports monic odd
models, while this curve has no rational Weierstrass point. Its completed
square model is

```
f = -x^4 + 2*x^2 + x
h = x^3 + x + 1
F = 4*f + h^2 = x^6 - 2*x^4 + 2*x^3 + 9*x^2 + 6*x + 1
```

Exact Sage computations give `F.roots(QQ) == []` and `odd_model(F) is None`.
For this degree-six model the points at infinity are not branch points, so
there is no rational branch point to move to infinity to obtain an odd
model over Q. The current rforest implementation therefore declines every
job. Hybrid `rforest` uses the existing native even-degree algorithm, with
Sage fallback for primes the native algorithm declines. Supporting this
example in the simultaneous algorithm requires an even-degree extension.

## Checks performed

- Below 1,000,000 there are 78,496 odd good-prime jobs (the bad primes are
  2 and 97). Strict `rforest-only` supplied **zero matrices**, declining
  every job as expected. This fast coverage check does not compute the
  Frobenius matrices and is not a speed measurement.
- At primes 23, 31, 37, 41 and 101, precision two, the hybrid and legacy
  4-by-4 matrices agree entry for entry. All five lifted characteristic
  polynomials also agree with independent PARI `lfungenus2`/`lfuneuler`.
- Both `CurveLFunction` and `StandardLFunction` produce exactly the same
  first 128 Dirichlet coefficients with `HYPELLFROB_BACKEND=legacy` and
  `HYPELLFROB_BACKEND=rforest`. This checks fallback, including the bad
  factors at 2 and 97, rather than rforest arithmetic.
- The three conductor-388 regression tests passed in 56.28 seconds. These
  cover the new fallback check and the two existing tests of the factor
  at 2, including a translated model.

A preliminary six-job legacy batch through prime 1009 exceeded its
120-second timeout. Its incomplete result is recorded in the raw report;
it supplies no timing comparison. The smaller completed sample runs were
concurrent with other checks, so their timings are not used to estimate a
speed ratio. No full million-bound even-degree fallback computation was
attempted.

## Independent L-value check

The full reference calculation uses PARI good H1 factors, the package's
automatically detected bad factors, and the package's standard factors,
Dirichlet coefficients, values and functional-equation checks. This is the
same good-factor protocol as the saved `conductor_388_lvalue.json` report;
it does not exercise rforest. All 655 good-prime reference factors were
computed afresh.

| Requested digits | Dirichlet coefficients | L''(std,1) | Functional-equation error |
| ---: | ---: | ---: | ---: |
| 8 | 4,928 | 2.6170256284 | 2.65e-23 |
| 15 | 4,928 | 2.61702562836175902 | 6.62e-24 |

Both values reproduce the [saved conductor-388 result](conductor_388_lvalue.json).
PARI warns that its extra guard precision requires more coefficients
(requests range from 9,764 to 17,998). As in that saved report, the
4,928-coefficient value agreement and small functional-equation errors are
numerical checks, not a rigorous accuracy certificate.

The standard conductor is 37,636. Its bad factors are

```
P_2(T)  = 1 + T + 2*T^2 - 4*T^3
P_97(T) = 1 - 83*T - 1261*T^2 - 9409*T^3
```

## Reproduction

From the `experimental/rforest` worktree, with the native driver built:

```sh
HYPELLFROB_EULER=$PWD/hypellfrob-threaded/build/euler \
  sage -python -m pytest tests/test_average_frobenius.py \
    tests/test_lfunction.py -k conductor_388 -q
```

The tested executable is from commit `434e654`. The
[raw report](frobenius_rforest_conductor388_20261101.json) records its hash,
all sample matrices, exact coverage counts, the timeout, and reference
L-value results. The main branch is unchanged.

To reproduce the independent full L-value check (this uses PARI good
factors, with no native rforest computation):

```python
from sage.all import QQ, PolynomialRing, pari, prime_range
from hyperell_regulator import HyperellipticCurveWithReduction
from hyperell_regulator.lfunction import StandardLFunction

x = PolynomialRing(QQ, "x").gen()
C = HyperellipticCurveWithReduction(-x**4 + 2*x**2 + x, x**3 + x + 1)
L = StandardLFunction(C)
n = 4928
bad = set(C.bad_primes())
for p in prime_range(3, n + 1):
    if p not in bad:
        e = L._pari_L.lfuneuler(p)
        L._h1cache[p] = L._RT(list(e.denominator().Vecrev()))
print(L.value(1, 2, n, 15))
print(L.check_functional_equation(n, 15))
```
