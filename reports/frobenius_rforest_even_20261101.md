# Even-degree rforest support, 1 November 2026

The `experimental/rforest` branch now accelerates even-degree models,
including **388.a.776.1**, using the actual rforest C library. For conductor
388, a fresh 15-digit standard L-value calculation with 20,000 coefficients
took **2.89 seconds**. A 30-digit calculation with 100,000 coefficients took
**12.40 seconds**; all 100,000 coefficients and the resulting value matched
an independent-factor reference exactly.

Below 1,000,000, **all 78,485 Euler factors supplied by strict rforest mode
matched independent smalljac results**, with zero mismatches. The remaining
11 small-prime factors were checked separately through the Python fallback
and also matched. This accounts for all **78,496 odd good primes** for the
curve. The direct native run took **187.16 seconds** with
**3.27 GiB** peak resident memory.

The [earlier restriction to odd models](frobenius_rforest_conductor388_20261101.md)
is superseded by this implementation. Main remains at `30f22b4`; the work
is isolated in `build/rforest-worktree` on `experimental/rforest`.

## Method

The backend uses direct integer reductions over the original monic even
model. It computes the affine `(2g+1)`-dimensional Frobenius action, verifies
its residue identity, and restricts to the residue-free projective basis
used by the existing reference implementation. Horizontal/diagonal path
ordering and a bounded initial prefix avoid the even-degree singularity.
For this curve, common-factor cancellation reduces the largest transition
degree from 22 to 7, and the longest initial prefix has two blocks.

No branch point or extension field is needed. The proposed quadratic-root
method is not universally available: this sextic is irreducible modulo 11,
so it has no root in the unramified quadratic extension at that prime.
The direct method also retains one integer polynomial across all primes,
which is what allows the remainder forest to share products. The
[algorithm notes](../hypellfrob-threaded/EVEN-DEGREE.md) give the path
inequalities, prefix construction, precision conditions and residue formula.

The Python caller additionally handles nonmonic integral even models by
monicizing a common quadratic twist and removing its local character from
the characteristic polynomial. Six cases, including negative leading
coefficients and a good prime dividing the original leading coefficient,
were checked against Sage. The native matrix interface still requires a
monic model. `legacy` and the original `remainder` mode retain their
existing algorithms; selecting `rforest` enables the new paths.

## Correctness and native performance

The completed-square model used throughout is

```
f = -x^4 + 2*x^2 + x,   h = x^3 + x + 1
F = 4*f + h^2 = x^6 - 2*x^4 + 2*x^3 + 9*x^2 + 6*x + 1.
```

The good odd primes below one million exclude 97; 2 is also bad. Precision
is the package's reconstruction bound (4 at 3, 3 at 5/7/11, and 2 thereafter).
Strict mode supplies every good prime at least 41. Its eleven declines are
3, 5, 7, 11, 13, 17, 19, 23, 29, 31 and 37. The first six use Sage fallback
in Python; hybrid native mode supplies the other five with the reference
kernel. Separate small-prime comparisons cover all eleven.

At the million bound the independent oracle compares reconstructed characteristic
polynomials, **not matrix entries**. Full matrix-entry comparisons with
the existing even-degree reference cover genera 1--3 and precisions 1--4.
An additional 100,000-bound check verifies every coefficient of all 9,579
native characteristic polynomials modulo the requested p-power before
applying the usual Weil-bound lifting.
They include dense and negative coefficients, zero constant terms, mixed
precisions, duplicate and unsorted jobs, shifted starts, both polynomial
evaluators, subtree boundaries and worker counts 1/4/16.

A fresh alternating comparison at bound 200 includes every native fallback
job, so the two processes return identical output for all 44 jobs,
including all 38 matrices supplied by the reference:

| Backend | Workers | Median wall time | Median CPU time | Repetitions |
| --- | ---: | ---: | ---: | ---: |
| Legacy | 8 | 11.543 s | 78.013 s | 3 |
| rforest, including native fallback | 8 | 0.322 s | 0.855 s | 3 |

This is about **35.9x faster** in wall time at this small bound.
These native timings include startup and output and exclude Sage fallback,
polynomial lifting and L-function evaluation. Larger-bound timings below
are single runs, and the shared machine shows substantial timing variation;
no speed ratio is extrapolated to a million-bound reference run.

Exploratory 100,000-bound runs all matched the independent oracle:

| Workers | Kappa | Wall seconds | Peak RSS MiB |
| ---: | ---: | ---: | ---: |
| 8 | 6 | 6.53 | 455.7 |
| 8 | 7 | 8.64 | 410.9 |
| 16 | 6 | 9.15 | 755.7 |
| 16 | 7 | 13.26 | 691.5 |

Eight workers and kappa 6 were fastest among these tested settings. The
automatic forest split already chooses kappa 6 at 100,000 and 1,000,000.
An 8-GiB address-space cap applied to all measured native processes.
The larger 1,000,000 run used eight workers and the automatic setting.

Smalljac 4.1.3 with ff_poly 1.2.7 was built separately under `/tmp`, solely
as an independent oracle. Its own coarse timer reported about 20 seconds
for all 78,496 factors using 16 processes. It was faster here for Euler
factors alone; rforest also computes full p-adic matrices. No globally
optimal Euler-only performance is claimed. Neither smalljac nor ff_poly is
required by, or bundled with, the new production backend.

Oracle sources: [official smalljac distribution](https://math.mit.edu/~drew/smalljac.html)
and Kedlaya--Sutherland, [Computing L-series of hyperelliptic curves](https://arxiv.org/abs/0801.2778).
Archive and data hashes, all repetitions, raw metrics, source/binary hashes,
and the exact comparisons are in the [raw report](frobenius_rforest_even_20261101.json).
The native build uses GCC/G++ 15.2.0, private thread-safe static NTL 11.6.0,
system GMP, `-O2` for C++17 and `-O3` for the vendored C11 rforest code.

## End-to-end standard L-value

The standard conductor is 37,636. Bad factors are computed by the existing
package code from the original `(f,h)` model. With eight workers:

| Coefficients | Digits | Fresh-object wall time | L''(std,1) | FE error |
| ---: | ---: | ---: | --- | ---: |
| 4,928 | 15 | 0.97 s, median of 3 | 2.61702562836175902 | 6.62e-24 |
| 20,000 | 15 | 2.89 s | 2.61702562836175902 | 6.31e-30 |
| 100,000 | 30 | 12.40 s | 2.61702562836175901615851458820520 | 6.84e-49 |

Every row starts with a new L-function object and no H1 cache. Timings
include setup, good-factor computation, lifting, coefficient expansion,
L-function construction, value evaluation and the FE check; Sage import
and subsequent coefficient-hash verification are excluded. Recorded native
calls confirm the use of `euler 8 rforest`. The 30-digit reference used
independently generated smalljac good factors, with no native calls. Its
100,000 exact Dirichlet coefficients have the same hash as the rforest
run, and its value is identical. That reference's timing excludes oracle
generation, so it is used for correctness rather than timing comparison.

As in the historical report, 4,928 coefficients trigger PARI's additional
guard-precision requests. The 20,000- and 100,000-coefficient runs produce
no coefficient-count warnings. Numerical values and FE checks are not
interval-arithmetic accuracy certificates. The prior PARI good-factor pass
at 4,928 coefficients took 507.4 seconds, but that is a historical single
run and is not used as a fresh paired speedup estimate.

## Validation and use

The 159-case targeted suite passed, including 35 new even-degree cases;
six additional nonmonic cases passed. After the API/configuration changes,
22 API/build/model checks and all 19 Frobenius doctests passed. Eight
model-dispatch/API cases were also rerun after preserving legacy dispatch
for nonmonic input. The binary
self-check requires genuine even-degree support and agrees at 1/4/16
workers. AddressSanitizer with leak detection passed four additional
cases and compared 2,256 full matrices, including a dense batch and
precision-four genus-three matrices, with no reported errors or leaks.
Source-distribution checks include the new sources, headers, algorithm
notes, benchmarks and licenses and exclude native build products.

From the experimental worktree, enable the backend and choose workers:

```sh
export HYPELLFROB_BACKEND=rforest
export HYPELLFROB_THREADS=8
```

An explicit `nthreads` argument overrides the environment setting. The
backend remains opt-in. For the full 15-digit example:

```python
from sage.all import QQ, PolynomialRing
from hyperell_regulator import HyperellipticCurveWithReduction
from hyperell_regulator.lfunction import StandardLFunction
x = PolynomialRing(QQ, "x").gen()
C = HyperellipticCurveWithReduction(-x**4+2*x**2+x, x**3+x+1)
L = StandardLFunction(C)
print(L.value(1, 2, num_coeffs=20000, digits=15))
print(L.check_functional_equation(num_coeffs=20000, digits=15))
```

After building the official test-only smalljac/ffpoly tools, the independent
million-bound comparison can be reproduced with:

```sh
lpdata /tmp/conductor388 'x^6-2*x^4+2*x^3+9*x^2+6*x+1' 1000000 1
sage -python hypellfrob-threaded/benchmark_rforest.py \
    --binary hypellfrob-threaded/build/euler --bound 1000000 \
    --workers 8 --coefficients 1 6 9 2 -2 0 1 \
    --euler-reference /tmp/conductor388_lpdata.txt \
    --output /tmp/conductor388-million.json
```

For a fresh small-bound comparison including all native fallback jobs:

```sh
sage -python hypellfrob-threaded/benchmark_rforest.py \
    --binary hypellfrob-threaded/build/euler --bound 200 \
    --workers 8 --repeats 3 --hybrid --coefficients 1 6 9 2 -2 0 1 \
    --output /tmp/conductor388-paired.json
```
