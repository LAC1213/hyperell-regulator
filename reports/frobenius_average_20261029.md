# Harvey 2013 Frobenius implementation and comparison, 29 October 2026

The new C++ backend implements the simultaneous reduction and accumulating
remainder tree algorithm in [arXiv:1210.8239v3](https://arxiv.org/pdf/1210.8239).
All tested full Frobenius matrices agree exactly with the existing Harvey
2007 implementation modulo the requested power of the prime. Independent
Sage Euler factors and both package L-function coefficient paths agree too.
**A significant practical speedup has not been achieved.** The new backend
remains experimental and opt-in; the existing backend remains the default.

These are historical timings of the first implementation. A subsequent
complexity audit found that it recomputed `root_modulus / 2` at every
product-tree node, introducing quadratic bit-work despite correct arithmetic.
It also retained right-child products after their final use. Those findings
mean these timings cannot establish the practical crossover of Harvey's
algorithm. The corrected large-bound experiment is recorded in the
[1,000,000-bound report](frobenius_average_1000000_20261029.md).

The implementation, arithmetic/integration tests, and `benchmark.py` are
archived on `experimental/harvey-remainder-tree`. The historical
[raw timings](frobenius_average_20261029.json) are retained here.

The test suite covers genus 1, 2, and 3; p-adic precisions 1, 2, 3, and 4;
negative coefficients; zero constant terms; randomized integral curves;
unsorted and duplicate jobs; mixed precisions; parallel trees; bad primes;
and unsupported models. The native comparisons use `remainder-only`, which
cannot hide a failed computation through fallback. The L-function test
compares Dirichlet coefficients through index 250 for both `CurveLFunction`
and `StandardLFunction`, on the conductor-249 curve. The 65-test arithmetic,
L-function, and build suite passed. A subsequently added strict-mode build
check regression also passed, alongside the 13 existing build tests.

Timings below are medians of three runs on this machine (Linux x86-64,
16 logical CPUs), for `Q = x^5 + x^3 + 1`. Both algorithms used the same
binary, g++ 15.2.0 with `-O2 -std=c++17`, a checksum-verified private NTL
11.6.0 built with threads enabled, and system GMP. Run order alternated
between algorithms. These are local wall-clock measurements on a shared
machine, not dedicated-hardware performance estimates.

| Prime bound | Jobs | Workers | Existing (s) | Remainder (s) | Existing / remainder |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 1,000 | 165 | 1 | 0.394 | 0.332 | 1.19 |
| 1,000 | 165 | 4 | 0.081 | 0.087 | 0.93 |
| 5,000 | 666 | 1 | 1.561 | 1.500 | 1.04 |
| 5,000 | 666 | 4 | 0.588 | 0.614 | 0.96 |
| 20,000 | 2,259 | 1 | 9.946 | 13.107 | 0.76 |
| 20,000 | 2,259 | 4 | 3.909 | 5.749 | 0.68 |

Each job requests the package's reconstruction precision, using its exact
integer equivalent in the benchmark. Jobs include every odd good prime
below the bound. Timings include process startup, output, and the existing
native algorithm for primes that the remainder backend declines. Both
backends decline five small primes that require Sage fallback. At bound
20,000, the strict remainder path computes 2,249 matrices, and native
fallback supplies another five. All returned entries agreed in every run,
including the separate strict-mode comparison.

These measurements exclude Sage characteristic-polynomial lifting, Sage
fallback for the five small primes, and L-function evaluation. They measure
the native Frobenius stage rather than end-to-end L-value computation.
They provide no evidence for an end-to-end L-function speedup. The paper's
asymptotic theorem does not guarantee a crossover at these bounds: this
implementation constructs several matrix-product trees at elevated working
precision, whereas the existing single-prime implementation is already
optimized. Further work would need to reduce those costs and measure larger
bounds before changing the default.

Reproduce with a freshly built accelerator:

```sh
sage -python hypellfrob-threaded/benchmark.py \
    --binary hypellfrob-threaded/build/euler \
    --bounds 1000 5000 20000 --workers 1 4 --repeats 3 \
    --output /tmp/frobenius-benchmark.json

HYPELLFROB_EULER="$PWD/hypellfrob-threaded/build/euler" \
    sage -python -m pytest tests/test_average_frobenius.py \
    tests/test_lfunction.py tests/test_native_build.py -q

sage -python -m hyperell_regulator.frobenius --check
```

Select the new implementation with
`frobenius_polynomials(f, primes, nthreads=4, backend="remainder")`, or set
`HYPELLFROB_BACKEND=remainder` for the existing L-function API. Models outside
the paper's scope retain the existing native/Sage fallback behavior.
