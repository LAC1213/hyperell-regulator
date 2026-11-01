# rforest Frobenius backend, 1 November 2026

The `experimental/rforest` branch uses the actual upstream rforest C
library for Harvey's simultaneous Frobenius algorithm. On
`Q = x^5 + x^3 + 1`, fresh paired measurements at bound 100,000 show a
**2.06x wall-time speedup** over the existing backend.
At bound 1,000,000, all **78,485** eligible full Frobenius matrices
matched the saved legacy results exactly: **1,255,760 entries, zero mismatches**.
The million-bound rforest run completed in 305.2 seconds.

This code is isolated from `main` in `build/rforest-worktree`, on branch
`experimental/rforest`. The previous C++ tree is retained as `remainder`;
`legacy` remains the default. Select the new backend with
`HYPELLFROB_BACKEND=rforest` for both L-function classes, or
`frobenius_polynomials(f, primes, backend="rforest")`.

## Measurements

| Bound | Model | Backend | Workers | Wall seconds | CPU seconds | Peak RSS MiB |
| ---: | --- | --- | ---: | ---: | ---: | ---: |
| 100,000 | `x^5+x^3+1` | legacy | 16 | 38.99 | 301.93 | 150.0 |
| 100,000 | same | rforest-only | 16 | 18.89 | 83.31 | 496.5 |
| 1,000,000 | same | rforest-only | 16 | 305.20 | 1290.25 | 3277.1 |
| 20,000 | conductor-249 odd model | legacy | 16 | 4.22 | 27.44 | 122.3 |
| 20,000 | same | rforest-only | 16 | 3.09 | 10.68 | 277.5 |

The 100,000 rows are medians of three runs per backend, alternating run
order, with peak RSS taken over all three. Rforest used
72.4% less CPU time. The other bounds are single runs.
The conductor-249 model is
`x^5 + 13*x^4 + 72*x^3 + 240*x^2 + 384*x + 256`, the integral monic odd
model of `(x^3+1)^2 + 4*(x^2+x)` used by the L-function tests. Its 2,250
strict-mode matrices all matched legacy.

The 100,000 comparison has 9,589 odd good-prime jobs and 9,579 strict
rforest matrices, all identical to legacy. At 1,000,000 there are 78,495
jobs, excluding bad primes 53 and 61. Precision is the package's Kedlaya
reconstruction bound: one job at precision 4, three at precision 3, and
78,491 at precision 2. Rforest declines ten small primes. Legacy supplies
five of these (17, 19, 23, 29, 31); both paths decline 3, 5, 7, 11, 13.
Hybrid `rforest` retains native fallback, and Python retains Sage fallback.
The saved million-bound input hash was verified before comparison.

The million-bound oracle is the saved October 29 legacy output. There is
**no fresh million-bound legacy timing** here. The earlier four-worker
751.86-second legacy result is historical and is not used to calculate
the speedup above. Worker count and machine load differ. The older
881.97-second, 13-GiB custom-tree result is also historical.

All timing rows cover native startup, Frobenius computation and output,
excluding Sage fallback, characteristic-polynomial lifting and L-value
evaluation. They do not measure an end-to-end L-value speedup. Strict
`rforest-only` prohibits fallback, so correctness cannot be hidden by
calling the existing implementation. The few extra legacy small-prime jobs
have negligible cost but their coverage difference is explicit above.

This was a shared 16-logical-CPU Intel Core Ultra 7 265H Linux machine.
CPU frequency, host load and scheduling vary; these are not dedicated
hardware measurements. Both paths are in the same executable and use
GCC/G++ 15.2.0, private thread-safe static NTL 11.6.0 and system GMP.
C++ uses `-O2 -std=c++17`; rforest C uses
`-O3 -std=gnu11 -fomit-frame-pointer -funroll-loops`.
An 8-GiB address-space limit applied to all measured native processes.
The unchanged legacy arithmetic objects were checked against the main build.
The [raw report](frobenius_rforest_20261101.json) includes all repetitions,
CPU times, memory, input/binary/source hashes, and exploratory tuning runs.

## What made the implementation faster

- Vendored rforest, pinned to upstream commit
  `3103d396c67cb1685131b1f11e84975cca335bdf`, supplies GMP remainder forests
  and Harvey's FFT integer matrix arithmetic.
- Exact common polynomial factors and integer contents are cancelled
  before evaluation. The largest precision-two transition degree on the
  test curve falls from 18 to 6. Exact remaining denominator valuations
  determine the additional working precision.
- Exact forward differences evaluate consecutive polynomial matrices by
  addition. A four-worker 100,000-bound comparison gave 29.95 seconds
  versus 33.45 for the original evaluator, with identical matrices.
- FFT plans and controls are local to each worker and reused across
  numerator/denominator and reduction-pair calls. Diagnostic shared memory
  counters are disabled. A leaked temporary vector in upstream `mproduct`
  is freed.
- The automatic forest split grows with `log(log(bound))`, keeping the
  number of subtrees logarithmic rather than fixing their leaf count.
  Exploratory four-worker tests of kappa 0 through 9 favoured kappa 6 at
  100,000. Disabling FFT increased that run from 32.94 to 49.80 seconds.
- Larger reduction pairs run first. A 100,000-bound worker sweep gave
  55.78, 25.26, 19.69 and 16.15 seconds at 2, 4, 8 and 16 workers,
  respectively; all matrices matched. Single exploratory timings guided
  tuning; the repeated paired measurements supply the speedup claim.

The algorithm still uses more memory than legacy, and the best worker and
forest settings depend on curve, precision and hardware. No global optimum
or universal crossover is claimed. The original paper's average polynomial
time bound describes scaling, whereas the forest, FFT, cancellation and
scheduling address the practical constants that limited the earlier tree.

## Validation and reproduction

The targeted suite passed **123 tests**, covering full matrices for genera
1--3 and precisions 1--4, independent Sage Euler factors, both L-function
coefficient paths, zero constant terms, negative coefficients, transformed
even models, mixed/duplicate/unsorted jobs, small and bad primes, subtree
boundaries, both evaluators and FFT controls. The executable self-check
also passed with 1, 4 and 16 workers. Source-distribution inspection verified
all vendored C sources, headers and license notices and excluded build products.
An AddressSanitizer build with leak detection passed four additional native
cases (including genus three, zero constant term, and 16 workers), comparing
3,006 full matrices against legacy with no memory errors or reported leaks.

```sh
# Run from the experimental worktree with the native executable built.
sage -python -m pytest tests/test_average_frobenius.py \
    tests/test_lfunction.py tests/test_native_build.py -q
sage -python -m hyperell_regulator.frobenius --check \
    --binary hypellfrob-threaded/build/euler

# Fresh paired native runs, including matrix comparison.
sage -python hypellfrob-threaded/benchmark_rforest.py \
    --binary hypellfrob-threaded/build/euler --bound 100000 \
    --workers 16 --repeats 3 --output /tmp/rforest-paired.json

# Fresh million-bound legacy/rforest comparison; takes longer.
sage -python hypellfrob-threaded/benchmark_rforest.py \
    --binary hypellfrob-threaded/build/euler --bound 1000000 \
    --workers 16 --output /tmp/rforest-million.json
```

The last command computes a new oracle rather than requiring local saved
output. `--reference` together with `--reference-report` instead uses a
saved legacy output and verifies its original input hash, as in the
million-bound result above. See the [native guide](../hypellfrob-threaded/README.md#rforest-backend)
for build instructions, strict/hybrid modes and tuning controls.
