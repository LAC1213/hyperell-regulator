# hyperell-frobenius

Frobenius matrices of hyperelliptic curves over `F_p`, for one curve at many
primes, batched and run across threads.

This is **David Harvey's `hypellfrob`**. The mathematics and essentially all
of the code are his:

> Copyright (C) 2007, 2008, David Harvey
> `hypellfrob` version 2.1.1, as distributed with Sage in
> `sage/schemes/hyperelliptic_curves/hypellfrob`
>
> David Harvey, *Kedlaya's algorithm in larger characteristic*,
> International Mathematics Research Notices, 2007.
> [arXiv:math/0610973](https://arxiv.org/abs/math/0610973)

`hypellfrob` is free software under the **GNU General Public License,
version 2 or later**; the full text is in [COPYING](COPYING), and this
package is distributed under the same terms. Every file keeps its original
copyright header, and the files in `src/` carry notices saying whether and
how they were modified, as the GPL requires.

The copies here came from a working tree of `hypellfrob-fast`, which applies
three performance changes to the dyadic-evaluation recursion in
`recurrences_ntl.cpp`; its output is bit-identical to Harvey's original.
`even_degree.*` is not Harvey's: it is a direct implementation of Kedlaya's
algorithm, linear in `p` rather than `O(sqrt p)`, which handles the even
degree models `matrix()` rejects.

## Why it is here

An `L`-function needs the Euler factor at every good prime below a bound,
which for a genus 2 curve is thousands of `p`-adic Frobenius computations.
Sage will do them one at a time; `euler` takes the whole list as a single
job, so the setup is paid once and the primes run concurrently. On this
package's test curve that is about 5.6x faster at 2000 primes, and 2258
primes below 20000 take under four seconds.

Nothing here is required. `hyperell_regulator.frobenius` falls back to Sage
for any prime this program does not supply, so results never depend on
whether it is built -- only the speed does.

## Building

```
make                       # if NTL is under /usr/local
make NTL_INC=$HOME/ntl/include NTL_LIB=$HOME/ntl/src NTL_LINK=$HOME/ntl/src/ntl.a
make check
```

The package installer downloads the pinned, checksum-verified NTL 11.6.0 GitHub release,
builds it with `NTL_THREADS=on` and `NTL_THREAD_BOOST=on`, and statically links
it into `hyperell_regulator/bin/euler`. It discovers Sage's GMP prefix
(including `SAGE_LOCAL`); `GMP_INC`/`GMP_LIB` override those paths.
`HYPERELL_REGULATOR_BUILD_JOBS` sets build parallelism (default: logical CPU count).
Installation needs network access, Perl, Make, a C11 compiler, and a C++17 compiler.
Sage's NTL is never changed. Package builds ignore local `config.mk` and
`NTL_INC`/`NTL_LIB`/`NTL_LINK`; these remain available for manual builds.
The build compares serial, four-worker, and CPU-count results before installing.

For manual builds, copy `config.mk.example` to `config.mk` and edit it. Threads need NTL
built with `NTL_THREADS` and `NTL_THREAD_BOOST`; without them pass `1` for
the thread count, which still batches.

`make check` only shows that the driver runs. What checks the arithmetic is
`tests/test_lfunction.py`, which compares every Euler factor against Sage's
own, prime by prime.

## Using it

```
build/euler [nthreads] [rforest|rforest-only|legacy] [kappa]
```

reads from stdin

```
1 0 0 1 0 1        # Q's coefficients, ascending: Q = x^5 + x^3 + 1
11 3               # one line per job: p and the p-adic precision N
13 3
```

and writes one line per job, in order:

```
<p> 1 <dim> <entries...>      row-major, the 2g x 2g matrix mod p^N
<p> 0                         declined: bad reduction, or a precondition
                              of the fast path failed
```

In `legacy` mode, odd degree goes through `hypellfrob::matrix()`, which is `O(sqrt p)` and
needs `Q` monic with `p > (2g+1)(2N-1)`; even degree through
`matrix_reference()`, which is correct but linear in `p`.

## rforest backend

The default `rforest` backend uses the actual
[rforest C library](https://math.mit.edu/~drew/rforest.html) for the
accumulating products in Harvey's 2012 algorithm. It computes full
Frobenius matrices modulo `p^N`, including the extra precision needed for
denominator cancellation. It is not restricted to Hasse--Witt matrices
modulo `p`.

```python
frobenius_polynomials(f, primes, nthreads=4, backend="rforest")
```

`rforest` is the default for the Python API, both L-function classes, and
the native driver. Select `backend="legacy"`, set
`HYPELLFROB_BACKEND=legacy`, or run `build/euler 4 legacy` to use the
per-prime algorithm explicitly. The default `build/euler 4` uses legacy
fallback for unsupported jobs; `rforest-only` declines those
jobs so arithmetic comparisons cannot pass through fallback.
Rforest now handles even degree models directly, including conductor 388,
without a rational Weierstrass point. Shifted horizontal/diagonal paths
avoid the even-degree singularity, and a residue check removes the extra
affine class. The result is a full matrix in the same projective basis as
the existing reference.
Python also monicizes nonmonic even input using a common quadratic twist
and removes its local character from the characteristic polynomial.

`HYPELLFROB_THREADS` sets the worker count for Python and both L-function
classes when no explicit `nthreads` is supplied. Eight workers with the
automatic forest split was fastest in the conductor-388 100,000-bound
sweep on the test machine. Choose workers for the available hardware;
more simultaneous trees use more memory and can be slower.

The pinned C sources and their licenses are vendored under
[`rforest`](rforest/README.hyperell.md). They build with the
driver, without a separate rforest installation. FFT state is local to
each worker and reused across reduction pairs. Exact polynomial factor
cancellation reduces both matrix degree and required precision; forward
differences replace repeated polynomial multiplications during evaluation.
The largest transition degree on `y^2 = x^5 + x^3 + 1` at precision two
drops from 18 to 6. Larger pairs run first to balance the workers.

The forest splits its work into subtrees to bound temporary storage.
The optional native `kappa` argument controls this split; `-1` chooses
automatically, and larger values generally use less tree memory but more
prefix updates. The automatic setting grows with `log(log(bound))`,
preserving the average polynomial time bound. A thread-safe NTL build
with matching headers is required, as for the other native backends.
Rforest's FFT implementation requires 64-bit limbs and a 64-bit platform.

Reproduce timings and full-matrix comparisons with:

```sh
sage -python hyperell-frobenius/benchmark_rforest.py \
    --binary hyperell-frobenius/build/euler --bound 100000 \
    --workers 4 --repeats 3 --output /tmp/rforest-benchmark.json
```

Use `--bound 1000000` for the million-bound experiment. The script
alternates run order, measures wall time, CPU time and peak RSS, and
compares every supplied matrix entry in strict mode. Timings cover the
native Frobenius computation; they exclude Sage characteristic-polynomial
lifting and L-function evaluation. `--reference` and `--reference-report`
allow comparison against a previously saved legacy run, without claiming
a fresh legacy timing. `--evaluator powers` and `--no-fft` expose the
evaluation and FFT controls used in performance tuning.
`--euler-reference file` instead accepts a smalljac `lpdata` CSV and
compares every supplied Euler factor independently. This works for even
models too, and omits the legacy run; it does not claim full matrix-entry
comparison against that oracle. No smalljac dependency is needed to build
or use the backend.
`--hybrid` includes the legacy fallback jobs when timing `rforest`; with a
legacy matrix reference, it requires identical output for every job.

## Correctness and archived comparisons

`sage -python -m hyperell_regulator.frobenius --check` compares the legacy
and strict rforest matrices with one, four, and CPU-count workers, for both
odd and even degree. Arithmetic and L-function tests live in
`tests/test_average_frobenius.py` and `tests/test_lfunction.py`.

The earlier custom accumulating-tree backend and its benchmark driver
remain on `experimental/harvey-remainder-tree` and `experimental/rforest`.
They are not included in the main-branch implementation. Historical
measurements remain in Git history. The current comparison utility is
`benchmark_rforest.py`.

The caller only sends an even degree model when it has to. `4f + h^2` has
even degree whenever `h` has degree `g+1`, but if the curve has a rational
Weierstrass point `r`, then `u^{\deg F} F(r + 1/u)` moves it to infinity and
has odd degree -- an isomorphism over `Q`, so the Frobenius characteristic
polynomial is unchanged -- and the fast path applies after monicising. A
curve with no rational Weierstrass point uses the even-degree forest in
`rforest` mode and the reference routine in `legacy` mode.

`hyperell_regulator.frobenius` drives this and lifts the matrices to
characteristic polynomials; it looks for the binary at `$HYPELLFROB_EULER`,
then at `hyperell_regulator/bin/euler`, then at
`hyperell-frobenius/build/euler` beside the package.
