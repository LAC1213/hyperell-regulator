# hyperell-frobenius

Full Frobenius matrices of hyperelliptic curves over `F_p`, computed for
many primes simultaneously using the bundled [rforest C library](rforest/README.hyperell.md).

The C++ implementation follows David Harvey, *Counting points on
hyperelliptic curves in average polynomial time*,
[arXiv:1210.8239](https://arxiv.org/abs/1210.8239), Sections 3–5.
It computes matrices modulo `p^N`, including the extra precision required
for denominator cancellation. The implementation extends the reduction
paths to even degree and restricts affine cohomology to the kernel of the
residue at infinity.

Rforest is the only native implementation. Unsupported jobs are declined;
`hyperell_regulator.frobenius` handles those primes with Sage. The earlier
per-prime implementations remain available in Git history.

The C++ code is distributed under GPL version 2 or later; see [COPYING](COPYING).
The bundled C library retains its own copyright and license notices.

## Building

```sh
make                       # if NTL is under /usr/local
make NTL_INC=$HOME/ntl/include NTL_LIB=$HOME/ntl/src NTL_LINK=$HOME/ntl/src/ntl.a
make check
```

Package installation downloads checksum-verified NTL 11.6.0, builds a private
copy with `NTL_THREADS=on` and `NTL_THREAD_BOOST=on`, and statically links it
into `hyperell_regulator/bin/euler`. Sage's NTL is untouched. GMP comes from
Sage's prefix; `GMP_INC` and `GMP_LIB` override its paths.
`HYPERELL_REGULATOR_BUILD_JOBS` controls build parallelism, which defaults
to the logical CPU count. Installation requires network access, Perl,
Make, a C11 compiler, and a C++17 compiler.

For manual builds, copy `config.mk.example` to `config.mk` and edit it.
Package builds ignore those local settings. Manual builds also require a
thread-safe NTL with matching headers, because NTL arithmetic runs across
independent reduction forests. The FFT code requires a 64-bit platform
with 64-bit GMP limbs.

Both `make check` and the installer compare full odd- and even-degree
matrices against fixed verified references using one, four, and CPU-count
workers. These checks run without Sage in build isolation. Additional
tests compare odd-degree matrices with Sage, all characteristic coefficients
modulo `p^N` with independent Sage polynomials, and L-function coefficients
with Sage and PARI.

## Using it

```text
build/euler [nthreads] [rforest|rforest-only] [kappa]
```

Both names select the same implementation and always decline unsupported
jobs. `rforest-only` remains an alias for existing native callers.
Removed backend names are rejected.

The driver reads coefficients of `Q` in ascending order, then jobs `p N`:

```text
1 0 0 1 0 1
101 3
103 3
```

It writes one row per job, in the input order:

```text
<p> 1 <dim> <entries...>      row-major 2g-by-2g matrix modulo p^N
<p> 0                         declined: bad reduction or an unmet precondition
```

## rforest backend

```python
frobenius_polynomials(f, primes, nthreads=4)
```

The Python API and both L-function classes use rforest automatically.
The optional `backend="rforest"` keyword and `HYPELLFROB_BACKEND=rforest`
setting remain accepted; other values raise an error.
`HYPELLFROB_THREADS` chooses workers unless `nthreads` is supplied.

Monic even degree models, including conductor 388, are handled directly
without a rational branch point or field extension. Shifted horizontal and
diagonal paths avoid singular reductions; the residue invariant is checked
before returning a projective matrix. Python monicizes nonmonic even models
using a common quadratic twist, then removes its local character from the
characteristic polynomial.

When a rational Weierstrass point exists, Python first moves it to infinity
to obtain a smaller odd degree model, isomorphic over `Q`. Primes where the
transformed model has bad reduction fall back to Sage on the original model.

The C library builds with the driver; no separate installation is needed.
FFT state is local to each worker and reused across reduction pairs.
Polynomial factor cancellation reduces transition degrees and the required
precision. Forward differences accelerate polynomial evaluation; larger
reduction pairs run first to balance workers.

The optional `kappa` argument controls subtree splitting: `-1` selects it
automatically; larger values generally reduce temporary tree storage at
the cost of more prefix updates. Choose workers for available hardware;
more concurrent forests consume more memory and can be slower.

## Correctness and timing

Check the installed driver:

```sh
sage -python -m hyperell_regulator.frobenius --check
```

Benchmark against independent Sage Euler factors at a modest bound:

```sh
sage -python hyperell-frobenius/benchmark_rforest.py \
    --binary hyperell-frobenius/build/euler --bound 1000 \
    --sage-reference --workers 4 --repeats 3 \
    --output /tmp/rforest-benchmark.json
```

The utility requires one oracle: `--sage-reference`, `--euler-reference file`
for a smalljac `lpdata` CSV, or `--reference file` for saved full native
matrices. `--reference-report` can verify the input hash of a saved matrix
reference. No smalljac dependency is required to build or use the backend.

For large bounds such as `1000000`, use an existing reference file to avoid
recomputing independent factors. Reported wall time, CPU time, and peak RSS
cover native computation only; reference generation, polynomial lifting,
and L-function evaluation are excluded. Every supplied matrix or Euler
factor is compared. Unsupported jobs are counted as declined.
`--evaluator powers`, `--no-fft`, and `--kappas` expose tuning controls.

The earlier custom accumulating-tree implementation is archived on
`experimental/harvey-remainder-tree` and `experimental/rforest`. Historical
measurements and the removed per-prime sources remain in Git history.

The Python wrapper searches for `euler` at `$HYPELLFROB_EULER`, then
`hyperell_regulator/bin/euler`, then `hyperell-frobenius/build/euler`.
