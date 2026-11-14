# hyperell-frobenius

Full Frobenius matrices of hyperelliptic curves over `F_p`, computed for
many primes simultaneously using the bundled [rforest C library](rforest/README.hyperell.md).

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
build/euler [nthreads] [rforest [kappa]]
```

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
`HYPELLFROB_THREADS` chooses workers unless `nthreads` is supplied.

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

The Python wrapper searches for `euler` at `$HYPELLFROB_EULER`, then
`hyperell_regulator/bin/euler`, then `hyperell-frobenius/build/euler`.
