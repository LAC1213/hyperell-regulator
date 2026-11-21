# hyperell-frobenius

Frobenius matrices of hyperelliptic curves over `F_p`, computed simultaneously
for many primes using [rforest](rforest/README.hyperell.md).

## Build

```sh
make
make check
```

Manual builds require thread-safe NTL headers and libraries. Copy
`config.mk.example` to `config.mk`, or pass paths explicitly:

```sh
make NTL_INC=$HOME/ntl/include NTL_LIB=$HOME/ntl/src NTL_LINK=$HOME/ntl/src/ntl.a
```

Package installation builds a private NTL 11.6.0 with `NTL_THREADS=on` and
`NTL_THREAD_BOOST=on`, statically linked into `hyperell_regulator/bin/euler`.
GMP comes from Sage; `GMP_INC` and `GMP_LIB` override its paths.
`HYPERELL_REGULATOR_BUILD_JOBS` sets build parallelism, defaulting to the
logical CPU count. The build requires network access, Perl, Make, C11 and
C++17 compilers, and a 64-bit platform with 64-bit GMP limbs.

`make check` and package installation compare odd and even degree matrices
with reference values using one, four and CPU-count workers.

## Driver

```text
build/euler [nthreads] [rforest [kappa]]
```

Input consists of the coefficients of `Q` in ascending degree, followed by
jobs `p N`:

```text
1 0 0 1 0 1
101 3
103 3
```

Output has one row per job, in input order:

```text
<p> 1 <dim> <entries...>    row-major 2g-by-2g matrix modulo p^N
<p> 0                       declined job
```

`kappa=-1` chooses subtree splitting automatically. Larger values reduce
tree storage at the cost of more prefix updates. Sparse requests (including
missing high primes after a cached range) use balanced modular products for
long gaps, rather than constructing one growing exact integer product.

`make check-arithmetic` compares FFT matrix multiplication and sparse
remainder forests with independent, step-by-step GMP arithmetic.

## Python

```python
from hyperell_regulator.frobenius import frobenius_polynomials

frobenius_polynomials(f, primes, nthreads=4)
```

`HYPELLFROB_THREADS` sets workers unless `nthreads` is supplied.
The wrapper searches for `euler` at `$HYPELLFROB_EULER`, then
`hyperell_regulator/bin/euler`, then `hyperell-frobenius/build/euler`.
Unsupported jobs use Sage.

Check the installed executable or benchmark against Sage:

```sh
sage -python -m hyperell_regulator.frobenius --check
sage -python hyperell-frobenius/benchmark_rforest.py \
    --binary hyperell-frobenius/build/euler --bound 1000 \
    --sage-reference --workers 4 --repeats 3 \
    --output /tmp/rforest-benchmark.json
```
