# hyperell-regulator

A SageMath package for numerical tests of Beilinson's conjecture for the
standard motive of genus 2 curves over the rationals. This is the rank-5
submotive of the wedge-square of h^1(C) given by the kernel of the cup
product. This packages computes regulator
integrals, periods, standard L-values, and reduction/component data.
Regulators, periods and L-values support both float64 and arbitrary precision.

[Documentation](https://lac1213.github.io/hyperell-regulator/)

The [interactive cut-system explorer](https://lac1213.github.io/hyperell-regulator/explorer/)
shows sixteen genus 2 curves of small conductor and thirty-two torsion functions, with
clickable lifted paths, ordinary/logarithmic integrals and regulator pairings
displayed to ten significant digits. Its data is computed locally with Sage;
the viewer runs entirely on GitHub Pages. See the
[export and regeneration guide](docs/explorer-guide.rst).

## Installation

Install into an existing SageMath environment:

```sh
sage -pip install git+https://github.com/LAC1213/hyperell-regulator.git
```

From a checkout, use `sage -pip install .`, or
`sage -pip install -e '.[test]'` for development and tests.
The patched `sage_cluster_pictures` dependency is installed automatically.

Installation downloads the [NTL 11.6.0 GitHub release](https://github.com/libntl/ntl/releases/tag/v11.6.0),
verifies its SHA-256, and builds a private
static copy with `NTL_THREADS=on` and `NTL_THREAD_BOOST=on` for the Frobenius
accelerator. Sage's installed NTL is untouched; GMP comes from Sage. This
requires network access, a C++17 compiler, Perl, and Make. Build parallelism
defaults to the logical CPU count (`HYPERELL_REGULATOR_BUILD_JOBS` overrides it).
Set `GMP_INC` and `GMP_LIB` if GMP is outside Sage's prefix.
If the build is unavailable, the package uses Sage's implementation.

Euler-factor computation defaults to one worker per logical CPU; use
`frobenius_polynomials(f, primes, nthreads=4)` to override the worker count.
The build checks serial, four-worker, and CPU-count results against each other.
To require a working accelerator during installation and then check it:

```sh
HYPERELL_REGULATOR_REQUIRE_EULER=1 sage -pip install -v git+https://github.com/LAC1213/hyperell-regulator.git
sage -python -m hyperell_regulator.frobenius --check
```

See [hypellfrob-threaded](hypellfrob-threaded/README.md) for manual build instructions.

## Python API

Run in Sage or with `sage -python`:

```python
from sage.all import QQ, PolynomialRing
from hyperell_regulator.regulator import (
    periods_and_regulator_stokes, c_invariants,
)

x = PolynomialRing(QQ, "x").gen()
F = x**5 + (x + 1)**2              # curve y^2 = F(x)
phi = (x + 1, 1)                  # function (x + 1) + y

periods, regulator = periods_and_regulator_stokes(F, phi)  # float64
periods, regulator = periods_and_regulator_stokes(F, phi, prec=30)
c_plus, c_minus = c_invariants(F, prec=30)
```

Use exact rational coefficients for arbitrary precision. `prec` is a decimal
accuracy target and convergence failures raise exceptions. 

For `c_invariants`, the finite branch points are ordered by increasing argument
about a perturbed centroid, starting at the root of largest absolute value
(ties by real, then imaginary part). Period paths share the regulator's
semilunar detours and embeddedness checks. The period matrix contains full
cycle integrals, twice the integrals along the lifted paths; `c_plus` and
`c_minus` are signed pivoted minors of its real and imaginary parts.
For genus two, using full periods multiplies each minor by four relative to
half-periods on the same paths. The new paths can also change signs and rational
normalisations. Existing saved reports use their original conventions;
experiment caches are invalidated automatically when rerun.

## Beilinson experiments

```sh
# Fetch curve data and test conductors below 1000; float64 by default.
sage -python -m hyperell_regulator.experiments.beilinson 1000 --fetch

# Reuse cached data and request 30-digit regulators and L-values.
sage -python -m hyperell_regulator.experiments.beilinson 250 --prec 30 --digits 30 --edge-seconds 120

# Test the parameter n=2 in the one-parameter family, without an L-value.
sage -python -m hyperell_regulator.experiments.family_n 2 --no-lvalue
```

Both experiments support `--prec`, `--digits`, `--no-lvalue`, `--verbose`,
`--edge-seconds`, `--force`, and `--out`; use `--help` for all options.
The default edge limit is 15 seconds. For `family_n`, `--max-conductor`
bounds the standard motive's conductor, rather than the curve conductor.

Reports compare `q = det / (2*pi*c_plus*c_minus*L''(std,1))` with a rational
guess. Data and results are cached under
`lmfdb_cache/` and `results/` in the working directory; set
`HYPERELL_REGULATOR_DATA` to choose another location.

The [saved 100-digit survey](reports/README.md) includes LaTeX tables,
full-precision values, exact functions, and regeneration instructions.

## Tests

Install the `test` extra as shown above, then run:

```sh
sage -python -m pytest
sage -t hyperell_regulator
```
