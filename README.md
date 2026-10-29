# hyperell-regulator

A SageMath package for numerical tests of Beilinson's conjecture for the
standard motive of genus 2 curves over the rationals. This is the rank-5
submotive of the wedge-square of h^1(C) given by the kernel of the cup
product. The package computes regulator
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

Installation builds an optional Frobenius accelerator using a private,
checksum-verified NTL 11.6.0. It needs network access, a C++17 compiler, Perl,
and Make; if unavailable, computations fall back to Sage.
See the [installation guide](docs/installation.rst) for build settings and
[manual native build instructions](hypellfrob-threaded/README.md).

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

`c_invariants` uses signed minors of **full cycle periods** along validated
semilunar paths. See the [period conventions](docs/quickstart.rst#period-conventions)
when comparing saved reports or other normalizations.

## Beilinson experiments

```sh
# Fetch curve data and test conductors below 1000; float64 by default.
sage -python -m hyperell_regulator.experiments.beilinson 1000 --fetch

# Reuse cached data and request 30-digit regulators and L-values.
sage -python -m hyperell_regulator.experiments.beilinson 250 --prec 30 --digits 30 --edge-seconds 120

# Test the parameter n=2 in the one-parameter family, without an L-value.
sage -python -m hyperell_regulator.experiments.family_n 2 --no-lvalue
```

Both commands use float64 regulators and a 15-second edge limit by default.
Use `--help` or the [experiment guide](docs/experiments.rst) for precision,
budget, conductor, and cache options.

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
