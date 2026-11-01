# hyperell-regulator

A SageMath package for numerical tests of Beilinson's conjecture for the
standard motive of genus 2 curves over the rationals. This is the rank-5
submotive of the wedge-square of h^1(C) given by the kernel of the cup
product. The package computes regulator
integrals, periods, standard L-values, and reduction/component data.
Regulators, periods and L-values support both float64 and arbitrary precision.

[Documentation](https://lac1213.github.io/hyperell-regulator/)

The [interactive cut-system explorer](https://lac1213.github.io/hyperell-regulator/explorer/)
shows the successful 100-digit survey curves and their torsion functions, with
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
checksum-verified NTL 11.6.0. It needs network access, C11 and C++17 compilers,
Perl, and Make; if unavailable, computations fall back to Sage.
See the [installation guide](docs/installation.rst) for build settings and
[manual native build instructions](hyperell-frobenius/README.md).

The native accelerator uses rforest for simultaneous Frobenius computations,
including even degree models such as conductor 388. It is the only native
backend; unsupported primes fall back to Sage. No backend setting is needed.
`HYPELLFROB_THREADS` controls workers for L-function computations; see the
[backend guide](hyperell-frobenius/README.md#rforest-backend) for tuning and
comparison commands. The bundled C library is in `hyperell-frobenius/rforest`.

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

`divisor_class_order` computes exact orders on odd and even degree models,
including curves without a rational Weierstrass point:

```python
from hyperell_regulator.components import HyperellipticCurveWithReduction, divisor_class_order

C = HyperellipticCurveWithReduction(x + 2*x**2 - x**4, 1 + x + x**3)
divisor_class_order(C, (0, 0), ("infinity", 0))   # 7
divisor_class_order(C, (0, 0), ("infinity", -1))  # 21
```

Here infinity labels are the limits of `y/x**(g+1)` on the original
`y**2 + h*y = f` model. The polynomial API `torsion_order(F, P, Q)` uses
coordinates on `y**2 = F`; its optional `multiple=N` finds the exact order
dividing `N` without the default 24-step search bound.

## Beilinson experiments

```sh
# Reuse cached data and request 100-digit regulators and L-values.
sage -python -m hyperell_regulator.experiments.beilinson 1300 --edge-seconds 600 --digits 100 --reg-digits 100

# Test the parameter n=2 in the one-parameter family, without an L-value.
sage -python -m hyperell_regulator.experiments.family_n 2 --no-lvalue
```

Both commands use float64 regulators and a 15-second edge limit by default.
Use `--help` or the [experiment guide](docs/experiments.rst) for precision,
budget, conductor, and cache options.

Reports use determinant weights `(1, 2 Re(z), |z|^2)` and compare
`q = pi*c_plus*c_minus*L''(std,1)/D` with a rational guess. Data and results are cached under
`lmfdb_cache/` and `results/` in the working directory; set
`HYPERELL_REGULATOR_DATA` to choose another location.

The [saved 100-digit survey for conductor < 1300](reports/README.md) includes LaTeX tables,
full-precision values, exact functions, and regeneration instructions.

## Tests

Install the `test` extra as shown above, then run:

```sh
sage -python -m pytest
sage -t hyperell_regulator
```
