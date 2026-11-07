# hyperell-regulator

This repository contains a SageMath package for numerical tests of Beilinson's conjecture for the
standard motive of genus-two curves over the rationals. This is the rank-5
submotive of the wedge-square of h^1(C) given by the kernel of the cup
product. The package computes regulator
integrals, periods, standard L-values, and reduction data for hyperelliptic curves.
Regulators, periods and L-values can be computed to arbitrary precision.
The complete documentation can be found [here](https://lac1213.github.io/hyperell-regulator/).

The [interactive cut-system explorer](https://lac1213.github.io/hyperell-regulator/explorer/)
shows the paths of integration and functions used for some genus-two curves of small conductor.

## Installation

Install into an existing SageMath environment:

```sh
sage -pip install git+https://github.com/LAC1213/hyperell-regulator.git
```

From a checkout, use `sage -pip install .`, or
`sage -pip install -e '.[test]'` for development and tests.
The patched `sage_cluster_pictures` dependency is installed automatically.

Installation builds a C++ program for fast computation of
Frobenius eigenvalues which uses [rforest](https://math.mit.edu/~drew/rforest.html) and implements the algorithm of David Harvey and Andrew V. Sutherland.

## Examples

Given a function `phi` with divisor of the form `N P - N Q` on a genus-two curve `C`, we consider the integrals of `log |phi| w`, where w is a (1,1)-form on `C`.
For example you can run the following in Sage or with `sage -python`:

```python
from sage.all import QQ, PolynomialRing
from hyperell_regulator.regulator import (
    periods_and_regulator_stokes, c_invariants,
)

x = PolynomialRing(QQ, "x").gen()
F = x**5 + (x + 1)**2              # curve y^2 = F(x)
phi = (x + 1, 1)                  # function (x + 1) + y with divisor 5 (0, -1) - 5 oo

periods, regulator = periods_and_regulator_stokes(F, phi, prec=30)
c_plus, c_minus = c_invariants(F, prec=30)
```

`prec` is the number of decimal digits to compute.

`divisor_class_order` computes exact orders of divisors of the form `P - Q`:

```python
from hyperell_regulator.components import HyperellipticCurveWithReduction, divisor_class_order

C = HyperellipticCurveWithReduction(x + 2*x**2 - x**4, 1 + x + x**3)
divisor_class_order(C, (0, 0), ("infinity", 0))   # 7
divisor_class_order(C, (0, 0), ("infinity", -1))  # 21
```

Here infinity labels are the limits of `y/x**(g+1)` on the
`y**2 + h*y = f` model.

For the curve `y^2 + (x^3 + 1)y = x^2 + x` (LMFDB 249.a.249.1),
compute the standard L-function's second derivative `L''(std, 1)` to
30 decimal digits:

```python
from sage.all import QQ, PolynomialRing
from hyperell_regulator import HyperellipticCurveWithReduction, standard_lvalue

x = PolynomialRing(QQ, "x").gen()
C = HyperellipticCurveWithReduction(x**2 + x, x**3 + 1)
result = standard_lvalue(C, digits=30, verbose=False)
print(result["value"])
```

## Beilinson experiments

```sh
# Reuse cached data and request 100-digit regulators and L-values.
sage -python -m hyperell_regulator.experiments.beilinson 1300 --edge-seconds 600 --digits 100 --reg-digits 100

# Test the parameter n=2 in the one-parameter family.
sage -python -m hyperell_regulator.experiments.family_n 2
```

Both commands use float64 regulators and a 15-second edge limit by default.
Use `--help` or the [experiment guide](docs/experiments.rst) for precision, conductor, and cache options.

Data and results are cached under
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
