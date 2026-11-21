# hyperell-regulator

This repository contains a SageMath package for numerical tests of Beilinson's conjecture for the
standard motive of genus-two curves over the rationals. This is the rank-5
submotive of the wedge-square of h^1(C) given by the kernel of the cup
product. The package computes regulator
integrals, periods, standard L-values, and reduction data for hyperelliptic curves.
Regulators, periods and L-values can be computed to arbitrary precision.
The complete documentation can be found [here](https://lac1213.github.io/hyperell-regulator/).

The [interactive cut-system explorer](https://lac1213.github.io/hyperell-regulator/explorer/)
shows one two-sheet branch diagram per curve and its logarithmic integrals for some genus-two curves of small conductor.

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

## Beilinson experiment

Search for pairs of integral motivic divisors and compare their regulator
determinants with the standard L-value.

```sh
# Known rational points on cached curves below conductor 1000
sage -python -m hyperell_regulator.experiments.beilinson 1000 --digits 20 --reg-digits 20

# Search the geometric torsion packet when the rational packet has fewer than two integral directions
sage -python -m hyperell_regulator.experiments.beilinson 1000 --packet geometric --digits 20 --reg-digits 20

# Certified regulator and conditional L-value enclosures
sage -python -m hyperell_regulator.experiments.beilinson 250 --certified --digits 100 --edge-seconds 600
```

The survey selects LMFDB rational endomorphism algebras `Q` and `CM`.
Use `--fetch` to refresh the curve cache and `--packet weierstrass` to
use the six geometric Weierstrass points. See the
[experiment guide](docs/experiments.rst) for precision and cache options.
Data are cached in `lmfdb_cache/` and `results/` under the working directory
or under `HYPERELL_REGULATOR_DATA`.

The [saved surveys](reports/README.md) give numerical candidates,
integration models, and integrality checks.

## Package examples

The numerical integrator computes periods and integrals weighted by `log |x-a|`. A minimum spanning tree joins the branch values and logarithmic singularities and the two-dimensional regulator integral is computed as a bilinear expression of line integrals along the edges of the tree.
For example you can run the following in Sage or with `sage -python`:

```python
from sage.all import QQ, PolynomialRing
from hyperell_regulator.regulator import (
    periods_and_regulator_stokes, regulators_stokes, c_invariants,
)

x = PolynomialRing(QQ, "x").gen()
F = x**5 + (x + 1)**2              # curve y^2 = F(x)
periods, regulator = periods_and_regulator_stokes(F, a=0, prec=30)
periods, columns = regulators_stokes(F, [-1, 0], prec=30)
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

`TorsionPacketIntegrality(C, points)` checks the sufficient multidegree-zero
criterion for `q_* B(sum(a_i P_i))` and computes integer kernels of linear
combinations. See the [integrality guide](docs/integrality.rst) for point
conventions, unresolved reductions, and the archived multidegree survey.

`ProjectedTorsionPacketIntegrality(C, points)` solves the projected
local boundary constraints and returns a basis of integral combinations.
It also accepts algebraic points over an absolute number field.
[Projected boundary examples](docs/integrality.rst#projected-boundaries-and-integral-classes)

`MotivicPacket(C)` gives the six geometric Weierstrass points and accepts
degree-zero divisors by coefficient row or point dictionary. Use
`D.determinant(E, digits=20)` for a pair of divisor classes,
`D.boundary(prime)` for the local graph pairing, and
`packet.integral_basis()` for integral divisor representatives.
See the [divisor and stable model API](docs/integrality.rst#divisor-classes-and-stable-models).

`geometric_torsion_packet(C, basepoint=O)` computes the complete geometric
hyperelliptic torsion packet using Poonen's external PARI/GP backend. See
the [geometric torsion guide](docs/geometric-torsion.rst) for installation,
supported base points, exact algebraic coordinates, and new classes on
249.a.6723.1.

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

Supply `certified=True` to the L-value and regulator/period methods to return
Sage Arb balls with a checked absolute radius below `10^-digits` (`prec` for
regulators). Certified L-values are conditional on the analytic
continuation and functional-equation of the L-function. Build their optional native
backend with `sage -python -m hyperell_regulator.build_certified`.
See the [ball arithmetic guide](docs/ball-arithmetic.rst) for requirements and
supported evaluations.

## Tests

Install the `test` extra as shown above, then run:

```sh
sage -python -m pytest
sage -t hyperell_regulator
```

Build and preview the documentation locally with the package installed.

```sh
sage -pip install -r docs/requirements.txt
sage -python -m sphinx -W --keep-going -b html docs docs/_build/html
sage -python -m http.server 8000 --bind 127.0.0.1 --directory docs/_build/html
```
