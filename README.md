# hyperell-regulator

Beilinson's conjecture for the square of a genus 2 curve.

For a genus 2 curve C over ℚ, we construct divisors D₁, D₂ supported on a
torsion packet, prove that the classes α(D₁), α(D₂) ∈ H³_𝓜(C², ℚ(2)) are
integral, and compute their regulator matrix R(D₁,D₂). We then check
numerically that

    π · c⁺(C) · c⁻(C) · L″(M_std(C), 1) ≈ q · det R(D₁,D₂)

for a small rational number q. Here M_std(C) is the rank-5 standard motive,
the kernel of the cup product on ∧²h¹(C), and c⁺(C), c⁻(C) are Deligne
periods. All quantities can be computed to arbitrary precision.

- [Interactive explorer](https://lac1213.github.io/hyperell-regulator/explorer/)
- [Documentation](https://lac1213.github.io/hyperell-regulator/)

## Installation

In a SageMath environment:

```sh
sage -pip install git+https://github.com/LAC1213/hyperell-regulator.git
```

From a checkout, use `sage -pip install -e '.[test]'`. Installation builds
a C++ program computing Frobenius polynomials with
[rforest](https://math.mit.edu/~drew/rforest.html) (Harvey–Sutherland).

## 1. The curve and its torsion packet

`load_curve` fetches a curve from the LMFDB, or reads it from the local cache.
Curves can also be built directly with `HyperellipticCurveWithReduction(f, h)`.

```python
from hyperell_regulator import load_curve, geometric_torsion_packet

C = load_curve("249.a.249.1")        # y^2 + (x^3 + 1)y = x^2 + x
T = geometric_torsion_packet(C)      # torsion packet containing the Weierstrass points

T.size                               # 10
T.points(rational_abscissae=True)    # members of the packet with rational x coordinate: [(-1, 0), ('infinity', 0), ('infinity', -1), (0, 0), (0, -1)]
T.order((0, 0))                      # 7
```

The packet is computed by Poonen's PARI/GP code.
Install it once with `install_poonen_backend()`. Points at infinity are
labelled by lim y/x³. See the [geometric torsion guide](https://lac1213.github.io/hyperell-regulator/geometric-torsion.html).

## 2. Integral motivic classes

A degree-zero divisor D supported on a torsion packet gives a class α(D). We prove that
α(D) is integral by solving the Laplace equation on the dual graph of the
special fibre at each bad prime, as in section 1 of the explorer.
`integral_basis` returns a basis of the integral divisors over ℚ.

```python
from hyperell_regulator import MotivicPacket

P = MotivicPacket.from_geometric_packet(T, rational_abscissae=True)
basis, complete = P.integral_basis()
[D.coefficients for D in basis]      # [(1, 0, 0, -1/2, -1/2), (0, 1/2, 1/2, -1/2, -1/2)]
complete                             # True: the integral rank is exactly 2

P.integrality.check(basis[0].coefficients)["q_integral"]   # True
P.boundary(basis[0], 3).projected    # (0): no obstruction at p = 3
```

See the [integrality guide](https://lac1213.github.io/hyperell-regulator/integrality.html) for the local models.

## 3. The regulator

When the integral rank is at least two, we compute det R(D₁,D₂). As in
section 2 of the explorer, we cut C along the preimage of a tree joining the
critical values and logarithmic centers, which turns the regulator integral
into a bilinear expression of line integrals along the edges of the tree.

```python
from sage.all import pi
from hyperell_regulator import standard_lvalue
from hyperell_regulator.regulator import c_invariants

D1, D2 = (P.divisor([2*a for a in D.coefficients]) for D in basis)   # clear denominators
det = P.determinant(D1, D2, digits=20)            # 133.95092605065823947547511

c_plus, c_minus = c_invariants(C.F, prec=20)      # -6.4459258435624592, -4.5685032931214575
L = standard_lvalue(C, digits=20, verbose=False)["value"]   # 5.0676270403470587136421

(pi * c_plus * c_minus * L / det).n()             # 3.50000000000000
```

So q = 7/2 for this curve.
Pass `certified=True` for Arb balls with radius below 10⁻ᵈⁱᵍⁱᵗˢ; certified
L-values assume the functional equation
([ball arithmetic guide](https://lac1213.github.io/hyperell-regulator/ball-arithmetic.html)).

## 4. All in one go

`beilinson_experiment` runs steps 1–3, searches for a nonzero determinant
and recognises q:

```python
from hyperell_regulator import beilinson_experiment

result = beilinson_experiment("249.a.249.1", digits=20, reg_digits=20)
result["comparison"]["rational"]     # -7/2
```

It first tries the rational points and falls back to the geometric packet
when they give fewer than two integral directions. Results are cached in
`lmfdb_cache/` and `results/` under the working directory or
`HYPERELL_REGULATOR_DATA`. See the [experiment guide](https://lac1213.github.io/hyperell-regulator/experiments.html).

## Console examples

One curve:

```console
$ sage -python -m hyperell_regulator.experiments.beilinson --label 249.a.249.1 --digits 20 --reg-digits 20
curve               cond status                      integral    determinant            q   rational rel. error  torsion divisors
249.a.249.1          249 ok                          integral     -133.95093         -3.5       -7/2 1.2468504e-19  D1 = 2*(-1, 0) - oo(Y=-1) - oo(Y=0) | D2 = (0, -1) + (0, 0) - oo(Y=-1) - oo(Y=0)
```

All cached curves of conductor below 300:

```console
$ sage -python -m hyperell_regulator.experiments.beilinson 300 --digits 20 --reg-digits 20
8 cached curves with conductor < 300; geometric packet, Galois traces over Q
curve               cond status                      integral    determinant            q   rational rel. error  torsion divisors
169.a.169.1          169 ok                          integral     -145.90597       -2.375      -19/8 7.0979887e-19  D1 = (-1, 0) + (-1, 1) - oo(Y=-1) - oo(Y=0) | D2 = (0, -1) + (0, 0) - oo(Y=-1) - oo(Y=0)
249.a.249.1          249 ok                          integral     -133.95093         -3.5       -7/2 1.2468504e-19  D1 = 2*(-1, 0) - oo(Y=-1) - oo(Y=0) | D2 = (0, -1) + (0, 0) - oo(Y=-1) - oo(Y=0)
249.a.6723.1         249 ok                          integral     -66.975463        -1.75       -7/4 5.9879455e-19  D1 = 2*(-1, 0) - oo(Y=-1) - oo(Y=0) | D2 = 4*(0, -2) + 4*(0, 1) - oo(Y=-1) - oo(Y=0) - 6*(3, -14)
  ***   Warning: unknown valuation of conductor at 2.
256.a.512.1          256 ok                          integral     -96.994895            5          5 2.973402e-18  D1 = (0, -1) + (0, 0) - (1, -1) - (1, 0) | D2 = (1/2, -1/2) - oo(Y=0)
277.a.277.1          277 ok                          integral     -143.44654        -3.75      -15/4 2.3135959e-19  D1 = 2*(-1, 0) - oo(Y=-1) - oo(Y=0) | D2 = (0, -1) + (0, 0) - oo(Y=-1) - oo(Y=0)
  277.a.277.2 Fewer than two integral directions found on any rational packet
  277.a.277.2 Loading or computing the Poonen torsion packet
  277.a.277.2 10 geometric torsion points in the cached packet
277.a.277.2          277 insufficient_integral_rank         -              -            -          -          -  -
  Geometric packet integral rank 1
295.a.295.1          295 ok                          integral      144.58555         -3.5       -7/2 2.5838663e-19  D1 = (0, -1) + (0, 0) - 2*(1, -1) | D2 = oo(Y=-1) - 2*(1, -1) + oo(Y=0)
  295.a.295.2 Fewer than two integral directions found on any rational packet
  295.a.295.2 Loading or computing the Poonen torsion packet
  295.a.295.2 6 geometric torsion points in the cached packet
295.a.295.2          295 insufficient_integral_rank         -              -            -          -          -  -
  Geometric packet integral rank 0

insufficient_integral_rank: 2; ok: 6
q = pi*c+*c-*L''(std,1) / D; rel. error = abs(q/rational - 1).
Results: results/beilinson_results.json
```

Here `determinant` is det R(D₁,D₂), `oo(Y=η)` is the point at infinity with
η = lim y/x³, and `D` in the last lines is the determinant. Add `--verbose`
for per-edge progress, `--certified` for rigorous enclosures, and `--fetch`
to refresh the LMFDB cache.

## Tests and documentation

```sh
sage -python -m pytest
sage -t hyperell_regulator

sage -pip install -r docs/requirements.txt
sage -python -m sphinx -W --keep-going -b html docs docs/_build/html
```
