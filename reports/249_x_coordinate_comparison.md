# Direct x-coordinate Stokes integration, conductor 249

This initial comparison uses semilunar cuts. The package now defaults to
minimum spanning trees; the verification driver now uses MST cuts. Original
timings below are historical. See the [paired MST comparison](mst_cut_comparison.md)
and [MST full-cycle periods](mst_periods.md).

The proposed simplification works. The numerical integrals need only the
degree-two hyperelliptic map. Functions with divisor `n(P-Q)` remain relevant
to constructing and certifying the algebraic cycles, but need not be used
as the covers for numerical integration. The package now uses this direct method throughout, with all logarithmic
weights sharing one two-sheet cut and no general chart machinery.

## Cut and Stokes formula

On `Y^2 = F(x)`, choose an embedded semilunar chain through all finite branch
values of x, **including 0 as a marked vertex**, starting at 0 and ending
with a ray to infinity. The zero and infinity are needed for the logarithm,
even when x is unramified there. Thus a chain connecting just the branch
points is not quite sufficient. For several `log(x-a)` weights, include every center a on the same cut.
Their logarithm jumps can differ from edge to edge.

The complement is simply connected and avoids the branch values. Its two
lifts are discs, and `log x` is holomorphic on each. Consequently

\[
g_a(z)=\int^z\log x\,\omega_a,\qquad
d(g_a\overline{\omega_b})=\log x\,\omega_a\wedge\overline{\omega_b}.
\]

Apply Stokes on both discs and glue the boundaries using the existing
Stokes assembly, retaining each edge's logarithm jump, a multiple of `2 pi i`. At a regular zero,
the two lifted endpoints are distinct and the first edge has identity
sheet gluing. At the two infinities of an even-degree model, each sheet
keeps its own infinity branch. Neither feature requires total ramification.

There is a convention issue with “take the real part.” Set

\[
W=\omega_1\wedge\overline{\omega_2}
 -\overline{\omega_1}\wedge\omega_2,\qquad \overline W=-W.
\]

For any chosen branch on the cut complement,

\[
\int_C\log|x|\,W
=i\,\operatorname{Re}\!\left(\frac1i\int_C\log x\,W\right).
\]

Thus take the real part **after dividing by i**, equivalently the imaginary
part of the complex integral and then multiply by i. The real part of the
unnormalized integral instead comes from `arg x`.

For `omega_1=dx/Y`, `omega_2=x dx/Y`, the package's returned regulator vector is

\[
r=\int_{\mathbb C}\frac{\log|x|}{|F(x)|}
 (1,\operatorname{Re}x,|x|^2)\,dA,
\qquad \int_C\log|x|\,W=-8i\,r_{12}.
\]

The determinant convention doubles the middle entry, so this is also
`-4i` times the middle entry of the `log|x|` matrix column.

## Numerical comparison

Both columns were computed together with `regulators_stokes(F, [-1, 0])`, at a
100-decimal-digit target, using the same integration coordinates and weights
as the saved high-degree runs. No torsion functions were constructed and no
L-values recomputed. The entire matrices, including the unweighted column,
were checked against the saved values.

| Curve | Saved cover degrees | Direct cover degrees | Largest absolute matrix difference | Direct time | Saved time |
|---|---:|---:|---:|---:|---:|
| 249.a.249.1 | 7, 14 | one shared degree-2 cut | 2.29e-100 | 8.78 s | 635.6 s |
| 249.a.6723.1 | 7, 28 | one shared degree-2 cut | 2.29e-100 | 10.41 s | 1041.0 s |

These are observed timings against previous runs, not a controlled benchmark.
The determinants and normalized quotients are

| Curve | D | pi c+ c- L'' / D |
|---|---:|---:|
| 249.a.249.1 | 33.48773151266455986457282727876349 | 28 (within 1.0e-99) |
| 249.a.6723.1 | 4.72042171469714946744437436356583 | 24.82978584930216841732078053849348 |

Full values, cut vertices, timings and error comparisons are in
[249_x_coordinate_comparison.json](249_x_coordinate_comparison.json).
The precision targets and comparisons are numerical, not interval-certified bounds.

For the first curve, the saved integration coordinate, denoted u in the
survey, has

\[
V^2=u^6+4u^5+4u^4+2u^3+1,
\]

and the direct `log|u|` integral of W is
`+3.33836038756391700965650136046851201029472208231958... i`.
The original LMFDB equation is `y^2+(x^3+1)y=x^2+x`.
Putting `Y=2y+x^3+1`, its polynomial is
`Y^2=x^6+2x^3+4x^2+4x+1`, and `u=1/x`, `V=u^3 Y`.
The two differentials pull back to `-u du/V` and `-du/V`, so W is preserved
but `log|x|=-log|u|`. In that **original coordinate**, the answer is therefore

\[
\boxed{\int_C\log|x|
 (\omega_1\wedge\overline{\omega_2}
 -\overline{\omega_1}\wedge\omega_2)
=-3.33836038756391700965650136046851201029472208231958\ldots\,i.}
\]

An additional direct 30-digit run on the original polynomial confirmed this
sign and value. Its regulator vector was
`(-1.5476867914348250215920716112, 0.41729504844548962620706267006,
1.1415824404631117209412054750)`.

## Reproduction

```sh
sage -python -m hyperell_regulator.experiments.x_coordinate_comparison
sage -python -m pytest -q tests/test_x_coordinate_stokes.py
```

The experiment reads the integral example from
`reports/beilinson_100_digits.json` and the second example from
`results/249_rational_approximations/input.json`.

For a single direct integral:

```python
from sage.all import QQ, PolynomialRing
from hyperell_regulator.regulator import periods_and_regulator_stokes

x = PolynomialRing(QQ, "x").gen()
F = (x**3 + 1)**2 + 4*(x**2 + x)  # original 249.a.249.1 coordinate
periods, r = periods_and_regulator_stokes(F, a=0, prec=100)
print(-8*r[1])  # imaginary coefficient of the requested integral
```
