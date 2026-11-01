# Even-degree simultaneous Frobenius reductions

The `rforest` backend computes full p-adic Frobenius matrices for monic
squarefree `Q` of degree `d = 2g+2`, using the same integer reduction
identities as the odd-degree implementation. No branch point or field
extension is required. The `remainder` backend retains its original
odd-degree restriction; `legacy` retains the direct even-degree reference.

## Paths avoiding the horizontal singularity

Write `W(s,t)` for the span of `x^(s+i) dx/y^(2t+1)`, `0 <= i < d`.
The horizontal (H), diagonal (D), and vertical (V) maps lower `(s,t)` by
`(1,0)`, `(1,1)`, and `(0,1)`, respectively. H has denominator

```
d*(2t-1) - 2s = -2R,   R = s - d*(t-1/2).
```

In even degree, `R=0` can occur at integer indices. Applying H there would
multiply the scalar denominator by zero. D and V have no such singularity;
V requires `Q(0)` to be a unit at the prime in question.

For a reduction pair `(a,b)`, with `b` odd, block `r` begins at

```
s_r = a*(2r+1)-1,   t_r = (b*(2r+1)-1)/2.
C = d*b - 2*a,
R_r = -C*r + d-1-C/2.
```

It ends at the corresponding state for `r-1`. The implementation chooses:

- If `2a < b`, take `2a` D steps and `b-2a` V steps. There are no H steps.
- If `2a > b` and `C <= 0`, take `b` D steps followed by `2a-b` H steps.
  Here R stays strictly positive throughout H, including when `C=0`.
- If `2a > b` and `C > 0`, take `2a-b` H steps followed by `b` D steps.
  Start this sequence at `r > r0`, where
  `r0 = max(0, floor((d-1-C/2)/C))`. Then R is strictly negative throughout H.

These inequalities hold for every subsequent integer r. As `C` is a
positive even integer in the last case, `r0 <= g`. Instead of multiplying
the omitted singular blocks, construct a direct initial map from
`W(s_r0,t_r0)` to `W(-1,0)`: use V only while `t > s+1`, then D until
`t=0`, then H until `s=-1`. The last H denominators are `-d-2s`, which
cannot vanish. The first row in the final `W(-1,0)` coordinates is zero.

The bounded initial map is exact over Q, with integer content cancelled.
Primes dividing its denominator are declined; the usual native/Sage
fallback remains available. Every required pair must also satisfy the
existing p-admissibility checks, including `p > d*b+2a`, good reduction,
and the V-step unit condition. Those checks ensure that all forest
endpoints lie beyond their bounded prefixes.

Each generic block remains a polynomial matrix in r of degree at most
`max(2a,b)`. Exact common polynomial factors with its scalar denominator
are cancelled. Only nonconstant affine factors are retained for valuation
counting, so the `C=0` case never attempts an inverse of a zero slope.
The remaining factors' roots are counted in
`r0 < r <= (p-1)/2`, and the matrix/scalar forests start at `r0+1` with the
direct initial map. They run modulo `p^(mu+rho)` and cancel the denominator's
p-part exactly before inverting its unit part. The Frobenius expansion
carries the same explicit factor p as in the odd-degree case.

This preserves a polynomial transition independent of p and a bounded
initial map. In particular, the extension preserves sharing between primes
and does not replace the forest by separate per-prime extension-field jobs.
The Frobenius expansion uses the existing `p > d*(2mu-1)` precision bound;
the even-degree cohomology and precision analysis are compatible with the
standard `dx/y` basis at these primes.

## Removing the extra affine class

Before restriction, the affine minus-space has dimension `d-1 = 2g+1`, in
the monomial `x^i dx/y` basis. Let

```
P(u) = u^d Q(1/u),
rho_i = [u^(i-g)] P(u)^(-1/2).
```

Then `rho_i=0` for `i<g`, `rho_g=1`, and the smooth projective cohomology
is the kernel of this residue covector. The implementation checks

```
rho * M = p * rho   (mod p^mu)
```

for every affine Frobenius matrix. With `i,j != g`, it then returns

```
M_projective[i,j] = M[i,j] - rho_j*M[i,g].
```

This is the action on `e_j-rho_j*e_g`, expressed by dropping coordinate g.
It is the same projective basis as the existing even-degree reference
implementation, allowing full matrix-entry comparisons at all tested
precisions. It removes the extra eigenvalue p rather than assuming the
entire affine matrix describes the L-function.

## Nonmonic input through the Python API

For an integral even-degree polynomial f with leading coefficient c, the
common model

```
Q(X) = c^(d-1) * f(X/c)
```

is monic and integral. It is isomorphic to the c-twist of the input curve.
The Python caller removes this twist from each lifted characteristic
polynomial using the local character `(c/p)`. Primes dividing c, or
singular on the transformed model, use the original-model Sage fallback.
Thus the transformation also preserves shared integer coefficients across
primes. The native matrix interface itself continues to require monic Q.

## Why a quadratic branch extension is not used

A root over the unramified quadratic extension is available only when
Q modulo p has a factor of degree one or two. For conductor 388,
`x^6-2*x^4+2*x^3+9*x^2+6*x+1` is irreducible modulo 11, so it has no root
over `F_(11^2)` or the unramified quadratic p-adic extension. Moreover,
choosing a different p-adic root at each prime changes the lifted
polynomial, removing the common integer transition used by the forest.
The direct even-degree paths avoid both issues.

## References and validation

- David Harvey, [Counting points on hyperelliptic curves in average
  polynomial time](https://arxiv.org/abs/1210.8239), for the shared
  Frobenius expansion, integer maps, and accumulating products.
- Michael Harrison, [An extension of Kedlaya's algorithm for hyperelliptic
  curves](https://arxiv.org/abs/1006.4206), for the even-degree affine
  cohomology, integrality analysis, and extra eigenvalue.

The path selection, shifted forest start, and direct bounded prefix above
are the even-degree adaptation implemented here. Tests compare full
matrices at precisions 1--4 in genera 1--3 with the existing reference;
independent Sage/PARI Euler factors, the residue identity, nonmonic twists,
subtree boundaries, mixed precision, and parallel runs supply further checks.
See the [conductor-388 results](../reports/frobenius_rforest_even_20261101.md).
