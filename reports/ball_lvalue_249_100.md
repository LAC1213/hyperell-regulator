# Conditional 100-digit standard L-value: 249.a.249.1

For the curve `y^2 + (x^3+1)y = x^2+x`, the computation encloses
`L''(M_std,1)`, with standard conductor `N=62001`. Rounded to **100 decimal
places**, its value is

```text
5.0676270403470587136420759409704921633659419584973778429481177000287672193848807749605780500214686347
```

The full stored interval has radius less than **5.228 × 10^-120** and 398 bits
of achieved relative accuracy. Both interval endpoints round to the displayed
decimal. The radius refers to the full interval, not to the rounded display.

This is a rigorous propagated numerical enclosure **conditional on the stated
analytic inputs**: the correct motivic Euler and gamma factors, and the
holomorphy, finite-order, functional-equation and coefficient-growth hypotheses
of the analytic algorithm. This computation does not prove those inputs or an
exact Beilinson relation.

- [Exact dyadic intervals, assumptions, checksums and provenance](ball_lvalue_249_100.json).
- [Complete exact Euler input, gzip compressed](ball_lvalue_249_100_euler.txt.gz).
- [Backend output](ball_lvalue_249_100_run.txt).
- [Implementation and investigation notes](../docs/ball-arithmetic.rst).

The computation evaluates `L(std,2)` and uses
`L''(std,1) = sqrt(N)/(4*pi^2) * L(std,2)`. The JSON also retains the real and
imaginary enclosures for `L(std,2)` and the normalization factor. It records
42861 prime Euler factors through coefficient index 517454, a 600-bit internal
target, 704-bit working arithmetic, and an explicit 333-bit acceptance floor.
Acceptance uses the final achieved bound, including coefficient tails and
Fourier/interpolation errors. Rank and zero counting are not requested.

The optional backend is Costa–Platt's `lfunctions`, revision
`42014e1204f91e93f72b46945c4f5701b3ed7b1b`, with both repository patches:

- [Values-only entry point](../hyperell_regulator/experiments/lfunctions_values_only.patch).
- [Fourier error-bound corrections](../hyperell_regulator/experiments/lfunctions_error_bounds.patch).

The latter fixes the destination index of two error terms and the
frequency-dependent coefficient cutoff used in the truncation bound. It also
certifies that the floating-point bucket choices are within the Taylor bound's
half-grid radius. The bridge requires a revision symbol from this correction;
it cannot silently link an uncorrected backend. The earlier
`ball_lvalue_249_80.json` record is withdrawn because it used the uncorrected
error propagation.

The exact factor construction uses the rationally isomorphic odd model
`Y^2=X^5+13X^4+72X^3+240X^2+384X+256`. The coordinate identity was checked
exactly. Native matrices were reconstructed by the trace identities for the
first two characteristic-polynomial coefficients and the proven Weil bound.
The new factors matched all 22044 cached factors in the overlap. Direct point
counts checked the factor at 2 and the nine small primes declined by the native
backend. Independent Sage computations checked four large primes.

After applying both patches and building the bridge as described in the
implementation notes, reproduce the value from the saved input, from the
repository root:

```sh
gzip -dc reports/ball_lvalue_249_100_euler.txt.gz > /tmp/ball-249-euler.txt
/tmp/ball-lvalue-probe 600 62001 /tmp/ball-249-euler.txt 333 517454
sage -python -m pytest -q tests/test_saved_ball_lvalue_100.py
```

The artifact tests check all prime positions and input hashes, the final radius
and decimal rounding with exact rational arithmetic, and the normalization
against independent rational bounds for pi from Machin's identity. Reproducing
the analytic enclosure additionally requires running the patched backend.
