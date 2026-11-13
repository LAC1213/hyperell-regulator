# Minimum spanning tree cuts versus semilunar paths

This is an archived comparison made before the semilunar backend was removed.
Only MST integration is now available. Full cycle periods also now use MST
edges; see [MST period normalizations](mst_periods.md).

The default regulator cut is now the Euclidean minimum spanning tree of
all finite branch values and selected logarithmic centers. A rightward ray
from the rightmost vertex connects the tree to infinity. This ray cannot
cross the finite tree, since every finite point has smaller or equal real
part. The former selectable semilunar backend has since been removed.

The segment substitutions, quadrature, two-sheet Stokes sum and second-moment
formula are shared by both methods. At a junction the tree boundary follows
the cyclic order of incident edges; each bank is traversed exactly once.
For an edge oriented towards infinity, the logarithm jump is `-2 pi i` if its
child subtree contains the logarithmic center and zero otherwise. The sheet
gluing switches signs precisely when that subtree has an odd number of branch
values. Tests check these two independent descriptions against the contour
labels. During this archived comparison, full cycle periods and their signed
minors retained the semilunar cycle basis, preserving its normalization.
Current full periods use MST edges, with normalization changes recorded in
[MST period normalizations](mst_periods.md).

## Timing and accuracy

The comparison uses all 22 successful saved survey curves, plus the
nonintegral 249.a.6723.1 example. Each method runs three times in the same
process, with method order alternating between repetitions. The table gives
median wall times. Timings include cut construction, topology, quadrature
and assembly; they exclude parsing the reference model and recomputing L-values.
Runs used a single process with BLAS/OMP threads limited to one.

At a 30-digit target:

| Curve | Semilunar seconds | MST seconds | Speedup | Edges including the ray |
|---|---:|---:|---:|---:|
| 249.a.249.1 | 1.417 | 1.386 | 1.02x | 7 → 7 |
| 249.a.6723.1 | 1.215 | 0.447 | 2.72x | 7 → 7 |
| 389.a.389.1 | 4.704 | 0.721 | 6.53x | 30 → 6 |
| 603.a.603.1 | 50.574 | 0.886 | 57.07x | 6 → 6 |
| 763.a.763.1 | 4.470 | 0.556 | 8.04x | 6 → 6 |
| 997.a.997.1 | 2.116 | 0.456 | 4.64x | 29 → 5 |
| 1051.b.1051.2 | 7.333 | 0.525 | 13.97x | 5 → 5 |
| 1147.a.35557.2 | 2.107 | 0.430 | 4.90x | 29 → 5 |

The median per-curve speedup is **1.68x**.
Summing the per-curve medians gives **91.49 s → 16.69 s** (5.48x).
The aggregate gain is strongly influenced by conductor 603; it is not a
prediction of the speedup for a typical curve. Finite edges total
193 → 121.

Conductor 603 illustrates a gain beyond reducing the number of edges:
it still has six edges including the ray, but quadrature evaluations fall
from 139,987 to
5,253. Both finite edge geometry and ray
attachment change. Length minimization alone does not guarantee a minimum
quadrature cost.

At a 100-digit target, also with three repetitions:

| Curve | Semilunar seconds | MST seconds | Speedup |
|---|---:|---:|---:|
| 249.a.249.1 | 4.406 | 3.201 | 1.38x |
| 997.a.997.1 | 8.079 | 1.838 | 4.40x |
| 1147.a.35557.2 | 9.248 | 2.699 | 3.43x |
| 249.a.6723.1 | 9.415 | 3.575 | 2.63x |

All paired matrices agree exactly after rounding to the requested output
precision, at both 30 and 100 digits. Each run is additionally checked against
the saved 100-digit matrix and determinant, with a scaled tolerance of
`10^(-digits+3)`. These numerical checks are not interval-certified bounds.

[Full timings, matrices, evaluation counts and source hashes](mst_cut_comparison.json)
retain every repetition. [Geometry comparison](mst_cut_comparison.svg) shows
both constructions with identical coordinate scales for each curve.

## Archived evidence and current verification

The saved JSON contains all measured repetitions, matrices, evaluation counts
and hashes of the implementation used at the time. The SVG preserves the
geometry comparison. The two-backend benchmark and rendering drivers were
retired with the semilunar backend.

Recheck the current MST integral values without overwriting this archive:

```sh
sage -python -m hyperell_regulator.experiments.x_coordinate_comparison --output /tmp/249-mst-verification.json
sage -python -m hyperell_regulator.experiments.mst_periods
```

The conductor-249 integral check includes the nonintegral example when
`results/249_rational_approximations/input.json` is available.
Explorer exports retain explicit tree endpoints and the cyclic boundary word.
