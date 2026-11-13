# Full cycle periods on the minimum spanning tree

Full cycle periods now use the Euclidean MST of the finite branch values,
with each row equal to twice the integral along one lifted edge. The
logarithmic Stokes integrator also uses MST cuts; the semilunar backend and
its selector have been removed.

For odd degree, all `2g` finite edges give an integral cycle basis. For even
degree, the `2g+1` doubled edges have one relation. Each leaf edge has
coefficient `+1` or `-1` in this relation, so omitting one gives an integral
basis. The implementation omits the leaf of largest modulus, breaking ties
by real and imaginary parts. Tree ties use a geometric ordering of vertices.
The topology is chosen in float64; integration endpoints retain guarded
precision.

The regression check reconstructs the change of basis against frozen
period matrices from the former implementation for five models (real and
complex branch values, odd and even degree, including conductor 249 and
277). The change matrices have integer entries, determinant `+1` or `-1`,
and reconstruction errors below `1e-27`. Before the engines were combined,
independent SciPy quadrature agreed with guarded quadrature on these paths.
Both output modes now use the guarded engine, and the saved matrices below
remain regression references for the common implementation.

## Normalization comparison at 100 digits

The package defines `c_plus` and `c_minus` as signed pivoted maximal minors
of the real and imaginary period matrices. Although the cycle bases are
integrally equivalent, these selected minors can change by rational factors
and signs. They are not a basis-independent normalization. Their definition
has been retained, and the table records the change explicitly.

All 22 successful saved survey curves and the additional nonintegral
conductor-249 example were checked at a 100-digit target. Each period-minor
ratio agrees with the rational factor shown to a scaled tolerance of
`1e-97`. For previously resolved rational comparisons, the new quotient
agrees with the old quotient multiplied by both factors to the same
relative tolerance. The saved regulator determinants and L-values are used
unchanged. These numerical checks are not interval-certified bounds.

| Curve | New/old c+ | New/old c− | New π c+ c− L″ / D |
|---|---:|---:|---:|
| 249.a.249.1 | 1 | 1/2 | 14 |
| 277.a.277.1 | -1 | 2 | 30 |
| 295.a.295.1 | 1 | 1/2 | -14 |
| 349.a.349.1 | -1 | -1/4 | 13/2 |
| 388.a.776.1 | -1 | 1 | 7 |
| 389.a.389.1 | 1 | -1 | 80 |
| 389.a.389.2 | -1 | 1 | -20 |
| 523.a.523.1 | -1 | 1 | -20 |
| 603.a.603.1 | 1 | 1 | 40960 |
| 603.a.603.2 | 1 | 2 | 10 |
| 691.a.691.1 | 1 | -1 | -16 |
| 709.a.709.1 | 1 | 2 | -16 |
| 763.a.763.1 | 1 | -1 | 40 |
| 909.a.909.1 | 1 | -1 | -8 |
| 925.a.925.1 | -1 | -1 | -4 |
| 997.a.997.1 | 1 | -1 | 64 |
| 997.a.997.2 | 1 | -1/2 | -16 |
| 1051.b.1051.1 | 1 | 1 | 16 |
| 1051.b.1051.2 | 1 | -1 | -16 |
| 1116.a.214272.1 | 1 | 2 | 27.6222341473748850… |
| 1123.a.1123.1 | -1 | -1 | 48 |
| 1147.a.35557.2 | -1 | 1 | 64 |
| 249.a.6723.1 | -1 | -1 | 24.8297858493021684… |

For `249.a.249.1`, `c_plus` is unchanged and `c_minus` is halved, so the
quotient is now **14**, formerly 28. The regulator integral itself is
unchanged. Archived survey tables and explorer captions retain their
original period normalization; new computations use the MST basis and the
experiment cache version has been advanced.

[Full period matrices, minors, quotients and source hashes](mst_periods.json)
are saved separately from the historical survey. Reproduce with:

```sh
sage -python -m hyperell_regulator.experiments.mst_periods
```

The additional nonintegral example is included when
`results/249_rational_approximations/input.json` is available.
