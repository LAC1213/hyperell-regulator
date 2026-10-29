# Paper factors with full semilunar cycle periods

Recomputed on 29 October 2026, using the manuscript models and the root of largest absolute value as the starting point. Period minors were evaluated at 50 and 90 decimal digits.

The manuscript convention is `q = 2*pi*c_plus*c_minus*L_second_derivative/D`, the reciprocal of the experiment ratio. Thus `q_new = q_old * (c_plus_new*c_minus_new)/(c_plus_old*c_minus_old)`.

| Curve | Paper q | Period-product multiplier | Updated q |
|---|---:|---:|---:|
| `249.a.249.1` | -7 | 8 | -56 |
| `277.a.277.1` | -15/2 | -8 | 60 |
| `295.a.295.1` | 7 | -16 | -112 |
| `349.a.349.1` | -13/2 | -8 | 52 |
| `388.a.776.1` | -7/2 | 8 | -28 |
| `389.a.389.1` | -20 | 16 | -320 |
| `389.a.389.2` | -10 | -8 | 80 |
| `523.a.523.1` | 10 | 8 | 80 |
| `709.a.709.1` | -8 | 4 | -32 |
| `763.a.763.1` | -10 | 16 | -160 |
| `925.a.925.1` | -2 | 8 | -16 |
| `997.a.997.2` | -16 | -8 | 128 |
| `1051.b.1051.1` | 8 | 8 | 64 |
| `1051.b.1051.2` | 8 | 8 | 64 |
| `1123.a.1123.1` | 24 | 8 | 192 |
| `1147.a.35557.2` | -16 | 16 | -256 |

All sixteen updated factors are integers.
The inferred period-product multipliers agree with the 90-digit computations to absolute error below `1e-85`; the individual minors agree between the two precisions to relative error below `1e-47`.

This calculation reuses the manuscript rational factors and the 110-digit period minors from the earlier audit. It does not freshly recompute regulators or L-values. Saved numerical identity checks are included where available; `388.a.776.1` rescales the manuscript claim because its L-value was not independently reproduced by that audit. These are numerical checks, without interval certificates.

The check also exposed quadrature-sized differences between mathematically equal pivots for `277.a.277.1`. Period-minor pivot ties now use quadrature accuracy, giving consistent signs under precision refinement.

Full values and provenance: [JSON evidence](paper_period_normalisation_20261029.json). The [input snapshot](paper_period_normalisation_input.json) contains the manuscript models, original factors and earlier period minors.

Reproduce from the repository root with `sage -python -m reports.check_paper_period_normalisation`.
