# Beilinson survey: conductor < 1300, 100 digits

This directory contains the current 100-digit Beilinson survey:

- [Full saved data and integration models](beilinson_100_digits.json).
- [Three LaTeX tables: torsion/integrality, models, and numerical data](beilinson_100_digits.tex).

The snapshot has 154 curve records, 22 successful
regulators and 22 numerical comparisons. Unsuccessful attempts and their
reasons are retained in the JSON. The three tables display completed
comparisons. This is a survey of locally cached geometrically simple
curves, without a claim of completeness over all curves.
The conductor-1116 record is retained in the JSON, but omitted from the LaTeX
tables and explorer.
The LaTeX tables and explorer include 21 successful regulators.

Regulators and L-values target 100 decimal digits. The determinant uses weights
`(1, 2 Re(z), |z|^2)`, full cycle periods, and
`q = pi*c_plus*c_minus*L_second_derivative/D`. The comparison error is
`abs(q/rational - 1)`. An undecided component verdict is displayed as
`Undecided` in the integrality column. Precision targets are not certified
error bounds.

Period normalizations follow the cycle basis used in each record. The saved
pre-MST snapshot uses its original cycle basis; current calculations use MST
cycles. Selected signed minors and quotients can therefore differ by rational
factors. See [the MST period comparison](mst_periods.md) for all successful
saved models, including the change from 28 to 14 for conductor 249.

The numerical table displays 10 significant digits. The JSON preserves the
complete saved decimal strings, matrices, precision settings, component
verdicts and integration models. Table points use the displayed integration
models. Reciprocal coordinate changes are marked by a dagger in the model
table, with their shift and scale given explicitly. The JSON retains both
the original points and their integration-model images.

After the calculation has finished, regenerate without repeating integrations:

```sh
sage -python -m hyperell_regulator.experiments.beilinson_report
```

Source: `results/beilinson_results.json` (read without modification).
Source SHA-256: `c04e1821eb7d9b279b7771e82fbf6f91d218efef47e479844fb250c918598d41`.
