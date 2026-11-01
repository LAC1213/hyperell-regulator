# Beilinson survey: conductor < 1300, 100 digits

This directory contains only the current 100-digit Beilinson survey:

- [Full saved data and exact functions](beilinson_100_digits.json).
- [LaTeX numerical and algebraic tables](beilinson_100_digits.tex).

The snapshot has 154 curve records, 22 successful
regulators and 22 numerical comparisons. Unsuccessful attempts and their
reasons are retained in the JSON; the tables list their status counts and any
numerical failures. This is a survey of locally cached geometrically simple
curves, without a claim of completeness over all curves.
The conductor-1116 record is retained in the JSON, but omitted from the LaTeX
tables and explorer. Table status counts use the same exclusion.
The LaTeX tables and explorer include 21 successful regulators.

Regulators and L-values target 100 decimal digits. The determinant uses weights
`(1, 2 Re(z), |z|^2)`, full cycle periods, and
`q = pi*c_plus*c_minus*L_second_derivative/D`. The comparison error is
`abs(q/rational - 1)`. An undecided component verdict remains undecided and is
marked by a dagger. Precision targets are not certified error bounds.

The numerical table displays 10 significant digits. The JSON preserves the
complete saved decimal strings, matrices, precision settings, component
verdicts and exact functions. Points are in the original cached model;
functions are on the displayed integration model, with its coordinate change.

After the calculation has finished, regenerate without repeating integrations:

```sh
sage -python -m hyperell_regulator.experiments.beilinson_report
```

The standalone LaTeX document uses `geometry`, `amsmath`, `amssymb`, `booktabs`
and `longtable`. No PDF is generated.

Source: `results/beilinson_results.json` (read without modification).
Source SHA-256: `c04e1821eb7d9b279b7771e82fbf6f91d218efef47e479844fb250c918598d41`.
