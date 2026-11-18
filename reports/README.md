# Beilinson survey: conductor < 1300, 100 digits

- [Saved data and integration models](beilinson_100_digits.json).
- [LaTeX tables: torsion and integrality, models, numerical data](beilinson_100_digits.tex).

The survey contains 154 cached geometrically simple curves, with 22 successful
regulators and comparisons. Failed attempts and their reasons are retained.
The tables and explorer omit conductor 1116 since the integrality in that case fails.

The numerical target is 100 decimal digits. We use full cycle periods on the
minimum spanning tree basis and determinant weights `(1, 2 Re(z), |z|^2)`, with
`q = pi*c_plus*c_minus*L_second_derivative/D`. The relative error is
`abs(q/rational - 1)`; accuracy targets are not certified bounds.
An undecided component test is displayed as `Undecided`.

Tables display 10 significant digits. The JSON retains full decimal strings,
matrices, component verdicts and exact integration models. Points are given
in the displayed model; a dagger marks a reciprocal coordinate change.

Recompute the survey in SageMath:

```sh
sage -python -m hyperell_regulator.experiments.beilinson 1300 --fetch --force --digits 100 --reg-digits 100
sage -python -m hyperell_regulator.experiments.beilinson_report
```

Fetching requires network access. Cached L-values are reused when their
precision and coefficient counts suffice.

Rebuild the tables from the saved JSON:

```sh
sage -python -m hyperell_regulator.experiments.beilinson_report --input reports/beilinson_100_digits.json
```

Source: `results/beilinson_results.json`.
SHA-256: `0148c0633c63ad769b7432fb57095844240b22c4dc68ad5141059e5ba1295514`.
