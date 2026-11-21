# Beilinson survey: conductor < 1300, 100 digits

- [Saved data and integration models](beilinson_100_digits.json).
- [LaTeX tables: torsion and integrality, models, numerical data](beilinson_100_digits.tex).

Arithmetic checks for all 154 entries are available as
[CSV](beilinson_arithmetic_100.csv) and [JSON](beilinson_arithmetic_100.json):
good ordinary reduction at 2 and 3, and surjectivity on Jacobian 3-torsion.
The JSON retains exact point counts, Frobenius witnesses and an LMFDB metadata
snapshot. All 106 positive surjectivity answers have Frobenius certificates;
43 negative answers follow from rational 3-torsion or rational endomorphisms.
The other five negative answers are LMFDB classifications, whose mod-3 image
computations use probabilistic methods. Recheck the data without network access:

```sh
sage -python reports/verify_beilinson_arithmetic_100.py
```

The survey contains 154 cached geometrically simple curves, with 22 successful
regulators and comparisons. Failed attempts and their reasons are retained.
The tables and explorer omit conductor 1116 since the integrality in that case fails.

The numerical target is 100 decimal digits. We use full cycle periods on the
minimum spanning tree basis and determinant weights `(1, 2 Re(z), |z|^2)`, with
`q = pi*c_plus*c_minus*L_second_derivative/D`. The relative error is
`abs(q/rational - 1)`; accuracy targets are not certified bounds.
An undecided component test is displayed as `Undecided`.

A separate [certified 100-digit regulator/period example](certified_regulator_249_100.json)
stores exact interval endpoints and a comparison bound for 249.a.249.1 in the
model stated in that file. Its L-value is conditional on the recorded analytic
hypotheses; the interval comparison does not prove exact rational equality.
See the [ball arithmetic guide](../docs/ball-arithmetic.rst) to reproduce such
computations with `certified=True` or `--certified`.

The [certified tree result](certified_tree_249_100.json) records the same
100-digit matrix using finite tree edges and one infinity ray. It includes
the boundary word, sheet permutations, integer logarithm turns and work
counts. Every matrix entry overlaps the independent parallel-ray enclosure.

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
