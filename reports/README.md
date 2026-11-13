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

Full cycle periods use the MST cycle basis recorded in each successful row.
The current code computes D and the signed period minors at a 100-digit target.
L-values use the precision-validated coefficient cache. The quotient uses
these ingredients in the same integration model as its displayed points.

The numerical table displays 10 significant digits. The JSON preserves the
complete saved decimal strings, matrices, precision settings, component
verdicts and integration models. Table points use the displayed integration
models. Reciprocal coordinate changes are marked by a dagger in the model
table, with their shift and scale given explicitly. The JSON retains both
the original points and their integration-model images.

From the repository root in SageMath, recompute the survey (network access is
needed to fetch curve data; cached L-values are reused after validation):

```sh
sage -python -m hyperell_regulator.experiments.beilinson 1300 --fetch --force --digits 100 --reg-digits 100 --edge-seconds 600
sage -python -m hyperell_regulator.experiments.beilinson_report
```

The bound is strict: conductor < 1300. A fresh cache also recomputes L-values
and can take substantially longer. Component/reduction computations require
the package's Sage dependencies. The command retains failed attempts and
does not claim that the locally cached curves form a complete list.

To rebuild the tables from the versioned JSON without repeating numerical
calculations or requiring local curve caches:

```sh
sage -python -m hyperell_regulator.experiments.beilinson_report --input reports/beilinson_100_digits.json
```

To typeset the survey, install the LaTeX packages in its preamble, including
`arydshln`. Compile into a separate directory to keep this folder tidy.
Run the final command twice to resolve cross-references:

```sh
mkdir -p /tmp/beilinson-survey
pdflatex -interaction=nonstopmode -halt-on-error -output-directory=/tmp/beilinson-survey reports/beilinson_100_digits.tex
```

Source: `results/beilinson_results.json` (read without modification).
Source SHA-256: `0148c0633c63ad769b7432fb57095844240b22c4dc68ad5141059e5ba1295514`.
