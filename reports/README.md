# Saved Beilinson survey

`beilinson_100_digits.tex` contains all 16 entries marked `ok` in the saved
100-digit survey. The numerical table and the accompanying algebraic table
are linked by curve label. Decimal values are abbreviated for readability;
`beilinson_100_digits.json` retains the full saved values and exact functions.

The error column is the absolute difference `abs(q - rational)`, computed
from the full decimal ratio. Three successful regulators have no L-value:
the two conductor-603 curves have an unsupported Frobenius action, and
1717.a.1717.2 exhausted the PARI stack. Their comparison cells are blank.
The omitted rows were not marked `ok`.

Points use the original cached curve coordinates. Functions use the displayed
integration model `v^2=G(u)`; the table gives its coordinate shift and scale.
The scalar normalization of each exact function is the package's normalization.
The normalized regulator includes the corresponding norm correction.

The saved TeX and JSON are retained as an archived snapshot.

The TeX uses `amsmath`, `amssymb`, `booktabs`, `longtable`, and `geometry`.
It is a standalone landscape document. Compile with `pdflatex` twice where
those LaTeX packages are installed; this workspace lacks a LaTeX format.

Source: `results/beilinson_results.json` (kept locally, excluded from Git).
Source SHA-256: `9179b3d7e7155adf507acae7b58da77a8c932454c6ef72e9086437e07340db31`.
The source was read without modification; numerical failures at 100 digits
were not rerun or suppressed.

## Integer quotient check

With weights `(1, 2 Re(z), |z|^2)` and `q = pi*c_plus*c_minus*L''/D`,
all 20 completed comparisons in the cached conductor-below-2000 survey and
all 7 family comparisons for `n=-4,...,2` give integers. The numerical
results and scope are saved in [beilinson_integer_factors.json](beilinson_integer_factors.json).
That saved survey excluded undecided component verdicts; the runner now
includes them by default.
`997.a.997.1` remains unresolved after a time-limited 15-digit retry and a
float64 sheet-closure failure; no quotient is claimed for that curve.

## Manuscript tables with full cycle periods

`paper_full_periods.tex` contains replacements for the manuscript's torsion,
model and numerical tables using the current `c_invariants` convention.
The script produces this single LaTeX file, with all three table environments.
The determinant uses weights `(1, 2 Re(z), |z|^2)`. All sixteen transported
paper factors `q = pi*c_plus*c_minus*L_second_derivative/D` are integers;
all sixteen now have numerical comparisons. The component condition at 2
for `388.a.776.1` remains undecided. Its freshly computed L-value and
functional-equation checks are saved in [conductor_388_lvalue.json](conductor_388_lvalue.json).

Regenerate from the repository root:

```sh
sage -python -m reports.remake_paper_tables
# Choose another output file:
sage -python -m reports.remake_paper_tables --out /tmp/paper_tables.tex
```

The script reads the portable `paper_full_periods_input.json`, so it does not
need the attachment, LMFDB caches or local experiment results. It computes
periods at 50 and 110 digits, doubles the saved determinants for the new
middle-row weight, and reuses the saved L-values.
The table captions identify the lower precision available for three rows and
the remaining integrality claim for `388.a.776.1`. The archived survey
above retains its historical normalization. The LaTeX uses `amsmath`, `booktabs` and `hyperref`.
