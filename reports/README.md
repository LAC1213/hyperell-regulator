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

## Manuscript tables with full cycle periods

`paper_full_periods.tex` contains replacements for the manuscript's torsion,
model and numerical tables using the current `c_invariants` convention.
The script produces this single LaTeX file, with all three table environments.
All sixteen updated paper factors
`q = 2*pi*c_plus*c_minus*L_second_derivative/D` are integers.

Regenerate from the repository root:

```sh
sage -python -m reports.remake_paper_tables
# Choose another output file:
sage -python -m reports.remake_paper_tables --out /tmp/paper_tables.tex
```

The script reads the portable `paper_full_periods_input.json`, so it does not
need the attachment, LMFDB caches or local experiment results. It computes
periods at 50 and 110 digits and reuses the saved determinants and L-values.
The table captions identify the lower precision available for two rows and
the unverified L-value/integrality claim for `388.a.776.1`; unavailable errors
are left blank. The LaTeX uses `amsmath`, `booktabs` and `hyperref`.
