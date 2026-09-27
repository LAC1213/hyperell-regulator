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

Regenerate without any LMFDB cache or numerical integration:

```sh
sage -python -m hyperell_regulator.experiments.latex_report reports/beilinson_100_digits.json --out reports/beilinson_100_digits.tex
```

The TeX uses `amsmath`, `amssymb`, `booktabs`, `longtable`, and `geometry`.
It is a standalone landscape document. Compile with `pdflatex` twice where
those LaTeX packages are installed; this workspace lacks a LaTeX format.

Source: `results/beilinson_results.json` (kept locally, excluded from Git).
Source SHA-256: `9179b3d7e7155adf507acae7b58da77a8c932454c6ef72e9086437e07340db31`.
The source was read without modification; numerical failures at 100 digits
were not rerun or suppressed.
