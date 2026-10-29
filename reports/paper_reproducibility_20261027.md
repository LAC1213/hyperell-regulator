# Reproducibility of the manuscript's final three tables

Checked on 27 October 2026 against checkout `c036a1b68b839b7919cf9c7b26949d7a1fc34b77`.
The input is the supplied manuscript, *On Beilinson's conjecture for the square of a genus 2 curve*.
Its SHA-256, package source hashes, input cache hashes, environment, exact checks and full numerical results are in [the JSON evidence](paper_reproducibility_20261027.json).
The package implementation was left unchanged.

The complete sixteen-row, approximately 100-digit claim was **not fully reproduced** by the current experiment. Nine rows passed at 100 digits. Five more passed at 30 digits, and `277.a.277.1` passed at 15 digits. `388.a.776.1` is skipped by the experiment because the component test at 2 is undecided; its standard L-value also requires an externally supplied bad Euler factor. A separate numerical run recovered its displayed determinant at a 15-digit target.

All sixteen displayed determinants match after accounting for coordinate and column changes. All sixteen models, thirty-two points, thirty-two torsion orders relative to the stated base point, and thirty-two displayed period values match. The fifteen supported L-values match their displayed values. The remaining L-value, for `388.a.776.1`, was not reproduced.

The first pass used `--force --prec 100 --digits 100 --edge-seconds 180`, with one selected curve per subprocess and eight concurrent subprocesses. `523.a.523.1` was retried with a 360-second edge limit and succeeded in about 451 seconds overall. Wall time limits depend on the machine and concurrent load; a timeout here leaves higher-precision reproduction unresolved.

The installed environment was SageMath 10.7, Python 3.11.14, NumPy 2.4.3, SciPy 1.17.1 and PARI 2.17.3. The manuscript's historical notes mention SageMath 10.9.

The experiment reused the existing LMFDB curve and component inputs and its saved L-value cache, copied into private data directories. Every regulator was freshly integrated. Independent fresh component checks were then run for the actual manuscript point/base pairs. Periods on the manuscript models were separately recomputed with a 110-digit target. For L-values, thirteen existing 100-digit entries passed the current coefficient-count requirements and their saved functional-equation checks; `277.a.277.1` and `389.a.389.1` were independently reevaluated at 30 digits. This is **not an independent recomputation of every 100-digit L-value**.

| Curve | Paper q | 100-digit experiment | Best completed identity check |
|---|---:|---|---|
| `249.a.249.1` | `-7` | Timeout: edge 5 (degree 14) | 30 digits; 6.29E-31 |
| `277.a.277.1` | `-15/2` | Timeout: edge 0, degree 15 | 15 digits; 1.94E-15 |
| `295.a.295.1` | `7` | Timeout: edge 0 (degree 14) | 30 digits; 9.40E-31 |
| `349.a.349.1` | `-13/2` | Timeout: edge 4 (degree 13) | 30 digits; 3.86E-31 |
| `388.a.776.1` | `-7/2` | Skipped: undecided at 2 | D only: printed digits match at 15-digit target |
| `389.a.389.1` | `-20` | Sheet closure error | 30 digits; 1.66E-30 |
| `389.a.389.2` | `-10` | Timeout: edge 2 (degree 5) | 30 digits; 4.06E-31 |
| `523.a.523.1` | `10` | Pass (360 s retry) | 100 digits; 1.85E-101 |
| `709.a.709.1` | `-8` | Pass | 100 digits; 1.19E-101 |
| `763.a.763.1` | `-10` | Pass | 100 digits; 3.39E-102 |
| `925.a.925.1` | `-2` | Pass | 100 digits; 3.27E-101 |
| `997.a.997.2` | `-16` | Pass | 100 digits; 8.98E-101 |
| `1051.b.1051.1` | `8` | Pass | 100 digits; 3.34E-104 |
| `1051.b.1051.2` | `8` | Pass | 100 digits; 8.96E-101 |
| `1123.a.1123.1` | `24` | Pass | 100 digits; 3.03E-101 |
| `1147.a.35557.2` | `-16` | Pass | 100 digits; 9.50E-101 |

The errors above were calculated from fresh determinants and the independently recomputed periods on the manuscript models, using the manuscript's rational factor. The successful 100-digit discrepancies are all below `1e-100`. The historical error mantissas need not repeat exactly. These are numerical comparisons, without rigorous interval enclosures.

The package reports `r = D/(2*pi*c_plus*c_minus*L_second_derivative)`; the manuscript uses its reciprocal `q`. The experiment can select a different base point, coordinate model or logarithmic column order. The checker transforms the determinant before comparing it with the manuscript. In particular, the reciprocal moment transformation has determinant -1; changes of the logarithmic basis contribute a further sign. Period minors were recomputed directly on the manuscript models rather than inferred from those signs.

For `389.a.389.1`, the 100-digit run ended with sheet-boundary discrepancy approximately `2.6078e-30`, against tolerance `1e-97`; its 30-digit run succeeded. For `277.a.277.1`, the 30-digit run also failed sheet closure, with discrepancy approximately `1.8347e-17`, against tolerance `1e-27`; its 15-digit run recovered the printed determinant and the factor `-15/2` with relative discrepancy approximately `1.94e-15`.

The fresh component tests returned True for both manuscript point/base pairs on each of fifteen curves. For `388.a.776.1`, both pairs returned `{2: undecided, 97: True}`. Its separate regulator-only runs do not settle integrality. The manuscript's manual argument at 2 and supplied Euler factor are additional inputs that the unmodified experiment does not automatically use.

Three caption corrections are needed in the torsion/integrality table:

- The third column's orders are relative to `O^+`, as specified in the surrounding model description. The caption says `O^-`. For example, `(0,1)` on `277.a.277.1` has order 3 relative to `O^+`, and order 5 relative to `O^-`.
- The fourth column records integrality, although the caption describes another order column.
- The error formula should use `2*pi`, consistently with the real determinant and real period columns and with the displayed identity in the body. Its extra `i` would make the reported tiny errors impossible under those real conventions.

The experiment checks geometrical simplicity through the cached `is_simple_geom` flag. It does not independently establish the stronger endomorphism-ring assertion `End(Jac(C)_Qbar) = Z` in the table caption.

Reproduce the primary run from the repository root with [the checker](reproduce_paper_tables.py):

```sh
export DOT_SAGE=/tmp/hyperell-paper-audit-sage
sage -python -m reports.reproduce_paper_tables \
  --paper '/path/to/the/supplied/manuscript.tex' \
  --prec 100 --digits 100 --edge-seconds 180 --jobs 8
sage -python -m reports.reproduce_paper_tables --components-only
sage -python -m reports.reproduce_paper_tables --summarize
```

In this workspace, Sage is `/home/lambert/anaconda3/envs/sage/bin/sage`. The default output directory is `results/paper_reproducibility_20261027`. To run a separate lower-precision survey, use `--prec 30 --digits 30 --out-dir results/paper_reproduction_30`. Use a distinct output directory to retain the 100-digit evidence.

The longer retry was:

```sh
HYPERELL_REGULATOR_DATA="$PWD/results/paper_reproducibility_20261027/523.a.523.1" \
sage -python -m hyperell_regulator.experiments.beilinson 1148 \
  --force --prec 100 --digits 100 --edge-seconds 360 --verbose \
  --out results/paper_reproducibility_20261027/523.a.523.1/retry_100.json
```

Full logs and intermediate results remain under `results/paper_reproducibility_20261027/`; the portable JSON evidence includes the completed runs and additional checks.
