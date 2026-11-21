# rforest

Source: <https://github.com/edgarcosta/rforest>,
commit `3103d396c67cb1685131b1f11e84975cca335bdf`.
Archive SHA-256: `3e94db94c2cc518fa86633d136607ea2af1f4d327f0f5cbcc9920d474ce55754`.

Andrew Sutherland's rforest library includes David Harvey's `zz` FFT code.
`LICENSE` covers rforest (MIT); `COPYING` covers `zz` (BSD 2-clause).
Both notices are included in `hyperell_regulator/RFOREST-COPYING.txt`.

Local changes:

- `hwmpz.c`, `hwmpz.h`: thread-local FFT state and reference-counted plans;
  require 64-bit `unsigned long` and GMP limbs. Matrix multiplication sizes
  transforms from the sum of operand bit lengths and streams left-hand
  matrix rows while reusing right-hand transforms to limit peak memory.
- `hwmem.h`, `zzmem.h`: disable diagnostic allocation counters.
- `rforest.c`: free the temporary GMP vector in `mproduct()` and evaluate
  polynomial matrices by exact forward differences. Long gaps between
  requested endpoints use balanced products reduced modulo the remaining
  moduli, including the initial gap in a precision retry with cached factors.
- `rforest.c`, `rforest.h`: configure thread-local controls through a C function.

The parent Makefile builds with C11, GMP and `-O3`. The C++ caller cancels
common polynomial factors, computes denominator valuations and schedules
reduction pairs across workers. For comparisons,
`HYPELLFROB_RFOREST_EVAL=powers` selects evaluation by powers and
`HYPELLFROB_RFOREST_FFT=0` disables FFT arithmetic.

References:

- Harvey and Sutherland, [Computing Hasse-Witt matrices of hyperelliptic
  curves in average polynomial time](https://arxiv.org/abs/1402.3246).
- Harvey and Sutherland, [Computing Hasse-Witt matrices of hyperelliptic
  curves in average polynomial time, II](https://arxiv.org/abs/1410.5222).
