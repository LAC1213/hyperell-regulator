# Vendored rforest

Source: https://github.com/edgarcosta/rforest

Pinned upstream commit: `3103d396c67cb1685131b1f11e84975cca335bdf`.
The maintained source is the rforest 1.0 C library distributed by Andrew
Sutherland, with its supporting David Harvey `zz` FFT implementation.

SHA-256 of the upstream commit archive:
`3e94db94c2cc518fa86633d136607ea2af1f4d327f0f5cbcc9920d474ce55754`.

`LICENSE` covers the rforest sources (MIT). `COPYING` contains the BSD
2-clause notice for the supporting `zz` sources. Both notices are also
shipped with binary distributions in `hyperell_regulator/RFOREST-COPYING.txt`.

Local changes:

- `hwmpz.c` and `hwmpz.h`: make FFT state and controls thread-local. Retain
  FFT plans through reference-counted setup/clear calls so a worker can
  reuse them across the numerator and denominator and successive pairs.
  Reject platforms without 64-bit `unsigned long` and GMP limbs at compile
  time, since the upstream arithmetic assumes these sizes.
- `hwmem.h` and `zzmem.h`: disable optional diagnostic allocation counters,
  which otherwise introduce shared writable state and counting overhead.
  Process-level peak memory is measured externally in the benchmarks.
- `rforest.c`: free the temporary initialized GMP vector in `mproduct()`.
- `rforest.c`: evaluate consecutive polynomial matrices using exact forward
  differences. `HYPELLFROB_RFOREST_EVAL=powers` selects the original evaluator
  for independent comparisons; both paths use the same forest arithmetic.
- `rforest.c` and `rforest.h`: configure worker-local controls through a C
  function, avoiding incompatible C/C++ thread-local wrapper symbols on macOS.

The FFT routines, crossover tables, and forest algorithm are upstream code.
The parent Makefile builds the C sources with `-O3`, C11, and GMP, without
an additional download or a Python wrapper. The C++ caller cancels common
polynomial factors before constructing a forest, computes the remaining
denominator valuations, and schedules reduction pairs across workers.
`HYPELLFROB_RFOREST_FFT=0` disables the FFT path for comparisons.

Algorithm references:

- Harvey and Sutherland, *Computing Hasse-Witt matrices of hyperelliptic
  curves in average polynomial time*, https://arxiv.org/abs/1402.3246.
- Harvey and Sutherland, *Computing Hasse-Witt matrices of hyperelliptic
  curves in average polynomial time, II*, https://arxiv.org/abs/1410.5222.
