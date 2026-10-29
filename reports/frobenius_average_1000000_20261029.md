# Frobenius comparison below 1,000,000, 29 October 2026

The corrected remainder-tree backend completed the million-bound test:
**all 78,485 supported Frobenius matrices matched the existing backend
exactly**, comprising 1,255,760 entries. Ten sampled Euler factors,
including the factor at 999983, also agreed with Sage independently.
The corrected implementation used 12.8% less total CPU time, but took
17.3% longer in wall time with four workers and used substantially more RAM.

This extends the [smaller-bound benchmark](frobenius_average_20261029.md)
for `Q = x^5 + x^3 + 1`. There are 78,495 odd good primes below 1,000,000,
excluding the bad primes 53 and 61. The package's reconstruction bound
requests one job at precision 4, three at precision 3, and 78,491 at precision 2.

| Backend | Workers | Wall time | Total CPU time | Peak resident memory |
| --- | ---: | ---: | ---: | ---: |
| Existing | 4 | 751.86 s (12m 32s) | 3001.82 s | 131.5 MiB |
| Corrected remainder tree | 4 | 881.97 s (14m 42s) | 2617.57 s | 13.0 GiB |

These are single runs on a shared 16-logical-CPU Linux machine, not medians
from dedicated hardware. Both binaries use g++ 15.2.0, `-O2 -std=c++17`,
private thread-safe NTL 11.6.0, and system GMP. The existing per-prime
algorithm was unchanged between builds. The [raw results](frobenius_average_1000000_20261029.json)
record executable and input hashes, all attempts, resource limits, memory
samples, CPU times, and the full comparison counts.

The remainder run used `remainder-only`, so fallback could not conceal a
failure. It declined ten exceptional primes. Five (17, 19, 23, 29, 31)
are handled by the existing native algorithm; both native paths decline
3, 5, 7, 11, and 13. Timing excludes Sage fallback, characteristic-polynomial
lifting, and L-function evaluation. The corrected run's additional five
declines save negligible work, but mean these are not timings of complete
end-to-end L-value computations.

The remainder run averaged about three cores of CPU use, while the baseline
averaged about four. Samples show the largest reduction trees continuing
after other workers finished. This poorer parallel utilization accompanies
the higher wall time despite lower total CPU use. No process swap was
reported in the corrected run's periodic samples. A 20 GiB virtual-memory
limit applied from startup was not exhausted.

The first implementation failed the four-worker test after 660.08 seconds,
with `out of memory` (return code 134, SIGABRT) after reaching 19.3 GiB
resident memory. A 20 GiB virtual-memory cap had been applied during that
run to protect the shared machine from exhausting RAM and swap. It returned
no matrices, so it supplied no completed arithmetic comparison. A subsequent
one-worker retry was stopped when the complexity audit found the defect
below. Its interrupted timing is not used for performance conclusions.

The audit found two implementation problems:

- `reduce_centered()` computed `root_modulus / 2` at every product node.
  There are O(N) nodes and the root modulus has O(N) bits at fixed curve
  and precision. Repeating that division/copy introduced quadratic bit-work.
  The corrected version computes the threshold once per tree.
- Right-child products were retained after their last use. The corrected
  version frees them once their parent has been built; remainder descent
  needs only left-sibling products. This reduces retained integer payload.

Both changes pass all 66 targeted arithmetic, L-function, and build tests,
including genera 1--3 and precisions 1--4. The full million-bound comparison
then passed with zero mismatches. The corrected binary is now the source
checkout's default executable; the initial executable is archived under
`hypellfrob-threaded/build/euler-before-threshold-fix`.

Harvey's [Theorem 1](https://arxiv.org/pdf/1210.8239) is an asymptotic bit-complexity
result. For fixed curve and genus it gives total cost N times polylogarithmic
factors, improving on the roughly N^(3/2), up to logarithmic factors, cost
of separate square-root-in-p computations. It specifies no practical
crossover at N=1,000,000. The initial implementation's quadratic overhead
also meant its timings could not fairly test the intended complexity.
Correct arithmetic alone did not establish an efficient implementation.
The corrected measurements still establish no wall-time speedup on this
machine at this bound; improving scheduling and memory use remains necessary
before selecting this backend automatically.

Independent Sage checks used primes 104773, 224797, 350429, 479957, 611993,
746807, 882439, 951277, 993197, and 999983; all ten Euler factors agreed.

The benchmark utility can repeat matrix comparisons at this bound:

```sh
sage -python hypellfrob-threaded/benchmark.py \
    --binary hypellfrob-threaded/build/euler \
    --bounds 1000000 --workers 4 --repeats 1 --timeout 7200 \
    --output /tmp/frobenius-million.json
```

That utility times the hybrid remainder mode and performs an additional
strict run to verify coverage, so its complete protocol takes longer than
a single remainder run. The raw results above time strict mode directly.
