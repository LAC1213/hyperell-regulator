"""Reproducible native matrix comparison and timing (run with sage -python).

Example:
  sage -python hypellfrob-threaded/benchmark.py --binary hypellfrob-threaded/build/euler \
      --bounds 1000 5000 20000 --workers 1 4 --repeats 3 --output /tmp/frobenius.json

Times include process startup, every requested prime, and native fallback for
the remainder backend's exceptional primes. They exclude Sage lifting and
L-function evaluation. Correctness compares every matrix entry modulo p^N.
"""
import argparse
import json
import os
import platform
from pathlib import Path
import statistics
import subprocess
import time

from sage.all import PolynomialRing, ZZ, prime_range


def precision(p, genus):
    # Exact integer form of p^N >= 2 binom(2g,g) p^(g/2).
    import math
    rhs = 4 * math.comb(2*genus, genus)**2 * p**genus
    n = 1
    while p**(2*n) < rhs:
        n += 1
    return n


def run(binary, workers, backend, batch, timeout):
    start = time.perf_counter()
    result = subprocess.run([binary, str(workers), backend], input=batch,
                            text=True, capture_output=True, check=True,
                            timeout=timeout)
    elapsed = time.perf_counter()-start
    if result.stderr:
        raise RuntimeError(result.stderr)
    return elapsed, [list(map(int, row.split())) for row in result.stdout.splitlines()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True)
    parser.add_argument("--coefficients", nargs="+", type=int, default=[1, 0, 0, 1, 0, 1])
    parser.add_argument("--bounds", nargs="+", type=int, default=[1000, 5000, 20000])
    parser.add_argument("--workers", nargs="+", type=int, default=[1, 4])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=600)
    parser.add_argument("--output")
    args = parser.parse_args()
    if args.repeats < 1 or any(w < 1 for w in args.workers):
        parser.error("repeats and worker counts must be positive")
    binary = str(Path(args.binary).resolve())
    q = PolynomialRing(ZZ, "x")(args.coefficients)
    if q.degree() < 3 or q.degree()%2 == 0 or not q.is_monic() or not q.is_squarefree():
        parser.error("the benchmark needs a squarefree monic odd degree polynomial")
    genus = (q.degree()-1)//2
    discriminant = q.discriminant()
    report = {
        "binary": binary, "platform": platform.platform(),
        "cpu_count": os.cpu_count(), "coefficients": args.coefficients,
        "repeats": args.repeats,
        "comparison": "full matrices modulo p^N; identical input jobs and worker counts",
        "measurements": [],
    }
    for bound in args.bounds:
        jobs = [(int(p), precision(int(p), genus)) for p in prime_range(3, bound)
                if discriminant % p]
        batch = " ".join(map(str, args.coefficients))+"\n"
        batch += "".join("%d %d\n" % job for job in jobs)
        for workers in args.workers:
            timings = {"legacy": [], "remainder": []}
            reference = None
            for repeat in range(args.repeats):
                order = ["legacy", "remainder"]
                if repeat%2:
                    order.reverse()
                for backend in order:
                    elapsed, rows = run(binary, workers, backend, batch, args.timeout)
                    assert len(rows) == len(jobs)
                    assert [row[0] for row in rows] == [p for p, n in jobs]
                    if reference is None:
                        reference = rows
                    assert rows == reference, (bound, workers, backend, "matrix mismatch")
                    timings[backend].append(elapsed)
            _, strict = run(binary, workers, "remainder-only", batch, args.timeout)
            assert len(strict) == len(reference)
            covered = 0
            for row, expected in zip(strict, reference):
                if row[1] == 1:
                    assert row == expected
                    covered += 1
            if jobs and covered == 0:
                raise RuntimeError("new backend handled no primes; increase the bound")
            medians = {k: statistics.median(v) for k, v in timings.items()}
            measurement = {
                "bound": bound, "jobs": len(jobs), "workers": workers,
                "remainder_primes": covered,
                "native_matrices": sum(row[1] == 1 for row in reference),
                "seconds": timings, "median_seconds": medians,
                "speedup": medians["legacy"]/medians["remainder"],
                "all_matrices_equal": True,
            }
            report["measurements"].append(measurement)
            print(json.dumps(measurement), flush=True)
            if args.output:
                Path(args.output).write_text(json.dumps(report, indent=2)+"\n")
    if not args.output:
        print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
