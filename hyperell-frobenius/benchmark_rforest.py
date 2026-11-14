"""Time native rforest runs against independent Sage or saved references.

Run with sage -python. Uses wait4 for per-process CPU time and peak RSS;
wall time is measured with perf_counter. Unsupported jobs are declined.
Reference generation is excluded from native timings.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time

from sage.all import GF, HyperellipticCurve, PolynomialRing, ZZ, prime_range
from sage.all import matrix as sage_matrix
def precision(p, genus):
    # Exact integer form of the package's Kedlaya reconstruction bound.
    import math
    rhs = 4 * math.comb(2*genus, genus)**2 * p**genus
    n = 1
    while p**(2*n) < rhs:
        n += 1
    return n


def parse_rows(text, jobs, genus):
    rows = [list(map(int, line.split())) for line in text.splitlines()]
    if len(rows) != len(jobs):
        raise RuntimeError("wrong number of native output rows")
    for (p, n), row in zip(jobs, rows):
        if len(row) < 2 or row[0] != p or row[1] not in (0, 1):
            raise RuntimeError("invalid native output row")
        if row[1] and (row[2] != 2*genus or len(row) != 3+(2*genus)**2 or
                       any(v < 0 or v >= p**n for v in row[3:])):
            raise RuntimeError("invalid Frobenius matrix")
    return rows


def run(binary, workers, backend, kappa, batch, timeout, memory_mib, environment):
    command = [binary, str(workers), backend, str(kappa)]

    def memory_limit():
        if memory_mib:
            limit = memory_mib*1024**2
            resource.setrlimit(resource.RLIMIT_AS, (limit, limit))

    with tempfile.TemporaryDirectory(prefix="rforest-benchmark-") as directory:
        directory = Path(directory)
        (directory/"input").write_text(batch)
        with (directory/"input").open() as inp, (directory/"output").open("w") as out, \
                (directory/"errors").open("w") as err:
            start = time.perf_counter()
            process = subprocess.Popen(command, stdin=inp, stdout=out, stderr=err,
                                       env=environment, preexec_fn=memory_limit)
            next_update = start+30
            try:
                while True:
                    pid, status, usage = os.wait4(process.pid, os.WNOHANG)
                    now = time.perf_counter()
                    if pid:
                        process.returncode = os.waitstatus_to_exitcode(status)
                        break
                    if now-start > timeout:
                        raise TimeoutError("native benchmark exceeded timeout")
                    if now >= next_update:
                        print(json.dumps({"event": "progress", "backend": backend,
                                          "workers": workers, "kappa": kappa,
                                          "elapsed_seconds": now-start}), flush=True)
                        next_update = now+30
                    time.sleep(0.1)
            except BaseException:
                process.kill()
                _, status, _ = os.wait4(process.pid, 0)
                process.returncode = os.waitstatus_to_exitcode(status)
                raise
            elapsed = time.perf_counter()-start
        stderr = (directory/"errors").read_text()
        if process.returncode:
            raise RuntimeError("native benchmark failed (%d): %s" % (process.returncode, stderr))
        output = (directory/"output").read_text()
    rss = usage.ru_maxrss if sys.platform != "darwin" else usage.ru_maxrss/1024
    return output, {"elapsed_seconds": elapsed, "user_seconds": usage.ru_utime,
                    "system_seconds": usage.ru_stime, "peak_rss_kib": rss,
                    "stderr": stderr}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", required=True)
    parser.add_argument("--bound", type=int, default=1000000)
    parser.add_argument("--workers", nargs="+", type=int, default=[4])
    parser.add_argument("--kappas", nargs="+", type=int, default=[-1])
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--coefficients", nargs="+", type=int, default=[1, 0, 0, 1, 0, 1])
    references = parser.add_mutually_exclusive_group(required=True)
    references.add_argument("--reference", help="saved native matrix output")
    references.add_argument("--euler-reference", help="smalljac lpdata CSV (p,a1,...,ag)")
    references.add_argument("--sage-reference", action="store_true",
                            help="compute independent Sage Euler factors before timing")
    parser.add_argument("--reference-report", help="report recording the reference's input hash")
    parser.add_argument("--evaluator", choices=["differences", "powers"], default="differences")
    parser.add_argument("--no-fft", action="store_true")
    parser.add_argument("--timeout", type=float, default=7200)
    parser.add_argument("--memory-mib", type=int, default=8192)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.reference_report and not args.reference:
        parser.error("--reference-report requires --reference")
    if args.repeats < 1 or args.bound < 3 or any(w < 1 for w in args.workers) or \
            any(k < -1 or k > 30 for k in args.kappas) or args.memory_mib < 0:
        parser.error("invalid benchmark bounds, worker count, or kappa")
    q = PolynomialRing(ZZ, "x")(args.coefficients)
    if q.degree() < 3 or not q.is_monic() or not q.is_squarefree():
        parser.error("a squarefree monic model of degree at least three is required")
    genus = (q.degree()-1)//2
    discriminant = q.discriminant()

    def input_batch(bound):
        jobs = [(int(p), precision(int(p), genus)) for p in prime_range(3, bound)
                if discriminant % p]
        batch = " ".join(map(str, args.coefficients))+"\n"
        batch += "".join("%d %d\n" % job for job in jobs)
        return jobs, batch

    jobs, batch = input_batch(args.bound)
    binary = str(Path(args.binary).resolve())
    environment = dict(os.environ, HYPELLFROB_RFOREST_EVAL=args.evaluator,
                       HYPELLFROB_RFOREST_FFT="0" if args.no_fft else "1")
    environment.pop("HYPELLFROB_PROFILE", None)
    report = {"bound": args.bound, "jobs": len(jobs), "coefficients": args.coefficients,
              "binary": binary, "binary_sha256": hashlib.sha256(Path(binary).read_bytes()).hexdigest(),
              "input_sha256": hashlib.sha256(batch.encode()).hexdigest(),
              "platform": platform.platform(), "cpu_count": os.cpu_count(),
              "evaluator": args.evaluator, "fft": not args.no_fft, "repeats": args.repeats,
              "memory_limit_mib": args.memory_mib, "runs": [], "summaries": [], "status": "running"}
    forest_backend = "rforest"
    report["strict_rforest"] = True
    expected = None
    euler_expected = None
    if args.euler_reference or args.sage_reference:
        sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
        from hyperell_regulator.frobenius import _lift
    if args.euler_reference:
        # smalljac writes the independent L-polynomial coefficients
        # 1+a1*T+...+ag*T^g+...+p^g*T^(2g). Complete the functional
        # equation and compare to the lifted native characteristic polynomial.
        reference_text = Path(args.euler_reference).read_text()
        euler_expected = {}
        for line in reference_text.splitlines():
            if not line.strip() or line.startswith("["):
                continue
            fields = line.split(",")
            if "?" in fields:
                continue
            values = list(map(int, fields))
            p, a = values[0], [1]+values[1:]
            if len(a) != genus+1:
                raise RuntimeError("wrong genus in Euler-factor reference at p=%d" % p)
            euler_expected[p] = [p**(genus-i)*a[i] for i in range(genus)]+list(reversed(a))
        if any(p not in euler_expected for p, _ in jobs):
            raise RuntimeError("Euler-factor reference does not cover every job")
        report["reference"] = str(Path(args.euler_reference).resolve())
        report["reference_sha256"] = hashlib.sha256(reference_text.encode()).hexdigest()
        report["reference_header"] = reference_text.splitlines()[0]
    if args.sage_reference:
        euler_expected = {}
        print(json.dumps({"event": "sage_reference", "jobs": len(jobs)}), flush=True)
        for p, _ in jobs:
            h = HyperellipticCurve(q.change_ring(GF(p)))
            euler_expected[p] = [int(c) for c in h.frobenius_polynomial()]
        report["reference"] = "Sage Frobenius polynomials"
    report["comparison"] = "euler_factors" if euler_expected is not None else "full_matrices"
    report["full_characteristic_moduli_checked"] = euler_expected is not None
    if args.reference:
        reference_text = Path(args.reference).read_text()
        report["reference"] = str(Path(args.reference).resolve())
        report["reference_sha256"] = hashlib.sha256(reference_text.encode()).hexdigest()
        if args.reference_report:
            previous = json.loads(Path(args.reference_report).read_text())
            assert previous["coefficients"] == args.coefficients
            ref_jobs, ref_batch = input_batch(previous["bound"])
            assert hashlib.sha256(ref_batch.encode()).hexdigest() == previous["input_sha256"]
            reference_rows = parse_rows(reference_text, ref_jobs, genus)
        else:
            reference_rows = [list(map(int, line.split())) for line in reference_text.splitlines()]
        by_prime = {row[0]: row for row in reference_rows}
        expected = [by_prime[p] for p, _ in jobs]
        parse_rows("\n".join(" ".join(map(str, row)) for row in expected), jobs, genus)

    def save():
        Path(args.output).write_text(json.dumps(report, indent=2)+"\n")

    save()
    for workers in args.workers:
        for repeat in range(args.repeats):
            variants = [(forest_backend, k) for k in args.kappas]
            if repeat%2:
                variants.reverse()
            for backend, kappa in variants:
                print(json.dumps({"event": "starting", "backend": backend,
                                  "bound": args.bound, "workers": workers, "kappa": kappa,
                                  "repeat": repeat}), flush=True)
                output, metrics = run(binary, workers, backend, kappa, batch, args.timeout,
                                      args.memory_mib, environment)
                rows = parse_rows(output, jobs, genus)
                if expected is None and euler_expected is None:
                    raise RuntimeError("no reference matrices available")
                compared = 0
                for i, ((p, n), row) in enumerate(zip(jobs, rows)):
                    if row[1]:
                        if euler_expected is not None:
                            dim = 2*genus
                            matrix = [row[3+j*dim:3+(j+1)*dim] for j in range(dim)]
                            raw = sage_matrix(ZZ, matrix).charpoly().list()
                            if any((a-b) % (p**n) for a, b in zip(raw, euler_expected[p])):
                                raise AssertionError("full characteristic polynomial mismatch at p=%d" % p)
                            if _lift(matrix, p, genus, n).list() != euler_expected[p]:
                                raise AssertionError("Euler-factor mismatch at p=%d" % p)
                        elif row != expected[i]:
                            raise AssertionError("matrix mismatch at p=%d" % p)
                        compared += int(row[1] == 1)
                    elif p > q.degree()*(4*n-1) and q[0]%p:
                        raise AssertionError("rforest unexpectedly declined p=%d" % p)
                if not compared:
                    raise AssertionError("rforest supplied no matrices")
                record = {"backend": backend, "workers": workers, "kappa": kappa,
                          "repeat": repeat, "declined_jobs": len(jobs)-compared,
                          "matrices_compared": compared if euler_expected is None else 0,
                          "euler_factors_compared": compared if euler_expected is not None else 0,
                          "all_matrices_equal": True if euler_expected is None else None,
                          "all_euler_factors_equal": True if euler_expected is not None else None,
                          "native_output_sha256": hashlib.sha256(output.encode()).hexdigest(),
                          **metrics}
                report["runs"].append(record)
                print(json.dumps({"event": "completed", **record}), flush=True)
                save()
        for backend, kappa in [(forest_backend, k) for k in args.kappas]:
            selected = [r for r in report["runs"] if r["workers"] == workers and
                        r["backend"] == backend and r["kappa"] == kappa]
            report["summaries"].append({"backend": backend, "workers": workers,
                "kappa": kappa, "median_seconds": statistics.median(r["elapsed_seconds"] for r in selected),
                "median_cpu_seconds": statistics.median(r["user_seconds"]+r["system_seconds"] for r in selected),
                "peak_rss_kib": max(r["peak_rss_kib"] for r in selected)})
    report["status"] = "complete"
    save()


if __name__ == "__main__":
    main()
