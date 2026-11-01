"""Time native rforest/legacy runs and compare full Frobenius matrices.

Run with sage -python. Uses wait4 for per-process CPU time and peak RSS;
wall time is measured with perf_counter. The strict rforest mode cannot
silently pass a comparison by falling back to the old implementation.
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

from sage.all import PolynomialRing, ZZ, prime_range
from benchmark import precision


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
    command = [binary, str(workers), backend]
    if backend.startswith("rforest"):
        command.append(str(kappa))

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
    parser.add_argument("--reference", help="saved legacy output; omit legacy timing")
    parser.add_argument("--reference-report", help="report recording the reference's input hash")
    parser.add_argument("--evaluator", choices=["differences", "powers"], default="differences")
    parser.add_argument("--no-fft", action="store_true")
    parser.add_argument("--timeout", type=float, default=7200)
    parser.add_argument("--memory-mib", type=int, default=8192)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    if args.repeats < 1 or args.bound < 3 or any(w < 1 for w in args.workers) or \
            any(k < -1 or k > 30 for k in args.kappas) or args.memory_mib < 0:
        parser.error("invalid benchmark bounds, worker count, or kappa")
    q = PolynomialRing(ZZ, "x")(args.coefficients)
    if q.degree() < 3 or q.degree()%2 == 0 or not q.is_monic() or not q.is_squarefree():
        parser.error("a squarefree monic odd degree model is required")
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
    expected = None
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
            variants = ([] if args.reference else [("legacy", -1)])
            variants += [("rforest-only", k) for k in args.kappas]
            if repeat%2:
                variants.reverse()
            for backend, kappa in variants:
                print(json.dumps({"event": "starting", "backend": backend,
                                  "bound": args.bound, "workers": workers, "kappa": kappa,
                                  "repeat": repeat}), flush=True)
                output, metrics = run(binary, workers, backend, kappa, batch, args.timeout,
                                      args.memory_mib, environment)
                rows = parse_rows(output, jobs, genus)
                if backend == "legacy" and expected is None:
                    expected = rows
                if expected is None:
                    raise RuntimeError("no reference matrices available")
                compared = 0
                for (p, n), row, ref in zip(jobs, rows, expected):
                    if backend == "legacy" or row[1]:
                        if row != ref:
                            raise AssertionError("matrix mismatch at p=%d" % p)
                        compared += int(row[1] == 1)
                    elif p > q.degree()*(4*n-1) and q[0]%p:
                        raise AssertionError("rforest unexpectedly declined p=%d" % p)
                if backend != "legacy" and not compared:
                    raise AssertionError("rforest supplied no matrices")
                record = {"backend": backend, "workers": workers, "kappa": kappa,
                          "repeat": repeat, "matrices_compared": compared,
                          "all_matrices_equal": True, **metrics}
                report["runs"].append(record)
                print(json.dumps({"event": "completed", **record}), flush=True)
                save()
        for backend, kappa in (([("legacy", -1)] if not args.reference else [])+
                               [("rforest-only", k) for k in args.kappas]):
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
