"""Precompute the manuscript's 16 curves / 32 functions for GitHub Pages.

    sage -python -m hyperell_regulator.experiments.paper_explorer --jobs 4

Input is a portable snapshot of the manuscript tables, not the experiment's
automatically selected triples. Completed curves are resumable. No L-values
are computed: the paper's displayed numbers are identified as reference data.
"""

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time


def save(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(".tmp")
    temporary.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    temporary.replace(path)


def algorithm_sources():
    """Hashes make resumed and published calculations traceable to the code."""
    package = Path(__file__).resolve().parents[1]
    paths = sorted((package/"regulator").glob("*.py")) + [Path(__file__).resolve()]
    return {str(path.relative_to(package.parent)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in paths}


def function_at_marked_point(F, point, order):
    """Construct phi with its pole specifically at O+ on even-degree models."""
    from sage.all import QQ
    from hyperell_regulator.regulator.functions import _field, _places, _in_x

    K, L = _field(F)
    infinities = _places(K, L, "infinity")
    if F.degree() % 2 == 0:
        x, y = L(K.gen()), L.gen()
        infinities = [q for q in infinities if (y/x**(F.degree()//2) - 1).valuation(q) > 0]
    if len(infinities) != 1:
        raise ValueError("Cannot identify the manuscript's O+ base point")
    p = _places(K, L, tuple(QQ(v) for v in point))[0]
    q = infinities[0]
    basis = (order*q.divisor() - order*p.divisor()).basis_function_space()
    if len(basis) != 1:
        raise ValueError("The marked torsion divisor is not principal")
    A, B = basis[0].list()
    return _in_x(A), _in_x(B)


def compute_curve(entry, out, prec, edge_seconds):
    from sage.all import PolynomialRing, QQ, RealField, matrix
    from hyperell_regulator.regulator.cover import Cover
    from hyperell_regulator.regulator.diagnostics import IntegrationMonitor
    from hyperell_regulator.regulator.explorer import cut_system_data, parse_expression

    S = PolynomialRing(QQ, "x")
    F = S(parse_expression(entry["F"], S.fraction_field()))
    label = entry["label"]
    source_hashes = algorithm_sources()
    functions = []
    R = RealField(int(3.322*(prec+15)))
    for name, point, n in zip(("P", "Q"), entry["points"], entry["orders"]):
        filename = "%s-%s.json" % (label, name)
        path = out / "data" / filename
        if path.is_file():
            data = json.loads(path.read_text())
            if (data["precision"] != prec or data.get("paper", {}).get("entry") != entry
                    or data.get("source_sha256") != source_hashes):
                raise ValueError("Cached inputs, precision or source changed; choose a new output directory")
        else:
            print("%s phi_%s: N=%d, target=%d" % (label, name, n, prec), flush=True)
            phi = function_at_marked_point(F, point, n)
            with IntegrationMonitor(max_edge_seconds=edge_seconds, max_edge_steps=1000000):
                data = cut_system_data(F, phi, prec)
            cov = Cover(F, *phi)
            num, den = cov.norm_exact()
            c = (S.gen() - QQ(point[0]))**n * den/num
            if c.numerator().degree() > 0 or c.denominator().degree() > 0:
                raise ArithmeticError("Norm normalization is not constant")
            logc = R(abs(QQ(c))).log()
            for pair in data["pairings"]:
                pair["normalized_regulator"] = str((2*R(pair["regulator"]) + logc*R(pair["period"]))/n)
            data["paper"] = dict(entry=entry, function=name, point=point, order=n,
                                 pole="O+" if F.degree() % 2 == 0 else "O", norm_c=str(c),
                                 normalization="normalized_regulator = (2*regulator + log|c|*period)/N; (x-x(P))^N = c*Norm(phi)")
            data["source_sha256"] = source_hashes
            save(path, data)
        functions.append(dict(name=name, order=n, point=point, file="data/"+filename,
                              warnings=data["warnings"], data=data))
    p, q = [f.pop("data") for f in functions]
    agreement = max(abs(R(a["period"])-R(b["period"])) for a, b in zip(p["pairings"], q["pairings"]))
    scale = max(R(1), *(abs(R(a["period"])) for a in p["pairings"]))
    if agreement/scale > R("1e-10"):
        raise ArithmeticError("Ordinary periods from the two covers disagree")
    rows = [[R(a["normalized_regulator"]), R(b["normalized_regulator"]), R(a["period"])]
            for a, b in zip(p["pairings"], q["pairings"])]
    determinant = matrix(R, rows).determinant()
    from decimal import Decimal, localcontext
    with localcontext() as context:
        context.prec = 80
        printed = Decimal(entry["numerical"]["det"])
        printed_match = abs(Decimal(str(determinant)) - printed) <= Decimal(5).scaleb(printed.as_tuple().exponent-1)
    if not printed_match:
        raise ArithmeticError("Fresh determinant %s disagrees with paper %s" % (determinant, printed))
    return dict(**entry, functions=functions, determinant=str(determinant),
                normalized_matrix=[[str(v) for v in row] for row in rows],
                validation=dict(period_agreement=str(agreement), determinant_rounds_as_printed=printed_match))


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, default=Path("reports/paper_explorer_inputs.json"))
    parser.add_argument("--out", type=Path, default=Path("docs/_extra/explorer"))
    parser.add_argument("--prec", type=int, default=15)
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--edge-seconds", type=float, default=120)
    parser.add_argument("--timeout", type=int, default=1200)
    parser.add_argument("--worker", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args(argv)
    source = json.loads(args.inputs.read_text())
    if args.worker is not None:
        entry = source["curves"][args.worker]
        result = compute_curve(entry, args.out, args.prec, args.edge_seconds)
        save(args.out / "data" / (entry["label"] + "-summary.json"), result)
        print("Validated " + entry["label"], flush=True)
        return
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    start = time.monotonic()
    def run(index):
        label = source["curves"][index]["label"]
        log = args.out / "data" / (label + ".log")
        log.parent.mkdir(parents=True, exist_ok=True)
        command = [sys.executable, "-m", "hyperell_regulator.experiments.paper_explorer",
                   "--inputs", str(args.inputs), "--out", str(args.out), "--prec", str(args.prec),
                   "--edge-seconds", str(args.edge_seconds), "--worker", str(index)]
        environment = dict(os.environ, OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
        with log.open("w") as stream:
            result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT,
                                    timeout=args.timeout, env=environment, check=False)
        if result.returncode:
            raise RuntimeError(label + ": " + log.read_text()[-2500:])
        log.unlink()
        print("Completed %s (%.1fs elapsed)" % (label, time.monotonic()-start), flush=True)
        return json.loads((args.out / "data" / (label + "-summary.json")).read_text())
    completed, errors = {}, []
    with ThreadPoolExecutor(max_workers=args.jobs) as executor:
        futures = {executor.submit(run, i): i for i in range(len(source["curves"]))}
        for future in as_completed(futures):
            try:
                result = future.result()
                completed[result["label"]] = result
            except Exception as exc:
                errors.append(str(exc))
                print(str(exc), file=sys.stderr, flush=True)
    if errors:
        raise SystemExit("%d curves need attention; successful exports are saved for resuming" % len(errors))
    manifest = dict(schema_version=1, display_digits=10, computation_digits=args.prec,
                    source=source["source"], input_sha256=hashlib.sha256(args.inputs.read_bytes()).hexdigest(),
                    source_sha256=algorithm_sources(),
                    curves=[completed[e["label"]] for e in source["curves"]])
    save(args.out / "catalog.json", manifest)
    print("Saved catalog: %d curves, %d functions" % (len(completed), 2*len(completed)), flush=True)


if __name__ == "__main__":
    main()
