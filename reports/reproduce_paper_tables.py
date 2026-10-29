"""Audit the three tables in the manuscript's final section with this checkout.

Run with Sage's Python; see paper_reproducibility_20261027.md for the invocation.
Only input curve/component data and the explicitly identified L-value cache
are reused. Regulator output is always forced into a separate data directory.
"""
import argparse
import concurrent.futures
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys


def load(path):
    return json.loads(Path(path).read_text())


def save(path, data):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(data, indent=2, default=str) + "\n")


def table(source, label):
    end = source.index(r"\label{" + label + "}")
    start = source.rindex(r"\begin{tabular}", 0, end)
    return source[start:source.index(r"\end{tabular}", start)]


def entries(body):
    for row in re.split(r"\\\\\s*\n", body):
        match = re.search(r"\\texttt\{([^}]+)\}", row)
        if match:
            yield match[1], row.split("&")[1:]


def clean(cell):
    return cell.strip().strip("$").strip()


def rational(text):
    text = clean(text)
    return re.sub(r"\\frac\{([^}]+)\}\{([^}]+)\}", r"\1/\2", text)


def parse_paper(path):
    source = "\n".join(line for line in Path(path).read_text().splitlines()
                       if not line.lstrip().startswith("%"))
    rows = {}
    for label, cells in entries(table(source, "tab:torsion-integrality")):
        points = re.findall(r"\(([^,()]+),([^,()]+)\)", cells[0])
        orders = re.search(r"\\begin\{gathered\}(.*?)\\end\{gathered\}", cells[1])[1]
        rows[label] = dict(label=label,
                           points=[[rational(a), rational(b)] for a, b in points],
                           orders=[int(v) for v in orders.split(r"\\")],
                           q=rational(cells[3]), reported_error=clean(cells[4]))
    for label, cells in entries(table(source, "tab:curve-models")):
        rows[label]["polynomial_latex"] = clean(cells[0])
        rows[label]["reciprocal"] = label in {
            "709.a.709.1", "763.a.763.1", "1123.a.1123.1", "1147.a.35557.2"}
    for label, cells in entries(table(source, "tab:data")):
        rows[label]["numerical"] = dict(zip(
            ("det", "c_plus", "c_minus", "L_second_derivative"),
            [clean(v) for v in cells[:4]]))
    return list(rows.values())


def polynomial(latex):
    from sage.all import PolynomialRing, QQ
    x = PolynomialRing(QQ, "x").gen()
    text = re.sub(r"[\s{}]", "", latex)
    terms = re.findall(r"[+-]?[^+-]+", text)
    out = x.parent()(0)
    for term in terms:
        m = re.fullmatch(r"([+-]?)(\d*)(x(?:\^(\d+))?)?", term)
        if m is None or not (m[2] or m[3]):
            raise ValueError("unrecognized polynomial term: " + term)
        coefficient = int(m[2] or 1) * (-1 if m[1] == "-" else 1)
        out += coefficient * x**(int(m[4] or 1) if m[3] else 0)
    return out


def exact_checks(expected, row):
    from sage.all import QQ, ZZ
    from hyperell_regulator.components import curve_from_row
    from hyperell_regulator.regulator.functions import _field, _places
    C = curve_from_row(row)
    G = polynomial(expected["polynomial_latex"])
    original = G.parent()(C.F)
    transformed = G.parent()(list(reversed(original.list() + [0]*(7-len(original.list())))))
    model_ok = G == (transformed if expected["reciprocal"] else original)
    K, L = _field(G)
    infinities = _places(K, L, "infinity")
    x, y = L(K.gen()), L.gen()
    if G.degree() == 6:
        plus = [q for q in infinities if (y/x**3 - 1).valuation(q) > 0]
        assert len(plus) == 1
        base = plus[0]
    else:
        assert len(infinities) == 1
        base = infinities[0]
    checks = []
    for point, claimed in zip(expected["points"], expected["orders"]):
        point = tuple(QQ(a) for a in point)
        on_curve = point[1]**2 == G(point[0])
        item = dict(point=point, on_curve=bool(on_curve), claimed_order=claimed)
        if on_curve:
            p = _places(K, L, point)[0]
            # The cached finite group exponent is a proven multiple in rank 0.
            exponent = ZZ(1)
            for n in row["mw_invs"]:
                exponent = exponent.lcm(ZZ(n))
            for n in exponent.divisors():
                if len((n*base.divisor()-n*p.divisor()).basis_function_space()) == 1:
                    item["actual_order_to_plus_infinity"] = int(n)
                    break
        checks.append(item)
    from hyperell_regulator.regulator import c_invariants
    cp, cm = c_invariants(G, prec=110)
    return dict(model_matches_coordinate_change=bool(model_ok), points=checks,
                period_target_digits=110,
                paper_periods={"c_plus": cp.str(truncate=False),
                               "c_minus": cm.str(truncate=False)})


def worker(args):
    from hyperell_regulator.experiments import beilinson as b
    expected = load(args.worker)
    row = load(Path(os.environ["HYPERELL_REGULATOR_DATA"])/"lmfdb_cache/g2c.json")[0]
    try:
        checks = exact_checks(expected, row)
    except Exception as exc:
        checks = {"error": repr(exc)}
    save(Path(args.worker).with_name("exact.json"), checks)
    b.run(1148, force=True, reg_digits=args.prec, digits=args.digits,
          max_edge_seconds=args.edge_seconds, verbose=True,
          out=str(Path(args.worker).with_name("experiment.json")))


def component_checks(expected, row):
    """Recompute the two manuscript point/base pairs on an affine model."""
    from sage.all import QQ
    from hyperell_regulator.components import (
        affine_model, curve_from_row, same_component_all_primes)
    C = curve_from_row(row)
    G = polynomial(expected["polynomial_latex"])
    points = []
    for u, v in expected["points"]:
        u, v = QQ(u), QQ(v)
        if expected["reciprocal"]:
            point = (("infinity", (v-C.h[3])/2) if u == 0 else
                     (1/u, (v/u**3-C.h(1/u))/2))
        else:
            point = (u, (v-C.h(u))/2)
        points.append(point)
    if expected["reciprocal"]:
        base = (QQ(0), -QQ(C.h[0])/2)
    else:
        base = ("infinity", (QQ(1 if G.degree() == 6 else 0)-QQ(C.h[3]))/2)
    D, images, shift = affine_model(C, points+[base])
    pairs = []
    for point, image in zip(points, images[:2]):
        verdict, per = same_component_all_primes(D, image, images[2])
        pairs.append(dict(point=point, base=base, same_component_everywhere=verdict,
                          per_prime={str(k): v for k, v in per.items()}))
    return dict(pairs=pairs)


def compare(expected, checks, result):
    """Convert changes of pole, logarithmic columns and reciprocal coordinate."""
    from sage.all import RealField, QQ
    R = RealField(420)
    answer = dict(label=expected["label"], status=result["status"])
    if "paper_periods" in checks:
        answer["periods_round_as_printed"] = {
            k: rounds_as_printed(v, expected["numerical"][k])
            for k, v in checks["paper_periods"].items()}
    if result["status"] != "ok":
        answer["error"] = result.get("error")
        return answer

    def paper_x(point):
        x = None if point[0] == "infinity" else QQ(point[0])
        if expected["reciprocal"]:
            return QQ(0) if x is None else (None if x == 0 else 1/x)
        return x

    marked = [QQ(p[0]) for p in expected["points"]]

    def basis(point):
        x = paper_x(point)
        if x is None:
            return [0, 0]
        if x not in marked:
            raise ValueError("experiment selected an additional x-coordinate")
        return [int(x == p) for p in marked]

    t = result["triple"]
    base = basis(t["R"])
    columns = [[a-b for a, b in zip(basis(t[k]), base)] for k in ("P", "Q")]
    log_det = columns[0][0]*columns[1][1]-columns[0][1]*columns[1][0]
    moment_det = -1 if (result.get("shift") is not None) != expected["reciprocal"] else 1
    assert abs(log_det) == 1 and result.get("pole", "infinity") == "infinity"
    determinant = R(result["det"])/(log_det*moment_det)
    answer.update(column_change_determinant=log_det, moment_change_determinant=moment_det,
                  paper_determinant=str(determinant),
                  determinant_rounds_as_printed=rounds_as_printed(determinant, expected["numerical"]["det"]))
    if result.get("L_second_derivative") is not None and "paper_periods" in checks:
        L = R(result["L_second_derivative"])
        cp, cm = [R(checks["paper_periods"][k]) for k in ("c_plus", "c_minus")]
        norm = 2*R.pi()*cp*cm*L
        answer.update(paper_q=str(norm/determinant),
                      relative_error_for_paper_q=str(abs(norm/(R(QQ(expected["q"]))*determinant)-1)),
                      lvalue_rounds_as_printed=rounds_as_printed(L, expected["numerical"]["L_second_derivative"]),
                      L_cached=result.get("L_cached"))
    return answer


def rounds_as_printed(value, printed):
    from decimal import Decimal, localcontext
    with localcontext() as ctx:
        ctx.prec = 160
        target = Decimal(printed)
        return abs(Decimal(str(value))-target) <= Decimal(5).scaleb(target.as_tuple().exponent-1)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--paper", type=Path)
    p.add_argument("--out-dir", type=Path, default=Path("results/paper_reproducibility_20261027"))
    p.add_argument("--prec", type=int, default=100)
    p.add_argument("--digits", type=int, default=100)
    p.add_argument("--edge-seconds", type=float, default=180)
    p.add_argument("--jobs", type=int, default=8)
    p.add_argument("--worker", type=Path, help=argparse.SUPPRESS)
    p.add_argument("--summarize", action="store_true")
    p.add_argument("--components-only", action="store_true")
    args = p.parse_args()
    if args.components_only:
        rows = {r["label"]: r for r in load("lmfdb_cache/g2c.json")}
        completed = {}
        for expected in load(args.out_dir/"expected.json"):
            label = expected["label"]
            try:
                completed[label] = component_checks(expected, rows[label])
            except Exception as exc:
                completed[label] = {"error": repr(exc)}
            save(args.out_dir/"fresh_components.json", completed)
            print(label, completed[label], flush=True)
        return
    if args.summarize:
        root = args.out_dir
        expected = {r["label"]: r for r in load(root/"expected.json")}
        items = load(root/"audit.json")
        checks = (load(root/"exact_periods.json") if (root/"exact_periods.json").exists()
                  else {item["label"]: item["exact"] for item in items})
        comparisons = []
        for item in items:
            label = item["label"]
            result = item["experiment"]
            retry = root/label/"retry_100.json"
            if retry.exists() and load(retry)[0].get("status") == "ok":
                result = load(retry)[0]
            comparison = compare(expected[label], checks[label], result)
            comparison["initial_status"] = item["experiment"]["status"]
            comparison["edge_seconds"] = result["settings"]["max_edge_seconds"]
            comparisons.append(comparison)
        save(root/"comparisons.json", comparisons)
        return
    if args.worker:
        worker(args)
        return
    if args.paper is None:
        p.error("--paper is required")
    expected = parse_paper(args.paper)
    root = args.out_dir.resolve()
    save(root/"expected.json", expected)
    save(root/"provenance.json", dict(
        commit=subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
        paper_sha256=hashlib.sha256(args.paper.read_bytes()).hexdigest(),
        source_sha256={str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                       for path in sorted(Path("hyperell_regulator").rglob("*.py"))},
        input_sha256={str(path): hashlib.sha256(Path(path).read_bytes()).hexdigest()
                      for path in ("lmfdb_cache/g2c.json", "results/same_component_results.json",
                                   "results/lvalue_cache.json")},
        regulator_digits=args.prec, lvalue_digits=args.digits,
        edge_seconds=args.edge_seconds, jobs=args.jobs,
        reused_inputs=["lmfdb_cache/g2c.json", "results/same_component_results.json",
                       "results/lvalue_cache.json"],
        cached_regulators_reused=False))
    curves = {r["label"]: r for r in load("lmfdb_cache/g2c.json")}
    components = {r["label"]: r for r in load("results/same_component_results.json")}

    def run_one(row):
        label = row["label"]
        work = root/label
        save(work/"expected.json", row)
        save(work/"lmfdb_cache/g2c.json", [curves[label]])
        save(work/"results/same_component_results.json", [components[label]])
        shutil.copyfile("results/lvalue_cache.json", work/"results/lvalue_cache.json")
        env = dict(os.environ, HYPERELL_REGULATOR_DATA=str(work),
                   OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
        with (work/"run.log").open("w") as log:
            command = [sys.executable, "-u", "-m", "reports.reproduce_paper_tables",
                       "--worker", str(work/"expected.json"), "--prec", str(args.prec),
                       "--digits", str(args.digits), "--edge-seconds", str(args.edge_seconds)]
            status = subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT).returncode
        print(label, "exit", status, flush=True)
        return dict(label=label, exit_code=status,
                    exact=load(work/"exact.json") if (work/"exact.json").exists() else None,
                    experiment=load(work/"experiment.json")[0] if (work/"experiment.json").exists() else None)

    completed = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.jobs) as pool:
        for future in concurrent.futures.as_completed([pool.submit(run_one, r) for r in expected]):
            completed.append(future.result())
            save(root/"audit.json", sorted(completed, key=lambda r: list(curves).index(r["label"])))


if __name__ == "__main__":
    main()
