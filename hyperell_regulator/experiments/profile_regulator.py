"""Bounded float64 diagnostics for both covers of conductor 249.

Run in Sage's Python with ``python -m hyperell_regulator.experiments.profile_regulator``.
Every integration edge has a wall-clock and a work limit. Progress is printed
while an edge runs, including its active phase and continuation counters.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np

from ..regulator import function_with_divisor, stokes
from ..regulator.cover import Cover
from ..regulator.diagnostics import IntegrationMonitor
from ..regulator.polynomials import ring


def _progress(event):
    state = event.get("state", {})
    print(json.dumps({k: event.get(k) for k in
                      ("event", "degree", "index", "seconds", "steps", "phase")}
                     | {"state": state, "counters": event.get("counters", {}),
                        "status": event.get("status")}), flush=True)


def profile_249(orders=(7, 14), max_edge_seconds=15.0):
    x = ring().gen()
    F = x**6 + 4*x**5 + 4*x**4 + 2*x**3 + 1
    report = {"curve": "249.a.249.1", "covers": []}
    for order in orders:
        point = {7: (-1, 0), 14: (0, -1)}[order]
        cov = Cover(F, *function_with_divisor(F, point, "infinity", order))
        z, kinds, direction = cov.cut_path()
        start = time.perf_counter()
        upper, lower, lam, jumps = stokes.cut_data(cov, z, direction)
        row = {"degree": order, "cut_seconds": time.perf_counter() - start,
               "nodes": [repr(complex(v)) for v in z], "failures": []}
        hold = stokes.hold_radii(z)
        periods, logs = [], []
        with IntegrationMonitor(max_edge_seconds=max_edge_seconds,
                                max_edge_steps=100000, callback=_progress) as monitor:
            for i in range(len(z)):
                try:
                    if i + 1 < len(z):
                        _, leg = stokes.edge_leg(
                            cov, z[i], z[i+1], lam[i], kinds[i], kinds[i+1],
                            hold_a=hold[i], hold_b=hold[i+1])
                    else:
                        _, leg = stokes.ray_leg(
                            cov, z[i], direction, lam[i], kind_in=kinds[i],
                            span=stokes.cut_scale(z)[0], hold=hold[i])
                    periods.append(leg.dI)
                    logs.append(leg.dL)
                except (RuntimeError, ArithmeticError, TimeoutError) as exc:
                    row["failures"].append({"edge": i, "error": str(exc)})
                    print(json.dumps(row["failures"][-1]), flush=True)
            row["edges"] = monitor.records
        if not row["failures"]:
            u, ell = np.stack(periods, axis=1), np.stack(logs, axis=1)
            for name, aa, bb in (("period", u, u), ("log", ell+jumps*u, ell)):
                defect = max(abs(sum(s*(aa[a, i, m] if s > 0 else bb[a, i, m])
                                     for i, m, s in walk))
                             for a in range(cov.g)
                             for walk in stokes.boundary_walks(upper, lower))
                row[name + "_closure"] = float(defect)
        report["covers"].append(row)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orders", nargs="+", type=int, choices=(7, 14), default=[7, 14])
    parser.add_argument("--max-edge-seconds", type=float, default=15.0)
    parser.add_argument("--output", type=Path, default=Path("results/float64_steps_249.json"))
    args = parser.parse_args()
    report = profile_249(args.orders, args.max_edge_seconds)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print("Saved", args.output, flush=True)
    if any(c["failures"] for c in report["covers"]):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
