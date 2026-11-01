"""Single 100-digit comparison for the conductor-555 nonintegral triple."""

from sage.all import QQ, PolynomialRing
from hyperell_regulator.components import HyperellipticCurveWithReduction
from hyperell_regulator.experiments.beilinson import _analyse_triple
from hyperell_regulator.regulator import IntegrationMonitor

x = PolynomialRing(QQ, "x").gen()
C = HyperellipticCurveWithReduction(
    3*x**5 - 2*x**4 - 4*x**3 + x**2 + x,
    x + 1,
)

triple = {
    # Internal keys P,Q are the two zeros; R is their common pole.
    "P": (QQ(-1)/3, QQ(-1)/3),
    "Q": (1, -1),
    "R": ("infinity", 0),
    "order_P": 2,
    "order_Q": 2,
    "verdict_P": False,
    "verdict_Q": False,
    "pairwise_nonconjugate": True,
}

digits = 100


def report_edge(event):
    detail = event["phase"] or event["edge"]
    if event["event"] == "end":
        detail = "%s, %s" % (event["status"], detail)
    print("  555.a.8325.1 edge %d %s: %.1fs, %s" %
          (event["index"], event["event"], event["seconds"], detail),
          flush=True)


with IntegrationMonitor(max_edge_seconds=300, max_edge_steps=1_000_000,
                        callback=report_edge, progress_interval=5):
    result = _analyse_triple(
        {"label": "555.a.8325.1", "cond": 555},
        C, triple, n_triples=1,
        reg_digits=digits,
        digits=digits,
        lindep_bound=digits // 2,
        verbose=True,
    )

if result["status"] != "ok":
    raise RuntimeError(result.get("error", result["status"]))

print(result["comparison"])
