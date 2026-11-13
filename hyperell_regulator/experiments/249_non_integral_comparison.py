"""Single 100-digit comparison for the conductor-249 nonintegral triple."""
from sage.all import QQ, PolynomialRing
from hyperell_regulator.components import HyperellipticCurveWithReduction
from hyperell_regulator.experiments.beilinson import _analyse_triple
from hyperell_regulator.regulator import IntegrationMonitor

x = PolynomialRing(QQ, "x").gen()
C = HyperellipticCurveWithReduction(
    -x**5 + x**3 + x**2 + 3*x + 2,
    x**3 + 1,
)

triple = {
    # Internal keys P,Q are the two zeros; R is their common pole.
    "P": (QQ(-1), QQ(0)),
    "Q": (QQ(0), QQ(1)),
    "R": ("infinity", 0),
    "order_P": 7,
    "order_Q": 28,
    "verdict_P": True,
    "verdict_Q": False,
    "pairwise_nonconjugate": True,
}

digits = 100


def report_edge(event):
    detail = event["phase"] or event["edge"]
    if event["event"] == "end":
        detail = "%s, %s" % (event["status"], detail)
    print("  249.a.6723.1 edge %d %s: %.1fs, %s" %
          (event["index"], event["event"], event["seconds"], detail),
          flush=True)


# Both logarithmic columns share a degree-two x-cut.
with IntegrationMonitor(max_edge_seconds=120, max_edge_steps=1_000_000,
                        callback=report_edge, progress_interval=5):
    result = _analyse_triple(
        {"label": "249.a.6723.1", "cond": 249},
        C, triple, n_triples=1,
        reg_digits=digits,
        digits=digits,
        lindep_bound=digits // 2,
        verbose=True,
    )

if result["status"] != "ok":
    raise RuntimeError(result.get("error", result["status"]))

print(result["comparison"])
