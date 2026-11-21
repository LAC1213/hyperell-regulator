"""Explicit local standard factors for the CM-labelled survey below 10000.

Coefficients are ascending, in the Euler polynomial variable T.  The
``StandardLFunction.bad_factors`` API takes wedge-square factors, so the
helper multiplies each standard polynomial by (1-p*T).  Wild conductor
exponents are supplied separately.  See reports/beilinson_cm_bad_factors.rst
for the local-representation derivation, including what ``CM`` means here.

Overrides are guarded by the exact cached equation y^2+h(x)y=f(x), rather
than applied to arbitrary models carrying the same label.
"""

import hashlib
import json

from sage.all import PolynomialRing, QQ


LOCAL_DATA = {
    "324.a.648.1": {
        "model": [[0, 0, 1, 2, 2, 1], [1, 1, 0, 1]],
        "standard_factors": {"3": [1, -3]},
        "conductor_exponents": {"3": 8},
        "newform": "18.2.c.a",
        "derivation": "Principal series at 3, character conductors (0,2).",
    },
    "784.c.614656.1": {
        "model": [[0, -1, -9, -13, -4, 1], []],
        "standard_factors": {"2": [1, 4, 8, 8]},
        "conductor_exponents": {"2": 2},
        "newform": "28.2.e.a",
        "derivation": "Tame supercuspidal at 2; nebentypus value has order 3.",
    },
    "1296.a.20736.1": {
        "model": [[0, -1, 4, -3, -1, 1], []],
        "standard_factors": {"2": [1, 4, 8, 8], "3": [1, -3]},
        "conductor_exponents": {"2": 2, "3": 8},
        "newform": "36.2.e.a",
        "derivation": "Tame supercuspidal at 2; principal series (0,2) at 3.",
    },
    "2704.a.43264.1": {
        "model": [[0, -1, -5, -5, 0, 1], []],
        "standard_factors": {"2": [1, 0, 0, 8]},
        "conductor_exponents": {"2": 2},
        "newform": "52.2.h.a",
        "derivation": "Tame supercuspidal at 2; nebentypus value has order 6.",
    },
    "3969.c.35721.1": {
        "model": [[0, -1, 0, 4, -5, 1], [0, 1, 1]],
        "standard_factors": {"3": [1, -3]},
        "conductor_exponents": {"3": 8},
        "newform": "63.2.i.a",
        "derivation": "Principal series at 3, character conductors (0,2).",
    },
    "4096.b.65536.1": {
        "model": [[0, -1, 0, 0, 0, 1], []],
        "standard_factors": {"2": [1, 2]},
        "conductor_exponents": {"2": 12},
        "newform": "64.2.b.a",
        "derivation": "Bolza CM decomposition: weight-3 level-8 form and three quadratic Tate twists.",
    },
    "8649.c.700569.1": {
        "model": [[0, -1, 4, 13, 9, 1], [0, 1, 1]],
        "standard_factors": {"3": [1, -3]},
        "conductor_exponents": {"3": 4},
        "newform": "93.2.g.b",
        "derivation": "Tame principal series at 3, character conductors (0,1).",
    },
}


def local_data_metadata(label):
    """JSON-safe local data and fingerprint used to reject stale caches."""
    if label not in LOCAL_DATA:
        return None
    data = json.loads(json.dumps(LOCAL_DATA[label]))
    data["source"] = "https://www.lmfdb.org/Genus2Curve/Q/" + label.replace(".", "/")
    digest = hashlib.sha256(json.dumps(data, sort_keys=True).encode()).hexdigest()
    return dict(data, revision=1, sha256=digest)


def lfunction_options(label, C):
    """Translate the guarded standard data to the wedge-square constructor API."""
    data = LOCAL_DATA.get(label)
    if data is None:
        return {}
    R = PolynomialRing(QQ, "T")
    if R(C.f.list()) != R(data["model"][0]) or R(C.h.list()) != R(data["model"][1]):
        raise ValueError("manual local data for %s require its recorded model" % label)
    T = R.gen()
    return {
        "bad_factors": {int(p): (1-int(p)*T)*R(a)
                        for p, a in data["standard_factors"].items()},
        "conductor_exponents": {int(p): a for p, a in data["conductor_exponents"].items()},
    }
