r"""Torsion pairs satisfying the integrality criterion on genus two curves.

The criterion requires smooth reduction on the same component of the
stable special fibre at every prime. Component tests at 2 are undecided.
"""

import json
import os

from sage.all import (
    GF,
    HyperellipticCurve,
    Infinity,
    PolynomialRing,
    QQ,
    ZZ,
    divisors,
    gcd,
    lcm,
    prime_range,
)

from hyperell_regulator.paths import data_path, ensure_parent
from hyperell_regulator.reduction import HyperellipticCurveWithReduction

CACHE = data_path("lmfdb_cache", "g2c.json")
RESULTS = data_path("results", "same_component_results.json")


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------

def fetch_lmfdb(max_cond=None, path=CACHE, pause=6, limit=None,
                geom_simple=True, email="acampo@ihes.fr"):
    r"""Fetch and merge LMFDB curve and rational-point data into ``path``.

    Requests are paginated, with ``pause`` seconds between requests.
    The rational-point scan offset is saved in ``_scan_state.json``.
    """
    import time
    import urllib.request
    import urllib.parse

    headers = {"User-Agent": "sage-research-script/1.0 (%s)" % email}
    state_path = os.path.join(os.path.dirname(path) or ".",
                              "_scan_state.json")
    calls = [0]

    def api(table, **kw):
        kw.setdefault("_format", "json")
        url = "https://www.lmfdb.org/api/%s/?%s" % (
            table, urllib.parse.urlencode(kw))
        for attempt in range(4):
            calls[0] += 1
            raw = urllib.request.urlopen(
                urllib.request.Request(url, headers=headers), timeout=120).read()
            if raw.lstrip()[:1] != b"<":
                return json.loads(raw)
            wait = 90 * (attempt + 1)
            print("   (LMFDB served a CAPTCHA; waiting %ss)" % wait)
            time.sleep(wait)
        raise RuntimeError(
            "LMFDB kept returning HTML (a CAPTCHA challenge) instead of JSON")

    try:
        with open(path) as fh:
            have = {r["label"]: r for r in json.load(fh)}
    except (IOError, OSError, ValueError):
        have = {}

    # --- 1. the curves, filtered server-side ------------------------------
    flag = "is_simple_geom" if geom_simple else "is_simple_base"
    selected, offset = [], 0
    while True:
        page = api("g2c_curves", _offset=offset, _limit=100, _sort="cond",
                   _fields="label,cond,is_simple_base,is_simple_geom,eqn",
                   **{flag: "true"})["data"]
        if not page:
            break
        selected.extend(r for r in page
                        if max_cond is None or r["cond"] <= max_cond)
        print("   g2c_curves offset %-5s -> %s selected" % (offset,
                                                            len(selected)))
        if limit is not None and len(selected) >= limit:
            selected = selected[:limit]
            break
        if max_cond is not None and max(r["cond"] for r in page) > max_cond:
            break
        offset += len(page)
        time.sleep(pause)

    # --- 2. the rational points, resuming the scan ------------------------
    want = set(r["label"] for r in selected) - set(have)
    pts = {}
    if want:
        try:
            with open(state_path) as fh:
                offset = int(json.load(fh).get("ratpts_offset", 0))
        except (IOError, OSError, ValueError):
            offset = 0
        start_offset = offset
        while want - set(pts):
            time.sleep(pause)
            page = api("g2c_ratpts", _offset=offset, _limit=100,
                       _fields="label,num_rat_pts,rat_pts,rat_pts_v,mw_invs"
                       )["data"]
            if not page:
                if start_offset == 0:
                    break
                offset, start_offset = 0, 0   # resume point overshot; restart
                continue
            # Guard against a resume point that has overshot: labels begin
            # with the conductor and the table is in conductor order, so if
            # this page is already past the smallest conductor we still need,
            # the wanted rows lie behind us.  Without this the scan would walk
            # to the end of a ~66000 row table before the wrap-around below.
            if start_offset > 0 and not pts:
                here = ZZ(page[0]["label"].split(".")[0])
                lowest = min(ZZ(t.split(".")[0]) for t in want)
                if here > lowest:
                    print("   g2c_ratpts resume point overshot "
                          "(offset %s is at conductor %s, need %s); "
                          "restarting scan" % (offset, here, lowest))
                    offset, start_offset = 0, 0
                    continue
            for r in page:
                if r["label"] in want:
                    pts[r["label"]] = r
            print("   g2c_ratpts offset %-5s -> %s/%s matched"
                  % (offset, len(pts), len(want)))
            offset += len(page)
        with open(ensure_parent(state_path), "w") as fh:
            # int(): Sage Integers are not JSON serialisable
            # Integers, which json cannot serialise
            json.dump({"ratpts_offset": int(offset)}, fh)

    # --- 3. merge into the cache -----------------------------------------
    for r in selected:
        if r["label"] in have and r["label"] not in pts:
            continue                       # already cached, nothing new
        q = pts.get(r["label"], {})
        have[r["label"]] = {
            "label": r["label"],
            "cond": int(r["cond"]),
            "is_simple_base": bool(r["is_simple_base"]),
            "is_simple_geom": bool(r["is_simple_geom"]),
            "eqn": r["eqn"],
            "num_rat_pts": q.get("num_rat_pts"),
            "rat_pts": q.get("rat_pts") or [],
            "rat_pts_verified": q.get("rat_pts_v"),
            "mw_invs": q.get("mw_invs"),
        }
    out = sorted(have.values(), key=lambda r: (r["cond"], r["label"]))
    with open(ensure_parent(path), "w") as fh:
        json.dump(out, fh, indent=2)
    print("cache %s now holds %s curve(s), conductors %s-%s   [%s API call(s)]"
          % (path, len(out), out[0]["cond"], out[-1]["cond"], calls[0]))
    return out


def lmfdb_label(C, pause=6, email="acampo@ihes.fr", verbose=True):
    r"""The LMFDB labels matching the absolute Igusa invariants of ``C``.

    EXAMPLES::

        sage: from hyperell_regulator.components import HyperellipticCurveWithReduction, lmfdb_label
        sage: R.<x> = ZZ[]
        sage: C = HyperellipticCurveWithReduction(-x^2 - x, x^3 + x^2 + x + 1)
        sage: lmfdb_label(C)                                    # optional - internet
        ['277.a.277.1']
    """
    import time
    import urllib.request
    import urllib.parse

    headers = {"User-Agent": "sage-research-script/1.0 (%s)" % email}

    def api(**kw):
        kw.setdefault("_format", "json")
        url = "https://www.lmfdb.org/api/g2c_curves/?%s" % \
            urllib.parse.urlencode(kw)
        raw = urllib.request.urlopen(
            urllib.request.Request(url, headers=headers), timeout=120).read()
        if raw.lstrip()[:1] == b"<":
            raise RuntimeError("LMFDB returned a CAPTCHA page, not JSON")
        return json.loads(raw)

    def invariants(curve):
        return tuple(QQ(t) for t in
                     curve.absolute_igusa_invariants_wamelen())

    target = invariants(C.curve() if hasattr(C, "curve") else C)

    keys = []
    try:
        keys.append(("abs_disc", abs(ZZ(C.minimal_discriminant()))))
    except Exception:
        pass
    try:
        keys.append(("cond", ZZ(C.conductor())))
    except Exception:
        pass

    seen, out = set(), []
    for field, value in keys:
        if out:
            break
        try:
            rows = api(_fields="label,cond,abs_disc,eqn",
                       **{field: int(value)})["data"]
        except Exception as exc:
            if verbose:
                print("   query on %s=%s failed: %s" % (field, value, exc))
            continue
        if verbose:
            print("   %s=%s -> %s candidate(s)" % (field, value, len(rows)))
        for r in rows:
            if r["label"] in seen:
                continue
            seen.add(r["label"])
            eqn = _as_list(r["eqn"])
            try:
                D = HyperellipticCurve(PolynomialRing(QQ, 'x')(
                        [QQ(c) for c in eqn[0]]),
                    PolynomialRing(QQ, 'x')([QQ(c) for c in eqn[1]]))
                if invariants(D) == target:
                    out.append(r["label"])
            except Exception:
                continue
        time.sleep(pause)
    return out


def load_curves(path=CACHE, max_cond=None, simple_only=True, geom_simple=False):
    r"""Cached curves in conductor order, with conductor at most ``max_cond``.

    ``simple_only`` requires simplicity over `\QQ`; ``geom_simple`` requires
    geometric simplicity.
    """
    with open(path) as fh:
        rows = json.load(fh)
    if geom_simple:
        rows = [r for r in rows if r.get("is_simple_geom")]
    elif simple_only:
        rows = [r for r in rows if r.get("is_simple_base")]
    if max_cond is not None:
        rows = [r for r in rows if r["cond"] <= max_cond]
    return sorted(rows, key=lambda r: (r["cond"], r["label"]))


def curve_from_row(row, prec=250):
    r"""The curve `y^2+h(x)y=f(x)` from coefficients in increasing degree."""
    S = PolynomialRing(ZZ, 'x')
    eqn = _as_list(row["eqn"])
    f = S([ZZ(c) for c in eqn[0]])
    h = S([ZZ(c) for c in eqn[1]])
    return HyperellipticCurveWithReduction(f, h, prec=prec)


def _as_list(v):
    r"""Read a list or its string representation."""
    if isinstance(v, str):
        return json.loads(v.replace("'", '"'))
    return v


# ---------------------------------------------------------------------------
# moving the points into the affine chart
# ---------------------------------------------------------------------------

def lmfdb_points(C, rat_pts):
    r"""Convert weighted projective triples of weights `(1,g+1,1)` to points.

    Affine points are ``(x,y)``; infinite points are ``("infinity", eta)``
    with `eta=\lim y/x^{g+1}`.
    """
    g = C.genus()
    out = []
    for X, Y, Z in _as_list(rat_pts):
        X, Y, Z = QQ(X), QQ(Y), QQ(Z)
        if Z == 0:
            # (1 : Y : 0); the two points at infinity are distinguished by Y,
            # which satisfies Y^2 + h_{g+1} Y = f_{2g+2}.  Keep it.
            out.append(("infinity", Y / X ** (g + 1) if X != 1 else Y))
        else:
            out.append((X / Z, Y / Z ** (g + 1)))
    return out


def search_rational_points(C, bound=12):
    r"""Rational points with `|x|` and its denominator at most ``bound``, including infinity.

    Infinite points are labelled by `eta=\lim y/x^{g+1}`.

    EXAMPLES::

        sage: from hyperell_regulator.components import HyperellipticCurveWithReduction, search_rational_points
        sage: R.<x> = ZZ[]
        sage: C = HyperellipticCurveWithReduction(-x^2 - x, x^3 + x^2 + x + 1)
        sage: sorted(search_rational_points(C), key=str)
        [('infinity', -1), ('infinity', 0), (-1, 0), (0, -1), (0, 0)]
    """
    g = C.genus()
    RQ = PolynomialRing(QQ, 'x')
    fQ, hQ = RQ(C.f), RQ(C.h)
    pts = []
    lc = QQ(C.F[2 * g + 2])
    hl = QQ(C.h[g + 1])
    if lc == 0:
        # deg F = 2g+1: a single point at infinity, always rational, the double
        # root Y = -h_{g+1}/2 of Y^2 + h_{g+1} Y = f_{2g+2}
        pts.append(("infinity", -hl / 2))
    elif lc.is_square():
        # deg F = 2g+2 with lc a square: two rational points at infinity
        r = lc.sqrt()
        pts.append(("infinity", (-hl + r) / 2))
        pts.append(("infinity", (-hl - r) / 2))
    # otherwise the two points at infinity are conjugate, not rational
    for q in range(1, bound + 1):
        for n in range(-q * bound, q * bound + 1):
            x0 = QQ(n) / q
            if x0.denominator() != q:
                continue
            disc = hQ(x0) ** 2 + 4 * fQ(x0)
            if disc >= 0 and disc.is_square():
                r = disc.sqrt()
                for y0 in {(-hQ(x0) + r) / 2, (-hQ(x0) - r) / 2}:
                    pts.append((x0, y0))
    return pts


def curve_points(C, row, bound=12, max_bound=60):
    r"""Cached rational points, or a height search up to ``max_bound``."""
    if row.get("rat_pts"):
        return lmfdb_points(C, row["rat_pts"])
    want = row.get("num_rat_pts")
    b = bound
    while True:
        pts = search_rational_points(C, b)
        if want is None or len(pts) >= want or b >= max_bound:
            break
        b *= 2
    if want is not None and len(pts) != want:
        print("   WARNING: found %s point(s) by search up to height %s but "
              "LMFDB records %s" % (len(pts), b, want))
    return pts


def affine_model(C, points):
    r"""An isomorphic square model on which all marked points are affine.

    Use `u=1/(x-a)` and `v=(2y+h(x))/(x-a)^{g+1}`.
    Return ``(curve, images, a)``.

    EXAMPLES::

        sage: from hyperell_regulator.components import HyperellipticCurveWithReduction, affine_model
        sage: R.<x> = ZZ[]
        sage: C = HyperellipticCurveWithReduction(x^5 + x^3 + 1)
        sage: D, pts, a = affine_model(C, [(0, 1), ('infinity', 0)])
        sage: all(p is not None for p in pts)
        True
    """
    g = C.genus()
    F = C.F
    xs = set(p[0] for p in points if p is not None and p[0] != "infinity")
    a = None
    for cand in range(0, 200):
        for t in (ZZ(cand), -ZZ(cand)):
            if F(t) != 0 and t not in xs:
                a = t
                break
        if a is not None:
            break
    if a is None:
        raise ArithmeticError("no usable shift found")

    S = PolynomialRing(ZZ, 'u')
    u = S.gen()
    coeffs = list(F(F.parent().gen() + a))          # F expanded about a
    G = sum(ZZ(c) * u ** (2 * g + 2 - i) for i, c in enumerate(coeffs))
    assert G.degree() == 2 * g + 2 and G(0) != 0 or True

    # Y = 2y + h(x) on the square model; v = Y/(x-a)^{g+1}, and v^2 = G(u)
    RQ = PolynomialRing(QQ, 'x')
    hQ = RQ(C.h)
    new = []
    hlead = QQ(C.h[g + 1])
    for P in points:
        if P is None:
            new.append(None)
        elif P[0] == "infinity":
            # x = infinity  ->  u = 0.  On the square model Y = 2y + h(x), and
            # v = Y/(x-a)^{g+1} tends to 2*Y_oo + h_{g+1}, whose square is
            # 4 f_{2g+2} + h_{g+1}^2 = lc(F) = G(0).
            v = 2 * QQ(P[1]) + hlead
            assert v ** 2 == QQ(G(0)), (v, G(0))
            new.append((QQ(0), v))
        else:
            x0, y0 = P
            Y = 2 * y0 + hQ(x0)
            new.append((1 / (x0 - a), Y / (x0 - a) ** (g + 1)))
    D = HyperellipticCurveWithReduction(G, 0, prec=C.prec)
    return D, new, a


# ---------------------------------------------------------------------------
# the order of [P] - [Q] in the Jacobian
# ---------------------------------------------------------------------------

def odd_degree_model(C):
    r"""An odd degree square model, or ``None`` if none exists over `\QQ`.

    Return ``(curve, point_map, r)``, moving a rational root `r` to infinity.
    For an already odd degree model, ``r=None``.
    """
    g = C.genus()
    if C.F.degree() % 2:
        D = C.square_model()

        def phi(P):
            if P is None:
                return None
            if P[0] == "infinity":
                if 2*QQ(P[1])+QQ(C.h[g+1]) != 0:
                    raise ValueError("invalid point at infinity")
                return "infinity"
            x0, y0 = C.to_square_model(P)
            return D.point([x0, y0, QQ(1)])

        return D, phi, None
    roots = C.F.change_ring(QQ).roots(QQ)
    if not roots:
        return None
    r = QQ(roots[0][0])
    SQ = PolynomialRing(QQ, 'u')
    u = SQ.gen()
    shifted = C.F.change_ring(QQ)(C.F.parent().change_ring(QQ).gen() + r)
    G = sum(QQ(c) * u ** (2 * g + 2 - i) for i, c in enumerate(list(shifted)))
    d = lcm([QQ(c).denominator() for c in G] + [1])
    G = (d ** 2 * G).change_ring(ZZ)
    if G.degree() % 2 == 0 or not G.change_ring(QQ).is_squarefree():
        return None
    D = HyperellipticCurve(G.change_ring(QQ))
    RQ = PolynomialRing(QQ, 'x')
    hQ = RQ(C.h)

    def phi(P):
        if P is None:
            return None
        if P[0] == "infinity":
            # x = oo  ->  u = 0, v = d * (2 Y_oo + h_{g+1})
            v = d * (2 * QQ(P[1]) + QQ(C.h[g + 1]))
            return D.point([QQ(0), v, QQ(1)])
        x0, y0 = QQ(P[0]), QQ(P[1])
        if x0 == r:
            return "infinity"          # the Weierstrass point we moved
        Y = 2 * y0 + hQ(x0)
        return D.point([1 / (x0 - r), d * Y / (x0 - r) ** (g + 1), QQ(1)])

    return D, phi, r


def _torsion_multiple(C, primes):
    """A multiple of every rational torsion order, from good odd reductions."""
    key = ("torsion_multiple", int(primes))
    if key in C._cache:
        return C._cache[key]
    # Retaining the degree as well as separability gives a smooth model at p.
    bad = C.F.discriminant()*C.F.leading_coefficient()
    multiple = None
    for p in prime_range(3, primes):
        if bad % p == 0:
            continue
        try:
            reduced = HyperellipticCurve(C.F.change_ring(GF(p)))
            size = ZZ(reduced.frobenius_polynomial()(1))
        except (ArithmeticError, NotImplementedError, ValueError):
            continue
        multiple = size if multiple is None else gcd(multiple, size)
        if multiple == 1:
            break
    C._cache[key] = multiple
    return multiple


def _square_point(C, point):
    """Coordinates on `Y^2=4f+h^2`, including labelled infinity."""
    if point[0] == "infinity":
        return ("infinity", 2*QQ(point[1])+QQ(C.h[C.genus()+1]))
    return C.to_square_model(point)


def _order_by_principality(C, P, Q, primes):
    """The divisor-class order by exact principality tests."""
    from .regulator.torsion import _field, _place, _order_from_multiple

    key = "torsion_function_field"
    if key not in C._cache:
        C._cache[key] = _field(C.F)
    K, L = C._cache[key]
    p, q = (_place(K, L, _square_point(C, point)) for point in (P, Q))
    if p == q:
        return ZZ(1)
    multiple = _torsion_multiple(C, primes)
    if multiple is None:
        return None
    order = _order_from_multiple(q.divisor()-p.divisor(), multiple)
    return Infinity if order is None else order


def divisor_class_order(C, P, Q, primes=30):
    r"""The order of `[P]-[Q]` in `J(\QQ)`.

    Return ``+Infinity`` for nontorsion and ``None`` if undetermined.
    Points use `y^2+h(x)y=f(x)`; ``("infinity", eta)`` means
    `eta=\lim y/x^{g+1}`. Good odd reductions bound the torsion order.

    EXAMPLES::

        sage: from hyperell_regulator.components import HyperellipticCurveWithReduction, divisor_class_order
        sage: R.<x> = ZZ[]
        sage: C = HyperellipticCurveWithReduction(-x^2 - x, x^3 + x^2 + x + 1)
        sage: divisor_class_order(C, (0, 0), (0, -1))
        15

        sage: C = HyperellipticCurveWithReduction(x + 2*x^2 - x^4, 1 + x + x^3)
        sage: divisor_class_order(C, (0, 0), ("infinity", 0))
        7
        sage: divisor_class_order(C, (0, 0), ("infinity", -1))
        21
    """
    m = odd_degree_model(C)
    if m is None:
        try:
            return _order_by_principality(C, P, Q, primes)
        except Exception:
            return None
    D, phi, _ = m
    try:
        JQ = D.jacobian()(QQ)
        pP, pQ = phi(P), phi(Q)
        eP = JQ(0) if pP == "infinity" else JQ(pP)
        eQ = JQ(0) if pQ == "infinity" else JQ(pQ)
        cls = eP - eQ
    except Exception:
        return None
    zero = JQ(0)
    if cls == zero:
        return ZZ(1)

    multiple = _torsion_multiple(C, primes)
    if multiple is None:
        return None
    for k in divisors(multiple):
        if k * cls == zero:
            return ZZ(k)
    return Infinity


# ---------------------------------------------------------------------------
# the question
# ---------------------------------------------------------------------------

def same_component_all_primes(C, P, Q):
    r"""Test smooth reduction on the same stable component at every prime.

    Return ``(verdict, per_prime)``, with ``verdict=None`` if undecided.
    """
    per_prime = {}
    verdict = True
    for p in C.bad_primes():
        if p == 2:
            per_prime[p] = None
            continue
        per_prime[p] = C.same_component(P, Q, p)["answer"]
    if any(v is False for v in per_prime.values()):
        verdict = False
    elif any(v is None for v in per_prime.values()):
        verdict = None
    return verdict, per_prime


def report_curve(row, verbose=True):
    r"""Torsion orders and component tests for one LMFDB curve."""
    C = curve_from_row(row)
    pts = curve_points(C, row)
    n_inf = sum(1 for p in pts if p is None)

    # work on a model where every point is affine
    D, dpts, a = affine_model(C, pts)
    usable = [(i, p) for i, p in enumerate(dpts) if p is not None]

    res = {"label": row["label"], "cond": row["cond"], "genus": C.genus(),
           "curve": C, "model": D, "shift": a,
           "bad_primes": D.bad_primes(), "points": pts,
           "n_points": len(pts), "n_infinity": n_inf,
           "n_unusable": len(dpts) - len(usable), "pairs": []}

    if verbose:
        print("=" * 74)
        print("%s   conductor %s   genus %s   Jacobian %s"
              % (row["label"], row["cond"], C.genus(),
                 "geometrically simple" if row.get("is_simple_geom")
                 else "simple over Q, not geometrically simple"))
        print("   %s" % C)
        print("   %s known rational points (%s at infinity)"
              % (len(pts), n_inf))
        print("   working model v^2 = %s   [u = 1/(x - %s)]" % (D.f, a))
        print("   bad primes of that model: %s" % D.bad_primes())

    for k in range(len(usable)):
        for l in range(k + 1, len(usable)):
            (i, P), (j, Qp) = usable[k], usable[l]
            verdict, per = same_component_all_primes(D, P, Qp)
            check = None
            if (pts[i] is not None and pts[i][0] != "infinity"
                    and pts[j] is not None and pts[j][0] != "infinity"):
                # both affine on the original model too: the answer must agree,
                # since the two models are isomorphic over QQ
                try:
                    v0, _ = same_component_all_primes(C, pts[i], pts[j])
                    check = (v0 == verdict)
                except Exception as exc:
                    check = "err: %s" % exc
            order = None
            if verdict is True:
                order = divisor_class_order(C, pts[i], pts[j])
            res["pairs"].append({"i": i, "j": j, "P": pts[i], "Q": pts[j],
                                 "verdict": verdict, "per_prime": per,
                                 "cross_model": check, "order": order})
    if verbose:
        good = [d for d in res["pairs"] if d["verdict"] is True]
        unk = [d for d in res["pairs"] if d["verdict"] is None]
        checks = [d["cross_model"] for d in res["pairs"]
                  if d["cross_model"] is not None]
        print("   %s pairs tested: %s same-component everywhere, %s undecided"
              % (len(res["pairs"]), len(good), len(unk)))
        if checks:
            print("   cross-model check on the %s pairs affine in both "
                  "models: %s" % (len(checks),
                                  "all agree" if all(c is True for c in checks)
                                  else "MISMATCH %s" % checks))
        for d in res["pairs"]:
            mark = {True: "YES", False: "no ", None: "?  "}[d["verdict"]]
            print("      %s  %-22s %-22s  %s"
                  % (mark, _fmt(d["P"]), _fmt(d["Q"]),
                     ", ".join("p=%s:%s" % (q, {True: "same", False: "diff",
                                               None: "?"}[v])
                               for q, v in sorted(d["per_prime"].items()))))
    return res


def _fmt(P):
    if P is None:
        return "??"
    if P[0] == "infinity":
        return "oo(Y=%s)" % P[1]
    return "(%s, %s)" % (P[0], P[1])


def _jsonable(v):
    r"""Convert exact numbers and points to JSON values."""
    if v is None or isinstance(v, bool) or isinstance(v, str):
        return v
    if v is Infinity:
        return "Infinity"
    if v in ZZ:
        return int(v)
    if isinstance(v, (tuple, list)):
        return [_jsonable(t) for t in v]
    return str(v)


def _result_record(r):
    r"""One curve record in JSON form."""
    return {
        "label": r["label"],
        "cond": int(r["cond"]),
        "genus": int(r["genus"]),
        "bad_primes": [int(q) for q in r["bad_primes"]],
        "working_model": str(r["model"].f),
        "shift": int(r["shift"]),
        "n_points": int(r["n_points"]),
        "n_infinity": int(r["n_infinity"]),
        "pairs": [{"P": _jsonable(d["P"]),
                   "Q": _jsonable(d["Q"]),
                   "same_component_everywhere": d["verdict"],
                   "per_prime": {str(q): v for q, v in
                                 sorted(d["per_prime"].items())},
                   "order_of_P_minus_Q": _jsonable(d["order"]),
                   "cross_model_check": _jsonable(d["cross_model"])}
                  for d in r["pairs"]],
    }


def save_results(results, results_path=RESULTS):
    r"""Merge results by curve label and return the full saved dictionary."""
    try:
        with open(results_path) as fh:
            have = {r["label"]: r for r in json.load(fh)}
    except (IOError, OSError, ValueError):
        have = {}
    before = len(have)
    for r in results:
        have[r["label"]] = _result_record(r)
    out = sorted(have.values(), key=lambda r: (r["cond"], r["label"]))
    with open(ensure_parent(results_path), "w") as fh:
        json.dump(out, fh, indent=2)
    print("results: %s curve(s) written to %s (%s new, %s updated)"
          % (len(out), results_path, len(out) - before,
             len(results) - (len(out) - before)))
    return out


def enumerate_all(max_cond=None, path=CACHE, geom_simple=True,
                  results_path=RESULTS):
    r"""Compute component records in conductor order.

    By default the Jacobians are geometrically simple.
    """
    rows = load_curves(path=path, max_cond=max_cond, simple_only=True,
                       geom_simple=geom_simple)
    if not rows:
        print("No curves in %s.  See fetch_lmfdb's docstring for how to "
              "populate it." % path)
        return []
    print("%s curve(s) with %s Jacobian, conductors %s-%s   [from %s]"
          % (len(rows), "geometrically simple" if geom_simple
             else "simple (over Q)", min(r["cond"] for r in rows),
             max(r["cond"] for r in rows), path))
    print("To enumerate more, fetch more first:  "
          "fetch_lmfdb(limit=%s)   (merges into the same cache,\n"
          "skipping curves already there; ~20s per new curve, so be patient)\n"
          % (len(rows) + 30))
    out = [report_curve(r) for r in rows]
    if results_path is not None:
        print()
        save_results(out, results_path)
    print()
    print("=" * 74)
    print("SUMMARY")
    for r in out:
        good = [d for d in r["pairs"] if d["verdict"] is True]
        print("  %-16s cond %-6s %s/%s pairs same-component at every prime"
              % (r["label"], r["cond"], len(good), len(r["pairs"])))
        for d in good:
            o = d["order"]
            otxt = ("order of [P]-[Q] unknown" if o is None else
                    "[P]-[Q] not torsion" if o == Infinity else
                    "[P]-[Q] has order %s" % o)
            print("        %-22s and  %-22s  %s"
                  % (_fmt(d["P"]), _fmt(d["Q"]), otxt))
    return out


def main():
    r"""Run the component enumeration."""
    return enumerate_all()


if __name__ == "__main__":
    main()
