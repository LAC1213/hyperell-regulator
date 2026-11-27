r"""LMFDB curves, rational points and exact torsion orders for genus 2."""

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


# ---------------------------------------------------------------------------
# data
# ---------------------------------------------------------------------------

def fetch_lmfdb(max_cond=None, path=CACHE, pause=6, limit=None,
                geom_simple=True, email="acampo@ihes.fr", end_q_z=False, end_algs=None,
                on_curve=None, _cancel_event=None):
    r"""Fetch and merge LMFDB curve and rational-point data into ``path``.

    Requests are paginated with ``pause`` seconds between them; the scan
    offset is saved in ``_scan_state.json``. ``end_q_z=True`` selects
    End_Q(J)=Z and ``end_algs`` selects LMFDB algebra labels. ``on_curve``
    is called with each curve as soon as it is cached.
    """
    import time
    import urllib.request
    import urllib.parse

    headers = {"User-Agent": "sage-research-script/1.0 (%s)" % email}
    state_path = os.path.join(os.path.dirname(path) or ".",
                              "_scan_state.json")
    calls = [0]

    def check_cancelled():
        if _cancel_event is not None and _cancel_event.is_set():
            raise InterruptedError('LMFDB fetch cancelled')

    def sleep(seconds):
        if _cancel_event is None:
            time.sleep(seconds)
        elif _cancel_event.wait(seconds):
            check_cancelled()

    def api(table, **kw):
        kw.setdefault("_format", "json")
        url = "https://www.lmfdb.org/api/%s/?%s" % (
            table, urllib.parse.urlencode(kw))
        for attempt in range(4):
            check_cancelled()
            calls[0] += 1
            raw = urllib.request.urlopen(
                urllib.request.Request(url, headers=headers), timeout=120).read()
            check_cancelled()
            if raw.lstrip()[:1] != b"<":
                return json.loads(raw)
            wait = 90 * (attempt + 1)
            print("   (LMFDB served a CAPTCHA; waiting %ss)" % wait)
            sleep(wait)
        raise RuntimeError(
            "LMFDB kept returning HTML (a CAPTCHA challenge) instead of JSON")

    try:
        with open(path) as fh:
            have = {r["label"]: r for r in json.load(fh)}
    except (IOError, OSError, ValueError):
        have = {}

    cached_labels, emitted, pts = set(have), set(), {}

    def merge(r):
        previous = have.get(r['label'], {})
        q = pts.get(r['label'], {})
        have[r['label']] = {
            **previous, 'label': r['label'], 'cond': int(r['cond']),
            'is_simple_base': bool(r['is_simple_base']),
            'is_simple_geom': bool(r['is_simple_geom']),
            'end_alg': r.get('end_alg'), 'eqn': r['eqn'],
            'num_rat_pts': q.get('num_rat_pts', previous.get('num_rat_pts')),
            'rat_pts': q.get('rat_pts', previous.get('rat_pts')) or [],
            'rat_pts_verified': q.get('rat_pts_v', previous.get('rat_pts_verified')),
            'mw_invs': q.get('mw_invs', previous.get('mw_invs')),
        }

    def save():
        # Readers in the analysis thread must always see complete JSON.
        import tempfile
        from pathlib import Path
        check_cancelled()
        destination = Path(ensure_parent(path))
        out = sorted(have.values(), key=lambda r: (r['cond'], r['label']))
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(mode='w', dir=destination.parent,
                                             prefix=destination.name+'.', delete=False) as fh:
                temporary = fh.name
                json.dump(out, fh, indent=2)
            os.replace(temporary, destination)
        finally:
            if temporary is not None and os.path.exists(temporary):
                os.unlink(temporary)
        return out

    def publish(rows):
        if on_curve is None or not rows:
            return
        for r in rows:
            merge(r)
        save()
        chosen = {r['label'] for r in selected}
        for r in rows:
            label = r['label']
            if label in chosen and label not in emitted:
                check_cancelled()
                on_curve(dict(have[label]))
                emitted.add(label)

    # --- 1. the curves, filtered server-side ------------------------------
    if end_algs is None and end_q_z:
        end_algs = ('Q',)
    # The supported rational endomorphism fields imply Q-simplicity,
    # but do not imply geometric simplicity.
    # Include endomorphism metadata for old cached rows that are now excluded.
    flag = "is_simple_geom" if geom_simple and end_algs is None else "is_simple_base"
    selected, catalog, offset = [], [], 0
    while True:
        page = api("g2c_curves", _offset=offset, _limit=100, _sort="cond",
                   _fields="label,cond,is_simple_base,is_simple_geom,eqn,end_alg",
                   **{flag: "true"})["data"]
        if not page:
            break
        bounded = [r for r in page if max_cond is None or r["cond"] <= max_cond]
        catalog.extend(bounded)
        selected.extend(r for r in bounded if end_algs is None or r.get('end_alg') in end_algs)
        if limit is not None:
            selected = selected[:limit]
        # Existing point data is already usable as soon as metadata refreshes.
        # New curves wait for their rational-point data before publication.
        publish([r for r in bounded if r['label'] in cached_labels])
        print("   g2c_curves offset %-5s -> %s selected" % (offset,
                                                            len(selected)))
        if limit is not None and len(selected) >= limit:
            break
        if max_cond is not None and max(r["cond"] for r in page) > max_cond:
            break
        offset += len(page)
        sleep(pause)

    # --- 2. the rational points, resuming the scan ------------------------
    want = set(r["label"] for r in selected) - cached_labels
    by_label = {r['label']: r for r in selected}
    if want:
        try:
            with open(state_path) as fh:
                offset = int(json.load(fh).get("ratpts_offset", 0))
        except (IOError, OSError, ValueError):
            offset = 0
        start_offset = offset
        while want - set(pts):
            sleep(pause)
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
                here = int(page[0]["label"].split(".")[0])
                lowest = min(int(t.split(".")[0]) for t in want)
                if here > lowest:
                    print("   g2c_ratpts resume point overshot "
                          "(offset %s is at conductor %s, need %s); "
                          "restarting scan" % (offset, here, lowest))
                    offset, start_offset = 0, 0
                    continue
            for r in page:
                if r["label"] in want:
                    pts[r["label"]] = r
            publish([by_label[r['label']] for r in page if r['label'] in want])
            print("   g2c_ratpts offset %-5s -> %s/%s matched"
                  % (offset, len(pts), len(want)))
            offset += len(page)
        check_cancelled()
        with open(ensure_parent(state_path), "w") as fh:
            # int(): Sage Integers are not JSON serialisable
            # Integers, which json cannot serialise
            json.dump({"ratpts_offset": int(offset)}, fh)

    # --- 3. merge into the cache -----------------------------------------
    chosen = {r['label'] for r in selected}
    for r in catalog:
        if r['label'] not in chosen and r['label'] not in have:
            continue
        merge(r)
    out = save()
    if on_curve is not None:
        for r in sorted(selected, key=lambda r: (r['cond'], r['label'])):
            if r['label'] not in emitted:
                check_cancelled()
                on_curve(dict(have[r['label']]))
                emitted.add(r['label'])
    span = ', conductors %s-%s' % (out[0]['cond'], out[-1]['cond']) if out else ''
    print("cache %s now holds %s curve(s)%s   [%s API call(s)]"
          % (path, len(out), span, calls[0]), flush=True)
    return out


def load_curves(path=CACHE, max_cond=None, simple_only=True, geom_simple=False, end_q_z=False,
                end_algs=None):
    r"""Cached curves in conductor order, with conductor at most ``max_cond``.

    ``simple_only`` requires simplicity over `\QQ`; ``geom_simple`` requires
    geometric simplicity. ``end_q_z`` selects End_Q(J)=Z independently of
    those two flags, using the cached ``end_alg`` field.
    ``end_algs`` selects a specified collection of rational algebra labels.
    """
    with open(path) as fh:
        rows = json.load(fh)
    if end_algs is not None:
        rows = [r for r in rows if r.get('end_alg') in end_algs]
    elif end_q_z:
        rows = [r for r in rows if end_q_is_z(r)]
    elif geom_simple:
        rows = [r for r in rows if r.get("is_simple_geom")]
    elif simple_only:
        rows = [r for r in rows if r.get("is_simple_base")]
    if max_cond is not None:
        rows = [r for r in rows if r["cond"] <= max_cond]
    return sorted(rows, key=lambda r: (r["cond"], r["label"]))


def end_q_is_z(row):
    """Whether the cached rational endomorphism algebra proves End_Q(J)=Z.

    LMFDB's ``end_alg='Q'`` specifies End_Q(J) tensor Q = Q. The only
    integral order in Q containing 1 is Z. Missing metadata is not evidence.
    """
    return row.get('end_alg') == 'Q'


def curve_from_row(row, prec=250):
    r"""The curve `y^2+h(x)y=f(x)` from coefficients in increasing degree."""
    S = PolynomialRing(ZZ, 'x')
    eqn = _as_list(row["eqn"])
    f = S([ZZ(c) for c in eqn[0]])
    h = S([ZZ(c) for c in eqn[1]])
    C = HyperellipticCurveWithReduction(f, h, prec=prec)
    C.lmfdb_record = dict(row)
    C.lmfdb_label = row.get('label')
    return C


def load_curve(label, *, cache=None, fetch=True, refresh=False, prec=250):
    r"""Load one LMFDB genus 2 curve over QQ with reduction methods.

    Uncached labels are fetched individually and merged into the cache.
    ``fetch=False`` requires cached data; ``refresh=True`` refetches. The
    curve carries ``lmfdb_label`` and ``lmfdb_record``.
    """
    import re
    import tempfile
    from pathlib import Path
    from urllib.parse import urlencode
    from urllib.request import Request, urlopen

    if not isinstance(label, str) or not re.fullmatch(r'[1-9][0-9]*\.[a-z]+\.[1-9][0-9]*\.[1-9][0-9]*', label):
        raise ValueError('expected an LMFDB genus 2 curve label, e.g. 394.a.394.1')
    destination = Path(cache if cache is not None else data_path('lmfdb_cache', 'g2c.json'))
    try:
        rows = json.loads(destination.read_text())
    except FileNotFoundError:
        rows = []
    have = {r['label']: r for r in rows}
    if label in have and not refresh:
        return curve_from_row(have[label], prec=prec)
    if not fetch:
        raise ValueError('%s is not available without fetching from LMFDB' % label)

    def query(table, fields):
        params = dict(label='s'+label, _format='json', _fields=fields, _limit=2)
        request = Request('https://www.lmfdb.org/api/%s/?%s' % (table, urlencode(params)),
                          headers={'User-Agent': 'hyperell-regulator/0.1 (genus-two research)'})
        with urlopen(request, timeout=120) as response:
            raw = response.read()
        if raw.lstrip().startswith(b'<'):
            raise RuntimeError('LMFDB returned HTML instead of JSON (possibly a CAPTCHA)')
        records = json.loads(raw)['data']
        if len(records) > 1 or any(r.get('label') != label for r in records):
            raise ValueError('LMFDB returned an unexpected curve for label '+label)
        return records

    curves = query('g2c_curves', 'label,cond,eqn,end_alg,is_simple_base,is_simple_geom')
    if not curves:
        raise ValueError('LMFDB genus 2 curve not found: '+label)
    row = dict(have.get(label, {}), **curves[0])
    points = query('g2c_ratpts', 'label,num_rat_pts,rat_pts,rat_pts_v,mw_invs')
    if points:
        for key in ('num_rat_pts', 'rat_pts', 'mw_invs'):
            row[key] = points[0].get(key)
        row['rat_pts_verified'] = points[0].get('rat_pts_v')
    C = curve_from_row(row, prec=prec)
    if C.genus() != 2 or int(row['cond']) != int(label.split('.')[0]):
        raise ValueError('LMFDB returned inconsistent genus 2 curve data for '+label)
    have[label] = row
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', dir=destination.parent, delete=False) as stream:
            temporary = Path(stream.name)
            json.dump(sorted(have.values(), key=lambda r: (int(r['cond']), r['label'])), stream, indent=2)
            stream.write('\n')
        temporary.replace(destination)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    return C


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
# the order of P - Q in the Jacobian
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
    r"""The order of `P-Q` in `J(\QQ)`.

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
