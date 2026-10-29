"""Embedded semilunar paths shared by regulator cuts and period cycles.

Geometry is chosen in float64. Labels retain the indices of the original
points, so arbitrary-precision integrations never round their endpoints.
"""

import numpy as np


def angular_order(points, start, eps=1e-13):
    """Increasing argument about a perturbed centroid, starting at ``start``."""
    v = np.asarray(points, complex)
    c = v.mean() + eps * (np.pi + np.e * 1j)
    order = sorted(range(len(v)), key=lambda i: np.angle(v[i] - c))
    k = order.index(start)
    return order[k:] + order[:k], c, k


def semilunar_path(points, start, eps=1e-13, radius=1.7, arc=24,
                   fan=64, include_ray=True):
    """Return validated nodes, original-point labels, and an optional ray.

    Follow the angular order, rotating to the requested first point. If the
    reintroduced closing edge crosses another point, route it outside the
    configuration. Ordinary detour corners have label ``None``.
    """
    v = np.asarray(points, complex)
    order, c, k = angular_order(v, start, eps)
    nodes, labels = list(v[order]), list(order)

    def validate():
        if include_ray:
            return _pick_ray(nodes, v, c, fan)
        try:
            check_path_embedded(nodes, crit=v, include_ray=False)
            return True
        except NotImplementedError:
            return None

    direction = validate()
    if direction is None:
        if k == 0:
            raise NotImplementedError("no embedded semilunar path")
        j = len(v) - 1 - k
        R = radius * max(np.abs(v - c))
        bend = _detour(nodes[j], nodes[j + 1], c, R, arc)
        nodes[j + 1:j + 1] = bend
        labels[j + 1:j + 1] = [None] * len(bend)
        direction = validate()
        if direction is None:
            raise NotImplementedError("no embedded semilunar path even with a detour")
    return np.array(nodes), labels, direction if include_ray else None


def _pick_ray(nodes, crit, c, fan):
    a = np.asarray(nodes, complex)
    last = a[-1]
    scale = max(np.abs(a - c).max(), 1.0)
    if abs(last - c) > 0:
        d = (last - c) / abs(last - c)
        try:
            check_path_embedded(a, d, crit)
            return d
        except NotImplementedError:
            pass
    best, room = None, -1.0
    for d in [np.exp(2j * np.pi * m / fan) for m in range(fan)]:
        try:
            check_path_embedded(a, d, crit)
        except NotImplementedError:
            continue
        probe = last + d * np.linspace(0.02, 6.0, 60) * scale
        gap = min(_dist_to_path(p, a[:-1]) for p in probe)
        if gap > room:
            best, room = d, gap
    return best


def _dist_to_path(p, nodes):
    best = np.inf
    for u, v in zip(nodes[:-1], nodes[1:]):
        t = np.clip(((p - u) * np.conj(v - u)).real / abs(v - u) ** 2, 0., 1.)
        best = min(best, abs(u + t * (v - u) - p))
    return best


def _detour(a, b, c, R, n):
    ta, tb = np.angle(a - c), np.angle(b - c)
    sweep = np.mod(tb - ta, 2 * np.pi)
    return list(c + R * np.exp(1j * (ta + sweep * np.linspace(0., 1., n))))


def check_path_embedded(e, ray=None, crit=None, tol=1e-9, include_ray=True):
    """Reject crossings, overlaps, or paths through other marked points.

    By default also validate a terminal ray, as required by regulator cuts.
    Period paths use ``include_ray=False`` and check all finite branch points.

    EXAMPLES::

        sage: from hyperell_regulator.regulator.cover import check_path_embedded
        sage: check_path_embedded([-2., -1., 0., 1., 2.])
        sage: check_path_embedded([0., 1., 2., -2., -1.])
        Traceback (most recent call last):
        ...
        NotImplementedError: ...
    """
    e = np.asarray(e, complex)
    scale = max(np.abs(e).max(), 1.0)
    pts = list(e)
    if include_ray:
        if ray is None:
            ray = e[-1] - e.mean()
            if abs(ray) < tol * scale:
                raise NotImplementedError("the last critical value sits at the centroid")
        if not abs(ray):
            raise NotImplementedError("the terminal ray has zero direction")
        pts.append(e[-1] + ray / abs(ray) * 10.0 * scale)
    if crit is None:
        crit = e
    crit = np.asarray(crit, complex)
    segs = list(zip(pts[:-1], pts[1:]))

    def side(a, b, c):
        return ((b - a).conjugate() * (c - a)).imag

    for i, (a, b) in enumerate(segs):
        if a == b:
            raise NotImplementedError("the path has a zero-length segment")
        for c in crit:
            if abs(c - a) < tol * scale or abs(c - b) < tol * scale:
                continue
            t = ((c - a).conjugate() * (b - a)).real / abs(b - a) ** 2
            if -tol < t < 1 + tol and abs(a + t * (b - a) - c) < tol * scale:
                raise NotImplementedError(
                    "the critical value %r lies on the segment %r--%r" % (c, a, b))
        for c, d in segs[i + 2:]:
            if (side(a, b, c) * side(a, b, d) < 0
                    and side(c, d, a) * side(c, d, b) < 0):
                raise NotImplementedError(
                    "the segments %r--%r and %r--%r cross" % (a, b, c, d))
            if (abs(side(a, b, c)) < tol * scale
                    and abs(side(a, b, d)) < tol * scale):
                u = (b - a) / abs(b - a)
                t0, t1 = ((c - a) / u).real, ((d - a) / u).real
                if min(max(t0, t1), abs(b - a)) - max(min(t0, t1), 0.) > tol * scale:
                    raise NotImplementedError(
                        "the segments %r--%r and %r--%r overlap" % (a, b, c, d))
