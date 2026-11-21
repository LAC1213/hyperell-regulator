r"""Regulator integrals by adaptive quadrature over `\CC`.

A smooth partition of unity separates clusters of nearby singularities.
"""

import scipy.integrate as spyint
import numpy as np
from scipy.special import psi

# Integrate close branch points in local boxes with a C^2 cutoff.
# The complementary weight covers the rest of the plane. Their sum is one,
# so quadrature sees no discontinuity at box boundaries.


def _find_clusters(points, threshold):
    """Index groups of at least two points joined by distances below ``threshold``."""
    n = len(points)
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        pi, pj = find(i), find(j)
        if pi != pj:
            parent[pi] = pj

    for i in range(n):
        for j in range(i + 1, n):
            if abs(points[i] - points[j]) < threshold:
                union(i, j)

    groups = {}
    for i in range(n):
        groups.setdefault(find(i), []).append(i)
    return [g for g in groups.values() if len(g) > 1]


def _smoothstep(u):
    """The quintic step from 0 to 1 with vanishing first and second endpoint derivatives."""
    u = np.clip(u, 0.0, 1.0)
    return 6 * u**5 - 15 * u**4 + 10 * u**3


def _make_boxes(points, cluster_threshold=0.4, max_half_side=0.6, buffer_factor=1.3):
    """Disjoint square neighborhoods of close singularities.

    Return ``(center, R, R_buf, x_breaks, y_breaks)`` for each cluster.
    ``R`` and ``R_buf`` are core and buffer half-sides.
    """
    pts_arr = np.array(points, dtype=complex)
    clusters = _find_clusters(pts_arr, cluster_threshold)
    boxes = []
    for idxs in clusters:
        pts = pts_arr[idxs]
        centroid = np.mean(pts)
        in_r = max(np.abs(pts - centroid)) * 3 + 0.05
        others = np.array([pts_arr[k] for k in range(len(pts_arr)) if k not in idxs])
        # leave room for the buffer_factor expansion when checking clearance
        clearance = (0.6 / buffer_factor) * min(np.abs(others - centroid)) if len(others) else max_half_side
        R = min(in_r, clearance, max_half_side)
        x_breaks = sorted(set(np.round((pts_arr[i] - centroid).real, 12) for i in idxs))
        y_breaks = sorted(set(np.round((pts_arr[i] - centroid).imag, 12) for i in idxs))
        boxes.append([centroid, R, buffer_factor * R, x_breaks, y_breaks])

    # safety pass: shrink any pair of boxes whose *buffers* would overlap
    for i in range(len(boxes)):
        for j in range(i + 1, len(boxes)):
            ci, Rbi = boxes[i][0], boxes[i][2]
            cj, Rbj = boxes[j][0], boxes[j][2]
            d = abs(ci - cj)
            if Rbi + Rbj > 0.9 * d:
                scale = 0.9 * d / (Rbi + Rbj)
                boxes[i][1] *= scale; boxes[i][2] *= scale
                boxes[j][1] *= scale; boxes[j][2] *= scale

    return [tuple(b) for b in boxes]


def _weight(z, box):
    """A `C^2` cutoff equal to 1 on the core square and 0 outside the buffer."""
    centroid, R, R_buf, _, _ = box
    d = z - centroid
    t = max(abs(d.real), abs(d.imag)) / R
    if t <= 1.0:
        return 1.0
    tb = R_buf / R
    if t >= tb:
        return 0.0
    return 1.0 - _smoothstep((t - 1.0) / (tb - 1.0))


def _total_weight(z, boxes):
    # boxes are kept disjoint (including buffers), so a simple sum
    # is safe and stays within [0, 1].
    return sum(_weight(z, b) for b in boxes) if boxes else 0.0


def _integrate_plane(f_corrected, boxes, local_opts=None, global_opts=None):
    r"""Integrate the corrected function over `\CC` using local and polar quadrature.

    Return ``(real_part, imag_part, summed_error_estimate)``.
    """

    if local_opts is None:
        local_opts = {'limit': 150, 'epsabs': 1e-10, 'epsrel': 1e-10}
    if global_opts is None:
        global_opts = {'limit': 150, 'epsabs': 1e-9, 'epsrel': 1e-9}

    total_re = 0.0
    total_im = 0.0
    total_err = 0.0

    # ---- local pieces: one smoothly-weighted square per cluster ----
    for box in boxes:
        centroid, R, R_buf, xb, yb = box
        g = lambda x, y: _weight(centroid + x + 1j * y, box) * f_corrected(centroid + x + 1j * y)
        g_re = lambda x, y: np.real(g(x, y))
        g_im = lambda x, y: np.imag(g(x, y))
        xpts = sorted(set([0.0, -R, R] + list(xb)))
        ypts = sorted(set([0.0, -R, R] + list(yb)))
        opts = [dict(local_opts, points=xpts), dict(local_opts, points=ypts)]
        vre, ere = spyint.nquad(g_re, [[-R_buf, R_buf], [-R_buf, R_buf]], opts=opts)
        vim, eim = spyint.nquad(g_im, [[-R_buf, R_buf], [-R_buf, R_buf]], opts=opts)
        total_re += vre
        total_im += vim
        total_err += ere + eim

    # ---- global piece: whole plane, weighted by (1 - total local weight) ----
    def weighted(z):
        w = _total_weight(z, boxes)
        return (1.0 - w) * f_corrected(z) if w < 1.0 else 0.0

    r = lambda s: s / (1 - s)
    integrand = lambda s, theta: weighted(r(s) * np.exp(theta * 1j)) * r(s) / (1 - s) ** 2
    integrand_re = lambda s, theta: np.real(integrand(s, theta))
    integrand_im = lambda s, theta: np.imag(integrand(s, theta))

    # Breakpoints are no longer load-bearing for correctness (the
    # integrand has no jump left to find), but still help nquad by
    # flagging where the weight varies fastest, along the real axis.
    s_breaks = set()
    for centroid, R, R_buf, xb, yb in boxes:
        if abs(centroid.imag) < 1e-9:
            for rr in (centroid.real - R_buf, centroid.real + R_buf):
                rr_abs = abs(rr)
                s_breaks.add(rr_abs / (1 + rr_abs))
    s_breaks = sorted(s_breaks)

    opts = [dict(global_opts, points=s_breaks),
            dict(global_opts, points=[0.0, np.pi, 2 * np.pi])]
    vre, ere = spyint.nquad(integrand_re, [[0, 1], [0, 2 * np.pi]], opts=opts)
    vim, eim = spyint.nquad(integrand_im, [[0, 1], [0, 2 * np.pi]], opts=opts)
    total_re += vre
    total_im += vim
    total_err += ere + eim

    return total_re, total_im, total_err


def _coeffs(f):
    r"""The coefficients of ``f`` in descending degree."""
    if hasattr(f, "parent") and hasattr(f, "list"):
        return [float(c) for c in reversed(f.list())] or [0.0]
    return list(f)


# a : a rational number

# Compute the integral G(z) log|z-a| /|f(z)| over the complex plane when f(a) = 0
def regulator_integral_0(f, a=1, G=lambda z: 1, tol=1e-12, cluster_threshold=0.4):
    f = _coeffs(f)
    deg = len(f) - 1
    assert deg == 5 or deg == 6
    assert np.polyval(f, a) == 0

    roots = np.roots(f)
    roots1 = [r for r in roots if abs(r - a) > tol]
    assert len(roots1) == len(roots) - 1
    fder = np.polyder(f)

    const_1 = [G(roots1[i]) * np.log(np.abs(roots1[i] - a)) / np.abs(np.polyval(fder, roots1[i]))
               for i in range(deg - 1)]
    const_2 = G(a) / np.abs(np.polyval(fder, a))

    bump_1 = lambda z, z0: np.exp(-np.abs(z - z0) ** 2) / np.abs(z - z0)
    bump_2 = lambda z, z0: np.log(np.abs(z - z0)) * np.exp(-np.abs(z - z0) ** 2) / np.abs(z - z0)

    correction_1 = lambda z: sum([const_1[i] * bump_1(z, roots1[i]) for i in range(deg - 1)])
    correction_2 = lambda z: const_2 * bump_2(z, a)

    f_corrected = lambda z: (G(z) * np.log(np.abs(z - a)) / np.abs(np.polyval(f, z))
                              - correction_1(z) - correction_2(z))

    singular_points = list(roots1) + [a]
    boxes = _make_boxes(singular_points, cluster_threshold=cluster_threshold)

    int_re, int_im, err = _integrate_plane(f_corrected, boxes)

    result = (int_re + 1j * int_im
              + sum([np.pi ** 1.5 * const_1[i] for i in range(deg - 1)])
              + (-1 / 2 * np.pi ** 1.5 * (np.log(4) - psi(1))) * const_2)

    return result, err, 0.0
# Note:
# int_\C 1/|z| e^{-|z|^2} dxdy = Pi^{3/2}
# = np.pi**1.5
# int_\C log|z|/|z| e^{-|z|^2} dxdy = -1/2*Pi^{3/2}*(log(4)-psi(1))
# = -1/2 * np.pi**1.5 * (np.log(4) - psi(1))


# Compute the integral G(z) log|z-a| /|f(z)| over the complex plane when f(a) != 0
def regulator_integral_1(f, a=1, G=lambda z: 1, cluster_threshold=0.4):
    f = _coeffs(f)
    deg = len(f) - 1
    assert deg == 5 or deg == 6
    assert np.polyval(f, a) != 0

    roots = np.roots(f)
    fder = np.polyder(f)

    const_1 = [G(roots[i]) * np.log(np.abs(roots[i] - a)) / np.abs(np.polyval(fder, roots[i]))
               for i in range(deg)]

    bump_1 = lambda z, z0: np.exp(-np.abs(z - z0) ** 2) / np.abs(z - z0)
    correction_1 = lambda z: sum([const_1[i] * bump_1(z, roots[i]) for i in range(deg)])
    f_corrected = lambda z: G(z) * np.log(np.abs(z - a)) / np.abs(np.polyval(f, z)) - correction_1(z)

    boxes = _make_boxes(roots, cluster_threshold=cluster_threshold)
    int_re, int_im, err = _integrate_plane(f_corrected, boxes)

    result = int_re + 1j * int_im + sum([np.pi ** 1.5 * const_1[i] for i in range(deg)])
    return result, err, 0.0
# Note:
# int_\C 1/|z| e^{-|z|^2} dxdy = Pi^{3/2}
# = np.pi**1.5


# Compute the integral G(z) /|f(z)| over the complex plane
def regulator_integral_2(f, G=lambda z: 1, cluster_threshold=0.4):
    f = _coeffs(f)
    deg = len(f) - 1
    assert deg == 5 or deg == 6

    roots = np.roots(f)
    fder = np.polyder(f)

    const_1 = [G(roots[i]) / np.abs(np.polyval(fder, roots[i])) for i in range(deg)]

    bump_1 = lambda z, z0: np.exp(-np.abs(z - z0) ** 2) / np.abs(z - z0)
    correction_1 = lambda z: sum([const_1[i] * bump_1(z, roots[i]) for i in range(deg)])
    f_corrected = lambda z: G(z) / np.abs(np.polyval(f, z)) - correction_1(z)

    boxes = _make_boxes(roots, cluster_threshold=cluster_threshold)
    int_re, int_im, err = _integrate_plane(f_corrected, boxes)

    result = int_re + 1j * int_im + sum([np.pi ** 1.5 * const_1[i] for i in range(deg)])
    return result, err, 0.0
# Note:
# int_\C 1/|z| e^{-|z|^2} dxdy = Pi^{3/2}
# = np.pi**1.5


# Compute the integral G(z) log|z-a| /|f(z)| over the complex plane
def regulator_integral(f, a=1, G=lambda z: 1, tol=1e-12, cluster_threshold=0.4):
    f = _coeffs(f)
    if np.polyval(f, a) == 0:
        return regulator_integral_0(f, a, G, tol, cluster_threshold)
    else:
        return regulator_integral_1(f, a, G, cluster_threshold)


def _weight_vector(integrate, f, **kwds):
    """Genus-two weights, in the order (1, Re(z), |z|^2)."""
    values = [integrate(f, G=weight, **kwds)
              for weight in (lambda z: 1, np.real, lambda z: np.abs(z)**2)]
    return [np.real(v) for v, _, _ in values], [err for _, err, _ in values]


def regulator_vector(f, a=1, tol=1e-12, cluster_threshold=0.4):
    r"""Integrals of `(1, \Re(z), |z|^2) \log|z-a| / |f(z)|` over the plane."""
    return _weight_vector(regulator_integral, _coeffs(f), a=a, tol=tol,
                          cluster_threshold=cluster_threshold)


def regulator_vector_2(f, cluster_threshold=0.4):
    r"""Integrals of `(1, \Re(z), |z|^2) / |f(z)|` over the plane."""
    return _weight_vector(regulator_integral_2, _coeffs(f),
                          cluster_threshold=cluster_threshold)


def _period_digits(epsabs, epsrel):
    """The decimal accuracy target corresponding to the integration tolerances."""
    tolerances = np.asarray([epsabs, epsrel], dtype=float)
    if not np.isfinite(tolerances).all() or not np.any(tolerances > 0):
        raise ValueError('period tolerances must be finite with at least one positive value')
    return max(15, int(np.ceil(-np.log10(min(tolerances[tolerances > 0])))))


def period_matrix(f, epsabs=1.49e-8, epsrel=1.49e-8, *, certified=False, prec=None):
    r"""Full cycle periods of `x^i dx/y` as a complex128 array.

    Rows are twice the lifted tree-edge integrals; in even degree omit the
    leaf edge of largest modulus. At least one tolerance must be positive;
    the internal target is at least 15 digits, with output limited to float64.
    """
    from .precision import period_matrix_mp

    if certified:
        return period_matrix_mp(f, _period_digits(epsabs, epsrel) if prec is None else prec,
                                certified=True)
    return np.asarray(period_matrix_mp(f, _period_digits(epsabs, epsrel)), complex)


def cplus(f, epsabs=1.49e-8, epsrel=1.49e-8, *, certified=False, prec=None):
    """The signed real minor of full cycle periods."""
    return _c_invariants(f, epsabs, epsrel, certified=certified, prec=prec)[0]


def cminus(f, epsabs=1.49e-8, epsrel=1.49e-8, *, certified=False, prec=None):
    """The signed imaginary minor of full cycle periods."""
    return _c_invariants(f, epsabs, epsrel, certified=certified, prec=prec)[1]


def _c_invariants(f, epsabs=1.49e-8, epsrel=1.49e-8, *, certified=False, prec=None):
    """Both signed minors from one period-matrix integration."""
    from .precision import c_invariants_mp

    if certified:
        return c_invariants_mp(f, _period_digits(epsabs, epsrel) if prec is None else prec,
                               certified=True)
    return tuple(float(v) for v in c_invariants_mp(f, _period_digits(epsabs, epsrel)))
