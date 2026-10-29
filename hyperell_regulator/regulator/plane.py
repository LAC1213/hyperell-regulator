r"""
Regulator integrals over the plane by adaptive quadrature.

The reference backend: each entry of the regulator matrix is computed as a
genuine two-dimensional integral over `\CC`, with the near-singular clusters
of branch points cut out by a smooth partition of unity (the "box trick"
below) so that adaptive quadrature sees no discontinuity.  Slow -- six 2D
quadratures for one `3\times3` determinant -- but it makes no assumption on
the position of the branch points, so it is what the faster backend
:mod:`hyperell_regulator.regulator.stokes` is checked against.

Use :func:`regulator_vector` for logarithmic weights and
:func:`regulator_vector_2` for unweighted periods. :func:`period_matrix`,
:func:`cplus` and :func:`cminus` compute full cycle periods and their minors.
"""

import scipy.integrate as spyint
import numpy as np
from scipy.special import psi

# Integrate close branch points in local boxes with a C^2 cutoff.
# The complementary weight covers the rest of the plane. Their sum is one,
# so quadrature sees no discontinuity at box boundaries.


def _find_clusters(points, threshold):
    """Union-find clustering of `points` (list/array of complex
    numbers) by pairwise distance < threshold. Returns a list of
    index-groups of size >= 2 -- singleton points don't need any
    special treatment."""
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
    """C^2 smoothstep on [0,1]: 0 at u=0, 1 at u=1, with zero first
    AND second derivative at both ends (so the weight below is C^2
    across both transition boundaries, not just C^0)."""
    u = np.clip(u, 0.0, 1.0)
    return 6 * u**5 - 15 * u**4 + 10 * u**3


def _make_boxes(points, cluster_threshold=0.4, max_half_side=0.6, buffer_factor=1.3):
    """Build one axis-aligned square box per cluster of close
    singular points among `points`. Returns a list of tuples
    (centroid, R, R_buf, x_breaks, y_breaks):
        R      -- half-side of the "core" box (weight == 1 inside)
        R_buf  -- half-side of the buffer box (weight == 0 outside);
                  R_buf = buffer_factor * R
        x_breaks / y_breaks -- coordinates of the cluster's own
                  points relative to the centroid, for quadrature
                  breakpoints inside the core box.
    Boxes (including their buffers) are kept disjoint from every
    other singular point and from each other.
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
    """Smooth cutoff weight for one box: 1 inside the core square,
    0 outside the buffer square, C^2 transition in between. Uses
    the L-infinity (Chebyshev) distance so the core/buffer regions
    are the axis-aligned squares used elsewhere."""
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
    """Integrate f_corrected (C -> C, already had all analytic bump
    corrections subtracted, so it's bounded on C) over the whole
    complex plane: smoothly-weighted local Cartesian quadrature
    inside each box's buffer region, plus the compactified-polar
    quadrature (centered at the origin) with the complementary
    smooth weight over the rest of the plane.
    Returns (real_part, imag_part, summed_error_estimate)."""

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
    r"""
    ``f`` as a list of coefficients in descending degree.

    This integrator is numpy's, so it wants coefficients rather than a
    polynomial; the package's entry points take Sage polynomials, and they
    are converted here, once, on the way in.
    """
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


# f : a list of rational numbers, corresponding to the polynomial f[0] * x**n + ... + f[n]
# Compute the period matrix of the hyperelliptic curve y^2 = f(x)
def period_matrix(f, epsabs=1.49e-8, epsrel=1.49e-8):
    r"""
    Full cycle periods along a validated semilunar branch-point chain.

    Start at the finite root of largest absolute value (ties by real then
    imaginary part), then follow increasing argument about a perturbed
    centroid. Detours use the same path builder and validation as regulator
    cuts. Each row is twice the integral along one lifted branch-to-branch
    path, for the differentials ``x^i dx/y``, ``0 <= i < g``.

    ``epsabs``/``epsrel`` are forwarded to SciPy's ``quad``. Its defaults
    give roughly eight digits; tighten them for more float64 accuracy.
    """
    from ._periods import period_paths, segment_differential

    f = _coeffs(f)
    g = (len(f) - 2) // 2
    roots = np.roots(f).astype(complex)
    periods = []
    for path in period_paths(roots, g, root_tol=128*np.finfo(float).eps):
        row = np.zeros(g, complex)
        previous_y = None
        for (a, ia), (b, ib) in zip(path[:-1], path[1:]):
            evaluate, previous_y = segment_differential(
                roots, complex(f[0]), a, b, ia, ib, g,
                sqrt=np.sqrt, sin=np.sin, cos=np.cos, pi=np.pi,
                sign_tol=128*np.finfo(float).eps, previous_y=previous_y)
            for i in range(g):
                re, _ = spyint.quad(lambda t: evaluate(t)[i].real, 0, 1,
                                    epsabs=epsabs, epsrel=epsrel)
                im, _ = spyint.quad(lambda t: evaluate(t)[i].imag, 0, 1,
                                    epsabs=epsabs, epsrel=epsrel)
                row[i] += re + 1j * im
        periods.append(2 * row)
    return np.asarray(periods, complex)


def _period_minor(periods, epsabs, epsrel):
    from ._periods import pivoted_minor

    return pivoted_minor(periods, max(epsabs, epsrel, 128*np.finfo(float).eps))


def cplus(f, epsabs=1.49e-8, epsrel=1.49e-8):
    """Signed real minor of full cycle periods, with deterministic pivot ties."""
    return _period_minor(period_matrix(f, epsabs, epsrel).real, epsabs, epsrel)


def cminus(f, epsabs=1.49e-8, epsrel=1.49e-8):
    """Signed imaginary minor of full cycle periods, with deterministic pivot ties."""
    return _period_minor(period_matrix(f, epsabs, epsrel).imag, epsabs, epsrel)


def _c_invariants(f, epsabs=1.49e-8, epsrel=1.49e-8):
    """Both signed minors from one period-matrix integration."""
    periods = period_matrix(f, epsabs, epsrel)
    return tuple(float(_period_minor(part, epsabs, epsrel))
                 for part in (periods.real, periods.imag))
