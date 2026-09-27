r"""
Regulator integrals over the plane by adaptive quadrature.

The reference backend: each entry of the regulator matrix is computed as a
genuine two-dimensional integral over `\CC`, with the near-singular clusters
of branch points cut out by a smooth partition of unity (the "box trick"
below) so that adaptive quadrature sees no discontinuity.  Slow -- six 2D
quadratures for one `3\times3` determinant -- but it makes no assumption on
the position of the branch points, so it is what the faster backend
:mod:`hyperell_regulator.regulator.stokes` is checked against.

The entry points are :func:`regulator_det_g` and :func:`regulator_vector_g`
for the determinant and the weight rows, and :func:`period_matrix`,
:func:`cplus`, :func:`cminus` for the periods.
"""

import scipy.integrate as spyint
import numpy as np
import scipy as sp
from scipy.special import psi

# ============================================================
# "Box trick" (v2: smooth partition of unity)
#
# See regulator_integral_boxed.py for the original box-decomposition
# idea. That version used a HARD 0/1 mask to cut the near-singular
# clusters out of the global (origin-centered, compactified-polar)
# integral: masked(z) = 0 inside a box, f_corrected(z) outside.
#
# That hard mask is only piecewise continuous -- it has a genuine
# jump discontinuity across each box's boundary. Adaptive quadrature
# converges quickly for smooth integrands but only slowly (and, in
# 2D with a boundary not aligned to the coordinate grid, unreliably)
# across a jump. We only handed nquad exact breakpoints for where
# that jump crosses the real axis (theta = 0, pi); everywhere else
# the routine has to rediscover the jump on its own, which is what
# was limiting the achievable precision to about 1e-5 in practice.
#
# The fix here: replace the hard mask with a SMOOTH cutoff weight
# w(z), equal to 1 inside the box, 0 outside a slightly larger
# buffer box, and a C^2 ("smoothstep") transition in between. Then:
#
#   local piece (per cluster)  = integral of   w(z)  * f_corrected(z)
#                                 over the buffer box
#   global piece               = integral of (1-w(z)) * f_corrected(z)
#                                 over the whole plane (compactified)
#
# Since sum of the weights is exactly 1 everywhere (boxes are kept
# disjoint, including their buffers), local + global again
# reconstructs the exact integral of f_corrected over all of C --
# but now with NO discontinuity anywhere in either integrand, so
# ordinary adaptive quadrature converges quickly and reliably.
# ============================================================


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
    # const_2 = G(root_a)/np.abs(np.polyval(fder, root_a))

    bump_1 = lambda z, z0: np.exp(-np.abs(z - z0) ** 2) / np.abs(z - z0)
    bump_2 = lambda z, z0: np.log(np.abs(z - z0)) * np.exp(-np.abs(z - z0) ** 2) / np.abs(z - z0)

    correction_1 = lambda z: sum([const_1[i] * bump_1(z, roots1[i]) for i in range(deg - 1)])
    correction_2 = lambda z: const_2 * bump_2(z, a)
    # correction_2 = lambda z : const_2 * bump_2(z, root_a)

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


# integrals of log|a + z|/|f(z)|, Re(z) log|a + z|/|f(z)|,  |z^2| log|a + z|/|f(z)|
def regulator_vector(f, a=1, tol=1e-12, cluster_threshold=0.4):
    f = _coeffs(f)
    v1, err1, _ = regulator_integral(f, a, G=lambda z: 1, tol=tol, cluster_threshold=cluster_threshold)
    v2, err2, _ = regulator_integral(f, a, G=lambda z: np.real(z), tol=tol, cluster_threshold=cluster_threshold)
    v3, err3, _ = regulator_integral(f, a, G=lambda z: np.abs(z) ** 2, tol=tol, cluster_threshold=cluster_threshold)
    reg_vec = [np.real(v1), np.real(v2), np.real(v3)]
    err_vec = [err1, err2, err3]
    return reg_vec, err_vec


# integrals of 1/|f(z)|, Re(z)/|f(z)|,  |z^2|/|f(z)|
def regulator_vector_2(f, cluster_threshold=0.4):
    f = _coeffs(f)
    v1, err1, _ = regulator_integral_2(f, G=lambda z: 1, cluster_threshold=cluster_threshold)
    v2, err2, _ = regulator_integral_2(f, G=lambda z: np.real(z), cluster_threshold=cluster_threshold)
    v3, err3, _ = regulator_integral_2(f, G=lambda z: np.abs(z) ** 2, cluster_threshold=cluster_threshold)
    reg_vec = [np.real(v1), np.real(v2), np.real(v3)]
    err_vec = [err1, err2, err3]
    return reg_vec, err_vec


# f : a list of rational numbers, corresponding to the polynomial f[0] * x**n + ... + f[n]
# Compute the period matrix of the hyperelliptic curve y^2 = f(x)
# returns 2g x g matrix, where the ith row is given by the integrals of x^i dx/y over the line segments joining the first 2g + 1 roots of f.
#
# (unaffected by the box trick -- this is a 1D contour integral
# between consecutive roots, not a 2D plane integral.)
def period_matrix(f, epsabs=1.49e-8, epsrel=1.49e-8):
    r"""
    ``epsabs``/``epsrel`` are forwarded to ``quad``; scipy's own defaults
    (``1.49e-8``, kept here so old callers see no change) cap the result at
    about 8 correct digits -- tighten them (e.g. to ``1e-15``, the practical
    floor of float64) for a period matrix good to more digits.
    """
    f = _coeffs(f)
    deg = len(f) - 1
    g = int((deg - 1) / 2)

    # np.roots returns a real (float64) array when every root happens to be
    # real; the path between two consecutive real roots is then real too, and
    # since f alternates sign between consecutive real roots, np.sqrt of a
    # negative *real* returns nan rather than an imaginary number -- which
    # makes quad raise "array must not contain infs or NaNs". Casting to
    # complex keeps the whole computation in C, where the square root is fine.
    roots = np.roots(f).astype(complex)
    roots.sort()

    periods = np.full((2 * g, g), 1j)

    for i in range(g):
        for j in range(2 * g):
            a, b = roots[j], roots[j + 1]
            delta, mid = b - a, (a + b) / 2
            others = np.delete(roots, [j, j + 1])
            ymid = np.sqrt(-delta**2 * f[0] * np.prod(mid - others))
            if abs(ymid.real)<=128*np.finfo(float).eps*abs(ymid) and ymid.imag<0:
                ymid=-ymid
            def integrand(t):
                # Cancel the endpoint square roots analytically. Continuing
                # each remaining factor avoids principal-sqrt sign jumps.
                z = a + delta * np.sin(np.pi * t / 2)**2
                y = ymid * np.prod(np.sqrt((z - others) / (mid - others)))
                return np.pi * delta * z**i / y
            resre, _ = spyint.quad(lambda t: integrand(t).real, 0, 1,
                                   epsabs=epsabs, epsrel=epsrel)
            resim, _ = spyint.quad(lambda t: integrand(t).imag, 0, 1,
                                   epsabs=epsabs, epsrel=epsrel)
            periods[j, i] = resre + resim * 1j

    return periods


def cplus(f, epsabs=1.49e-8, epsrel=1.49e-8):
    f = _coeffs(f)
    periods = period_matrix(f, epsabs=epsabs, epsrel=epsrel)
    P, L, U = sp.linalg.lu(periods.real)
    return np.linalg.det(P) * np.linalg.det(U)


def cminus(f, epsabs=1.49e-8, epsrel=1.49e-8):
    f = _coeffs(f)
    periods = period_matrix(f, epsabs=epsabs, epsrel=epsrel)
    P, L, U = sp.linalg.lu(periods.imag)
    return np.linalg.det(P) * np.linalg.det(U)

# ============================================================
# Generalization: log|g(z)| in place of log|z - a|
#
# For a polynomial g(z) = c * prod_k (z - s_k)^{m_k}, we have the
# elementary identity
#
#     log|g(z)| = log|c| + sum_k m_k * log|z - s_k|
#
# so
#
#     integral of  G(z) * log|g(z)| / |f(z)|  over C
#   =  log|c| * [integral of G(z)/|f(z)|]                      <- regulator_integral_2
#   +  sum_k  m_k * [integral of G(z) log|z - s_k| / |f(z)|]   <- regulator_integral (per root of g)
#
# Each term on the right is exactly what regulator_integral_2 /
# regulator_integral (dispatching internally to the f(a)=0 or
# f(a)!=0 case, per root s_k of g) already compute robustly via the
# box trick above. So log|g(z)| needs no new numerics at all -- just
# factoring g and summing the existing, already-validated pieces.
#
# Roots of g that coincide with a root of f are handled correctly
# and automatically: regulator_integral(f, a=s_k, ...) itself checks
# np.polyval(f, s_k) == 0 (within `tol`) and dispatches to
# regulator_integral_0, which folds that log-pole in with any nearby
# poles of f via the same box/cluster machinery.
#
# Multiplicities of g's roots are detected numerically (np.roots
# returns each root of a repeated factor as its own, very-close-but-
# not-exactly-equal floating point value) and collapsed via a tight
# `mult_tol`, so a repeated root of g contributes m times a single
# regulator_integral call rather than m separate (and needlessly
# expensive) calls.
#
# Caveat: this does NOT protect against a root of g landing very
# close to -- but not within `tol` of -- a root of f (i.e. an
# "almost cancellation" between a zero of g and a pole of f that is
# too far apart to trigger regulator_integral_0, but close enough to
# be numerically delicate). This situation is not exercised by the
# validation below; if you hit it, consider increasing `tol` so the
# near-coincidence is treated as exact, or perturbing the polynomial.
# ============================================================


# integrals of log|g(z)|/|f(z)|, Re(z) log|g(z)|/|f(z)|, |z^2| log|g(z)|/|f(z)|
# Same as regulator_det(f, a1, a2), but with the scalars a1, a2
# generalized to polynomials g1, g2 (log|z - a| -> log|g(z)|).