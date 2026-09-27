r"""
Integration of `\omega_a`, `\log\varphi\,\omega_a` and the second moments
along the edges of the cut graph.

An edge is a lift of a straight segment in the `w`-plane, and is integrated in
`w` -- except near its ends, which are ramification points of `\varphi`, where
`dx/dw` blows up.  There the edge is integrated in a local coordinate of the
curve instead, along a straight path, which is legitimate because the
integrands are holomorphic and only the homotopy class matters:

======================  ===============================================
target of the leg       chart
======================  ===============================================
`y \neq 0`, `x` finite  `x`; `y' = F'/2y`, `w = \varphi(x,y)`
`y = 0`                 `y`; `x' = 2y/F'`, regular at the branch point
`x = \infty`            `\tau`, `x = \tau^{-2}`, `y = \eta\tau^{-(2g+1)}`
======================  ===============================================

At `P` and at `Q`, where `\varphi` vanishes or blows up to order `N`,
`\lambda = \log\varphi` diverges logarithmically.  Those legs run down their
straight path exponentially, `c = c_* + (c_0 - c_*)e^{-\xi}`, which makes
`\lambda` linear in `\xi` while `\omega_a` decays, so the truncated tail is
negligible.  `w` itself is never carried as a state there: near `P` it
underflows and `d\lambda = dw/w` becomes noise, so it is recovered from the
norm, `w = \operatorname{N}(\varphi)/(A - By)`, and near `Q` from the
`\tau`-chart's `\hat w = w\tau^N`.
"""

import numpy as np
from scipy.integrate import solve_ivp

from . import walks
from .cover import Cover, _shift_poly, _vanishes
from .diagnostics import checkpoint, phase, traced

RTOL, ATOL = 1e-12, 1e-14
# How far the tolerance may be relaxed when 1e-12 is below what the
# right-hand side can actually deliver: see :func:`_solve`.
RTOL_FLOOR = 1e-8
XI = 45.0


class LegFailure(RuntimeError):
    """The chart leg strayed off its branch; the caller shortens it."""


class QuadratureFailure(ArithmeticError):
    """The integral did not meet its tolerance within the refinement budget.

    This is not a chart failure: shortening an endpoint tail does not
    resolve a nearby interior singularity of the analytically continued
    integrand. In particular it must escape the stand-off retry loop.
    """


@traced("ode_solve", fields=("span",))
def _solve(rhs, span, y0, **kwds):
    r"""
    ``solve_ivp`` with the tolerance relaxed until the right-hand side can
    meet it.

    ``rtol = 1e-12`` is below the accuracy of the right-hand side on a cover
    with a wide dynamic range, and there the integration does not fail
    gracefully: it reports *"required step size is less than spacing between
    numbers"*, which reads like a singularity on the path and is nothing of
    the kind.  On the degree-6 model of LMFDB 249.a.249.1 the fibre over an
    edge midpoint carries `|y|` up to `1.2\cdot10^3` and `|F'(x)|` up to
    `8\cdot10^5` while `y\,d\varphi/dx` is `O(10^2)`, so four or five digits
    go in cancellation and the rates are good to about `10^{-11}`.  Asked for
    `10^{-12}` the solver dies after 446 evaluations having advanced
    `2\cdot10^{-7}` of the way; asked for `10^{-10}` it completes, and for
    `10^{-8}` it completes in 266.  Every stand-off then fails for the same
    reason, whatever the geometry, and the caller reports that it can find
    none -- which is how a tolerance that is merely too tight came to look
    like a cut system that could not be integrated.

    So a failure of exactly that kind is retried a hundredfold looser, down
    to :data:`RTOL_FLOOR`; anything else is raised at once.  The tolerance
    that succeeded is returned with the solution, since it caps what the leg
    is worth.
    """
    rtol, atol = kwds.pop("rtol", RTOL), kwds.pop("atol", ATOL)
    def counted_rhs(t, y):
        checkpoint("ode_evaluations", t=t)
        return rhs(t, y)

    while True:
        checkpoint("ode_attempts", rtol=rtol, atol=atol)
        sol = solve_ivp(counted_rhs, span, y0, rtol=rtol, atol=atol, **kwds)
        if sol.success:
            return sol, rtol
        if "step size" not in sol.message or rtol >= RTOL_FLOOR:
            raise LegFailure(sol.message)
        rtol, atol = rtol * 100.0, atol * 100.0


class Leg:
    r"""
    What one leg contributes, per sheet: `dI[a] = \int\omega_a`,
    `dL[a] = \int\lambda\,\omega_a`, and `k[a][b] = \int(\int\omega_a)
    \bar\omega_b`, or ``None`` when the moments were not asked for.
    """

    __slots__ = ("dI", "dL", "k", "frac", "charts")

    def __init__(self, dI, dL, k):
        self.dI, self.dL, self.k = dI, dL, k
        self.frac = None        # the stand-off :func:`side` settled on
        self.charts = None      # verified preflight chart choices, if present

    def reversed(self):
        k = None
        if self.k is not None:
            k = -self.k + np.einsum("a...,b...->ab...", self.dI, np.conj(self.dI))
        return Leg(-self.dI, -self.dL, k)

    def __add__(self, o):
        k = None
        if self.k is not None and o.k is not None:
            k = (self.k + np.einsum("a...,b...->ab...", self.dI, np.conj(o.dI))
                 + o.k)
        return Leg(self.dI + o.dI, self.dL + o.dL, k)

    @staticmethod
    def stack(legs):
        """One leg per sheet into a single leg over all of them."""
        dI = np.stack([L.dI for L in legs], axis=-1)
        dL = np.stack([L.dL for L in legs], axis=-1)
        k = (None if legs[0].k is None
             else np.stack([L.k for L in legs], axis=-1))
        return Leg(dI, dL, k)


def _pack(g, n, p0, lam0, moments):
    nb = 1 + g + 1 + g + (g * g if moments else 0)
    y0 = np.zeros(nb * n, complex)
    y0[:n] = p0
    y0[(1 + g) * n:(2 + g) * n] = lam0
    return np.concatenate([y0.real, y0.imag]), nb


def _run(rhs, span, g, n, p0, lam0, moments):
    y0, nb = _pack(g, n, p0, lam0, moments)
    sol, _ = _solve(rhs, span, y0, method="DOP853")
    y = sol.y[:nb * n, -1] + 1j * sol.y[nb * n:, -1]
    b = [y[i * n:(i + 1) * n] for i in range(nb)]
    dI = np.stack(b[1:1 + g])
    dL = np.stack(b[2 + g:2 + 2 * g])
    k = None
    if moments:
        k = np.stack([np.stack(b[2 + 2 * g + a * g:2 + 2 * g + (a + 1) * g])
                      for a in range(g)])
    return b[0], b[1 + g], Leg(dI, dL, k)


def _rates(g, I, dI, lam, dlam, moments):
    out = list(dI) + [dlam] + [lam * d for d in dI]
    if moments:
        out += [I[a] * np.conj(dI[b]) for a in range(g) for b in range(g)]
    return out


def _wrap(g, n, moments, body):
    """Split the real/imaginary packing off the chart-specific body."""
    nb = 1 + g + 1 + g + (g * g if moments else 0)

    def rhs(s, yr):
        y = yr[:nb * n] + 1j * yr[nb * n:]
        b = [y[i * n:(i + 1) * n] for i in range(nb)]
        I = b[1:1 + g]
        dp, dI, dlam = body(s, b[0], I, b[1 + g])
        d = np.concatenate([dp] + _rates(g, I, dI, b[1 + g], dlam, moments))
        return np.concatenate([d.real, d.imag])

    return rhs


# --------------------------------------------------------------- the charts

def _tanh_sinh(level, lam=np.pi / 2, reach=24.0):
    r"""
    Nodes in `(0,1)` and weights for `\int_0^1 f(t)\,dt`, by the
    double-exponential change of variable `t = (1 + \tanh(\lambda\sinh\tau))/2`.

    The leg stops a stand-off short of a critical value, and just beyond that
    end the integrand grows like `(w-z)^{1/e-1}`: `\alpha = 1/2` at a simple
    branch point and `1 - 1/N` at the `N`-fold point `P`.  Gauss-Legendre is
    spectral on an analytic integrand and only algebraic on that one, which is
    why raising its order does nothing -- 8, 16 and 24 all leave the closing
    defect at `8\cdot10^{-6}` -- and why no mesh drawn by hand fixed it
    either.  The double-exponential substitution pushes the singularity to
    infinity and clusters the nodes at the ends doubly exponentially, which is
    what Molin and Neurohr integrate these periods with (Math. Comp. 88,
    2019, §6).

    ``reach`` is where `\cosh^{-2}` has underflowed, so the nodes past it
    contribute nothing.
    """
    h = 1.0 / (1 << int(level))
    kmax = int(np.ceil(np.arcsinh(reach / lam) / h))
    tau = np.arange(-kmax, kmax + 1) * h
    sh = lam * np.sinh(tau)
    u = np.tanh(sh)
    w = h * lam * np.cosh(tau) / np.cosh(sh) ** 2
    nodes = (1.0 + u) / 2.0
    # The affine map can round an interior tanh value to exactly 1.
    keep = (nodes > 0.0) & (nodes < 1.0)
    return nodes[keep], w[keep] / 2.0


@traced("walk_quadrature", fields=("level", "logarithmic"))
def _walk_quadrature(cov, coord, rate, x, y, lam, charts, level,
                     cap=20000, wof=None, logarithmic=False, max_level=10,
                     chosen=None, breaks=None):
    """Refine the quadrature independently of the sheet-tracking steps.

    A well-resolved fibre does not imply a well-resolved integral.  On the
    degree-14 cover the level-5 rule misses 3e-5 on one finite edge, even
    though the full preimage sums to zero at every quadrature node.
    Nearby critical values off the path can require further levels: the
    degree-5 conductor-277 ray needs level 9, despite successful tracking.
    ``max_level`` bounds quadrature work separately from chart retries;
    the active IntegrationMonitor also enforces its time and step budgets.
    """
    if int(max_level) <= int(level):
        raise ValueError("max_level must allow at least one refinement")
    end = -np.log1p(-np.nextafter(1.0, 0.0)) if logarithmic else 1.0
    if chosen is None:
        with phase("choose_charts"):
            chosen = walks.choose(cov, coord, x, y, charts, end=end)
    x, y, lam = [np.atleast_1d(np.asarray(v, complex)) for v in (x, y, lam)]
    n = len(x)
    samples = [{} for _ in range(n)]
    active = np.arange(n)
    # Each sheet retains its own verified samples and last converged result.
    # A difficult sheet must not force every other sheet onto a finer rule.
    result = [np.empty(n, complex) for _ in range(3)]
    result += [np.empty((cov.g, n), complex) for _ in range(2)]
    previous = None
    last_error = None
    for refinement in range(int(level), int(max_level) + 1):
        nodes, wts = _tanh_sinh(refinement)
        if breaks is not None:
            intervals = list(zip(breaks[:-1], breaks[1:]))
            wts = np.concatenate([(b-a)*wts for a,b in intervals])
            nodes = np.concatenate([a+(b-a)*nodes for a,b in intervals])
        if logarithmic:
            wts = wts / (1.0 - nodes)
            nodes = -np.log1p(-nodes)
        sheet_ids = active.tolist()
        with phase("quadrature_level", level=refinement, nodes=len(nodes),
                   active_sheets=sheet_ids):
            out = walks.leg(cov, coord, rate, x[active], y[active], lam[active],
                            nodes, wts, charts=charts, cap=cap, wof=wof,
                            end=float(nodes[-1]) if logarithmic else end,
                            chosen=[chosen[i] for i in active],
                            samples=[samples[i] for i in active],
                            sheet_ids=sheet_ids)
        converged = np.zeros(len(active), dtype=bool)
        if previous is not None:
            converged[:] = True
            last_error = []
            for value, old in zip(out[3:], previous):
                error = np.max(np.abs(value - old[:, active]), axis=0)
                scale = np.maximum(1.0, np.max(np.abs(value), axis=0))
                converged &= error <= 1e-11 * scale
                last_error.append(float(np.max(error / scale)))
            checkpoint(level=refinement, scaled_quadrature_errors=last_error)
        for full, value in zip(result, out):
            full[..., active] = value
        checkpoint("converged_sheets", amount=int(np.count_nonzero(converged)),
                   level=refinement, converged_sheets=active[converged].tolist())
        active = active[~converged]
        if not len(active):
            return tuple(result)
        # Copies retain the previous level until every comparison is made;
        # result itself is updated only at the active original sheet indices.
        previous = [value.copy() for value in result[3:]]
    raise QuadratureFailure(
        "walk quadrature did not converge by level %d on sheets %s; "
        "last scaled period/log errors=%s (tolerance 1e-11)" %
        (max_level, active.tolist(), last_error))



def _pinch_detour(cov, wa, wb):
    """Avoid projection collisions inside disks containing no branch value.

    At B=0, two distinct curve points share x. The map phi is unramified,
    so moving the base path inside a disk free of critical values and zero
    preserves both its lifted endpoints and the integrals of omega/log(phi).
    This keeps the existing root-separation test meaningful.
    """
    if getattr(cov, "exact", None) is None or cov.pure_x or len(cov.b) < 2:
        return None
    if not hasattr(cov, "_ordinary_pinches"):
        from sage.all import QQ, CC, PolynomialRing
        S = PolynomialRing(QQ, ["x", "w"])
        T = PolynomialRing(QQ, "W")
        pinch, _ = cov._pinch_polynomials(S, T)
        if pinch is None or pinch.degree() < 1:
            cov._ordinary_pinches = np.array([], complex)
        else:
            branch = cov._branch_polynomial()
            pinch = pinch // pinch.gcd(branch)
            cov._ordinary_pinches = np.array(
                [complex(z) for z in pinch.roots(CC, multiplicities=False)])
    points = cov._ordinary_pinches
    if not len(points) or wa == wb:
        return None
    wa, wb = complex(wa), complex(wb)
    length = abs(wb-wa)
    direction = (wb-wa)/length
    if not hasattr(cov, "_pinch_obstacles"):
        cov._pinch_obstacles = np.r_[cov.critical_values(), 0j]
    obstacles = cov._pinch_obstacles
    bends = []
    for i, z in enumerate(points):
        position = (z-wa)/direction
        if not 0 < position.real < length:
            continue
        distances = [position.real, length-position.real,
                     float(np.min(abs(obstacles-z)))]
        if len(points) > 1:
            distances.append(float(np.min(abs(np.delete(points,i)-z))))
        # Even the rectangle's farthest corner is inside half the
        # distance to every excluded value. Use that clearance to keep
        # colliding x-roots well separated along the detour.
        radius = 0.25*min(distances)
        if radius <= 0 or abs(position.imag) >= radius/2:
            continue
        center = z-1j*position.imag*direction
        normal = (1j if position.imag <= 0 else -1j)*direction
        bends.append((position.real, [center-radius*direction,
                     center-radius*direction+radius*normal,
                     center+radius*direction+radius*normal,
                     center+radius*direction]))
    if not bends:
        return None
    path = [wa]
    for _, vertices in sorted(bends, key=lambda item:item[0]):
        path.extend(vertices)
    return path+[wb]


@traced("leg_in_w", fields=("wa", "wb", "integrate", "moments"))
def leg_in_w(cov, wa, wb, x0, y0, lam0, moments=False, level=5, cap=20000,
             toward=None, integrate=True, _avoid_pinches=True, _prepared=None,
             _charts=None):
    r"""
    Track the fibre along a straight w-segment and integrate its lifts.

    Newton solves the cover equation at each sample, with adaptive step
    and root-separation checks in :mod:`.walks`. Near a critical endpoint,
    w=z+s^e desingularizes the path parameter. Nested quadrature rules
    resolve the integral independently of the continuation mesh.

    ``moments`` wants the running primitive at every node, which a panel rule
    does not give, and keeps the ODE; that leaves the slow form of the
    regulator an independent check on the fast one.
    """
    if moments:
        return _leg_in_w_ode(cov, wa, wb, x0, y0, lam0, moments)
    path = _pinch_detour(cov, wa, wb) if _avoid_pinches else None
    if path is not None:
        checkpoint("pinch_detours", segments=len(path)-1)
        # Center the continuation polynomial at the projection collision.
        # A multiple zero of B otherwise loses root separation by cancellation
        # even on the safe detour (notably x=-1 on conductor 295).
        if not hasattr(cov, "_pinch_x_values"):
            from sage.all import CC
            F, a, ad, b, bd = cov.exact
            squarefree = b // b.gcd(b.derivative())
            cov._pinch_x_values = [(complex(a(r)/ad(r)), complex(r))
                                  for r in squarefree.roots(CC, multiplicities=False)]
            cov._pinch_charts = {}
        centered = []
        for j in range(1, len(path)-1, 4):
            value = (path[j]+path[j+3])/2
            _, center = min(cov._pinch_x_values, key=lambda item:abs(item[0]-value))
            if center.imag == 0:
                if center not in cov._pinch_charts:
                    cov._pinch_charts[center] = walks.Chart(cov, "x", center)
                centered.append(cov._pinch_charts[center])
        target = toward[0] if toward is not None else wb
        charts = centered + walks.charts_over(cov, target)
        x, y, lam, total = x0, y0, lam0, None
        prepared = []
        for j, (a, b) in enumerate(zip(path, path[1:])):
            x, y, lam, part = leg_in_w(
                cov, a, b, x, y, lam, level=min(level, 3), cap=cap, toward=toward,
                integrate=integrate, _avoid_pinches=False,
                _prepared=None if _prepared is None else _prepared[j],
                _charts=charts)
            prepared.append(part.charts)
            total = part if total is None else total+part
        total.charts = prepared
        return x, y, lam, total
    g = cov.g
    x = np.atleast_1d(np.asarray(x0, complex)).copy()
    y = np.atleast_1d(np.asarray(y0, complex)).copy()
    lam = np.atleast_1d(np.asarray(lam0, complex)).copy()
    n = len(x)
    dw = complex(wb) - complex(wa)
    if dw == 0:
        return x, y, lam, Leg(np.zeros((g, n), complex),
                              np.zeros((g, n), complex), None)

    # In w the integrand is not analytic on the leg: it stops a stand-off
    # short of a critical value, and just beyond that end x^a/(y dphi/dx)
    # grows like (w-z)^{1/e-1} -- at the N-fold point P, with e = N = 7, that
    # is w^{-6/7}.  Gauss-Legendre is spectral on an analytic integrand and
    # only algebraic on that one, which is why raising its order does
    # nothing: 8, 16 and 24 all leave the closing defect at 8e-6.  In the
    # local uniformizer s, w = z + s^e, the integrand is regular -- omega_a
    # is a holomorphic differential and s a uniformizer, so omega_a/ds is
    # holomorphic -- and the same rule is spectral again.  The leg runs along
    # a ray from z, so a straight path in s covers the same segment in w.
    coord = rate = None
    breaks = None
    if toward is not None and toward[1] > 1:
        zc, e = complex(toward[0]), int(toward[1])
        ratio = (complex(wb) - zc) / (complex(wa) - zc)
        if abs(ratio.imag) <= 1e-9 * abs(ratio) and ratio.real > 0:
            sa = (complex(wa) - zc) ** (1.0 / e)
            sb = sa * ratio.real ** (1.0 / e)
            ds = sb - sa
            # Preserve the small endpoint when |sa| >> |sb|. Computing
            # sa+t*(sb-sa) loses sb through cancellation on long edges.
            parameter = lambda t: (1.0-t)*sa + t*sb
            coord = lambda t: zc + parameter(t) ** e
            rate = lambda t: e * parameter(t) ** (e - 1) * ds
            if ratio.real < 1e-6:
                # A long leg can span many radial scales. Exponential
                # spacing resolves each scale without squeezing the last
                # part into a tiny boundary layer of the quadrature rule.
                # This is the same straight w-path, with positive radius.
                log_ratio = np.log(ratio.real)
                radial = lambda t: (complex(wa)-zc)*np.exp(t*log_ratio)
                coord = lambda t: zc + radial(t)
                rate = lambda t: log_ratio*radial(t)
                # Nearby branch values (and the logarithm at zero) create
                # narrow peaks even after radial scaling. Put their closest
                # projections at panel boundaries, where tanh-sinh clusters.
                projected = np.real((np.r_[cov.critical_values(), 0j]-zc)
                                    / (complex(wa)-zc))
                cuts = [float(np.log(r)/log_ratio) for r in projected
                        if ratio.real < r < 1.0]
                breaks = sorted(set([0.0, *cuts, 1.0]))
                level = min(level, 3)
    if coord is None:
        coord = lambda t: (1.0-t)*complex(wa) + t*complex(wb)
        rate = lambda t: dw
    zc = complex(toward[0]) if toward is not None else complex(wb)
    try:
        charts = walks.charts_over(cov, zc) if _charts is None else _charts
        if not integrate:
            # Test the proposed hand-over before paying for its integral.
            # side() may reject several stand-offs while finding a valid
            # endpoint chart; those rejected paths need only continuation.
            ends = []
            for i, (xx, yy) in enumerate(zip(x, y)):
                with phase("preflight", sheet=i):
                    ends.append(walks.race(cov, coord, xx, yy, charts))
            placeholder = Leg(np.zeros((g, n), complex), np.zeros((g, n), complex), None)
            placeholder.charts = [v[3] for v in ends]
            return (np.array([v[0] for v in ends]),
                    np.array([v[1] for v in ends]),
                    lam + np.log(complex(wb) / complex(wa)), placeholder)
        x, y, lam, I, L = _walk_quadrature(
            cov, coord, rate, x, y, lam, charts, level, cap=cap,
            chosen=_prepared, breaks=breaks)
    except RuntimeError as exc:
        raise LegFailure(str(exc))
    return x, y, lam, Leg(I, L, None)


def _leg_in_w_ode(cov, wa, wb, x0, y0, lam0, moments=False, **_):
    r"""
    The same leg by marching `dx/dw = 1/(d\varphi/dx)`, kept for the moments.

    Faster than the tracker and without its guarantee; see :func:`leg_in_w`.
    """
    g, n = cov.g, len(np.atleast_1d(x0))
    dw = wb - wa
    # y is a state as well: (x, y) both move, and y is not a function of
    # (x, w) when B = 0
    nb = 1 + g + 1 + g + (g * g if moments else 0)

    def rhs(s, yr):
        z = yr[:(nb + 1) * n] + 1j * yr[(nb + 1) * n:]
        b = [z[i * n:(i + 1) * n] for i in range(nb + 1)]
        x, y = b[0], b[nb]
        I = b[1:1 + g]
        ydphi = cov.y_dphi_dx(x, y)
        dxdw = y * dw / ydphi
        dydw = np.polyval(cov.Fp, x) * dw / (2 * ydphi)
        d = np.concatenate([dxdw] + _rates(
            g, I, [x ** a * dw / ydphi for a in range(g)], b[1 + g],
            np.full(n, dw / (wa + s * dw)), moments) + [dydw])
        return np.concatenate([d.real, d.imag])

    y0v = np.zeros((nb + 1) * n, complex)
    y0v[:n] = x0
    y0v[(1 + g) * n:(2 + g) * n] = lam0
    y0v[nb * n:] = y0
    sol, _ = _solve(rhs, (0.0, 1.0), np.concatenate([y0v.real, y0v.imag]),
                    method="DOP853")
    z = sol.y[:(nb + 1) * n, -1] + 1j * sol.y[(nb + 1) * n:, -1]
    b = [z[i * n:(i + 1) * n] for i in range(nb + 1)]
    k = None
    if moments:
        k = np.stack([np.stack(b[2 + 2 * g + a * g:2 + 2 * g + (a + 1) * g])
                      for a in range(g)])
    return b[0], b[nb], b[1 + g], Leg(np.stack(b[1:1 + g]),
                                      np.stack(b[2 + g:2 + 2 * g]), k)


@traced("endpoint_x", fields=("xend", "exp"))
def leg_in_x(cov, x0, y0, xend, lam0, moments=False, exp=False):
    r"""
    Straight in `x` down to `x_{\mathrm{end}}`, working in `u = x -
    x_{\mathrm{end}}`.

    `y` is carried as `\log y`, whose rate `F'(u)/2F(u)` never divides by `y`
    and never chooses a branch: `y = \exp(\log y)` is continuous along the
    path by construction and satisfies `y^2 = F` to the solver's accuracy.

    Neither of the two obvious alternatives will do.  Integrating `y` itself
    is what the `w`-leg does, and it only satisfies `y^2 = F(x)` to about
    `10^{-14}`; carrying that on while `F(x)` shrinks past it leaves `x` as
    noise and the rates that divide by it explode.  Taking
    `y = y_0\sqrt{F(x)/F(x_0)}` with the principal root is exact on the
    curve, but only continuous while the path stays clear of the zeros of
    `F`: at the far node of the degree-14 cover on LMFDB 249.a.249.1 the
    ramification point has `y = 0.057`, the straight path passes close
    enough to a zero of `F` for the root to jump, and one of the two sheets
    that meet there comes back with the wrong sign -- two per cent of the
    answer, and the only two sheets whose boundaries then fail to close.
    """
    g = cov.g
    x0, xend = complex(x0), complex(xend)
    cs = cov.shifted(xend)
    u0 = x0 - xend
    y0 = complex(y0)
    nb = 1 + g + 1 + g + (g * g if moments else 0)

    def coord(s):
        if exp:
            return u0 * np.exp(-s), -u0 * np.exp(-s)
        return u0 * (1 - s), -u0

    def rhs(s, yr):
        z = yr[:nb + 1] + 1j * yr[nb + 1:]
        b = [z[i:i + 1] for i in range(nb + 1)]
        ly = b[nb]
        I = b[1:1 + g]
        u, du = coord(s)
        y = np.exp(ly)
        dly = np.atleast_1d(
            np.polyval(cs.Fp, u) / (2 * np.polyval(cs.F, u)) * du)
        dy = y * dly
        dlam = cs.dlogphi(u, y[0], du, dy[0])
        d = np.concatenate([np.atleast_1d(du)] + _rates(
            g, I, [(xend + u) ** a * du / y for a in range(g)], b[1 + g],
            np.atleast_1d(dlam), moments) + [dly])
        return np.concatenate([d.real, d.imag])

    y0v = np.zeros(nb + 1, complex)
    y0v[0] = u0
    y0v[1 + g] = complex(lam0)
    y0v[nb] = np.log(y0)
    span = (0.0, XI) if exp else (0.0, 1.0)
    sol, _ = _solve(rhs, span, np.concatenate([y0v.real, y0v.imag]),
                    method="DOP853")
    z = sol.y[:nb + 1, -1] + 1j * sol.y[nb + 1:, -1]
    b = [z[i:i + 1] for i in range(nb + 1)]
    k = None
    if moments:
        k = np.stack([np.stack(b[2 + 2 * g + a * g:2 + 2 * g + (a + 1) * g])
                      for a in range(g)])
    return b[0], b[1 + g], Leg(np.stack(b[1:1 + g]),
                               np.stack(b[2 + g:2 + 2 * g]), k)


@traced("endpoint_P", fields=("xend",))
def leg_into_P(cov, x0, y0, xend, lam0, moments=False, track=0):
    r"""
    Into the totally ramified point `P`, along the sheet rather than across it.

    The other chart legs follow a path laid down in the plane; this one
    follows the curve.  `\varphi` has a zero of order `N` at `P`, so down the
    ray `w = w_0e^{-Ns}` the sheet satisfies

    .. math:: \frac{du}{ds} = \frac{-N}{(\log\varphi)'(u)},
              \qquad \lambda = \lambda_0 - Ns, \qquad u = x - x_P ,

    which near `P` is `du/ds = -u`: the whole fibre slides into `P`
    exponentially, and `\lambda` needs no integrating at all.

    Why not the `x`-chart, which does this job at every other node: the `N`
    sheets leave `P` as `x - x_P \sim c\,w^{1/N}\zeta_N^k`, so the stand-off
    where the `w`-leg hands over -- as close to `P` as the neighbouring
    critical values allow -- is still `|w|^{1/N}` away in `x`, which at
    `N = 15` is two thirds of a unit.  A straight path in `x` from there is
    not short, and on the degree-15 cover of LMFDB 277.a.277.1 it lands in
    the wrong homotopy class: the legs come back wrong in the first digit and
    the `N` of them fail to sum to zero.

    Two things make this stable where the obvious forms are not.
    `(\log\varphi)'` is taken as `\sum_j m_j/(u + x_P - r_j)` over the zeros
    and poles of `\operatorname{N}(\varphi)`
    (:meth:`~.cover.Cover.dlog_norm_about`), so the `N`-fold zero contributes
    exactly `N/u` rather than a quotient that cancels to nothing below
    `|u| = \varepsilon^{1/N}`.  And `y` is carried as a state: deriving it as
    `y_0\sqrt{F/F_0}` is what the straight-path legs do, but this path bends,
    and where it bends around a zero of `F` the principal root jumps and the
    solver stalls on the discontinuity.

    With ``track`` the leg also returns `(s, u, y)` sampled at that many
    points for diagnostics.
    """
    g, n = cov.g, len(np.atleast_1d(x0))
    N = cov.N
    xend = complex(xend)
    cs = cov.shifted(xend)
    dlogN = cov.dlog_norm_about(xend)
    nb = 1 + g + 1 + g + (g * g if moments else 0)

    def rhs(s, yr):
        z = yr[:(nb + 1) * n] + 1j * yr[(nb + 1) * n:]
        b = [z[i * n:(i + 1) * n] for i in range(nb + 1)]
        u, y = b[0], b[nb]
        I = b[1:1 + g]
        # log(phi) = log N(phi) - log(A - By), the first term by its poles
        big = cs.A(u) - cs.B(u) * y
        dydu = np.polyval(cs.Fp, u) / (2 * y)
        dbig = (cs.dA(u) - cs.dB(u) * y) - cs.B(u) * dydu
        duds = -N / (dlogN(u) - dbig / big)
        d = np.concatenate([duds] + _rates(
            g, I, [(xend + u) ** a * duds / y for a in range(g)], b[1 + g],
            np.full(n, -float(N)), moments) + [dydu * duds])
        return np.concatenate([d.real, d.imag])

    y0v = np.zeros((nb + 1) * n, complex)
    y0v[:n] = np.atleast_1d(x0) - xend
    y0v[(1 + g) * n:(2 + g) * n] = lam0
    y0v[nb * n:] = np.atleast_1d(y0)
    grid = np.linspace(0.0, XI, track) if track else None
    sol, _ = _solve(rhs, (0.0, XI), np.concatenate([y0v.real, y0v.imag]),
                    method="DOP853", t_eval=grid)
    z = sol.y[:(nb + 1) * n, -1] + 1j * sol.y[(nb + 1) * n:, -1]
    b = [z[i * n:(i + 1) * n] for i in range(nb + 1)]
    k = None
    if moments:
        k = np.stack([np.stack(b[2 + 2 * g + a * g:2 + 2 * g + (a + 1) * g])
                      for a in range(g)])
    leg = Leg(np.stack(b[1:1 + g]), np.stack(b[2 + g:2 + 2 * g]), k)
    if not track:
        return b[0], b[1 + g], leg
    zz = sol.y[:(nb + 1) * n, :] + 1j * sol.y[(nb + 1) * n:, :]
    return b[0], b[1 + g], leg, (sol.t, zz[:n, :], zz[nb * n:, :])


def _weierstrass_chart(cov, seed):
    """Center the y-chart at a root of F, retaining exact rational zeros.

    At a Weierstrass zero of phi, A and B have common factors in x-x_P.
    Centering at a rounded polynomial root replaces those factors by tiny
    nonzero constants; dlog(phi) then has the wrong valuation as y tends
    to zero. Rational branch points are therefore shifted over QQ before
    conversion to float64. Other branch points retain the numerical chart.
    """
    if not hasattr(cov, "_weierstrass_charts"):
        roots = np.roots(cov.F)
        exact = {}
        if cov.exact is not None:
            from sage.all import QQ

            for root in cov.exact[0].roots(QQ, multiplicities=False):
                i = int(np.argmin(np.abs(roots - complex(root))))
                roots[i] = complex(root)
                exact[i] = root
        cov._weierstrass_charts = (roots, exact, {})
    roots, exact, charts = cov._weierstrass_charts
    i = int(np.argmin(np.abs(roots - complex(seed))))
    if i not in charts:
        if i in exact:
            F, a, ad, b, bd = cov.exact
            shifted = F.parent().gen() + exact[i]
            charts[i] = Cover(F(shifted), a(shifted) / ad(shifted),
                              b(shifted) / bd(shifted))
        else:
            charts[i] = cov.shifted(roots[i])
    return roots[i], charts[i]


@traced("endpoint_y", fields=("xend", "yend", "exp"))
def leg_in_y(cov, x0, y0, xend, yend, lam0, moments=False, exp=False):
    r"""
    Straight in `y`; the chart that works at a Weierstrass point.

    `x` is *not* integrated either, for the same reason.  `x_{\mathrm{end}}` is
    a simple root of `F`, so with `u = x - x_{\mathrm{end}}` and
    `F(x_{\mathrm{end}}+u) = u\,G(u)`, `G(0) = F'(x_{\mathrm{end}}) \neq 0`, the
    point on the curve over `y` is the solution of

    .. math:: u = y^2 / G(u),

    which a fixed point iteration finds to full relative accuracy however
    small `y` is.  If it fails to converge the leg is refused and the caller
    shortens it, which makes `y_0` smaller and the iteration easier.
    """
    g = cov.g
    y0, yend = complex(y0), complex(yend)
    xend, cs = _weierstrass_chart(cov, xend)
    G = cs.F[:-1]                       # F(xend + u) = u * G(u)
    dy0 = yend - y0

    def solve_u(y):
        w = y * y
        u = w / np.polyval(G, 0.0)
        for _ in range(60):
            checkpoint("endpoint_iterations", chart="y")
            nxt = w / np.polyval(G, u)
            if abs(nxt - u) <= 1e-15 * abs(nxt):
                return nxt
            u = nxt
        raise LegFailure("the fibre point over y = %r did not converge" % (y,))

    def coord(s):
        if exp:
            return yend + (y0 - yend) * np.exp(-s), -(y0 - yend) * np.exp(-s)
        return y0 + s * dy0, dy0

    def body(s, _u, I, lam):
        y, dy = coord(s)
        u = solve_u(y)
        Fp = np.polyval(cs.Fp, u)
        dx = 2 * y / Fp * dy
        dlam = cs.dlogphi(u, y, dx, dy)
        return (np.atleast_1d(dx),
                [np.atleast_1d(2 * (xend + u) ** a * dy / Fp)
                 for a in range(g)],
                np.atleast_1d(dlam))

    span = (0.0, XI) if exp else (0.0, 1.0)
    return _run(_wrap(g, 1, moments, body), span, g, 1,
                np.atleast_1d(complex(x0) - xend),
                np.atleast_1d(complex(lam0)), moments)


@traced("endpoint_tau", fields=("branch",))
def leg_in_tau(cov, tau0, eta0, lam0, moments=False, branch=1):
    r"""
    Straight in `\tau` down to `\tau = 0`, that is to a point at infinity,
    run exponentially because `\lambda` diverges there.

    In odd degree `x = \tau^{-2}`, `y = \eta\tau^{-(2g+1)}` and

    .. math:: \eta^2 = \tilde F(\tau) = \tau^{4g+2}F(\tau^{-2}), \qquad
              \omega_a = \frac{-2\,\tau^{2g-2a-2}}{\eta}\,d\tau ;

    in even degree infinity is two points, `x = \tau^{-1}`,
    `y = \eta\tau^{-(g+1)}`, and

    .. math:: \eta^2 = \tilde F(\tau) = \tau^{2g+2}F(\tau^{-1}), \qquad
              \omega_a = \frac{-\tau^{g-1-a}}{\eta}\,d\tau ,

    with `\eta(0) = \pm\sqrt{f_{2g+2}}` telling the two apart.  Both are
    regular at `\tau = 0` (indexing `\omega_a = x^a\,dx/y`).

    `\varphi` has a pole of order `p` there, so what is carried is
    `\hat w = w\tau^{p}`, finite and nonzero, and
    `\lambda = \log\hat w - p\log\tau` is linear in the exponential
    parameter.  `p` is `N` when `Q` is the point in hand, but the chart is
    needed whenever a leg runs to infinity, and `p` may then be zero or
    negative -- `\varphi(\infty)` is just another critical value.  In even
    degree it is also where the two points part company, the leading terms of
    `A` and `B y` cancelling at one of them, so it is taken from
    :meth:`~.cover.Cover.pole_order_at_infinity` rather than from the degrees.
    """
    g = cov.g
    pA, pB = cov.tau_powers()
    ra, rad = cov.a[::-1], cov.ad[::-1]
    rb, rbd = cov.b[::-1], cov.bd[::-1]
    rF = cov.F[::-1]
    sq = 1 if cov.even else 2          # F~(tau) is F[::-1] at tau or tau^2

    def at(c, u):
        """``c`` reversed, evaluated at the chart's argument."""
        return np.polyval(c, u if sq == 1 else u * u)

    def ratio(num, den, u):
        """R(tau) and dR/dtau for R = num(tau^sq)/den(tau^sq)."""
        t = u if sq == 1 else u * u
        dt = 1.0 if sq == 1 else 2 * u
        n, d = np.polyval(num, t), np.polyval(den, t)
        dn = dt * np.polyval(np.polyder(num), t)
        dd = dt * np.polyval(np.polyder(den), t)
        return n / d, (dn * d - n * dd) / (d * d)

    F0 = at(rF, tau0)
    omega = (lambda tau, a: -tau ** (g - 1 - a)) if cov.even else (
        lambda tau, a: -2 * tau ** (2 * g - 2 * a - 2))
    pw0 = pA if pB is None else max(pA, pB)      # before any cancellation
    rn, rd = cov.norm_num[::-1], cov.norm_den[::-1]
    kN = sq * (len(cov.norm_den) - len(cov.norm_num))

    def body(s, _e, I, lam):
        tau = tau0 * np.exp(-s)
        dtau = -tau
        eta = eta0 * np.sqrt(at(rF, tau) / F0)                 # not integrated
        dF = (np.polyval(np.polyder(rF), tau) if sq == 1
              else 2 * tau * np.polyval(np.polyder(rF), tau * tau))
        deta = dF / (2 * eta) * dtau
        # both terms scaled by tau^pw0, the order the degrees give.  The
        # true order at this point may be smaller -- that is exactly the
        # cancellation below -- but the scaling has to be the one the lambda
        # rates are written for, or the two disagree by that difference.
        RA, dRA = ratio(ra, rad, tau)
        kA = pw0 - pA
        T1 = tau ** kA * RA
        dT1 = ((kA * tau ** (kA - 1) * RA if kA else 0.0)
               + tau ** kA * dRA) * dtau
        if pB is not None:
            RB, dRB = ratio(rb, rbd, tau)
            kB = pw0 - pB
            T2 = tau ** kB * RB * eta
            dT2 = (((kB * tau ** (kB - 1) * RB * eta if kB else 0.0)
                    + tau ** kB * dRB * eta) * dtau
                   + tau ** kB * RB * deta)
        else:
            T2, dT2 = 0.0 * eta, 0.0 * eta
        # phi tau^pw0 and (phi o iota) tau^pw0.  At a point where the
        # leading terms of A and B y cancel -- which is what makes the pole
        # order there smaller than the degrees say, and happens at one of the
        # two points at infinity in even degree -- the first of these is a
        # difference of two numbers of size tau^(-pw0) and carries no digits
        # at all.  The other one does, and phi = N(phi)/(A - By) recovers
        # phi from it and the norm, which is a ratio of polynomials.
        Tp, dTp = T1 + T2, dT1 + dT2
        Tm, dTm = T1 - T2, dT1 - dT2
        if abs(Tp) >= abs(Tm):
            dlam = dTp / Tp + pw0
        else:
            rat = ratio(rn, rd, tau)
            dlam = (-kN + rat[1] / rat[0] * dtau) - (pw0 + dTm / Tm)
        dI = [np.atleast_1d(omega(tau, a) / eta * dtau) for a in range(g)]
        return np.atleast_1d(deta), dI, np.atleast_1d(dlam)

    lam_hat = np.atleast_1d(complex(lam0))
    return _run(_wrap(g, 1, moments, body), (0.0, XI), g, 1,
                np.atleast_1d(complex(eta0)), lam_hat, moments)


def tau_start(cov, x, y):
    r"""
    `(\tau, \eta)` of a point near infinity, and which point it is.

    Odd degree: `x = \tau^{-2}`, `y = \eta\tau^{-(2g+1)}`, and either branch
    of `\tau = x^{-1/2}` will do -- `(\tau,\eta)` and `(-\tau,-\eta)` are
    the same point of the curve and the integrands are invariant under the
    swap, their exponent being even.  Even degree: `\tau = 1/x`,
    `\eta = y\tau^{-(g+1)}`, with nothing to choose; the sign of
    `\eta(0)` is what says which of the two points at infinity this is.

    Returns `(\tau, \eta, \mathrm{branch})`.
    """
    if cov.even:
        tau = 1.0 / complex(x)
        eta = complex(y) * tau ** (cov.g + 1)
        lead = np.sqrt(complex(cov.F[0]))
        return tau, eta, (1 if abs(eta - lead) <= abs(eta + lead) else -1)
    tau = 1.0 / np.sqrt(complex(x))
    return tau, y * tau ** (2 * cov.g + 1), 1


# ------------------------------------------------------- assembling an edge

def _deflate(c, m, seed):
    r"""
    The `m`-fold root of the polynomial ``c`` nearest ``seed``, from the
    derivative that has it simply, or ``None`` when there is none.

    ``np.roots`` places an `m`-fold root only to within `\varepsilon^{1/m}`:
    a hundredth of the scale for the `N`-fold point `P` of a degree-7 cover.
    A leg aimed there stops that far short, and the piece it drops turns with
    the sheet, by `\zeta_m` from one to the next.  The `(m-1)`-st derivative
    has the same root simply, so ``np.roots`` places *it* to the last digit.
    """
    for _ in range(m - 1):
        c = np.polyder(c)
    if len(c) < 2:
        return None
    r = np.roots(c)
    return r[int(np.argmin(np.abs(r - seed)))]


def ramifies(cov, x, y, xs=None, ys=None, tol=1e-3):
    r"""
    Whether the sheet through `(x,y)` ramifies over the value it lies over.

    Not "is it near a fibre point that does": on the degree-15 cover of LMFDB
    277.a.277.1 an unramified sheet passes within 0.003 of a branch point,
    and aiming its leg there moves it by a hundredth.  Nor "is
    `\mathcal{P}_x` small": that polynomial's coefficients run to `10^6` on
    the same cover, so no threshold on its value means anything.
    `d\varphi/dx` vanishes at the ramification points of the fibre and
    nowhere else on it, and measured against the rest of the fibre the two
    are four orders of magnitude apart.
    """
    if xs is None:
        xs, ys = cov.fibre(complex(cov.phi(np.array([x]), np.array([y]))[0]))
    ok = np.isfinite(xs)
    v = np.abs(cov.y_dphi_dx(xs[ok], ys[ok]))
    scale = float(np.median(v)) if len(v) else 0.0
    if not scale:
        return False
    here = abs(cov.y_dphi_dx(np.array([x]), np.array([y]))[0])
    return bool(here <= tol * scale)


def chart_at(cov, x, y, big=1e6, tol=128*np.finfo(float).eps):
    r"""
    Which of `\tau`, `y`, `x` is a local uniformizer at the point `(x,y)`:
    ``"tau"``, ``"y"`` or ``"x"``.

    The question is about the **point on the curve**, not about its image
    under `\varphi`.  On `y^2 = F(x)` the function `x - x_0` has a simple zero
    at every finite point except a Weierstrass one, where `F(x_0) = 0`, the
    two points over `x_0` have come together, and the zero is double; there
    `y` is the uniformizer instead, `x - x_0 = y^2/G(y)`.  At infinity
    neither is finite and `\tau` is the parameter.  That is the whole rule,
    and it depends on nothing but `(x,y)`.

    So the one thing to decide is whether `F(x_0)` vanishes, and that is
    asked of `F` itself rather than of `|y|`: :func:`~.cover._vanishes`
    measures the value against the size of the terms that made it, which is
    the only scale-free way to ask.  Reading it off `|y| < 10^{-6}` instead
    -- which is what this did -- is an absolute threshold on a quantity with
    no absolute scale. The relative threshold is at the roundoff scale:
    a loose geometric tolerance can mistake a nearby ordinary point for
    a Weierstrass point and move the endpoint of integration.
    """
    if not np.isfinite(x) or abs(x) > big:
        return "tau"
    return "y" if bool(_vanishes(cov.F, x, tol)) else "x"


def w_is_uniformizer(cov, xt, yt, x0, y0, xs=None, ys=None, tol=1e-3):
    r"""
    Whether `w - z_*` is a local uniformizer at `(x_t, y_t)`, the limit of the
    sheet that starts at `(x_0, y_0)` -- that is, whether the sheet stays
    unramified, so the `w`-leg may run to the end instead of handing over to
    a chart.

    Two things are asked, and the first is the one that matters.

    *Against the sheet's own start.*  `d\varphi/dx` collapses along a leg
    exactly when the leg runs into a ramification point, so the ratio of its
    value at the two ends says so directly.  Nothing else does: on
    `\varphi = x` over `y^2 = x(x^2-1)(x^2-4)` *every* point of the fibre
    above a critical value is a ramification point, the critical value is
    known to `2\cdot10^{-15}`, and so `y = \sqrt{F}` comes out as
    `1.1\cdot10^{-7}` rather than zero on both sheets alike.  A test relative
    to that fibre sees a perfectly ordinary point; against the midpoint,
    where `|y\,d\varphi/dx| = 1.19`, the ratio is `9.4\cdot10^{-8}` and the
    ramification is unmistakable.

    *Against the rest of the target fibre*, as :func:`ramifies` does, which
    catches a sheet whose `d\varphi/dx` was already small to begin with.

    This is deliberately not the negation of :func:`ramifies`.  That one asks
    "does this sheet ramify", and where it cannot tell it answers no, its
    caller then refining nothing and being no worse off.  Here the same
    ignorance must answer no as well, and for a sharper reason: a `w`-leg run
    into a ramification point does not come back.  Both are asked positively
    and both fall to the cautious side.
    """
    here = abs(cov.y_dphi_dx(np.array([xt]), np.array([yt]))[0])
    start = abs(cov.y_dphi_dx(np.array([x0]), np.array([y0]))[0])
    if not start or here <= tol * start:
        return False
    if xs is None:
        xs, ys = cov.fibre(complex(cov.phi(np.array([xt]), np.array([yt]))[0]))
    ok = np.isfinite(xs)
    v = np.abs(cov.y_dphi_dx(xs[ok], ys[ok]))
    scale = float(np.median(v)) if len(v) else 0.0
    return bool(scale and here > tol * scale)


def _on_curve(cov, x, yseed):
    r"""`(x, \pm\sqrt{F(x)})`, on the branch of ``yseed``."""
    y = np.sqrt(complex(np.polyval(cov.F, x)))
    return x, (y if abs(y - yseed) <= abs(y + yseed) else -y)


def _target(cov, x0, y0, zstar, kind="simple", slack=1e-2, fib=None):
    r"""
    The fibre point over `z_*` that the sheet through `(x_0,y_0)` runs to.

    The fibre locates a ramification point only as well as ``np.roots``
    locates a multiple root, so a multiple one is refined by :func:`_deflate`:
    on `\operatorname{N}(\varphi)`, which vanishes to order `N` at `x_P`, at
    `P`; on `\mathcal{P}_x(\cdot,z_*)`, which has a simple ramification
    point as a simple root, at an ordinary node.  A refinement is taken only
    when it lands within ``slack`` of where the fibre put it, which is what
    tells a sheet that ramifies over `z_*` from one that merely passes over
    it.
    """
    # Keep endpoint identification in the same centered polynomial used
    # for the initial fibre. Deflating the expanded derivative near a
    # high-order zero can otherwise select a different nearby endpoint.
    if cov.exact is not None:
        if not hasattr(cov, "_centered_fibre"):
            cov.fibre(zstar)
        centered = getattr(cov, "_centered_fibre", None)
        if centered is not None:
            center, local = centered
            local_fib = None if fib is None else (fib[0]-center, fib[1])
            xt, yt = _target(local, x0-center, y0, zstar, kind, slack, local_fib)
            return xt+center, yt
    xs, ys = cov.fibre(zstar) if fib is None else fib
    d = np.abs(xs - x0) + 0.1 * np.abs(ys - y0)
    j = int(np.argmin(d))
    xt, yt = xs[j], ys[j]
    if not np.isfinite(xt):
        return xt, yt
    if kind == "P":
        # every sheet runs into P, so there is nothing to decide here.  The
        # fibre's own N-fold cluster is spread by eps^(1/N) -- a tenth of the
        # scale at N = 15 -- so the cluster test below would not recognise it
        # and the unrefined point would be a tenth out
        zeros = cov.norm_factors()[0]
        if len(zeros) == 1 and zeros[0][1] == cov.N:
            # The exact squarefree norm knows this point. Deflating the
            # rounded N-fold polynomial can displace it by eps^(1/N).
            return _on_curve(cov, zeros[0][0], yt)
        r = _deflate(cov.norm_num, cov.N, xt)
        return (xt, yt) if r is None else _on_curve(cov, r, yt)
    # In degree three the two ramifying sheets dominate the median used
    # by ramifies(). Also compare with this sheet before the endpoint:
    # dphi/dx collapses along a ramifying lift, while its starting value is
    # nonzero. The derivative-root proximity check below remains necessary.
    here = abs(cov.y_dphi_dx(np.array([xt]), np.array([yt]))[0])
    start = abs(cov.y_dphi_dx(np.array([x0]), np.array([y0]))[0])
    approaching_ramification = start > 0 and here <= 1e-3 * start
    if not (ramifies(cov, xt, yt, xs, ys) or approaching_ramification):
        return xt, yt                  # this sheet does not ramify over z_*
    r = _deflate(np.polyder(cov.poly_in_x(zstar)), 1, xt)
    if r is None or abs(r - xt) > slack * (1.0 + abs(xt)):
        return xt, yt
    if cov.exact is not None and not cov.pure_x:
        # A double fibre root is sensitive to the rounded critical value.
        # Instead solve the exact critical-point equation (phi_x = 0),
        # shifted before float conversion so its simple root is well scaled.
        if not hasattr(cov, "_critical_x_polynomial"):
            F, a, ad, b, bd = cov.exact
            A, B = a/ad, b/bd
            q = (F*A.derivative()**2
                 - (F*B.derivative()+F.derivative()*B/2)**2).numerator()
            cov._critical_x_polynomial = q // q.gcd(q.derivative())
        from sage.all import QQ
        q = cov._critical_x_polynomial
        center = QQ(complex(r).real)
        shifted = q(q.parent().gen()+center)
        fine = _deflate(np.asarray([complex(c) for c in shifted.list()[::-1]]),
                        1, r-complex(center))
        if fine is not None:
            fine += complex(center)
            if abs(fine-r) <= slack*(1+abs(r)):
                r = fine
    return _on_curve(cov, r, yt)


@traced("end_leg", fields=("zstar", "kind", "at_infinity"))
def end_leg(cov, w0, x0, y0, lam0, zstar, kind, moments, big=1e6,
            at_infinity=False, fib=None):
    r"""
    One sheet's leg from `(x_0,y_0)` over `w_0` in to its limit over `z_*`.

    `kind` is ``"simple"`` when `\lambda` stays finite there, ``"P"`` when
    `\varphi \to 0` and ``"Q"`` when `\varphi \to \infty`; the last two run
    exponentially.

    A chart is needed only where `w` stops being a local coordinate, and that
    is exactly where the sheet *ramifies*.  Over an ordinary point of the
    fibre `w - z_*` is itself a uniformizer, `d\varphi/dx \neq 0`, and the
    `w`-leg simply runs on to `z_*` -- no chart, no stand-off, no hand-over,
    and measured against the chart route it agrees to `10^{-14}`.  At a
    simple critical value that is `N-2` of the `N` sheets, so the charts are
    for the two that do ramify.  Which they are is not a matter of distance:
    :func:`ramifies` asks whether `d\varphi/dx` vanishes there.

    For the sheets that do ramify the chart is the local uniformizer at the
    limit *point*, which :func:`chart_at` reads off `(x_t, y_t)` alone.
    """
    exp = kind in ("P", "Q")
    if kind == "Q":
        # phi -> infinity, so this sheet runs to a pole of phi, and which
        # point that is belongs to the *cover*, not to the sheet: in even
        # degree infinity is two points, phi has its pole at just one of
        # them, and every sheet runs to that one.
        #
        # Reading it off the sheet's own eta -- which is what this did --
        # cannot work.  phi has a pole of order N there, so phi ~ x^N and
        # |x| ~ |w|^(1/N): on the degree-6 model of LMFDB 249.a.249.1, with
        # N = 14, the fibre is still at |x| ~ 4 when |w| = 5e8, and no reach
        # brings eta near +-sqrt(f_{2g+2}).  Two of the fourteen sheets then
        # named the point where phi is *finite*, whereupon this looked for a
        # finite pole, found ad*bd constant, and declared that phi has no
        # pole at all.
        pole = next((b for b in ((1, -1) if cov.even else (1,))
                     if cov.pole_order_at_infinity(b) == cov.N), None)
        if pole is not None:
            tau, eta, _ = tau_start(cov, x0, y0)
            return leg_in_tau(cov, tau, eta, lam0, moments, branch=pole)[2]
        den = np.polymul(cov.ad, cov.bd)
        if len(den) < 2:
            raise LegFailure("phi has no pole: it cannot be N P - N Q")
        r = np.roots(den)
        xt = r[int(np.argmin(np.abs(r - x0)))]
        # 1/phi vanishes to order N at a finite Q, so N(phi)'s denominator
        # does too, and the same deflation places it to the last digit
        fine = _deflate(cov.norm_den, cov.N, xt)
        if fine is not None and abs(fine - xt) <= 1e-2 * (1.0 + abs(xt)):
            xt = fine
        xt, yt = _on_curve(cov, xt, y0)
    elif at_infinity:
        tau, eta, br = tau_start(cov, x0, y0)
        return leg_in_tau(cov, tau, eta, lam0, moments, branch=br)[2]
    else:
        xt, yt = _target(cov, x0, y0, zstar, kind, fib=fib)
    # a sheet running to infinity comes back as a huge root rather than as one
    # at infinity: the cover polynomial's leading coefficient vanishes there
    # in exact arithmetic but not in float64, so ``np.roots`` returns 1e16
    # instead of dropping the root.  Either way the tau-chart is what reaches
    # it -- in even degree, at whichever of the two points at infinity the
    # sheet's own eta says.
    chart = chart_at(cov, xt, yt, big)
    checkpoint(chart=chart, target_x=xt, target_y=yt)
    if chart == "tau":
        tau, eta, br = tau_start(cov, x0, y0)
        return leg_in_tau(cov, tau, eta, lam0, moments, branch=br)[2]
    if chart == "y":
        # leg_in_y centers at a root of F and retains exact rational roots.
        # Replacing an exact endpoint by np.roots(F) here destroys the
        # common zeros of A and B and corrupts the logarithmic integral.
        return leg_in_y(cov, x0, y0, xt, 0.0, lam0, moments, exp)[2]
    # An unramified sheet could be run straight in with the w-leg, w being a
    # uniformizer at its limit point (:func:`w_is_uniformizer`), and that is
    # what the ODE did.  The certified tracker cannot: its step is bounded by
    # the distance to the nearest value where *any* two roots collide, and
    # z_* is one, so rho is zero there whatever sheet is being followed.  The
    # short-cut measured as neutral -- same defect to 1e-15, same time -- so
    # it goes rather than the guarantee.
    _check_clear(cov, x0, xt)
    if not exp and abs(yt) > 0:
        # A short straight x-tail may still end on the wrong y-sheet:
        # the true lift can turn at a Weierstrass point before reaching
        # its critical endpoint. Continue sqrt(F) along this segment by
        # its linear root factors; each factor has argument change < pi.
        # If its sign disagrees with the target, move the handover closer.
        if not hasattr(cov, "_tail_branch_points"):
            cov._tail_branch_points = np.roots(cov.F)
        roots = cov._tail_branch_points
        continued_y = y0*np.prod(np.sqrt((xt-roots)/(x0-roots)))
        if abs(continued_y-yt) > abs(continued_y+yt):
            raise LegFailure("straight x-tail reaches the opposite y-sheet")
    return leg_in_x(cov, x0, y0, xt, lam0, moments, exp)[2]


def _check_clear(cov, x0, xt, room=0.25):
    r"""
    Refuse an `x`-chart leg that would pass a zero of `F`.

    The leg is a straight path in `x` and `y` follows it; that is only the
    continuation the integral wants while the path stays on one side of every
    branch point of `y`.  At the far node of the degree-14 cover on LMFDB
    249.a.249.1 the ramification point has `F(x) = 0.003`, a zero of `F`
    three thousandths away, and a leg starting a hundredth out winds around
    it: two per cent of the answer, on exactly the two sheets that meet
    there.  Raising here puts the caller's stand-off loop to work, which
    shortens the leg until it is clear -- the stand-off that :func:`side`
    would otherwise choose comes from the spacing of the *critical values*,
    and says nothing about how close a Weierstrass point happens to be.
    """
    r = np.roots(cov.F) if len(cov.F) > 1 else np.zeros(0, complex)
    if not len(r):
        return
    d = np.abs(r - xt)
    d = d[d > 1e-10 * max(1.0, abs(xt))]      # xt itself is the y-chart's job
    if len(d) and abs(x0 - xt) > room * d.min():
        raise LegFailure("the x-chart leg would pass a zero of F: %.2e away "
                         "against a reach of %.2e" % (d.min(), abs(x0 - xt)))


def _delta_data(cov, chart=None):
    r"""
    `\mathcal{P}` read as a polynomial in the fibre's coordinate, with the
    values of `w` over which its roots collide or escape.

    ``chart`` is ``None`` for `x` itself, or a constant `c` for
    `x' = 1/(x-c)`, in which case the coefficients are `\mathcal{P}`'s
    Taylor-shifted by `c` and reversed.

    A chart is needed because the bound's `M` is a bound on the *roots*,
    `\max_k(|a_k|/|a_N|)^{1/(k+1)}`, so it diverges wherever one escapes --
    and on an even degree model one does: infinity is two points, `\varphi`
    has its pole at one of them and a finite value at the other, so that
    value carries a preimage at infinity and `a_N` vanishes there.  Nor is
    `1/x` always the answer: on the degree-6 model of LMFDB 249.a.249.1 the
    node `w = 2` has `a_N(2) = 0` *and* `a_0(2) = 0`, a point at infinity
    and a point at `x = 0`, so `x` and `1/x` fail in opposite ways and only
    a `c` clear of the whole fibre works.
    """
    key = None if chart is None else complex(chart)
    if getattr(cov, "_dd", None) is None or not isinstance(cov._dd, dict):
        cov._dd = {}
    if key not in cov._dd:
        n = cov.N
        ps = (cov.p if chart is None
              else [_shift_poly(q, complex(chart)) for q in cov.p])
        aks = []
        for k in range(n + 1):
            col = [q[len(q) - 1 - k] if k < len(q) else 0.0 for q in ps]
            aks.append(np.array(col[::-1], complex))     # a_k(w), descending
        if chart is not None:
            aks = aks[::-1]
        top = np.trim_zeros(np.asarray(aks[n], complex), "f")
        if not len(top):
            top = np.zeros(1, complex)
        roots = np.roots(top) if len(top) > 1 else np.zeros(0, complex)
        bad = np.concatenate([np.asarray(cov.critical_values(), complex),
                              np.asarray(roots, complex)])
        cov._dd[key] = (aks, top, roots, bad)
    return cov._dd[key]


def _nice(z):
    """A nearby rational, so the certificate is built over QQ and cached."""
    from sage.all import QQ

    return QQ(int(round(z.real))) if abs(z.imag) < 1.0 else QQ(int(round(z.real)))


def _delta_data_s(cov, z, e, chart=None):
    r"""
    :func:`_delta_data` in the local uniformizer `s`, where `w = z + s^e`.

    Substituting into `\mathcal{P}` leaves a polynomial in `x` of the same
    degree whose coefficients are polynomials in `s`: writing
    `a_k(w) = Aw^2 + Bw + C`, the coefficient of `x^k` becomes

    .. math:: A\,s^{2e} + (2Az + B)\,s^{e} + (Az^2 + Bz + C).

    Its roots degenerate where `w` does, so the bad values of `s` are the
    `e`-th roots of `w_j - z` over every bad `w_j` -- which puts `s = 0`
    among them, `z` itself being one.  Bounding the step in `s` is not the
    same as bounding it in `w` and dividing by `|dw/ds|`: that rate collapses
    like `s^{e-1}` at the end of the leg, so the converted step runs away
    exactly where it must not.
    """
    key = (complex(z), int(e), None if chart is None else complex(chart))
    if getattr(cov, "_dds", None) is None:
        cov._dds = {}
    if key not in cov._dds:
        aks, _, _, badw = _delta_data(cov, chart)
        zc = complex(z)
        out = []
        for a in aks:
            A, B, C = (complex(a[0]) if len(a) > 2 else 0.0,
                       complex(a[-2]) if len(a) > 1 else 0.0, complex(a[-1]))
            col = np.zeros(2 * e + 1, complex)
            col[0] = A                                  # s^{2e}
            col[e] = 2.0 * A * zc + B                   # s^{e}
            col[2 * e] = A * zc * zc + B * zc + C       # s^0
            out.append(np.trim_zeros(col, "f"))
        top = out[cov.N]
        if not len(top):
            top = np.zeros(1, complex)
        roots = np.roots(top) if len(top) > 1 else np.zeros(0, complex)
        bad = []
        for wj in np.concatenate([badw, np.asarray(roots, complex) * 0.0]):
            r = complex(wj) - zc
            if r == 0:
                bad.append(0.0 + 0.0j)
                continue
            rad, th = abs(r) ** (1.0 / e), np.angle(r) / e
            bad += [rad * np.exp(1j * (th + 2.0 * np.pi * j / e))
                    for j in range(e)]
        bad += list(roots)
        cov._dds[key] = (out, top, roots, np.array(bad, complex))
    return cov._dds[key]


def certified_delta(cov, w, xs, eps, unif=None, chart=None):
    r"""
    How far the fibre may be moved in one step, from Sage's bound
    (:meth:`RiemannSurface._compute_delta`) written against `\mathcal{P}`.

    With ``unif = (z, e)`` the bound is taken in the local uniformizer
    `s`, `w = z + s^e`, and returned as a step in `s`; ``w`` is then the
    value of `s`.  Without it the bound and the step are in `w`.

    With `\rho` half the distance to the nearest value where the roots
    collide or escape, `Y` the largest `|\partial_s/\partial_x|` on the
    fibre, and `M` a bound on the roots read off the coefficients,

    .. math:: \delta = \rho\,
              \frac{\sqrt{(\rho Y - \varepsilon)^2 + 4\varepsilon M}
                    - (\rho Y + \varepsilon)}{2M - 2\rho Y},

    small enough that Newton from the previous fibre cannot leave the
    `\varepsilon` ball around it, `\varepsilon` being the least distance
    between the fibre's points divided by three.
    """
    if unif is None:
        aks, top, a0roots, bad = _delta_data(cov, chart)
        wv, scale = complex(w), 1.0
    else:
        zc, e = complex(unif[0]), int(unif[1])
        aks, top, a0roots, bad = _delta_data_s(cov, zc, e, chart)
        wv = complex(w)
        # |ds/dx| against |dw/dx|: dw = e s^{e-1} ds
        scale = e * abs(wv) ** (e - 1) if e > 1 else 1.0
    if not len(bad):
        return None
    rho = float(np.abs(bad - wv).min()) / 2.0
    if not rho:
        return None
    if not np.isfinite(eps):
        # a single distinct root: there is no other for Newton to be
        # confused with, so nothing bounds the step but the structure
        # itself.  This is Sage's own fallback, half the distance to the
        # branch locus, and it is what the whole formula degenerates to.
        return rho
    wq = complex(unif[0]) + wv ** int(unif[1]) if unif is not None else wv
    px, pw = cov.Px(xs, wq), cov.Pw(xs, wq)
    ok = np.abs(px) > 0.0
    if not np.any(ok):
        return None
    # dtau/dw = -tau^2 dx/dw, so the rate in the reciprocal chart carries
    # a factor 1/x^2 -- the whole point being that it is *small* where x is
    # large, which is what keeps the bound finite as a root escapes
    rate = (1.0 if chart is None else
            1.0 / (np.asarray(xs, complex)[ok] - complex(chart)) ** 2)
    Y = float((np.abs(pw[ok] / px[ok]) * np.abs(rate)).max()) * scale
    r = abs(wv) + rho
    upper = [float(np.polyval(np.abs(a), r)) for a in aks[:-1]]
    if len(a0roots):
        low = abs(complex(top[0])) * float(
            np.prod(np.abs(a0roots - wv) - rho)) / 2.0
    else:
        low = abs(complex(top[-1])) / 2.0
    if not np.isfinite(low) or low <= 0.0:
        return None
    M = 2.0 * max((upper[k] / low) ** (1.0 / (k + 1))
                  for k in range(len(upper)))
    den, root = 2.0 * M - 2.0 * rho * Y, (rho * Y - eps) ** 2 + 4.0 * eps * M
    if den <= 0.0 or root < 0.0:
        return None
    out = rho * (np.sqrt(root) - (rho * Y + eps)) / den
    return float(out) if np.isfinite(out) and out > 0.0 else None


@traced("certified_reach", fields=("m", "zend", "cap"))
def certified_reach(cov, m, zend, cap=0.08, floor=1e-9, tries=400):
    r"""
    The stand-off to approach `z_{\mathrm{end}}` with, as a fraction of
    `|m - z_{\mathrm{end}}|`, from the certified step.

    Stepping by :func:`certified_delta` and re-solving the fibre each time, no
    sheet can be exchanged for another: the step is bounded so that every
    point stays inside its own `\varepsilon` ball, `\varepsilon` being the
    least distance between the fibre's points divided by three.  Where the
    stepping stalls is where the sheets stop being certifiably apart, and
    *that* is how close a leg may be taken before a chart has to finish the
    job -- not eight per cent of an edge, and not a quarter of the way to the
    nearest other node, neither of which knows anything about the fibre.

    Capped at ``cap``, so this can only bring the hand-over in, never push it
    out past what the caller already trusts.
    """
    span = abs(complex(zend) - complex(m))
    if not span or cov.exact is None:
        return None
    t = 0.0
    xs, _ = cov.fibre(complex(m))
    for _ in range(tries):
        checkpoint("certified_attempts", t=t)
        w = complex(m) + t * (complex(zend) - complex(m))
        fin = np.isfinite(xs)
        if int(fin.sum()) < 2:
            break
        d = np.abs(xs[fin][:, None] - xs[fin][None, :])
        np.fill_diagonal(d, np.inf)
        step = certified_delta(cov, w, xs[fin], float(d.min()) / 3.0)
        if step is None or step / span < floor:
            break
        t = min(1.0, t + step / span)
        if t >= 1.0:
            break
        xs, _ = cov.fibre(complex(m) + t * (complex(zend) - complex(m)))
    return max(min(1.0 - t, cap), floor) if t > 0.0 else None


@traced("side", fields=("m", "zend", "kind"))
def side(cov, m, x, y, lam, zend, kind, moments, frac=0.08, hold=None):
    r"""
    Half an edge, from the reference point `m` out to `z_{\mathrm{end}}`.

    The chart legs are only valid near their target, so the stand-off is pulled
    in until every sheet's leg lands, which also keeps the straight chart path
    short enough to stay on its branch.  Continuation first tests the
    proposed hand-over and its endpoint legs; only a usable hand-over is
    integrated, so rejected stand-offs do not repeat the quadrature.

    ``hold`` caps how far from `z_{\mathrm{end}}` the chart leg may start, in
    absolute terms.  A fraction of the edge is not a cap at all when the edges
    are of wildly different lengths: on the degree-15 cover of LMFDB
    277.a.277.1 one edge is `8\cdot10^4` long, so eight per cent of it starts
    the chart leg three thousand units from its target and the leg is simply
    wrong -- silently, since nothing about it fails.  Entered from four units
    away the same leg is right to `10^{-11}`.
    """
    if kind == "none":
        out = leg_in_w(cov, m, zend, x, y, lam, moments)[3]
        out.frac = 1.0
        return out
    last, certified = None, None
    if hold is not None and abs(m - zend) > 0:
        frac = min(frac, hold / abs(m - zend))
    # A relative cutoff rejects long edges before any attempt when hold
    # already makes frac tiny. Bound retries, and stop only when the actual
    # target displacement is no longer representable.
    for attempt in range(16):
        checkpoint("standoff_attempts", fraction=frac)
        w0 = zend + frac * (m - zend)
        if w0 == zend or w0 == m:
            last = LegFailure("endpoint stand-off is not representable")
            break
        try:
            # Numerical critical values must not be rationalized as if they
            # were exact fibres: that turns a ramified endpoint into e=1.
            e = cov.endpoint_ramification(zend, kind)
            checkpoint(ramification=e)
            x0, y0, lam0, L1 = leg_in_w(cov, m, w0, x, y, lam, moments,
                                        toward=(zend, e), integrate=moments)
            if kind == "P":
                xt, yt = _target(cov, x0[0], y0[0], zend, "P")
                if np.isfinite(xt) and abs(yt) > 1e-6 * max(1.0, abs(xt)):
                    # the whole fibre slides into P together
                    tail = leg_into_P(cov, x0, y0, xt, lam0, moments)[2]
                    if not moments:
                        L1 = leg_in_w(cov, m, w0, x, y, lam,
                                      toward=(zend, e), _prepared=L1.charts)[3]
                    out = L1 + tail
                    out.frac = frac
                    return out
            # which sheets run to a point at infinity is a matter of
            # counting, not of distance: the fibre over z_* says how many of
            # its points are there, and those are the sheets whose x has run
            # furthest by the stand-off.  Distance cannot answer it -- every
            # finite fibre point is nearer to them than infinity is -- and
            # magnitude alone cannot either, since the stand-off may leave
            # them at |x| = 6 while a finite point sits at 1.6.
            fib = cov.fibre(complex(zend))
            far = int(np.count_nonzero(~np.isfinite(fib[0])))
            gone = set(np.argsort(-np.abs(x0))[:far].tolist()) if far else set()
            legs = []
            for i in range(len(x0)):
                with phase("endpoint_tail", sheet=i):
                    legs.append(end_leg(cov, w0, x0[i], y0[i], lam0[i], zend,
                                        kind, moments, at_infinity=(i in gone), fib=fib))
            if not moments:
                L1 = leg_in_w(cov, m, w0, x, y, lam, toward=(zend, e),
                              _prepared=L1.charts)[3]
            out = L1 + Leg.stack([Leg(L.dI[:, 0], L.dL[:, 0],
                                       None if L.k is None else L.k[:, :, 0])
                                   for L in legs])
            out.frac = frac
            return out
        except (LegFailure, FloatingPointError) as exc:
            checkpoint("failed_standoffs", failure=str(exc), fraction=frac)
            last = exc
            if certified is None:
                # the fixed fraction did not work; ask the certified tracker
                # where the sheets actually stop being apart, rather than
                # walking in by thirds and hoping
                certified = certified_reach(cov, m, zend, cap=frac) or False
                if certified:
                    frac = certified
                    continue
            frac /= 3.0
    raise RuntimeError("no usable stand-off approaching %r (%s)" % (zend, last))


@traced("finite", fields=("za", "zb", "kind_a", "kind_b", "moments"), is_edge=True)
def edge_leg(cov, za, zb, lam_ref, kind_a, kind_b, moments=False,
             hold_a=None, hold_b=None):
    r"""
    The lifts of `[z_a, z_b]`, oriented `z_a \to z_b`, indexed by the fibre at
    the midpoint, which is returned alongside.

    `lam_ref` is `\log\varphi` at the midpoint on the side of the cut in hand.
    """
    m = (za + zb) / 2
    x, y = cov.fibre(m)
    lam = np.full(len(x), complex(lam_ref))
    return (x, y), (side(cov, m, x, y, lam, za, kind_a, moments,
                         hold=hold_a).reversed()
                    + side(cov, m, x, y, lam, zb, kind_b, moments, hold=hold_b))


def _odd_infinity_region(cov):
    """An exterior disk containing no finite branch point or zero of phi.

    Coefficient bounds give |F/(lc(F)*x^deg(F))-1| <= 1/4 and the
    same bound for Norm(phi) throughout |x| >= radius. Thus the infinity
    square root has one analytic branch in tau=1/sqrt(x), and phi has
    no zero there. Polynomial A,B have no finite poles.
    """
    radius = 1.0
    for _ in range(64):
        if all(sum(abs(c[j]/c[0])*radius**(-j)
                   for j in range(1, len(c))) <= 0.25
               for c in (cov.F, cov.norm_num)):
            break
        radius *= 2.0
    else:
        raise LegFailure("could not bound the odd infinity chart")
    # An upper bound for |phi| on |x|=radius, on either sheet.
    bound = (np.polyval(abs(cov.a), radius)/abs(cov.ad[0])
             + np.polyval(abs(cov.b), radius)/abs(cov.bd[0])
             * np.sqrt(np.polyval(abs(cov.F), radius)))
    if not np.isfinite(bound):
        raise LegFailure("odd infinity chart bound overflowed")
    return radius, float(bound)


@traced("odd_infinity_ray")
def _odd_infinity_ray(cov, zr, d, span, x, y, lam, level):
    """Track a finite ray segment, then integrate in the true uniformizer.

    Once |w| exceeds the boundary bound, a lifted outward ray cannot cross
    |x|=radius. Its tail and the straight tau tail lie in the same infinity
    chart. This avoids separating sheets in 1/x, which is ramified at an
    odd-degree infinity and makes their x-projections nearly coincide.
    """
    radius, bound = _odd_infinity_region(cov)
    distance = span * abs(d)
    # For every later point |zr + length*d/|d|| >= length-|zr| > bound.
    # All poles are at infinity, so every lift eventually lies outside the
    # disk; none can cross its boundary on this whole final ray.
    factor = max(1.0, (2*bound + 2*abs(zr) + distance)/distance)
    log_factor = np.log(factor)
    sigma0 = 1.0/(span*complex(d))
    coord = lambda t: sigma0*np.exp(-log_factor*t)
    rate = lambda t: -log_factor*coord(t)
    wof = lambda t: complex(zr)+1.0/coord(t)
    checkpoint(infinity_radius=radius, boundary_phi_bound=bound,
               ray_length_factor=factor)
    xf, yf, lf, I, L = _walk_quadrature(
        cov, coord, rate, x, y, lam,
        walks.charts_over(cov, None, sigma=True, zr=complex(zr)),
        level, wof=wof)
    if np.any(abs(xf) <= radius):
        raise LegFailure("ray handover did not enter the infinity chart")
    tails = []
    for i in range(len(xf)):
        tau, eta, branch = tau_start(cov, xf[i], yf[i])
        with phase("endpoint_tail", sheet=i):
            tail = leg_in_tau(cov, tau, eta, lf[i], branch=branch)[2]
        tails.append(Leg(tail.dI[:, 0], tail.dL[:, 0], None))
    return Leg(I, L, None) + Leg.stack(tails)


@traced("ray", fields=("zr", "d", "span", "moments"), is_edge=True)
def ray_leg(cov, zr, d, lam_ref, moments=False, reach=40.0,
            kind_in="simple", span=1.0, hold=None, level=5):
    r"""
    The lifts of the ray `z_r \to \infty`, indexed by the fibre at
    `z_r + \mathrm{span}\,d`.

    The inward piece is an ordinary edge.  The outward piece runs to `Q`
    in `\sigma = 1/(w - z_r)`: the ray is then the straight segment
    from `1/(\mathrm{span}\,d)` to `\sigma = 0`, and `\sigma = 0` is `Q`.
    On odd-degree polynomial covers, the outward ray is tracked only to
    an exterior region bounded by the coefficients of F and Norm(phi),
    then finished in the true infinity uniformizer tau=1/sqrt(x). Tracking
    1/x all the way to zero would artificially coalesce opposite tau sheets.
    The remaining covers use uniformizer quadrature with logarithmic tracking
    so every step remains resolvable near the endpoint.  ``span`` is the size of the cut system
    (:func:`~.stokes.cut_scale`), and the reference point must be the one the
    racetrack labelled the sheets at.

    What this replaces was a straight `w`-leg out to ``reach`` spans followed
    by a chart at `Q` per sheet.  Neither half of that was satisfactory.  The
    reach is a guess -- `\varphi` has a pole of order `N`, so
    `|x| \sim |w|^{1/N}`, and on the degree-14 cover of LMFDB 249.a.249.1 the
    fibre is still at `|x| \approx 4` when `|w| = 5\cdot10^{8}`, which is to
    say that no reach arrives.  The chart then had to be told *which* point
    at infinity each sheet was running to, a question the sheet cannot answer
    from that distance, and the ray's sheets summed to `1.8` where a full
    preimage must sum to zero.  In `\sigma` there is nothing to guess: the
    endpoint is a point, and the quadrature approaches it in a regular
    coordinate raced for like any other.

    ``reach`` is accepted and ignored; ``moments`` keeps the old route, a
    panel rule having no running primitive to offer.
    """
    ref = zr + span * d
    x, y = cov.fibre(ref)
    lam = np.full(len(x), complex(lam_ref))
    inward = side(cov, ref, x, y, lam, zr, kind_in, moments, hold=hold)
    if moments:
        far = zr + reach * span * d
        x1, y1, lam1, out1 = _leg_in_w_ode(cov, ref, far, x, y, lam, moments)
        tails = []
        for i in range(len(x1)):
            with phase("endpoint_tail", sheet=i):
                tails.append(end_leg(cov, far, x1[i], y1[i], lam1[i], np.inf,
                                     "Q", moments))
        out = out1 + Leg.stack([Leg(L.dI[:, 0], L.dL[:, 0],
                                    None if L.k is None else L.k[:, :, 0])
                                for L in tails])
        return (x, y), inward.reversed() + out

    if not cov.even and len(cov.ad) == len(cov.bd) == 1:
        out = _odd_infinity_ray(cov, zr, d, span, x, y, lam, level)
        return (x, y), inward.reversed() + out

    # sigma = s^e with s linear, e the ramification of phi at Q: the
    # uniformizer substitution that makes omega_q/ds regular there, exactly
    # as the legs into P and into the branch points make it regular at theirs
    pls = cov.places_over(None)
    e = max(int(q["e"]) for q in pls) if pls else cov.N
    sa = (1.0 / (span * complex(d))) ** (1.0 / e)
    # Keep the same quadrature in the uniformizer s = sa*(1-t), but
    # track in r = -log(1-t).  Near t=1 a single floating-point ulp can
    # move a sheet by half its separation; halving that step cannot resolve
    # it.  In r, relative changes of s remain representable all the way to
    # the last quadrature node.  The endpoint itself is never sampled.
    coord = lambda r: (sa * np.exp(-r)) ** e
    rate = lambda r: -e * coord(r)
    wof = lambda r: complex(zr) + 1.0 / coord(r)
    try:
        _, _, _, I, L = _walk_quadrature(
            cov, coord, rate, x, y, lam,
            walks.charts_over(cov, None, sigma=True, zr=complex(zr)), level,
            wof=wof, logarithmic=True)
    except RuntimeError as exc:
        raise LegFailure(str(exc))
    return (x, y), inward.reversed() + Leg(I, L, None)
