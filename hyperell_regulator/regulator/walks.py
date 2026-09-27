r"""Adaptive Newton continuation of preimages of phi.

Each sheet is tracked in x-c or 1/(x-c). Endpoint place data guides the
choice: on an even-degree model the inverse coordinate is a uniformizer at
infinity. At a finite Weierstrass point, and at odd-degree infinity, these
x-coordinates are ramified; the endpoint integrator in :mod:`._legs` uses
y or a local infinity parameter there.

Newton corrections and step displacements are checked against numerical
root separation. These are adaptive numerical checks, not interval
certificates; monodromy, quadrature convergence and boundary closure provide
additional independent checks.

For inverse coordinates, eta=u^(g+1)y and the holomorphic differentials are
-((cu+1)^q u^(g-1-q)/eta) du. On the curve, the cover polynomial satisfies
P_w=-2*b*bd*ad^2*y. Cancelling eta in this identity before evaluating the
inverse differential removes its apparent singularity at a Weierstrass
crossing while retaining the homogeneous scaling needed near infinity.
"""

import numpy as np

from .polynomials import evaluate as _eval

from .diagnostics import checkpoint
from .cover import _shift_poly

__all__ = ["Chart", "Walker", "charts_over", "walker", "race",
           "choose", "leg"]


class Chart:
    r"""
    A Moebius coordinate `u` on the `x`-line, with the cover polynomial,
    `F`, and the holomorphic differentials all read in it.

    ``kind`` is ``"x"`` for `u = x - c` and ``"inv"`` for `u = 1/(x - c)`.
    The second is what makes a point at infinity ordinary: it sits at
    `u = 0`, and the cover polynomial's coefficients are simply reversed.
    """

    __slots__ = ("kind", "c", "cw", "dcw", "Fu", "g", "n", "sigma", "zr",
                 "cov", "local", "sheet_a", "sheet_b", "sheet_d", "_fibres",
                 "sheet_t", "sheet_t_power", "branch_points")

    def __init__(self, cov, kind, c, sigma=False, zr=0.0):
        self.cov, self.kind, self.c, self.g = cov, kind, complex(c), cov.g
        self.sigma, self.zr = bool(sigma), complex(zr)
        self.local = cov.shifted(self.c) if kind == "x" else cov.reciprocal(self.c)
        self._fibres = {}
        # B == 0 makes P a perfect square -- both points over an x are in the
        # fibre and x does not separate them -- and Newton on a double root
        # converges linearly to half the digits.  Its square root
        # a - ad*w has the same roots, simple, and the sign of eta carries
        # the two sheets apart, exactly as Cover.fibre does it.
        pw = ([np.asarray(cov.a, complex), -np.asarray(cov.ad, complex)]
              if cov.pure_x else cov.p)
        if kind == "x":
            ps = ([np.asarray(self.local.a, complex), -np.asarray(self.local.ad, complex)]
                  if cov.pure_x else self.local.p)
        else:
            ps = [np.asarray(_shift_poly(q, self.c), complex) for q in pw]
        self.n = max(len(q) for q in ps) - 1
        # a_k(w), the coefficient of (x - c)^k, as a polynomial in w
        aks = []
        for k in range(self.n + 1):
            col = [q[len(q) - 1 - k] if k < len(q) else 0.0 for q in ps]
            aks.append(np.array(col[::-1], complex))
        if kind == "inv":
            aks = aks[::-1]
        # cw[k, j] is the coefficient of u^k w^j; every evaluation below is
        # then one matrix-vector product rather than a loop over the sheets'
        # worth of np.polyval calls
        d = max(len(a) for a in aks)
        self.cw = np.zeros((self.n + 1, d), complex)
        for k, a in enumerate(aks):
            self.cw[k, :len(a)] = a[::-1]
        # A Moebius coordinate on the w-line as well: sigma = 1/(w - z_r)
        # turns sigma^D P(x, z_r + 1/sigma) into a polynomial in sigma, which
        # is the shifted coefficients read backwards -- the same trick the
        # chart at x = infinity plays on the other variable.  It is what
        # makes w = infinity, which is Q, an ordinary endpoint, and the shift
        # by z_r is what keeps the straight path in sigma the image of the
        # *ray* rather than of some other path out of the configuration.
        if self.sigma:
            self.cw = np.array([_shift_poly(r[::-1], self.zr)
                                for r in self.cw], complex)
        self.dcw = self.cw[:, 1:] * np.arange(1, d)
        # eta^2 = Fu(u): F(c+u) in the first chart, u^{2g+2}F(c+1/u) in the
        # second, which is that polynomial's coefficients padded to degree
        # 2g+2 and reversed
        Fs = np.asarray(_shift_poly(cov.F, self.c), complex)
        if kind == "x":
            self.Fu = self.local.F
        else:
            asc = np.zeros(2 * cov.g + 3, complex)
            asc[:len(Fs)] = Fs[::-1]
            self.Fu = asc

        # A + B*y = w also identifies the sign of eta.  Continuity of
        # sqrt(F) alone fails when an unramified walk crosses y=0: x turns
        # round there and the correct square root changes sign.
        polys = [np.polymul(cov.a, cov.bd), np.polymul(cov.b, cov.ad),
                 np.polymul(cov.ad, cov.bd)]
        if kind == "x":
            local = self.local
            polys = [np.polymul(local.a, local.bd),
                     np.polymul(local.b, local.ad),
                     np.polymul(local.ad, local.bd)]
        else:
            polys = [np.asarray(_shift_poly(q, self.c), complex) for q in polys]
        if kind == "inv":
            degree = max(len(polys[0]) - 1,
                         len(polys[1]) - 1 + self.g + 1,
                         len(polys[2]) - 1)
            polys = [np.pad(q[::-1], (0, degree - offset - len(q) + 1))
                     for q, offset in zip(polys, (0, self.g + 1, 0))]
        self.sheet_a, self.sheet_b, self.sheet_d = polys
        self.branch_points = np.roots(self.Fu)
        # On the curve P_w = -2 b bd ad^2 y. Homogenization gives
        # P_w/eta = -2 u^(n-g-1) (b bd ad^2)(c+1/u). Cancel eta
        # symbolically, before evaluation at a Weierstrass point.
        if cov.exact is not None:
            F, a, ad, b, bd = cov.exact
            p0 = b*b*ad*ad*F - bd*bd*a*a
            p1, p2 = 2*bd*bd*ad*a, -(bd*ad)**2
            common = p0.gcd(p1).gcd(p2)
            quotient, remainder = (b*bd*ad*ad).quo_rem(common)
            if remainder:
                raise ArithmeticError("fibre cancellation does not divide P_w/eta")
            t = np.asarray([complex(c) for c in quotient.list()[::-1]])
        else:
            t = np.polymul(cov.b, np.polymul(cov.bd, np.polymul(cov.ad, cov.ad)))
        t = np.asarray(_shift_poly(t, self.c), complex)
        self.sheet_t = t[::-1]
        self.sheet_t_power = self.n - self.g - len(t)

    # ------------------------------------------------------------ geometry

    def _wpow(self, w, d):
        return complex(w) ** np.arange(d)

    def poly(self, w):
        r"""The polynomial in `u` whose roots are the fibre over ``w``."""
        return self.fibre_data(w)[0]

    def fibre_data(self, w):
        """Polynomial and derivative shared by sheets and nested rules.

        Roots are filled lazily by the walker.  The cache belongs to this
        path's chart and is discarded with it, not retained across curves.
        """
        key = complex(w)
        if key not in self._fibres:
            if len(self._fibres) >= 16384:
                self._fibres.clear()
            c = (self.cw @ self._wpow(key, self.cw.shape[1]))[::-1]
            self._fibres[key] = [c, np.polyder(c), None]
        return self._fibres[key]

    def to_u(self, x):
        r"""`u` at the point with `x`-coordinate ``x``; `0` at infinity."""
        with np.errstate(divide="ignore", invalid="ignore"):
            if self.kind == "x":
                return np.asarray(x, complex) - self.c
            u = 1.0 / (np.asarray(x, complex) - self.c)
        return np.where(np.isfinite(u), u, 0.0)

    def to_x(self, u):
        r"""The `x`-coordinate of the point at ``u``; `\infty` at `u = 0`."""
        u = np.asarray(u, complex)
        if self.kind == "x":
            return u + self.c
        with np.errstate(divide="ignore", invalid="ignore"):
            return self.c + 1.0 / u

    def eta_of(self, u, y):
        r"""`\eta` from a `y`: `y` itself, or `u^{g+1}y` at infinity."""
        if self.kind == "x":
            return np.asarray(y, complex)
        return np.asarray(u, complex) ** (self.g + 1) * np.asarray(y, complex)

    def y_of(self, u, eta):
        r"""The inverse of :meth:`eta_of`."""
        if self.kind == "x":
            return np.asarray(eta, complex)
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.asarray(eta, complex) / np.asarray(u, complex) ** (
                self.g + 1)

    def eta(self, u, near, zeta=None):
        r"""`\pm\sqrt{F_u(u)}`, identified by the cover equation.

        Continuity breaks ties at zeros of B, and selects the sheet when
        phi factors through x and both square roots have the same image.
        """
        r = np.sqrt(_eval(self.Fu, np.asarray(u, complex)))
        chosen = np.where(np.abs(r - near) <= np.abs(r + near), r, -r)
        if zeta is None or self.cov.pure_x:
            return chosen
        a, b, d = (_eval(q, u) for q in
                   (self.sheet_a, self.sheet_b, self.sheet_d))
        if self.sigma:
            lhs, rhs = d * (1 + self.zr * zeta) - a * zeta, b * r * zeta
        else:
            lhs, rhs = d * zeta - a, b * r
        plus, minus = np.abs(lhs - rhs), np.abs(lhs + rhs)
        return np.where(plus < minus, r,
                        np.where(minus < plus, -r, chosen))

    # --------------------------------------------------------- differentials

    def refine(self, u, eta, zeta, rounds=4):
        r"""
        Polish the point by Newton on `\varphi = \zeta` *along the curve*,
        which the cover polynomial cannot do near the `N`-fold point.

        `\mathcal{P}(\cdot,w)` has an `N`-fold root at `x_P`, and evaluating
        it there is the subtraction of two numbers of size `1` to get one of
        size `|w|`: on the degree-7 cover of LMFDB 249.a.249.1 the root can
        then be located to no better than `10^{-4}` of the distance to its
        neighbours, and a leg tracked that way closes to `6\cdot10^{-6}`
        where the same leg marched as an ODE closes to `10^{-12}`.

        `\varphi = \operatorname{N}(\varphi)/(A - By)` has that zero
        factored out -- :meth:`~.cover.Cover.phi_stable` -- so the residual
        is computed to full relative accuracy however small it is.  The
        step is taken as `y(\varphi - \zeta)/(y\,d\varphi/dx)` rather than
        `(\varphi - \zeta)/(d\varphi/dx)`, the denominator being the form
        that stays finite at a Weierstrass point.

        In a `\sigma` chart the equation solved is `1/(\varphi - z_r) =
        \sigma` instead, which is the same relative accuracy read at the
        other end: `\varphi` runs to infinity at `Q` and its reciprocal
        runs to zero, and it is the reciprocal that the path measures.

        Returns ``(u, eta, residual)``; the residual is the last correction,
        and is ``inf`` when this chart cannot do it, which sends the caller
        back to the polynomial's own answer.
        """
        if self.kind == "inv" and self.sigma:
            # Keep the homogeneous Newton equation on the infinity ray;
            # reconstructing phi here can overflow even when sigma is small.
            return u, eta, np.inf
        cov, z = self.local, complex(zeta)
        x, y, prev, m = complex(u), complex(eta), np.inf, np.inf
        for _ in range(rounds):
            den = complex(cov.y_dphi_dx(x, y))
            if den == 0 or not np.isfinite(den):
                return u, eta, np.inf
            val = complex(cov.phi_stable(x, y))
            if self.sigma:
                v = val - self.zr
                if v == 0 or not np.isfinite(v):
                    return u, eta, np.inf
                # d(1/v)/dx = -v'/v^2, so the step carries a factor v^2
                step = -y * (1.0 / v - z) * v * v / den
            else:
                step = y * (val - z) / den
            if not np.isfinite(step):
                return u, eta, np.inf
            x = x - step
            r = np.sqrt(complex(_eval(cov.F, x)))
            y = r if abs(r - y) <= abs(r + y) else -r
            m = abs(step)
            if m <= 1e-15 * (1.0 + abs(x)) or m >= 0.9 * prev:
                break
            prev = m
        return x, y, m

    def integrand(self, u, eta, zeta):
        r"""
        `\omega_q/d\zeta` at the point `u` over `\zeta`, as the pair
        `(x^q)_{q<g}` and the denominator they are divided by.

        `\zeta` is `w`, or `\sigma = 1/(w - z_r)` in a chart built with
        ``sigma=True``; the two differ by `d\sigma/dw = -\sigma^2`, which is
        applied to the denominator rather than to the differential.

        The two halves come back separately because the caller has a third
        factor, `d\zeta/dt`, that cancels against this one.  Approaching `Q`
        along the ray, `\sigma` vanishes to order `e` in the uniformizer, so
        the denominator is `10^{-208}` at the last quadrature node of a
        degree-14 cover and `d\zeta/dt` is `10^{-208}` too, while their
        ratio is `O(1)`: dividing first and multiplying afterwards would
        pass through `10^{208}` for no reason and overflow outright on a
        cover of higher degree.

        Both halves are taken in this chart, and neither is ever the ratio
        of two numbers that have separately blown up:

        .. math:: \frac{\omega_q}{dw} = \frac{x^q}{\eta\,dw/du},
                  \qquad \eta\,\frac{dw}{du} = -\eta\,
                  \frac{\mathcal{P}_u}{\mathcal{P}_w}.

        In the chart at a finite `x` the denominator is instead read off
        :meth:`~.cover.Cover.y_dphi_dx`, which is the same product
        `y\,d\varphi/dx` written as `y(A' + B'y) + BF'/2`.  They agree, but
        the second is right at a Weierstrass point and the first is not:
        there `\eta = \sqrt{F(x)}` cancels to zero while `dw/du` stays
        finite, so their product reads as `0` where the true value is
        `BF'/2`.  On the degree-6 model of LMFDB 249.a.249.1 the sheet
        ending at `x = -1` over `w = 2` is exactly this, and it turned the
        whole leg into ``nan``.

        In reciprocal coordinates, the transformed cover supplies the same
        norm-based derivative on finite legs. On rays in sigma, use the
        homogeneous polynomial identity below: its factors stay scaled even
        when phi itself would overflow.
        """
        u = np.asarray(u, complex)
        g = self.g
        with np.errstate(divide="ignore", invalid="ignore"):
            if self.kind == "x":
                x = u + self.c
                den = self.local.y_dphi_dx(u, eta)
                if self.sigma:
                    den = -den * complex(zeta) * complex(zeta)
                return np.array([x ** q for q in range(g)]), den
            t = self.c * u + 1.0
            num = np.array([-t ** q * u ** (g - 1 - q) for q in range(g)])
            # In reciprocal coordinates the same stable norm identity is
            # available. Use it where representable; the homogeneous form
            # below remains necessary at extreme infinity (phi overflows).
            if not self.sigma:
                with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                    stable = self.local.y_dphi_dx(u, eta)
                if np.all(np.isfinite(stable)) and np.all(stable != 0):
                    return num, stable
            c, derivative, _ = self.fibre_data(zeta)
            Pu = _eval(derivative, u)
            if self.cov.pure_x:
                dc = (self.dcw @ self._wpow(zeta, self.dcw.shape[1]))[::-1]
                den = -eta * Pu / _eval(dc, u)
            else:
                # For sigma=1/(w-zr), (sigma^2 P)_sigma/eta = +2 T.
                # This stays scaled at infinity and removes the artificial
                # 0/0 at eta=0 without evaluating huge rational A and B.
                T = _eval(self.sheet_t, u) * u ** self.sheet_t_power
                den = (-1 if self.sigma else 1) * Pu / (2 * T)
            t = self.c * u + 1.0
            return (np.array([-t ** q * u ** (g - 1 - q) for q in range(g)]),
                    den)

    def dw_du(self, u, w):
        r"""`-\mathcal{P}_u/\mathcal{P}_w`, by implicit differentiation of
        the cover polynomial read in this chart."""
        u = np.asarray(u, complex)
        c = self.poly(w)
        Pu = _eval(np.polyder(c), u) if len(c) > 1 else np.zeros_like(u)
        dc = (self.dcw @ self._wpow(w, self.dcw.shape[1]))[::-1]
        with np.errstate(divide="ignore", invalid="ignore"):
            return -Pu / _eval(dc, u)

    def __repr__(self):
        center = "%g" % self.c.real if self.c.imag == 0 else str(self.c)
        return ("u = x - %s" % center if self.kind == "x"
                else "u = 1/(x - %s)" % center)


def charts_over(cov, z, sigma=False, zr=0.0):
    r"""
    The coordinates to race for a walk ending over ``z``, one per place.

    Sage's ``local_uniformizer`` at each place of `\varphi^{-1}(z)` says
    which coordinate is the right one *there*: `1/x` at a place at infinity,
    `x - x_0` at a finite one.  They are collected here and raced rather than
    chosen, because which sheet ends at which place is not known before the
    walk -- that is what the walk decides.

    Both `x` and `1/x` are offered, with `1/x` first when the target
    fibre contains a point at infinity.  Finite poles use centered
    coordinates so that their small displacements remain representable.  ``z = None`` asks for the places over `\varphi =
    \infty`, and ``sigma`` reads the cover polynomial in `1/w`, which is
    what the ray leg wants at both ends.

    EXAMPLES::

        sage: from hyperell_regulator.regulator.cover import Cover
        sage: from hyperell_regulator.regulator.polynomials import ring
        sage: from hyperell_regulator.regulator.walks import charts_over
        sage: x = ring().gen()
        sage: [repr(c) for c in charts_over(Cover(x^5 + (1+x)^2, x + 1, 1), 2)]
        ['u = x - 0', 'u = 1/(x - 0)']
    """
    if z is None:
        poles = cov.norm_factors()[1]
        if poles:
            # At a finite pole xQ, x rounds to xQ long before the ray's
            # last node.  Only the centered coordinate retains x-xQ;
            # evaluating its differential must use shifted coefficients too.
            return [Chart(cov, "x", r, sigma, zr) for r, _ in poles]
    out = [Chart(cov, "x", 0.0, sigma, zr),
           Chart(cov, "inv", 0.0, sigma, zr)]
    try:
        pls = cov.places_over(z) or []
    except Exception:                                          # noqa: BLE001
        return out
    if z is None and pls and all(q["infinite"] for q in pls):
        # All poles lie at infinity: x is outside the endpoint chart and
        # cannot win a useful trial. In even degree 1/x is a uniformizer;
        # in odd degree it has order two (handled by the ray substitution).
        return [out[1]]
    if any(q["infinite"] for q in pls):
        # The trial ends short of the target, so x can still be finite
        # while its derivative has already lost digits by cancellation.
        # Prefer the chart supplied by the known point at infinity; a
        # failed walk there can still fall back to a finite chart.
        out.reverse()
    seen = {0.0}
    if z is not None:
        # A high-order zero can dominate conditioning even when the
        # target is another value. Offer its centered coordinate on every
        # finite leg, preserving common factors of A, B and F exactly.
        for root, multiplicity in cov.norm_factors()[0]:
            if complex(root).imag == 0 and complex(root).real not in seen:
                seen.add(complex(root).real)
                out.insert(0, Chart(cov, "x", root, sigma, zr))
    for q in pls:
        if q["infinite"]:
            continue                      # already offered as 1/x
        # a uniformizer x - x0 names the place's x-coordinate; a shift there
        # puts the place at u = 0, which is where Newton is best behaved
        try:
            r = cov.exact[0].parent()(q["uniformizer"].numerator())
            for x0 in r.roots(multiplicities=False):
                key = complex(x0).real
                if key not in seen and abs(complex(x0)) < 1e6:
                    seen.add(key)
                    # Prefer the endpoint coordinate; global x can lose
                    # digits through cancellation near a high-order zero.
                    out.insert(0, Chart(cov, "x", complex(x0), sigma, zr))
        except Exception:                                      # noqa: BLE001
            pass
    return out


def _newton(c, seed, rounds=8, tol=1e-14, mile=1e-8, derivative=None):
    r"""
    Newton on the cover polynomial ``c``, already read in a chart.

    Returns the root, an iteration count for the step control to read, and
    the size of the last correction, which is what the root is located to.

    Three things here are not the textbook loop, and all are about the
    difference between *finding* a root and *polishing* it.

    Convergence is declared at ``tol`` *relative to the root* or, at any
    scale, as soon as the step stops contracting.  How well a polynomial's root can be located is set
    by the conditioning of its evaluation, and no constant knows that
    number: on the pinch model the floor is `3\cdot10^{-14}` against a
    tolerance of `2.6\cdot10^{-14}`, while on the degree-7 cover of LMFDB
    249.a.249.1 approaching its sevenfold point the same polynomial can only
    be evaluated to `5\cdot10^{-10}`, a factor of twenty thousand apart.
    Both stalled a walk that asked for more than was there.  Nothing is lost
    by not asking: the caller is handed the last step and rejects the point
    unless it is located far better than its neighbours are far away, which
    is the only sense in which the root has been identified at all.

    At a genuine double root, where stopping early would be wrong, Newton
    still halves the step every time and never looks stagnant.

    The count returned is the iteration that first got within ``mile``, not
    the one that finished.  The step control grows the step when Newton
    converges in four iterations and holds it at five to six, so counting
    the polishing makes an easy step look like a hard one: on that same
    degree-7 cover the two extra iterations it takes to reach the floor
    pinned the count at six for a whole leg, the step stayed at
    `9\cdot10^{-7}`, and the walk ran out of its budget four thousand steps
    later still three-quarters of the way along.
    """
    d = derivative if derivative is not None else np.polyder(c)
    u = complex(seed)
    prev, k, m = np.inf, 0, np.inf
    for j in range(1, rounds + 1):
        checkpoint("newton_iterations")
        den = _eval(d, u)
        if den == 0:
            return u, rounds + 1, m
        step = _eval(c, u) / den
        u -= step
        # relative to the root, not to one: approaching Q the fibre sits at
        # |u| = 10^-23, where a tolerance of 1e-14 * (1 + |u|) is satisfied
        # by the first step whatever it is, and Newton returns after one
        # iteration having placed the root to a part in twenty
        m, sc = abs(step), abs(u) or 1.0
        if k == 0 and m <= mile * sc:
            k = j
        if m <= tol * sc or m >= 0.9 * prev:
            return u, k or j, m
        prev = m
    return u, rounds + 1, m


class WalkBudgetExceeded(RuntimeError):
    """A resumable walk reached its cumulative step allowance."""


class Walker:
    r"""
    One point of the fibre, moving along the path in a fixed chart.

    The step is halved whenever Newton is slow or the correction moves
    further than ``sep`` of the distance to the nearest other root, which is
    what stops the walk changing sheets; it is raised by half again whenever
    Newton converges in four iterations or fewer.  A walk that cannot keep
    to that raises :exc:`RuntimeError`, and it is by racing those failures
    that the right chart is found.

    The step is kept across calls to :meth:`to`, so a leg that stops at a
    hundred quadrature nodes does not restart the control a hundred times.
    """

    __slots__ = ("chart", "coord", "u", "eta", "t", "h", "steps", "sep")

    def __init__(self, chart, coord, u, eta, sep=0.25):
        self.chart, self.coord, self.sep = chart, coord, sep
        self.u, self.eta = complex(u), complex(eta)
        self.t, self.h, self.steps = 0.0, 1.0 / 16, 0

    def to(self, t1, budget=4096):
        r"""
        Move to the parameter ``t1``, in as many steps as it takes.

        ``self.h`` is the control step and survives the stop.  Clipping it
        to land on a quadrature node is not information about the tracking
        and must not be mistaken for any: tanh-sinh piles its nodes up
        exponentially at the ends, and a walker that took the clipped step
        for its state spent six steps per node climbing back out of
        `10^{-10}` instead of one, or -- worse, on a rejection -- kept half
        of it and never recovered.  A *halving* is information, and that is
        what the control step follows.
        """
        t1 = float(t1)
        while self.t < t1:
            rest = t1 - self.t
            h, cut = min(self.h, rest), False
            while True:
                if self.steps >= budget:
                    if cut:
                        self.h = h
                    raise WalkBudgetExceeded("budget")
                self.steps += 1
                checkpoint("walk_attempts", chart=repr(self.chart),
                           t=self.t, h=h, walker_steps=self.steps)
                zeta = self.coord(self.t + h)
                data = self.chart.fibre_data(zeta)
                c, dc, roots = data
                un, k, res = _newton(c, self.u, derivative=dc)
                if roots is None:
                    roots = data[2] = np.roots(c)
                d = np.abs(roots - un)
                # A single distinct x-root still cannot jump through a
                # pole of the coordinate: there is no neighbouring root to
                # limit its movement, so use the coordinate's own scale.
                gap = (np.partition(d, 1)[1] if len(d) > 1
                       else max(1.0, abs(self.u)))
                # the polynomial finds which root this is; phi_stable then
                # says where it is, to a relative accuracy the polynomial
                # cannot reach near an N-fold zero
                etan = complex(self.chart.eta(un, self.eta, zeta))
                un2, etan2, res2 = self.chart.refine(un, etan, zeta)
                refined = (np.isfinite(res2) and res2 <= 0.01 * gap
                           and abs(un2 - un) <= 0.01 * gap)
                if refined:
                    # Polynomial corrections can round to zero despite a
                    # poor true phi residual. Prefer the norm-based solve,
                    # provided it stays in this root's separation ball.
                    un, etan, res = un2, etan2, res2
                    if res2 <= 1e-14 * (1.0 + abs(un)):
                        # This solve converged in at most four polishing
                        # iterations. A stalled expanded-polynomial solve
                        # must not keep an already tiny control step frozen.
                        k = min(k, 4)
                        checkpoint("stable_newton_steps")
                continuous = True
                if abs(etan-self.eta) > abs(etan+self.eta):
                    # Opposite y can project almost onto the same x at B=0.
                    # In a branch-point-free u-disk, bound the argument of
                    # F(u)/F(u_old) by the sum of the root-factor bounds.
                    # If it is < pi, its continued square root cannot move
                    # into the opposite half-plane. Reject this sheet jump.
                    distances = np.abs(self.chart.branch_points-self.u)
                    movement = abs(un-self.u)
                    if len(distances) and movement < np.min(distances):
                        angle = np.sum(np.arcsin(movement/distances))
                        continuous = angle >= np.pi
                # the point must have moved less than a quarter of the way
                # to its nearest neighbour -- that is what keeps it on its
                # own sheet -- and must be pinned a hundred times finer than
                # that, or it has not been identified at all
                if (continuous and (k <= 6 or refined) and np.isfinite(un) and np.isfinite(etan)
                        and res <= 0.01 * gap
                        and abs(un - self.u) <= self.sep * max(gap, 1e-300)):
                    checkpoint("walk_accepted")
                    break
                checkpoint("walk_rejected", newton_count=k,
                           sheet_jump=not continuous,
                           relative_correction=float(res / gap),
                           relative_movement=float(abs(un - self.u) / gap))
                h, cut = h * 0.5, True
                # Reject a step that cannot advance the parameter.  In
                # particular, never snap to an unchecked Newton root when
                # the remaining interval is only a few ulps wide.  Rays
                # use a logarithmic parameter to avoid this loss of
                # resolution near their ramified endpoint.
                if h < 1e-10 * rest or self.t + h == self.t:
                    raise RuntimeError("stalled at %r" % (zeta,))
            self.eta, self.u, self.t = etan, un, self.t + h
            if cut:
                self.h = h
            elif h >= self.h and k <= 4:
                self.h = h * 1.5
        return self.u, self.eta

    def xy(self):
        r"""The point in `(x, y)`, which is `(\infty, \infty)` at `u = 0`
        in the chart at infinity."""
        return (complex(self.chart.to_x(self.u)),
                complex(self.chart.y_of(self.u, self.eta)))


def walker(cov, coord, x0, y0, chart, sep=0.25):
    r"""A :class:`Walker` started at the point `(x_0, y_0)` of the fibre."""
    u = complex(chart.to_u(x0))
    return Walker(chart, coord, u, complex(chart.eta_of(u, y0)), sep)


def race(cov, coord, x0, y0, charts, budgets=(64, 256, 1024, 4096), big=1e8,
         end=1.0, _fallback=True):
    r"""
    Walk the point in every chart and keep whichever finishes first.

    Trials receive increasing cumulative step budgets, in a fixed chart
    order. Accepted points are retained when a budget expires; increasing
    the budget resumes the same walk. The answer does not depend on timing. Returns
    ``(x, y, steps, chart)``.

    Finishing is not quite enough.  A chart whose `u` has run to ``big`` has
    not carried the point to the endpoint, it has followed it out of its own
    domain: approaching `Q` along the ray every sheet leaves the `x`-line,
    and `x` tracks them perfectly well out to `10^{24}` -- distinct roots,
    Newton converging -- while `y\,d\varphi/dx` there is `x^{N+g}`, which on
    a degree-14 cover is `10^{384}` and simply does not exist.  Such a walk
    is kept only as a last resort, after every chart has been tried.
    """
    last, fallback = None, None
    pending = [(ch, walker(cov, coord, x0, y0, ch)) for ch in charts]
    for budget in tuple(budgets) + (None,):
        if budget is None:
            # Every chart finished badly or not at all.  Fewest steps is a
            # proxy for best conditioned and it fails here: running into Q
            # the x-chart follows the sheet out of its own domain in three
            # hundred steps while 1/x, which is right, needs five thousand.
            # So before settling for the bad answer, give the others room.
            if fallback is None:
                break
            budget = max(budgets) * 8
        remaining = []
        for ch, wk in pending:
            checkpoint("chart_trials", chart=repr(ch), trial_budget=budget)
            try:
                wk.to(end, budget)
            except WalkBudgetExceeded as exc:
                last = exc
                remaining.append((ch, wk))
                continue
            except RuntimeError as exc:
                last = exc
                continue
            x, y = wk.xy()
            if abs(wk.u) <= big:
                return x, y, wk.steps, ch
            if fallback is None:
                fallback = (x, y, wk.steps, ch)
        pending = remaining
        if not pending:
            break
    if _fallback and len(charts) == 1 and charts[0].kind == "inv":
        # The endpoint can be infinity while this particular inverse chart
        # has a pole along the path (x=c). A translated inverse is still an
        # endpoint uniformizer. Try it only after the first chart fails.
        ch = charts[0]
        for shift in (1j, -1j):
            alternative = Chart(cov, "inv", ch.c + shift, ch.sigma, ch.zr)
            try:
                return race(cov, coord, x0, y0, [alternative], budgets, big,
                            end=end, _fallback=False)
            except RuntimeError as exc:
                last = exc
    if fallback is not None:
        return fallback
    raise RuntimeError("no chart could walk the point: %s" % last)


def choose(cov, coord, x0, y0, charts, **kwds):
    r"""
    The winning chart for each sheet, by racing a trial walk of the path.

    A chart that fails anywhere on the path fails the whole walk, so the
    winner is good for the path entire and not only at its far end.
    """
    x0 = np.atleast_1d(np.asarray(x0, complex))
    y0 = np.atleast_1d(np.asarray(y0, complex))
    return [race(cov, coord, x0[i], y0[i], charts, **kwds)[3]
            for i in range(len(x0))]


def leg(cov, coord, rate, x0, y0, lam0, nodes, wts, charts=None,
        cap=20000, wof=None, end=1.0, chosen=None, samples=None,
        sheet_ids=None, **kwds):
    r"""
    `\bigl(\int\omega_q, \int\lambda\,\omega_q\bigr)` along ``coord``,
    one sheet at a time and each in its own chart.

    ``nodes, wts`` are the quadrature rule on `[0,end]`; ``coord`` and ``rate``
    parametrise the path in the charts' variable -- `w`, or `\sigma = 1/w`
    -- and its derivative.  ``wof`` gives `w` itself when that is not
    ``coord``, `\lambda = \log\varphi` being a branch of `\log w` whichever
    variable the path is drawn in.  Returns ``(x, y, lam, I, L)`` with the
    fibre at the far end.

    The walk visits the nodes in order, taking as many steps between two of
    them as the tracking needs, so the mesh that resolves the *fibre* never
    dictates the mesh that resolves the *integrand*.  ``chosen`` reuses a
    chart per sheet from a previous trial walk when refining the quadrature.
    ``samples`` reuses verified nodes between nested rules; it belongs to
    this same path, logarithm branch, and sheet ordering only.
    """
    g = cov.g
    x0 = np.atleast_1d(np.asarray(x0, complex))
    y0 = np.atleast_1d(np.asarray(y0, complex))
    n = len(x0)
    if charts is None:
        charts = charts_over(cov, coord(end))
    if wof is None:
        wof = coord
    won = (choose(cov, coord, x0, y0, charts, end=end, **kwds)
           if chosen is None else chosen)
    I = np.zeros((g, n), complex)
    L = np.zeros((g, n), complex)
    xe, ye = np.zeros(n, complex), np.zeros(n, complex)
    lam0 = np.atleast_1d(np.asarray(lam0, complex))
    lam = lam0.copy()
    try:
        wend = complex(wof(end))
        wend = wend if np.isfinite(wend) and wend != 0 else None
    except ZeroDivisionError:
        wend = None
    for i in range(n):
        sheet = i if sheet_ids is None else sheet_ids[i]
        checkpoint("quadrature_sheets", sheet=sheet)
        ch = won[i]
        wk = walker(cov, coord, x0[i], y0[i], ch)
        li, prev = complex(lam0[i]), 0.0
        cache = {} if samples is None else samples[i]
        for t, wt in zip(nodes, wts):
            if t in cache:
                checkpoint("cached_nodes")
                u, eta, li, f = cache[t]
                # This point was already tracked and checked on this same
                # path at a coarser level of the nested quadrature rule.
                wk.u, wk.eta, wk.t = u, eta, float(t)
            else:
                checkpoint("quadrature_nodes", sheet=sheet)
                u, eta = wk.to(t, cap)
                li += np.log(wof(t) / wof(prev))
                num, den = ch.integrand(u, eta, coord(t))
                with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
                    f = num.ravel() * (rate(t) / den)
                if not np.all(np.isfinite(f)):
                    raise RuntimeError("non-finite differential in %r at %r" % (ch, t))
                cache[t] = (u, eta, li, f)
            prev = t
            I[:, i] += wt * f
            L[:, i] += wt * li * f
        wk.to(end, cap)
        xe[i], ye[i] = wk.xy()
        # lambda = log phi is infinite at Q, where the ray ends; the caller
        # does not use it there, and the last node's value is what there is
        lam[i] = li + (np.log(wend / wof(prev)) if wend else 0.0)
    return xe, ye, lam, I, L
