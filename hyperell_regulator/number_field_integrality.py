r"""Sufficient multidegree integrality over a number field.

At every bad prime ideal we follow K-rational sections through the
hyperelliptic charts and successive blowups at singular reductions, with
exact arithmetic in K. A blowup of ``uv=pi^n`` at a section through its
center is ``UV=pi^(n-2)`` in the uniformizer chart (Stacks Project 0CDB).
Smooth special-fibre points then identify components of a regular model.
Sections still singular after ``max_blowups`` remain unknown.
"""

from sage.all import (
    HyperellipticCurve, Infinity, PolynomialRing, QQ, ZZ,
    gcd, lcm, matrix, prime_range,
)

from .integrality import TorsionPacketIntegrality


def _curve_polynomials(C, K):
    R = PolynomialRing(K, 'x')
    f, h = (R(C.f), R(C.h)) if hasattr(C, 'f') else map(R, C.hyperelliptic_polynomials())
    F = h*h+4*f
    if F.degree() not in (5, 6) or not F.is_squarefree():
        raise ValueError('a smooth genus 2 hyperelliptic curve is required')
    if f.degree() > 6 or h.degree() > 3:
        raise NotImplementedError('use a hyperelliptic model with deg(f)<=6 and deg(h)<=3')
    return f, h


def _binary_discriminant(f, h):
    F = h*h+4*f
    return F.discriminant()*(F[5]**2 if F.degree() == 5 else 1)


def _factors(f, h):
    """Polynomial-in-x coefficients of all irreducible monic y-factors."""
    R = f.parent()
    S = PolynomialRing(R.fraction_field(), 'Y')
    Y = S.gen()
    return [[R(a) for a in g] for g, _ in (Y*Y+S(h)*Y-S(f)).factor()]


def _value(factor, x, y):
    return sum(a(x)*y**j for j, a in enumerate(factor))


def _infinity_factor(factor):
    """Closure of an affine component in t=1/x, z=y/x^3."""
    R = factor[0].parent()
    t = R.gen()
    weight = max(i+3*j for j, a in enumerate(factor) for i, c in enumerate(a) if c)
    result = [sum(c*t**(weight-i-3*j) for i, c in enumerate(a) if c)
              for j, a in enumerate(factor)]
    common = min(a.valuation() for a in result if a)
    return [R(a//t**common) for a in result]


class _LocalSections:
    """Exact section reductions on one regular model over O_(K,q)."""

    def __init__(self, f, h, prime, *, max_blowups):
        self.f, self.h, self.prime = f, h, prime
        self.K = f.base_ring()
        self.k = prime.residue_field()
        self.reduction = self.k.reduction_map()
        self.lift = self.k.lift_map()
        self.R = PolynomialRing(self.k, 'x')
        self.max_blowups = max_blowups
        self.uniformizer = None
        self.charts = {}
        self.main = _factors(self.reduce_poly(f), self.reduce_poly(h))
        self.main_infinity = [_infinity_factor(g) for g in self.main]
        self.inf_f = f.parent()([f[6-i] for i in range(7)])
        self.inf_h = h.parent()([h[3-i] for i in range(4)])

    def valuation(self, a):
        return Infinity if a == 0 else a.valuation(self.prime)

    def reduce_poly(self, f):
        return self.R([self.reduction(a) for a in f])

    def chart(self, patch, path, f, h):
        key = patch, path
        if key not in self.charts:
            fb, hb = self.reduce_poly(f), self.reduce_poly(h)
            factors = ((self.main if patch == 'finite' else self.main_infinity)
                       if not path else _factors(fb, hb))
            self.charts[key] = fb, hb, factors
        return self.charts[key]

    def point(self, P):
        """Return a component label or unknown, plus a reproducible path."""
        if P[0] == 'infinity':
            patch, a, b = 'infinity', self.K(0), P[1]
        elif self.valuation(P[0]) < 0:
            patch, a, b = 'infinity', 1/P[0], P[1]/P[0]**3
        else:
            patch, (a, b) = 'finite', P
        f, h = (self.f, self.h) if patch == 'finite' else (self.inf_f, self.inf_h)
        path, steps = (), []
        for depth in range(self.max_blowups+1):
            if self.valuation(a) < 0 or self.valuation(b) < 0:
                return None, 'the section is not integral in the selected chart', steps
            if b*b+h(a)*b != f(a):
                raise ArithmeticError('the transformed section left the exact curve equation')
            fb, hb, factors = self.chart(patch, path, f, h)
            ab, bb = self.reduction(a), self.reduction(b)
            smooth = (2*bb+hb(ab) != 0 or hb.derivative()(ab)*bb-fb.derivative()(ab) != 0)
            steps.append(dict(patch=patch, residue_x=str(ab), residue_y=str(bb),
                              reduced_f=str(fb), reduced_h=str(hb), smooth=bool(smooth)))
            if smooth:
                hits = [j for j, g in enumerate(factors) if _value(g, ab, bb) == 0]
                if len(hits) != 1:
                    raise ArithmeticError('a smooth section must lie on a unique component')
                label = ('main:%s' % hits[0] if not path else
                         '%s:%s:component:%s' % (patch, '/'.join(path), hits[0]))
                return label, None, steps
            if depth == self.max_blowups:
                return None, 'the blowup budget was exhausted', steps
            if self.uniformizer is None:
                self.uniformizer = self.K.uniformizer(self.prime)
                if self.valuation(self.uniformizer) != 1:
                    raise ArithmeticError('the local uniformizer does not have valuation one')
            pi = self.uniformizer
            A, B = self.K(self.lift(ab)), self.K(self.lift(bb))
            x = f.parent().gen()
            shifted_f, shifted_h = f(A+pi*x), h(A+pi*x)
            # The section through the centre forces the constant term
            # to be divisible by pi^2. The vanishing Jacobian at the
            # centre forces both linear terms to be divisible by pi^2
            # as well. The Y^2 coefficient has valuation exactly two,
            # so this is the saturated equation of the pi-chart.
            f = (shifted_f-B*shifted_h-B*B)/pi**2
            h = (shifted_h+2*B)/pi
            if any(self.valuation(c) < 0 for c in list(f)+list(h)):
                raise ArithmeticError('the uniformizer blowup chart is not integral')
            tangent_discriminant = self.reduction(h[1]**2+4*f[2])
            centre_name = 'node' if tangent_discriminant else 'centre'
            path += ('%s(%s,%s)' % (centre_name, ab, bb),)
            a, b = (a-A)/pi, (b-B)/pi
        raise AssertionError('unreachable blowup state')


class NumberFieldTorsionPacketIntegrality(TorsionPacketIntegrality):
    r"""Multidegree criterion for K-rational sections of a genus 2 curve.

    ``field`` is an absolute number field containing the curve coefficients
    and all point coordinates. Points use the original ``y^2+h*y=f`` model,
    with ``('infinity', eta)`` for ``eta=lim(y/x^3)``.

    ``torsion_multiples`` (one integer, or one per point) annihilate each
    ``P_i-P_0``; we verify them and reduce to exact orders. If omitted, we
    bound torsion over K using good unramified reductions.

    We check every prime above a bad rational prime, or, without reduction
    data, every prime dividing the discriminant of an integral model.
    ``check``, ``combine`` and ``kernel`` behave as in
    :class:`TorsionPacketIntegrality`.
    """

    def __init__(self, C, points, *, field, torsion_multiples=None,
                 torsion_primes=30, max_blowups=32):
        if field is QQ or not field.is_absolute():
            raise ValueError('field must be an absolute number field')
        if max_blowups < 0:
            raise ValueError('max_blowups must be nonnegative')
        self.curve, self.field = C, field
        self.f, self.h = _curve_polynomials(C, field)
        self.points = tuple(self._normalize_point(P) for P in points)
        if not self.points or len(set(self.points)) != len(self.points):
            raise ValueError('the torsion packet must be nonempty with distinct points')
        n = len(self.points)
        self.torsion_evidence = []
        if torsion_multiples is None:
            common = self._torsion_multiple(torsion_primes)
            if common is None:
                raise NotImplementedError('no suitable good unramified reductions; supply torsion_multiples')
            multiples = [common]*n
        elif torsion_multiples in ZZ:
            multiples = [ZZ(torsion_multiples)]*n
        else:
            multiples = list(torsion_multiples)
        if len(multiples) != n or any(N not in ZZ or N <= 0 for N in multiples):
            raise ValueError('positive integer torsion multiples are required for every point')
        self.orders = self._verify_orders(multiples)

        # Scaling y alone gives a globally integral model without changing
        # its x-map; keep the original coordinates in the public interface.
        coefficients = list(self.f)+list(self.h)
        self.y_scale = (ZZ(1) if all(c.is_integral() for c in coefficients) else
                        lcm([c.denominator() for c in coefficients]))
        f, h = self.y_scale**2*self.f, self.y_scale*self.h
        local_points = [(P[0], self.y_scale*P[1]) for P in self.points]
        rational_curve = all(c in QQ for c in coefficients)
        if rational_curve and hasattr(C, 'bad_primes'):
            primes = [q for p in C.bad_primes() for q in field.primes_above(p)]
            # Scaling may introduce singularities at additional model primes.
            primes += [q for q, _ in field.ideal(self.y_scale).factor() if q not in primes]
        else:
            primes = [q for q, _ in field.ideal(2*_binary_discriminant(f, h)).factor()]
        self.bad_primes = tuple(sorted(primes, key=lambda q: (int(q.smallest_integer()), str(q))))
        self.local_data = {}
        rows = [[ZZ(1)]*n]
        self.constraint_labels = [('degree', None)]
        for q in self.bad_primes:
            labels, reasons, paths = [None]*n, {}, {}
            local = _LocalSections(f, h, q, max_blowups=max_blowups)
            for i, P in enumerate(local_points):
                label, reason, steps = local.point(P)
                labels[i], paths[i] = label, steps
                if reason is not None:
                    reasons[i] = reason
            for label in sorted(set(labels)-{None}):
                rows.append([ZZ(a == label) for a in labels])
                self.constraint_labels.append((q, label))
            for i, label in enumerate(labels):
                if label is None:
                    rows.append([ZZ(i == j) for j in range(n)])
                    self.constraint_labels.append((q, 'unresolved point %s' % i))
            self.local_data[q] = dict(
                components=labels, unresolved=reasons, paths=paths,
                rational_prime=int(q.smallest_integer()),
                ramification_index=int(q.ramification_index()),
                residue_degree=int(q.residue_class_degree()),
                uniformizer=local.uniformizer,
            )
        self.constraints = matrix(ZZ, rows)
        self.complete = all(not d['unresolved'] for d in self.local_data.values())

    def _normalize_point(self, P):
        K = self.field
        if len(P) != 2:
            raise ValueError('use (x,y) or ("infinity", eta) coordinates')
        if P[0] == 'infinity':
            eta = K(P[1])
            if eta*eta+self.h[3]*eta != self.f[6]:
                raise ValueError('invalid point at infinity')
            return 'infinity', eta
        a, b = K(P[0]), K(P[1])
        if b*b+self.h(a)*b != self.f(a):
            raise ValueError('point is not on the original curve over the chosen field')
        return a, b

    def _torsion_multiple(self, bound):
        """Good unramified odd reductions give injective torsion reduction."""
        K = self.field
        discriminant = _binary_discriminant(self.f, self.h)
        multiple = None
        for p in prime_range(3, bound+1):
            for q in K.primes_above(p):
                if q.ramification_index() != 1 or q.norm() > 200:
                    continue
                if discriminant.valuation(q) != 0 or any(c.valuation(q) < 0 for c in list(self.f)+list(self.h) if c):
                    continue
                k = q.residue_field()
                R = PolynomialRing(k, 'x')
                red = k.reduction_map()
                D = HyperellipticCurve(R([red(c) for c in self.f]), R([red(c) for c in self.h]))
                size = ZZ(D.frobenius_polynomial()(1))
                multiple = size if multiple is None else gcd(multiple, size)
                self.torsion_evidence.append(dict(prime=q, cardinality=size))
                if len(self.torsion_evidence) >= 3 or multiple == 1:
                    return multiple
        return multiple

    def _verify_orders(self, multiples):
        F = self.h*self.h+4*self.f

        def weierstrass(P):
            return (F.degree() == 5 if P[0] == 'infinity' else
                    2*P[1]+self.h(P[0]) == 0)

        O = self.points[0]
        rational_curve = (hasattr(self.curve, 'bad_primes') and
                          all(c in QQ for c in list(self.f)+list(self.h)))
        rational_base = all(a in QQ for a in O if a != 'infinity')
        base = None
        orders = [ZZ(1)]
        for P, N in zip(self.points[1:], multiples[1:]):
            if weierstrass(P) and weierstrass(O):
                # The x-map gives 2(P-O) principal. A degree-one function
                # cannot exist on a genus-two curve, so the order is 2.
                if N % 2:
                    raise ValueError('the supplied torsion multiple does not annihilate a Weierstrass difference')
                orders.append(ZZ(2))
                continue
            if rational_curve and rational_base and all(a in QQ for a in P if a != 'infinity'):
                from .components import divisor_class_order
                def rational_point(A):
                    return ('infinity', QQ(A[1])) if A[0] == 'infinity' else tuple(map(QQ, A))
                order = divisor_class_order(self.curve, rational_point(P), rational_point(O))
                if order is not None:
                    if order is Infinity or N % order:
                        raise ValueError('the supplied torsion multiple does not annihilate a rational difference')
                    orders.append(ZZ(order))
                    continue
            if base is None:
                base = self._odd_model_jacobian()
            J, point_class = base
            difference = point_class(P)-point_class(O)
            order = ZZ(N)
            if order*difference != J(0):
                raise ValueError('the supplied torsion multiple does not annihilate %s - basepoint' % (P,))
            for p in order.prime_divisors():
                while order % p == 0 and (order//p)*difference == J(0):
                    order //= p
            orders.append(order)
        return orders

    def _odd_model_jacobian(self):
        """Use Mumford arithmetic after adjoining one branch point if needed."""
        K = self.field
        F = self.h*self.h+4*self.f
        root = None
        if F.degree() == 6:
            factor = min((g for g, _ in F.factor()), key=lambda g: g.degree())
            if factor.degree() == 1:
                root = -factor[0]/factor[1]
            else:
                K = K.extension(factor.monic(), 'w')
                root = K.gen()
        R = PolynomialRing(K, 'u')
        u = R.gen()
        F = R(F.list())
        if root is not None:
            shifted = F(u+root)
            F = R([shifted[6-i] for i in range(7)])
        C = HyperellipticCurve(F)
        J = C.jacobian()(K)
        h = R(self.h.list())

        def point_class(P):
            x, y = P
            if x == 'infinity':
                if root is None:
                    return J(0)
                return J(C.point([K(0), 2*K(y)+h[3], K(1)]))
            x, y = K(x), K(y)
            if root is not None and x == root:
                return J(0)
            Y = 2*y+h(x)
            if root is not None:
                x, Y = 1/(x-root), Y/(x-root)**3
            return J(C.point([x, Y, K(1)]))

        return J, point_class
