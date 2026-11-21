r"""Sufficient integrality criterion for sums of torsion-packet classes.

For a degree-zero divisor ``D = sum(n_i P_i)`` supported on a rational
torsion packet, test ``sum_{P_i -> Z} n_i = 0`` on every component at every
bad prime. This implements the supplied cocycle/multidegree criterion for
the projection ``q_* B(D)`` to the standard motive, subject to its geometric
proof. It does not compute the exact Chow-valued boundary obstruction.

Smooth points on semistable stable components are labelled, as are smooth
integral affine points on the original model at 2. These smooth patches
persist in a regular proper model. Other points (including contracted
chains) and unavailable nonsemistable fibres remain unknown. Kernel computations force their
coefficients to zero, yielding a certified sublattice when data are partial.
"""

from collections.abc import Mapping

from sage.all import GF, Infinity, PolynomialRing, QQ, ZZ, identity_matrix, lcm, matrix, vector

from .components import affine_model, divisor_class_order


def _point(C, P):
    """Canonical rational coordinates on the original model, with validation."""
    if isinstance(P, (tuple, list)) and len(P) == 2 and P[0] == 'infinity':
        eta = QQ(P[1])
        g = C.genus()
        if eta**2 + C.h[g+1]*eta != C.f[2*g+2]:
            raise ValueError('invalid point at infinity')
        return ('infinity', eta)
    return tuple(C._coords(P))


def _integer_matrix(A):
    """Clear row denominators without changing the integral kernel."""
    rows = [list(lcm([a.denominator() for a in row]) * row) for row in A.rows()]
    return matrix(ZZ, len(rows), A.ncols(), rows)


def _smooth_affine_components(C, points, p):
    r"""Label smooth integral points on an integral affine model at ``p``.

    Factor the monic quadratic in y over F_p(x), and use the Jacobian
    criterion at each reduced point. Monicity and Gauss's lemma ensure the
    factors have polynomial coefficients, and that these are all affine
    components. Resolving a proper compactification leaves this smooth
    open unchanged, even if the whole curve has wild reduction.
    """
    p = ZZ(p)
    coefficients = list(C.f) + list(C.h)
    if any(QQ(a).denominator() % p == 0 for a in coefficients):
        raise NotImplementedError('the original affine model is not integral at %s' % p)
    R = PolynomialRing(GF(p), 'x')
    x = R.gen()
    f, h = R(C.f.list()), R(C.h.list())
    S = PolynomialRing(R.fraction_field(), 'y')
    y = S.gen()
    factors = [a for a, _ in (y*y + S(h)*y - S(f)).factor()]
    # Conversion also guards against unexpected denominators in factors.
    polynomials = [[R(a) for a in factor] for factor in factors]
    labels, reasons = [], {}
    for i, P in enumerate(points):
        reason = None
        if P[0] == 'infinity':
            reason = 'the point at infinity is outside the original affine patch'
        elif any(QQ(a).denominator() % p == 0 for a in P):
            reason = 'the point is not integral on the original affine patch'
        else:
            a, b = GF(p)(P[0]), GF(p)(P[1])
            if 2*b + h(a) == 0 and h.derivative()(a)*b - f.derivative()(a) == 0:
                reason = 'the point reduces to a singular point of the affine fibre'
            else:
                hits = [j for j, factor in enumerate(polynomials)
                        if sum(c(a)*b**k for k, c in enumerate(factor)) == 0]
                if len(hits) != 1:
                    raise ArithmeticError('a smooth point must lie on exactly one affine component')
                labels.append('affine_%s' % hits[0])
                continue
        labels.append(None)
        reasons[i] = reason
    evidence = dict(method='smooth original affine patch', prime=int(p),
                    reduced_f=str(f), reduced_h=str(h),
                    component_factors=[str(a) for a in factors])
    return labels, reasons, evidence


class TorsionPacketIntegrality:
    r"""Multidegrees and balanced combinations on one verified torsion packet.

    ``points`` are rational points on ``y^2+h*y=f``; infinity is represented
    by ``('infinity', eta)`` with ``eta=lim(y/x^(g+1))``. All differences
    from the first point are checked for torsion by exact principality or
    Jacobian arithmetic. A non-torsion difference raises ``ValueError``;
    an undetermined order raises ``NotImplementedError``.

    ``check(D)`` accepts a coefficient vector in point order, or a mapping
    from points to rational coefficients. ``kernel()`` computes a saturated
    integer lattice of degree-zero divisors. ``kernel(generators)`` instead
    finds integer combinations of the given degree-zero divisors, for
    example the columns ``N_i(P_i-Q_i)`` representing ``alpha_phi_i``.

    A failed sufficient criterion never asserts nonintegrality. A nonzero
    formal divisor need not give a nonzero or independent motivic class.
    No assertion about the additional Tate/mean-log condition for ``p'_*``
    is made here.
    """

    def __init__(self, C, points, torsion_primes=30):
        self.curve = C
        self.points = tuple(_point(C, P) for P in points)
        if not self.points:
            raise ValueError('a torsion packet must contain at least one point')
        if len(set(self.points)) != len(self.points):
            raise ValueError('torsion packet points must be distinct')
        self.orders = [ZZ(1)]
        for P in self.points[1:]:
            order = divisor_class_order(C, P, self.points[0], primes=torsion_primes)
            if order is Infinity:
                raise ValueError('points do not form a torsion packet: %s' % (P,))
            if order is None:
                raise NotImplementedError('torsion order is undetermined for %s' % (P,))
            self.orders.append(ZZ(order))

        # The point reduction API is affine. This is a Q-isomorphism and
        # bad_primes() uses the minimal discriminant, not the chart discriminant.
        self.model, images, self.shift = affine_model(C, self.points)
        self.local_data = {}
        n = len(self.points)
        rows = [[ZZ(1)] * n]       # degree zero is always required
        self.constraint_labels = [('degree', None)]
        for p in self.model.bad_primes():
            labels, reasons = [None] * n, {}
            evidence = None
            try:
                if p == 2:
                    labels, reasons, evidence = _smooth_affine_components(C, self.points, p)
                elif self.model.is_semistable_at(p) is not True:
                    raise NotImplementedError('regular components of nonsemistable fibres are unavailable')
                else:
                    for i, P in enumerate(images):
                        try:
                            red = self.model.reduce_point(p, P)
                            if red.is_smooth:
                                labels[i] = red.component.label
                            else:
                                reasons[i] = red.reason
                        except (NotImplementedError, ArithmeticError) as exc:
                            reasons[i] = str(exc)
            except (NotImplementedError, ArithmeticError) as exc:
                reasons = {i: str(exc) for i in range(n)}
            for label in sorted(set(labels) - {None}):
                rows.append([ZZ(label == lab) for lab in labels])
                self.constraint_labels.append((int(p), label))
            for i, label in enumerate(labels):
                if label is None:
                    rows.append([ZZ(i == j) for j in range(n)])
                    self.constraint_labels.append((int(p), 'unresolved point %s' % i))
            self.local_data[int(p)] = dict(components=labels, unresolved=reasons)
            if evidence is not None:
                self.local_data[int(p)]['evidence'] = evidence
        self.constraints = matrix(ZZ, rows)
        self.complete = all(not data['unresolved'] for data in self.local_data.values())

    def _normalize_point(self, point):
        return _point(self.curve, point)

    def divisor(self, coefficients):
        """Validate a rational degree-zero divisor and return its vector."""
        if isinstance(coefficients, Mapping):
            values = [QQ(0)] * len(self.points)
            for P, a in coefficients.items():
                point = self._normalize_point(P)
                if point not in self.points:
                    raise ValueError('divisor contains a point outside the torsion packet')
                values[self.points.index(point)] += QQ(a)
        else:
            values = list(coefficients)
        if len(values) != len(self.points):
            raise ValueError('one coefficient per packet point is required')
        D = vector(QQ, values)
        if sum(D) != 0:
            raise ValueError('the divisor must have degree zero')
        return D

    def combine(self, divisors, coefficients):
        r"""The divisor of ``sum_j coefficients[j] B(divisors[j])``."""
        divisors, coefficients = list(divisors), list(coefficients)
        if len(divisors) != len(coefficients):
            raise ValueError('one coefficient per divisor is required')
        D = vector(QQ, len(self.points))
        for E, a in zip(divisors, coefficients):
            D += QQ(a) * self.divisor(E)
        return D

    def check(self, coefficients):
        r"""Test multidegree zero; only ``q_integral=True`` is a certificate.

        ``criterion_satisfied`` is true, false, or unknown. Local
        multidegrees omit unresolved points, whose nonzero coefficients
        make that local answer unknown. Known nonzero multidegrees in a
        fully resolved fibre give false, meaning only that this sufficient
        criterion fails. Missing components have multidegree zero.
        """
        D = self.divisor(coefficients)
        local = {}
        for p, data in self.local_data.items():
            degrees, unresolved = {}, {}
            for i, (a, label) in enumerate(zip(D, data['components'])):
                if label is None:
                    if a:
                        unresolved[i] = data['unresolved'][i]
                else:
                    degrees[label] = degrees.get(label, QQ(0)) + a
            answer = None if unresolved else all(a == 0 for a in degrees.values())
            local[p] = dict(multidegrees=degrees, unresolved=unresolved,
                            criterion_satisfied=answer)
        answers = [d['criterion_satisfied'] for d in local.values()]
        answer = False if False in answers else (None if None in answers else True)
        return dict(divisor=D, criterion_satisfied=answer,
                    q_integral=True if answer is True else None, local=local)

    def kernel(self, generators=None):
        r"""A saturated integer lattice of combinations certified by the criterion.

        Basis vectors give coefficients of ``generators`` (point divisors
        when omitted); ``divisors`` gives their actual point coefficients.
        ``complete=False`` means unresolved support was forced to zero.
        Formal relations giving the zero divisor are retained and identified
        by their zero rows in ``divisors``.
        """
        n = len(self.points)
        G = (identity_matrix(QQ, n) if generators is None else
             matrix(QQ, [self.divisor(D) for D in generators]).transpose())
        if G.nrows() != n:
            # The empty generator list still has n rows.
            if G.ncols() == 0:
                G = matrix(QQ, n, 0)
            else:
                raise ValueError('invalid generators')
        A = _integer_matrix(self.constraints.change_ring(QQ) * G)
        basis = A.right_kernel().basis_matrix()
        return dict(basis=basis, divisors=basis * G.transpose(),
                    constraints=A, complete=self.complete)


def torsion_packets(C, points, torsion_primes=30):
    """Partition points by verified torsion differences, retaining unknown pairs.

    Return ``(packets, unresolved)``; singleton packets are included. The
    latter lists comparisons with undetermined order, so a partial torsion
    computation cannot silently be presented as a complete partition.
    """
    packets, unresolved = [], []
    for P in dict.fromkeys(_point(C, P) for P in points):
        for packet in packets:
            order = divisor_class_order(C, P, packet[0], primes=torsion_primes)
            if order is None:
                unresolved.append((P, packet[0]))
            elif order is not Infinity:
                packet.append(P)
                break
        else:
            packets.append([P])
    return packets, unresolved
