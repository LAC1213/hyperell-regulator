r"""Weierstrass torsion-packet divisors, regulators and finite-place boundaries."""

from dataclasses import dataclass
from functools import cached_property
from collections.abc import Mapping

from sage.all import (
    Matrix, QQ, QQbar, RealField, lcm, number_field_elements_from_algebraics,
)

from .number_field_integrality import _curve_polynomials
from ._galois import galois_trace_basis
from .projected_integrality import ProjectedTorsionPacketIntegrality
from .regulator import regulators_stokes


DETERMINANT_WEIGHTS = (1, 2, 1)


def integration_chart(F, abscissae, shift=None):
    """Choose a common real chart and retain the centres at infinity."""
    xs = tuple(abscissae)
    if shift is None and 'infinity' not in xs:
        return F, xs, None
    if shift is None:
        avoid = [complex(x) for x in xs if x != 'infinity']
        avoid += [complex(x) for x in F.roots(QQbar, multiplicities=False)]
        choices = [QQ(a) for a in range(-4, 5) if F(a) != 0 and a not in xs]
        if not choices:
            shift = next(QQ(j) for j in range(5, 6+F.degree()+len(xs))
                         if F(j) != 0 and QQ(j) not in xs)
        else:
            shift = max(choices, key=lambda a: min(abs(complex(a)-z) for z in avoid)/(1+abs(a)))
    shift = QQ(shift)
    x = F.parent().gen()
    f = F(x+shift)
    G = F.parent()([f[6-i] for i in range(7)])
    return G, tuple(0 if a == 'infinity' else
                    'infinity' if a == shift else 1/(a-shift) for a in xs), shift


class DivisorRegulators:
    """Evaluate linear combinations of logarithmic moments on one common cut.

    A point at infinity contributes zero in this chart. Each row represents
    a degree-zero divisor, so the otherwise arbitrary logarithmic
    constant has no effect on the determinant with the period column.
    """

    def __init__(self, F, abscissae, *, digits=None, certified=False):
        self.abscissae = tuple(abscissae)
        self.certified = certified
        self.digits = (15 if digits is None and certified else digits)
        self._column_cache = {}
        self.polynomial = F
        self._period = None

    def _integrate(self, support):
        xs = list(dict.fromkeys(self.abscissae[i] for i in support
                                if self.abscissae[i] != 'infinity'))
        if self._period is not None and all(x in self._column_cache for x in xs):
            return
        digits = self.digits
        if self.certified:
            from .regulator.balls import digits_value
            digits = digits_value(digits)+12
        options = {'certified': True} if self.certified else {}
        xs = list(dict.fromkeys([*self._column_cache, *xs]))
        per, raw = regulators_stokes(self.polynomial, xs, prec=digits, **options)
        # Whenever the cut changes, recompute all columns used together.
        ring = (per[0].parent() if self.certified else
                RealField(int(3.33*self.digits)+24) if self.digits is not None else float)
        self._period = tuple(ring(w*p) for w, p in zip(DETERMINANT_WEIGHTS, per))
        self._column_cache = {x: tuple(ring(w*p) for w, p in zip(DETERMINANT_WEIGHTS, a))
                              for x, a in zip(xs, raw)}
        self._ring = ring

    def columns(self, divisors):
        """Return the period column and regulator columns for exact coefficient rows."""
        rows = [tuple(map(QQ, d)) for d in divisors]
        if any(len(d) != len(self.abscissae) or sum(d) for d in rows):
            raise ValueError('one coefficient per abscissa and total degree zero are required')
        self._integrate(i for i in range(len(self.abscissae)) if any(d[i] for d in rows))
        cols = []
        for d in rows:
            terms = [(self._ring(a), self._column_cache[x]) for a, x in zip(d, self.abscissae)
                     if a and x != 'infinity']
            cols.append(tuple(sum((a*c[j] for a, c in terms), self._ring(0)) for j in range(3)))
        return self._period, cols

    def determinant(self, left, right):
        """Return the determinant and matrix in the Beilinson normalization."""
        p, (a, b) = self.columns((left, right))
        return self.determinant_of_columns(a, b, p)

    def determinant_of_columns(self, left, right, period):
        """Combine previously evaluated columns without another integration."""
        rows = [list(r) for r in zip(left, right, period)]
        det = Matrix(period[0].parent(), rows).determinant() if self.certified else Matrix(rows).determinant()
        if self.certified:
            from .regulator.balls import narrow, CertificationError
            if not narrow(det, self.digits):
                raise CertificationError('determinant did not achieve the requested absolute radius')
        return det, rows


@dataclass(frozen=True)
class MotivicDivisor:
    """A degree-zero formal divisor representing the standard projection of B."""

    packet: object
    coefficients: tuple

    def regulator(self, **options):
        """Return the real regulator column in the chosen chart."""
        return self.packet.regulators([self], **options)[1][0]

    def determinant(self, other, **options):
        """Return the determinant with another class and the period column."""
        return self.packet.determinant(self, other, **options)

    def boundary(self, prime):
        """Return the local cohomological boundary in graph coordinates."""
        return self.packet.boundary(self, prime)


def _divisor_evaluator(F, abscissae, divisors, *, digits=None, certified=False, shift=None):
    xs = tuple(dict.fromkeys(abscissae))
    groups = [tuple(i for i, a in enumerate(abscissae) if a == x) for x in xs]
    rows = [[sum(d[i] for i in g) for g in groups] for d in divisors]
    support = [i for i in range(len(xs)) if any(d[i] for d in rows)]
    F, centres, chart = integration_chart(F, [xs[i] for i in support], shift)
    evaluator = DivisorRegulators(F, centres, digits=digits, certified=certified)
    return evaluator, [[d[i] for i in support] for d in rows], chart


class MotivicPacket:
    r"""Degree-zero divisor representatives of standard motivic classes.

    The default packet consists of the six geometric Weierstrass points.
    An explicit torsion packet may include other points. Coordinates are
    exact algebraic numbers on the original rational model. Torsion of
    supplied points is verified before using their motivic classes.
    """

    def __init__(self, C, points=None, *, torsion_multiples=None, max_model_degree=None):
        self.curve = C
        self.f, self.h = _curve_polynomials(C, QQ)
        self.F = self.h*self.h+4*self.f
        self._weierstrass = points is None
        if points is None:
            points = [(x, -self.h(x)/2) for x in sorted(
                self.F.roots(QQbar, multiplicities=False), key=lambda x: (x.real(), x.imag()))]
            if self.F.degree() == 5:
                points.append(('infinity', -self.h[3]/2))
            torsion_multiples = 2
        self.points = tuple((x if x == 'infinity' else QQbar(x), QQbar(y)) for x, y in points)
        if not self.points or len(set(self.points)) != len(self.points):
            raise ValueError('the packet must contain distinct points')
        for x, y in self.points:
            f, h = (self.f[6], self.h[3]) if x == 'infinity' else (self.f(x), self.h(x))
            if y*y+h*y != f:
                raise ValueError('a point is off the curve')
        self.abscissae = tuple(P[0] for P in self.points)
        self._torsion_multiples = torsion_multiples
        self._max_model_degree = max_model_degree

    @classmethod
    def from_geometric_packet(cls, packet, *, rational_abscissae=False):
        """Use the computed Weierstrass torsion packet including nonbranch points."""
        N = 2*lcm([n for g, n in packet.factors
                   if not rational_abscissae or g is None or g.degree() == 1])
        return cls(packet.curve, packet.points(rational_abscissae=rational_abscissae),
                   torsion_multiples=N)

    @cached_property
    def integrality(self):
        """The number-field arithmetic of the packet and all its bad places."""
        flat = [v for x, y in self.points for v in ((y,) if x == 'infinity' else (x, y))]
        K, entries, _ = number_field_elements_from_algebraics(flat, embedded=True)
        entries = iter(entries)
        points = [(x if x == 'infinity' else next(entries), next(entries))
                  for x, _ in self.points]
        return ProjectedTorsionPacketIntegrality(self.curve, points, field=K,
                                                 torsion_multiples=self._torsion_multiples,
                                                 use_stable_extension=True,
                                                 max_model_degree=self._max_model_degree)

    def divisor(self, coefficients):
        """Build a degree-zero class from a point mapping or a coefficient row."""
        if isinstance(coefficients, MotivicDivisor):
            if coefficients.packet is not self:
                raise ValueError('the divisor belongs to a different packet')
            coefficients = coefficients.coefficients
        if isinstance(coefficients, Mapping):
            D = [QQ(0)]*len(self.points)
            for (x, y), a in coefficients.items():
                P = (x if x == 'infinity' else QQbar(x), QQbar(y))
                if P not in self.points:
                    raise ValueError('the support lies outside the packet')
                D[self.points.index(P)] += QQ(a)
        else:
            D = tuple(map(QQ, coefficients))
        if len(D) != len(self.points) or sum(D):
            raise ValueError('one coefficient per point and total degree zero are required')
        return MotivicDivisor(self, tuple(D))

    def _prepare_regulators(self, divisors, *, digits=None, certified=False, shift=None):
        if not self._weierstrass:
            # Construction on explicit points requires torsion verification.
            self.integrality
        divisors = [self.divisor(d).coefficients for d in divisors]
        return _divisor_evaluator(self.F, self.abscissae, divisors, digits=digits,
                                  certified=certified, shift=shift)

    def regulators(self, divisors, **options):
        """Return the period column, regulator columns and common chart shift."""
        evaluator, divisors, chart = self._prepare_regulators(divisors, **options)
        p, cols = evaluator.columns(divisors)
        return p, cols, chart

    def determinant(self, left, right, **options):
        """Compute the regulator determinant of two divisors on this packet."""
        evaluator, divisors, _ = self._prepare_regulators((left, right), **options)
        return evaluator.determinant(*divisors)[0]

    def integral_basis(self):
        """Return integral Galois-trace directions over QQ and completeness.

        Solve the local equations in the span of degree-zero Galois
        traces, modulo the symmetric ordinate identifications. Algebraic
        anti-invariant directions have zero trace and are excluded.
        Regulator determinants test independence of the resulting images
        in motivic cohomology. For the unrestricted number-field kernel,
        use ``packet.integrality.kernel()``.
        """
        traces = galois_trace_basis(self.integrality.abscissae)
        generators = [self.integrality.abscissa_divisor(v) for v in traces]
        result = self.integrality.kernel(generators)
        divisors = tuple(self.divisor(self.integrality.abscissa_divisor(v))
                         for v in result['abscissa_basis'])
        return divisors, result['complete']

    def boundary(self, divisor, prime):
        """Compute the full graph boundary functional on a supported model."""
        D = self.divisor(divisor)
        return self.integrality.graph_boundary(D.coefficients, prime)

    def stable_model(self, prime):
        """Return the local model with marks in the arithmetic abscissa order."""
        return self.integrality.stable_model(prime)
