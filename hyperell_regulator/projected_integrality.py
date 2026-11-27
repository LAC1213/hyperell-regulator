r"""Exact local constraints for the standard projection of torsion-packet classes.

At semistable places we solve the Laplace equation on the graph of the
special fibre. At 2 we use verified auxiliary models for the wild equations
in the small-conductor survey. Elsewhere the multidegree criterion applies.
"""

from sage.all import (
    NumberField, PolynomialRing, QQ, QQbar, identity_matrix, matrix,
    number_field_elements_from_algebraics,
)

from ._dyadic_projected import dyadic_constraints
from ._projected_reduction import LocalConstraints, component_constraints, semistable_constraints
from .integrality import _integer_matrix
from .number_field_integrality import (
    NumberFieldTorsionPacketIntegrality, _LocalSections, _curve_polynomials,
)


class ProjectedTorsionPacketIntegrality(NumberFieldTorsionPacketIntegrality):
    r"""Decide integrality of :math:`q_*B(D)` and solve for integral combinations.

    Points lie in one torsion packet on ``y^2+h*y=f``, with
    ``('infinity', eta)`` for ``eta=lim(y/x^3)``. Pass ``field`` for
    nonrational coordinates, or use :meth:`from_abscissae`.

    ``check(D)['q_integral']`` is true, false, or ``None``; false means a
    nonzero projected boundary at some place. ``kernel(generators)`` returns
    a saturated lattice of integral combinations of the generators (points
    by default). Its dimension need not equal the motivic rank.
    """

    def __init__(self, C, points, *, field=QQ, torsion_multiples=None,
                 torsion_primes=30, max_blowups=32, use_stable_extension=False,
                 max_model_degree=None):
        if field is QQ:
            R = PolynomialRing(QQ, 'a')
            field = NumberField(R.gen(), 'a')
        super().__init__(C, points, field=field, torsion_multiples=torsion_multiples,
                         torsion_primes=torsion_primes, max_blowups=max_blowups)
        # On the standard motive the simultaneous hyperelliptic involution
        # acts trivially. Pair the two representatives only if both were
        # verified to lie in this torsion packet.
        groups = {}
        for i, (x, _) in enumerate(self.points):
            groups.setdefault(x, []).append(i)
        self.abscissae = tuple(groups)
        self._groups = tuple(tuple(g) for g in groups.values())
        self._max_model_degree = max_model_degree
        self._use_stable_extension = use_stable_extension
        self._max_blowups = max_blowups
        n = len(self.points)
        self._orbit_sum = matrix(QQ, [[int(i in g) for i in range(n)] for g in self._groups])
        self.projected_local_data = {}
        f, h = self.y_scale**2*self.f, self.y_scale*self.h
        points = [(x, self.y_scale*y) for x, y in self.points]
        rational_model = all(a in QQ for a in list(f)+list(h))
        for prime in self.bad_primes:
            p = int(prime.smallest_integer())
            stable = C.local_data(p)['stable_type'] if rational_model and hasattr(C, 'local_data') else None
            stable = stable.split(',')[0] if stable else None
            if stable == 'I':
                data = LocalConstraints('potential_good_reduction',
                                        matrix(QQ, 0, len(self._groups)), True)
            else:
                local = _LocalSections(f, h, prime, max_blowups=max_blowups)
                data = self._projected_constraints(local, points, stable)
            self.projected_local_data[prime] = data
        solution = self.kernel()
        self.constraints = solution['constraints']
        self.necessary_constraints = solution['necessary_constraints']
        self.complete = solution['complete']
        self.constraint_labels = [('degree', None)]
        for prime, data in self.projected_local_data.items():
            self.constraint_labels.extend((prime, (data.method, i))
                                          for i in range(data.matrix.nrows()))
            self.constraint_labels.extend((prime, ('unresolved', i)) for i in data.unresolved)

    def _projected_constraints(self, local, points, stable):
        try:
            if local.k.characteristic() == 2:
                try:
                    return dyadic_constraints(local, points, self._groups, stable)
                except NotImplementedError:
                    pass
            return semistable_constraints(local, points, self._groups)
        except (NotImplementedError, ArithmeticError) as error:
            if self._use_stable_extension and local.k.characteristic() != 2:
                from ._projected_reduction import graph_constraints
                from .stable_model import _tame_model

                groups = tuple(tuple(points[i] for i in g) for g in self._groups)
                try:
                    model = _tame_model(local.f, local.h, local.prime, groups,
                                        self._max_model_degree)
                    return graph_constraints(model, 'stable_extension')
                except (NotImplementedError, ArithmeticError) as stable_error:
                    error = stable_error
            labels = self.local_data[local.prime]['components']
            fallback = component_constraints(self._groups, labels)
            fallback.evidence['projected_model_unavailable'] = str(error)
            return fallback

    def stable_model(self, prime):
        """Return the marked local model in abscissa order."""
        if prime in QQ:
            primes = self.field.primes_above(prime)
            if len(primes) != 1:
                raise ValueError('specify a prime ideal of the packet field')
            prime = primes[0]
        data = self.projected_local_data.get(prime)
        if data is not None and data.evidence and 'model' in data.evidence:
            return data.evidence['model']
        from ._projected_reduction import nodal_model
        from .stable_model import _dyadic_model, _tame_model

        f, h = self.y_scale**2*self.f, self.y_scale*self.h
        points = [(x, self.y_scale*y) for x, y in self.points]
        local = _LocalSections(f, h, prime, max_blowups=self._max_blowups)
        if local.k.characteristic() == 2:
            return _dyadic_model(self.curve, local, points, self._groups, self._max_model_degree)
        try:
            return nodal_model(local, points, self._groups)
        except NotImplementedError:
            pass
        groups = tuple(tuple(points[i] for i in g) for g in self._groups)
        return _tame_model(f, h, prime, groups, self._max_model_degree)

    def graph_boundary(self, coefficients, prime):
        """Evaluate the cohomological boundary of a point divisor at a place."""
        return self.stable_model(prime).boundary(self._orbit_sum*self.divisor(coefficients))

    @classmethod
    def from_abscissae(cls, C, abscissae, **kwargs):
        """Use both ordinates above each exact abscissa of a rational model.

        Abscissae may be rational or ``QQbar``; use ``None`` or
        ``'infinity'`` for infinity. First-occurrence order is kept.
        """
        f, h = _curve_polynomials(C, QQ)
        F = h*h+4*f
        xs = tuple(dict.fromkeys(None if x is None or x == 'infinity' else QQbar(x)
                                 for x in abscissae))
        if not xs:
            raise ValueError('at least one abscissa is required')
        sizes, coords = [], []
        for x in xs:
            r, delta = (h[3], F[6]) if x is None else (h(x), F(x))
            s = QQbar(delta).sqrt()
            if x is not None:
                coords.append(x)
            coords.append((-r+s)/2)
            if s:
                coords.append((-r-s)/2)
            sizes.append(2 if s else 1)
        field, elements, _ = number_field_elements_from_algebraics(coords, embedded=True)
        elements = iter(elements)
        points = []
        for x, size in zip(xs, sizes):
            a = 'infinity' if x is None else next(elements)
            points.extend((a, next(elements)) for _ in range(size))
        return cls(C, points, field=field, **kwargs)

    def abscissa_divisor(self, coefficients):
        r"""Lift coefficients in ``abscissae`` order to a degree-zero divisor."""
        coefficients = list(coefficients)
        if len(coefficients) != len(self._groups):
            raise ValueError('one coefficient per abscissa is required')
        values = [QQ(0)]*len(self.points)
        for a, group in zip(coefficients, self._groups):
            for i in group:
                values[i] = QQ(a)/len(group)
        return self.divisor(values)

    @staticmethod
    def _applicable_constraints(data, coefficients):
        evidence = data.evidence or {}
        if 'galois_trace_matrix' in evidence and all(
                all(coefficients.row(i) == coefficients.row(group[0]) for i in group)
                for group in evidence['galois_trace_orbits']):
            return evidence['galois_trace_matrix'], ()
        return data.matrix, data.unresolved

    def check(self, coefficients):
        """Evaluate the projected local conditions on a degree-zero divisor."""
        D = self.divisor(coefficients)
        a = self._orbit_sum*D
        local = {}
        for prime, data in self.projected_local_data.items():
            constraints, unresolved = self._applicable_constraints(data, matrix(QQ, [a]).transpose())
            values = constraints*a
            unsupported = [i for i in unresolved if a[i]]
            verdict = None if unsupported else (True if not values else
                                                False if data.exact else None)
            local[prime] = dict(q_integral=verdict, values=values, method=data.method,
                                exact=data.exact, unresolved=unsupported, evidence=data.evidence)
        answers = [d['q_integral'] for d in local.values()]
        verdict = False if False in answers else None if None in answers else True
        return dict(divisor=D, q_integral=verdict, local=local)

    def kernel(self, generators=None):
        r"""Solve the local conditions and return integral divisor combinations.

        Rows of ``basis`` are coefficients of the generators and rows of
        ``divisors`` the resulting point divisors. ``abscissa_basis`` gives
        independent directions after identifying the two ordinates over each
        abscissa.
        """
        n = len(self.points)
        G = (identity_matrix(QQ, n) if generators is None else
             matrix(QQ, [self.divisor(D) for D in generators]).transpose())
        if G.ncols() == 0:
            G = matrix(QQ, n, 0)
        m = G.ncols()
        degree = [sum(g) for g in G.columns()]
        sufficient, necessary = [degree], [degree]
        T = self._orbit_sum*G
        for data in self.projected_local_data.values():
            constraints, unresolved = self._applicable_constraints(data, T)
            rows = (constraints*T).rows()
            unknown = [T.row(i) for i in unresolved if T.row(i)]
            sufficient.extend(rows+unknown)
            if data.exact and not unknown:
                necessary.extend(rows)
        A = _integer_matrix(matrix(QQ, len(sufficient), m, sufficient))
        N = _integer_matrix(matrix(QQ, len(necessary), m, necessary))
        basis = A.right_kernel().basis_matrix()
        abscissa_basis = (basis*T.transpose()).row_space().basis_matrix()
        return dict(basis=basis, divisors=basis*G.transpose(),
                    abscissa_basis=abscissa_basis, constraints=A,
                    necessary_constraints=N, complete=A.rank() == N.rank())
