"""Construct split stable reduction data and its marked dual graph."""

from dataclasses import dataclass
from functools import cached_property

from sage.all import NumberField, PolynomialRing, QQ

from .graph_boundary import FibreGraph
from .number_field_integrality import _LocalSections, _curve_polynomials


@dataclass
class StableModel:
    """Special-fibre charts, node thicknesses and specialization of marked points.

    ``graph`` retains the vertices needed to specialize sections.
    ``stable_graph`` contracts the unstable rational components after
    forgetting the marks. The component equations describe the
    normalizations of the components indexed by ``graph``.
    """

    field: object
    prime: object
    graph: FibreGraph
    specializations: tuple
    component_equations: tuple = ()
    charts: tuple = ()
    residue_extension_degree: int = 1
    length_normalization: str = 'uniformizer'
    node_coordinates: tuple = ()

    @cached_property
    def stabilization(self):
        return self.graph.stabilize()

    @property
    def stable_graph(self):
        return self.stabilization[0]

    @property
    def stable_component_equations(self):
        if not self.component_equations:
            return ()
        return tuple(self.component_equations[i] for i in self.stabilization[1])

    def boundary(self, coefficients):
        """Compute the special-fibre boundary of the marked divisor."""
        coeffs = tuple(map(QQ, coefficients))
        if len(coeffs) != len(self.specializations) or sum(coeffs):
            raise ValueError('one coefficient per mark and total degree zero are required')
        d = [QQ(0)]*len(self.graph.genera)
        for a, positions in zip(coeffs, self.specializations):
            if not a:
                continue
            if positions is None:
                raise NotImplementedError('a point in the divisor has unresolved specialization')
            for i in positions:
                d[i] += a/len(positions)
        return self.graph.boundary(d)


def _extend(field, prime, polynomial, name):
    relative_name = name+'0' if field.variable_name() == name else name
    L = field.extension(polynomial, relative_name)
    K = L.absolute_field(name)
    _, to_K = K.structure()
    embed = lambda a: to_K(L(a))
    q = K.ideal([embed(g) for g in prime.gens()]).prime_factors()[0]
    return K, q, embed, to_K(L.gen())


def _unramified_quadratic(field, prime, max_degree):
    if max_degree is not None and 2*field.degree() > max_degree:
        raise NotImplementedError('residue splitting exceeds the requested extension degree')
    k = prime.residue_field()
    lift = k.lift_map()
    T = PolynomialRing(field, 'z')
    modulus = T([field(lift(c)) for c in k['z'].irreducible_element(2)])
    K, q, embed, _ = _extend(field, prime, modulus, 'z')
    if (q.ramification_index() != prime.ramification_index() or
            q.residue_class_degree() != 2*prime.residue_class_degree()):
        raise ArithmeticError('the model residue extension did not split as expected')
    return K, q, embed


def _split_branches(f, h, prime, groups, max_degree):
    K = f.base_ring()
    while True:
        irreducible = next((g for g, _ in (h*h+4*f).factor() if g.degree() > 1), None)
        if irreducible is None:
            return f, h, prime, groups
        if max_degree is not None and K.degree()*irreducible.degree() > max_degree:
            raise NotImplementedError('branch splitting exceeds the requested extension degree')
        K, prime, embed, _ = _extend(K, prime, irreducible.monic(), 'b')
        R = PolynomialRing(K, 'x')
        f, h = [R([embed(c) for c in g]) for g in (f, h)]
        groups = tuple(tuple((x if x == 'infinity' else embed(x), embed(y)) for x, y in g)
                       for g in groups)


def _odd_chart(f, h, groups):
    F = h*h+4*f
    if F.degree() == 5:
        return F, tuple(tuple(P[0] for P in g) for g in groups), None
    root = F.roots(multiplicities=False)[0]
    x = F.parent().gen()
    shifted = F(x+root)
    abscissae = tuple(tuple(0 if a == 'infinity' else
                            'infinity' if a == root else 1/(a-root)
                            for a, _ in g) for g in groups)
    return F.parent()([shifted[6-i] for i in range(7)]), abscissae, root


def _disc_tree(abscissae, valuation, make_vertex):
    """Stable tree of the marked projective line, rooted at the infinity mark."""
    positions = {}

    def visit(xs):
        centre = xs[0]
        depth = min(valuation(x-y) for i, x in enumerate(xs) for y in xs[:i])
        parts = []
        for x in xs:
            for part in parts:
                if x == part[0] or valuation(x-part[0]) > depth:
                    part.append(x)
                    break
            else:
                parts.append([x])
        vertex = make_vertex(centre, depth, xs, parts)
        for part in parts:
            if len(part) == 1:
                positions[part[0]] = vertex
            else:
                child, child_depth = visit(part)
                vertex['children'].append((child, child_depth-depth))
        return vertex, depth

    root, _ = visit(abscissae)
    positions['infinity'] = root
    return positions


def _tame_model(f, h, prime, groups, max_degree):
    f, h, prime, groups = _split_branches(f, h, prime, groups, max_degree)
    K = f.base_ring()
    F, marks, reciprocal_centre = _odd_chart(f, h, groups)
    roots = tuple(F.roots(multiplicities=False))
    xs = list(dict.fromkeys([*roots, *(x for g in marks for x in g if x != 'infinity')]))
    val = lambda x: x.valuation(prime)
    # A single quadratic ramification makes every valuation of F at a
    # vertex even and every annulus branched at its ends divisible by two.
    T = PolynomialRing(K, 's')
    if max_degree is not None and 2*K.degree() > max_degree:
        raise NotImplementedError('ramification exceeds the requested extension degree')
    K, q, embed, pi = _extend(K, prime, T.gen()**2-K.uniformizer(prime), 's')
    if q.ramification_index() != 2*prime.ramification_index():
        raise ArithmeticError('the quadratic model extension did not ramify')
    k = q.residue_field()
    red = k.reduction_map()
    S = PolynomialRing(k, 't')
    t = S.gen()
    genera, edges, equations, charts = [], [], [], []
    need_unramified = False

    def make_vertex(centre, depth, disk, parts):
        nonlocal need_unramified
        inside = [r for r in roots if r in disk]
        odd = [part[0] for part in parts if sum(r in part for r in inside) % 2]
        # Remove the square factors on the residual double cover.
        nu = val(F.leading_coefficient())
        unit = embed(F.leading_coefficient())/pi**(2*nu)
        for r in roots:
            n = depth if r in inside else min(depth, val(centre-r))
            # For a root inside the disc the leading t coefficient is 1.
            if r not in inside:
                unit *= embed(centre-r)/pi**(2*n)
            nu += n
        g = S(red(unit))
        for a in odd:
            g *= t-red(embed(a-centre)/pi**(2*depth))
        branch_count = len(odd)+len(inside)%2
        vertices = tuple(range(len(genera), len(genera)+(1 if branch_count else 2)))
        if branch_count:
            genera.append((branch_count-2)//2)
            equations.append(g)
        else:
            genera.extend([0, 0])
            equations.extend([None, None])
        # The two branches over an unbranched node are defined over at
        # most the quadratic residue extension, since the base is split.
        need_unramified |= len(inside)%2 == 0 and not g.leading_coefficient().is_square()
        need_unramified |= any(
            sum(r in part for r in inside)%2 == 0
            and not g(red(embed(part[0]-centre)/pi**(2*depth))).is_square()
            for part in parts if len(part)>1)
        chart = dict(vertices=vertices, centre=embed(centre), scale=pi**(2*depth),
                     y_scale=pi**nu, residual_polynomial=g)
        charts.append(chart)
        return dict(vertices=vertices, odd=len(inside)%2, children=[], chart=chart)

    positions = _disc_tree(xs, val, make_vertex)

    reciprocal_centre = None if reciprocal_centre is None else embed(reciprocal_centre)
    specializations = []
    for group in marks:
        v = positions[group[0]]['vertices']
        # Two conjugate sections lie on opposite sheets over a split
        # component. A singleton on a split chart needs sheet gluing.
        specializations.append(v if len(v) == 1 or len(group) == 2 else None)
    if need_unramified:
        K, q, embed = _unramified_quadratic(K, q, max_degree)
        lift = k.lift_map()
        red = q.residue_field().reduction_map()
        S = PolynomialRing(q.residue_field(), 't')

        def change(g):
            return S([red(embed(lift(c))) for c in g])

        equations = [None if g is None else change(g) for g in equations]
        for chart in charts:
            for name in ('centre', 'scale', 'y_scale'):
                chart[name] = embed(chart[name])
            chart['residual_polynomial'] = change(chart['residual_polynomial'])
        if reciprocal_centre is not None:
            reciprocal_centre = embed(reciprocal_centre)
    charts[0]['reciprocal_centre'] = reciprocal_centre
    nodes = []

    def glue(v):
        parent = v['chart']
        for child, length in v['children']:
            chart = child['chart']
            a, b = v['vertices'], child['vertices']
            t = red((chart['centre']-parent['centre'])/parent['scale'])
            if child['odd']:
                edges.append((a[0], b[0], length))
                nodes.append(((t, S.base_ring()(0)), ('infinity', S.base_ring()(0))))
            else:
                # Each matching of branches differs by the involution on
                # the subtree below this edge of the marked base tree.
                y = parent['residual_polynomial'](t).sqrt()
                eta = chart['residual_polynomial'].leading_coefficient().sqrt()
                for i, sign in enumerate((1, -1)):
                    edges.append((a[i % len(a)], b[i % len(b)], 2*length))
                    nodes.append(((t, sign*y), ('infinity', sign*eta)))
            glue(child)

    glue(positions['infinity'])
    graph = FibreGraph(genera, edges)
    if sum(genera)+graph.cycles.nrows() != 2:
        raise ArithmeticError('the constructed double cover has wrong arithmetic genus')
    return StableModel(K, q, graph, tuple(specializations), tuple(equations), tuple(charts),
                       2 if need_unramified else 1, node_coordinates=tuple(nodes))


def _dyadic_good_model(local, mark_count, max_degree):
    """Search for a smooth chart by a square completion and disc rescaling."""
    from itertools import product

    K, q, k = local.K, local.prime, local.k
    R = local.f.parent()
    T = PolynomialRing(K, 's')
    x, s = R.gen(), T.gen()
    for reciprocal, f, h in ((False, local.f, local.h), (True, local.inf_f, local.inf_h)):
        reduced_h = local.reduce_poly(h)
        if not reduced_h:
            continue
        for root in reduced_h.roots(multiplicities=False):
            c = K(local.lift(root))
            f0, h0 = f(x+c), h(x+c)
            # Lifts differing by a multiple of the residue characteristic
            # can give distinct square completions at a finer scale.
            for terms in product((-1, 0, 1), repeat=4):
                a = R(terms)
                f1, h1 = f0-a*a-h0*a, h0+2*a
                if not f1 or not h1:
                    continue
                j = min(i for i in range(h1.degree()+1) if h1[i] and not h1[i].valuation(q))
                if not j:
                    continue
                r = min(b.valuation(q)/QQ(2*j-i) for i, b in enumerate(f1)
                        if b and i < 2*j) if any(f1[i] for i in range(2*j)) else None
                if r is None or r <= 0:
                    continue
                if any(b.valuation(q)+(i-2*j)*r < 0 for i, b in enumerate(f1) if b):
                    continue
                if any(b.valuation(q)+(i-j)*r < 0 for i, b in enumerate(h1) if b):
                    continue
                degree = r.denominator()
                if max_degree is not None and K.degree()*degree > max_degree:
                    continue
                pi = K.uniformizer(q)
                if degree == 1:
                    L, prime, embed, rho = K, q, K, pi**r.numerator()
                else:
                    L, prime, embed, rho = _extend(K, q, s**degree-pi**r.numerator(), 's')
                S = PolynomialRing(L, 'u')
                g = S([embed(b)*rho**(i-2*j) for i, b in enumerate(f1)])
                h_new = S([embed(b)*rho**(i-j) for i, b in enumerate(h1)])
                ext = _LocalSections(g, h_new, prime, max_blowups=local.max_blowups)
                from ._projected_reduction import nodal_model

                try:
                    model = nodal_model(ext, [], [])
                except NotImplementedError:
                    continue
                if model.graph.genera == (2,):
                    return StableModel(L, prime, model.graph, ((0,),)*mark_count,
                                       charts=(*model.charts, dict(reciprocal=reciprocal,
                                           centre=embed(c), shift=S([embed(b) for b in a]),
                                           x_scale=rho, y_scale=rho**j, f=g, h=h_new)))
    raise NotImplementedError('no smooth dyadic chart found by square completion')


def _dyadic_model(C, local, points, groups, max_degree):
    from ._projected_reduction import nodal_model

    try:
        return nodal_model(local, points, groups)
    except NotImplementedError:
        pass
    stable = C.local_data(2).get('stable_type') if hasattr(C, 'local_data') else None
    stable = stable.split(',')[0] if stable else None
    if stable in (None, 'I'):
        return _dyadic_good_model(local, len(groups), max_degree)
    from ._dyadic_projected import dyadic_constraints

    constraints = dyadic_constraints(local, points, groups, stable)
    if 'model' not in constraints.evidence:
        raise NotImplementedError('only the projected boundary is available for this wild fibre')
    return constraints.evidence['model']


def stable_model(C, prime, *, field=QQ, points=(), max_degree=None, max_blowups=32):
    """Construct the stable special fibre over a finite extension at a place.

    Supply a rational prime over QQ or a prime ideal of an absolute
    number field containing the marked points. Points have the coordinates
    of the original model, including ``('infinity', eta)``. The odd
    residue characteristic construction uses the split marked double
    cover. At two, the supported models are reduced nodal equations,
    verified smooth disc charts, and the verified compact-type charts.
    """
    if field is QQ:
        R = PolynomialRing(QQ, 'a')
        field = NumberField(R.gen(), 'a')
    if prime in QQ:
        prime, = field.primes_above(prime)
    f, h = _curve_polynomials(C, field)
    groups = tuple(((x if x == 'infinity' else field(x), field(y)),) for x, y in points)
    for ((x, y),) in groups:
        value, linear = (f[6], h[3]) if x == 'infinity' else (f(x), h(x))
        if y*y+linear*y != value:
            raise ValueError('a marked point is off the curve')
    coefficients = list(f)+list(h)
    exponent = max([0, *(-c.valuation(prime) for c in coefficients if c)])
    if exponent:
        scale = field.uniformizer(prime)**exponent
        f, h = scale**2*f, scale*h
        groups = tuple(((x, scale*y),) for ((x, y),) in groups)
    local = _LocalSections(f, h, prime, max_blowups=max_blowups)
    marks = [g[0] for g in groups]
    indices = [(i,) for i in range(len(groups))]
    if int(prime.smallest_integer()) != 2:
        from ._projected_reduction import nodal_model

        try:
            model = nodal_model(local, marks, indices)
            if not model.graph.edges:
                return model
        except (NotImplementedError, ArithmeticError):
            pass
        return _tame_model(f, h, prime, groups, max_degree)
    return _dyadic_model(C, local, marks, indices, max_degree)
