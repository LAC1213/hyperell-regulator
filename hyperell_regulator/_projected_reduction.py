"""Exact local graph constraints for projected torsion-packet classes."""

from dataclasses import dataclass

from sage.all import PolynomialRing, QQ, lcm, matrix, vector

from .graph_boundary import FibreGraph
from .number_field_integrality import _LocalSections, _factors, _value


@dataclass
class LocalConstraints:
    method: str
    matrix: object
    exact: bool
    unresolved: tuple = ()
    evidence: object = None


def component_constraints(groups, labels, method='multidegree'):
    """A sufficient constraint on all points whose components are known."""
    components = sorted({labels[i] for group in groups for i in group}-{None})
    rows = [[QQ(sum(labels[i] == c for i in group))/len(group)
             for group in groups] for c in components]
    unknown = tuple(j for j, group in enumerate(groups) if any(labels[i] is None for i in group))
    return LocalConstraints(method, matrix(QQ, len(rows), len(groups), rows),
                            False, unknown, dict(components=labels))


def node_thickness(local, f, h, a, b):
    """Certify the critical value of the local equation using Newton lifting."""
    K = local.K
    x, y = K(local.lift(a)), K(local.lift(b))
    for _ in range(8):
        H = matrix(K, [[h.derivative(2)(x)*y-f.derivative(2)(x), h.derivative()(x)],
                       [h.derivative()(x), 2]])
        if local.valuation(H.det()) != 0:
            raise NotImplementedError('the special fibre has a singularity other than a node')
        gradient = vector(K, [h.derivative()(x)*y-f.derivative()(x), 2*y+h(x)])
        value = y*y+h(x)*y-f(x)
        # The critical point agrees to precision v(gradient). The value
        # is then determined to precision at least 2*v(gradient).
        if value and local.valuation(value) < 2*min(map(local.valuation, gradient)):
            return int(local.valuation(value))
        x, y = vector(K, [x, y])-H.solve_right(gradient)
    raise ArithmeticError('could not certify the thickness of an ordinary node')


def _singular_abscissae(f, h):
    if f.base_ring().characteristic() == 2:
        return h.gcd(h.derivative()**2*f-f.derivative()**2)
    discriminant = h*h+4*f
    return discriminant.gcd(discriminant.derivative())


def _split_node_field(local, points):
    f, h = map(local.reduce_poly, (local.f, local.h))
    singular = _singular_abscissae(f, h)
    if singular.is_zero():
        raise NotImplementedError('the fibre is singular along a component')
    degree = lcm([g.degree() for g, _ in singular.factor()])
    if degree % 2 and len(local.main) == 1:
        k2 = local.k.extension(2, 'b')
        if len(_factors(f.change_ring(k2), h.change_ring(k2))) == 2:
            degree *= 2
    local, points = _extend_residue(local, points, degree)
    # Even when the singular point itself is rational, the two tangent
    # lines can be exchanged by Frobenius.
    k = local.k
    T = PolynomialRing(k, 'w')
    w = T.gen()
    for f, h, infinity in ((local.f, local.h, False), (local.inf_f, local.inf_h, True)):
        f, h = map(local.reduce_poly, (f, h))
        singular = _singular_abscissae(f, h)
        roots = [k(0)] if infinity and not singular(0) else (
            [] if infinity else singular.roots(multiplicities=False))
        x = f.parent().gen()
        for a in roots:
            b = f(a).sqrt() if k.characteristic() == 2 else -h(a)/2
            g, j = f(x+a), h(x+a)
            if not (w*w+j[1]*w+b*j[2]-g[2]).roots():
                return _extend_residue(local, points, 2)
    return local, points


def _extend_residue(local, points, degree):
    if degree == 1:
        return local, points
    k = local.k
    R = PolynomialRing(local.K, 'z')
    modulus = k['z'].irreducible_element(degree)
    L = local.K.extension(R([local.K(local.lift(c)) for c in modulus]), 'z')
    K = L.absolute_field('a')
    _, to_K = K.structure()

    def embed(a):
        return to_K(L(a))

    prime, = K.ideal([embed(a) for a in local.prime.gens()]).prime_factors()
    if (prime.ramification_index() != local.prime.ramification_index()
            or prime.residue_class_degree() != degree*local.prime.residue_class_degree()):
        raise ArithmeticError('the node splitting extension is not unramified')
    S = PolynomialRing(K, 'x')
    split = _LocalSections(S([embed(a) for a in local.f]), S([embed(a) for a in local.h]),
                           prime, max_blowups=local.max_blowups)
    return split, [(x if x == 'infinity' else embed(x), embed(y)) for x, y in points]


def _nodes(local):
    nodes = {}
    count = len(local.main)
    for patch, f, h, factors in (
        ('finite', local.f, local.h, local.main),
        ('infinity', local.inf_f, local.inf_h, local.main_infinity),
    ):
        fb, hb = map(local.reduce_poly, (f, h))
        singular = _singular_abscissae(fb, hb)
        if not singular:
            raise NotImplementedError('the fibre is singular along a component')
        xs = (singular.roots(multiplicities=False) if patch == 'finite' else
              [local.k(0)] if singular(0) == 0 else [])
        for a in xs:
            b = fb(a).sqrt() if local.k.characteristic() == 2 else -hb(a)/2
            n = node_thickness(local, f, h, a, b)
            if count == 1:
                ends = [0, 0]
            else:
                ends = [i for i, g in enumerate(factors) if _value(g, a, b) == 0]
                if len(ends) != 2:
                    raise NotImplementedError('the node does not join two distinct components')
            nodes[patch, str(a), str(b)] = (*ends, n)
    return nodes


def _section_positions(reductions, paths):
    labels = [q[0] for q in reductions]
    if any(label is None for label in labels):
        return None
    if all(label.startswith('main:') for label in labels):
        return tuple(int(label.split(':')[-1]) for label in labels)
    keys = [(steps[0]['patch'], steps[0]['residue_x'], steps[0]['residue_y'])
            for _, _, steps in reductions]
    depths = {len(steps)-1 for _, _, steps in reductions}
    if len(set(keys)) != 1 or len(depths) != 1 or keys[0] not in paths:
        return None
    depth, = depths
    path = paths[keys[0]]
    if 2*depth > len(path)-1:
        return None
    if path[0] != path[-1] and 2*depth != len(path)-1 and len(set(labels)) != 2:
        return None
    # Reflection on a loop preserves the boundary signature. On an
    # edge between distinct vertices the conjugates give both positions.
    return (path[depth], path[-1-depth])


def nodal_model(local, points, groups):
    from .stable_model import StableModel
    local, points = _split_node_field(local, points)
    nodes = _nodes(local)
    count = len(local.main)
    if count == 1 and len(nodes) <= 2:
        genera = [2-len(nodes)]
    elif count == 2 and len(nodes) == 3 and all(len(g) == 2 for g in local.main):
        genera = [0, 0]
    else:
        raise NotImplementedError('unsupported nodal genus-two fibre')
    paths, edges = {}, []
    for key, (a, b, length) in nodes.items():
        start = len(genera)
        path = [a, *range(start, start+length-1), b]
        paths[key] = path
        genera.extend([0]*(length-1))
        edges.extend((i, j, 1) for i, j in zip(path, path[1:]))
    reductions = [local.point(P) for P in points]
    marks = tuple(_section_positions([reductions[i] for i in group], paths)
                  for group in groups)
    charts = (dict(f=local.reduce_poly(local.f), h=local.reduce_poly(local.h),
                   nodes=nodes, paths=paths),)
    return StableModel(local.K, local.prime, FibreGraph(genera, edges), marks, charts=charts)


def graph_constraints(model, method):
    weights = model.graph.vertex_boundaries
    columns = [sum((weights[i] for i in positions))/len(positions) if positions else
               weights[0] for positions in model.specializations]
    unknown = tuple(i for i, col in enumerate(model.specializations) if col is None)
    A = matrix(QQ, len(columns), len(weights[0]), columns).transpose()
    return LocalConstraints(method, A, True, unknown, dict(model=model))


def semistable_constraints(local, points, groups):
    model = nodal_model(local, points, groups)
    out = graph_constraints(model, 'semistable')
    out.evidence.update(model.charts[0])
    return out
