"""Geometric packet kernels at simple odd nodes without a global compositum."""

from sage.all import GF, NumberField, QQ, ZZ, lcm, matrix, vector

from .._galois import galois_trace_basis
from ..motivic import MotivicPacket


def _colliding_block(branch, p, node):
    """Return the size of a collision block and the partner map for size two."""
    degree = branch.degree()
    # There is an inertia transposition at the simple node. Its conjugates
    # give a graph whose connected components are cliques: the edge
    # transpositions generate the symmetric group on each component.
    # Transitivity makes their sizes equal, so prime degree gives one block.
    if ZZ(degree).is_prime():
        return int(degree), None

    field = NumberField(branch, 'a')
    a = field.gen()
    # Odd thickness gives a unique ramified quadratic factor at the node.
    primes = [q for q in field.primes_above(p) if (a-int(node)).valuation(q) > 0]
    if len(primes) != 1:
        return None
    prime, = primes
    factors = []
    for factor, _ in branch.change_ring(field).factor():
        value = factor(a)
        if value and value.valuation(prime) > 0:
            factors.append(factor)
    if len(factors) != 1:
        return None

    # This factor contains the other root specializing to the node.
    # Over Q(a) it consists of that root's orbit under the stabilizer
    # of a, i.e. all the other roots in a's clique.
    partner, = factors
    if partner.degree() == 1:
        return 2, branch.parent()((-partner[0]).list())
    return int(partner.degree()+1), None


def _node_collision(packet, p, disc):
    """Return the collision and its evidence when the prime can be handled here."""
    F = packet.F
    if p == 2 or F.leading_coefficient().valuation(p):
        return None
    k = GF(p)
    reduced = F.change_ring(k)
    repeated = reduced.gcd(reduced.derivative())
    if repeated.degree() != 1:
        return None
    node = repeated.roots(multiplicities=False)[0]
    branches = []
    for g, _ in packet.factors:
        if g is None:
            continue
        if any(c.denominator() % p == 0 for c in g):
            return None
        if g.change_ring(k)(node):
            continue
        if F % g:
            return None
        branches.append(g)
    thickness = int(disc.valuation(p))
    data = dict(method='simple_node_galois_pairs', exact=True, unresolved=[],
                node=str(node), node_thickness=thickness)
    if len(branches) == 1:
        branch, = branches
        if thickness % 2 == 0 or branch.change_ring(k) % repeated**2:
            return None
        try:
            block = _colliding_block(branch, p, node)
        except (ArithmeticError, NotImplementedError, ValueError):
            return None
        if block is None:
            return None
        size, partner = block
        data.update(branch_factor=str(branch), collision_block_size=size)
        return (branch, size, partner, None), data
    if len(branches) == 2:
        linear, branch = sorted(branches, key=lambda g: g.degree())
        if linear.degree() != 1:
            return None
        # The rational root is fixed. Conjugating the other root
        # gives every pair with this root and a root of its factor.
        data['branch_factors'] = [str(linear), str(branch)]
        return (branch, None, None, -linear[0]), data
    return None


def simple_nodal_kernel(packet):
    """Use single-node constraints, stopping if they already force rank zero.

    A colliding pair of Weierstrass points reduces to the midpoint of the
    loop, whose length is twice the valuation of the root difference.
    The projected boundary is a nonzero multiple of the sum of the
    coefficients on that pair, for any positive thickness.
    """
    if not packet.complete or any(c not in ZZ for c in list(packet.f)+list(packet.h)):
        return None
    disc = packet.F.discriminant()
    collisions, places, pending = [], {}, []
    for p in packet.curve.bad_primes():
        resolved = _node_collision(packet, p, disc)
        if resolved is None:
            try:
                stable = packet.curve.local_data(p).get('stable_type')
            except (ArithmeticError, NotImplementedError, ValueError):
                stable = None
            if stable and stable.split(',')[0] == 'I':
                places[str(p)] = dict(method='potential_good_reduction', exact=True, unresolved=[])
            else:
                pending.append(p)
        else:
            collision, data = resolved
            collisions.append(collision)
            places[str(p)] = data
    if pending and not collisions and len(packet.factors) > 1:
        return None

    from ._beilinson_divisors import _kernel_record
    from .geometric_torsion import encode_point

    points = packet.points()
    groups = {}
    for i, (x, _) in enumerate(points):
        groups.setdefault(x, []).append(i)
    n = len(groups)
    constraints = [[1]*n]
    for branch, size, partner, fixed_root in collisions:
        roots = {x: i for i, x in enumerate(groups) if x != 'infinity' and not branch(x)}
        if fixed_root is not None:
            fixed_index = next(i for i, x in enumerate(groups) if x == fixed_root)
            supports = [(fixed_index, i) for i in roots.values()]
        elif size >= 3:
            # Vanishing of all pair sums in a clique of size at least
            # three is equivalent to vanishing of each coefficient.
            supports = [(i,) for i in roots.values()]
        elif partner is None:
            supports = [tuple(roots.values())]
        else:
            mates = {i: roots[partner(x)] for x, i in roots.items()}
            supports = [(i, j) for i, j in mates.items() if i < j]
        constraints.extend([int(i in support) for i in range(n)] for support in supports)
    A = matrix(QQ, constraints)
    orbit_sum = matrix(QQ, [[int(i in group) for i in range(len(points))]
                            for group in groups.values()])

    def lift(v):
        return tuple(vector(QQ, [a/len(g) for a, g in zip(v, groups.values())])*orbit_sum)

    traces = galois_trace_basis(tuple(groups))
    candidates = (A*traces.transpose()).right_kernel().basis_matrix()*traces
    integral = [lift(v) for v in candidates]
    complete, upper = True, len(integral)
    solver = None
    if pending and integral:
        # Every integral divisor lies in the surviving span, so the
        # other points need never enter a compositum for local arithmetic.
        support = [i for i in range(len(points)) if any(row[i] for row in integral)]
        if len(support) == len(points):
            return None
        try:
            subpacket = MotivicPacket(packet.curve, [points[i] for i in support],
                                     torsion_multiples=2*lcm(order for _, order in packet.factors))
            solver = subpacket.integrality
            restricted = solver.kernel([[row[i] for i in support] for row in integral])
        except (ArithmeticError, NotImplementedError, ValueError):
            return None
        upper -= restricted['necessary_constraints'].rank()
        integral = list(restricted['basis']*matrix(QQ, integral))
        complete = restricted['complete']
        for place, data in solver.projected_local_data.items():
            if place.smallest_integer() in pending:
                places[str(place)] = dict(method=data.method, exact=data.exact,
                                          unresolved=list(data.unresolved), support_indices=support)
    else:
        for p in pending:
            places[str(p)] = dict(method='skipped_after_zero_kernel', exact=False,
                                  unresolved=list(range(n)))

    def check(D):
        if A*(orbit_sum*vector(QQ, D)):
            return False
        return True if solver is None else solver.check([D[i] for i in support])['q_integral']

    result = _kernel_record([encode_point(P) for P in points], integral, complete, orbit_sum,
                            lift, check, places, traces)
    result['rank_upper_bound'] = upper
    # The full packet contains distinct Weierstrass points, hence a
    # difference of order two. If ord(2P-kappa)=n, adjoining that order-two
    # difference makes the exponent exactly 2*lcm(n), including odd n.
    result.update(packet_torsion_order=str(2*lcm(n for _, n in packet.factors)),
                  packet_torsion_order_exact=True)
    return result
