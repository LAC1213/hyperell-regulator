r"""Certified embedded-tree cuts, analytic continuation and one infinity ray.

The rotation system is exact. At a finite vertex the boundary traverses the
clockwise sector between consecutive incident edges. At infinity it makes
one counterclockwise turn. These homotopies determine integer logarithm turns
and square-root signs without a sampled contour or floating-point unwrap.
"""
from contextlib import nullcontext
from functools import cmp_to_key

from sage.all import QQ, QQbar, ComplexField

from ._tree import minimum_spanning_tree, orient_to_root
from .balls import (CertificationError, digits_value, fields, narrow, principal,
                    rational, _model, edge_integrals, _assemble_walks)
from .diagnostics import _active, checkpoint


def _log(z, C):
    """Principal value at an exact point, choosing +pi on the negative axis."""
    w = C(z)
    if w.real() < 0:
        side = -1 if z.imag() < 0 else 1
        return (-w).log()+side*C.gen(0)*C.pi()
    return w.log()


def _sqrt(z, C):
    w = C(z)
    if w.real() < 0:
        side = -1 if z.imag() < 0 else 1
        return side*C.gen(0)*(-w).sqrt()
    return w.sqrt()


def _cyclic_order(directions):
    """Sort exact nonzero directions counterclockwise, starting at +Re."""
    def half(z):
        return 0 if z.imag() > 0 or (z.imag() == 0 and z.real() > 0) else 1
    def compare(a, b):
        u, v = directions[a], directions[b]
        if half(u) != half(v):
            return -1 if half(u) < half(v) else 1
        cross = (u.conjugate()*v).imag()
        if cross == 0:
            raise CertificationError('incident tree edges overlap')
        return -1 if cross > 0 else 1
    return sorted(directions, key=cmp_to_key(compare))


def _separated(z, a, b, c, d):
    for part in (lambda w: w.real(), lambda w: w.imag()):
        ab, cd = [part(z[i]) for i in (a, b)], [part(z[i]) for i in (c, d)]
        if max(v.upper() for v in ab) < min(v.lower() for v in cd):
            return True
        if max(v.upper() for v in cd) < min(v.lower() for v in ab):
            return True
    def cross(i, j, k):
        return ((z[j]-z[i]).conjugate()*(z[k]-z[i])).imag()
    u, v = cross(a, b, c), cross(a, b, d)
    if (u > 0 and v > 0) or (u < 0 and v < 0):
        return True
    u, v = cross(c, d, a), cross(c, d, b)
    return (u > 0 and v > 0) or (u < 0 and v < 0)


def _validate_embedding(nodes, edges, C):
    """Verify every geometric property used by the rotation system."""
    z = list(map(C, nodes))
    for i, (a, b) in enumerate(edges):
        for j, point in enumerate(nodes):
            if j in (a, b):
                continue
            position = (z[j]-z[a])/(z[b]-z[a])
            if not position.imag().contains_zero() or position.real() < 0 or position.real() > 1:
                continue
            exact = (point-nodes[a])/(nodes[b]-nodes[a])
            if exact.imag() == 0 and 0 <= exact.real() <= 1:
                raise CertificationError('tree edge passes through another cut vertex')
        for c, d in edges[:i]:
            if len({a, b, c, d}) == 4 and not _separated(z, a, b, c, d):
                raise CertificationError('could not certify disjoint tree edges')


def _precision_tree(nodes, C):
    P = ComplexField(C.precision())
    z = list(map(P, nodes))
    parent = list(range(len(z)))
    def find(i):
        while parent[i] != i:
            i = parent[i]
        return i
    result = []
    for _, a, b in sorted((abs(z[a]-z[b]), a, b) for a in range(len(z)) for b in range(a)):
        i, j = find(a), find(b)
        if i != j:
            parent[i] = j
            result.append((a, b))
    return result


class TreeCut:
    """Exact vertices, certified planar embedding and its boundary rotation."""
    def __init__(self, F, centers, digits):
        self.F, self.centers = F, centers
        self.roots = F.roots(QQbar, multiplicities=False)
        self.nodes = list(self.roots)
        for c in centers:
            if QQbar(c) not in self.nodes:
                self.nodes.append(QQbar(c))
        C, R = fields(digits)
        try:
            edges = minimum_spanning_tree(self.nodes)
            _validate_embedding(self.nodes, edges, C)
        except (ValueError, CertificationError):
            edges = _precision_tree(self.nodes, C)
            _validate_embedding(self.nodes, edges, C)
        root = max(range(len(self.nodes)), key=lambda i: (self.nodes[i].real(), self.nodes[i].imag()))
        self.edges = orient_to_root(len(self.nodes), edges, root)
        self.ray_vertex, self.direction = root, QQbar(1)
        # All finite edges are in Re(x)<=Re(root); the positive horizontal
        # ray therefore meets the tree only at root. This also makes its
        # regularized square-root factors stay in the right half-plane.
        z = list(map(C, self.nodes))
        width = max(w.real().upper() for w in z)-min(w.real().lower() for w in z)
        height = max(w.imag().upper() for w in z)-min(w.imag().lower() for w in z)
        self.span = QQ(max(1, (width+height).ceil()))
        self.references = [(self.nodes[a]+self.nodes[b])/2 if b is not None else self.nodes[a]+self.span
                           for a, b in self.edges]
        self.labels = list(range(len(self.roots)))+[None]*(len(self.nodes)-len(self.roots))
        self.g = (F.degree()-1)//2
        infinity = len(self.nodes)
        neighbors = [{} for _ in range(infinity+1)]
        lookup = {}
        for i, (a, b) in enumerate(self.edges):
            if b is None:
                b = infinity
                delta = self.direction
            else:
                delta = self.nodes[b]-self.nodes[a]
            neighbors[a][b], neighbors[b][a] = delta, -delta
            lookup[a, b], lookup[b, a] = (i, 1), (i, -1)
        self.cyclic_orders = [_cyclic_order(row) for row in neighbors]
        leaf = next(i for i in range(infinity) if len(neighbors[i]) == 1)
        first = (leaf, self.cyclic_orders[leaf][0])
        a, b = first
        self.boundary_word, self.banks = [], []
        while True:
            self.boundary_word.append(lookup[a, b])
            self.banks.append((a, b))
            row = self.cyclic_orders[b]
            a, b = b, row[(row.index(a)-1) % len(row)]
            if (a, b) == first:
                break
            if len(self.banks) > 2*len(self.edges):
                raise CertificationError('tree boundary does not close')
        if len(self.banks) != 2*len(self.edges) or len(set(self.boundary_word)) != len(self.banks):
            raise CertificationError('tree boundary omits or repeats a bank')
        self.turns = self._turns(C, R)
        parents = dict(self.edges)
        for i, (a, _) in enumerate(self.edges):
            for node in range(len(self.nodes)):
                ancestor = node
                while ancestor is not None and ancestor != a:
                    ancestor = parents[ancestor]
                expected = -1 if ancestor == a else 0
                if self.turns[i, 1][node]-self.turns[i, -1][node] != expected:
                    raise CertificationError('logarithm jump disagrees with tree component')

    def _turns(self, C, R):
        angles = [[_log(ref-c, C).imag() for c in self.nodes] for ref in self.references]
        turns, current = {}, [0]*len(self.nodes)
        pi2 = 2*R.pi()
        for k, (i, side) in enumerate(self.boundary_word):
            turns[i, side] = list(current)
            j, _ = self.boundary_word[(k+1) % len(self.banks)]
            vertex = self.banks[k][1]
            if vertex == len(self.nodes):
                current = [n+1 for n in current]
                continue
            v, ri, rj = self.nodes[vertex], self.references[i], self.references[j]
            for c_index, c in enumerate(self.nodes):
                if c_index == vertex:
                    ratio = (rj-v)/(ri-v)
                    angle = _log(ratio, C).imag()
                    # Clockwise exterior sector, a full turn at a leaf.
                    gap = -angle if ratio.imag() < 0 else pi2-angle
                    increment = -gap
                else:
                    increment = (_log((v-c)/(ri-c), C).imag()
                                 +_log((rj-c)/(v-c), C).imag())
                # This expression is an integer by the straight-segment /
                # vertex-sector homotopy. Isolate that integer with balls.
                delta = (angles[i][c_index]+increment-angles[j][c_index])/pi2
                try:
                    lo, hi = delta.lower().ceil(), delta.upper().floor()
                    if lo != hi:
                        raise ValueError('turn interval has no unique integer')
                    change = int(lo)
                except ValueError as exc:
                    raise CertificationError('could not isolate logarithm turn') from exc
                current[c_index] += change
        if any(current):
            raise CertificationError('tree logarithm continuation does not close')
        return turns


def _sheet_sign(ratio):
    """Isolate +1 or -1 using the identity ratio^2=1."""
    if not ratio.imag().contains_zero() or (not ratio.real() > 0 and not ratio.real() < 0):
        raise CertificationError('could not certify reference sheet sign')
    return -1 if ratio.real() < 0 else 1


def _ray_tail(cut, digits):
    r"""Validated reciprocal tail B -> infinity, with v=u^4.

    x=z+D/u^8 gives root factors (D+(z-r)*u^8)/(B-r), all
    in the right half-plane along [0,1]. Endpoint powers are nonnegative;
    the only logarithmic singularity is an explicit -8*log(u).
    """
    C, R = fields(digits)
    z_exact = cut.nodes[cut.ray_vertex]
    B = cut.references[-1]
    D, z = C(cut.span), C(z_exact)
    roots = [(C(z_exact-r), C(B-r)) for r in cut.roots]
    logs = [(C(z_exact-c), C(B-c), _log(B-c,C)) for c in cut.centers]
    yB = _sqrt(QQbar(cut.F.leading_coefficient()), C)
    for r in cut.roots:
        yB *= _sqrt(B-r, C)
    def differential(u, analytic, j):
        h = u**8
        X = D+z*h
        y = yB
        for a, b in roots:
            y *= principal((D+a*h)/b, analytic, 'sqrt')
        return 8*D*u**(4*(cut.F.degree()-2-2*j)-1)*X**j/y
    def logarithm(u, analytic, record, omit=False):
        a,b,seed = record
        value = seed+principal((D+a*u**8)/b,analytic,'log')
        return value if omit else value-8*principal(u,analytic,'log')
    eps, tolerance = QQ(1)/10**((digits+23)//4), R(10)**(-digits-5)
    endpoint = C(R(eps/2).add_error(R(eps/2)))
    def integrate(j, record=None):
        bound = R(1) if record is None else (abs(logarithm(endpoint,False,record,True))
                                             +8*(1-R(eps).log()))
        error = R(eps)*abs(differential(endpoint,False,j))*bound
        def evaluate(u,analytic):
            checkpoint('quadrature_nodes',phase='certified tree ray tail')
            value = differential(u,analytic,j)
            return value if record is None else value*logarithm(u,analytic,record)
        value = C.integral(evaluate,eps,1,abs_tol=tolerance,rel_tol=tolerance).add_error(error)
        if not narrow(value,digits):
            raise CertificationError('tree ray tail did not achieve requested radius')
        return value
    return dict(periods=[integrate(j) for j in range(cut.g)],
                logarithmic=[[integrate(j,record) for j in range(cut.g)] for record in logs],
                reference_y=yB, reference_logs=[record[2] for record in logs])


def _ray(cut,digits):
    """Join a regularized finite segment and a certified reciprocal tail."""
    z,B = cut.nodes[cut.ray_vertex],cut.references[-1]
    first = edge_integrals(cut.F,z,B,centers=cut.centers,powers=range(cut.g),digits=digits+3)
    tail = _ray_tail(cut,digits+3)
    C = tail['reference_y'].parent()
    mid = (z+B)/2
    # Continue the finite kernel's midpoint sheet to B. There are no
    # intervening cut vertices; every ratio has positive real part.
    yend = C(first['reference_y'])*C(2).sqrt()**int(z in cut.roots)
    for r in cut.roots:
        if r != z:
            yend *= (C(B-r)/C(mid-r)).sqrt()
    sign = _sheet_sign(yend/tail['reference_y'])
    tail['periods'] = [sign*a+b for a,b in zip(first['periods'],tail['periods'])]
    tail['logarithmic'] = [[sign*a+b for a,b in zip(left,right)]
                           for left,right in zip(first['logarithmic'],tail['logarithmic'])]
    if not all(narrow(v,digits) for v in tail['periods']+sum(tail['logarithmic'],[])):
        raise CertificationError('complete tree ray did not achieve requested radius')
    return tail


def compute_tree(F, centers=(), prec=15, *, details=None):
    digits = digits_value(prec)
    F, centers = _model(F), list(dict.fromkeys(rational(c) for c in centers))
    cut = TreeCut(F, centers, digits)
    for guard in (10, 25, 50):
        target = digits+guard
        C, R = fields(target)
        results = []
        for i, (a, b) in enumerate(cut.edges):
            mon = _active.get()
            with mon.edge('certified tree finite' if b is not None else 'certified tree ray',
                          degree=2, precision=target) if mon else nullcontext():
                results.append(edge_integrals(F, cut.nodes[a], cut.nodes[b], centers=centers,
                                              powers=range(cut.g), digits=target)
                               if b is not None else _ray(cut, target))
        upper, lower = [], []
        for i, ref in enumerate(cut.references):
            base = _sqrt(QQbar(F.leading_coefficient()), C)
            for r in cut.roots:
                base *= _sqrt(ref-r, C)
            ratio = base/results[i]['reference_y']
            seed = int(_sheet_sign(ratio) < 0)
            for side, output in ((1, upper), (-1, lower)):
                sign = (seed+sum(cut.turns[i, side][:len(cut.roots)])) % 2
                output.append([sign, 1-sign])
        cut.upper, cut.lower = upper, lower
        cut.walks = [[(i, (upper if side > 0 else lower)[i][sheet], side)
                      for i, side in cut.boundary_word] for sheet in (0, 1)]
        u = [[[row['periods'][a], -row['periods'][a]] for row in results] for a in range(cut.g)]
        per = _assemble_walks(C, R, cut.g, cut.walks, u, u, [C(0)]*len(results))
        cols, all_jumps, all_logs = [], [], []
        for k, c in enumerate(centers):
            index = cut.nodes.index(QQbar(c))
            jumps = [2*C.gen(0)*C.pi()*(cut.turns[i, 1][index]-cut.turns[i, -1][index])
                     for i in range(len(results))]
            shifts = [2*C.gen(0)*C.pi()*cut.turns[i, -1][index] for i in range(len(results))]
            ell = [[[row['logarithmic'][k][a]+shifts[i]*row['periods'][a],
                     -(row['logarithmic'][k][a]+shifts[i]*row['periods'][a])]
                    for i, row in enumerate(results)] for a in range(cut.g)]
            cols.append(_assemble_walks(C, R, cut.g, cut.walks, u, ell, jumps))
            all_jumps.append(jumps)
            all_logs.append([row['reference_logs'][k]+shifts[i] for i, row in enumerate(results)])
        cut.jumps, cut.logs = all_jumps, all_logs
        if all(narrow(v, digits) for v in per+sum(cols, [])):
            if details is not None:
                details.update(certified=True, cut_kind='tree', cut=cut,
                               exact_vertices=cut.nodes, exact_direction=cut.direction,
                               working_bits=C.precision())
            return per, cols
    raise CertificationError('assembled tree regulator did not achieve requested radius')
