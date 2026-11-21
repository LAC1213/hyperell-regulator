"""The cohomological boundary of torsion-packet classes on a fibre graph."""

from dataclasses import dataclass
from functools import cached_property

from sage.all import QQ, identity_matrix, matrix, vector


@dataclass(frozen=True)
class GraphBoundary:
    """Coordinates of the boundary functional on special-fibre cohomology.

    The potential is normalized at vertex zero. The coordinates pair with
    the diagonal mixed classes on positive-genus components and with
    tensors of graph cycles, respectively. The latter pairing has the
    opposite sign to ``cycle_form`` in the product blowup convention.
    ``projected`` gives the class modulo the diagonal Tate direction.
    """

    potential: object
    component_potentials: object
    cycle_form: object
    tate_signature: object
    projected: object


class FibreGraph:
    """A connected metric dual graph with the genus of each normalization.

    Edges are triples ``(tail, head, length)``. Loops and parallel edges
    are allowed. All arithmetic is rational. In a regular model edges
    have unit length. Contracting a chain adds its lengths.
    """

    def __init__(self, genera, edges):
        self.genera = tuple(map(int, genera))
        self.edges = tuple((int(i), int(j), QQ(n)) for i, j, n in edges)
        m = len(self.genera)
        if not m or any(g < 0 for g in self.genera):
            raise ValueError('the graph needs vertices with nonnegative genera')
        if any(not (0 <= i < m and 0 <= j < m and n > 0) for i, j, n in self.edges):
            raise ValueError('each edge needs valid ends and positive length')
        # Also check connectivity before making a boundary calculation.
        if self.laplacian.rank() != m-1:
            raise ValueError('the fibre graph must be connected')

    @cached_property
    def incidence(self):
        M = matrix(QQ, len(self.genera), len(self.edges))
        for k, (i, j, _) in enumerate(self.edges):
            M[i, k] += 1
            M[j, k] -= 1
        M.set_immutable()
        return M

    @cached_property
    def cycles(self):
        C = self.incidence.right_kernel().basis_matrix()
        C.set_immutable()
        return C

    @cached_property
    def laplacian(self):
        D = matrix.diagonal([1/n for _, _, n in self.edges])
        L = self.incidence*D*self.incidence.transpose()
        L.set_immutable()
        return L

    @cached_property
    def _inverse(self):
        keep = list(range(1, len(self.genera)))
        return self.laplacian.matrix_from_rows_and_columns(keep, keep).inverse()

    def _signature(self, potential):
        u = vector(QQ, [value for g, value in zip(self.genera, potential) if g])
        D = matrix.diagonal([n*(potential[i]+potential[j]) for i, j, n in self.edges])
        Q = self.cycles*D*self.cycles.transpose()
        return u, Q

    @cached_property
    def tate_signature(self):
        u, Q = self._signature([QQ(1)]*len(self.genera))
        out = vector(QQ, list(u)+Q.list())
        out.set_immutable()
        return out

    @cached_property
    def _quotient(self):
        tate = self.tate_signature
        return (matrix(QQ, [tate]).right_kernel().basis_matrix() if tate
                else identity_matrix(QQ, len(tate)))

    def boundary(self, multidegree):
        """Solve the vertical-divisor equation and return the boundary."""
        d = vector(QQ, multidegree)
        if len(d) != len(self.genera) or sum(d):
            raise ValueError('the graph divisor must have degree zero')
        u = vector(QQ, [0]+list(self._inverse*vector(QQ, d[1:])))
        g, Q = self._signature(u)
        return GraphBoundary(u, g, Q, self.tate_signature,
                             self._quotient*vector(QQ, list(g)+Q.list()))

    @cached_property
    def vertex_boundaries(self):
        n = len(self.genera)
        return tuple(self.boundary([int(i == j)-int(i == 0) for i in range(n)]).projected
                     for j in range(n))

    def stabilize(self):
        """Contract unmarked unstable rational vertices and return their survivors."""
        vertices = set(range(len(self.genera)))
        edges = list(self.edges)
        while True:
            for v in sorted(vertices):
                if self.genera[v]:
                    continue
                incident = [(k, b if a == v else a, n)
                            for k, (a, b, n) in enumerate(edges) if v in (a, b)]
                valence = sum(2 if w == v else 1 for _, w, _ in incident)
                if valence not in (1, 2) or any(w == v for _, w, _ in incident):
                    continue
                for k, _, _ in reversed(incident):
                    del edges[k]
                if valence == 2:
                    (_, a, s), (_, b, t) = incident
                    edges.append((a, b, s+t))
                vertices.remove(v)
                break
            else:
                break
        vertices = tuple(sorted(vertices))
        index = {v: i for i, v in enumerate(vertices)}
        return FibreGraph([self.genera[v] for v in vertices],
                          [(index[a], index[b], n) for a, b, n in edges]), vertices
