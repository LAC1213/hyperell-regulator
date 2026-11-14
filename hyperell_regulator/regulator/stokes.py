r"""Direct Stokes integration on the two sheets of y^2=F(x).

A minimum spanning tree contains all branch values and logarithmic zeros,
with an outward ray at its rightmost vertex.
Finite edges use x=a+(b-a) sin(pi*t/2)^2; the ray uses a
fixed reciprocal substitution. Root factors cancel endpoint singularities
analytically. There are no general covers, charts or Newton continuation.

The imaginary symmetric part of each boundary sum uses the closed formula
for second moments. Several log|x-a| weights share one set of edge periods.
"""
import numpy as np
from contextlib import nullcontext
from sage.all import QQ
from ._tree import minimum_spanning_tree, orient_to_root, tree_contour
from ._periods import segment_differential
from .polynomials import polynomial
from .quadrature import Arithmetic, quadrature
from .diagnostics import _active, checkpoint


def cut_scale(z):
    """Largest and smallest separations between cut nodes."""
    a = np.asarray(z)
    d = np.abs(a[:, None] - a[None, :])
    span = float(d.max()) or 1.0
    gap = float(d[d > 0].min()) if np.any(d > 0) else span
    return span, gap


class Cut:
    """A labelled tree cut on the x-sphere, independent of torsion orders."""
    def __init__(self, F, centers=(), prec=15):
        self.ar = Arithmetic(prec)
        C, R = self.ar.C, self.ar.R
        self.F = polynomial(F)
        if self.F is None:
            raise ValueError('F must have exact rational coefficients')
        if self.F.degree() < 5 or self.F.gcd(self.F.derivative()).degree() > 0:
            raise ValueError('F must be squarefree and have degree at least five')
        self.g = int((self.F.degree()-1)//2)
        self.roots = list(self.F.roots(C, multiplicities=False))
        # Preserve exact rational coincidences with logarithmic zeros.
        for r in self.F.roots(QQ, multiplicities=False):
            i = min(range(len(self.roots)), key=lambda i: abs(self.roots[i]-C(r)))
            self.roots[i] = C(r)
        self.centers = list(dict.fromkeys(QQ(a) for a in centers))
        points = list(self.roots)
        labels = list(range(len(points)))
        for a in self.centers:
            if C(a) not in points:
                points.append(C(a))
                labels.append(None)
        self.nodes, self.labels = points, labels
        root = max(range(len(points)), key=lambda i: (float(points[i].real()), float(points[i].imag())))
        self.edges = orient_to_root(len(points), minimum_spanning_tree(points), root)
        self.direction = C(1)
        self.ray_vertex = self.edges[-1][0]
        self.span = R(cut_scale(self.nodes)[0])
        self.kinds = ['branch' if i is not None else 'regular' for i in self.labels]
        self.references = [(self.nodes[a]+self.nodes[b])/2 if b is not None else
                           self.nodes[a]+self.span*self.direction for a, b in self.edges]
        self.reference_y = []
        for i in range(len(self.edges)-1):
            evaluate, _ = self.segment(i)
            a, b, _, _ = self.edge_points(i)
            self.reference_y.append((b-a)*R.pi()/2/evaluate(R('.5'))[0])
        z = self.references[-1]
        self.reference_y.append((C(self.F.leading_coefficient())*prod(z-r for r in self.roots)).sqrt())
        self.upper, self.lower, self.logs, self.jumps = self.topology()
        self.walks = [[(i, (self.upper if side > 0 else self.lower)[i][sheet], side)
                       for i, side in self.boundary_word] for sheet in (0, 1)]

    def edge_points(self, i, end=None):
        a, b = self.edges[i]
        return (self.nodes[a], self.nodes[b] if end is None else end,
                self.labels[a], self.labels[b] if end is None else None)

    def segment(self, i, end=None):
        a, b, ia, ib = self.edge_points(i, end)
        ar = self.ar
        return segment_differential(self.roots, ar.C(self.F.leading_coefficient()), a, b,
                                    ia, ib, self.g, sqrt=lambda z: z.sqrt(),
                                    sin=lambda z: z.sin(), cos=lambda z: z.cos(),
                                    pi=ar.R.pi(), sign_tol=ar.root_tol)

    def topology(self):
        """Read bank signs and logarithm jumps from a small contour of the cut.

        Each root factor has its argument unwrapped separately. This tracks
        only square-root signs; no algebraic fibres or Newton solves occur.
        """
        z = np.asarray([complex(v) for v in self.nodes])
        path, marks, self.boundary_word = tree_contour(
            z, self.edges, complex(self.direction), self.references,
            .02*cut_scale(z)[1], float(self.span))
        self.seed = complex(path[0])
        roots = np.asarray([complex(r) for r in self.roots])
        angles = np.unwrap(np.angle(path[:, None]-roots), axis=0)
        if np.max(abs(angles[-1]-angles[0])) > 1e-8:
            raise RuntimeError('cut contour has nonzero winding around a branch value')
        phases = np.exp(.5j*angles.sum(axis=1))*np.sqrt(complex(self.F.leading_coefficient()))
        upper, lower = [], []
        for i, ref in enumerate(self.references):
            for side, output in ((1, upper), (-1, lower)):
                k = marks[i, side]
                # Continue the bank sample straight to the reference point.
                factor = np.prod(np.sqrt((complex(ref)-roots)/(path[k]-roots)))
                y = phases[k]*np.prod(np.sqrt(abs(path[k]-roots)))*factor
                j = 0 if abs(y-complex(self.reference_y[i])) <= abs(y+complex(self.reference_y[i])) else 1
                output.append([j, 1-j])
        logs, jumps = [], []
        for a in self.centers:
            theta = np.unwrap(np.angle(path-float(a)))
            values, differences = [], []
            for i, ref in enumerate(self.references):
                angle = np.angle(complex(ref)-float(a))
                turns = {s: round((theta[marks[i, s]]-angle)/(2*np.pi)) for s in (1, -1)}
                values.append((ref-self.ar.C(a)).log()+self.ar.C.gen()*2*self.ar.R.pi()*turns[-1])
                differences.append(self.ar.C.gen()*2*self.ar.R.pi()*(turns[1]-turns[-1]))
            logs.append(values)
            jumps.append(differences)
        return upper, lower, logs, jumps


def prod(values):
    result = 1
    for v in values:
        result *= v
    return result


def _finite(cut, i, end=None, reference_logs=None):
    ar, g = cut.ar, cut.g
    a, b, _, _ = cut.edge_points(i, end)
    mid = (a+b)/2
    evaluate, yend = cut.segment(i, end)
    logs = [row[i] for row in cut.logs] if reference_logs is None else reference_logs
    def sample(t):
        values = evaluate(t)
        s, c = (ar.R.pi()*t/2).sin(), (ar.R.pi()*t/2).cos()
        weighted = []
        for center, lam in zip(cut.centers, logs):
            center = ar.C(center)
            if center == a:
                ll = lam+ar.R(2).log()+2*s.log()
            elif center == b:
                ll = lam+ar.R(2).log()+2*c.log()
            else:
                ll = lam+(((a-center)*c*c+(b-center)*s*s)/(mid-center)).log()
            weighted.extend(ll*v for v in values)
        return values+weighted
    return quadrature(ar, sample, g*(1+len(cut.centers))), yend


def _ray(cut):
    ar, g = cut.ar, cut.g
    C, R = ar.C, ar.R
    z, B = cut.nodes[cut.ray_vertex], cut.references[-1]
    D, yB = cut.span*cut.direction, cut.reference_y[-1]
    mid = (z+B)/2
    # Logs on this initial segment are continued back from the ray reference.
    logs = [lam[-1]+((mid-C(a))/(B-C(a))).log() for a, lam in zip(cut.centers, cut.logs)]
    first, yend = _finite(cut, len(cut.edges)-1, B, logs)
    sign = 1 if abs(yend-yB) <= abs(yend+yB) else -1
    first = [sign*v for v in first]
    degree = cut.F.degree()
    def sample(t):
        v = 1-t
        h = v*v
        X = B*h+D*(1-h)
        y = yB*prod(((X-r*h)/(B-r)).sqrt() for r in cut.roots)
        values = [2*D*v**(degree-3-2*j)*X**j/y for j in range(g)]
        weighted = []
        for a, lam in zip(cut.centers, cut.logs):
            ll = lam[-1]+((X-C(a)*h)/(B-C(a))).log()-2*v.log()
            weighted.extend(ll*r for r in values)
        return values+weighted
    tail = quadrature(ar, sample, len(first))
    return [a+b for a, b in zip(first, tail)]


def _sum(cut, u, ell, jumps):
    """Boundary primitives plus the closed, symmetrized second moments."""
    ar, g, n = cut.ar, cut.g, len(cut.edges)
    H = [[ar.C(0) for _ in range(g)] for _ in range(g)]
    forward, backward = {}, {}
    for k, walk in enumerate(cut.walks):
        for j, (i, sheet, side) in enumerate(walk):
            (forward if side > 0 else backward)[i, sheet] = k, j
    for a in range(g):
        U = []
        for walk in cut.walks:
            row = [ar.C(0)]
            for i, s, side in walk:
                row.append(row[-1]+(ell[a][i][s]+jumps[i]*u[a][i][s] if side > 0 else -ell[a][i][s]))
            U.append(row)
        scale = max([ar.R(1)]+[abs(v) for row in ell[a] for v in row])
        defect = max(abs(row[-1]) for row in U)/scale
        checkpoint(boundary_defect=str(defect))
        if defect > ar.R(10)**(-ar.digits+3):
            raise RuntimeError('sheet boundary does not close: %s' % defect)
        for i in range(n):
            for s in (0, 1):
                k, j = forward[i, s]
                l, jj = backward[i, s]
                constant = U[k][j]-U[l][jj]+ell[a][i][s]
                for b in range(g):
                    H[a][b] += constant*u[b][i][s].conjugate()
    result = []
    for a in range(g):
        for b in range(a, g):
            moment = sum(jumps[i].imag()/2*(u[a][i][s]*u[b][i][s].conjugate()).real()
                         for i in range(n) for s in (0, 1))
            result.append(ar.output(-((H[a][b].imag()+H[b][a].imag())/2+moment)/4))
    return result


def compute(F, centers=(0,), prec=15, *, details=None):
    """Return periods and log|x-a| vectors, sharing all integrations."""
    cut = Cut(F, centers, prec)
    ar, g = cut.ar, cut.g
    u = [[] for _ in range(g)]
    ell = [[[] for _ in range(g)] for _ in cut.centers]
    for i in range(len(cut.edges)):
        mon = _active.get()
        with mon.edge('finite' if i+1 < len(cut.edges) else 'ray', degree=2, precision=prec) if mon else nullcontext():
            values = _finite(cut, i)[0] if i+1 < len(cut.edges) else _ray(cut)
            for a in range(g):
                u[a].append([values[a], -values[a]])
                for k in range(len(cut.centers)):
                    v = values[(k+1)*g+a]
                    ell[k][a].append([v, -v])
    zero = [ar.C(0)]*len(cut.edges)
    periods = _sum(cut, u, u, zero)
    regulators = [_sum(cut, u, e, jump) for e, jump in zip(ell, cut.jumps)]
    if details is not None:
        details.update(cut=cut, u=u, ell=ell)
    return periods, regulators
