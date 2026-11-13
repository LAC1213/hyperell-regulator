r"""Direct Stokes integration on the two sheets of y^2=F(x).

One cut contains the branch values and all logarithmic zeros, followed by
an outward ray. Finite edges use x=a+(b-a) sin(pi*t/2)^2; the ray uses a
fixed reciprocal substitution. Root factors cancel endpoint singularities
analytically. There are no general covers, charts or Newton continuation.

The imaginary symmetric part of each boundary sum uses the closed formula
for second moments. Several log|x-a| weights share one set of edge periods.
"""
import numpy as np
from contextlib import nullcontext
from sage.all import QQ
from ._paths import semilunar_path
from ._periods import segment_differential
from .polynomials import polynomial
from .quadrature import Arithmetic, quadrature
from .diagnostics import _active, checkpoint


def _line(a, b, n):
    return list(a + (b - a) * np.linspace(0, 1, n)[1:])


def _arc(c, rad, th0, sweep, n):
    return list(c + rad * np.exp(1j * (th0 + sweep * np.linspace(0, 1, n)[1:])))


def _miter(v, eps, ea, eb):
    r"""
    Where the two offset lines meet at a corner concave for the side in hand.

    The offset is always `\varepsilon` to the left of the direction of travel,
    so on the return pass it is `i\varepsilon e_a` for the *reversed* direction.
    Putting an arc here instead would carry the path round the vertex and pick
    up its monodromy spuriously.
    """
    rhs = 1j * eps * (eb - ea)
    A = np.array([[ea.real, -eb.real], [ea.imag, -eb.imag]])
    s = np.linalg.solve(A, np.array([rhs.real, rhs.imag]))[0]
    return v + 1j * eps * ea + s * ea


def cut_scale(z):
    """Largest and smallest separations between cut nodes."""
    a = np.asarray(z)
    d = np.abs(a[:, None] - a[None, :])
    span = float(d.max()) or 1.0
    gap = float(d[d > 0].min()) if np.any(d > 0) else span
    return span, gap


def racetrack(z, d, eps=0.02, M=60.0, dens=30.0):
    """Boundary of a narrow neighborhood of the cut, with bank references."""
    span, gap = cut_scale(z)
    eps, M = eps * gap, M * span          # see cut_scale
    pts = list(z) + [z[-1] + M * d]
    r = len(z)
    e = [(pts[k + 1] - pts[k]) / abs(pts[k + 1] - pts[k]) for k in range(r)]
    ref = [(z[k] + z[k + 1]) / 2 for k in range(r - 1)] + [z[-1] + span * d]

    path, marks = [pts[0] + 1j * eps * e[0]], {}
    for s in (+1, -1):
        for k in (range(r) if s > 0 else range(r - 1, -1, -1)):
            end = (pts[k + 1] if s > 0 else pts[k]) + 1j * s * eps * e[k]
            mid = ref[k] + 1j * s * eps * e[k]
            path += _line(path[-1], mid,
                          max(8, int(dens * abs(mid - path[-1]) / span)))
            marks[(k, s)] = len(path) - 1
            path += _line(mid, end, max(8, int(dens * abs(end - mid) / span)))
            nxt = k + 1 if s > 0 else k - 1
            if 0 <= nxt < r:
                v = pts[k + 1] if s > 0 else pts[k]
                ea, eb = (e[k], e[nxt]) if s > 0 else (-e[k], -e[nxt])
                angle = np.angle(eb / ea)
                if abs(angle) < 1e-12:
                    pass              # collinear: the offsets already meet
                elif s * angle < 0:
                    path += _arc(v, eps, np.angle(path[-1] - v), angle, 60)
                else:
                    path += _line(path[-1], _miter(v, eps, ea, eb), 6)
        if s > 0:
            R = abs(path[-1] - z[-1])
            th0 = np.angle(path[-1] - z[-1])
            th1 = np.angle(pts[-1] - 1j * eps * e[-1] - z[-1])
            path += _arc(z[-1], R, th0, np.mod(th1 - th0, 2 * np.pi), 900)
    v = pts[0]
    th0, th1 = np.angle(path[-1] - v), np.angle(1j * e[0])
    sweep = np.mod(th1 - th0, 2 * np.pi)
    if np.mod(np.angle(e[0]) - th0, 2 * np.pi) < sweep:
        sweep -= 2 * np.pi
    path += _arc(v, eps, th0, sweep, 60)
    return np.array(path), marks, ref


class Cut:
    """A labelled semilunar cut on the x-sphere, independent of torsion orders."""
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
        start = points.index(C(self.centers[0])) if self.centers else max(range(len(points)), key=lambda i: abs(points[i]))
        nodes, indices, direction = semilunar_path(points, start)
        self.nodes = [C(z) if i is None else points[i] for z, i in zip(nodes, indices)]
        self.labels = [None if i is None else labels[i] for i in indices]
        self.direction = C(direction)
        self.span = R(cut_scale(nodes)[0])
        self.kinds = ['branch' if i is not None else 'regular' for i in self.labels]
        self.references = [(a+b)/2 for a, b in zip(self.nodes[:-1], self.nodes[1:])]+[self.nodes[-1]+self.span*self.direction]
        self.reference_y = []
        for i in range(len(self.nodes)-1):
            evaluate, _ = self.segment(i)
            a, b = self.nodes[i:i+2]
            self.reference_y.append((b-a)*R.pi()/2/evaluate(R('.5'))[0])
        z = self.references[-1]
        self.reference_y.append((C(self.F.leading_coefficient())*prod(z-r for r in self.roots)).sqrt())
        self.upper, self.lower, self.logs, self.jumps = self.topology()
        n = len(self.nodes)
        self.walks = [[(i, self.upper[i][s], 1) for i in range(n)]+
                      [(i, self.lower[i][s], -1) for i in reversed(range(n))] for s in (0, 1)]

    def segment(self, i, end=None):
        a = self.nodes[i]
        b = self.nodes[i+1] if end is None else end
        ib = self.labels[i+1] if end is None else None
        ar = self.ar
        return segment_differential(self.roots, ar.C(self.F.leading_coefficient()), a, b,
                                    self.labels[i], ib, self.g, sqrt=lambda z: z.sqrt(),
                                    sin=lambda z: z.sin(), cos=lambda z: z.cos(),
                                    pi=ar.R.pi(), sign_tol=ar.root_tol)

    def topology(self):
        """Read bank signs and logarithm jumps from a small contour of the cut.

        Each root factor has its argument unwrapped separately. This tracks
        only square-root signs; no algebraic fibres or Newton solves occur.
        """
        z = np.asarray([complex(v) for v in self.nodes])
        path, marks, _ = racetrack(z, complex(self.direction))
        roots = np.asarray([complex(r) for r in self.roots])
        angles = np.unwrap(np.angle(path[:, None]-roots), axis=0)
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
    a = cut.nodes[i]
    b = cut.nodes[i+1] if end is None else end
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
    z, B = cut.nodes[-1], cut.references[-1]
    D, yB = cut.span*cut.direction, cut.reference_y[-1]
    mid = (z+B)/2
    # Logs on this initial segment are continued back from the ray reference.
    logs = [lam[-1]+((mid-C(a))/(B-C(a))).log() for a, lam in zip(cut.centers, cut.logs)]
    first, yend = _finite(cut, len(cut.nodes)-1, B, logs)
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
    ar, g, n = cut.ar, cut.g, len(cut.nodes)
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
    for i in range(len(cut.nodes)):
        mon = _active.get()
        with mon.edge('finite' if i+1 < len(cut.nodes) else 'ray', degree=2, precision=prec) if mon else nullcontext():
            values = _finite(cut, i)[0] if i+1 < len(cut.nodes) else _ray(cut)
            for a in range(g):
                u[a].append([values[a], -values[a]])
                for k in range(len(cut.centers)):
                    v = values[(k+1)*g+a]
                    ell[k][a].append([v, -v])
    zero = [ar.C(0)]*len(cut.nodes)
    periods = _sum(cut, u, u, zero)
    regulators = [_sum(cut, u, e, jump) for e, jump in zip(ell, cut.jumps)]
    if details is not None:
        details.update(cut=cut, u=u, ell=ell)
    return periods, regulators
