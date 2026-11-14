"""Guarded tanh-sinh quadrature for explicit, single-valued differentials."""
import math
from sage.all import ComplexField, RealField
from .diagnostics import checkpoint, phase


class QuadratureFailure(RuntimeError):
    """Successive quadrature levels did not meet the requested accuracy."""


class Arithmetic:
    def __init__(self, digits):
        if int(digits) != digits or digits < 5:
            raise ValueError('prec must be an integer of at least 5 decimal digits')
        self.digits = int(digits)
        self.bits = int(math.ceil((digits+30)*math.log2(10)))
        self.C, self.R = ComplexField(self.bits), RealField(self.bits)
        self.output = RealField(int(math.ceil(digits*math.log2(10))))
        self.tol = self.R(10)**(-digits-3)
        self.root_tol = self.R(10)**(-digits-20)
        self._rules = {}

    def rule(self, level):
        if level not in self._rules:
            R = self.R
            h, lam = R(2)**(-level), R.pi()/2
            reach = (self.digits+10)*R(10).log()/2
            count = int((reach/lam).arcsinh()/h)+1
            nodes = []
            for k in range(-count, count+1):
                t = h*k
                sh = lam*t.sinh()
                q = (2*abs(sh)).exp()
                v = 1/(1+q)
                node = v if sh < 0 else 1-v
                if 0 < node < 1:
                    nodes.append((node, h*lam*t.cosh()*2*q/(1+q)**2))
            self._rules[level] = nodes
        return self._rules[level]


def quadrature(ar, evaluate, dimension, max_level=15):
    """Integrate on (0,1), caching samples and checking successive levels."""
    cache, previous, error = {}, None, None
    for level in range(3, max_level+1):
        with phase('quadrature', level=level, precision=ar.digits):
            total = [ar.C(0) for _ in range(dimension)]
            for t, weight in ar.rule(level):
                if t not in cache:
                    checkpoint('quadrature_nodes')
                    cache[t] = evaluate(t)
                for i, value in enumerate(cache[t]):
                    total[i] += weight*value
            if previous is not None:
                error = max(abs(a-b)/max(1, abs(a)) for a, b in zip(total, previous))
                checkpoint(scaled_quadrature_error=str(error), level=level)
                if error <= ar.tol:
                    return total
            previous = total
    raise QuadratureFailure('quadrature did not converge: error=%s, tolerance=%s' % (error, ar.tol))
