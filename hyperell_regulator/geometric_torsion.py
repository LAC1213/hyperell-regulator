r"""Complete genus-two hyperelliptic torsion packets via Poonen's algorithm.

Uses the author's PARI/GP implementation of ``Computing torsion points on
curves`` (2001), with the 2024 verification update, as an external backend.
Its successful result includes exact verification and a p-adic completeness
argument, rather than a bounded search for torsion orders. Download the
backend separately with :func:`install_poonen_backend`.

The supported packet is the packet containing the Weierstrass points. Any
base point in this packet, including a non-Weierstrass algebraic point, can
be selected. A base point outside it is explicitly unsupported: returning
this packet for such a base point would be incorrect.
"""

import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import time

from sage.all import (
    HyperellipticCurve, PolynomialRing, QQ, QQbar, ZZ, lcm,
    number_field_elements_from_algebraics,
)

from . import _subprocess
from .paths import data_path, ensure_parent

SOURCE_URL = 'https://math.mit.edu/~poonen/papers/torsion.gp'
PAPER_URL = 'https://math.mit.edu/~poonen/papers/curvetorsion.pdf'
SOURCE_SHA256 = '668d836ed87de83a23ea260aa3a7b14c994d5851eef4ee28fd99333341f86df6'
RESULT_MARKER = 'HYPERELL_TORSION_RESULT='


def install_poonen_backend(path=None, *, source_path=None):
    """Install the pinned upstream GP file in the local results cache.

    With ``source_path``, copy an already downloaded file without network
    access. The SHA256 is checked before any destination file is written.
    The upstream source is not redistributed as part of this package.
    """
    if source_path is None:
        import urllib.request
        with urllib.request.urlopen(SOURCE_URL, timeout=60) as response:
            content = response.read()
    else:
        content = Path(source_path).read_bytes()
    if hashlib.sha256(content).hexdigest() != SOURCE_SHA256:
        raise ValueError('Poonen backend checksum differs from the pinned upstream version')
    destination = Path(path or data_path('results', 'poonen', 'torsion.gp'))
    ensure_parent(destination)
    destination.write_bytes(content)
    return destination


def _polynomials(C):
    R = PolynomialRing(QQ, 'x')
    f, h = (R(C.f), R(C.h)) if hasattr(C, 'f') else map(R, C.hyperelliptic_polynomials())
    F = 4*f + h*h
    if F.degree() not in (5, 6) or not F.is_squarefree():
        raise ValueError('a smooth genus-two curve over QQ is required')
    return f, h, F


def _gp_script(F, backend):
    """Input contains only validated polynomial coefficients and a quoted path."""
    # Integer coefficients avoid parsing expressions supplied by callers.
    polynomial = '+'.join('(%s)*x^%s' % (ZZ(a), i) for i, a in enumerate(F))
    return '\n'.join([
        'read(%s);' % json.dumps(str(Path(backend).resolve())),
        'hrcurve=%s;' % polynomial,
        '{iferr(',
        'hrresult=torsion(hrcurve);',
        'if(type(hrresult)!="t_VEC",error("torsion computation did not complete"));',
        'setupsextic(hrcurve);hrrows=List();',
        'for(hri=1,#hrresult[2],hrpoly=hrresult[2][hri];',
        'if(hrpoly==infty,',
        'if(poldegree(hrcurve)==5,hrorder=1,',
        'setupsextic(x^6*subst(hrcurve,x,1/x));hrorder=findorder(xpoint(0)));',
        'listput(hrrows,[[],hrorder]);setupsextic(hrcurve),',
        'hrpoly=hrpoly/content(hrpoly);if(pollead(hrpoly)<0,hrpoly=-hrpoly);',
        'hrorder=if(hrcurve%hrpoly==0,1,numberfieldtesttorsion(hrpoly));',
        'if(type(hrorder)!="t_INT",error("exact torsion verification failed"));',
        'listput(hrrows,[Vecrev(hrpoly),hrorder])));',
        'print("%s",[hrresult[1],Vec(hrrows)]),' % RESULT_MARKER,
        'hrerror,print("HYPERELL_TORSION_ERROR=",hrerror);quit(1))}',
        'quit;', '',
    ])


def _validate_result(F, data):
    """Check irreducibility, counts, Weierstrass coverage and order witnesses."""
    R = F.parent()
    if not isinstance(data, list) or len(data) != 2:
        raise ValueError('malformed torsion result')
    count, rows = data
    if type(count) is not int or not isinstance(rows, list):
        raise ValueError('malformed torsion count or rows')
    found, branch_product, total, infinite = set(), R(1), 0, False
    for entry in rows:
        if not isinstance(entry, list) or len(entry) != 2:
            raise ValueError('malformed torsion factor')
        coefficients, order = entry
        if type(order) is not int or order < 1 or order > 1000:
            raise ValueError('invalid exact doubled-divisor order')
        if not isinstance(coefficients, list) or any(type(a) is not int for a in coefficients):
            raise ValueError('invalid torsion polynomial coefficients')
        if not coefficients:
            if infinite:
                raise ValueError('duplicate point at infinity')
            infinite = True
            if F.degree() == 5 and order != 1:
                raise ValueError('the Weierstrass point at infinity has a nontrivial doubled order')
            total += 1 if F.degree() == 5 else 2
            continue
        g = R(coefficients).monic()
        if g.degree() < 1 or not g.is_irreducible() or g in found:
            raise ValueError('duplicate or reducible torsion factor')
        found.add(g)
        if F % g == 0:
            if order != 1:
                raise ValueError('Weierstrass factor has a nontrivial doubled order')
            branch_product *= g
            total += g.degree()
        else:
            total += 2*g.degree()
    if branch_product.monic() != F.monic() or (F.degree() == 5 and not infinite):
        raise ValueError('torsion result omits Weierstrass points')
    if total != count:
        raise ValueError('torsion packet size disagrees with its factors')
    return data


class GeometricTorsionPacket:
    """An exactly verified complete torsion packet, with QQbar coordinates.

    ``points()`` returns coordinates on the original ``y^2+h*y=f`` model.
    The polynomial factors compactly represent all geometric points, not
    just rational x-coordinates. ``order(P)`` computes the exact order of
    ``P-basepoint`` by Mumford arithmetic over an explicit number field.
    """

    def __init__(self, C, data, *, basepoint=None, metadata=None):
        self.curve = C
        self.f, self.h, self.F = _polynomials(C)
        _validate_result(self.F, data)
        self.size = int(data[0])
        self.factors = tuple((self.F.parent()(a).monic() if a else None, ZZ(n))
                             for a, n in data[1])
        self.metadata = dict(metadata or {})
        self.complete = self.metadata.get('complete') is True
        self._points = None
        if basepoint is None:
            basepoint = (('infinity', -QQbar(self.h[3])/2) if self.F.degree() == 5
                         else (self.F.roots(QQbar, multiplicities=False)[0], None))
            if basepoint[1] is None:
                basepoint = (basepoint[0], -self.h(basepoint[0])/2)
        self.basepoint = self._point(basepoint)
        if not self.contains(self.basepoint):
            raise NotImplementedError('base point is outside the hyperelliptic torsion packet; '
                                      'general non-Weierstrass packets are not implemented')

    def _point(self, P):
        if isinstance(P, str) and P == 'infinity':
            if self.F.degree() != 5:
                raise ValueError('specify which point at infinity')
            P = ('infinity', -self.h[3]/2)
        if len(P) != 2:
            raise ValueError('use (x,y) or ("infinity", eta) coordinates')
        if P[0] == 'infinity':
            eta = QQbar(P[1])
            if eta*eta + self.h[3]*eta != self.f[6]:
                raise ValueError('invalid point at infinity')
            return ('infinity', eta)
        a, b = map(QQbar, P)
        if b*b + self.h(a)*b != self.f(a):
            raise ValueError('point is not on the original curve')
        return a, b

    def contains(self, P):
        """Membership in this packet, after validating the curve coordinates."""
        P = self._point(P)
        if P[0] == 'infinity':
            return any(g is None for g, _ in self.factors)
        return any(g is not None and g(P[0]) == 0 for g, _ in self.factors)

    def points(self, *, rational_abscissae=False):
        """Exact QQbar points, optionally restricted to rational abscissae and infinity."""
        if rational_abscissae:
            return self._points_from_factors((g, n) for g, n in self.factors
                                             if g is None or g.degree() == 1)
        if self._points is None:
            points = self._points_from_factors(self.factors)
            if len(points) != self.size:
                raise ArithmeticError('QQbar enumeration disagrees with the certified packet size')
            self._points = tuple(points)
        return list(self._points)

    def _points_from_factors(self, factors):
        points = []
        for g, _ in factors:
            abscissae = ('infinity',) if g is None else g.roots(QQbar, multiplicities=False)
            for a in abscissae:
                h, F = (self.h[3], self.F[6]) if a == 'infinity' else (self.h(a), self.F(a))
                r = QQbar(F).sqrt()
                points.append((a, (-h+r)/2))
                if r:
                    points.append((a, (-h-r)/2))
        return points

    def _doubled_order(self, P):
        for g, n in self.factors:
            if (g is None and P[0] == 'infinity') or (g is not None and P[0] != 'infinity' and g(P[0]) == 0):
                return n
        raise ValueError('point is outside this torsion packet')

    def order(self, P, *, max_field_degree=64):
        """Exact order of P-basepoint, using a proved multiple from the backend.

        The backend verifies the order of ``2P-canonical``. Thus twice the
        lcm of the two doubled orders annihilates P-basepoint. This is a
        rigorous bound; no finite-prime bound on rational torsion is used
        for algebraic points. Large composita raise an explicit error.
        """
        P, O = self._point(P), self.basepoint
        if not self.contains(P):
            raise ValueError('point is outside this torsion packet')
        if P == O:
            return ZZ(1)
        w = None if self.F.degree() == 5 else self.F.roots(QQbar, multiplicities=False)[0]
        values = [a for point in (P, O) for a in point if a != 'infinity']
        if w is not None:
            values.append(w)
        K, elements, embedding = number_field_elements_from_algebraics(values, embedded=True)
        if K.degree() > max_field_degree:
            raise NotImplementedError('exact order needs a number field of degree %s' % K.degree())
        field_values = iter(elements)
        mapped = []
        for point in (P, O):
            mapped.append(('infinity', next(field_values)) if point[0] == 'infinity'
                          else (next(field_values), next(field_values)))
        R = PolynomialRing(K, 'u')
        u = R.gen()
        if w is None:
            G = R(self.F)
        else:
            root = next(field_values)
            shifted = R(self.F)(u+root)
            G = R([shifted[6-i] for i in range(7)])
        D = HyperellipticCurve(G)
        J = D.jacobian()(K)

        def cls(point):
            if point[0] == 'infinity':
                if w is None:
                    return J(0)
                x0, y0 = K(0), 2*point[1]+K(self.h[3])
            else:
                a, b = point
                Y = 2*b + R(self.h)(a)
                if w is None:
                    x0, y0 = a, Y
                elif a == root:
                    return J(0)
                else:
                    x0, y0 = 1/(a-root), Y/(a-root)**3
            return J(D.point([x0, y0, K(1)]))

        difference = cls(mapped[0])-cls(mapped[1])
        order = 2*lcm(self._doubled_order(P), self._doubled_order(O))
        if order*difference != J(0):
            raise ArithmeticError('the backend torsion multiple failed exact Jacobian verification')
        for p in order.prime_divisors():
            while order % p == 0 and (order//p)*difference == J(0):
                order //= p
        return order


def geometric_torsion_packet(C, basepoint=None, *, backend=None, gp=None,
                             timeout=300, cache=True, progress=None):
    """Compute a complete geometric packet; timeout/failure never means empty.

    ``backend`` is the pinned upstream ``torsion.gp`` file; it defaults to
    ``HYPERELL_POONEN_GP`` or the installed local cache. ``gp`` can select
    a PARI/GP executable. Successful, validated results are cached by curve
    and backend checksum. A timeout raises ``subprocess.TimeoutExpired``.
    """
    _, _, F = _polynomials(C)
    denominator = lcm([a.denominator() for a in F])
    integral = (denominator**2*F).change_ring(ZZ)
    backend = Path(backend or os.environ.get('HYPERELL_POONEN_GP') or
                   data_path('results', 'poonen', 'torsion.gp'))
    if not backend.is_file():
        raise FileNotFoundError('install the external GP backend with install_poonen_backend() '
                                'or supply backend=/path/to/torsion.gp')
    digest = hashlib.sha256(backend.read_bytes()).hexdigest()
    if digest != SOURCE_SHA256:
        raise ValueError('GP backend is not the pinned, inspected upstream version')
    key = hashlib.sha256(json.dumps([list(map(str, integral)), digest]).encode()).hexdigest()
    cache_file = Path(data_path('results', 'poonen', key+'.json'))
    if cache and cache_file.is_file():
        saved = json.loads(cache_file.read_text())
        if (saved.get('source_sha256') == digest and saved.get('curve') == list(map(str, integral))
                and saved.get('complete') is True):
            return GeometricTorsionPacket(C, saved['data'], basepoint=basepoint,
                                          metadata=dict(saved, cached=True))
    executable = gp or shutil.which('gp')
    if not executable:
        raise FileNotFoundError('PARI/GP executable not found; supply gp=/path/to/gp')
    started = time.monotonic()
    completed = _subprocess.run([str(executable), '-q', '-f', '-s', '256M'],
                               input=_gp_script(integral, backend), timeout=timeout,
                               progress=progress, label='Poonen geometric torsion packet')
    lines = [line[len(RESULT_MARKER):] for line in completed.stdout.splitlines()
             if line.startswith(RESULT_MARKER)]
    if completed.returncode or completed.stderr.strip() or 'SUCCESS!' not in completed.stdout or len(lines) != 1:
        raise ArithmeticError('Poonen backend did not certify a complete packet: '+
                              (completed.stderr or completed.stdout)[-1800:])
    data = _validate_result(F, json.loads(lines[0]))
    saved = dict(data=data, curve=list(map(str, integral)), source=SOURCE_URL,
                 paper=PAPER_URL, source_sha256=digest, seconds=time.monotonic()-started,
                 complete=True, transcript_sha256=hashlib.sha256(completed.stdout.encode()).hexdigest())
    packet = GeometricTorsionPacket(C, data, basepoint=basepoint, metadata=dict(saved, cached=False))
    if cache:
        ensure_parent(cache_file)
        with tempfile.NamedTemporaryFile(mode='w', dir=cache_file.parent, delete=False) as fh:
            temporary = fh.name
            json.dump(saved, fh, indent=2)
            fh.write('\n')
        os.replace(temporary, cache_file)
        cache_file.with_suffix('.log').write_text(completed.stdout)
    return packet
