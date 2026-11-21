"""Find new geometric torsion classes and compare their regulator moments.

Default examples are LMFDB 249.a.6723.1 and 256.a.512.1. Results retain
complete packet certificates and exact number-field coordinates. Numerical
regulators test new directions at the chosen complex embeddings; they do
not certify integrality over number fields or descent to Q.
"""

import argparse
import json
from pathlib import Path

import numpy as np
from sage.all import PolynomialRing, QQ, QQbar, matrix, number_field_elements_from_algebraics

from ..components import CACHE, curve_from_row, curve_points
from ..geometric_torsion import geometric_torsion_packet
from ..integrality import TorsionPacketIntegrality
from ..paths import data_path, ensure_parent
from ..regulator import regulators_stokes


def _ordered_roots(polynomial):
    """Fix the embedding convention by exact lexicographic real/imag order."""
    return sorted(polynomial.roots(QQbar, multiplicities=False),
                  key=lambda z: (z.real(), z.imag()))


def encode_point(point):
    """Exact coordinates in a number field, with an exact embedding index."""
    infinite = point[0] == 'infinity'
    values = [point[1]] if infinite else list(point)
    K, elements, embedding = number_field_elements_from_algebraics(values, embedded=True)
    if K is QQ:
        return dict(field_polynomial=None, root_index=None,
                    coordinates=['infinity', str(elements[0])] if infinite else list(map(str, elements)))
    roots = _ordered_roots(K.polynomial())
    root_index = roots.index(QQbar(embedding(K.gen())))
    coefficients = [[str(a) for a in element.list()] for element in elements]
    return dict(field_polynomial=list(map(str, K.polynomial())), root_index=root_index,
                coordinates=['infinity', coefficients[0]] if infinite else coefficients)


def decode_point(record):
    """Recover the exact QQbar point from an encoded record."""
    coordinates = record['coordinates']
    infinite = coordinates[0] == 'infinity'
    if record['field_polynomial'] is None:
        return (('infinity', QQbar(QQ(coordinates[1]))) if infinite else
                tuple(QQbar(QQ(a)) for a in coordinates))
    R = PolynomialRing(QQ, 'a')
    root = _ordered_roots(R(record['field_polynomial']))[record['root_index']]
    values = [R(a)(root) for a in (coordinates[1:] if infinite else coordinates)]
    return ('infinity', values[0]) if infinite else tuple(values)


def regulator_comparison(packet, rational_points, *, digits=25):
    """Compare new B(P-O) with previously certified rational pair classes.

    The base O must be a rational Weierstrass point. Put O at infinity,
    compute logarithmic moment vectors on one common cut, and test linear
    independence modulo the period vector using 3x3 determinants. These
    are numerical tests, not interval-certified nonvanishing proofs.
    """
    C, O = packet.curve, packet.basepoint
    R = PolynomialRing(QQ, 'u')
    u = R.gen()
    if O[0] == 'infinity' and packet.F.degree() == 5:
        G = R(packet.F)
        coordinate = lambda P: None if P[0] == 'infinity' else P[0]
    elif O[0] != 'infinity' and O[0] in QQ and packet.F(O[0]) == 0:
        root = QQ(O[0])
        shifted = R(packet.F)(u+root)
        G = R([shifted[6-i] for i in range(7)])
        coordinate = lambda P: (QQbar(0) if P[0] == 'infinity' else
                                None if P[0] == root else 1/(P[0]-root))
    else:
        raise NotImplementedError('regulator comparison needs a rational Weierstrass base')
    old = [packet._point(P) for P in rational_points]
    new = [P for P in packet.points() if P not in old]
    rational_test = TorsionPacketIntegrality(C, rational_points)
    pairs = []
    for i, P in enumerate(old):
        for j, Q in enumerate(old[:i]):
            D = [0]*len(old)
            D[i], D[j] = 1, -1
            if coordinate(P) != coordinate(Q) and rational_test.check(D)['q_integral'] is True:
                pairs.append((i, j))
    centers = list(dict.fromkeys(coordinate(P) for P in old+new if coordinate(P) is not None))
    S, columns = regulators_stokes(G, centers, prec=digits)
    lookup = dict(zip(centers, columns))
    zero = [v*0 for v in S]

    def column(P):
        center = coordinate(P)
        return zero if center is None else lookup[center]

    old_columns = [[a-b for a, b in zip(column(old[i]), column(old[j]))] for i, j in pairs]
    new_columns = [column(P) for P in new]
    def rank(columns):
        array = np.asarray([list(S)]+columns, dtype=float)
        return int(np.linalg.matrix_rank(array, tol=1e-10))-1
    determinants = [[str(matrix([col, new_col, S]).determinant()) for col in old_columns]
                    for new_col in new_columns]
    return dict(model=list(map(str, G)), digits=digits, certified=False,
                normalization='logarithmic moments (1, Re(u), |u|^2), compared modulo period vector',
                period_vector=list(map(str, S)),
                old_certified_pairs=[list(pair) for pair in pairs],
                old_vectors=[list(map(str, col)) for col in old_columns],
                new_vectors=[list(map(str, col)) for col in new_columns],
                determinants_old_new_period=determinants,
                old_regulator_rank=rank(old_columns),
                combined_regulator_rank=rank(old_columns+new_columns),
                new_points=[encode_point(P) for P in new])


def analyse(row, *, backend=None, timeout=300, digits=25, with_regulator=True):
    """Exact full packet and new algebraic classes for one cached LMFDB curve."""
    C = curve_from_row(row)
    old = curve_points(C, row)
    base = (-1, 0) if row['label'] == '249.a.6723.1' else None
    packet = geometric_torsion_packet(C, basepoint=base, backend=backend, timeout=timeout,
                                      progress=lambda s: print(s, flush=True))
    known = [packet._point(P) for P in old]
    points = packet.points()
    result = dict(label=row['label'], input_curve=row, size=packet.size,
                  basepoint=encode_point(packet.basepoint), complete=packet.complete,
                  backend=packet.metadata,
                  factors=[dict(polynomial=list(map(str, g)) if g is not None else None,
                                doubled_order=int(n)) for g, n in packet.factors],
                  points=[dict(point=encode_point(P), order=int(packet.order(P)),
                               previously_known=P in known) for P in points],
                  new_points=sum(P not in known for P in points),
                  integrality_of_new_classes='not established over number fields',
                  descent_to_Q='not attempted')
    if with_regulator:
        result['regulator_comparison'] = regulator_comparison(packet, old, digits=digits)
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--label', action='append')
    parser.add_argument('--cache', default=CACHE)
    parser.add_argument('--backend')
    parser.add_argument('--timeout', type=float, default=300)
    parser.add_argument('--digits', type=int, default=25)
    parser.add_argument('--no-regulator', action='store_true')
    parser.add_argument('--out', default=data_path('results', 'geometric_torsion.json'))
    args = parser.parse_args(argv)
    source = json.loads(Path(args.cache).read_text())
    rows = source['input_curves'] if isinstance(source, dict) else source
    labels = args.label or ['249.a.6723.1', '256.a.512.1']
    by_label = {row['label']: row for row in rows}
    results = []
    for label in labels:
        print('Computing '+label, flush=True)
        result = analyse(by_label[label], backend=args.backend, timeout=args.timeout,
                         digits=args.digits, with_regulator=not args.no_regulator)
        results.append(result)
        Path(ensure_parent(args.out)).write_text(json.dumps(
            dict(version=1, input_curves=[by_label[label] for label in labels], curves=results),
            indent=2)+'\n')
        print('%s: %s points, %s new' % (label, result['size'], result['new_points']), flush=True)


if __name__ == '__main__':
    main()
