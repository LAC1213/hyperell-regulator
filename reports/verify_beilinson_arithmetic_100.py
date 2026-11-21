"""Verify the saved local computations and mod-3 surjectivity witnesses.

Run with ``sage -python reports/verify_beilinson_arithmetic_100.py``.
The negative torsion/endomorphism arguments use the saved algebraic LMFDB
data. Five negative image classifications are database results, not proved
by this verifier. Positive answers use Tests 4.4 and 4.5 of
https://arxiv.org/html/2301.02222v2#S4, with ell=3.
"""
import hashlib
import csv
import json
from pathlib import Path

from sage.all import GF, QQ, ZZ, HyperellipticCurve, PolynomialRing


ROOT = Path(__file__).resolve().parents[1]
R = PolynomialRing(QQ, 'x')
S = PolynomialRing(GF(3), 'X')
# Reductions of the two sets in Lemma 2.4 modulo 3.
C1920 = {(a, b) for a in range(3) for b in range(3)} - {(2, 0)}
C720 = {(0, 0), (0, 1), (0, 2), (1, 0), (1, 1)}


def counts(f, h, p):
    values = []
    for degree in (1, 2):
        field = GF(p**degree, 'a')
        T = PolynomialRing(field, 'x')
        g, j = T(f), T(h)
        affine = sum(1 for x in field for y in field if y*y+j(x)*y == g(x))
        infinity = sum(1 for y in field if y*y+field(h[3])*y == field(f[6]))
        values.append(affine+infinity)
    return values


def euler_coefficients(n1, n2, p):
    a1 = n1-p-1
    numerator = n2-p*p-1+a1*a1
    assert numerator % 2 == 0
    a2 = numerator//2
    return [1, a1, a2, p*a1, p*p]


def verify():
    data = json.loads((ROOT/'reports/beilinson_arithmetic_100.json').read_text())
    report = ROOT/data['report_path']
    assert hashlib.sha256(report.read_bytes()).hexdigest() == data['report_sha256']
    assert [r['label'] for r in data['records']] == [r['label'] for r in json.loads(report.read_text())]
    metadata = {r['label']: r for r in data['lmfdb_records']}
    certified_full = negative_arguments = database_only = 0
    for row in data['records']:
        raw = metadata[row['label']]
        eqn = json.loads(raw['eqn']) if isinstance(raw['eqn'], str) else raw['eqn']
        f, h = map(R, eqn)
        F = 4*f+h*h
        D = ZZ(F.discriminant()*(F.leading_coefficient()**2 if F.degree() == 5 else 1)/4096)
        assert D == row['discriminant'] and abs(D) == raw['abs_disc']
        for p in (2, 3):
            local = row['local'][str(p)]
            assert local['good'] == (D % p != 0)
            if not local['good']:
                assert local['status'] == 'bad' and not local['good_ordinary']
                continue
            n1, n2 = counts(f, h, p)
            # Compare independent direct enumeration with Sage point counting.
            T = PolynomialRing(GF(p), 'x')
            assert HyperellipticCurve(T(f), T(h)).count_points(2) == [n1, n2]
            assert local['counts'] == [n1, n2]
            coefficients = euler_coefficients(n1, n2, p)
            assert coefficients == local['euler_coefficients']
            assert local['ordinary'] == local['good_ordinary'] == bool(coefficients[2] % p)
            assert local['status'] == ('ordinary' if local['ordinary'] else 'nonordinary')
        image = row['galois_mod3']
        assert row['all_three'] == (image['surjective'] and all(
            row['local'][p]['good_ordinary'] for p in ('2', '3')))
        labels = [g for g in raw['modell_images'] if g.startswith('3.')]
        assert image['surjective'] == (not labels)
        if image['surjective']:
            witnesses = image['certificate']['witnesses']
            checked = set()
            for condition, witness in witnesses.items():
                p = witness['prime']
                assert p != 3 and D % p != 0
                T = PolynomialRing(GF(p), 'x')
                n1, n2 = map(int, HyperellipticCurve(T(f), T(h)).count_points(2))
                coefficients = list(reversed(euler_coefficients(n1, n2, p)))
                assert coefficients == witness['frobenius_coefficients']
                a, b = -coefficients[3], coefficients[2]
                pair = (a*a*pow(p, -1, 3) % 3, b*pow(p, -1, 3) % 3)
                polynomial = S(coefficients)
                tests = dict(exceptional_1920=pair not in C1920,
                             exceptional_720=pair not in C720,
                             irreducible=polynomial.is_irreducible(),
                             simple_linear_nonzero_trace=(a % 3 != 0 and any(
                                 g.degree() == 1 and k == 1 for g, k in polynomial.factor())))
                assert tests[condition]
                checked.add(condition)
            assert checked == set(tests)
            # The remaining exceptional test is automatic because ell != 7.
            # The similitude character is the surjective cyclotomic character mod 3.
            certified_full += 1
        elif image['basis'] == 'Nonzero rational 3-torsion':
            assert raw['torsion_order'] % 3 == 0
            negative_arguments += 1
        elif image['basis'] == 'Nontrivial rational endomorphisms':
            assert raw['end_alg'] not in ('Q', None)
            negative_arguments += 1
        else:
            assert labels and not image['certified_without_probabilistic_image_classification']
            database_only += 1
    with (ROOT/'reports/beilinson_arithmetic_100.csv').open() as fh:
        table = list(csv.DictReader(fh))
    assert [r['label'] for r in table] == [r['label'] for r in data['records']]
    for entry, row in zip(table, data['records']):
        for p in ('2', '3'):
            assert entry['reduction_at_'+p] == row['local'][p]['status']
            assert entry['good_ordinary_at_'+p] == str(row['local'][p]['good_ordinary'])
        assert entry['mod3_surjective'] == str(row['galois_mod3']['surjective'])
        assert entry['all_three'] == str(row['all_three'])
    print('Verified %d local records, %d surjectivity certificates, %d negative arguments; '
          '%d further negative answers use LMFDB classifications.' %
          (len(data['records']), certified_full, negative_arguments, database_only))


if __name__ == '__main__':
    verify()
