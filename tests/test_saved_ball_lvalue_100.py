"""Exact-arithmetic checks of the saved conditional 100-digit L-value enclosure.

These check the artifact, its normalization and its inputs. Reproducing the
analytic L-value bound additionally requires the documented patched backend.
"""

from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).parents[1]


@pytest.fixture(scope='module')
def certificate():
    return json.loads((ROOT / 'reports/ball_lvalue_249_100.json').read_text())


def interval(row):
    scale = Fraction(2) ** row['exponent']
    return int(row['lower']) * scale, int(row['upper']) * scale


def rounded(value, places):
    scale = 10**places
    n = (value * scale + Fraction(1, 2)).__floor__()
    return str(n // scale) + '.' + str(n % scale).zfill(places)


def test_certified_width_and_decimal_rounding(certificate):
    assert certificate['status'] == 'conditional_enclosure'
    assert certificate['run']['exit_code'] == 0
    assert certificate['run']['achieved_relative_bits'] >= 333
    lo, hi = interval(certificate['intervals']['L_second_derivative'])
    assert 0 < lo < hi
    assert (hi-lo)/2 < Fraction(1, 10**100)
    assert (hi-lo)/(hi+lo) < Fraction(1, 10**100)
    assert rounded(lo, 100) == rounded(hi, 100) == certificate['rounded_decimal_100_places']
    ilo, ihi = interval(certificate['intervals']['L_at_2_imag'])
    assert ilo <= 0 <= ihi


def arctan_reciprocal_bounds(q, terms=200):
    total = sum((Fraction((-1)**k, (2*k+1) * q**(2*k+1)) for k in range(terms)), Fraction(0))
    next_term = Fraction((-1)**terms, (2*terms+1) * q**(2*terms+1))
    return min(total, total+next_term), max(total, total+next_term)


def test_normalization_with_independent_rational_pi_bounds(certificate):
    # Machin's identity and alternating-series bounds give a rational pi
    # enclosure much narrower than the backend's working precision.
    a, b = arctan_reciprocal_bounds(5), arctan_reciprocal_bounds(239)
    pi_lo, pi_hi = 16*a[0] - 4*b[1], 16*a[1] - 4*b[0]
    expected_lo, expected_hi = 249/(4*pi_hi**2), 249/(4*pi_lo**2)
    scale_lo, scale_hi = interval(certificate['intervals']['normalization_scale'])
    assert scale_lo <= expected_lo <= expected_hi <= scale_hi
    llo, lhi = interval(certificate['intervals']['L_at_2_real'])
    dlo, dhi = interval(certificate['intervals']['L_second_derivative'])
    assert dlo <= expected_lo*llo <= expected_hi*lhi <= dhi


def test_complete_exact_euler_input_and_source_hashes(certificate):
    data = certificate['euler_input']
    compressed = (ROOT / data['path']).read_bytes()
    raw = gzip.decompress(compressed)
    assert hashlib.sha256(compressed).hexdigest() == data['compressed_sha256']
    assert hashlib.sha256(raw).hexdigest() == data['sha256']
    cutoff = data['coefficient_cutoff']
    sieve = bytearray(b'\x01') * (cutoff+1)
    sieve[:2] = b'\x00\x00'
    for p in range(2, int(cutoff**0.5)+1):
        if sieve[p]:
            sieve[p*p::p] = b'\x00' * ((cutoff-p*p)//p+1)
    primes = [p for p, yes in enumerate(sieve) if yes]
    rows = [list(map(int, row.split())) for row in raw.splitlines()]
    assert [row[0] for row in rows] == primes
    assert len(rows) == data['prime_count']
    for p, *c in rows:
        assert len(c) == 6 and c[0] == 1
        if p not in (3, 83):
            assert c[5] == -p**5 and c[4] == -p**3*c[1] and c[3] == -p*c[2]
            assert abs(c[1]) <= 5*p and abs(c[2]) <= 10*p*p
    for source in certificate['run']['source_hashes']:
        # Historical input generators are archived when the public API evolves.
        path = source.get('archived_path', source['path'])
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == source['sha256']
