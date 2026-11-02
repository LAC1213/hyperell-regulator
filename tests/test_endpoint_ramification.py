"""Endpoint ramification and coordinate identities for float64 continuation."""

from functools import lru_cache
from unittest.mock import patch

import numpy as np
import pytest
from sage.all import QQ, PolynomialRing

from hyperell_regulator.regulator import function_with_divisor
from hyperell_regulator.regulator.cover import Cover

R = PolynomialRing(QQ, "x")
x = R.gen()
F249 = x**6 + 4*x**5 + 4*x**4 + 2*x**3 + 1


@lru_cache(maxsize=1)
def _cover249():
    return Cover(F249, *function_with_divisor(F249, (0, -1), "infinity", 14))


def test_float_and_high_precision_vertices_share_the_exact_branch_polynomial():
    cov = Cover(x**5+x**2+2*x+1, x+1, 1)
    with patch.object(cov, 'exact_polynomial', wraps=cov.exact_polynomial) as build:
        low = cov.critical_values()
        high = cov.critical_values_hp(160)
        for value in low:
            if abs(value)>1e-8:
                assert cov.endpoint_ramification(value, 'simple') == 2
        assert build.call_count == 1
    assert len(low) == len(high)
    for value in high:
        assert min(abs(complex(value)-z) for z in low) < 1e-12


def test_weierstrass_zero_of_B_remains_a_true_branch_value():
    cov = Cover(x*(x**4+x+1), x+1, x)
    assert cov._branch_polynomial()(1) == 0
    assert 1 in cov.critical_values()
    assert cov.endpoint_ramification(1) == 2


def test_rational_endpoint_is_not_replaced_by_its_rounded_regular_fibre():
    value = QQ(1)/3 + QQ(1)/10**30
    cov = Cover(x*(x**4+x+1), x+value, x)
    # These share a float64 representation but define different exact fibres.
    assert max(p['e'] for p in cov.places_over(float(value))) == 1
    assert max(p['e'] for p in cov.places_over(value)) == 2
    assert cov.endpoint_ramification(complex(value)) == 2
    with pytest.raises(NotImplementedError):
        cov.endpoint_ramification(np.nextafter(float(value), np.inf))


def test_approximate_real_critical_values_keep_ramification():
    cov = _cover249()
    nodes, kinds, _ = cov.cut_path()
    for z, kind in zip(nodes, kinds):
        expected = cov.N if kind == "P" else 2 if kind == "simple" else 1
        assert cov.endpoint_ramification(z, kind) == expected
    # Rationalizing this approximation changes the fibre to a regular one.
    z = -1.8817501026561414
    assert max(p["e"] for p in cov.places_over(z)) == 1
    assert cov.endpoint_ramification(z, "simple") == 2


def test_exact_higher_ramification_is_not_inferred_from_nearby_float():
    cov = Cover(x**5 + x + 1, x**3 + 2, 0)
    assert cov.endpoint_ramification(2, "simple") == 3
    try:
        cov.endpoint_ramification(np.nextafter(2.0, 3.0), "simple")
    except NotImplementedError:
        pass
    else:
        raise AssertionError("a nearby regular fibre was treated as the exact critical fibre")


def test_mixed_ramification_uses_a_common_exponent():
    # Over w=2: x=0 gives two points of index three, and x=1 is a
    # Weierstrass point of index two. A common substitution needs s^6.
    cov = Cover((x - 1)*(x**4 + x + 2), x**3*(x - 1) + 2, 0)
    assert sorted(p["e"] for p in cov.places_over(2)) == [2, 3]
    assert cov.endpoint_ramification(2, "simple") == 6


def test_inexact_endpoint_requires_ramification_data():
    cov = Cover(x**5 + x + 1, x + 1, 1)
    cov.exact = None
    try:
        cov.endpoint_ramification(2.0, "simple")
    except NotImplementedError as exc:
        assert "exact coefficients or monodromy" in str(exc)
    else:
        raise AssertionError("inexact data was used as an exact certificate")
    # Exact place data may be supplied independently of the coefficients.
    cov.places_over = lambda z: [{"e": 2}, {"e": 3}]
    assert cov.endpoint_ramification(None, "Q") == 6


def test_reciprocal_coordinate_preserves_cover_and_norm():
    cov = _cover249()
    for center in (0, 2, 0.25 + 0.5j):
        local = cov.reciprocal(center)
        assert local.N == cov.N
        for u in (0.35 + 0.17j, -0.72 + 0.11j):
            xx = center + 1/u
            yy = np.sqrt(np.polyval(cov.F, xx))
            eta = u**(cov.g + 1) * yy
            w = complex(cov.phi_stable(xx, yy))
            assert abs(np.polyval(local.F, u) - eta*eta) < 1e-11 * max(1, abs(eta*eta))
            assert abs(complex(local.phi_stable(u, eta)) - w) < 1e-9 * max(1, abs(w))
            assert abs(local.norm_value(u) / cov.norm_value(xx) - 1) < 1e-10
            lhs = np.polyval(local.poly_in_x(w), u)
            scale = np.polyval(abs(local.poly_in_x(w)), abs(u))
            assert abs(lhs) < 1e-12 * scale
    inv = cov.reciprocal()
    numerator, denominator = inv.norm_exact()
    assert numerator == 4 and denominator == x**14
    # The unreduced expression u^42/u^56 would underflow here.
    assert abs(inv.norm_value(1e-20) / 4e280 - 1) < 1e-14


def test_reciprocal_differential_is_regular_at_weierstrass_point():
    cov = _cover249()
    inv = cov.reciprocal()
    u, eta = -1.0, 0.0
    expected = -u**(cov.g - 1) * cov.y_dphi_dx(-1.0, 0.0)
    assert expected == -12
    assert abs(inv.y_dphi_dx(u, eta) - expected) < 1e-12
