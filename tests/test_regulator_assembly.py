"""Shared Stokes assembly, with independent plane values and call-count checks."""

from unittest.mock import patch

import numpy as np

from hyperell_regulator.regulator import periods_and_regulator_stokes
from hyperell_regulator.regulator import stokes
from hyperell_regulator.regulator.polynomials import ring

x = ring().gen()
F = x**5 + x**2 + x
PHI = (x, 0)
PERIODS = np.array([9.22548668, -1.22108503, 8.98356248])
REGULATOR = np.array([-5.80825813, 0.45209334, 5.53626738])


def test_combined_stokes_integrates_each_edge_once():
    # These values come from the independent plane quadrature, not from a
    # second Stokes assembly that could share the same algebraic mistake.
    with patch.object(stokes, "edge_data", wraps=stokes.edge_data) as edges:
        per, reg = periods_and_regulator_stokes(F, PHI)
    assert edges.call_count == 1
    np.testing.assert_allclose(per, PERIODS, rtol=1e-7, atol=1e-8)
    np.testing.assert_allclose(reg, REGULATOR, rtol=1e-7, atol=1e-8)


def test_combined_stokes_preserves_integrated_moments():
    with patch.object(stokes, "edge_data", wraps=stokes.edge_data) as edges:
        per, reg = periods_and_regulator_stokes(F, PHI, fast=False)
    assert edges.call_count == 1
    assert edges.call_args.args[-1] is True
    np.testing.assert_allclose(per, PERIODS, rtol=1e-7, atol=1e-8)
    np.testing.assert_allclose(reg, REGULATOR, rtol=1e-7, atol=1e-8)


def test_combined_stokes_preserves_precision_dispatch():
    from hyperell_regulator.regulator import precision

    expected_per, expected_reg = [1, 2, 3], [4, 5, 6]
    with patch.object(precision, "periods_and_regulator_stokes_mp",
                      return_value=(expected_per, expected_reg)) as combined:
        assert periods_and_regulator_stokes(
            F, PHI, prec=30, check_tol=1e-8) == (expected_per, expected_reg)
    combined.assert_called_once_with(F, PHI, 30, check_tol=1e-8)
    try:
        periods_and_regulator_stokes(F, PHI, fast=False, prec=30)
    except ValueError as exc:
        assert "only for the fast version" in str(exc)
    else:
        raise AssertionError("arbitrary precision must still reject fast=False")


def test_beilinson_determinant_integrates_two_covers():
    from hyperell_regulator.experiments.beilinson import beilinson_determinant_stokes

    # Both points are Weierstrass, with phi=x and x-1 of order two.  Their
    # norms are squares, so the normalized columns equal the raw regulators.
    f = x * (x - 1) * (x**3 + x + 1)
    p_per = np.array([2., -1., 3.])
    q_per = np.array([2. + 1e-12, -1., 3.])
    p_reg, q_reg = np.array([1., 2., 3.]), np.array([4., 5., 7.])
    # Mock the expensive integration only: divisor construction, normalization
    # of the columns, reuse of P's periods, and the determinant remain real.
    from hyperell_regulator import regulator
    with patch.object(regulator, "periods_and_regulator_stokes",
                      side_effect=[(p_per, p_reg), (q_per, q_reg)]) as integrate:
        det, matrix, error = beilinson_determinant_stokes(
            f, (0, 0), 2, (1, 0), 2)
    assert integrate.call_count == 2
    expected = np.column_stack([p_reg, q_reg, p_per])
    np.testing.assert_allclose(np.asarray(matrix, float), expected, rtol=0, atol=0)
    assert abs(float(det) - np.linalg.det(expected)) < 1e-12
    assert error is None
