"""Logarithmic regulator weights with finite common poles."""
import numpy as np
from hyperell_regulator.motivic import DivisorRegulators
from hyperell_regulator.regulator import IntegrationMonitor, regulator_stokes
from hyperell_regulator.regulator.polynomials import ring


def test_finite_pole_column_equals_difference_of_infinity_columns():
    x=ring().gen();F=x*(x-1)*(x**3+x+1)
    with IntegrationMonitor(max_edge_seconds=15) as monitor:
        _, (divisor,) = DivisorRegulators(F, [0, 1]).columns([[1, -1]])
        finite = regulator_stokes(F, 0, pole=1)
    np.testing.assert_allclose(np.asarray(divisor, float), np.asarray(finite)*[1, 2, 1], atol=1e-8, rtol=0)
    assert all(e['status']=='ok' for e in monitor.records)
