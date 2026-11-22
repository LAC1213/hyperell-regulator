"""Logarithmic regulator weights with finite common poles."""
import numpy as np
from hyperell_regulator.experiments import beilinson as b
from hyperell_regulator.regulator import IntegrationMonitor
from hyperell_regulator.regulator.polynomials import ring


# Exact LMFDB equations y^2 + h(x)y = f(x), with coefficients in
# ascending degree. Keep the regression inputs independent of local caches.
CURVE_249 = {"label": "249.a.6723.1",
             "eqn": [[2, 3, 1, 1, 0, -1], [1, 0, 0, 1]]}
CURVE_353 = {"label": "353.a.353.1",
             "eqn": [[0, 0, 1], [1, 1, 0, 1]]}


def test_finite_pole_column_equals_difference_of_infinity_columns():
    x=ring().gen();F=x*(x-1)*(x**3+x+1)
    with IntegrationMonitor(max_edge_seconds=15) as monitor:
        finite=b.regulator_column_stokes(F,(0,0),pole=(1,0))
        p=b.regulator_column_stokes(F,(0,0))
        r=b.regulator_column_stokes(F,(1,0))
    np.testing.assert_allclose(np.asarray(finite,float),np.asarray(p,float)-np.asarray(r,float),atol=1e-8,rtol=0)
    assert all(e['status']=='ok' for e in monitor.records)
