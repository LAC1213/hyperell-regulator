"""Regulator determinants from common edge periods."""
from unittest.mock import patch
import numpy as np
from hyperell_regulator.regulator import IntegrationMonitor
from hyperell_regulator.regulator.polynomials import ring
from hyperell_regulator.motivic import DivisorRegulators


def test_shared_columns_integrate_each_edge_once():
    from hyperell_regulator.regulator.stokes import Cut
    x=ring().gen();F=x**6+4*x**5+4*x**4+2*x**3+1
    with IntegrationMonitor(max_edge_seconds=15) as monitor:
        det,M=DivisorRegulators(F,[-1,0,'infinity']).determinant([1,0,-1],[0,1,-1])
    assert len(monitor.records)==len(Cut(F,[-1,0]).nodes)
    assert all(e['degree']==2 for e in monitor.records)
    assert abs(float(det)-33.48773151266456)<1e-12


def test_numerical_determinant_uses_divisor_coefficients():
    from hyperell_regulator import motivic
    per=np.array([2.,-1.,3.]);cols=[np.array([1.,2.,3.]),np.array([4.,5.,7.])]
    x=ring().gen();F=x**5+x+1
    with patch.object(motivic,'regulators_stokes',return_value=(per,cols)) as integrate:
        det,M=DivisorRegulators(F,[0,1,'infinity']).determinant([2,-1,-1],[-1,3,-2])
    integrate.assert_called_once_with(F,[0,1],prec=None)
    expected=np.array([1,2,1])[:,None]*np.column_stack([2*cols[0]-cols[1],-cols[0]+3*cols[1],per])
    np.testing.assert_array_equal(np.array(M),expected)
    assert abs(float(det)-np.linalg.det(expected))<1e-12
