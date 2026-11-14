"""Accuracy checks and work budgets for explicit quadrature."""
import pytest
from hyperell_regulator.regulator.quadrature import Arithmetic, quadrature, QuadratureFailure
from hyperell_regulator.regulator import IntegrationMonitor, IntegrationLimitError, regulator_stokes
from hyperell_regulator.regulator.polynomials import ring


def test_singular_endpoint_integrals_and_cached_refinement():
    ar=Arithmetic(30);calls=[]
    def sample(t):
        calls.append(t)
        return [ar.C(t.sqrt()),ar.C(t.log())]
    result=quadrature(ar,sample,2)
    assert abs(result[0]-ar.R(2)/3)<ar.R('1e-29')
    assert abs(result[1]+1)<ar.R('1e-29')
    assert len(calls)==len(set(calls))


def test_nearby_complex_poles_require_refinement_and_failure_is_explicit():
    ar=Arithmetic(20);w=ar.R('.012');c=ar.R('.43')
    sample=lambda t:[ar.C(1)/(1+((t-c)/w)**2)]
    expected=w*((1-c)/w).arctan()+w*(c/w).arctan()
    with pytest.raises(QuadratureFailure):
        quadrature(ar,sample,1,max_level=4)
    assert abs(quadrature(ar,sample,1)[0]-expected)<ar.R('1e-20')


def test_rules_stay_strictly_inside_interval():
    ar=Arithmetic(100)
    assert all(0<t<1 and w>0 for t,w in ar.rule(5))


def test_budget_interrupts_actual_integration():
    x=ring().gen()
    with pytest.raises(IntegrationLimitError):
        with IntegrationMonitor(max_edge_steps=2,max_edge_seconds=None) as monitor:
            regulator_stokes(x**5+x+1)
    assert monitor.records[0]['status']=='limit'
    assert monitor.records[0]['steps']==3
