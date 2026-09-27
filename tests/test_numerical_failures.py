"""Regressions for endpoint coordinates and false polynomial convergence."""
import numpy as np
from unittest.mock import patch
from hyperell_regulator.regulator import IntegrationMonitor, function_with_divisor
from hyperell_regulator.regulator import _legs
from hyperell_regulator.regulator.cover import Cover
from hyperell_regulator.regulator.polynomials import ring
from hyperell_regulator.regulator.walks import Chart, charts_over


def test_centered_high_order_zero_keeps_exact_factors():
    from sage.all import ComplexField
    x=ring().gen()
    F=x**6+2*x**5+3*x**4+4*x**3-x*x-2*x+1
    cov=Cover(F,*function_with_divisor(F,(-1,0),'infinity',15))
    chart=charts_over(cov,0)[0]
    assert chart.kind=='x' and chart.c==-1
    local=chart.local
    assert local.exact is not None
    num,den=local.norm_exact()
    assert num==num.leading_coefficient()*x**15 and den.degree()==0
    CC=ComplexField(120)
    for u in (1e-3+2e-3j,1e-6+2e-6j):
        z=CC(u);fq,a,ad,b,bd=local.exact;eta=fq(z).sqrt()
        expected=a(z)/ad(z)+b(z)/bd(z)*eta
        actual=local.phi_stable(u,complex(eta))
        assert abs(complex(actual)/complex(expected)-1)<1e-12


def test_cancelled_finite_pole_inverse_differential():
    x=ring().gen();F=x*(x-1)*(x**3+x+1)
    cov=Cover(F,*function_with_divisor(F,(0,0),(1,0),2))
    for sigma in (False,True):
        ch=Chart(cov,'inv',2,sigma=sigma)
        for xx in (0.4+0.2j,1.3+0.1j):
            y=np.sqrt(complex(F(xx)));u=complex(ch.to_u(xx));eta=complex(ch.eta_of(u,y))
            w=complex(cov.phi_stable(xx,y));z=1/w if sigma else w
            num,den=ch.integrand(u,eta,z)
            expected=np.array([1,xx])/cov.y_dphi_dx(xx,y)
            actual=num/den*(-z*z if sigma else 1)
            np.testing.assert_allclose(actual,expected,atol=1e-12,rtol=1e-12)


def test_long_edge_small_relative_standoff_is_attempted():
    # A nearby endpoint on a very long edge needs a small relative fraction.
    class Cov:
        def endpoint_ramification(self,*args):return 2
        def fibre(self,z):return np.array([1]),np.array([1])
    zero=_legs.Leg(np.zeros((2,1),complex),np.zeros((2,1),complex),None)
    def leg(cov,wa,wb,x,y,lam,*args,**kwargs):
        assert abs(wb-2)<1e-3 and wb!=2
        return x,y,lam,zero
    with patch.object(_legs,'leg_in_w',side_effect=leg) as calls:
        with patch.object(_legs,'end_leg',return_value=zero):
            out=_legs.side(Cov(),1e9,np.array([1]),np.array([1]),np.array([0]),2,'simple',False,hold=1e-4)
    assert calls.call_count==2 and out.frac<1e-7


def test_conductor_691_failed_edge_converges_and_trace_closes():
    x=ring().gen();F=4*x**5-4*x**3-3*x*x+2*x+1
    cov=Cover(F,*function_with_divisor(F,(0,-1),'infinity',8))
    a=0.05639794211475057+0.0029352394840220578j;b=-0.11897341953289166
    with IntegrationMonitor(max_edge_seconds=15) as monitor:
        _,leg=_legs.edge_leg(cov,a,b,np.log((a+b)/2),'simple','simple')
    assert monitor.records[0]['status']=='ok'
    assert np.max(abs(leg.dI.sum(axis=1)))<1e-9
    assert np.max(abs(leg.dL.sum(axis=1)))<1e-8


def test_uniformizer_interpolation_preserves_small_endpoint():
    class Cov:
        g=1
        def critical_values(self): return np.array([],complex)
    def quadrature(cov,coord,rate,x,y,lam,charts,level,**kwargs):
        assert abs(coord(1)/1e-8-1)<1e-14
        assert coord(0)==1e24
        assert np.isfinite(rate(1))
        return x,y,lam,np.zeros((1,1)),np.zeros((1,1))
    with patch.object(_legs.walks,'charts_over',return_value=[]):
        with patch.object(_legs,'_walk_quadrature',side_effect=quadrature):
            for toward in (None,(0,2)):
                _legs.leg_in_w(Cov(),1e24,1e-8,[0],[1],[0],toward=toward)
