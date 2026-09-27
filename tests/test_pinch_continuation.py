"""Projection collisions are ordinary points on the smooth curve."""
import numpy as np
from hyperell_regulator.regulator import IntegrationMonitor
from hyperell_regulator.regulator import _legs
from hyperell_regulator.regulator.cover import Cover
from hyperell_regulator.regulator.polynomials import ring


def test_unramified_pinch_agrees_with_straight_curve_ode():
    x=ring().gen()
    # At x=0 the two points y=+/-1 both map to w=2; their derivatives
    # are 3 and 1. They must remain distinct throughout continuation.
    cov=Cover(x**5+x+1,2+2*x,x)
    a,b=1.95,2.05
    xs,ys=cov.fibre(a);lam=np.full(cov.N,np.log(a),complex)
    with IntegrationMonitor(max_edge_seconds=15) as monitor:
        with monitor.edge('pinch',degree=cov.N):
            xe,ye,le,actual=_legs.leg_in_w(cov,a,b,xs,ys,lam)
        with monitor.edge('independent ODE',degree=cov.N):
            xo,yo,lo,reference=_legs._leg_in_w_ode(cov,a,b,xs,ys,lam)
        with monitor.edge('reverse pinch',degree=cov.N):
            xr,yr,lr,reverse=_legs.leg_in_w(cov,b,a,xe,ye,le)
    assert monitor.records[0]['counters']['pinch_detours']==1
    np.testing.assert_allclose(xe,xo,atol=1e-10,rtol=0)
    np.testing.assert_allclose(ye,yo,atol=1e-10,rtol=0)
    np.testing.assert_allclose(actual.dI,reference.dI,atol=1e-10,rtol=0)
    np.testing.assert_allclose(actual.dL,reference.dL,atol=1e-10,rtol=0)
    np.testing.assert_allclose(xr,xs,atol=1e-10,rtol=0)
    np.testing.assert_allclose(yr,ys,atol=1e-10,rtol=0)
    np.testing.assert_allclose(actual.dI+reverse.dI,0,atol=1e-10,rtol=0)
    np.testing.assert_allclose(actual.dL+reverse.dL,0,atol=1e-10,rtol=0)


def test_pinch_detour_does_not_cross_critical_value_or_zero():
    x=ring().gen();cov=Cover(x**5+x+1,2+2*x,x)
    path=_legs._pinch_detour(cov,1.95,2.05)
    assert path[0]==1.95 and path[-1]==2.05
    assert min(abs(z-2) for z in path)>0
    radius=max(abs(z-2) for z in path[1:-1])
    assert radius<min(abs(z-2) for z in np.r_[cov.critical_values(),0j])
    assert _legs._pinch_detour(cov,1.95+1j,2.05+1j) is None


def test_conductor_295_regulator_converges_with_bounded_edges():
    from unittest.mock import patch
    from hyperell_regulator.experiments import beilinson as b
    x=ring().gen();F=x**6+2*x**3-4*x*x+1
    periods=[];original=b._regulator_column_and_periods
    def column(*args,**kwargs):
        col,per=original(*args,**kwargs)
        periods.append(np.asarray(per,float))
        return col,per
    with IntegrationMonitor(max_edge_seconds=15,max_edge_steps=100000) as monitor:
        with patch.object(b,'_regulator_column_and_periods',column):
            det,matrix,_=b.beilinson_determinant_stokes(F,(0,1),14,(1,0),7)
    assert len(monitor.records)==10
    assert all(e['status']=='ok' for e in monitor.records)
    np.testing.assert_allclose(periods[0],periods[1],atol=1e-9,rtol=0)
    assert abs(float(det)-18.0731942157733)<1e-8
    for e in monitor.records:
        for key in ('trace_periods','trace_log_periods'):
            assert max(abs(complex(*v)) for v in e[key])<1e-9
    print('295 det',float(det),'independent period difference',
          np.max(abs(periods[0]-periods[1])),
          'slowest edge',max(e['seconds'] for e in monitor.records),flush=True)


def test_ramified_pinch_is_not_treated_as_an_ordinary_collision():
    x=ring().gen();cov=Cover(x**5+x+1,2+2*x*x,x*x)
    # Here both points over x=0 have ramification index two over w=2.
    # Detouring them could change sheet labels and must not be automatic.
    assert _legs._pinch_detour(cov,1.95,2.05) is None


def test_stable_newton_recovers_a_small_control_step():
    from unittest.mock import patch
    from hyperell_regulator.regulator import walks
    x=ring().gen();cov=Cover(x**5+x+1,x,0)
    chart=walks.Chart(cov,'x',0)
    start=0.2;yy=np.sqrt(complex(cov.F[-1]+start+start**5))
    wk=walks.walker(cov,lambda t:start+0.1*t,start,yy,chart)
    wk.h=1e-6
    polynomial_newton=walks._newton
    def stalled_polishing(*args,**kwargs):
        u,_,res=polynomial_newton(*args,**kwargs)
        return u,9,res
    with patch.object(walks,'_newton',stalled_polishing):
        wk.to(1,budget=100)
    assert abs(wk.xy()[0]-0.3)<1e-12
    assert wk.steps<100


def test_preflight_chart_choices_are_reused_without_skipping_nodes():
    from unittest.mock import patch
    from hyperell_regulator.regulator import walks
    x=ring().gen();cov=Cover(x**5+x+1,x,0)
    a,b=0.2,0.3;xx,yy=cov.fibre(a);lam=np.full(cov.N,np.log(a),complex)
    _,_,_,preflight=_legs.leg_in_w(cov,a,b,xx,yy,lam,integrate=False)
    with patch.object(walks,'choose',side_effect=AssertionError('repeated chart search')):
        xe,ye,_,leg=_legs.leg_in_w(cov,a,b,xx,yy,lam,_prepared=preflight.charts)
    np.testing.assert_allclose(xe,b,atol=1e-12,rtol=0)
    np.testing.assert_allclose(ye*ye,complex((x**5+x+1)(b)),atol=1e-12,rtol=0)
    assert np.all(np.isfinite(leg.dI)) and np.max(abs(leg.dI))>0


def test_295_preflight_preserves_y_sheet_at_tangent_projection_collision():
    from hyperell_regulator.regulator import function_with_divisor
    x=ring().gen(); F=x**6+2*x**3-4*x*x+1
    cov=Cover(F,*function_with_divisor(F,(0,1),'infinity',14))
    end=0.0065551065034296085; a=end/2; b=end+0.08*(a-end)
    xx,yy=cov.fibre(a); lam=np.full(cov.N,np.log(a),complex)
    with IntegrationMonitor(max_edge_seconds=15) as monitor:
        with monitor.edge('295 preflight',degree=14):
            actual=_legs.leg_in_w(cov,a,b,xx,yy,lam,
                                 toward=(end,2),integrate=False)
        with monitor.edge('295 straight ODE',degree=14):
            reference=_legs._leg_in_w_ode(cov,a,b,xx,yy,lam)
    np.testing.assert_allclose(actual[0],reference[0],atol=1e-8,rtol=0)
    # Previously the preflight jumped to the opposite y-sheet (error 4.2).
    np.testing.assert_allclose(actual[1],reference[1],atol=1e-8,rtol=0)
