"""Regression checks for conditioned continuation and long finite legs."""
import json
from unittest.mock import patch
import numpy as np
from hyperell_regulator.experiments import beilinson as b
from hyperell_regulator.regulator import IntegrationMonitor, function_with_divisor, _legs
from hyperell_regulator.regulator.cover import Cover
from hyperell_regulator.regulator.polynomials import ring


def _curve389():
    x=ring().gen()
    return x**6+4*x**5-6*x**4-32*x**3+x*x+64*x+28


def test_nearby_weierstrass_point_does_not_change_endpoint_chart():
    F=_curve389()
    cov=Cover(F,*function_with_divisor(F,(-2,0),'infinity',10))
    xt=-2.1220562189328414
    yt=0.02507982151474817
    assert _legs.chart_at(cov,xt,yt)=='x'
    assert _legs.chart_at(cov,-2,0)=='y'
    # A branch-point-free straight x-tail cannot change y sign. Place
    # the start close enough to isolate this check from the disk guard.
    x0=xt+1e-6
    local=cov.shifted(-2)
    y0=-np.sqrt(complex(np.polyval(local.F,x0+2)))
    w0=complex(local.phi_stable(x0+2,y0))
    with IntegrationMonitor() as monitor:
        with monitor.edge('wrong-sheet endpoint',degree=10):
            try:
                _legs.end_leg(cov,w0,x0,y0,
                    np.log(w0),2.7426115808649313e-5,
                    'simple',False)
            except _legs.LegFailure as exc:
                assert 'opposite y-sheet' in str(exc)
            else:
                raise AssertionError('accepted a tail ending on the opposite sheet')


def _check_curve(label,F,P,nP,Q,nQ,expected,edges):
    periods=[]
    original=b._regulator_column_and_periods
    def column(*args,**kw):
        col,per=original(*args,**kw)
        periods.append(np.array(per,float))
        return col,per
    def progress(event):
        if event['event'] != 'start':
            print('EDGE',label,event['index'],event['event'],
                  round(event['seconds'],2),event['phase'],flush=True)
    with IntegrationMonitor(max_edge_seconds=15,max_edge_steps=100000,
                            callback=progress,progress_interval=3) as monitor:
        with patch.object(b,'_regulator_column_and_periods',column):
            det,M,_=b.beilinson_determinant_stokes(F,P,nP,Q,nQ)
    assert len(monitor.records)==edges
    assert all(e['status']=='ok' and e['seconds']<15 for e in monitor.records)
    np.testing.assert_allclose(periods[0],periods[1],atol=1e-9,rtol=1e-9)
    np.testing.assert_allclose(float(det),expected,atol=0,rtol=1e-7)
    traces={k:max(abs(complex(*v)) for e in monitor.records for v in e[k])
            for k in ('trace_periods','trace_log_periods')}
    assert max(traces.values())<1e-8
    print('CURVE_RESULT',json.dumps(dict(label=label,det=float(det),matrix=np.array(M,float).tolist(),
          period_difference=float(np.max(abs(periods[0]-periods[1]))),traces=traces,
          edges=edges,slowest_edge=max(e['seconds'] for e in monitor.records),
          integration_seconds=sum(e['seconds'] for e in monitor.records))),flush=True)


def test_conductor_389_full_regulator_with_trace_checks():
    _check_curve('389.a.389.1',_curve389(),(-2,0),10,(2,0),5,
                 1.95775597024077,58)


def test_conductor_603_full_regulator_with_trace_checks():
    x=ring().gen()
    F=64*x**5+72704*x**4-59392*x**3+114688*x*x+16384*x
    # The model normalization scales v; phi_P=272*x^2-128*x+256+v.
    _check_curve('603.a.603.1',F,(4,-4096),5,(0,0),2,
                 -3.427815790547324e-11,10)
