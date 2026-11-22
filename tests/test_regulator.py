"""Regulators compared with plane integrals and survey values."""
import json
from pathlib import Path
import numpy as np
import pytest
from sage.all import QQ, RealField
from hyperell_regulator.regulator import (regulators_stokes, periods_and_regulator_stokes,
    periods_stokes, regulator_stokes, IntegrationMonitor, plot_cut_system)
from hyperell_regulator.regulator.polynomials import ring


@pytest.mark.parametrize('expression,periods,regulator', [
 ('x^5+x^2+x',[9.22548668,-1.22108503,8.98356248],[-5.80825813,.45209334,5.53626738]),
 ('x^5-5*x^3+4*x',[3.35349283,0,6.70698566],[-.89959203,0,6.44811227]),
 ('x^5+(x+1)^2',[6.05371274,-1.63861909,8.17585620],[-1.301698964,.50760058,5.770923704]),
])
def test_matches_independent_plane_integrals(expression, periods, regulator):
    F=ring()(expression)
    with IntegrationMonitor(max_edge_seconds=15) as monitor:
        p,r=periods_and_regulator_stokes(F)
    np.testing.assert_allclose(p,periods,atol=2e-8,rtol=0)
    np.testing.assert_allclose(r,regulator,atol=2e-8,rtol=0)
    assert all(v['degree']==2 and v['status']=='ok' for v in monitor.records)


def test_entire_saved_survey_at_twenty_digits():
    from hyperell_regulator.motivic import DivisorRegulators
    from hyperell_regulator.experiments.geometric_torsion import decode_point
    from hyperell_regulator.experiments._survey_divisors import decimal_value
    rows=json.loads((Path(__file__).parents[1]/'reports/beilinson_100_digits.json').read_text())
    R=RealField(150)
    successful=[r for r in rows if r.get('status')=='ok']
    assert len(successful)==41
    for row in successful:
        F=ring()(row['model'].replace('u','x'))
        data=row['divisor_data']
        points=data['support']
        abscissae=[decode_point(P['point'])[0] for P in points]
        divisors=[[P['coefficients'][name] for P in points] for name in data['matrix_divisor_ids']]
        det,M=DivisorRegulators(F,abscissae,digits=20).determinant(*divisors)
        error=max(abs(R(a)-R(decimal_value(b)))/max(1,abs(R(decimal_value(b))))
                  for rr,ss in zip(M,row['matrix']) for a,b in zip(rr,ss))
        assert error<R('1e-16'),(row['label'],error)
        assert abs(R(det)/R(decimal_value(row['det']))-1)<R('1e-16'),row['label']


def test_precision_growth_and_rational_logarithmic_center():
    x=ring().gen();q=QQ(1)/3;F=(x-q)**5+(x-q+1)**2
    p20,r20=periods_and_regulator_stokes(F,q,prec=20)
    p40,r40=periods_and_regulator_stokes(F,q,prec=40)
    R=r40[0].parent()
    for lo,hi in ((p20,p40),(r20,r40)):
        assert max(abs(R(a)-b) for a,b in zip(lo,hi))<R('1e-19')
    _,r=periods_and_regulator_stokes(x**5+(x+1)**2,0,prec=40)
    # x=u+q transforms the differential basis, not the scalar logarithm.
    assert abs(r40[0]-r[0])<R('1e-37')
    assert abs(r40[1]-(r[1]+q*r[0]))<R('1e-37')
    assert abs(r40[2]-(r[2]+2*q*r[1]+q*q*r[0]))<R('1e-37')


def test_reciprocal_curve_controls_infinity_and_log_sign():
    x=ring().gen();F=x**6+1
    p,r=periods_and_regulator_stokes(F,0,prec=30)
    R=r[0].parent()
    assert abs(p[0]-p[2])<R('1e-28')
    assert abs(r[0]+r[2])<R('1e-28')
    assert abs(r[1])<R('1e-28')


def test_period_only_and_shared_weights_agree_with_separate_cuts():
    x=ring().gen();F=x**6+4*x**5+4*x**4+2*x**3+1
    p,columns=regulators_stokes(F,[-1,0],prec=25)
    R=p[0].parent()
    assert max(abs(a-b) for a,b in zip(p,periods_stokes(F,prec=25)))<R('1e-23')
    for a,col in zip((-1,0),columns):
        assert max(abs(v-w) for v,w in zip(col,regulator_stokes(F,a,prec=25)))<R('1e-23')


def test_plot_and_legacy_polynomial_coefficient_input(tmp_path):
    F=ring().gen()**5+ring().gen()**2+ring().gen()
    p,r=periods_and_regulator_stokes([1,0,0,1,1,0])
    np.testing.assert_allclose(p,periods_stokes(F),atol=1e-12)
    filename=tmp_path/'cut.png'
    fig,_=plot_cut_system(F,[0,1],path=filename)
    assert filename.stat().st_size>1000
    import matplotlib.pyplot as plt
    plt.close(fig)
