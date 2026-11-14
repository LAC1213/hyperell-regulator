"""Precision-preserving reports, caches and cycle periods."""
import json
from unittest.mock import patch
from sage.all import QQ, RealField
from hyperell_regulator.regulator import c_invariants
from hyperell_regulator.regulator.polynomials import ring
from hyperell_regulator.experiments import beilinson as b

def test_report_preserves_decimal_digits_and_cli_precision():
    x=ring().gen();F=x**5+x+1;R=RealField(150)
    value=R('1.23456789012345678901234567890123456789')
    triple=dict(P=(0,1),Q=(1,1),R=('infinity',1),order_P=5,order_Q=5)
    with patch.object(b,'_model_for_triple',return_value=(F,[(0,1),(1,1)],None,'infinity')),patch.object(b,'_resolve_triple_orders'),patch.object(b,'beilinson_determinant_stokes',return_value=(value,[[value]*3]*3,None)):
        out=b._analyse_triple({'label':'test','cond':1},None,triple,1,
            verbose=False,reg_digits=30,with_periods=False,with_lvalue=False)
    saved=json.loads(json.dumps(out,default=str))
    assert isinstance(saved['det'],str)
    assert abs(R(saved['det'])-value)<R('1e-30')
    assert abs(R(saved['matrix'][0][0])-value)<R('1e-30')
    with patch.object(b,'run',return_value=[]) as run:
        b.main(['250','--prec','30','--edge-seconds','45'])
    assert run.call_args.kwargs['reg_digits']==30
    assert run.call_args.kwargs['max_edge_seconds']==45


def test_lvalue_cache_retains_precision_and_rejects_old_float_records():
    from tempfile import TemporaryDirectory
    from pathlib import Path
    from unittest.mock import Mock
    R=RealField(140);value=R('2.345678901234567890123456789012345')
    lf=Mock()
    lf._coeffs_needed.return_value=100
    lf.conductor.return_value=249
    lf.check_functional_equation.return_value=R('1e-40')
    lf.value.return_value=value
    with TemporaryDirectory() as d:
        path=Path(d)/'L.json'
        path.write_text(json.dumps({'curve':dict(L_second_derivative=float(value),
            L_pari_precision=True,L_digits=30,L_num_coeffs=100)}))
        with patch.object(b,'StandardLFunction',return_value=lf):
            result=b.standard_lvalue_cached('curve',None,digits=30,path=path)
            assert lf.value.call_count==1
            assert abs(R(result['L_second_derivative'])-value)<R('1e-35')
            cached=b.standard_lvalue_cached('curve',None,digits=30,path=path)
            assert cached['L_cached'] and lf.value.call_count==1


def test_709_periods_continue_across_principal_sqrt_cut():
    x=ring().gen();F=4*x**5-7*x**4+4*x
    cp,cm=c_invariants(F,prec=35)
    # Full periods on the MST cycle basis.
    assert abs(cp-cp.parent()('-4.5902904081140363700932752793993016'))<cp.parent()('1e-32')
    for lo,hi in zip(c_invariants(F),(cp,cm)):
        assert abs(lo-float(hi))<1e-10
