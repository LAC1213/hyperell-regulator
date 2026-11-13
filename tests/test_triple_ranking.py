"""Cut conditioning is independent of torsion orders."""
from unittest.mock import patch
from hyperell_regulator.experiments import beilinson as b
from hyperell_regulator.regulator.polynomials import ring


def test_torsion_orders_do_not_change_cut_score():
    x=ring().gen();F=x**6+4*x**5+4*x**4+2*x**3+1
    low=dict(order_P=7,order_Q=14)
    high=dict(order_P=700,order_Q=1400)
    with patch.object(b,'_model_for_triple',return_value=(F,[(-1,0),(0,-1)],None,'infinity')),patch.object(b,'_resolve_triple_orders'):
        assert b.triple_conditioning(None,low)==b.triple_conditioning(None,high)
