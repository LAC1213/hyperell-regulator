"""Full-precision map derivatives and continuation through cached samples."""
import pytest
from sage.all import QQ, PolynomialRing

from hyperell_regulator.regulator.cover import Cover
from hyperell_regulator.regulator._mp_continuation import Arithmetic, Chart, Walker
from hyperell_regulator.regulator._mp_stokes import Integrator


@pytest.mark.parametrize('mapping', ['affine', 'polynomial', 'constant', 'rational', 'pure'])
def test_mp_map_and_derivative_match_exact_rational_functions(mapping):
    x=PolynomialRing(QQ,'x').gen();ar=Arithmetic(100)
    F=x**5+x+1
    A,B={'affine': (x-1,0), 'polynomial': (3*x*x-1,1),
         'constant': ((x**3+1)/3,QQ(2)/3),
         'rational': ((x*x+1)/(x-1),(2*x+1)/(x+2)),
         'pure': (x**3+2,0)}[mapping]
    cov=Cover(F,A,B)
    for kind,center in [('x',0),('x',1),('inv',2)]:
        for sigma in (False,True):
            ch=Chart(cov,ar,kind,center,sigma=sigma,zr=ar.C('.125','.25'))
            u=ar.C('.25','.4');eta=ch.F(u).sqrt()
            f,a,ad,b,bd=ch.local.exact
            a,b=a/ad,b/bd
            expected=a(u)+b(u)*eta
            derivative=eta*a.derivative()(u)+eta**2*b.derivative()(u)+b(u)*f.derivative()(u)/2
            value,den=ch.values(u,eta)
            assert abs(value-expected)<ar.R('1e-110')*max(1,abs(expected))
            assert abs(den-derivative)<ar.R('1e-110')*max(1,abs(derivative))


def test_mp_continuation_preserves_monodromy_after_restoring_samples():
    x=PolynomialRing(QQ,'x').gen();ar=Arithmetic(100)
    engine=Integrator(Cover(x**5+2,x**3+2,0),ar)
    coord=lambda t:2+ar.R('.125')*(2*ar.C.gen()*ar.R.pi()*t).exp()
    points=engine.fibre(coord(0))
    walk=Walker(engine.chart,coord,points)
    half=list(walk.to(ar.R('.5')))
    walk.to(ar.R(1));expected=walk.global_points()
    rotation=(2*ar.C.gen()*ar.R.pi()/3).exp()
    for p,q in zip(points,expected):
        assert abs(q[0]-rotation*p[0])<ar.R('1e-110')
        assert abs(q[1]**2-(x**5+2)(q[0]))<ar.R('1e-110')
    # A cached sample can restore either backward or forward in the path.
    walk.restore(ar.R('.5'),half)
    walk.to(ar.R(1))
    for p,q in zip(expected,walk.global_points()):
        for v,w in zip(p,q):assert abs(v-w)<ar.R('1e-110')
    restarted=Walker(engine.chart,coord,points)
    restarted.to(ar.R('.25'))
    restarted.restore(ar.R('.5'),half)
    restarted.to(ar.R(1))
    for p,q in zip(expected,restarted.global_points()):
        for v,w in zip(p,q):assert abs(v-w)<ar.R('1e-110')


def test_mp_shifted_high_order_zero_retains_relative_accuracy():
    from hyperell_regulator.regulator import function_with_divisor

    x=PolynomialRing(QQ,'x').gen()
    F=x**6+2*x**5+3*x**4+4*x**3-x*x-2*x+1
    cov=Cover(F,*function_with_divisor(F,(-1,0),'infinity',15))
    ar=Arithmetic(100);ch=Chart(cov,ar,center=-1)
    reference=Arithmetic(180);hi=Chart(cov,reference,center=-1)
    u=ar.C('1e-10','2e-11');v=reference.C(u)
    value,den=ch.values(u,ch.F(u).sqrt())
    f,a,ad,b,bd=hi.local.exact
    a,b=a/ad,b/bd;eta=hi.F(v).sqrt()
    expected=a(v)+b(v)*eta
    derivative=eta*a.derivative()(v)+eta**2*b.derivative()(v)+b(v)*f.derivative()(v)/2
    assert abs(value/ar.C(expected)-1)<ar.R('1e-110')
    assert abs(den/ar.C(derivative)-1)<ar.R('1e-110')
