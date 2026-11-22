"""The two elliptic discs and rational traces on y^2=x^5-x at 2."""

from sage.all import NumberField, PolynomialRing, QQ

from hyperell_regulator import HyperellipticCurveWithReduction, MotivicPacket
from hyperell_regulator._dyadic_projected import _mean_clipped_depth
from hyperell_regulator.geometric_torsion import GeometricTorsionPacket


DATA = [22, [[[-1, 1], 1], [[0, 1], 1], [[1, 1], 1], [[1, 0, 1], 1], [[], 1],
             [[1, 2, 2, -2, 1], 3], [[1, -2, 2, 2, 1], 3]]]


def test_newton_average_accounts_for_an_exact_centre_and_different_depths():
    R=PolynomialRing(QQ, 'z');z=R.gen()
    K=NumberField(z, 'a');prime=K.primes_above(2)[0]
    # Distances from the centre 1 to the roots 1 and 3 are infinity and 1.
    assert _mean_clipped_depth((z-1)*(z-3), K(1), prime, QQ(2), QQ(1)/2) == QQ(3)/2
    assert _mean_clipped_depth((z-1)*(z-3), K(0), prime, QQ(2), QQ(1)/2) == QQ(1)/2


def test_bolza_full_galois_trace_basis_is_integral_and_complete():
    R=PolynomialRing(QQ, 'x');x=R.gen()
    C=HyperellipticCurveWithReduction(x**5-x)
    geometric=GeometricTorsionPacket(C, DATA, metadata={'complete':True})
    packet=MotivicPacket.from_geometric_packet(geometric)
    basis,complete=packet.integral_basis()
    assert complete and len(basis)==6
    solver=packet.integrality
    data,=solver.projected_local_data.values()
    assert data.method=='bolza_two_elliptic_discs'
    assert data.evidence['join_depth']==QQ(1)/2
    assert data.evidence['model'].graph.stabilize()[0].genera==(1,1)
    assert data.evidence['galois_trace_depths']==[QQ(1)/2]*len(solver.abscissae)
    assert all(solver.check(D.coefficients)['q_integral'] is True for D in basis)
    # All six Weierstrass points retract to the bridge midpoint.
    for i,P in enumerate(packet.points):
        if P[0]=='infinity' or not (2*P[1]+packet.h(P[0])):
            assert data.evidence['weights'][i]==QQ(1)/2
    # A single unresolved nonrational point is not replaced by a full orbit.
    unknown=data.unresolved[0]
    coefficients=[0]*len(packet.points)
    coefficients[solver._groups[unknown][0]]=1
    coefficients[next(i for i,P in enumerate(packet.points) if P[0]=='infinity')]=-1
    assert solver.check(coefficients)['q_integral'] is None
    assert not solver.kernel()['complete']
