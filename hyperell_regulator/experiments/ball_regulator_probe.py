"""Experimental certified finite-edge integrals, not a full regulator backend.

Use Sage's Arb integrator on regularized differentials, with explicit bounds
for the two omitted logarithmic endpoint tails. Sheet choices are local to
each edge: global Stokes topology and the infinity ray are not certified here.
"""

from hyperell_regulator.regulator.balls import edge_integrals


def standard_second_from_value_at_two(value, conductor, *, sign=1):
    r"""Propagate an L(2) ball to L''(1), conditional on the functional equation.

    For gamma shifts (-1,-1,0,0,1) and Lambda(s)=sign*Lambda(3-s),
    L''(1) = sign*sqrt(N)/(4*pi^2)*L(2), provided Lambda is holomorphic
    at 1 and 2. This does not certify the supplied L(2) or the assumptions.
    """
    from sage.all import ZZ
    from sage.rings.real_arb import RealBall

    if not isinstance(value, RealBall):
        raise TypeError("value must be a real ball enclosure of L(2)")
    N = ZZ(conductor)
    if N <= 0 or sign not in (-1, 1):
        raise ValueError("conductor must be positive and sign must be +1 or -1")
    R = value.parent()
    return sign * R(N).sqrt() * value / (4 * R.pi()**2)
