from sage.all import *

from regulator_integral_numpy import *

def family1(t):
    return [1, 2*t, t**2 + 2*t, 2*t**2 - 2, t**2 + 2*t, 2*t, 1]

def companion_curves(f):
    R = PolynomialRing(QQ, 'x')
    x = R.gen()
    g = f[0] * x**3 + f[1] * x**2 + (f[2] - 3*f[0])*x + f[3] - 2*f[1]
    g1 = R(g(1/x - 2)*x**3)
    g2 = R(g(1/x + 2)*x**3)
    S = PolynomialRing(QQ, ['x', 'y', 'z'])
    x, y, z = S.gens()
    an = list(g1)
    c1 = y**2*z - sum([an[i] * x**i * z**(3-i) for i in range(4)])
    E1 = EllipticCurve(c1, [0,1,0])
    an = list(g2)
    c2 = y**2*z - sum([an[i] * x**i * z**(3-i) for i in range(4)])
    E2 = EllipticCurve(c2, [0,1,0])
    return (E1, E2)

def tensor_euler_factor(E1, E2, p, bad_factors = None):
    """
    Return Euler factor at prime p for L(E1 ⊗ E2, s)
    bad_factors: dict {p: polynomial in T} overriding default
    """
    R = PolynomialRing(QQ, 'T')
    T = R.gen()

    # Custom bad factor override
    if bad_factors and p in bad_factors:
        return bad_factors[p]

    a1 = E1.ap(p)
    a2 = E2.ap(p)
    
    # The following gets the correct factors if both E1 and E2 
    # have either good or semistable reduction at p
    if E1.has_good_reduction(p):
        if E2.has_good_reduction(p):
            return (
                1
                - a1*a2*T
                + (a1**2 + a2**2 - 2*p)*p*T**2
                - p**2*a1*a2*T**3
                + p**4*T**4
            )
        else:
            return 1 - a1*a2*T + p*T**2
    else:
        if E2.has_good_reduction(p):
            return 1 - a1*a2*T + p*T**2
        else:
            return (1 - a1*a2*T)*(1 - a1*a2*p*T)
            
def tensor_dirichlet_series(E1, E2, N, bad_factors=None):
    """
    Compute Dirichlet coefficients a_n for n ≤ N
    using Euler product expansion.
    """
    A = [0]*(N+1)
    A[1] = 1

    for n in range(2, N+1):
        fac = factor(n)
        coeff = 1

        for p, e in fac:
            R = PowerSeriesRing(QQ, 'x', default_prec=e+1)
            x = R.gen()

            P = tensor_euler_factor(E1, E2, p, bad_factors)

            # substitute T = x
            P_series = R(P(x))

            inv = 1 / P_series
            coeff *= inv[e]

        A[n] = coeff

    return A

def tensor_lfun(E1, E2):
    assert is_squarefree(E1.conductor()) and is_squarefree(E2.conductor()), "tensor L-function only implemented for semistable reduction"
    
    cond = prod([p**2 for (p, _) in factor(E1.conductor()*E2.conductor())])
    L = Dokchitser(cond, [-1, 0, 0, 1], 3, 1)
    N = L.num_coeffs()
    L.init_coeffs(tensor_dirichlet_series(E1, E2, N)[1:])
    return L

def palindrome_regulator(f):
    integrand = lambda z : (1 - np.abs(z**2)) * np.log(np.abs(z))
    return hyperell_integral(f, lambda z : integrand(z))

maxt = 30

with open('palindromic_beilinson.csv','w') as file:
    for t in range(-maxt,maxt + 1):
        if t == -2 or t == 0:
            continue
        f = family1(t)
        E1, E2 = companion_curves(f)
        if not is_squarefree(E1.conductor()) or not is_squarefree(E2.conductor()):
            continue
        print(f"t = {t} ")
        print("companion curves:")
        print(E1.minimal_model())
        print(E2.minimal_model())
        L = tensor_lfun(E1, E2)
        print("L-fun check:", L.check_functional_equation())
        l2 = L(2)
        l1 = L.derivative(1,1)
        print("L(2) = ", l2)
        print("L'(1) = ", l1)
        reg, err = palindrome_regulator(f)
        print("reg, err = ", reg, err)
        rel = list(gp.lindep([RR(reg/pi), l1], 8))
        print("lindep [reg/pi, L'(1)]: ", rel)
        err2 = RR(reg/pi) * Integer(rel[0]) + RR(l1) * Integer(rel[1])
        print("lindep err: ", err2)
        print(60*"=")
        
        file.write(f"{t}, {RR(l1)}, {RR(l2)}, {RR(reg/pi)}, {Integer(rel[0])}, {Integer(rel[1])}, {max(err,err2)}\n")

