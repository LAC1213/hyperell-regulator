from mpmath import mp, mpf, mpc, quad, sqrt, conj, log, pi, exp, sign, arg, inf
import numpy as np

from sage.interfaces.magma_free import magma_free
from sage.interfaces.gp import *
from sage.rings import real_mpfr
 
import requests
import pickle 

import homology

# set arbitrary precision
mp.dps = 1000

# --- helpers ---
def mp_polyval(f, x):
    """Evaluate polynomial with coefficients f (descending) at x (mpf/mpc)."""
    y = mpf('0')
    for c in f:
        y = y*x + mpf(c)
    return y

def lift_sqrt_along_path(f, xs, y0=None):
    ys = [None]*len(xs)
    val0 = mp_polyval(f, xs[0])
    ys[0] = sqrt(val0) if y0 is None else y0
    for k in range(1, len(xs)):
        val = mp_polyval(f, xs[k])
        r1 = sqrt(val)
        r2 = -r1
        ys[k] = r1 if abs(r1 - ys[k-1]) < abs(r2 - ys[k-1]) else r2
    return ys

# --- cycle integrals with tanh-sinh ---
def cycle_integrals(f, omega, basepoints, cycles, npts=200):
    integrals = [mpc(0) for _ in cycles]
    for j, loop_indices in enumerate(cycles):
        yprev = None
        for k in range(len(loop_indices)-1):
            z1 = basepoints[loop_indices[k]]
            z2 = basepoints[loop_indices[k+1]]
            xs = [z1 + (z2-z1)*mpf(t)/(npts-1) for t in range(npts)]
            ys = lift_sqrt_along_path(f, xs, yprev)
            xbranch = lambda t : z1 + (z2 - z1)*t
            def ybranch(t):
                yplus = sqrt(mp_polyval(f, xbranch(t)))
                m = min(int(t*npts), npts - 1)
                if abs(ys[m] - yplus) < abs(ys[m] + yplus):
                    return yplus
                else:
                    return -yplus
            yprev = ybranch(1)
            integrand = lambda t: omega(xbranch(t), ybranch(t))*(z2 - z1)
            integrals[j] += quad(integrand, [0,1], method='tanh-sinh')
    return integrals

def bilinear_form(B, v, w):
    s = mpc(0)
    for i in range(len(v)):
        for j in range(len(w)):
            s += mpc(B[i,j])*v[i]*w[j]
    return s

def hermitian_form(B, v, w):
    s = mpc(0)
    for i in range(len(v)):
        for j in range(len(w)):
            s += mpc(B[i,j])*v[i]*conj(w[j])
    return s

def partial_period(f, a, omega):
    path = lambda t : a - t/(1 - t)
    dpath = lambda t : -1/(1 - t)**2
    npts = 200
    xs = [path((k + 1)/(npts + 1)) for k in range(npts)]
    ys = lift_sqrt_along_path(f, xs)
    def ybranch(t):
        yplus = sqrt(mp_polyval(f, path(t)))
        m = min(int(t*npts), npts - 1)
        if abs(ys[m] - yplus) < abs(ys[m] + yplus):
            return yplus
        else:
            return -yplus
    integrand = lambda t : omega(path(t),ybranch(t)) * dpath(t)
    return quad(integrand, [0, 1], method='tanh-sinh')

def poly_to_str(coeffs, var="x"):
    """
    Convert a list of coefficients [a_d, a_{d-1}, ..., a_0]
    to a human-readable polynomial string.
    """
    terms = []
    degree = len(coeffs) - 1
    
    for i, a in enumerate(coeffs):
        power = degree - i
        if a == 0:
            continue  # skip zero coefficients

        # Determine sign
        if a < 0:
            sign = " - " if terms else "-"
            a = -a
        else:
            sign = " + " if terms else ""

        # Determine coefficient string
        if a == 1 and power != 0:
            coeff_str = ""
        else:
            coeff_str = str(a)

        # Determine power string
        if power == 0:
            power_str = ""
        elif power == 1:
            power_str = var
        else:
            power_str = f"{var}^{power}"

        term = f"{sign}{coeff_str}{power_str}"
        terms.append(term)

    return "".join(terms) if terms else "0"

class HyperEllCurve:
    
    def __init__(self, f, realroots = False, inithomology = True):
        self.f = f
        self.genus = (len(f) - 2)//2
        if inithomology:
            self.basepoints, self.cycles = homology.homology_basis(self.f, [], realroots = realroots)
            self.intersection_matrix = homology.cycles_intersection(self.f, self.basepoints, self.cycles)
            self.basepoints = [mpc(x) for x in self.basepoints]
        else:
            self.basepoints = None
            self.cycles = None
            self.intersection_matrix = None
        self._periods = None
        self._periods2 = None
    
    def _str_homology(self):
        str = "Homology basepoints:"
        str += f"\n {self.basepoints}"
        str += "\nHomology cycles:"
        for c in self.cycles:
            str += f"\n {c}"
        str += "\nIntersection matrix:"
        str += f"\n{np.matrix(self.intersection_matrix)}"
        return str
    
    def __str__(self):
        str = f"Genus {self.genus} curve: y^2 = " + poly_to_str(self.f)
        if not (self.cycles is None):
            str += "\n" + self._str_homology()
        if not (self._periods is None):
            str += "\nPeriod matrix:"
            for p in self._periods:
                str += f"\n {p}"
        return str
    
    def init_homology(self, avoid, realroots = False):
        self.basepoints, self.cycles = homology.homology_basis(self.f, avoid, realroots = realroots)
        self.intersection_matrix = homology.cycles_intersection(self.f, self.basepoints, self.cycles)
        self.basepoints = [mpc(x) for x in self.basepoints]
           
    def plot_cycles(self):
        homology.plot_cycles(self.f, self.basepoints, self.cycles)

    def periods(self):
        if self._periods is None:
            self._periods = [cycle_integrals(self.f, lambda x,y: x**i/y, self.basepoints, self.cycles) for i in range(self.genus)]
        return self._periods    
    
    def periods2(self):
        if self._periods2 is None:
            self._periods2 = [cycle_integrals(self.f, lambda x,y: x**(self.genus + i)/y, self.basepoints, self.cycles) for i in range(self.genus)]
        return self._periods2    
    
    def derhamrational_test(self):
        print("The following should be rational numbers: ")
        B = np.linalg.inv(self.intersection_matrix)
        for p1 in self.periods():
            for p2 in self.periods2():
                print(bilinear_form(B, p1, p2)/2/pi/mpc(1j))
                
    def jacobian_volume(self):
        ps = self.periods()
        B = np.linalg.inv(self.intersection_matrix)
        gram = [[hermitian_form(B, p1, p2) for p1 in ps] for p2 in ps]
        G = mp.matrix(gram)
        return mp.det(G)
    
    @classmethod
    def g2fromLMFDB(cls, label, realroots = False, inithomology = True):
        url = "https://www.lmfdb.org/api/g2c_curves/"
        params = {
            "_format": "json",
            "label": f"{label}"
        }

        response = requests.get(url, params)
        data = response.json()
        curve = data["data"][0]

        # Hyperelliptic equation coefficients
        eqn = eval(curve["eqn"])
        f = eqn[0]
        f.reverse()
        g = eqn[1]
        g.reverse()
        
        h = np.polyadd(np.polymul(f,[4]),np.polymul(g, g))
        return cls([int(c) for c in h], realroots, inithomology)


class PointedHyperEllCurve(HyperEllCurve):
    
    def __init__(self, f, points, realroots = False, inithomology = True):
        if isinstance(f, HyperEllCurve):
            self.f = f.f
            self.genus = f.genus
            self.basepoints = f.basepoints
            self.cycles = f.cycles
            self.intersection_matrix = f.intersection_matrix
            self._periods = f._periods
            self._periods2 = f._periods2
        else:
            super().__init__(f, realroots = realroots, inithomology=False)
            if inithomology:
                self.basepoints, self.cycles = homology.homology_basis(self.f, points, realroots = realroots)
                self.intersection_matrix = homology.cycles_intersection(self.f, self.basepoints, self.cycles)
                self.basepoints = [mpc(x) for x in self.basepoints]
        self.points = [mpc(p) for p in points]
        self._logperiods = None
        self._partialperiods = None

    def __str__(self):
        str = f"Pointed (x in {self.points}) genus {self.genus} curve: y^2 = " + poly_to_str(self.f)
        if not (self.cycles is None):
            str += "\n" + self._str_homology()
        if not (self._periods is None):
            str += "\nPeriod matrix:"
            for p in self._periods:
                str += f"\n {p}"
        if not (self._logperiods is None):
            str += "\nLog-Period matrix:"
            for p in self._logperiods:
                str += f"\n {p}"
        if not (self._partialperiods is None):
            str += "\nPartial periods:"
            for p in self._partialperiods:
                str += f"\n {p}"
        return str
    
    def save(self, name):
        with open(f"{name}.pickle", 'wb') as file:
            pickle.dump(self, file)
    
    @classmethod
    def read(cls, name):
        with open(f"{name}.pickle", 'rb') as file:
            return pickle.load(file)
    
    def logperiods(self):
        if self._logperiods is None:
            self._logperiods = []
            for a in self.points:
                self._logperiods.append([cycle_integrals(self.f, lambda x,y : log(x - a)*x**i/y, self.basepoints, self.cycles) for i in range(self.genus)])
        return self._logperiods
    
    def partialperiods(self):
        if self._partialperiods is None:
            self._partialperiods = []
            for a in self.points:
                self._partialperiods.append([partial_period(self.f, a, lambda x,y : x**i/y) for i in range(self.genus)])
        return self._partialperiods

    def beilinson_matrix(self):
        self.periods()
        self.logperiods()
        self.partialperiods()
        B = np.linalg.inv(self.intersection_matrix)
        
        entry = lambda v, w, i, j : hermitian_form(B, v[i], w[j]) + hermitian_form(B, v[j], w[i])
        pairs = [(0,0), (0,1), (1,1)]
        mat = [[entry(self._periods, self._periods, i, j).imag/2/pi for i,j in pairs]]
        for k in range(len(self.points)):
            bilinear_term = [entry(self._logperiods[k], self._periods, i, j)/2/pi/mpc(1j) for i,j in pairs]
            pp = self._partialperiods[k]
            membrane_term = [pp[i]*conj(pp[j]) + conj(pp[i])*pp[j] for i,j in pairs ]
            mat.append([(b + m).real for b, m in zip(bilinear_term, membrane_term)])
  
        return mp.matrix(mat)
    
    def test_beilinson(self, lvalue, lindep_digits, s = 1):
        assert s in [1,2], "s must be either 1 or 2" 
        vol = self.jacobian_volume().real
        beilinson = self.beilinson_matrix()
        print("Beilinson matrix: ")
        print(beilinson)
        beilinson_det = mp.det(beilinson).real
        
        if s == 1:
            print("L'(1): ", lvalue)
            print("Vol: ", vol)
            print("Reg: ", beilinson_det)
            for i in range(5, lindep_digits):
                a, b = gp.lindep([lvalue * vol * pi, beilinson_det], i)
                print("lindep [L'(1) * vol * pi, reg]: ", a, b)
                print("lindep error: ", int(a)* lvalue * vol * pi + int(b) * beilinson_det )
        else:
            print("L(2): ", lvalue)
            print("Vol: ", vol)
            print("Reg: ", beilinson_det)
            for i in range(5, lindep_digits):
                a, b = gp.lindep([lvalue * vol, pi**4 * beilinson_det], i)
                print("lindep [L(2) * vol, pi^4 * reg]: ", a, b)
                print("lindep error: ", int(a)* lvalue * vol + int(b) * pi**4 * beilinson_det )
                
def parse_mpc(s):
    """
    Parse a string like:
    '-1.91... + 4.80...e-1*1j'
    into mpmath.mpc, handling scientific notation.
    """
    s = s.replace(" ", "")  # remove all spaces

    # purely imaginary?
    if s.endswith("*1j"):
        s = s[:-3]  # remove '*1j'
        for i in range(1, len(s), 1):
            if s[i] in "+-" and s[i-1] != 'E':
                split_idx = i
                break
        else:
            # no '+' or '-', purely imaginary
            return mpc(mpf(0), mpf(s))

        real_part = s[:split_idx]
        imag_part = s[split_idx:]
        return mpc(mpf(real_part), mpf(imag_part))
    else:
        # purely real
        return mpc(mpf(s), mpf(0))

def test3125():
    f = [4, 0, 0, 0, 0, -1]
    C = PointedHyperEllCurve(f, [0, -1], inithomology=False)
    C.init_homology([0, -0.5, -1])
    print(C)

    lvalue = parse_mpc(magma_free("M := JacobiMotive([1/5,1/5,1/5], [3/5]); chi := Grossencharacter(M); L := LSeries(chi : Precision := 100); Evaluate(L,2);"))
    lvalue = lvalue.real

    gp("default(realprecision, 100)")
    ldir = mpf(str(gp("lfun(5,1)")))
    lvalue = lvalue * ldir

    #lvalue = mpf('0.35445162981482511890665852501155568022')
    C.test_beilinson(lvalue, 10, s = 2)
    
    C.save("Cond3125")

# https://www.lmfdb.org/Genus2Curve/Q/394/a/3152/1
def test394():
    C0 = HyperEllCurve.g2fromLMFDB("394.a.3152.1", inithomology = False)
    C = PointedHyperEllCurve(C0, [0, -1])
    print(C)
    C.init_homology([0, -0.5, -1])
    
    lvalue2 = mpf('0.85875464247178527993941189389477680680')
    lvalue1 = mpf('8.5704886281089963868641526976072444724') # lfuncheck = 1.9e-31
    
    C.test_beilinson(lvalue2, 20, s = 2)
    
    C.save("394.a.3152.1")

from regulator_integral_numpy import hyperell_integral, cplus, cminus

def testvolume():
    f = [-4,0, 0, 1, 2, 1]
    C = HyperEllCurve(f)
    print(C)
    C.plot_cycles()
    vol = C.jacobian_volume().real
    print("Cplus: ", cplus(f))
    print("Cminus: ", cminus(f))
    print("Volume: ", vol)
    
    print("lindep: ", )
    a, b = gp.lindep([vol, cplus(f)*cminus(f)], 8)
    print("lindep: ", a, b)
    print("evaluate: ", int(a)* vol + int(b)*cplus(f)*cminus(f) )


def testoldbeilinson():
    f = [-4, 0, 0, 1, 2, 1]
    C = PointedHyperEllCurve(f, [0, -1], inithomology=False)
    C.init_homology([0, -0.5, -1])
    C.plot_cycles()
    beilinson = C.beilinson_matrix()
    print(C)
    print("Beilinson matrix: ")
    print(beilinson)
    
    integrands = [lambda z : 1, lambda z : z.real, lambda z : np.abs(z)**2]
    old_beilinson = [[hyperell_integral(f, om)[0]/np.pi for om in integrands]]
    old_beilinson.append([hyperell_integral(f, lambda z : np.log(np.abs(z)) * om(z))[0]/np.pi for om in integrands])
    old_beilinson.append([hyperell_integral(f, lambda z : np.log(np.abs(z + 1)) * om(z))[0]/np.pi for om in integrands])
    print("Old Beilinson matrix: ")
    for row in old_beilinson:
        print(row)
        

#testoldbeilinson()

#test394()
#test3125()
#testvolume()


