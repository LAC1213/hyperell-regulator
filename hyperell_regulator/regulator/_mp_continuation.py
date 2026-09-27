"""Precision-controlled Newton continuation and nested quadrature.

Exact cover coefficients are converted directly to Sage complex fields.
Float64 supplies discrete cut topology, chart preflight, and initial root labels.
All accepted fibre points, differentials, quadrature sums and error tests
are evaluated at the requested working precision.
"""
import math
from sage.all import ComplexField, RealField, PolynomialRing, PowerSeriesRing, LaurentSeriesRing, QQ
from .cover import Cover
from .diagnostics import checkpoint, phase
from ._legs import LegFailure, QuadratureFailure


class Arithmetic:
    def __init__(self, digits):
        if int(digits) != digits or digits < 5:
            raise ValueError('prec must be an integer of at least 5 decimal digits')
        self.digits = int(digits)
        self.bits = int(math.ceil((self.digits+30)*math.log2(10)))
        self.C = ComplexField(self.bits)
        self.R = RealField(self.bits)
        self.output = RealField(int(math.ceil(self.digits*math.log2(10))))
        self.tol = self.R(10)**(-self.digits-3)
        self.newton_tol = self.R(10)**(-self.digits-20)
        self._rules = {}

    def rule(self, level):
        if level not in self._rules:
            R = self.R
            h = R(2)**(-level)
            lam = R.pi()/2
            reach = (self.digits+10)*R(10).log()/2
            count = int((reach/lam).arcsinh()/h)+1
            nodes = []
            for k in range(-count,count+1):
                t = h*k
                sh = lam*t.sinh()
                # Avoid subtracting two nearly equal numbers at the ends.
                q = (2*abs(sh)).exp()
                v = 1/(1+q)
                node = v if sh < 0 else 1-v
                weight = h*lam*t.cosh()*2*q/(1+q)**2
                if 0 < node < 1:
                    nodes.append((node,weight))
            self._rules[level] = nodes
        return self._rules[level]


class Chart:
    def __init__(self, cov, ar, kind='x', center=0, sigma=False, zr=0):
        if cov.exact is None:
            raise ValueError('arbitrary precision requires exact rational cover coefficients')
        self.ar,self.C,self.R = ar,ar.C,ar.R
        self.kind,self.center,self.sigma,self.zr = kind,self.C(center),sigma,self.C(zr)
        self.cov = cov
        if kind=='x':
            F,a,ad,b,bd=cov.exact
            shifted=F.parent().gen()+QQ(center)
            self.local=Cover(F(shifted),a(shifted)/ad(shifted),b(shifted)/bd(shifted))
        else:
            F,a,ad,b,bd=cov.exact
            u=F.parent().fraction_field().gen()
            inverse=QQ(center)+1/u
            self.local=Cover(F.parent()(u**(2*cov.g+2)*F(inverse)),
                a(inverse)/ad(inverse),b(inverse)/bd(inverse)/u**(cov.g+1))
        if self.local.exact is None:
            raise ValueError('continuation chart centers must be rational and real')
        F,a,ad,b,bd = self.local.exact
        S = PolynomialRing(self.C,'u')
        self.F,self.a,self.ad,self.b,self.bd = [S(q) for q in (F,a,ad,b,bd)]
        self.Fp = self.F.derivative()
        self.da,self.dad,self.db,self.dbd = [q.derivative() for q in (self.a,self.ad,self.b,self.bd)]
        # Use the exact, cancelled fibre polynomial, never float64 arrays.
        ps = [b*b*ad*ad*F-bd*bd*a*a,2*bd*bd*ad*a,-(bd*ad)**2]
        common = ps[0].gcd(ps[1]).gcd(ps[2])
        self.p = [S(q//common) for q in ps]
        self.roots = [self.C(r) for r in F.roots(self.C,multiplicities=False)]
        for r in F.roots(QQ,multiplicities=False):
            index=min(range(len(self.roots)),key=lambda i:abs(self.roots[i]-self.C(r)))
            self.roots[index]=self.C(r)
        self.norm_factors = self.local.norm_factors_hp(ar.bits+16)
        self.norm_factors = [[(self.C(r),int(m)) for r,m in part] for part in self.norm_factors]
        self.norm_lead = self.C(self.local.norm_leading())
        self.g,self.N,self.pure = cov.g,cov.N,cov.pure_x

    def norm(self,u):
        value = self.norm_lead
        for sign,part in ((1,self.norm_factors[0]),(-1,self.norm_factors[1])):
            for r,m in part:
                value *= (u-r)**(sign*m)
        return value

    def dlognorm(self,u):
        return sum(sign*m/(u-r) for sign,part in ((1,self.norm_factors[0]),(-1,self.norm_factors[1])) for r,m in part)

    def values(self,u,eta):
        a,ad,b,bd = self.a(u),self.ad(u),self.b(u),self.bd(u)
        if (self.F[0] and (not self.ad[0] or not self.bd[0])
                and abs(u)<self.R('1e-6')*min([self.R(1)]+[abs(r) for r in self.roots if r])):
            # Extend a cancelling branch analytically through infinity.
            # Near u=0, subtracting logarithmic derivatives of poles also
            # loses precision, so use the same local series in a small disk.
            key=1 if abs(eta-self.F[0].sqrt())<=abs(eta+self.F[0].sqrt()) else -1
            if not hasattr(self,'_regular_zero'):self._regular_zero={}
            if key not in self._regular_zero:
                depth=self.ar.digits+max(self.ad.degree(),self.bd.degree())+20
                exactF,a0,ad0,b0,bd0=self.local.exact
                base=QQ if exactF[0].is_square() else self.C
                S=LaurentSeriesRing(base,'s',default_prec=depth)
                root=S(PowerSeriesRing(base,'s',default_prec=depth)(exactF).sqrt())*key
                value=S(a0)/S(ad0)+S(b0)/S(bd0)*root
                self._regular_zero[key]=(None if value.valuation()<0 else
                    self.F.parent()([self.C(value[k]) for k in range(value.prec())]))
            regular=self._regular_zero[key]
            if regular is not None:return regular(u),eta*regular.derivative()(u)
            if not u:raise LegFailure('coordinate zero is a pole on this sheet')
        A,B = a/ad,b/bd
        da = eta*(self.da(u)*ad-a*self.dad(u))/ad**2
        db = eta**2*(self.db(u)*bd-b*self.dbd(u))/bd**2
        bf = B*self.Fp(u)/2
        small,big = A+B*eta,A-B*eta
        den = da+db+bf
        if not self.pure and big and abs(small) <= abs(big):
            value = self.norm(u)/big
            den = value*(eta*self.dlognorm(u)-(da-db-bf)/big)
        else:
            value = small
        return value,den

    def polynomial(self,z):
        if self.pure:
            return self.a*z-self.ad*(1+self.zr*z) if self.sigma else self.a-self.ad*z
        if self.sigma:
            v = 1+self.zr*z
            return self.p[0]*z*z+self.p[1]*z*v+self.p[2]*v*v
        return self.p[0]+self.p[1]*z+self.p[2]*z*z

    def eta(self,u,near,z):
        r = self.F(u).sqrt()
        if self.pure:
            return r if abs(r-near)<=abs(r+near) else -r
        if not self.ad(u) or not self.bd(u):
            return r if abs(r-near)<=abs(r+near) else -r
        A,B = self.a(u)/self.ad(u),self.b(u)/self.bd(u)
        lhs = (1+self.zr*z)-A*z if self.sigma else z-A
        rhs = B*r*(z if self.sigma else 1)
        return r if abs(lhs-rhs)<abs(lhs+rhs) else -r

    def solve(self,z,u,eta,poly=None,derivative=None):
        p = self.polynomial(z) if poly is None else poly
        dp = p.derivative() if derivative is None else derivative
        # The polynomial identifies the sheet; the factored norm polishes it.
        for k in range(10):
            den = dp(u)
            if not den: break
            step = p(u)/den
            u -= step
            if abs(step) < self.R('1e-8')*max(abs(u),self.R('1e-100000')): break
        eta = self.eta(u,eta,z)
        last = self.R(float('inf'))
        for k in range(16):
            value,den = self.values(u,eta)
            if not den: raise LegFailure('zero Newton derivative')
            if self.sigma:
                v = value-self.zr
                step = -eta*(1/v-z)*v*v/den
            else:
                step = eta*(value-z)/den
            u -= step
            r = self.F(u).sqrt()
            eta = r if abs(r-eta)<=abs(r+eta) else -r
            error = abs(step)
            if error <= self.ar.newton_tol*max(abs(u),self.R('1e-100000')):
                return (*self.refine_point(z,u,eta),error,k+1)
            if error >= last and error <= self.ar.tol*self.R('1e-8')*max(1,abs(u)):
                return (*self.refine_point(z,u,eta),error,k+1)
            last = error
        raise LegFailure('Newton did not attain the working precision')

    def refine_point(self,z,u,eta):
        # Near an unramified Weierstrass point, sqrt(F(u)) loses half
        # the coordinate digits. The fibre equation is linear in eta.
        if not self.pure and abs(eta)<self.R('1e-4') and self.ad(u) and self.bd(u):
            original=u
            p=self.polynomial(z);dp=p.derivative()
            for _ in range(4):
                if not dp(u):break
                correction=p(u)/dp(u)
                u-=correction
                if abs(correction)<=self.ar.newton_tol*max(1,abs(u)):break
            A,B=self.a(u)/self.ad(u),self.b(u)/self.bd(u)
            if B:
                w=self.zr+1/z if self.sigma else z
                refined=(w-A)/B
                if abs(refined**2-self.F(u))<=self.ar.newton_tol*max(1,abs(self.F(u))):
                    return u,refined
            return original,eta
        return u,eta

    def to_local(self,x,y):
        u = x-self.center if self.kind=='x' else 1/(x-self.center)
        return u,y if self.kind=='x' else y*u**(self.g+1)

    def to_global(self,u,eta):
        return (u+self.center,eta) if self.kind=='x' else (self.center+1/u,eta/u**(self.g+1))

    def differential(self,u,eta,z,rate):
        value,den = self.values(u,eta)
        if self.sigma: den *= -z*z
        if self.kind=='x':
            x = u+self.center
            return [x**a*rate/den for a in range(self.g)]
        t = 1+self.center*u
        return [-t**a*u**(self.g-1-a)*rate/den for a in range(self.g)]


class Walker:
    """A fibre with adaptive Newton steps, root separation and y-sheet checks."""
    def __init__(self,chart,coord,points):
        self.ch,self.coord,self.ar = chart,coord,chart.ar
        self.points = [chart.to_local(*p) for p in points]
        self.t,self.h = self.ar.R(0),self.ar.R(1)/16
        self.points = [chart.solve(coord(self.t),u,e)[:2] for u,e in self.points]
        # For pure x maps two curve points deliberately share each x root.
        self.partners = []
        for i,(u,e) in enumerate(self.points):
            self.partners.append({j for j,(v,f) in enumerate(self.points)
                                  if i==j or (chart.pure and abs(u-v)<self.ar.tol*max(1,abs(u)))})

    def to(self,end):
        R = self.ar.R
        while self.t < end:
            h = min(self.h,end-self.t)
            cut = False
            while True:
                checkpoint('walk_attempts',t=float(self.t),h=float(h),precision=self.ar.digits)
                z = self.coord(self.t+h)
                good = True
                try:
                    p = self.ch.polynomial(z)
                    dp = p.derivative()
                    new = [self.ch.solve(z,u,e,p,dp) for u,e in self.points]
                    for i,((u,e),(v,f,error,_)) in enumerate(zip(self.points,new)):
                        others = [abs(v-q[0]) for j,q in enumerate(new) if j not in self.partners[i]]
                        if others:
                            gap=min(others)
                        else:
                            # Taylor coefficients give a root-separation
                            # radius even when an ordinary lift crosses u=0.
                            # A bound proportional to |u| falsely stalls there.
                            shifted=p(p.parent().gen()+v)
                            linear=abs(shifted[1])
                            radii=[(linear/(2*abs(shifted[k])))**(R(1)/(k-1))/2
                                   for k in range(2,shifted.degree()+1) if shifted[k]]
                            gap=min(radii) if radii else R(1)+abs(v)
                        if abs(v-u)>gap/4 or error>gap/100:
                            good=False;break
                        if abs(f-e)>abs(f+e):
                            distances = [abs(r-u) for r in self.ch.roots]
                            movement = abs(v-u)
                            if distances and movement<min(distances):
                                if sum((movement/d).arcsin() for d in distances)<R.pi():
                                    good=False;break
                except (LegFailure,ZeroDivisionError,ValueError,OverflowError):
                    good = False
                if good: break
                checkpoint('walk_rejected')
                h /= 2;cut=True
                if h < R('1e-14')*(end-self.t) or self.t+h==self.t:
                    raise LegFailure('arbitrary-precision continuation stalled at %s' % z)
            self.points = [(q[0],q[1]) for q in new]
            self.t += h
            checkpoint('walk_accepted')
            if cut: self.h=h
            elif h>=self.h: self.h=h*R('1.5')
        return self.points

    def global_points(self):
        return [self.ch.to_global(*p) for p in self.points]


def quadrature(ar,evaluate,dimension,breaks=None,level=3,max_level=20):
    """Nested tanh-sinh with cached ordered nodes and an explicit error test."""
    C,R = ar.C,ar.R
    cuts = [R(0),R(1)] if breaks is None else sorted(set(map(R,breaks)))
    cache,previous = {},None
    last = None
    for lev in range(level,max_level+1):
        with phase('quadrature_level',level=lev,precision=ar.digits):
            samples = []
            for a,b in zip(cuts[:-1],cuts[1:]):
                for t,w in ar.rule(lev):
                    node = a+(b-a)*t
                    if a<node<b: samples.append((node,(b-a)*w))
            evaluate.reset()
            total = [C(0) for _ in range(dimension)]
            for t,w in samples:
                if t in cache:
                    value,state = cache[t]
                    evaluate.restore(t,state)
                else:
                    checkpoint('quadrature_nodes')
                    value,state = evaluate(t)
                    cache[t] = value,state
                for i,v in enumerate(value): total[i] += w*v
            if previous is not None:
                last = max(abs(a-b)/max(1,abs(a)) for a,b in zip(total,previous))
                checkpoint(scaled_quadrature_error=str(last),level=lev)
                if last <= ar.tol:
                    return total
            previous = total
    raise QuadratureFailure('MP quadrature did not converge by level %d; error=%s, tolerance=%s' % (max_level,last,ar.tol))


class FunctionSamples:
    def __init__(self,function): self.function=function
    def reset(self): pass
    def restore(self,t,state): pass
    def __call__(self,t): return self.function(t),None


class WalkSamples:
    def __init__(self,ch,coord,rate,points,lam,wof=None,logarithmic=False):
        self.ch,self.coord,self.rate,self.initial,self.lam=ch,coord,rate,points,lam
        self.wof=coord if wof is None else wof
        self.logarithmic=logarithmic
        self.w0=self.wof(self.ch.ar.R(0))
    def reset(self):
        self.walk=Walker(self.ch,self.coord,self.initial)
        self.prev=self.ch.ar.R(0);self.log=self.lam
    def restore(self,t,state):
        parameter,points,lam=state
        self.walk.t,self.walk.points=parameter,list(points)
        self.prev,self.log=parameter,lam
    def __call__(self,t):
        parameter=-(1-t).log() if self.logarithmic else t
        points=self.walk.to(parameter)
        self.log += (self.wof(parameter)/self.wof(self.prev)).log()
        self.prev=parameter
        scale=1/(1-t) if self.logarithmic else 1
        vals=[self.ch.differential(u,e,self.coord(parameter),self.rate(parameter)*scale) for u,e in points]
        flat=[v for row in vals for v in row]
        return flat+[self.log*v for v in flat],(parameter,list(points),self.log)
