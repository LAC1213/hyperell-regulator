"""Arbitrary-precision implementation of the float64 cut-and-lift algorithm."""
import numpy as np
from sage.all import PolynomialRing, QQ
from .cover import Cover
from . import _legs
from .diagnostics import checkpoint, phase, _active
from .stokes import cut_data, cut_scale, hold_radii, boundary_walks
from ._mp_continuation import Arithmetic, Chart, Walker, WalkSamples, FunctionSamples, quadrature


class Integrator:
    def __init__(self,cov,ar):
        self.cov,self.ar,self.C,self.R = cov,ar,ar.C,ar.R
        zeros=cov.norm_exact()[0].roots(QQ)
        self.center=zeros[0][0] if len(zeros)==1 else QQ(0)
        self.chart=Chart(cov,ar,center=self.center)
        self.global_chart=self.chart if self.center==0 else Chart(cov,ar)
        self.charts={}
        self.g,self.n=cov.g,cov.N
        self.zero=self.C(0)
        F,a,ad,b,bd=cov.exact
        A,B=a/ad,b/bd
        q=(A.derivative().numerator() if cov.pure_x else
           (F*A.derivative()**2-(F*B.derivative()+F.derivative()*B/2)**2).numerator())
        self.critical_points=[] if not q else [self.C(r) for r in (q//q.gcd(q.derivative())).roots(self.C,multiplicities=False)]
        self.critical_values=[self.C(v) for v in cov.critical_values_hp(ar.bits+32)]

    def zeros(self): return [[self.zero for _ in range(self.n)] for _ in range(self.g)]

    def fibre(self,w,reference=None):
        x,y=self.cov.fibre(complex(w) if reference is None else reference)
        points=[(self.C(a),self.C(b)) for a,b in zip(x,y)]
        walker=Walker(self.chart,lambda t:w,points)
        return walker.global_points()

    def unpack(self,values,n=None):
        n=self.n if n is None else n
        k=n*self.g
        return ([[values[i*self.g+a] for i in range(n)] for a in range(self.g)],
                [[values[k+i*self.g+a] for i in range(n)] for a in range(self.g)])

    def walk(self,ch,coord,rate,points,lam,wof=None,logarithmic=False,breaks=None):
        samples=WalkSamples(ch,coord,rate,points,lam,wof,logarithmic)
        values=quadrature(self.ar,samples,2*self.g*len(points),breaks=breaks)
        I,L=self.unpack(values,len(points))
        if logarithmic:
            return None,None,I,L
        samples.walk.to(self.R(1))
        endpoint=samples.walk.global_points()
        wof=coord if wof is None else wof
        endlam=samples.log+(wof(1)/wof(samples.prev)).log()
        return endpoint,endlam,I,L

    def finite(self,a,b,points,lam,toward=None,avoid=True):
        path=_legs._pinch_detour(self.cov,complex(a),complex(b)) if avoid else None
        if path is not None:
            checkpoint('pinch_detours',segments=len(path)-1)
            path=[a]+[self.C(z) for z in path[1:-1]]+[b]
            I,L=self.zeros(),self.zeros()
            for left,right in zip(path[:-1],path[1:]):
                points,lam,i,l=self.finite(left,right,points,lam,toward,False)
                I=self.add(I,i);L=self.add(L,l)
            return points,lam,I,L
        delta=b-a
        coord=lambda t:(1-t)*a+t*b
        rate=lambda t:delta
        cuts=None
        if toward is not None and toward[1]>1:
            z,e=toward
            ratio=(b-z)/(a-z)
            if abs(ratio.imag()) <= self.R('1e-9')*abs(ratio) and ratio.real()>0:
                ratio=ratio.real()
                sa=(a-z)**(self.R(1)/e)
                sb=sa*ratio**(self.R(1)/e)
                parameter=lambda t:(1-t)*sa+t*sb
                coord=lambda t:z+parameter(t)**e
                rate=lambda t:e*parameter(t)**(e-1)*(sb-sa)
                if ratio<self.R('1e-6'):
                    lr=ratio.log()
                    radial=lambda t:(a-z)*(lr*t).exp()
                    coord=lambda t:z+radial(t)
                    rate=lambda t:lr*radial(t)
                    projected=[((v-z)/(a-z)).real() for v in self.critical_values+[self.zero]]
                    cuts=sorted(set([self.R(0),self.R(1)]+[r.log()/lr for r in projected if ratio<r<1]))
        # The x coordinate is singular at a finite-value point at infinity.
        # Reuse the float64 preflight's chart choices, then integrate each
        # chart group entirely at the requested precision.
        xx=np.array([complex(p[0]) for p in points]);yy=np.array([complex(p[1]) for p in points])
        _,_,_,prepared=_legs.leg_in_w(self.cov,complex(a),complex(b),xx,yy,
            np.zeros(len(points),complex),toward=None if toward is None else (complex(toward[0]),toward[1]),
            integrate=False,_avoid_pinches=False)
        groups={}
        for j,ch in enumerate(prepared.charts):
            key=(ch.kind,QQ(ch.c.real)) if ch.kind=='inv' else ('x',self.center)
            groups.setdefault(key,[]).append(j)
        ends=[None]*len(points);I,L=self.zeros(),self.zeros()
        for key,indices in groups.items():
            if key not in self.charts:self.charts[key]=Chart(self.cov,self.ar,key[0],key[1])
            tail,ll,ii,logs=self.walk(self.charts[key],coord,rate,[points[j] for j in indices],lam,breaks=cuts)
            for k,j in enumerate(indices):
                ends[j]=tail[k]
                for q in range(self.g):I[q][j]=ii[q][k];L[q][j]=logs[q][k]
        return ends,ll,I,L

    @staticmethod
    def add(A,B,sign=1): return [[a+sign*b for a,b in zip(x,y)] for x,y in zip(A,B)]

    def target(self,seed,z,kind):
        x,y=seed
        chart=_legs.chart_at(self.cov,x,y)
        M=self.global_chart
        if chart=='tau': return None,None,'tau'
        if chart=='y':
            xt=min(M.roots,key=lambda r:abs(r-self.C(x)))
            return xt,self.zero,'y'
        if kind in ('P','Q'):
            roots=M.norm_factors[0 if kind=='P' else 1]
            xt=min(roots,key=lambda q:abs(q[0]-self.C(x)))[0]
        else:
            candidates=sorted(self.critical_points,key=lambda r:abs(r-self.C(x)))
            xt=None
            if candidates and abs(candidates[0]-self.C(x))<self.R('1e-5')*max(1,abs(x)):
                r=candidates[0];v=M.F(r).sqrt()
                v=v if abs(v-self.C(y))<abs(v+self.C(y)) else -v
                if abs(M.values(r,v)[0]-z)<self.R('1e-8')*max(1,abs(z)):
                    xt=r
            if xt is None:
                u,e=self.chart.to_local(self.C(x),self.C(y))
                u,e,_,_=self.chart.solve(z,u,e)
                return (*self.chart.to_global(u,e),'x')
        yt=M.F(xt).sqrt()
        yt=yt if abs(yt-self.C(y))<=abs(yt+self.C(y)) else -yt
        return xt,yt,'x'

    def endpoint(self,point,lam,xt,yt,chart,kind):
        """The same x/y/tau uniformizers as the float64 endpoint integrator."""
        x0,y0=point
        if chart=='tau': return self.tau(point,lam)
        if kind=='P' and chart=='x':
            # Follow phi=w0*exp(-N*r), the float64 P-uniformizer path.
            u,e=self.chart.to_local(x0,y0)
            w0=self.chart.values(u,e)[0]
            coord=lambda r:w0*(-self.R(self.n)*r).exp()
            rate=lambda r:-self.n*coord(r)
            _,_,I,L=self.walk(self.chart,coord,rate,[point],lam,logarithmic=True)
            return [row[0] for row in I]+[row[0] for row in L]
        M=self.global_chart
        S=M.F.parent();u=S.gen()
        F=M.F(u+xt)
        if chart=='y':
            F-=F[0]
        Fp=F.derivative()
        a,ad,b,bd=[q(u+xt) for q in (M.a,M.ad,M.b,M.bd)]
        factors=[]
        for sign,part in ((1,M.norm_factors[0]),(-1,M.norm_factors[1])):
            for root,m in part:
                shift=xt-root
                if abs(shift)<self.ar.newton_tol*max(1,abs(xt)): shift=self.zero
                factors.append((shift,sign*m))
        def phi(v,y):
            A,B=a(v)/ad(v),b(v)/bd(v)
            small,big=A+B*y,A-B*y
            if abs(small)<=abs(big) and big:
                norm=M.norm_lead
                for shift,m in factors: norm*=(v+shift)**m
                return norm/big
            return small
        u0=x0-xt
        phi0=phi(u0,y0)
        if chart=='y':
            def evaluate(t):
                y=y0*(1-t)
                v=y*y/Fp(0)
                for k in range(32):
                    checkpoint('endpoint_iterations',chart='y',precision=self.ar.digits)
                    step=(F(v)-y*y)/Fp(v)
                    v-=step
                    if abs(step)<=self.ar.newton_tol*max(abs(v),self.R('1e-100000')): break
                else: raise _legs.LegFailure('MP Weierstrass endpoint Newton failed')
                ll=lam+(phi(v,y)/phi0).log()
                values=[-2*(xt+v)**j*y0/Fp(v) for j in range(self.g)]
                return values+[ll*r for r in values]
        else:
            F0=F(u0)
            def evaluate(t):
                v=u0*(1-t)
                y=y0
                for root in M.roots:y*=((xt+v-root)/(x0-root)).sqrt()
                ll=lam+(phi(v,y)/phi0).log()
                values=[-(xt+v)**j*u0/y for j in range(self.g)]
                return values+[ll*r for r in values]
        return quadrature(self.ar,FunctionSamples(evaluate),2*self.g)

    def tau(self,point,lam):
        M=self.global_chart;x0,y0=point
        even=self.cov.even;g=self.g
        tau0=1/x0 if even else 1/x0.sqrt()
        eta0=y0*tau0**(g+1 if even else 2*g+1)
        sq=1 if even else 2
        S=M.F.parent();t=S.gen()
        rev=lambda q:S(list(q)[::-1])
        F=rev(M.F)
        arg=lambda tau:tau**sq
        F0=F(arg(tau0))
        ra,rad,rb,rbd=[rev(q) for q in (M.a,M.ad,M.b,M.bd)]
        pA=sq*(M.a.degree()-M.ad.degree())
        pB=None if not M.b else sq*(M.b.degree()-M.bd.degree())+(g+1 if even else 2*g+1)
        pw=max(pA,pB) if pB is not None else pA
        def value(tau,eta):
            A=tau**(pw-pA)*ra(arg(tau))/rad(arg(tau))
            B=0 if pB is None else tau**(pw-pB)*rb(arg(tau))/rbd(arg(tau))*eta
            return A+B,A-B
        phi0,big0=value(tau0,eta0)
        nn,nd=self.cov.norm_exact()
        rn,rd=rev(S(nn)),rev(S(nd))
        norm_power=sq*(nd.degree()-nn.degree())
        rn0=rn(arg(tau0))/rd(arg(tau0))
        use_norm=abs(phi0)<abs(big0)
        def evaluate(t):
            tau=tau0*(1-t)
            eta=eta0*(F(arg(tau))/F0).sqrt()
            small,big=value(tau,eta)
            if use_norm:
                ratio=rn(arg(tau))/rd(arg(tau))
                ll=lam+(norm_power+pw)*(1-t).log()+(ratio/rn0).log()+(big0/big).log()
            else:
                ll=lam-pw*(1-t).log()+(small/phi0).log()
            values=[(tau**(g-1-j) if even else 2*tau**(2*g-2*j-2))*tau0/eta for j in range(g)]
            return values+[ll*r for r in values]
        return quadrature(self.ar,FunctionSamples(evaluate),2*g)

    def side(self,m,z,kind,points,lam,hold=None):
        if kind=='none': return self.finite(m,z,points,lam)[2:]
        frac=min(self.R('.08'),self.R(hold)/abs(m-z)) if hold is not None else self.R('.08')
        startx=np.array([complex(p[0]) for p in points]);starty=np.array([complex(p[1]) for p in points])
        e=self.cov.endpoint_ramification(complex(z),kind)
        last=None
        for attempt in range(16):
            checkpoint('standoff_attempts',fraction=str(frac))
            w=z+frac*(m-z)
            try:
                with phase('endpoint_preflight'):
                    xx,yy,_,_= _legs.leg_in_w(self.cov,complex(m),complex(w),startx,starty,np.zeros(self.n,complex),toward=(complex(z),e),integrate=False)
                    fibre=self.cov.fibre(complex(z))
                    far=int(np.count_nonzero(~np.isfinite(fibre[0])))
                    gone=set(np.argsort(-abs(xx))[:far].tolist()) if far else set()
                    seeds=[(np.inf,np.inf) if j in gone else _legs._target(self.cov,a,b,complex(z),kind,fib=fibre)
                           for j,(a,b) in enumerate(zip(xx,yy))]
                    targets=[self.target(seed,z,kind) for seed in seeds]
                    if not(kind=='P' and all(t[2]=='x' for t in targets)):
                        for j,(a,b) in enumerate(zip(xx,yy)):
                            _legs.end_leg(self.cov,complex(w),a,b,0.,complex(z),kind,False,
                                          at_infinity=j in gone,fib=fibre)
                ends,ll,I,L=self.finite(m,w,points,lam,(z,e))
                if kind=='P' and all(target[2]=='x' for target in targets):
                    coord=lambda r:w*(-self.R(self.n)*r).exp()
                    rate=lambda r:-self.n*coord(r)
                    _,_,tailI,tailL=self.walk(self.chart,coord,rate,ends,ll,logarithmic=True)
                    return self.add(I,tailI),self.add(L,tailL)
                for j,(point,target) in enumerate(zip(ends,targets)):
                    tail=self.endpoint(point,ll,*target,kind)
                    for a in range(self.g): I[a][j]+=tail[a];L[a][j]+=tail[self.g+a]
                return I,L
            except _legs.LegFailure as exc:
                last=exc;frac/=3
        raise _legs.LegFailure('no MP endpoint handover: %s'%last)

    def ray(self,z,d,span,points,lam):
        cov=self.cov
        poles=cov.norm_exact()[1].roots(QQ)
        if poles:
            kind,center='x',poles[0][0]
        else:
            kind,center='inv',0
        ch=Chart(cov,self.ar,kind,center,sigma=True,zr=z)
        sigma0=1/(span*d)
        if not cov.even and len(cov.ad)==len(cov.bd)==1:
            radius,bound=_legs._odd_infinity_region(cov)
            factor=max(self.R(1),(2*self.R(bound)+2*abs(z)+span*abs(d))/(span*abs(d)))
            lr=factor.log()
            coord=lambda t:sigma0*(-lr*t).exp()
            rate=lambda t:-lr*coord(t)
            wof=lambda t:z+1/coord(t)
            ends,ll,I,L=self.walk(ch,coord,rate,points,lam,wof)
            if any(abs(x)<=radius for x,y in ends): raise _legs.LegFailure('MP ray outside infinity chart')
            for j,point in enumerate(ends):
                tail=self.tau(point,ll)
                for a in range(self.g):I[a][j]+=tail[a];L[a][j]+=tail[self.g+a]
            return I,L
        pls=cov.places_over(None)
        e=max(int(q['e']) for q in pls) if pls else cov.N
        coord=lambda r:sigma0*(-self.R(e)*r).exp()
        rate=lambda r:-e*coord(r)
        wof=lambda r:z+1/coord(r)
        return self.walk(ch,coord,rate,points,lam,wof,logarithmic=True)[2:]


def stokes_sum(ar,W,u,ell,dl,tolerance):
    C=ar.C;g=len(u);r=len(u[0]);n=len(u[0][0])
    H=[[C(0) for _ in range(g)] for _ in range(g)]
    forward,backward={},{}
    for k,walk in enumerate(W):
        for j,(i,s,sg) in enumerate(walk):(forward if sg>0 else backward)[i,s]=(k,j)
    for a in range(g):
        plus=[[ell[a][i][s]+dl*u[a][i][s] for s in range(n)] for i in range(r)]
        U=[]
        for walk in W:
            vals=[C(0)]
            for i,s,sg in walk:vals.append(vals[-1]+(plus[i][s] if sg>0 else -ell[a][i][s]))
            U.append(vals)
        scale=max([ar.R(1)]+[abs(v) for row in plus for v in row])
        defect=max(abs(row[-1]) for row in U)/scale
        checkpoint(boundary_defect=str(defect),precision=ar.digits)
        if defect>tolerance:raise RuntimeError('MP sheet boundaries do not close: %s (tolerance %s)'%(defect,tolerance))
        for i in range(r):
            for s in range(n):
                k,j=forward[i,s];l,jp=backward[i,s]
                c=U[k][j]-U[l][jp]+ell[a][i][s]
                for b in range(g):H[a][b]+=c*u[b][i][s].conjugate()
    return H


def compute(F,phi,prec,check_tol=1e-6):
    ar=Arithmetic(prec);C,R=ar.C,ar.R
    cov=Cover(F,*phi)
    engine=Integrator(cov,ar)
    z,kinds,d=cov.cut_path()
    upper,lower,lam_minus,dlam=cut_data(cov,z,d)
    W=boundary_walks(upper,lower)
    zh=[C(t) if k=='none' else min(engine.critical_values,key=lambda v:abs(v-C(t))) for t,k in zip(z,kinds)]
    dh=C(d);span=R(cut_scale(z)[0]);holds=hold_radii(z)
    g,n,r=cov.g,cov.N,len(z)
    u=[[] for _ in range(g)];ell=[[] for _ in range(g)]
    tolerance=min(R(check_tol),R(10)**(-prec+3))
    from contextlib import nullcontext
    for i in range(r):
        mon=_active.get()
        context=mon.edge('finite_mp' if i<r-1 else 'ray_mp',degree=n,precision=prec) if mon else nullcontext()
        with context:
            m=(zh[i]+zh[i+1])/2 if i<r-1 else zh[-1]+span*dh
            reference=(z[i]+z[i+1])/2 if i<r-1 else z[-1]+cut_scale(z)[0]*d
            points=engine.fibre(m,reference)
            turns=round((complex(lam_minus[i]).imag-np.angle(complex(m)))/(2*np.pi))
            lam=m.log()+C.gen()*2*R.pi()*turns
            Ia,La=engine.side(m,zh[i],kinds[i],points,lam,holds[i])
            if i<r-1:Ib,Lb=engine.side(m,zh[i+1],kinds[i+1],points,lam,holds[i+1])
            else:Ib,Lb=engine.ray(zh[-1],dh,span,points,lam)
            I,L=engine.add(Ib,Ia,-1),engine.add(Lb,La,-1)
            for vals in (I,L):
                scale=max([R(1)]+[abs(v) for row in vals for v in row])
                defect=max(abs(sum(row)) for row in vals)/scale
                if defect>tolerance:raise RuntimeError('MP edge trace does not vanish: %s'%defect)
            if mon:
                mon._current['trace_periods']=[str(sum(row)) for row in I]
                mon._current['trace_log_periods']=[str(sum(row)) for row in L]
            for a in range(g):u[a].append(I[a]);ell[a].append(L[a])
    dl=C.gen()*2*R.pi()*round(complex(dlam).imag/(2*np.pi))
    H0=stokes_sum(ar,W,u,u,C(0),tolerance)
    H=stokes_sum(ar,W,u,ell,dl,tolerance)
    per,reg=[],[]
    for a in range(g):
        for b in range(a,g):
            extra=dl.imag()/2*sum((u[a][i][s]*u[b][i][s].conjugate()).real() for i in range(r) for s in range(n))
            per.append(ar.output(-(H0[a][b].imag()+H0[b][a].imag())/8))
            reg.append(ar.output(-((H[a][b].imag()+H[b][a].imag())/2+extra)/4))
    return per,reg
