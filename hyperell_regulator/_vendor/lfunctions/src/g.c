#include <stdio.h> // must precede <gmp.h>: GMP only declares its FILE* I/O (mpz_inp_str) when <stdio.h> is seen first
#include <gmp.h>  // must precede FLINT headers so the GMP-interop decls (fmpz_set_mpz) are visible
#include "glfunc.h"
#include "glfunc_internals.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <flint/arb_poly.h>

#ifdef __cplusplus
extern "C"{
#endif

#define maxr MAX_DEGREE

// On-disk G-cache format.  The file opens with a self-describing header line
//   GCACHE <version> <degree> <gprec> <twomu_0> ... <twomu_{degree-1}>
// so a loaded file can be verified against the current request rather than
// trusted on the strength of its (lossy, mus-only) filename.  twomu_i is the
// i-th sorted analytic mu times two -- an exact integer for the half-integer
// mus the library enforces -- which also captures the normalisation (already
// folded into L->mus).  GCACHE_VERSION must be bumped whenever the on-disk
// layout or the printarb/printarf number format changes; a change to
// EXTRA_BITS or the gprec formula needs no bump, since the stored gprec is
// checked for sufficiency (a too-low cached gprec is rejected automatically).
#define GCACHE_MAGIC "GCACHE"
#define GCACHE_VERSION 1

// Outcome of validating a cached file against the current request.  GCACHE_STALE
// (foreign / old version / different L-function / insufficient precision) means
// recompute and overwrite the file; GCACHE_CORRUPT (header verified usable but
// the body then failed to parse) is a genuinely broken file and stays fatal.
typedef enum { GCACHE_OK, GCACHE_STALE, GCACHE_CORRUPT } gcache_status;


  static void printarf(FILE *fp,const arf_t x) {
    static int init;
    static fmpz_t m,e;

    if (!init) {
      fmpz_init(m); fmpz_init(e);
      init = 1;
    }
    arf_get_fmpz_2exp(m,e,x);
    fmpz_fprint(fp,m); fprintf(fp," "); fmpz_fprint(fp,e);
  }

  static void printarb(FILE *fp,const arb_t x) {
    static int init;
    static arf_t a;

    if (!init) {
      arf_init(a);
      init = 1;
    }

#if 1
    printarf(fp,arb_midref(x));
    fprintf(fp," ");
    arf_set_mag(a,arb_radref(x));
    printarf(fp,a);
#else
    arb_printd(x,15);
#endif
  }


  static void my_hurwitz_zeta(arb_t res,long s,const arb_t a,long prec) {
    static int init;
    static arb_t A,t;

    if (!init) {
      arb_init(A); arb_init(t);
      init = 1;
    }
    arb_set(A,a);
    arb_zero(res);
    while (!arb_is_positive(A)) {
      arb_pow_ui(t,A,s,prec);
      arb_inv(t,t,prec);
      arb_add(res,res,t,prec);
      arb_add_ui(A,A,1,prec);
    }
    arb_set_si(t,s);
    arb_hurwitz_zeta(t,t,A,prec);
    arb_add(res,res,t,prec);
  }

  // Maclaurin series of log(Gamma_R(s+n/2)/Gamma_R(n/2)),
  //   n positive or not a multiple of 4
  // log(Gamma_R(s+n/2)*s/2*(-n/4)!*(-Pi)^(n/4)),
  //   n a non-positive multiple of 4
  static void lgammaR(arb_poly_t res,long r,long n,long prec) {
    long i,m;
    static int init;
    static arb_t s,t;

    if (!init) {
      arb_init(s); arb_init(t);
      init = 1;
    }

    arb_poly_zero(res);
    arb_const_pi(s,prec);
    arb_log(s,s,prec);
    for (m=1;m<=r;m++) {
      if (!(n&3) && n <= 0) {
        if (m == 1) {
          arb_const_euler(t,prec);
          arb_add(s,s,t,prec);
        } else
          arb_zeta_ui(s,m,prec);
        for (i=n>>2;i<0;i++) {
          arb_si_pow_ui(t,i,m,prec);
          arb_inv(t,t,prec);
          arb_add(s,s,t,prec);
        }
      } else {
        arb_set_si(t,n);
        arb_mul_2exp_si(t,t,-2);
        if (m == 1) {
          arb_digamma(t,t,prec);
          arb_sub(s,s,t,prec);
        } else {
          arb_set_si(s,m);
          my_hurwitz_zeta(s,m,t,prec);
        }
      }
      arb_mul_2exp_si(s,s,-m);
      arb_div_si(s,s,(m&1)?-m:m,prec);
      arb_poly_set_coeff_arb(res,m,s);
    }
  }

  // Polar part of \exp((1/2-s)u)\prod_j\Gamma_\R(s+\mu_j) around s=-n/2
  // returned as a polynomial, with residue in highest degree term
  // returns the order of pole
  // caches each coeff as a poly in u for repeat evaluation
  static struct {
    arb_poly_t poly[maxr];
    long order;
  } *cache;
  static long mcache;
  static long polarpart(arb_poly_t res,long twomu[],long r,long n,
      arb_srcptr u,long prec) {
    long i,j,k,m;
    arb_ptr p;
    static int init;
    static arb_t c,t,pi;
    static arb_poly_t f,g;

    if (!init) {
      arb_init(c); arb_init(t); arb_init(pi);
      arb_poly_init(f); arb_poly_init(g);
      init = 1;
    }

    arb_poly_zero(res);
    if (n >= mcache) {
      k = 2*n;
      if (k < 1024) k = 1024;
      cache = realloc(cache,k*sizeof(*cache));
      while (mcache < k) {
        cache[mcache].order = -1;
        for (j=0;j<maxr;j++)
          arb_poly_init(cache[mcache].poly[j]);
        mcache++;
      }
    } else if (cache[n].order >= 0)
      goto computeres;

    m = 0;
    arb_one(c);
    arb_const_pi(pi,prec);
    for (j=0;j<r;j++)
      if ((k=n-twomu[j]) >= 0 && !(k&3)) {
        for (k>>=2;k>0;k--) {
          arb_mul(c,c,pi,prec);
          arb_div_si(c,c,-k,prec);
        }
        arb_mul_2exp_si(c,c,1);
        m++;
      } else {
        // c *= \gamma_R(-k/2)
        arb_set_si(t,-k);
        arb_mul_2exp_si(t,t,-2);
        arb_gamma(t,t,prec);
        arb_mul(c,c,t,prec);
        arb_set_si(t,k);
        arb_mul_2exp_si(t,t,-2);
        arb_pow(t,pi,t,prec);
        arb_mul(c,c,t,prec);
      }

    cache[n].order = m;
    if (!m) return 0;

    arb_poly_zero(f);
    for (j=0;j<r;j++) {
      lgammaR(g,m-1,twomu[j]-n,prec);
      arb_poly_add(f,f,g,prec);
    }
    arb_poly_exp_series(g,f,m,prec);
    for (j=0;j<m;j++)
      if ((p = arb_poly_get_coeff_ptr(g,j)))
        arb_mul(p,p,c,prec);

    for (i=0;i<m;i++) {
      arb_poly_zero(cache[n].poly[i]);
      for (j=0;j<m-i;j++) {
        arb_poly_get_coeff_arb(t,g,m-(i+1)-j);
        for (k=j;k>0;k--)
          arb_div_si(t,t,-k,prec);
        arb_poly_set_coeff_arb(cache[n].poly[i],j,t);
      }
    }

computeres:
    // scale by exp((n+1)*u/2)
    arb_mul_si(t,u,n+1,prec);
    arb_mul_2exp_si(t,t,-1);
    arb_exp(c,t,prec);
    for (i=0;i<cache[n].order;i++) {
      arb_poly_evaluate(t,cache[n].poly[i],u,prec);
      arb_mul(t,t,c,prec);
      arb_poly_set_coeff_arb(res,cache[n].order-1-i,t);
    }
    return cache[n].order;
  }

  static void myabs(arb_t y,const arb_t x,long prec) {
    static int init;
    static arf_t l,r;

    if (!init) {
      arf_init(l); arf_init(r);
      init = 1;
    }
    arb_abs(y,x);
    arb_get_interval_arf(l,r,y,prec);
    if (arf_sgn(l) < 0) {
      arf_zero(l);
      arb_set_interval_arf(y,l,r,prec);
    }
  }

  // compute, for l=1,...,r,
  // Res_{s=-n/2}\exp((1/2-s)u)\prod_{j=1}^r\Gamma_\R(s+\mu_j)
  //             * \prod_{j=1}^l(-s-\mu_j)
  // and maximum laurent coeff for each l
  static void residues(arb_t res[],arb_t maxc[],long twomu[],long r,
      long n,arb_srcptr u,long prec) {
    long i,j,m;
    static int init;
    static arb_t t;
    static arb_poly_t f,g;

    if (!init) {
      arb_init(t);
      arb_poly_init(f);
      arb_poly_init(g);
      init = 1;
    }

    m = polarpart(g,twomu,r,n,u,prec);
    if (!m) {
      for (j=0;j<r;j++) arb_zero(res[j]);
      return;
    }

    arb_poly_one(f); arb_poly_neg(f,f);
    arb_poly_shift_left(f,f,1);
    for (j=0;j<r;j++) {
      arb_poly_get_coeff_arb(res[j],g,m-1);
      arb_set(maxc[j],res[j]);
      for (i=m-2;i>=0;i--) {
        arb_poly_get_coeff_arb(t,g,i);
        arb_union(maxc[j],maxc[j],t,prec);
      }
      myabs(maxc[j],maxc[j],prec);
      arb_set_si(t,n-twomu[j]);
      arb_mul_2exp_si(t,t,-1);
      arb_poly_set_coeff_arb(f,0,t);
      arb_poly_mullow(g,g,f,m,prec);
    }
  }

  // Taylor polynomial of G around u
  static void gtaylor(arb_t res[],long twomu[],long r,
      arb_srcptr u,long k,long prec,long prec2) {
    long i,j,n;
    static int init;
    static arb_t c,t,thresh,temp[maxr],maxc[maxr];
    static arb_poly_t f,g;

    if (!init) {
      for (j=0;j<maxr;j++) {
        arb_init(temp[j]);
        arb_init(maxc[j]);
      }
      arb_init(c); arb_init(t); arb_init(thresh);
      arb_poly_init(f);
      arb_poly_init(g);
      init = 1;
    }

    // c = (2*Pi)^r*exp(2*u)
    arb_const_pi(t,prec2);
    arb_mul_2exp_si(t,t,1);
    arb_pow_ui(c,t,r,prec2);
    arb_mul_2exp_si(t,u,1);
    arb_exp(t,t,prec2);
    arb_mul(c,c,t,prec2);

    // thresh = 2^{-prec}
    arb_one(thresh);
    arb_mul_2exp_si(thresh,thresh,-prec-EXTRA_BITS);

    for (j=0;j<r;j++)
      arb_zero(res[j]);
    for (n=0;;n++) {
      residues(temp,maxc,twomu,r,n,u,prec2);
      for (j=0;j<r;j++)
        arb_add(res[j],res[j],temp[j],prec2);
      arb_one(t); arb_mul_2exp_si(t,t,-1);
      for (j=0;j<r;j++) {
        if (n < twomu[j]+4 || !arb_lt(maxc[j],thresh)) break;
        arb_mul_ui(t,t,n-twomu[j]-2,prec2);
        arb_mul_2exp_si(t,t,-1);
      }
      if (j == r && arb_gt(t,c)) {
        for (j=0;j<r;j++) {
          arb_neg(t,maxc[j]);
          arb_union(t,t,maxc[j],prec2);
          arb_add(res[j],res[j],t,prec2);
        }
        break;
      }
    }

    arb_poly_one(g);
    arb_poly_one(f); arb_poly_shift_left(f,f,1);
    for (j=0;j<r;j++) {
      for (i=0;i<j;i++) {
        arb_poly_get_coeff_arb(t,g,i);
        arb_submul(res[j],t,res[i],prec2);
      }
      arb_set_si(t,-twomu[j]-1);
      arb_mul_2exp_si(t,t,-1);
      arb_poly_set_coeff_arb(f,0,t);
      arb_poly_mul(g,g,f,prec2);
    }

    arb_set_si(t,-2);
    arb_poly_set_coeff_arb(f,0,t);
    if (r & 1) arb_neg(c,c);
    for (;j<k;j++) {
      arb_mul(res[j],c,res[j-r],prec2);
      for (i=0;i<j;i++) {
        arb_poly_get_coeff_arb(t,g,i);
        arb_submul(res[j],t,res[i],prec2);
      }
      arb_poly_mul(g,g,f,prec2);
    }

    arb_one(t);
    for (j=2;j<k;j++) {
      arb_div_ui(t,t,j,prec2);
      arb_mul(res[j],res[j],t,prec2);
    }
    for (j=0;j<k;j++)
      arb_trim(res[j],res[j]);
  }

  // compute k >= r such that |eps^k*G^{(k)}(u)/k!| < thresh for all u
  // replace thresh by an interval containing eps^k*G^{(k)}(u)/k!
  static long taylor_terms(arb_t thresh,long twomu[],long r,
      arb_srcptr eps,long prec) {
    long j,k;
    static int init;
    static arb_t a,b,t,x,exppi2,four_pir;

    if (!init) {
      arb_init(a); arb_init(b);
      arb_init(t); arb_init(x);
      arb_init(exppi2); arb_init(four_pir);
      init = 1;
    }
    arb_const_pi(t,prec);
    arb_mul_2exp_si(exppi2,t,-1);
    arb_exp(exppi2,exppi2,prec);
    arb_mul_ui(t,t,r,prec);
    arb_inv(t,t,prec);
    arb_mul_2exp_si(four_pir,t,2);

    // b = (2*exp(Pi/2)/(Pi^2*r))^(1/4)
    arb_mul_2exp_si(b,exppi2,1);
    arb_div_ui(b,b,r,prec);
    arb_const_pi(t,prec);
    arb_mul(t,t,t,prec);
    arb_div(b,b,t,prec);
    arb_root_ui(b,b,4,prec);

    // a = sqrt(2)/b
    arb_sqrt_ui(a,2,prec);
    arb_div(a,a,b,prec);

    arb_one(x);
    for (j=0;j<r;j++) {
      // t = 1+r+r/2*(mu_j-1/2)
      arb_set_ui(t,4+r*(twomu[j]+3));
      arb_mul_2exp_si(t,t,-2);
      arb_gamma(t,t,prec);
      arb_root_ui(t,t,r,prec);
      arb_mul(x,x,t,prec);
      arb_pow_ui(t,b,twomu[j],prec);
      arb_mul(t,t,a,prec);
      if (!twomu[j]) arb_mul(t,t,exppi2,prec);
      arb_mul(x,x,t,prec);
    }
    arb_set_ui(t,r+1);
    arb_gamma(t,t,prec);
    arb_div(x,x,t,prec);

    arb_mul(t,four_pir,eps,prec);
    arb_pow_ui(t,t,r,prec);
    arb_mul(t,t,four_pir,prec);
    arb_mul(x,x,t,prec);
    arb_const_pi(t,prec);
    arb_div(x,x,t,prec);

    for (k=r;!arb_lt(x,thresh);) {
      k++;
      arb_one(t);
      for (j=0;j<r;j++) {
        // 1+r*(mu-1/2)/(2*k)
        arb_mul_si(t,t,4*k+r*(twomu[j]-1),prec);
        arb_div_si(t,t,4*k,prec);
      }
      arb_root_ui(t,t,r,prec);
      arb_mul(t,t,four_pir,prec);
      arb_mul(t,t,eps,prec);
      arb_mul(x,x,t,prec);
    }

    arb_set(thresh,x);
    return k;
  }

  // return B such that |G(u)| <= 2^{-B}, or 0 if u <= 0
  static long error_bound(long twomu[],long r,double u) {
    long j,j1,j2;
    double y,g,nu;

    if (u <= 0) return 0;

    // find two smallest mu_j
    for (j1=0,j=1;j<r;j++)
      if (twomu[j] < twomu[j1])
        j1 = j;
    for (j2=!j1,j=1;j<r;j++)
      if (j != j1 && twomu[j] < twomu[j2])
        j2 = j;

    y = exp(2*u/r);
    g = M_LN2*r/2+u/2-M_PI*r*y;
    for (j=0;j<r;j++) {
      nu = (twomu[j]-2)*0.25;
      if (j == j1 || j == j2) nu += 0.25;
      g += nu*log(y+nu/M_PI);
    }

    return (long)(-g/M_LN2);
  }

  static int isprime(long n) {
    long k,b=(long)sqrtl((long double)n);
    for (k=2;k<=b;k++)
      if (n % k == 0) return 0;
    return 1;
  }

  // compute C such that |a_n| <= Cn assuming Ramanujan
  // \prod_{p<r^{1/\alpha}}\prod_{j<\frac{r-1}{p^\alpha-1}}\frac{r+j-1}{jp^\alpha}
  static void coeff_bound(arf_t C,long r,long prec) {
    long j,p;

    arf_one(C);
    for (p=2;p<r;p++)
      if (isprime(p))
        for (j=1;(p-1)*j<r-1;j++) {
          arf_mul_ui(C,C,r+j-1,prec,ARF_RND_UP);
          arf_div_ui(C,C,j*p,prec,ARF_RND_UP);
        }
  }

  // compute G data into L
  // if(op) then also write the data to fp (in the cache dierctory)
  static void computeall(Lfunc *L, double umin,double Binv,long prec, bool op, FILE *fp) 
  {
    long i, j, k, prec2, imin, imax;
    double delta;
    arb_t *g;
    arb_t u,eps,thresh;
    arf_t m;
    long twomu[maxr];

    arb_init(u); arb_init(eps);
    arb_init(thresh); arf_init(m);

    for(i = 0; i < (long)L->degree; i++)
      twomu[i]=L->mus[i]*2.0;

    // clear any cached polar parts from last run
    for (i=0;i<mcache;i++)
      cache[i].order = -1;

    delta = 2*M_PI*Binv;
    arb_one(thresh);
    arb_mul_2exp_si(thresh,thresh,-prec);
    imin = (long)floor(umin/delta);
    for (imax=imin;error_bound(twomu,L->degree,imax*delta)<prec;imax++);

    coeff_bound(m,L->degree,53);
    arb_init(L->C);
    arb_set_arf(L->C,m);
    arb_init(L->alpha);
    arb_set_ui(L->alpha,1);

    if(op) {
      // Self-describing header (see GCACHE_MAGIC): version, degree, the gprec
      // this grid was computed at, and the sorted analytic mus as exact halved
      // integers.  read_gfile reads and verifies this before the body.
      fprintf(fp, "%s %d %lu %ld", GCACHE_MAGIC, GCACHE_VERSION,
              (unsigned long)L->degree, prec);
      for(i = 0; i < (long)L->degree; i++)
        fprintf(fp, " %ld", twomu[i]);
      fprintf(fp, "\n");
      printarf(fp,m);
      fprintf(fp," 1\n");
    }
    prec2 = prec + (long)(exp(2*imax*delta/L->degree)*M_PI*L->degree/M_LN2) + 100;
    arb_const_pi(u,prec2); arb_set_d(eps,Binv); arb_mul(eps,eps,u,prec2);
    k = taylor_terms(thresh,twomu,L->degree,eps,prec);
    L->max_K=k;
    L->one_over_B=Binv;
    if(verbose) printf("1/B set to %f\n",Binv);
    L->low_i=imin;
    L->hi_i=imax;

    arb_init(L->eq59);
    arb_set(L->eq59,thresh);
    if(op) {
      fprintf(fp,"%.9f %ld %ld\n", Binv, imin, imax);
      arb_get_abs_ubound_arf(m,thresh,prec);
      fprintf(fp,"%ld ",k);
      printarf(fp,m);
      fprintf(fp,"\n");
      fflush(fp);
    }

    L->Gs=(arb_t **)malloc(sizeof(arb_t *)*k);
    if(!L->Gs)
    {
      fprintf(stderr,"Fatal error allocating memory in computeall. Exiting.\n");
      exit(0);
    }
    for(i=0;i<k;i++)
    {
      L->Gs[i]=(arb_t *)malloc(sizeof(arb_t)*(imax-imin+1));
      if(!L->Gs[i])
      {
        fprintf(stderr,"Fatal error allocating memory in computeall. Exiting.\n");
        exit(0);
      }
    }
    for(i=0;i<k;i++)
      for(j=0;j<=imax-imin;j++)
        arb_init(L->Gs[i][j]);


    g = calloc(k,sizeof(g[0]));
    for (i=0;i<k;i++)
      arb_init(g[i]);
    int64_t ii;
    for (i=imin,ii=0;i<=imax;i++,ii++) {
      arb_mul_si(u,eps,2*i,prec2);
      gtaylor(g,twomu,L->degree,u,k,prec,prec2);
      for (j=0;j<k;j++) {
        arb_set(L->Gs[j][ii],g[j]);
        if( op ) {
          fprintf(fp, "%ld %ld ", i, j);
          printarb(fp,g[j]);
          fprintf(fp,"\n");
          fflush(fp);
        }
      }
    }
    for (i=0;i<k;i++)
      arb_clear(g[i]);
    free(g);
    arb_clear(u); arb_clear(eps);
    arb_clear(thresh); arf_clear(m);

  }

  bool read_arb(arb_ptr res, FILE *infile)
  {
    static bool init=false;
    static fmpz_t a,b;
    static mpz_t x,e;
    static arb_t radius;
    if(!init)
    {
      init=true;
      fmpz_init(a);
      fmpz_init(b);
      mpz_init(x);
      mpz_init(e);
      arb_init(radius);
    }


    if(!mpz_inp_str(x,infile,10))
      return(false);
    if(!mpz_inp_str(e,infile,10))
      return(false);
    fmpz_set_mpz(a,x);
    fmpz_set_mpz(b,e);
    arb_set_fmpz_2exp(res,a,b);

    if(!mpz_inp_str(x,infile,10))
      return(false);
    if(!mpz_inp_str(e,infile,10))
      return(false);
    fmpz_set_mpz(a,x);
    fmpz_set_mpz(b,e);
    arb_set_fmpz_2exp(radius,a,b);

    arb_add_error(res,radius);

    return(true);
  }


  // each line consists of <i> <k> <m1> <e1> <m2> <e2>
  // so Gs[i][k]=m1*2^e1 +/- m2*2^e2
  bool read_Gs(FILE *infile, Lfunc *L)
  {

    //int64_t prec=L->gprec;
    int64_t fi, i, j;
    uint64_t fk;

    for(i=L->low_i,j=0;i<=L->hi_i;i++,j++)
      for(uint64_t k=0;k<L->max_K;k++)
      {
        if(fscanf(infile,"%" PRId64 " %" PRIu64 "", &fi, &fk) != 2)
          return(false);
        if(fi!=i)
          return(false);
        if(fk!=k)
          return(false);
        if(!read_arb(L->Gs[k][j],infile))
          return(false);
      }
    return(true);
  }

  // Read and validate the GCACHE header (see GCACHE_MAGIC).  Returns false --
  // so the caller recomputes (overwriting the file) rather than returning wrong
  // numbers -- when the file is foreign or from an older format (magic/version
  // mismatch), describes a different L-function (degree or mus mismatch), or
  // was computed at too low a precision for the current request
  // (cached gprec < req_gprec; recall hi_i and max_K grow with gprec).  This
  // runs before read_gfile allocates anything, so a rejected file leaks nothing.
  static bool read_gheader(FILE *infile, const Lfunc *L, int64_t req_gprec)
  {
    char magic[16];
    int version;
    unsigned long degree;
    int64_t cached_gprec;
    if(fscanf(infile, "%15s %d %lu %" PRId64, magic, &version, &degree, &cached_gprec) != 4)
      return false;
    if(strcmp(magic, GCACHE_MAGIC) != 0)
      return false;
    if(version != GCACHE_VERSION)
      return false;
    if(degree != L->degree)
      return false;
    for(uint64_t i = 0; i < L->degree; i++) {
      long twomu;
      if(fscanf(infile, "%ld", &twomu) != 1)
        return false;
      if(twomu != (long)round(L->mus[i] * 2.0)) // identity: exact halved mu
        return false;
    }
    if(cached_gprec < req_gprec) // sufficiency: cached grid is precise enough
      return false;
    return true;
  }

  // Read a file written previously by computeall, validating it against the
  // current request first (see read_gheader).  Returns GCACHE_STALE -- so the
  // caller recomputes and overwrites rather than returning wrong numbers -- when
  // the header does not match (foreign / old format / different L-function /
  // insufficient precision); GCACHE_CORRUPT when the header is usable but the
  // body then fails to parse (a genuinely broken or truncated file, fatal); and
  // GCACHE_OK when the file is usable and fully read.
  gcache_status read_gfile(FILE *infile, Lfunc *L, int64_t req_gprec)
  {
    if(!read_gheader(infile, L, req_gprec)) // verify before allocating anything
      return GCACHE_STALE;
    int64_t m,e,alpha;
    //double dalpha;
    // C = coeff_bound(degree, 53) is canonical in the degree, which the header
    // has just verified, so the value read below is safe to trust.
    arb_init(L->C);
    arb_init(L->alpha);
    if(fscanf(infile,"%" PRId64 " %" PRId64 " %" PRId64 "\n",&m,&e,&alpha)!=3)
      return GCACHE_CORRUPT;
    arb_set_si(L->C,m);
    arb_mul_2exp_si(L->C,L->C,e);
    arb_set_si(L->alpha,alpha);
    if(fscanf(infile,"%lf %" PRId64 " %" PRId64 "\n",&L->one_over_B,&L->low_i,&L->hi_i)!=3)
      return GCACHE_CORRUPT;
    if(fscanf(infile,"%" PRIu64 "",&L->max_K)!=1)
      return GCACHE_CORRUPT;

    fmpz_t a,b;
    mpz_t x,ee;
    fmpz_init(a);
    fmpz_init(b);
    mpz_init(x);
    mpz_init(ee);
    arb_init(L->eq59);
    if(!mpz_inp_str(x,infile,10))
      return GCACHE_CORRUPT;
    if(!mpz_inp_str(ee,infile,10))
      return GCACHE_CORRUPT;
    fmpz_set_mpz(a,x);
    fmpz_set_mpz(b,ee);
    arb_set_fmpz_2exp(L->eq59,a,b);
    fmpz_clear(a);
    fmpz_clear(b);
    mpz_clear(x);
    mpz_clear(ee);

    L->Gs=(arb_t **)malloc(sizeof(arb_t *)*L->max_K);
    if(!L->Gs)
      return GCACHE_CORRUPT;
    for(uint64_t k=0;k<L->max_K;k++)
    {
      L->Gs[k]=(arb_t *)malloc(sizeof(arb_t)*(L->hi_i-L->low_i+1));
      if(!L->Gs[k])
        return GCACHE_CORRUPT;
      for(int64_t j=0;j<=L->hi_i-L->low_i;j++)
        arb_init(L->Gs[k][j]);
    }
    return read_Gs(infile,L) ? GCACHE_OK : GCACHE_CORRUPT;

  }

  Lerror_t compute_g(Lfunc *L)
  {

    Lerror_t ecode=ERR_SUCCESS;
    bool op = false; // true if we are going to write a cache file
    FILE *ofile = NULL; // file to write to

    // "Default mode" means the caller left gprec at 0; we then both derive the
    // working gprec ourselves and may consult the on-disk cache.  Capture that
    // state before we overwrite L->gprec below.
    bool use_cache = (L->gprec == 0) && (L->target_prec == DEFAULT_TARGET_PREC)
                     && (L->cache_dir);

    // Derive the required gprec up front (previously done only on a cache
    // miss), so it is available to validate a cached file's sufficiency.  The
    // wprec floor matters: a caller who raises wprec for a large conductor
    // needs a higher-precision grid than a default run wrote.
    if(L->gprec == 0)
    {
      double gfac = 0.0;
      for(uint64_t d=0; d < L->degree; d++)
        gfac += lgamma(0.25 + L->mus[d]/2.0);
      gfac /= M_LN2;
      L->gprec = L->target_prec + ceil(gfac) + EXTRA_BITS;
      if(L->gprec < L->wprec)
        L->gprec = L->wprec;
    }
    if(verbose)
      printf("g precision set to %" PRId64 " bits\n", L->gprec);

    if(use_cache) {
      char fname[1337];
      char fname1[1024] = "";
      size_t off = 0;
      // Build the cache filename, bailing out if any snprintf hits an encoding
      // error (returns < 0) or would truncate. snprintf reports the would-be
      // length, so capture it in an int (off is size_t: adding a negative return
      // directly would wrap) and advance only once a write fully fit. A truncated
      // name would drop trailing mus / path chars and let distinct L-functions
      // collide on one cache file, so on any failure we skip the cache
      // (best-effort) and compute the G data fresh rather than read or write the
      // wrong data. These checks bound both writes for any mus and cache_dir --
      // nothing here relies on a cap on the number or magnitude of the mus.
      bool name_fits = true;
      for(uint64_t r=0; r<L->degree; r++) {
        int n = snprintf(fname1+off, sizeof(fname1)-off, "_%.1f", L->mus[r]);
        if(n < 0 || (size_t) n >= sizeof(fname1)-off) {
          name_fits = false;
          break;
        }
        off += (size_t) n;
      }
      if(name_fits) {
        int n = snprintf(fname, sizeof(fname), "%s/g%s", L->cache_dir, fname1);
        name_fits = (n >= 0) && ((size_t) n < sizeof(fname));
      }
      if(name_fits) {
        FILE *infile = fopen(fname, "r");
        if(infile) // a cache file with this name already exists
        {
          gcache_status st = read_gfile(infile, L, L->gprec); // validate + read
          fclose(infile);
          if(st == GCACHE_OK) // header verified and body read: usable as-is
            return ecode;
          if(st == GCACHE_CORRUPT) // header usable but body unparsable: broken file
            return ecode|ERR_G_INFILE;
          // GCACHE_STALE: foreign / old version / different L-function /
          // insufficient precision.  Fall through to recompute the G data and
          // overwrite the stale file rather than reuse it or abort.
        }
        // cache miss, or a stale file we are about to overwrite
        ofile = fopen(fname, "w"); // truncates any stale file
        if( !ofile ) {
          ecode |= ERR_G_OUTFILE; // couldn't open outfile. Not fatal
          op = false;
        } else {
          op = true; // file open ok so we can output
        }
      }
    }

    computeall(L, -32*M_LN2, (double)L->degree/512, L->gprec, op, ofile);
    if(op)
      fclose(ofile);
    return ecode;
  }

#ifdef __cplusplus
}
#endif
