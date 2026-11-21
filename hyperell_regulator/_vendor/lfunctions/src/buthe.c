#include <flint/acb_poly.h>
#include "glfunc_internals.h"

#ifdef __cplusplus
extern "C"{
#endif
#ifdef BUTHE


  // We use Buthe's method to verify our list of zeros
  // we set h=8
  // and use the zeros up to height 96/r to check the list to height buthe_b=64/r
  // This is probably huge overkill


  // setup Buthe zero check stuff
  void init_buthe(Lfunc *L, int64_t prec) {
    /*
       buthe_ints.out contains output from buthe_ints.gp
       each entry is a value of
       2\int\limits_0^\infty \frac{\exp(-(1/2+\mu_k)t)}{1-\exp(-2t)}\left[\frac{b}{\pi}-\frac{\sin(bt)}{\pi t\cosh(2t)} \dif t
       where the \mu_k are the (analytic) mus so non-negative 1/2 integers
       each "row" consists of 2*MAX_MUI_2+1 entries corresponding to \mu=0,0.5,1.0...
       a row for each degree (and therefore b) from 2
       The entries are integers which must be multiplied by 2^BUTHE_INT_SHIFT
       Once the entries in a row get small enough, we pad with zeros.
       If mu is outside MAX_MU_2 or the relevant entry is zero, use the nearest preceding non-zero entry uinioned with 0.
       */

#define BUTHE_INT_SHIFT (-20)
    uint64_t buthe_r_ints [(MAX_R-1)*(2*MAX_MUI_2+1)]=
#include "../gp/buthe_ints.out"

    arb_t one;
    arb_init(one);
    arb_set_ui(one,1);
    for(uint64_t i=0;i<(MAX_R-1)*(2*MAX_MUI_2+1);i++)
    {
      arb_init(L->buthe_ints[i]);
      arb_set_ui(L->buthe_ints[i],buthe_r_ints[i]);
      arb_add_error(L->buthe_ints[i],one);
      arb_mul_2exp_si(L->buthe_ints[i],L->buthe_ints[i],BUTHE_INT_SHIFT);
    }
    arb_init(L->buthe_Wf);
    arb_init(L->buthe_Winf);
    arb_init(L->buthe_Ws);
    arb_init(L->buthe_b); // will confirm RH in [0,b]
    arb_init(L->buthe_C);
    arb_init(L->buthe_h);


    arb_set_ui(L->buthe_C,L->degree);
    arb_div_ui(L->buthe_b,L->B,OUTPUT_RATIO,prec);
    //printf("Buthe b set to ");arb_printd(L->buthe_b,20);printf("\n");
    //
    // we set buthe_h as big as possible
    // we need b-a>5h/Pi were b-a =512/16/10
    // 
    arb_set_ui(L->buthe_h,BUTHE_H); // defined in buthe_ints.out
    //printf("Buthe h set to ");arb_printd(L->buthe_h,20);printf("\n");

    // check h<pi*2*b/5 will be OK so long as degree <=10
    arb_t tmp;
    arb_init(tmp);
    arb_mul(tmp,L->pi,L->buthe_b,prec);
    arb_div_ui(tmp,tmp,5,prec);
    arb_mul_2exp_si(tmp,tmp,1);
    if(verbose) {
      printf("Max allowed h is ");
      arb_printd(tmp,20);
      printf("\n");
      fflush(stdout);
    }
    arb_sub(tmp,tmp,L->buthe_h,prec);
    if(!arb_is_positive(tmp))
    {
      fprintf(stderr,"Buthe h is too large at ");
      arb_fprintd(stderr,L->buthe_h,20);
      fprintf(stderr,"for Buthe B = ");
      arb_fprintd(stderr,L->buthe_b,20);
      fprintf(stderr,". Exiting.\n");
      exit(0);
    }

    arb_clear(tmp);
    arb_clear(one);
  }

  // the terms making up wf are computed using the coefficients bm
  // whilst processing the Euler Factors. see use_inv_lpoly in coeff.c
  //
  // compute a term for wf
  // bm is real part of the algebraic coefficient for p^m in -L'/L(s)=sum_{n>1} Lambda(n) b(n)n^{-s}
  // subtract 2 Re(log(p)bm f(mlog p)/sqrt(p^m) from Wf
  // Re f(mlog(p))=sin(b mlog(p))/(Pi m log(p) cosh(h mlog(p)/2))
  // so we want 2 bm sin(b m log(p))/(Pi m cosh(h m log(p)/2))
  void wf1(Lfunc *L, uint64_t m, uint64_t pm, arb_t bm, int64_t prec)
  {
    if(arb_is_zero(bm)) // nothing to do
      return;

    arb_t s,tmp,tmp1,tmp2,logp,logpm;
    arb_init(s);
    arb_init(tmp);
    arb_init(tmp1);
    arb_init(tmp2);
    arb_init(logp);
    arb_init(logpm);


    arb_log_ui(logpm,pm,prec);

    arb_mul(tmp1,logpm,L->buthe_h,prec);
    arb_mul_2exp_si(tmp1,tmp1,-1);
    arb_cosh(tmp2,tmp1,prec); // cosh(hm/2 log p)
    arb_mul(tmp1,tmp2,L->pi,prec);
    arb_mul_ui(tmp2,tmp1,m,prec); // pi m cosh(.)
    arb_sqrt_ui(tmp,pm,prec);
    arb_mul(tmp1,tmp,tmp2,prec); // sqrt(p^m) pi m cosh(.)
    arb_mul(tmp,L->buthe_b,logpm,prec);
    arb_sin(s,tmp,prec);
    arb_mul(tmp,s,bm,prec); // b(p^m) sin(.)
    arb_div(tmp2,tmp,tmp1,prec);
    arb_mul_2exp_si(tmp2,tmp2,1);
    arb_add(L->buthe_Wf,L->buthe_Wf,tmp2,prec);

    arb_clear(s);
    arb_clear(tmp);
    arb_clear(tmp1);
    arb_clear(tmp2);
    arb_clear(logp);
    arb_clear(logpm);

    return;
  }


  // fp1 is f^-1, fp is f
  void wf(Lfunc *L, uint64_t p, acb_poly_t fp1, acb_poly_t fp, int64_t prec)
  {
    static bool init=false;
    static acb_t tmp1,acm;
    if(!init)
    {
      init=true;
      acb_init(tmp1);
      acb_init(acm);
    }

    uint64_t pm=p,m=1;
    //if(verbose) if(p==2) {printf("p=%lu\n",p);acb_poly_printd(fp,20);printf("\n--------\n");acb_poly_printd(fp1,20);printf("\n--------\n");}
    while(pm<=L->buthe_M)
    {
      acb_zero(acm);
      // this is lazy. if m-1 exceeds length of fp1, adjust i,j accordingly
      for(int64_t i=1,j=m-1;(i<acb_poly_length(fp))&&(j>=0);i++,j--)
      {
        if(j>=acb_poly_length(fp1)) continue; // and avoid this
        acb_mul(tmp1,acb_poly_get_coeff_ptr(fp1,j),acb_poly_get_coeff_ptr(fp,i),prec);
        acb_mul_si(tmp1,tmp1,i,prec);
        acb_add(acm,acm,tmp1,prec);
      }
      wf1(L,m,pm,acb_realref(acm),prec);
      pm*=p;
      m++;
    }
  }

  // We assume that L=exp(sum c(p^m) p^{-ms})
  // and that |c(p^m)|<=rp^(m/2)
  // sigma_0=1, sigma_1=3/2, C=r
  // see Lemma 3.4 of Buthe
  void buthe_Wf_error(Lfunc *L)
  {
    int64_t prec=L->wprec;
    uint64_t r=L->degree;
    arb_t tmp,tmp1,tmp2;
    arb_init(tmp);
    arb_init(tmp1);
    arb_init(tmp2);
    if(verbose)
      {
	printf("In buthe_Wf_error with M=%lu\n",L->buthe_M);
	printf("   and h = ");arb_printd(L->buthe_h,20);printf("\n");
      }
    
    arb_log_ui(tmp,L->buthe_M,prec);
    arb_sub_ui(tmp1,L->buthe_h,2,prec); // h-2 = 6
    if(verbose) {printf("h-2 = ");arb_printd(tmp1,10);printf("\n");}
    arb_mul_2exp_si(tmp1,tmp1,-1); // (h-2)/2 = 3
    arb_mul(tmp2,tmp,tmp1,prec); // (h-2)/2 log M
    arb_neg(tmp2,tmp2); // (2-h)/2 log M
    arb_exp(tmp,tmp2,prec); // M^((2-h)/2)
    if(verbose)
      {
	printf("M^((2-h)/2)=");
	arb_printd(tmp,20);
	printf("\n");
      }
    arb_div(tmp2,tmp,tmp1,prec); // 2M^()/(h-2)
    arb_mul_2exp_si(tmp2,tmp2,2); // 8 M^()/(h-2)
    arb_mul_ui(tmp,tmp2,r,prec); // 8 r M^()/(h-2)
    arb_div(tmp1,tmp,L->pi,prec); // /pi
    if(verbose){printf("error in Buthe Wf <= ");arb_printd(tmp1,20);printf("\n");}
    arb_add_error(L->buthe_Wf,tmp1);
    arb_clear(tmp);
    arb_clear(tmp1);
    arb_clear(tmp2);

  }

  void buthe_fhat(arb_t res, arb_t z, Lfunc *L, int64_t prec)
  {
    static bool init=false;
    static arb_t tmp,tmp1,tmp2,tmp3;
    if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
    }
    arb_div(tmp,L->pi,L->buthe_h,prec); // Pi/h
    arb_add(tmp1,z,L->buthe_b,prec); // z+b
    arb_mul(tmp2,tmp,tmp1,prec); // Pi/h z
    arb_exp(tmp1,tmp2,prec); // exp(Pi/h (z+b))
    arb_atan(tmp2,tmp1,prec); // atan(exp(Pi/h (z+b)))
    arb_sub(tmp1,z,L->buthe_b,prec); // z-b
    arb_mul(tmp3,tmp1,tmp,prec); // Pi/h(z-b)
    arb_exp(tmp,tmp3,prec); // exp(Pi/h(z-b))
    arb_atan(tmp3,tmp,prec); // atan(exp(Pi/h(z-b)))
    arb_sub(tmp,tmp2,tmp3,prec);
    arb_div(res,tmp,L->pi,prec);
    arb_mul_2exp_si(res,res,1); // 2/Pi*(atan - atan)
  }


  // sum over zeros for a dual l-function
  void buthe_Ws_dual(arb_t res, Lfunc *Lf, arb_t *zeros, int64_t prec)
  {
    static bool init=false;
    static arb_t tmp,tmp1;
    if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
    }
    for(uint64_t z=0;;z++)
    {
      if((z==MAX_ZEROS)||(arb_is_zero(zeros[z])))
        break;
      arb_set(tmp1,zeros[z]);
      buthe_fhat(tmp,tmp1,Lf,prec);
      if(verbose) {
        printf("Zero at ");arb_printd(tmp1,20);printf(" contributed ");arb_printd(tmp,20);printf("\n");
      }
      arb_add(res,res,tmp,prec);
    }
    arb_mul_2exp_si(res,res,1);
    if(Lf->rank>0)
    {
      arb_zero(tmp1);// we assume the rank is correct so zeros are at exactly 0
      buthe_fhat(tmp,tmp1,Lf,prec);
      arb_mul_ui(tmp1,tmp,Lf->rank,prec);
      arb_add(res,res,tmp1,prec);
    }
  }

  // sum over zeros for a non-dual l-function (will work with dual as well)
  void buthe_Ws_non_dual(arb_t res, Lfunc *Lf, arb_t *zeros, uint64_t side, int64_t prec)
  {
    static bool init=false;
    static arb_t tmp,tmp1;
    if(!init) {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
    }
    for(uint64_t z=0;;z++)
    {
      if((z==MAX_ZEROS)||(arb_is_zero(zeros[z])))
        break;
      arb_set(tmp1,zeros[z]);
      buthe_fhat(tmp,tmp1,Lf,prec);
      if(verbose) {
        printf("Zero at ");arb_printd(tmp1,20);printf(" contributed ");arb_printd(tmp,20);printf("\n");
      }
      arb_add(res,res,tmp,prec);
    }
    if((side==0)&&(Lf->rank>0)) // there are some central zeros so include them
    {
      arb_zero(tmp1); // we assume the rank is correct so zeros are at exactly 0
      buthe_fhat(tmp,tmp1,Lf,prec);
      arb_mul_ui(tmp1,tmp,Lf->rank,prec);
      arb_add(res,res,tmp1,prec);
    }
  }

  // add digamma(1/4+mu/2)
  void buthe_lgam1(arb_t res, double mu, int64_t prec)
  {
    static bool init=false;
    static arb_t s,tmp;
    if(!init)
    {
      arb_init(s);
      arb_init(tmp);
    }
    arb_set_d(s,1.0/4.0+(double)mu/2.0); // all normalised to (1-s)
    arb_digamma(tmp,s,prec);
    //if(verbose) {printf("lgam1(%f) returning ",mu);arb_printd(tmp,20);printf("\n");}
    arb_add(res,res,tmp,prec);
  }

  // add log N - r log pi
  void buthe_lgam2(arb_t res, uint64_t r, uint64_t N, arb_t logpi, int64_t prec)
  {
    static bool init=false;
    static arb_t tmp,tmp1,tmp2,tmp3;
    if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(tmp3);
    }
    arb_log_ui(tmp1,N,prec); // log N
    arb_mul_ui(tmp2,logpi,r,prec); // log Pi^r
    arb_sub(tmp3,tmp1,tmp2,prec); // log(N/Pi^r)
    if(verbose){printf("lgam2 adding ");arb_printd(tmp3,20);printf("\n");}
    arb_add(res,res,tmp3,prec);
  }

  void buthe_Winf(arb_t res, Lfunc *L, int64_t prec)
  {
    static bool init=false;
    static arb_t tmp,tmp1,logpi;
    if(!init) {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(logpi);
    }
    arb_log(logpi,L->pi,prec);
    arb_zero(res);
    for(uint64_t k=0;k<L->degree;k++)
      buthe_lgam1(res,L->mus[k],prec);
    buthe_lgam2(res,L->degree,L->conductor,logpi,prec);
    arb_div(tmp,L->buthe_b,L->pi,prec);
    arb_mul(res,res,tmp,prec);
    if(verbose) {printf("Winf before integral = ");arb_printd(res,20);printf("\n");}
    for(uint64_t k=0;k<L->degree;k++) {
      uint64_t ptr;
      bool approx;
      if(L->mus[k]>MAX_MU)
      {
        ptr=(L->degree-1)*(2*MAX_MUI_2+1)-1;
        approx=true;
      }
      else
      {
        ptr=(L->degree-2)*(2*MAX_MUI_2+1)+(uint64_t)(2.0*L->mus[k]);
        approx=arb_is_zero(L->buthe_ints[ptr]);
      }
      if(approx)
      {
        while(arb_is_zero(L->buthe_ints[ptr])) ptr--;
        arb_set(tmp,L->buthe_ints[ptr]);
        arb_zero(tmp1);
        arb_union(tmp,tmp,tmp1,prec);
      }
      else
        arb_set(tmp,L->buthe_ints[ptr]);
      //if(verbose){printf("Integral contributing ");arb_printd(buthe_ints[(Lf->r-2)*(2*MAX_MUI_2+1)+(uint64_t)(2.0*Lf->mus[k])],20);printf("\n");}
      arb_add(res,res,tmp,prec);
    }
    if(verbose) {printf("Winf = ");arb_printd(res,20);printf("\n");}
  }

  Lerror_t buthe_check_RH(Lfunc *L)
  {
    static bool init=false;
    static arb_t sum,two_zeros;
    if(!init) {
      init=true;
      arb_init(sum);
      arb_init(two_zeros);
      arb_set_ui(two_zeros,98);
      arb_div_ui(two_zeros,two_zeros,100,100);
    }
    int64_t prec=L->wprec;
    if(verbose)
    {printf("Going to use Weil-Barner to confirm list of zeros.\n");fflush(stdout);}
    arb_zero(L->buthe_Ws);
    if(L->self_dual==YES)
      buthe_Ws_dual(L->buthe_Ws,L,L->zeros[0],L->wprec);
    else // unknown or definately not
    {
      buthe_Ws_non_dual(L->buthe_Ws,L,L->zeros[0],0,L->wprec);
      buthe_Ws_non_dual(L->buthe_Ws,L,L->zeros[1],1,L->wprec);
    }

    buthe_Winf(L->buthe_Winf,L,L->wprec);
    //if(verbose){printf("Winf = ");arb_printd(L_func.buthe_Winf,20);printf("\n");}


    //printf("Random negation of Wf happening.\n");arb_neg(L->buthe_Wf,L->buthe_Wf);
    arb_add(sum,L->buthe_Wf,L->buthe_Winf,prec);
    arb_sub(sum,sum,L->buthe_Ws,prec);
    if(verbose)
    {
      printf("Ws = ");arb_printd(L->buthe_Ws,20);
      printf("\nWinf = ");arb_printd(L->buthe_Winf,20);
      printf("\nWf = ");arb_printd(L->buthe_Wf,20);
      printf("\nS = ");arb_printd(sum,20);printf("\n");
    }
    if(arb_is_negative(sum))
    {
      if(verbose) printf("Error in Weil-Barner check. Winf+Wf-Ws* must allow >=0. RH not confirmed.\n");
      return ERR_BUT_ERROR;
    }

    arb_sub(sum,sum,two_zeros,prec);

    if(!arb_is_negative(sum))
    {
      if(verbose) printf("Looks like we've missed a some (pairs of) zeros.\n");
      return ERR_RH_ERROR;
    }
    return ERR_SUCCESS;
  }

#ifdef __cplusplus
}
#endif

#endif
