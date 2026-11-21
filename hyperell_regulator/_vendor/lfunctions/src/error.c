#include "glfunc_internals.h"
#include <flint/acb.h>
#include <flint/acb_hypgeom.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
extern "C"{
#endif

// Lets the experimental bridge require this error-bound correction at link
// time, rather than accidentally accepting an unpatched shared library.
int Lfunc_error_bounds_revision(void) { return 1; }

// return Q(s) is analytic conductor defined pg 387 col 2
void Q(acb_t res, acb_t s, Lfunc *L, uint64_t conductor, uint64_t prec)
{
  //printf("In logQ with s=");acb_printd(s,10);printf("\n");
  arb_t two_pi,tmp1,tmp2;
  acb_t stmp1,stmp2,stmp3;
  arb_init(two_pi);arb_init(tmp1);arb_init(tmp2);
  acb_init(stmp1);acb_init(stmp2);acb_init(stmp3);
  arb_mul_2exp_si(two_pi,L->pi,1);

  acb_set_ui(stmp1,conductor);
  arb_set(acb_imagref(stmp2),acb_imagref(s));
  uint64_t j;
  for(j=0;j<L->degree;j++)
  {
    arb_set_d(tmp2,L->mus[j]);
    arb_add(acb_realref(stmp2),acb_realref(s),tmp2,prec);
    acb_mul(stmp3,stmp1,stmp2,prec);
    acb_div_arb(stmp1,stmp3,two_pi,prec);
  }
  //printf("Analytic conductor = ");acb_printd(stmp1,10);printf("\n");
  acb_set(res,stmp1);
  arb_clear(two_pi);arb_clear(tmp1);arb_clear(tmp2);
  acb_clear(stmp1);acb_clear(stmp2);acb_clear(stmp3);
}

// |Q(s)|
void absQ(arb_t res, acb_t s, Lfunc *L, uint64_t conductor, uint64_t prec)
{
  acb_t tmp1;
  acb_init(tmp1);
  Q(tmp1,s,L,conductor,prec);
  acb_abs(res,tmp1,prec);
  acb_clear(tmp1);
}

// returns |Pi^(-s/2) Gamma[s/2]|
void abs_gamma_r(arb_t res, acb_t s, Lfunc *L, uint64_t prec)
{
  acb_t s_by_2,lg;
  arb_t tmp,tmp0,log_pi;
  acb_init(s_by_2);
  acb_init(lg);
  arb_init(tmp);
  arb_init(tmp0);
  arb_init(log_pi);
  arb_log(log_pi,L->pi,prec);
  //printf("in abs_gamma_r with s = ");acb_printd(s,10);printf("\n");
  acb_mul_2exp_si(s_by_2,s,-1); // s/2
  acb_lgamma(lg,s_by_2,prec);
  //printf("loggamma(s/2) = ");acb_printd(lg,10);printf("\n");

  arb_mul(tmp,log_pi,acb_realref(s_by_2),prec);

  arb_sub(tmp0,acb_realref(lg),tmp,prec);
  arb_exp(res,tmp0,prec);
  //printf("gamma_r returning ");arb_printd(res,prec);printf("\n");
  acb_clear(s_by_2);
  acb_clear(lg);
  arb_clear(tmp);
  arb_clear(tmp0);
  arb_clear(log_pi);
}


// |gamma(s)| per top col 2 pg 387 without N^1/2(s-1/2)
// returns \prod |Gamma_R(s + u_j)|
void abs_gamma(arb_t res, acb_t s, Lfunc *L, int64_t prec)
{
  acb_t tmp1,tmp2,tmp3;
  acb_init(tmp1);
  acb_init(tmp2);
  acb_init(tmp3);
  arb_t tmp;
  arb_init(tmp);

  arb_set_ui(res,1);
  for(uint64_t j = 0; j < L->degree; j++)
  {
    arb_set_d(tmp,L->mus[j]);
    acb_add_arb(tmp1,s,tmp,prec);
    abs_gamma_r(tmp,tmp1,L,prec);
    arb_mul(res,res,tmp,prec);
  }
  //printf("abs_gamma returning ");arb_printd(res,10);printf("\n");
  acb_clear(tmp1);
  acb_clear(tmp2);
  acb_clear(tmp3);
  arb_clear(tmp);
}

// Lemma 5.7 of Artin's conjecture,....
// we compute this for the family by leaving out the conductor
void init_ftwiddle_error(Lfunc *L, int64_t prec)
{
  acb_t s,s_plus_mu;
  arb_t tmp1,tmp2,tmp3,tmp4,E,two_pi,four_by_pi2,beta;
  acb_init(s);
  acb_init(s_plus_mu);
  arb_init(tmp1);
  arb_init(tmp2);
  arb_init(tmp3);
  arb_init(tmp4);
  arb_init(E);
  arb_init(two_pi);
  arb_init(four_by_pi2);
  arb_init(beta);
  arb_const_pi(two_pi,prec);
  arb_inv(tmp1,two_pi,prec);
  arb_mul(four_by_pi2,tmp1,tmp1,prec);
  arb_mul_2exp_si(four_by_pi2,four_by_pi2,2); // 4/Pi^2
  arb_mul_2exp_si(two_pi,two_pi,1); // 2Pi

  arb_set_d(acb_realref(s),0.5);
  arb_set(acb_imagref(s),L->B); // 1/2+iB

  arb_set_d(tmp2,1.5);
  arb_zeta(tmp1,tmp2,prec);
  arb_pow_ui(tmp2,tmp1,L->degree,prec);

  abs_gamma(tmp1,s,L,prec);
  arb_mul(tmp3,tmp1,tmp2,prec);
  /*
     arb_set(acb_imagref(s_plus_mu),acb_imagref(s));
     for(uint64_t j=0;j<L->degree;j++)
     {
     arb_set_d(tmp1,L->mus[j]);
     arb_add(acb_realref(s_plus_mu),acb_realref(s),tmp1,prec);
     acb_abs(tmp2,s_plus_mu,prec);
     arb_div(tmp1,tmp2,two_pi,prec);
     arb_mul(E,E,tmp1,prec);
     }
     */
  absQ(tmp1,s,L,1,prec);
  arb_sqrt(tmp2,tmp1,prec);
  arb_mul(E,tmp2,tmp3,prec);

  // E without the N

  arb_mul_2exp_si(E,E,1);



  arb_mul_2exp_si(tmp1,two_pi,-3); //pi/4
  arb_mul_ui(beta,tmp1,L->degree,prec); // pi r/4
  uint64_t j;
  for(j=0;j<L->degree;j++)
  {
    arb_set_d(tmp1,0.5+L->mus[j]);
    arb_div(tmp2,tmp1,L->B,prec);
    arb_atan(tmp1,tmp2,prec);
    arb_sub(beta,beta,tmp1,prec);
  }

  arb_mul(tmp2,L->B,L->B,prec);
  arb_zero(tmp4);
  for(j=0;j<L->degree;j++)
  {
    arb_set_d(tmp1,(0.5+L->mus[j])*(0.5+L->mus[j]));
    arb_sub(tmp3,tmp2,tmp1,prec);
    arb_inv(tmp1,tmp3,prec);
    arb_add(tmp4,tmp4,tmp1,prec);
  }
  arb_mul(tmp1,tmp4,four_by_pi2,prec);
  arb_sub(beta,beta,tmp1,prec);

  arb_mul(tmp1,beta,L->B,prec);
  arb_neg(tmp1,tmp1);
  arb_exp(tmp2,tmp1,prec);
  arb_sub_ui(tmp1,tmp2,1,prec);
  arb_neg(tmp1,tmp1);
  arb_div(L->pre_ftwiddle_error,E,tmp1,prec);
  //printf("pre F twiddle error set to ");arb_printd(L->pre_ftwiddle_error,10);
  //printf("\n");

  acb_clear(s);
  acb_clear(s_plus_mu);
  arb_clear(tmp1);
  arb_clear(tmp2);
  arb_clear(tmp3);
  arb_clear(tmp4);
  arb_clear(two_pi);
  arb_clear(E);
  arb_clear(four_by_pi2);
  arb_clear(beta);

}

// multiply the error precomputed in L_family by sqrt{conductor} in Q(s)
void complete_ftwiddle_error(Lfunc *L, int64_t prec)
{
  arb_t tmp;
  arb_init(tmp);
  arb_sqrt_ui(tmp,L->conductor,prec);
  arb_mul(L->ftwiddle_error,L->pre_ftwiddle_error,tmp,prec);
  arb_clear(tmp);
}

// given n, what is the last a_m have we summed into res[n]
// we need to add the error for sum m>=m+1
// returns 0 if no coefficients have been added in
// then use F_hat_twiddle bound
uint64_t inv_m(uint64_t hi_i, uint64_t n, double one_over_B, double dc)
{
  return(exp(((double) hi_i-(double) n+0.5)*2.0*M_PI*one_over_B)*dc);
}

// see M_error1.tex Lemma 2
bool M_error_at_cutoff(arb_t res, arb_t x, Lfunc *L, uint64_t cutoff, int64_t prec)
{
  arb_t k1,tmp,tmp1,tmp2,tmp3,k0,rootNM2r,k2,k1k0,two_r;
  acb_t ctmp1,ctmp2,ctmp3;
  acb_init(ctmp1);
  acb_init(ctmp2);
  acb_init(ctmp3);
  arb_init(k1);
  arb_init(tmp);
  arb_init(tmp1);
  arb_init(tmp2);
  arb_init(tmp3);
  arb_init(two_r);
  arb_set_ui(tmp,2);
  arb_div_ui(two_r,tmp,L->degree,prec); // 2/r
  arb_init(k0);
  arb_init(rootNM2r); // (sqrt(N)/M)^(2/r)
  arb_init(k2);
  arb_mul_2exp_si(tmp,L->alpha,1);
  arb_add_ui(tmp1,tmp,1,prec);
  arb_mul_2exp_si(tmp1,tmp1,-2);
  arb_set_ui(tmp,L->conductor);
  arb_pow(tmp2,tmp,tmp1,prec); // N^((2alpha+1)/4)
  arb_mul(tmp,tmp2,L->C,prec); // CN...
  arb_mul_ui(tmp1,tmp,L->degree,prec); // rCN...
  arb_log_ui(tmp,2,prec);
  arb_set_si(tmp2,L->degree-2);
  arb_mul_2exp_si(tmp2,tmp2,-1);
  arb_mul(tmp3,tmp,tmp2,prec);
  arb_exp(tmp,tmp3,prec); // 2^((r-2)/2)
  arb_mul(k2,tmp,tmp1,prec); // 2^((r-2)/2) r C N....
  arb_init(k1k0);

  arb_set_d(tmp,0.5);
  arb_add(tmp1,tmp,L->alpha,prec);
  arb_add(tmp,tmp1,L->nu,prec);
  arb_mul_ui(k0,tmp,L->degree,prec);
  arb_mul_2exp_si(k0,k0,-1); // k0 = r(nu+alpha+1/2)/2
  arb_div_ui(tmp,x,L->degree,prec);
  arb_mul_2exp_si(tmp,tmp,1);
  arb_exp(tmp1,tmp,prec);
  arb_mul(tmp,tmp1,L->pi,prec);
  arb_mul_ui(k1,tmp,L->degree,prec); // k1 = pi r exp(2x/r)

  arb_log(tmp,k1,prec);
  arb_mul(tmp1,tmp,k0,prec);
  arb_neg(tmp1,tmp1);
  arb_exp(k1k0,tmp1,prec); // k1^(-k0)

  arb_sqrt_ui(tmp,L->conductor,prec);
  arb_div_ui(tmp1,tmp,cutoff,prec);
  arb_pow(rootNM2r,tmp1,two_r,prec); // (sqrt(N)/M)^(2/r)

  arb_set_d(tmp,-0.5);
  arb_add(tmp1,tmp,L->alpha,prec);
  arb_add(tmp,tmp1,L->nu,prec);
  arb_div(tmp1,L->pi,tmp,prec);
  arb_mul_2exp_si(tmp1,tmp1,1);
  if(!arb_ge(tmp1, rootNM2r)) {
    //printf("Failed 2piM^(2/r)... test in M_error.\n");
    //arb_printd(rootNM2r,20);printf("\n");
    return false;
  }

  arb_mul(tmp,L->nu,x,prec);
  arb_exp(tmp1,tmp,prec);
  arb_mul(tmp,tmp1,k2,prec);

  arb_div(acb_realref(ctmp1),k1,rootNM2r,prec);
  arb_set(acb_realref(ctmp2),k0);
  acb_hypgeom_gamma_upper(ctmp3,ctmp2,ctmp1,0,prec); // Gamma(k0,k1(M/sqrt(N))^(2/r))
  arb_mul(tmp1,tmp,acb_realref(ctmp3),prec);

  arb_mul(tmp,tmp1,k1k0,prec);

  arb_mul_ui(tmp1,rootNM2r,L->degree,prec);
  arb_div(tmp3,tmp1,k1,prec);
  for(uint64_t j=0; j < L->degree; j++)
  {
    arb_mul(tmp1,tmp3,L->nus[j],prec);
    arb_add_ui(tmp2,tmp1,1,prec);
    arb_pow(tmp1,tmp2,L->nus[j],prec);
    arb_mul(tmp,tmp,tmp1,prec);
  }
  arb_set(res,tmp);
  //printf("M_error with M=%" PRIu64 " x=",L->M);arb_printd(x,20);printf(" returning ");arb_printd(res,20);printf("\n");

  acb_clear(ctmp1);
  acb_clear(ctmp2);
  acb_clear(ctmp3);
  arb_clear(k1);
  arb_clear(tmp);
  arb_clear(tmp1);
  arb_clear(tmp2);
  arb_clear(tmp3);
  arb_clear(two_r);
  arb_clear(k0);
  arb_clear(rootNM2r); // (sqrt(N)/M)^(2/r)
  arb_clear(k2);
  arb_clear(k1k0);
  return true;
}

bool M_error(arb_t res, arb_t x, Lfunc *L, int64_t prec)
{
  return M_error_at_cutoff(res, x, L, L->M, prec);
}

#define BIG_ZETA (10)
#define BIG_ZETAD (1.001) // zeta([10,infinty]) in [1,1.001]
void djp_zeta(arb_t res, const arb_t x, slong prec)
{
  // big_zeta / big_zetad are process-lifetime cached constants; init once.
  // (The previous code never set init=true, so it re-ran arb_init every call,
  // leaking big_zetad's heap mantissa each time -- numerically a no-op since
  // the constants are recomputed identically.)
  static bool init=false;
  static arb_t big_zeta,big_zetad;
  if(!init)
    {
      arb_init(big_zeta);
      arb_init(big_zetad);
      arb_set_ui(big_zeta,1);
      arb_set_d(big_zetad,BIG_ZETAD);
      arb_union(big_zetad,big_zetad,big_zeta,prec);
      arb_set_ui(big_zeta,BIG_ZETA);
      init=true;
    }

  // z is pure scratch: compute x-big_zeta, test its sign, discard. Keep it a
  // plain local (init+clear each call) so its heap mantissa is not leaked.
  arb_t z;
  arb_init(z);
  arb_sub(z,x,big_zeta,prec);
  bool very_large_x = arb_is_positive(z);
  arb_clear(z);
  if(very_large_x) // very large x
    {
      arb_set(res,big_zetad);
      return;
    }
  arb_zeta(res,x,prec);
}
      
// Lemma 5 of ARB's g.pdf
void F_hat_twiddle_error(arb_t res, arb_t x, Lfunc *L)
{
  int64_t prec=L->wprec;
  arb_t tmp1,tmp2,tmp3,u,X,half_log_N,twoXbyr;
  arb_init(tmp1);
  arb_init(tmp2);
  arb_init(tmp3);
  arb_init(u);
  arb_init(X);
  arb_init(half_log_N);
  arb_init(twoXbyr);
  if(verbose){printf("x=");arb_printd(x,10);printf("\n");}
  arb_log_ui(half_log_N,L->conductor,prec);
  arb_mul_2exp_si(half_log_N,half_log_N,-1);
  // u = x-1/2 log N
  arb_sub(u,x,half_log_N,prec);
  if(verbose){printf("u=");arb_printd(u,10);printf("\n");}
  // X = Pi r exp(2u/r)
  arb_div_ui(tmp1,u,L->degree,prec);
  arb_mul_2exp_si(tmp1,tmp1,1);
  arb_exp(tmp2,tmp1,prec);
  arb_mul_ui(tmp1,tmp2,L->degree,prec);
  arb_mul(X,tmp1,L->pi,prec);
  // check X>r/2
  arb_set_d(tmp1,L->degree/2.0);
  arb_sub(tmp2,X,tmp1,prec);
  if(verbose){printf("X=");arb_printd(X,10);printf("\n");}
  if(!arb_is_positive(tmp2))
  {
    fprintf(stderr,"Fhat twiddle error failed X>r/2 test. Exiting.\n");
    exit(0);
  }
  // 2X/r
  arb_div_ui(twoXbyr,X,L->degree,prec);
  arb_mul_2exp_si(twoXbyr,twoXbyr,1);
  // arb_zeta won't do x that straddle an int. Will happen when
  // twoXbyr is large
  djp_zeta(tmp2,twoXbyr,prec);
  arb_pow_ui(tmp1,tmp2,L->degree,prec);
  arb_set_d(tmp2,0.5);
  arb_sub(tmp3,tmp2,twoXbyr,prec);
  arb_mul(tmp2,tmp3,L->arb_A,prec);
  arb_mul(tmp3,tmp2,L->pi,prec);
  arb_mul_2exp_si(tmp3,tmp3,1); // 2piA(1/2-2X/r)
  arb_exp(tmp2,tmp3,prec);
  arb_sub_ui(tmp3,tmp2,1,prec);
  arb_neg(tmp3,tmp3);
  arb_div(tmp2,tmp1,tmp3,prec);
  arb_mul(tmp1,u,L->nu,prec);
  arb_sub(tmp3,tmp1,X,prec);
  arb_exp(tmp1,tmp3,prec);
  arb_mul(tmp3,tmp1,tmp2,prec);
  arb_sqrt_ui(tmp1,2,prec);
  arb_pow_ui(tmp2,tmp1,L->degree,prec);
  arb_mul(res,tmp2,tmp3,prec);
  for(uint64_t j=0;j<L->degree;j++)
  {
    arb_mul_ui(tmp1,L->nus[j],L->degree,prec);
    arb_div(tmp2,tmp1,X,prec);
    arb_add_ui(tmp1,tmp2,1,prec);
    arb_pow(tmp2,tmp1,L->nus[j],prec);
    arb_mul(res,res,tmp2,prec);
  }
  //printf("F_hat_twiddle_error returning ");arb_printd(res,10);printf("\n");
  arb_clear(tmp1);
  arb_clear(tmp2);
  arb_clear(tmp3);
  arb_clear(u);
  arb_clear(X);
  arb_clear(half_log_N);
  arb_clear(twoXbyr);

}

// argument function that copes with interval starddling -ve real axis
void acb_reasonable_arg(arb_t out, acb_t in, arb_t pi, uint64_t prec)
{
  acb_t tmp;
  acb_init(tmp);

  if(arb_is_negative(acb_realref(in)) && arb_contains_zero(acb_imagref(in))) 
  {
    acb_neg(tmp, in);
    acb_arg(out, tmp, prec);
    arb_add(out,out,pi,prec);
  }
  else 
    acb_arg(out, in, prec);
  acb_clear(tmp);
}


// compute sqrt_sign from x and y such that x*sqrt_sign = conj(y*sqrt_sign)  
void fix_epsilon(acb_t res, acb_t x, acb_t y, arb_t pi, uint64_t prec)
{
  arb_t th1,th2;
  arb_t ep,two_pi;
  arb_init(two_pi);
  arb_mul_2exp_si(two_pi,pi,1);
  arb_init(th1);
  arb_init(th2);
  arb_init(ep);

  acb_reasonable_arg(th1,x,pi,prec);
  if(arb_is_negative(th1))
    arb_add(th1,th1,two_pi,prec);
  acb_reasonable_arg(th2,y,pi,prec);
  if(arb_is_negative(th2))
    arb_add(th2,th2,two_pi,prec);

  if(arb_gt(th1,th2))
  {
    arb_sub(ep,th1,th2,prec);
    arb_mul_2exp_si(ep,ep,-1);
    arb_add(ep,ep,th2,prec);
  }
  else
  {
    arb_sub(ep,th2,th1,prec);
    arb_mul_2exp_si(ep,ep,-1);
    arb_add(ep,ep,th1,prec);
  }
  arb_neg(ep,ep);

  arb_sin_cos(acb_imagref(res),acb_realref(res),ep,prec);
  arb_clear(two_pi);
  arb_clear(th1);
  arb_clear(th2);
  arb_clear(ep);

}


Lerror_t do_pre_iFFT_errors(Lfunc *L)
{
  int64_t prec=L->wprec;
  arb_t err,x,fhattwiddle,two_pi_A,tmp,th,th1,th2;
  acb_t ctmp1,ctmp2;
  arb_init(th);
  arb_init(th1);
  arb_init(th2);
  arb_init(fhattwiddle);
  arb_init(err);arb_init(x);
  arb_init(two_pi_A);
  arb_mul(two_pi_A,L->pi,L->arb_A,prec);
  arb_mul_2exp_si(two_pi_A,two_pi_A,1);
  arb_init(tmp);
  acb_init(ctmp1);
  acb_init(ctmp2);
  uint64_t i,n,j;
  uint64_t cutoff = L->M;
  const double grid_step = 2.0*M_PI*L->one_over_B;
  int64_t last_index = INT64_MIN;
  for(uint64_t m=1; m<=L->M; m++)
  {
    int64_t index = calc_m(m, grid_step, L->dc);
    if(index < last_index)
      return ERR_M_ERROR;
    last_index = index;
    // The bucket is proposed in double precision. Certify that its exact
    // displacement lies inside the Taylor remainder's half-grid radius.
    arb_mul_ui(tmp, L->one_over_root_N, m, prec);
    arb_log(th, tmp, prec);
    arb_mul_si(tmp, L->two_pi_by_B, index, prec);
    arb_sub(th1, th, tmp, prec);
    arb_abs(th1, th1);
    arb_mul_2exp_si(th2, L->two_pi_by_B, -1);
    if(!arb_le(th1, th2))
      return ERR_M_ERROR;
  }
  arb_zero(x);
  for(i=0;i<=L->fft_N/2;i++)
  {
    // Match the actual rounded index used by the convolution, rather than
    // an approximate inverse of that index. The included coefficients form
    // a prefix; this loop makes at most L->M decrements over the whole grid.
    while(cutoff && calc_m(cutoff, grid_step, L->dc)
                          + (int64_t)i > L->hi_i)
      cutoff--;
    if(cutoff==0)
      break;
    if(!M_error_at_cutoff(err,x,L,cutoff,prec))
      break; // bound this and subsequent whole Fourier terms below instead
    if(verbose&&((i%(L->fft_N/32))==0))
    {
      printf("M Error for n=%" PRIu64 " is ",i);
      arb_printd(err,10);
      printf("\n");
    }
    arb_add_error(acb_realref(L->res[i]),err);
    arb_add_error(acb_imagref(L->res[i]),err);
    arb_add_error(acb_realref(L->res[L->fft_N-1]),err);
    arb_add_error(acb_imagref(L->res[L->fft_N-1]),err);
    arb_add(x,x,L->two_pi_by_B,prec);
  }

  // error sum k\neq 0 F_hat(x+2\pi k A)
  arb_sub(err,two_pi_A,x,prec);
  F_hat_twiddle_error(fhattwiddle,err,L);
  arb_mul_2exp_si(fhattwiddle,fhattwiddle,1);
  if(verbose)
  {
    printf("F_hat_twiddle error at n=%" PRIu64 " = ",i);
    arb_printd(fhattwiddle,10);
    printf("\n");
  }

  arb_mul(err,L->eq59,L->sum_ans,prec);
  if(verbose){printf("Adding eq 5-9 error = ");arb_printd(err,10);printf("\n");}
  for(j=0;j<i;j++)
  {
    arb_add_error(acb_realref(L->res[j]),fhattwiddle);
    arb_add_error(acb_imagref(L->res[j]),fhattwiddle);
    arb_add_error(acb_realref(L->res[j]),err);
    arb_add_error(acb_imagref(L->res[j]),err);
  }

  // figure out epsilon
  acb_abs(tmp,L->res[0],prec);
  arb_set_d(th,1e-20);
  arb_sub(tmp,tmp,th,prec);

  if(!arb_is_positive(tmp)) // F_hat[0] too small to be useful
  {
    // first set the error terms for F_hat[-1]
    arb_add_error(acb_realref(L->res[L->fft_N-1]),fhattwiddle);
    arb_add_error(acb_imagref(L->res[L->fft_N-1]),fhattwiddle);
    arb_add_error(acb_realref(L->res[L->fft_N-1]),err);
    arb_add_error(acb_imagref(L->res[L->fft_N-1]),err);
    if(verbose)
    {
      printf("F_hat[1]=");acb_printd(L->res[1],30);printf("\n");
      printf("F_hat[-1]=");acb_printd(L->res[L->fft_N-1],30);printf("\n");
      printf("Using F_hat[1] and F_hat[-1] to fix epsilon.\n");
    }
    fix_epsilon(L->sqrt_sign,L->res[1],L->res[L->fft_N-1],L->pi,prec);
  }
  else // can use F_hat[0]
  {  
    if(verbose)
    {
      printf("Using F_hat[0] to fix epsilon.\n");
      printf("F_hat[0]=");acb_printd(L->res[0],10);printf("\n");
    }
    if((arb_is_negative(acb_realref(L->res[0])))&&(arb_contains_zero(acb_imagref(L->res[0]))))
    {
      if(verbose) printf("F_hat[0] straddles -ve real axis.\n");
      acb_neg(ctmp1,L->res[0]);
      acb_arg(th,ctmp1,prec);
      if(verbose) {printf("arg=");arb_printd(th,20);printf("\n");}
      arb_sin_cos(th1,th2,th,prec);
      arb_set(acb_realref(L->sqrt_sign),th2);
      arb_set(acb_imagref(L->sqrt_sign),th1);
    }
    else
    {
      if(verbose) printf("F_hat[0] doesn't straddle -ve real axis.\n");
      acb_arg(th,L->res[0],prec);
      //acb_arg(th2,L->res[L->fft_N/2],prec);
      //arb_set(th,th1); //arb_intersect(th,th1,th2,prec);
      arb_neg(th,th);
      arb_sin_cos(th1,th2,th,prec);
      arb_set(acb_realref(L->sqrt_sign),th2);
      arb_set(acb_imagref(L->sqrt_sign),th1);
      //fix_epsilon(ctmp1,L->res[1],L->res[L->fft_N-1],prec);
      //printf("Or we could have used ");acb_printd(ctmp1,10);printf("\n");
    }
  }
  //arb_one(acb_imagref(L->epsilon));arb_zero(acb_realref(L->epsilon));
  if(verbose)
  {
    printf("sqrt_sign set to ");
    acb_printd(L->sqrt_sign,10);
    printf("\n");
  }
  acb_sqr(L->sign,L->sqrt_sign,prec);
  for(n=0;n<i;n++)
    acb_mul(L->res[n],L->res[n],L->sqrt_sign,prec);

  // the rest of the vector we just approximate with F_hat_twiddle error

  // error sum F_hat(x+2\pi k A)
  F_hat_twiddle_error(fhattwiddle,x,L);
  // we have sum k\geq 0 F_hat(x+2\pi k A)
  // since x<\pi A we can just double this
  arb_mul_2exp_si(fhattwiddle,fhattwiddle,1);
  if(verbose) {
    printf("F_hat_twiddle error beyond n=%" PRIu64 " (x = ",i);
    arb_printd(x,10);
    printf(" ) = ");
    arb_printd(fhattwiddle,10);
    printf("\n");
  }
  for(;i<=L->fft_NN/2;i++)
  {
    acb_zero(L->res[i]);
    arb_add_error(acb_realref(L->res[i]),fhattwiddle);
    arb_add_error(acb_imagref(L->res[i]),fhattwiddle);
  }
  for(uint64_t n=L->fft_NN/2+1;n<L->fft_NN;n++)
    acb_conj(L->res[n],L->res[L->fft_NN-n]);


  arb_clear(th);
  arb_clear(th1);
  arb_clear(th2);
  arb_clear(fhattwiddle);
  arb_clear(err);arb_clear(x);
  arb_clear(two_pi_A);
  arb_clear(tmp);
  acb_clear(ctmp1);
  acb_clear(ctmp2);

  return ERR_SUCCESS;

}

#ifdef __cplusplus
}
#endif
