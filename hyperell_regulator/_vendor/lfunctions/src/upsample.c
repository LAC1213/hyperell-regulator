#include "glfunc.h"
#include "glfunc_internals.h"

// upsampling parameters from upsample_error_final.gp on lmfdb4 code/generic

/*
   uint64_t Ns[MAX_DEGREE-1]={1400,249,296,342,386,430,473,515};
   uint64_t Hs[MAX_DEGREE-1]={251658,31947,27352,24297,22077,20373,19016,17899};
   double errs[MAX_DEGREE-1]={5.5093645e-69,7.547808e-61,1.043349e-60,7.601862e-61,9.261187e-61,8.876395e-61,9.698565e-61,1.126559e-60};
   double D_errs[MAX_DEGREE-1]={8.354e-45,1.739e-49,2.025e-55,3.861e-62,1.889e-68,2.322e-75,5.289e-82,2.243e-88};
   */

#ifdef __cplusplus
extern "C"{
#endif

  int check_M(double M)
  {
    return M>0.5;
  }

  int check_H(double H, double h, double A, double *mus, uint64_t r) {
    double s=0.0;
    for(uint64_t d=0; d < r; d++)
      s+=mus[d];
    s+=(double)r/2.0;
    s/=M_PI;
    s=sqrt(s)*A*h/2.0;
    return H>s;
  }

  /*
     double zeta(double s)
     {
     const double precision = 1e-7;
     double output = 0.0;
     double calculation = 1.0;

     for (int denom = 2; abs(output - calculation) > precision; ++denom)
     {
     output += calculation;
     calculation = 1.0 / pow(denom, s);
     }

     return output;
     }

#define M_LN2l      0xb.17217f7d1cf79acp-4L
#define M_PIl       0xc.90fdaa22168c235p-2L

long double E(long double s, long double t, long double h, double *mus, uint64_t r, uint64_t N, long double T)
{
long double res=1.0;
for(uint64_t j=0;j<r;j++)
res*=expl((s+(long double)mus[j]-0.5l)/2.0l*logl(sqrtl((s+(long double)mus[j]-0.5l)*(s+(long double)mus[j]-0.5l)+(T+t)*(T+t))/(2.0l*M_PIl)));
res*=expl(-M_PIl*t*t/(h*h));
res*=expl(s/2.0l*logl(N));
res*=expl((long double)r/2.0l*M_LN2l);
res*=expl((long double)r*logl(zeta(s+0.5l)));
return res;
}
*/

int check_l(double M, double H, double A, uint64_t l)
{
  double tmp=M*A;
  if(H<tmp) tmp=H; // min(H,MA)
  return l<M_PI*tmp;
}
/*
// upsample error.
double upsample_error(long double M, long double H, long double h, long double A, double *mus, uint64_t r, uint64_t N, long double T, long double imz, uint64_t l)
{
if(!check_M(M)) return 0.0;
if(!check_H(H,h,A,mus,r)) return 0.0;
if(!check_l(M,H,A,l)) return 0.0;
uint64_t r0=0;
for(uint64_t d=0;d<r;d++)
if(mus[d]<0.5) r0++;
long double delta=0.0;
for(uint64_t j=0;j<r;j++)
delta+=(M+(long double)mus[j]-0.5)/((M+(long double)mus[j]-0.5)*(M+(long double)mus[j]-0.5)+T*T);
delta*=h*h;
delta/=4*M_PIl;
if(delta>=1.0) return 0.0;
long double E0=E(M,0,h,mus,r,N,T);
long double E1=E(1,H/A,h,mus,r,N,T);
long double E2=E(1,(H+1)/A,h,mus,r,N,T);
return 2.0l*coshl(M_PIl*A*imz)*expl((long double)l*logl(M_PIl*A))*(A*h*expl(M_PIl/(h*h)*(M*M+(delta*T)*(delta*T)/(long double)1-delta)-M_PIl*M*A)/((M_PIl*M*A-(long double)l)*sqrtl((long double)1-delta))*E0+expl(0.75l*(long double)r0*logl(3.0l))*E1*E1/((M_PIl*H-(long double)l)*(E1-E2)));
}
*/

void delta_term(arb_t delta,double M,double T, double mu, int64_t prec)
{
  arb_t tmp,tmp1,tmp2,T2;
  arb_init(tmp);
  arb_init(tmp1);
  arb_init(tmp2);
  arb_init(T2);
  arb_set_d(T2,T*T);
  arb_set_d(tmp,M-0.5); // M-0.5
  arb_set_d(tmp1,mu);
  arb_add(tmp2,tmp,tmp1,prec); // mu+M-0.5
  arb_mul(tmp,tmp2,tmp2,prec);
  arb_add(tmp1,tmp,T2,prec); //(mu+M-0.5)^2+T^2
  arb_div(tmp,tmp2,tmp1,prec);
  arb_add(delta,delta,tmp,prec);
  arb_clear(tmp);
  arb_clear(tmp1);
  arb_clear(tmp2);
  arb_clear(T2);
}

void arb_E(arb_t res, double sigma, arb_t t, uint64_t r, uint64_t N, double *mus, double T, arb_t pi, double h, int64_t prec)
{
  acb_t ctmp;
  arb_t tmp,tmp1,tmp2,tmp3;
  acb_init(ctmp);
  arb_init(tmp);
  arb_init(tmp1);
  arb_init(tmp2);
  arb_init(tmp3);
  arb_zero(tmp3); // noop?
  for(uint64_t j=0;j<r;j++)
  {
    arb_set_d(acb_realref(ctmp),sigma+mus[j]-0.5);
    arb_set(tmp2,acb_realref(ctmp));
    arb_set_d(acb_imagref(ctmp),T);
    arb_add(acb_imagref(ctmp),acb_imagref(ctmp),t,prec);
    acb_abs(tmp,ctmp,prec);
    arb_div(tmp1,tmp,pi,prec);
    arb_mul_2exp_si(tmp1,tmp1,-1);
    arb_log(tmp,tmp1,prec);
    arb_mul(tmp1,tmp,tmp2,prec);
    arb_add(tmp3,tmp3,tmp1,prec);
  } // tmp3 = log prod
  arb_set_d(tmp,h);
  arb_div(tmp1,t,tmp,prec);
  arb_mul(tmp,tmp1,tmp1,prec);
  arb_mul(tmp1,tmp,pi,prec); // pi t^2/h^2
  arb_sub(tmp,tmp3,tmp1,prec); // log exp()
  arb_log_ui(tmp1,N,prec);
  arb_set_d(tmp2,sigma/2.0);
  arb_mul(tmp3,tmp1,tmp2,prec); // log N^sig/2
  arb_add(tmp1,tmp,tmp3,prec); // log [N^ exp() prod]
  arb_set_d(tmp,sigma+0.5);
  arb_zeta(tmp2,tmp,prec);
  arb_log(tmp,tmp2,prec);
  arb_mul_ui(tmp2,tmp,r,prec);
  arb_add(tmp,tmp1,tmp2,prec); // log [zeta()^r N^r/2 exp() prod]
  arb_set_d(tmp2,(double) r/2.0);
  arb_log_ui(tmp1,2,prec);
  arb_mul(tmp3,tmp1,tmp2,prec);
  arb_add(tmp1,tmp,tmp3,prec);
  arb_exp(res,tmp1,prec);
  acb_clear(ctmp);
  arb_clear(tmp);
  arb_clear(tmp1);
  arb_clear(tmp2);
  arb_clear(tmp3);
}

Lerror_t arb_upsampling_error(
    arb_t res,
    double M,
    double H,
    double h,
    double A,
    double *mus,
    uint64_t r,
    uint64_t N,
    double T,
    arb_t imz,
    uint64_t l,
    arb_t pi,
    int64_t prec) {
  Lerror_t ecode=ERR_SUCCESS;
  if((!check_M(M))||(!check_H(H,h,A,mus,r))||(!check_l(M,H,A,l)))
  {
    arb_zero(res);
    return ecode|ERR_SPEC_VALUE;
  }
  uint64_t r0=0;
  for(uint64_t d=0;d<r;d++)
    if(mus[d]<0.5) r0++;
  arb_t delta,tmp,tmp1,tmp2,tmp3;
  arb_init(delta);
  arb_init(tmp);
  arb_init(tmp1);
  arb_init(tmp2);
  arb_init(tmp3);
  for(uint64_t j=0;j<r;j++)
    delta_term(delta,M,T,mus[j],prec);
  arb_set_d(tmp,h*h);
  arb_mul(delta,delta,tmp,prec);
  arb_div(delta,delta,pi,prec);
  arb_mul_2exp_si(delta,delta,-2);
  arb_sub_ui(tmp,delta,1,prec);
  if(!arb_is_negative(tmp))
    return ecode|ERR_SPEC_VALUE;

  arb_t E0,E1,E2;
  arb_init(E0);
  arb_init(E1);
  arb_init(E2);

  arb_zero(tmp);
  arb_E(E0,M,tmp,r,N,mus,T,pi,h,prec);
  arb_set_d(tmp,H);
  arb_set_d(tmp1,A);
  arb_div(tmp3,tmp,tmp1,prec); // H/A
  arb_E(E1,1.0,tmp3,r,N,mus,T,pi,h,prec);
  arb_add_ui(tmp,tmp,1,prec);
  arb_div(tmp3,tmp,tmp1,prec); // H/A
  arb_E(E2,1.0,tmp3,r,N,mus,T,pi,h,prec);

  arb_sub(tmp,E1,E2,prec);
  arb_mul(tmp1,E1,E1,prec);
  arb_div(tmp2,tmp1,tmp,prec); // E1^2/(E1-E2)

  arb_log_ui(tmp,3,prec);
  arb_set_d(tmp1,(double)r0*0.75);
  arb_mul(tmp3,tmp,tmp1,prec);
  arb_exp(tmp1,tmp3,prec); // 3^(3r0/4)
  arb_mul(tmp,tmp1,tmp2,prec); // 3^(3r0/4)*E1^2/(E1-E2)
  arb_set_d(tmp1,H);
  arb_mul(tmp2,tmp1,pi,prec);
  arb_sub_ui(tmp1,tmp2,l,prec);
  arb_t term2;
  arb_init(term2);
  arb_div(term2,tmp,tmp1,prec);

  arb_set_ui(tmp1,1);
  arb_t one_del;
  arb_init(one_del);
  arb_sub(one_del,tmp1,delta,prec);
  arb_sqrt(tmp1,one_del,prec);
  arb_set_d(tmp,M*A);
  arb_t piMA;
  arb_init(piMA);
  arb_mul(piMA,tmp,pi,prec);
  arb_sub_ui(tmp,piMA,l,prec);
  arb_mul(tmp3,tmp1,tmp,prec); // (piMA-l)sqrt(1-delta)

  arb_set_d(tmp,T*T);
  arb_mul(tmp2,tmp,delta,prec);
  arb_mul(tmp,tmp2,delta,prec);
  arb_div(tmp2,tmp,one_del,prec); // del^2T^2/(1-del)
  arb_set_d(tmp,M*M);
  arb_add(tmp1,tmp,tmp2,prec);
  arb_mul(tmp,tmp1,pi,prec);
  arb_set_d(tmp1,h*h);
  arb_div(tmp2,tmp,tmp1,prec);
  arb_sub(tmp,tmp2,piMA,prec);
  arb_exp(tmp1,tmp,prec);
  arb_div(tmp,tmp1,tmp3,prec);
  arb_mul(tmp1,tmp,E0,prec);
  arb_set_d(tmp,A*h);
  arb_t term1;
  arb_init(term1);
  arb_mul(term1,tmp,tmp1,prec);

  arb_add(tmp,term1,term2,prec);

  arb_set_d(tmp1,A);
  arb_mul(tmp1,tmp1,imz,prec);
  arb_mul(tmp2,tmp1,pi,prec);
  arb_cosh(tmp1,tmp2,prec);
  arb_mul(tmp3,tmp,tmp1,prec);

  arb_set_d(tmp,A);
  arb_mul(tmp1,tmp,pi,prec);
  arb_log(tmp,tmp1,prec);
  arb_mul_ui(tmp1,tmp,l,prec);
  arb_exp(tmp,tmp1,prec);
  arb_mul_2exp_si(tmp,tmp,1);

  arb_mul(res,tmp3,tmp,prec);



  arb_clear(term1);
  arb_clear(term2);
  arb_clear(one_del);
  arb_clear(piMA);
  arb_clear(E0);
  arb_clear(E1);
  arb_clear(E2);
  arb_clear(tmp);
  arb_clear(tmp1);
  arb_clear(tmp2);
  arb_clear(tmp3);
  arb_clear(delta);
  return ecode;
}

Lerror_t init_upsampling(Lfunc *L)
{
  Lerror_t ecode=ERR_SUCCESS;
  int64_t prec=L->wprec;

  // we allow for sampling a less than every point
  L->u_stride=1; // but actually sample at every point!
  arb_init(L->u_pi_A);
  arb_init(L->u_A);
  arb_mul_ui(L->u_A,L->arb_A,L->u_stride,prec);
  arb_init(L->u_one_over_A);
  arb_inv(L->u_one_over_A,L->u_A,prec);
  if(verbose){printf("A for upsampling set to ");arb_printd(L->u_A,10);printf("\n");}
  arb_mul(L->u_pi_A,L->pi,L->u_A,prec);

  double T=512.0/(double)L->degree/(double)OUTPUT_RATIO;
  double A=L->A*L->u_stride;
  double h=sqrt(1.0/A)*1.001;
  double H=ceil(A*A*h*h/2.0);
  double M=H/A;
  arb_init(L->upsampling_error);
  arb_t tmp,zero;
  arb_init(tmp);arb_init(zero);
  while(true)
  {
    ecode|=arb_upsampling_error(L->upsampling_error,M,H,h,A,L->mus,L->degree,L->conductor,T,zero,MAX_L,L->pi,L->wprec);
    if(fatal_error(ecode))
    {
      arb_clear(tmp);
      arb_clear(zero);
      return ecode;
    }
    arb_mul_2exp_si(tmp,L->upsampling_error,L->target_prec);
    arb_sub_ui(tmp,tmp,1,prec);
    if(arb_is_negative(tmp))
      break;
    h*=1.01;
    H=ceil(A*A*h*h/2.0);
    M=H/A;
  }
  arb_clear(tmp);
  arb_clear(zero);
  if(verbose){printf("Upsampling error set to ");arb_printd(L->upsampling_error,10);printf("\n");}


  arb_init(L->u_H);
  arb_set_d(L->u_H, h);
  L->u_N=H;
  if(verbose){printf("Upsampling H set to ");arb_printd(L->u_H,20);printf("\n");}
  arb_init(L->u_pi_by_H2); // -pi/H^2
  arb_inv(L->u_pi_by_H2,L->u_H,prec); // 1/H
  arb_mul(L->u_pi_by_H2,L->u_pi_by_H2,L->u_pi_by_H2,prec); // 1/H^2
  arb_mul(L->u_pi_by_H2,L->u_pi_by_H2,L->pi,prec);
  arb_neg(L->u_pi_by_H2,L->u_pi_by_H2); // -pi/H^2
  if (verbose){
    printf("-pi/H^2 = ");
    arb_printd(L->u_pi_by_H2,20);
    printf("\n");
  }

  if (verbose)
    printf("Upsampling N set to %" PRIu64 "\n",L->u_N);

  L->u_no_values=L->fft_NN/OUTPUT_RATIO+L->fft_NN/TURING_RATIO+L->u_N*4*L->u_stride+1;
  L->u_values[0]=(arb_t *)malloc(sizeof(arb_t)*L->u_no_values);
  L->u_values[1]=(arb_t *)malloc(sizeof(arb_t)*L->u_no_values);
  L->u_no_values_off=L->u_no_values-L->u_N*L->u_stride*2;
  L->u_values_off[0]=L->u_values[0]+L->u_N*L->u_stride*2;
  L->u_values_off[1]=L->u_values[1]+L->u_N*L->u_stride*2;
  for(uint64_t n=0;n<L->u_no_values;n++)
  {
    arb_init(L->u_values[0][n]);
    arb_init(L->u_values[1][n]);
  }


  return ERR_SUCCESS;
}

int64_t left_n(arb_ptr diff, arb_ptr t0, arb_t A, uint64_t prec)
{
  static arb_t tmp,tmp1;
  static fmpz_t fmpz_tmp;
  static bool init=false;
  if(!init)
  {
    init=true;
    arb_init(tmp);
    arb_init(tmp1);
    fmpz_init(fmpz_tmp);
  }
  arb_mul(tmp,t0,A,prec);
  arb_get_mid_arb(tmp1,tmp);
  arb_floor(tmp,tmp1,prec);
  if(arb_get_unique_fmpz(fmpz_tmp,tmp)==0)
  {
    //printf("t0 = ");arb_printd(t0,20);printf(" tmp = ");arb_printd(tmp,20);printf("\n");
    //printf("Error rounding to int in upsample routines.\n");
    return(BAD_64);
  }
  int64_t res=fmpz_get_si(fmpz_tmp); // this is going to be 1st point in left tail
  arb_set_si(tmp,res);
  arb_div(tmp,tmp,A,prec);
  arb_sub(diff,tmp,t0,prec); // will be -ve n pts to left of t0
  return(res);
}

// sinc x = sin(pi x)/(pi x)
// on entry sin_x = sin(pi x)
void sinc(arb_t res, arb_t x, arb_t sin_x, arb_t pi, uint64_t prec)
{
  static arb_t tmp1;
  static bool init=false;
  if(!init)
  {
    init=true;
    arb_init(tmp1);
  }
  arb_mul(tmp1,x,pi,prec);
  arb_div(res,sin_x,tmp1,prec);
}

//
void do_point(arb_ptr res, Lfunc *L, int64_t n, arb_t t, arb_t delta, arb_t sin_delta, arb_t A, uint64_t side, uint64_t prec)
{
  static arb_t tmp,tmp1,tmp2;
  static bool init=false;
  if(!init)
  {
    init=true;
    arb_init(tmp);
    arb_init(tmp1);
    arb_init(tmp2);
  }
  arb_mul(tmp,delta,delta,prec);
  arb_mul(tmp1,tmp,L->u_pi_by_H2,prec);
  arb_mul_ui(tmp,L->pi,L->degree,prec);
  arb_mul(tmp2,tmp,t,prec);
  arb_mul_2exp_si(tmp2,tmp2,-2); // pi r t/4
  arb_add(tmp,tmp1,tmp2,prec);
  arb_exp(tmp1,tmp,prec);
  arb_mul(tmp2,delta,A,prec);
  sinc(tmp,tmp2,sin_delta,L->pi,prec);
  arb_mul(tmp2,tmp1,tmp,prec);
  arb_mul(res,tmp2,L->u_values[side][n],prec);
}

// estimate f(t0) by upsampling off Lu->values
bool upsample_stride(arb_ptr res, arb_ptr t0, Lfunc *L, uint64_t side, uint64_t prec)
{
  static arb_t step,diff,this_diff,term,A,t,t1,t_delta,sin_diff,neg_sin_diff;
  static bool init=false;
  if(!init)
  {
    init=true;
    arb_init(step);
    arb_init(t);
    arb_init(t1);
    arb_init(t_delta);
    arb_init(A);
    arb_init(diff);
    arb_init(this_diff);
    arb_init(term);
    arb_init(sin_diff);
    arb_init(neg_sin_diff);
  }
  arb_mul_ui(step,L->one_over_A,L->u_stride,prec);
  arb_div_ui(A,L->arb_A,L->u_stride,prec);

  int64_t n=left_n(diff,t0,L->arb_A,prec);
  if(n==BAD_64)
    return(false);
  arb_mul(neg_sin_diff,diff,L->u_pi_A,prec);
  arb_sin(sin_diff,neg_sin_diff,prec);
  arb_neg(neg_sin_diff,sin_diff);
  //printf("Nearest n = %" PRId64 "\n",n);
  arb_mul_ui(t_delta,L->one_over_A,L->u_stride,prec);
  arb_mul_si(t1,L->one_over_A,n,prec);
  arb_set(t,t1);
  int64_t nn=n+L->u_N*L->u_stride*2,nn1=nn; // offset into values

  // do nearest point
  do_point(res,L,nn,t,diff,sin_diff,A,side,prec);

  arb_set(this_diff,diff);
  // do Lu->N-1 points to left of left
  for(uint64_t count=0;count<L->u_N-1;count++)
  {
    arb_sub(this_diff,this_diff,step,prec);
    arb_sub(t,t,t_delta,prec);
    nn-=L->u_stride;
    if(count&1)
      do_point(term,L,nn,t,this_diff,sin_diff,A,side,prec);
    else
      do_point(term,L,nn,t,this_diff,neg_sin_diff,A,side,prec);

    arb_add(res,res,term,prec);
  }

  arb_set(t,t1); // point to left again
  arb_set(this_diff,diff);
  nn=nn1;

  // do N points to right of left
  for(uint64_t count=0;count<L->u_N;count++)
  {
    arb_add(this_diff,this_diff,step,prec);
    arb_add(t,t,t_delta,prec);
    nn+=L->u_stride;
    if(count&1)
      do_point(term,L,nn,t,this_diff,sin_diff,A,side,prec);
    else
      do_point(term,L,nn,t,this_diff,neg_sin_diff,A,side,prec);
    arb_add(res,res,term,prec);
  }
  arb_add_error(res,L->upsampling_error);
  return(true);
}


#define STACK_SIZE (32)
// estimate f(t0) by upsampling off Lu->values
bool new_upsample_stride(arb_ptr res, arb_ptr t0, Lfunc *L, uint64_t side, uint64_t prec)
{
  static arb_t step,diff,this_diff,term,A,t,t1,t_delta,sin_diff,neg_sin_diff,res_stack[STACK_SIZE];
  static bool init=false;
  if(!init)
  {
    init=true;
    for(uint64_t i=0;i<STACK_SIZE;i++)
      arb_init(res_stack[i]);
    arb_init(step);
    arb_init(t);
    arb_init(t1);
    arb_init(t_delta);
    arb_init(A);
    arb_init(diff);
    arb_init(this_diff);
    arb_init(term);
    arb_init(sin_diff);
    arb_init(neg_sin_diff);
  }
  arb_mul_ui(step,L->one_over_A,L->u_stride,prec);
  arb_div_ui(A,L->arb_A,L->u_stride,prec);

  int64_t n=left_n(diff,t0,L->arb_A,prec);
  if(n==BAD_64)
    return false;
  arb_mul(neg_sin_diff,diff,L->u_pi_A,prec);
  arb_sin(sin_diff,neg_sin_diff,prec);
  arb_neg(neg_sin_diff,sin_diff);
  //printf("Nearest n = %" PRId64 "\n",n);
  arb_mul_ui(t_delta,L->one_over_A,L->u_stride,prec);
  arb_mul_si(t1,L->one_over_A,n,prec);
  arb_set(t,t1);
  int64_t nn=n+L->u_N*L->u_stride*2,nn1=nn; // offset into values

  // do nearest point
  do_point(res_stack[0],L,nn,t,diff,sin_diff,A,side,prec);
  if(verbose){printf("upsample after first pt = ");arb_printd(res_stack[0],20);printf("\n");}

  uint64_t res_ptr=1,res_count=1;

  arb_set(this_diff,diff);
  // do Lu->N-1 points to left of left
  for(uint64_t count=0;count<L->u_N-1;count++)
  {
    if(res_ptr>=STACK_SIZE)
      return false;
    arb_sub(this_diff,this_diff,step,prec);
    arb_sub(t,t,t_delta,prec);
    nn-=L->u_stride;
    if(count&1)
      do_point(res_stack[res_ptr++],L,nn,t,this_diff,sin_diff,A,side,prec);
    else
      do_point(res_stack[res_ptr++],L,nn,t,this_diff,neg_sin_diff,A,side,prec);
    res_count++;
    for(uint64_t i=res_count;!(i&1);i>>=1)
    {
      arb_add(res_stack[res_ptr-2],res_stack[res_ptr-2],res_stack[res_ptr-1],prec);
      res_ptr--;
    }
  }

  arb_set(t,t1); // point to left again
  arb_set(this_diff,diff);
  nn=nn1;

  // do N points to right of left
  for(uint64_t count=0;count<L->u_N;count++) {
    if(res_ptr>=STACK_SIZE)
      return false;
    arb_add(this_diff,this_diff,step,prec);
    arb_add(t,t,t_delta,prec);
    nn+=L->u_stride;
    if(count&1)
      do_point(res_stack[res_ptr++],L,nn,t,this_diff,sin_diff,A,side,prec);
    else
      do_point(res_stack[res_ptr++],L,nn,t,this_diff,neg_sin_diff,A,side,prec);
    res_count++;
    for(uint64_t i=res_count;!(i&1);i>>=1) {
      arb_add(res_stack[res_ptr-2],res_stack[res_ptr-2],res_stack[res_ptr-1],prec);
      res_ptr--;
    }
  }
  arb_zero(res);
  for(uint64_t i=0;i<res_ptr;i++)
    arb_add(res,res,res_stack[i],prec);
  arb_add_error(res,L->upsampling_error);
  return(true);
}

// delta=n/A-t
// sin_delta=sin(Pi*A*(n/A-t))
// compute differential at t
void do_point_f_dash(arb_ptr res, Lfunc *L, int64_t n, arb_t t, arb_t delta, arb_t sin_delta, arb_t cos_delta, uint64_t side, uint64_t prec) {
  static arb_t tmp,tmp1,tmp2,tmp3;
  static arb_t A_pi_delta;
  static bool init=false;
  if(!init)
  {
    init=true;
    arb_init(tmp);
    arb_init(tmp1);
    arb_init(tmp2);
    arb_init(tmp3);
    arb_init(A_pi_delta);
  }
  //printf("In do_point_f_dash with t= ");arb_printd(t,20);printf("\n");
  //printf("In do_point_f_dash with delta= ");arb_printd(delta,20);printf("\n");
  arb_mul(A_pi_delta,delta,L->u_pi_A,prec);
  arb_div(tmp,sin_delta,A_pi_delta,prec); // sin(A pi delta)/A pi delta
  arb_sub(tmp1,tmp,cos_delta,prec); // sin(A pi delta)/A pi delta-cos(A pi delta)
  arb_div(tmp,tmp1,delta,prec);
  arb_mul(tmp1,delta,delta,prec);
  arb_mul(tmp3,tmp1,L->u_pi_by_H2,prec);
  arb_mul_ui(tmp2,L->pi,L->degree,prec); // pi r
  arb_mul(tmp1,tmp2,t,prec); // pi r t
  arb_mul_2exp_si(tmp1,tmp1,-2); // (pi r t)/4
  arb_add(tmp2,tmp3,tmp1,prec); // (pi r t)/4 -pi(t-t0)^2/H^2
  arb_exp(tmp1,tmp2,prec);
  //printf("exp() =");arb_printd(tmp1,20);printf("\n)");
  arb_mul(tmp2,tmp1,tmp,prec);
  arb_mul(res,tmp2,L->u_values[side][n],prec);
  //printf("do_point_f_dash returning ");arb_printd(res,20);printf("\n");
}

bool f_dash(arb_ptr res, arb_ptr t0, Lfunc *L, uint64_t side, uint64_t prec)
{
  static arb_t step,diff,this_diff,term,A,t,t1,t_delta;
  //static arb_t t_error,diff_error;
  static arb_t sin_diff,neg_sin_diff,cos_diff,neg_cos_diff;
  static bool init=false;
  if(!init)
  {
    init=true;
    arb_init(step);
    arb_init(t);
    arb_init(t1);
    arb_init(t_delta);
    arb_init(A);
    arb_init(diff);
    arb_init(this_diff);
    arb_init(sin_diff);
    arb_init(neg_sin_diff);
    arb_init(cos_diff);
    arb_init(neg_cos_diff);
    arb_init(term);
    //arb_init(t_error);
    //arb_init(diff_error);
    //arb_mul_2exp_si(t_error,Lc->one_over_A,-1);
    //arb_set_d(diff_error,Lu->stride*2);
    //arb_inv(diff_error,diff_error,prec);
  }
  arb_mul_ui(step,L->one_over_A,L->u_stride,prec);
  arb_div_ui(A,L->arb_A,L->u_stride,prec);

  int64_t n=left_n(diff,t0,L->arb_A,prec);
  if(n==BAD_64)
    return(false);
  //printf("Nearest n = %" PRId64 " diff = ",n);arb_printd(diff,10);printf("\n");
  arb_mul(neg_sin_diff,diff,L->arb_A,prec);
  arb_sin_cos_pi(sin_diff,cos_diff,neg_sin_diff,prec);
  arb_neg(neg_sin_diff,sin_diff);
  arb_neg(neg_cos_diff,cos_diff);
  arb_mul_ui(t_delta,L->one_over_A,L->u_stride,prec);
  arb_mul_si(t1,L->one_over_A,n,prec);
  arb_set(t,t1);
  int64_t nn=n+L->u_N*L->u_stride*2,nn1=nn; // offset into values

  // do nearest point
  do_point_f_dash(res,L,nn,t,diff,sin_diff,cos_diff,side,prec);
  arb_set(this_diff,diff);
  // do Lu->N points to left of nearest
  for(uint64_t count=0;count<L->u_N;count++)
  {
    arb_sub(this_diff,this_diff,step,prec);
    arb_sub(t,t,t_delta,prec);
    nn-=L->u_stride;
    do_point_f_dash(term,L,nn,t,this_diff,(count&1)?sin_diff:neg_sin_diff,(count&1)?cos_diff:neg_cos_diff,side,prec);
    arb_add(res,res,term,prec);
  }

  arb_set(t,t1); // point to middle again
  arb_set(this_diff,diff);
  nn=nn1;

  for(uint64_t count=0;count<L->u_N;count++)
  {
    arb_add(this_diff,this_diff,step,prec);
    arb_add(t,t,t_delta,prec);
    nn+=L->u_stride;
    do_point_f_dash(term,L,nn,t,this_diff,(count&1)?sin_diff:neg_sin_diff,(count&1)?cos_diff:neg_cos_diff,side,prec);
    arb_add(res,res,term,prec);
  }
  // nothing rigorous about N-R
  //arb_add_error(res,Lu->upsampling_error);
  return true;
}


#define N_NEWTON_ITS (10)
// Newton iteration to isolate zeros
// using arb is overkill here as its not rigorous anyway
// the zero found will be checked rigorously later
bool newton(arb_ptr res, arb_ptr t0, Lfunc *L, uint64_t side, uint64_t prec)
{
  static bool init=false;
  static arb_t f0,fd,t,tmp1;
  if(!init)
  {
    init=true;
    arb_init(f0);
    arb_init(fd);
    arb_init(t);
    arb_init(tmp1);
  }
  arb_set(t,t0);
  for(uint64_t n=0;n<N_NEWTON_ITS;n++) {
    //printf("Newton t = ");arb_printd(t,60);printf("\n");
    if(!upsample_stride(f0,t,L,side,prec)) // f(t)
      return false;
    //printf("Newton f(t) = ");arb_printd(f0,60);printf("\n");
    if(!f_dash(fd,t,L,side,prec)) // f'(t)
      return false;
    //printf("Newton f'(t) = ");arb_printd(fd,60);printf("\n");

    arb_div(tmp1,f0,fd,prec);
    arb_sub(t,t,tmp1,prec); // t := t-f(t)/f'(t)
  }
  arb_set(res,t);
  return true;
}

#ifdef __cplusplus
}
#endif
