#include "glfunc.h"
#include "glfunc_internals.h"
#include <math.h>
#include <flint/acb_hypgeom.h>



#ifdef __cplusplus
extern "C"{
#endif
  /*
  // find the first point in the left tail for upsampling
  int64_t s_left_n(acb_ptr z, arb_t A, uint64_t prec) {
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
    arb_mul(tmp,acb_realref(z),A,prec);
    arb_get_mid_arb(tmp1,tmp);
    arb_floor(tmp,tmp1,prec);
    if(arb_get_unique_fmpz(fmpz_tmp,tmp)==0)
      return BAD_64;
    return fmpz_get_si(fmpz_tmp); // this is going to be 1st point in left tail
  }
  */

  // sinc x = sin(pi x)/(pi x)
  // on entry sin_x = sin(pi x)
  Lerror_t s_sinc(acb_t res, acb_t pi_x, acb_t sin_pi_x, int64_t prec)
  {
    (void)sin_pi_x;
    acb_sinc(res,pi_x,prec);
    return ERR_SUCCESS;
  }

  Lerror_t s_sinc_dash(acb_t res, acb_t pi_x, acb_t sin_pi_x, acb_t cos_pi_x, int64_t prec)
  {
    static bool init=false;
    static acb_t tmp1,tmp2,tmp3;
    if(!init)
      {
	init=true;
	acb_init(tmp1);
	acb_init(tmp2);
	acb_init(tmp3);
      }
    if(acb_contains_zero(pi_x))
    {
      // d/dz sinc(z) = -z/3 * 0F1(5/2; -z^2/4), an entire function.
      acb_t parameter, argument;
      acb_init(parameter); acb_init(argument);
      acb_set_ui(parameter,5); acb_mul_2exp_si(parameter,parameter,-1);
      acb_sqr(argument,pi_x,prec);
      acb_mul_2exp_si(argument,argument,-2); acb_neg(argument,argument);
      acb_hypgeom_0f1(res,parameter,argument,0,prec);
      acb_mul(res,res,pi_x,prec); acb_div_ui(res,res,3,prec); acb_neg(res,res);
      acb_clear(parameter); acb_clear(argument);
      return ERR_SUCCESS;
    }
    acb_div(tmp1,cos_pi_x,pi_x,prec);
    acb_div(tmp2,sin_pi_x,pi_x,prec);
    acb_div(tmp3,tmp2,pi_x,prec);
    acb_sub(res,tmp1,tmp3,prec);
    return ERR_SUCCESS;
  }

  void pi_r_z_by_4(acb_t res, acb_t z, Lfunc *L, int64_t prec)
  {
    static bool init=false;
    static arb_t rtmp;
    if(!init)
      {
	init=true;
	arb_init(rtmp);
      }
    arb_mul_ui(rtmp,L->pi,L->degree,prec); // pi r
    arb_mul_2exp_si(rtmp,rtmp,-2); // pi r/4
    acb_mul_arb(res,z,rtmp,prec); // pi r z /4
  }

  // -Pi(z-t0)^2/h^2
  void gaussian(arb_t res, acb_t z, arb_t pi_by_H2, int64_t prec)
  {
    static bool init=false;
    static arb_t tmp,tmp1;
    if(!init)
      {
	init=true;
	arb_init(tmp);
	arb_init(tmp1);
      }
    arb_sqr(tmp,acb_imagref(z),prec); // -(z-t0)^2
    arb_mul(res,tmp,pi_by_H2,prec); // pi(z-t0)^2/h^2
  }

  // return pi r/4 - 2 pi(z-2t0)/H^2
  void Dexp(acb_t res,Lfunc *L,acb_t z,arb_t pi_by_H2, int64_t prec)
  {
    arb_t rtmp;
    acb_t zm,ctmp;
    //printf("In Dexp with z = ");acb_printd(z,20);printf("\n");
    arb_init(rtmp);acb_init(zm);acb_init(ctmp);
    arb_mul_ui(rtmp,L->pi,L->degree,prec);
    arb_mul_2exp_si(rtmp,rtmp,-2); // Pi r/4
    arb_set(acb_imagref(zm),acb_imagref(z));
    arb_neg(acb_realref(zm),acb_realref(z)); // zm:=z-2t0
    acb_mul_arb(ctmp,zm,pi_by_H2,prec);
    acb_mul_2exp_si(ctmp,ctmp,1); // -2 pi/H^2 (z-2t0)
    acb_add_arb(res,ctmp,rtmp,prec); // -2pi/H^2 (z-2t0) + pi r /4
    arb_clear(rtmp);acb_clear(zm);acb_clear(ctmp);
  }
  
  
  // compute W(k/A)
  void W_k_A(arb_ptr res, Lfunc *L, int64_t k, int64_t prec, arb_t t0, arb_t pi_by_H2, arb_t A) {
    static bool init=false;
    static arb_t tmp,tmp1,tmp2,ka;
    if(!init)
    {
      init=true;
      arb_init(tmp);
      arb_init(tmp1);
      arb_init(tmp2);
      arb_init(ka);
    }
    //printf("t0 = ");arb_printd(t0,20);printf("\n");
    arb_set_si(tmp,k);
    // ka = // k/A
    arb_div(ka,tmp,A,prec);
    arb_sub(tmp,ka,t0,prec);
    arb_mul(tmp1,tmp,tmp,prec);
    arb_mul(tmp,tmp1,pi_by_H2,prec); // -pi (k/A-t0)^2/h^2
    arb_mul_ui(tmp1,ka,L->degree,prec);
    arb_mul(tmp2,tmp1,L->pi,prec);
    arb_mul_2exp_si(tmp2,tmp2,-2); // pi r k/(4A)
    arb_add(tmp1,tmp,tmp2,prec);
    // tmp = exp((k pi r)/(4 A) - (pi (k/A - t0)^2)/h^2)
    arb_exp(tmp,tmp1,prec);
    //printf("factor = ");arb_printd(tmp,20);printf("\n");
    if(k>=0)
      arb_mul(res,L->u_values_off[0][k],tmp,prec); // *Lambda(k/A)
    else
      arb_mul(res,L->u_values_off[1][-k],tmp,prec);
    //printf("W(%" PRId64 "/A) = W(",k);arb_printd(ka,20);printf(") = ");;arb_printd(res,20);printf("\n");
  }
  
  // W(n/A)sinc'(Pi(Az-k))
  // delta = Pi(Az-k)
  Lerror_t s_do_point_dash(acb_ptr res, Lfunc *L, int64_t n, acb_t z, acb_t delta, acb_t sin_delta, acb_t cos_delta, arb_t A, int64_t prec, arb_t pi_by_H2)
  {
    static acb_t tmp,tmp1;
    static arb_t tmp2;
    static bool init=false;
    if(!init)
    {
      init=true;
      acb_init(tmp);
      acb_init(tmp1);
      arb_init(tmp2);
    }
    Lerror_t ecode=s_sinc_dash(tmp,delta,sin_delta,cos_delta,prec);
    if(fatal_error(ecode))
      return ecode;
    W_k_A(tmp2,L,n,prec,acb_realref(z),pi_by_H2,A);
    acb_mul_arb(tmp1,tmp,tmp2,prec);
    acb_mul_arb(tmp,tmp1,A,prec);
    acb_mul_arb(res,tmp,L->pi,prec);
    return ecode;
  }

  
  // W(n/A)sinc(Pi(Az-n))
  Lerror_t s_do_point(acb_ptr res, Lfunc *L, int64_t n, acb_t z, acb_t delta, acb_t sin_delta, arb_t A, int64_t prec, arb_t pi_by_H2)
  {
    static acb_t tmp;
    static arb_t tmp2;
    static bool init=false;
    if(!init)
    {
      init=true;
      acb_init(tmp);
      arb_init(tmp2);
    }

    Lerror_t ecode=s_sinc(tmp,delta,sin_delta,prec);
    if(fatal_error(ecode))
      return ecode;
    W_k_A(tmp2,L,n,prec,acb_realref(z),pi_by_H2,A);
    acb_mul_arb(res,tmp,tmp2,prec);
    return ecode;
  }

  // estimate W(1/2+iz) by upsampling off Lu->u_values_off
  Lerror_t s_upsample_stride_dash(acb_ptr res, acb_ptr z, double t0, Lfunc *L, int64_t prec, arb_t pi_by_H2, arb_t err, int64_t N, uint64_t stride)
  {
    static arb_t A,pi; // the A for upsampling = usual A / stride
    static acb_t diff,this_diff,term,sin_diff,cos_diff;
    static bool init=false;
    if(!init)
    {
      init=true;
      arb_init(A);
      acb_init(diff);
      acb_init(this_diff);
      acb_init(term);
      acb_init(sin_diff);
      acb_init(cos_diff);
      arb_init(pi);
      arb_const_pi(pi,prec);
    }

    // sum runs for |k-t_0*A|<N
    int64_t k0=-(N-t0*L->A/stride-1);
    int64_t k1=(N+t0*L->A/stride-1);
    if(verbose)printf("Sampling for k=[%ld,%ld] and n=[%ld,%ld]\n",k0,k1,k0*stride,k1*stride);

    arb_div_ui(A,L->arb_A,stride,prec);
    if(verbose){printf("Special point upsampling A set to ");arb_printd(A,20);printf("\n");}
    acb_mul_arb(diff,z,A,prec);
    acb_sub_si(diff,diff,k0,prec);
    acb_mul_arb(diff,diff,L->pi,prec);
    acb_sin(sin_diff,diff,prec);
    acb_cos(cos_diff,diff,prec);
    acb_zero(res);
    int64_t n=k0*stride;
    for(int64_t k=k0;k<=k1;k++)
      {
	Lerror_t ecode=s_do_point_dash(term,L,n,z,diff,sin_diff,cos_diff,L->arb_A,prec,pi_by_H2);
	if(fatal_error(ecode))
	  return ecode;
	acb_add(res,res,term,prec);
	acb_neg(sin_diff,sin_diff);
	acb_neg(cos_diff,cos_diff);
	n+=stride;
	acb_sub_arb(diff,diff,pi,prec);
      }
    acb_div_ui(res,res,stride,prec); // because we use big A in (Pi A)^l
    arb_add_error(acb_realref(res),err);
    arb_add_error(acb_imagref(res),err);

    return ERR_SUCCESS;
  }

  
  // estimate W(1/2+iz) by upsampling off Lu->u_values_off
  Lerror_t s_upsample_stride(acb_ptr res, acb_ptr z, double t0, Lfunc *L, int64_t prec, arb_t pi_by_H2, arb_t err, int64_t N, uint64_t stride)
  {
    static arb_t A,pi; // the A for upsampling = usual A / stride
    static acb_t diff,this_diff,term,sin_diff,neg_sin_diff;
    static bool init=false;
    if(!init)
    {
      init=true;
      arb_init(A);
      acb_init(diff);
      acb_init(this_diff);
      acb_init(term);
      acb_init(sin_diff);
      acb_init(neg_sin_diff);
      arb_init(pi);
      arb_const_pi(pi,prec);
    }

    // sum runs for |k-t_0*A|<N
    int64_t k0=-(N-t0*L->A/stride-1);
    int64_t k1=(N+t0*L->A/stride-1);
    arb_div_ui(A,L->arb_A,stride,prec);
    acb_mul_arb(diff,z,A,prec);
    acb_sub_si(diff,diff,k0,prec);
    acb_mul_arb(diff,diff,L->pi,prec);
    acb_sin(sin_diff,diff,prec);

    acb_zero(res);
    int64_t n=k0*stride;
    for(int64_t k=k0;k<=k1;k++)
      {
	Lerror_t ecode=s_do_point(term,L,n,z,diff,sin_diff,L->arb_A,prec,pi_by_H2);
	if(verbose&&(k==0))
	  {
	    printf("W(0)sinc(stuff)  returned ");
	    acb_printd(term,20);
	    printf("\n");
	  }
	if(fatal_error(ecode))
	  return ecode;
	acb_add(res,res,term,prec);
	acb_neg(sin_diff,sin_diff);
	n+=stride;
	acb_sub_arb(diff,diff,pi,prec);
      }
    arb_add_error(acb_realref(res),err);
    arb_add_error(acb_imagref(res),err);

    return ERR_SUCCESS;
      
  }

  //1/gamma_r(s)
  void acb_rgamma_r(acb_t res, acb_t s, Lfunc *L, uint64_t prec)
  {
    acb_t s_by_2,lg,ctmp,ctmp1;
    arb_t log_pi;
    acb_init(s_by_2);
    acb_init(ctmp);
    acb_init(lg);
    acb_init(ctmp1);
    arb_init(log_pi);
    arb_log(log_pi,L->pi,prec);
    //printf("in abs_gamma_r with s = ");acb_printd(s,10);printf("\n");
    acb_mul_2exp_si(s_by_2,s,-1); // s/2

    acb_mul_arb(ctmp,s_by_2,log_pi,prec); 
    acb_exp(ctmp1,ctmp,prec);
    acb_rgamma(ctmp,s_by_2,prec);
    acb_mul(res,ctmp,ctmp1,prec);
    //printf("gamma_r returning ");acb_printd(res,20);printf("\n");
    acb_clear(s_by_2);
    acb_clear(lg);
    acb_clear(ctmp);
    acb_clear(ctmp1);
    arb_clear(log_pi);
  }


  // gamma(s) per l.pdf with epsilon and N^1/2(s-1/2)
  void spec_rgamma(acb_t res, acb_t s, Lfunc *L, int64_t prec)
  {
    acb_t tmp1,tmp2,tmp3;
    acb_init(tmp1);
    acb_init(tmp2);
    acb_init(tmp3);
    arb_t tmp;
    arb_init(tmp);

    acb_set_ui(res,1);
    for(uint64_t j=0;j<L->degree;j++) {
      arb_set_d(tmp,L->mus[j]);
      acb_add_arb(tmp1,s,tmp,prec);
      acb_rgamma_r(tmp2,tmp1,L,prec);
      acb_mul(res,res,tmp2,prec);
    }
    arb_set_d(tmp,0.5);
    acb_sub_arb(tmp1,s,tmp,prec); // s-1/2
    acb_mul_2exp_si(tmp1,tmp1,-1); // (s-1/2)/2
    arb_log_ui(tmp,L->conductor,prec);
    acb_mul_arb(tmp2,tmp1,tmp,prec);
    acb_neg(tmp2,tmp2);
    acb_exp(tmp1,tmp2,prec);
    //printf("N^(-(s-1/2)/2) = ");acb_printd(tmp1,20);printf("\n");
    acb_mul(res,res,tmp1,prec);
    acb_mul(res,res,L->sqrt_sign,prec);
    //printf("spec_gamma returning ");acb_printd(res,10);printf("\n");
    acb_clear(tmp1);
    acb_clear(tmp2);
    acb_clear(tmp3);
    arb_clear(tmp);
  }

  void Lam_to_L(acb_t res,acb_t Lam,acb_t an_s,Lfunc *L,int64_t prec)
  {
    acb_t ctmp;
    acb_init(ctmp);
    spec_rgamma(ctmp,an_s,L,prec);
    if(verbose){printf("going from Lambda to L by multiplying by ");acb_printd(ctmp,20);printf("\n");}
    acb_mul(res,Lam,ctmp,prec);
    acb_clear(ctmp);
  }


  void W_to_Lam(acb_t Lam,acb_t W,acb_t z,Lfunc *L,arb_t pi_by_H2,int64_t prec)
  {
    acb_t s,ctmp;
    arb_t tmp;
    acb_init(s);acb_init(ctmp);arb_init(tmp);
    acb_mul_arb(s,z,L->pi,prec);
    acb_mul_ui(ctmp,s,L->degree,prec);
    acb_mul_2exp_si(ctmp,ctmp,-2); // pi r z /4
    arb_mul(tmp,acb_imagref(z),acb_imagref(z),prec);
    arb_mul(tmp,tmp,pi_by_H2,prec);
    arb_sub(acb_realref(ctmp),acb_realref(ctmp),tmp,prec);
    acb_exp(s,ctmp,prec);
    if(verbose){printf("going from W to Lambda by dividing by ");acb_printd(s,20);printf("\n");}
    acb_div(Lam,W,s,prec);
    acb_clear(s);acb_clear(ctmp);arb_clear(tmp);

  }

  void Lam_dash_to_L_dash(acb_t res, acb_t Lam_dash, acb_t Lam,
                         acb_t s, Lfunc *L, int64_t prec)
  {
    /* L = R Lambda, hence L' = R (Lambda' + Lambda (log R)').
       The driver excludes gamma poles for first derivatives. */
    acb_t u, psi, dlog, scale;
    arb_t tmp, logpi;
    acb_init(u); acb_init(psi); acb_init(dlog); acb_init(scale);
    arb_init(tmp); arb_init(logpi);
    arb_log_ui(tmp, L->conductor, prec);
    arb_mul_2exp_si(tmp, tmp, -1);
    acb_set_arb(dlog, tmp); acb_neg(dlog, dlog);
    arb_log(logpi, L->pi, prec);
    for (uint64_t j=0; j<L->degree; j++) {
      arb_set_d(tmp, L->mus[j]);
      acb_add_arb(u, s, tmp, prec); acb_mul_2exp_si(u, u, -1);
      acb_digamma(psi, u, prec);
      acb_neg(psi, psi); acb_add_arb(psi, psi, logpi, prec);
      acb_mul_2exp_si(psi, psi, -1); acb_add(dlog, dlog, psi, prec);
    }
    spec_rgamma(scale, s, L, prec);
    acb_mul(u, Lam, dlog, prec); acb_add(u, u, Lam_dash, prec);
    acb_mul(res, scale, u, prec);
    acb_clear(u); acb_clear(psi); acb_clear(dlog); acb_clear(scale);
    arb_clear(tmp); arb_clear(logpi);
  }
  
    // go from W'(z) to Lam'(s)
  // W(z)=Lam(1/2+iz)exp(pi r t/4 - 
  void W_dash_to_Lam_dash(acb_t res,acb_t W_dash,acb_t W,acb_t z,Lfunc *L,
			  arb_t pi_by_H2,int64_t prec)
  {
    acb_t ctmp1,ctmp2,ctmp3,dexp,pir;
    arb_t G;
    arb_init(G);
    acb_init(pir);acb_init(ctmp1);acb_init(ctmp2);
    acb_init(ctmp3);acb_init(dexp);
    
    gaussian(G,z,pi_by_H2,prec); // Pi/h^2 (z-t0)^2
    //printf("Gaussian bit = ");arb_printd(G,20);printf("\n");
    pi_r_z_by_4(pir,z,L,prec); // Pi r z/4
    acb_sub_arb(ctmp1,pir,G,prec); // Pi r z/4 - Pi/h^2 (z-t0)^2
    acb_exp(ctmp2,ctmp1,prec); // exp ()
    //printf("exp(pi r z/4-pi(z-t0)^2/h^2) = ");acb_printd(ctmp2,20);printf("\n");
    
    Dexp(dexp,L,z,pi_by_H2,prec); // Pi r/4-2Pi(z-2t0)/h^2
    //printf("pi r/4-2 pi (z-2t0)/h^2 = ");acb_printd(dexp,20);printf("\n");
    acb_mul(ctmp1,dexp,W,prec); // *W
    acb_sub(ctmp3,W_dash,ctmp1,prec); // -W'
    acb_div(res,ctmp3,ctmp2,prec);
    acb_mul_onei(res,res);
    acb_neg(res,res);
    arb_clear(G);acb_clear(pir);
    acb_clear(ctmp1);acb_clear(ctmp2);acb_clear(ctmp3);acb_clear(dexp);
  }


  Lerror_t Lfunc_special_value(acb_t res, Lfunc_t LL, double alg_res, double alg_ims)
    {
      return Lfunc_special_value_choice(res,NULL,LL,alg_res,alg_ims,false,false);
    }
  
  // compute L(s) or Lam(s) into res where s is given in algebraic normalisation
  // if res_dash is given, compute L'(s) or Lam'(s) into res_dash
  Lerror_t Lfunc_special_value_choice(acb_t res, acb_t res_dash, Lfunc_t LL, double alg_res, double alg_ims, bool lam_p, bool do_dash)
  {
    double T=alg_ims;
    if(T<0.0)
      return ERR_SPEC_NZ; // need s in upper half plane
    
    Lerror_t ecode=ERR_SUCCESS;
    Lfunc *L=(Lfunc *) LL;
    int64_t prec=L->wprec;
    //printf("Algebraic s = %f + i%f\n",alg_res,alg_ims);

    arb_t arb_err,tmp;
    acb_t an_s; // analytic version of s
    acb_init(an_s);
    arb_init(arb_err);
    arb_init(tmp);
    arb_set_d(tmp, L->normalisation);
    arb_set_d(arb_err, alg_res);
    arb_sub(acb_realref(an_s),arb_err,tmp,prec); // an re(s) = alg re(s)- norm
    arb_set_d(acb_imagref(an_s),alg_ims);
    //printf("Analytic s = ");acb_printd(an_s,20);printf("\n");
    
    acb_t z;
    acb_init(z);
    arb_set(acb_realref(z),acb_imagref(an_s));
    arb_set_d(tmp,0.5);
    arb_sub(acb_imagref(z),tmp,acb_realref(an_s),prec);
    if(verbose) {printf("special value z = ");acb_printd(z,20);printf("\n");}    
    uint64_t stride=32; // should be computed dynamically
    double A=L->A/(double)stride;

    double h=sqrt(1.0/A)*1.001,best_h=h;
    double H=ceil(A*A*h*h/2.0),best_H=H;
    double M=H/A;
    double dimz=0.5-(alg_res - L->normalisation); // double approx to Im z
    if(verbose) printf("A = %f Im z = %f\n",A,dimz);
    ecode|=arb_upsampling_error(arb_err,M,H,h,A,L->mus,L->degree,L->conductor,T,acb_imagref(z),0,L->pi,prec); // first effort
    if(fatal_error(ecode))
    {
      arb_clear(tmp);
      acb_clear(z);
      acb_clear(an_s);
      arb_clear(arb_err);
      return ecode;
    }

    arb_t best_err;
    arb_init(best_err);
    arb_set(best_err,arb_err);
    while(true)
    {
      int64_t extra_bits=(int64_t)(M_PI*(dimz*dimz/(h*h)+fabs(dimz)*A)/M_LN2)+10;
      arb_mul_2exp_si(tmp,arb_err,L->target_prec+extra_bits);
      arb_sub_ui(tmp,tmp,1,prec);
      if(arb_is_negative(tmp)) // achieved target error
        break;
      h*=1.01;
      H=ceil(A*A*h*h/2.0);
      if(H*stride>L->u_no_values_off) // run out of data
	{
	  arb_set(arb_err,best_err);
	  h=best_h;H=best_H;
	  //ecode|=ERR_SPEC_PREC; // not necessarily, wait till end
	  break;
	}

      M=H/A;
      ecode|=arb_upsampling_error(arb_err,M,H,h,A,L->mus,L->degree,L->conductor,T,acb_imagref(z),0,L->pi,prec);
      if(fatal_error(ecode))
      {
        arb_clear(best_err);
        arb_clear(tmp);
        acb_clear(z);
        acb_clear(an_s);
        arb_clear(arb_err);
        return ecode;
      }
      arb_sub(tmp,best_err,arb_err,prec);
      if(arb_is_positive(tmp)) // found a better h,H so use them
      {
        best_h=h;
        best_H=H;
        arb_set(best_err,arb_err);
      }
    }
    arb_clear(best_err);

    arb_t arb_err_dash;
    arb_init(arb_err_dash);

    //printf("special value h set to %f\n",h);
    // use same parameters to compute error of differential
    if(do_dash)
      {
	ecode|=arb_upsampling_error(arb_err_dash,M,H,h,A,L->mus,L->degree,L->conductor,T,acb_imagref(z),1,L->pi,prec);
	if(fatal_error(ecode))
	  {
	    arb_clear(tmp);
	    acb_clear(z);
	    acb_clear(an_s);
	    arb_clear(arb_err);
	    arb_clear(arb_err_dash);
	    return ecode;
	  }
      }
    
    uint64_t iH=H;
    if(verbose) printf("H = %" PRIu64 " h = %f\n",iH,h);
    if(verbose) {printf("Upsample error set to ");arb_printd(arb_err,20);printf("\n");}

    arb_t pi_by_H2;
    arb_init(pi_by_H2);
    arb_set_d(tmp,h);
    arb_mul(tmp,tmp,tmp,prec);
    arb_div(pi_by_H2,L->pi,tmp,prec);
    arb_neg(pi_by_H2,pi_by_H2); // -Pi/h^2
    acb_t W;
    acb_init(W);
    ecode|=s_upsample_stride(W, z, T, L, prec, pi_by_H2, arb_err, iH, stride);
    if(fatal_error(ecode))
    {
      arb_clear(pi_by_H2);
      arb_clear(tmp);
      acb_clear(z);
      acb_clear(an_s);
      arb_clear(arb_err);
      arb_clear(arb_err_dash);
      acb_clear(W);
      return ecode;
    }
    //{printf("W(z) = ");acb_printd(W,20);printf("\n");}
    // go from W(z)->Lam(1/2+iz)
    acb_t Lam;
    acb_init(Lam);
    W_to_Lam(Lam,W,z,L,pi_by_H2,prec);
    if(verbose) {printf("Lambda(z) = ");acb_printd(Lam,20);printf("\n");}
    if(lam_p)
      acb_set(res,Lam);
    else
      // user wants L(s)
      // go from Lam->L by dividing out gamma(s)
      Lam_to_L(res,Lam,an_s,L,prec);
    /*
        // this code was checking that Lam was precise enough
	acb_abs(tmp,res,prec);
	arb_get_rad_arb(tmp,tmp);
	arb_mul_2exp_si(tmp,tmp,L->target_prec);
	arb_sub_ui(tmp,tmp,1,prec);
	if(!arb_is_negative(tmp))
	  ecode|=ERR_SPEC_PREC;
      }
    */
    
    if(!do_dash) // no need to do derivative
      {
	arb_clear(pi_by_H2);
	arb_clear(tmp);
	acb_clear(z);
	acb_clear(an_s);
	arb_clear(arb_err);
	arb_clear(arb_err_dash);
	acb_clear(W);
	acb_clear(Lam);
	return ecode;
      }

    // now do the derivative
    acb_t W_dash;
    acb_init(W_dash);
    ecode|=s_upsample_stride_dash(W_dash, z, T, L, prec, pi_by_H2, arb_err_dash, iH, stride);
    if(fatal_error(ecode))
      {
	arb_clear(pi_by_H2);
	arb_clear(tmp);
	acb_clear(z);
	acb_clear(an_s);
	arb_clear(arb_err);
	arb_clear(arb_err_dash);
	acb_clear(W);
	acb_clear(Lam);
	acb_clear(W_dash);
	return ecode;
      }
    
    // W_dash now contains W'(z)
    //printf("W'() = ");acb_printd(W_dash,20);printf("\n");
    acb_t Lam_dash;
    acb_init(Lam_dash);
    W_dash_to_Lam_dash(Lam_dash,W_dash,W,z,L,pi_by_H2,prec);
    if(lam_p)
      acb_set(res_dash,Lam_dash);
    else
      Lam_dash_to_L_dash(res_dash,Lam_dash,Lam,an_s,L,prec);
    
    acb_clear(Lam_dash);
    acb_clear(z);
    acb_clear(an_s);
    arb_clear(pi_by_H2);
    arb_clear(tmp);
    arb_clear(arb_err);
    arb_clear(arb_err_dash);
    acb_clear(W);
    acb_clear(Lam);
    acb_clear(W_dash);
    return ecode;
  }

#ifdef __cplusplus
}
#endif
