#include "glfunc.h"
#include "glfunc_internals.h"

#ifdef __cplusplus
extern "C"{
#endif

#define arb_cclear(x) if(x) arb_clear(x)
#define acb_cclear(x) if(x) acb_clear(x)
#define arf_cclear(x) if(x) arf_clear(x)

  void Lfunc_clear(Lfunc_t LL)
  {
    Lfunc *L=(Lfunc *) LL;

    free(L->mus);
    arb_cclear(L->zero_prec);
    arb_cclear(L->zero_error);
    arb_cclear(L->mu);
    arb_cclear(L->nu);
    if(L->nus)
      {
	for(uint64_t i=0;i<L->degree;i++)
	  arb_clear(L->nus[i]);
	free(L->nus);
      }
    arb_cclear(L->C);
    arb_cclear(L->alpha);
    arb_cclear(L->B);
    arb_cclear(L->two_pi_by_B);
    arb_cclear(L->pi);
    arb_cclear(L->eq59);
    if(L->Gs)
      {
	for(uint64_t k=0;k<L->max_K;k++)
	  if(L->Gs[k])
	    {
	      for(int64_t i=0;i<=L->hi_i-L->low_i;i++)
		arb_cclear(L->Gs[k][i]);
	      free(L->Gs[k]);
	    }
	free(L->Gs);
      }

    arb_cclear(L->arb_A);
    arb_cclear(L->one_over_A);
    for(uint64_t side=0;side<2;side++)
      if(L->zeros[side])
	{
	  for(uint64_t z=0;z<MAX_ZEROS;z++)
	    arb_cclear(L->zeros[side][z]);
	  free(L->zeros[side]);
	}
    arb_cclear(L->delta);
    arb_cclear(L->exp_delta);
    arb_cclear(L->pre_ftwiddle_error);
    arb_cclear(L->ftwiddle_error);
    
    #ifdef BUTHE
    arb_cclear(L->buthe_Wf);
    arb_cclear(L->buthe_Winf);
    arb_cclear(L->buthe_Ws);
    arb_cclear(L->buthe_b);
    arb_cclear(L->buthe_C);
    arb_cclear(L->buthe_h);
    for(uint64_t i=0;i<(MAX_R-1)*(2*MAX_MUI_2+1);i++)
      arb_cclear(L->buthe_ints[i]);
    #endif

    #ifdef TURING
    arb_cclear(L->imint);
    arb_cclear(L->X);
    #endif
    
    arb_cclear(L->one_over_root_N);
    arb_cclear(L->sum_ans);
    arb_cclear(L->u_H);
    arb_cclear(L->u_pi_by_H2);
    arb_cclear(L->u_A);
    arb_cclear(L->u_one_over_A);
    for(uint64_t side=0;side<2;side++)
      if(L->u_values[side])
	{
	  for(uint64_t i=0;i<L->u_no_values;i++)
	    arb_cclear(L->u_values[side][i]);
	  free(L->u_values[side]);
	}
    arb_cclear(L->u_pi_A);
    arb_cclear(L->upsampling_error);
    arb_cclear(L->Lam_d);
    arb_cclear(L->L_d);
    
    if(L->G)
      {
	for(uint64_t i=0;i<L->fft_N;i++)
	  acb_cclear(L->G[i]);
	free(L->G);
      }
    if(L->w)
      {
	for(uint64_t i=0;i<L->fft_N/2;i++)
	  acb_cclear(L->w[i]);
	free(L->w);
      }
    if(L->ww)
      {
	for(uint64_t i=0;i<L->fft_NN/2;i++)
	  acb_cclear(L->ww[i]);
	free(L->ww);
      }
    if(L->kres)
      {
	for(uint64_t i=0;i<L->fft_N;i++)
	  acb_cclear(L->kres[i]);
	free(L->kres);
      }
    if(L->skm)
      {
	for(uint64_t k=0;k<L->max_K;k++)
	  if(L->skm[k])
	    {
	      for(uint64_t i=0;i<L->fft_N;i++)
		acb_cclear(L->skm[k][i]);
	      free(L->skm[k]);
	    }
	free(L->skm);
      }
    if(L->res)
      {
	for(uint64_t i=0;i<L->fft_NN;i++)
	  acb_cclear(L->res[i]);
	free(L->res);
      }
    acb_cclear(L->sign);
    acb_cclear(L->sqrt_sign);
    if(L->ans)
      {
	for(uint64_t i=0;i<L->allocated_M;i++)
	  acb_cclear(L->ans[i]);
	free(L->ans);
      }

    arf_cclear(L->arf_A);
    arf_cclear(L->arf_one_over_A);

    free(L);
  }
#ifdef __cplusplus
}
#endif
