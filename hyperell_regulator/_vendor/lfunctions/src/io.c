#include "glfunc.h"
#include "glfunc_internals.h"

#ifdef __cplusplus
extern "C"{
#endif

acb_srcptr Lfunc_sign(Lfunc_t L)
{
  Lfunc *LL=(Lfunc *) L;
  return (acb_srcptr)LL->sign;
}

acb_srcptr Lfunc_sqrt_sign(Lfunc_t L)
{
  Lfunc *LL=(Lfunc *) L;
  return (acb_srcptr)LL->sqrt_sign;
}

arb_srcptr Lfunc_zeros(Lfunc_t L, uint64_t side)
{
  if(side>1)
    return (arb_srcptr) NULL;
  Lfunc *LL=(Lfunc *) L;
  return (arb_srcptr)LL->zeros[side];
}

int64_t Lfunc_rank(Lfunc_t L)
{
  Lfunc *LL=(Lfunc *) L;
  return LL->rank;
}

double arb_getd(arb_t x)
{
  double res=arf_get_d(arb_midref(x),ARF_RND_NEAR);
  return res;
}

double normalised(Lfunc *L, uint64_t side, uint64_t ptr, double t)
{
  static bool init=false;
  static acb_t s;
  static arb_t tmp1,tmp2;
  if(!init)
    {
      init=true;
      acb_init(s);
      arb_init(tmp1);
      arb_init(tmp2);
    }

  arb_set_d(acb_realref(s),0.5);
  arb_set_d(acb_imagref(s),t);
  abs_gamma(tmp1,s,L,100);
  arb_div(tmp2,L->u_values_off[side][ptr],tmp1,100);
  return arb_getd(tmp2);
}

Lplot_t *Lfunc_plot_data(Lfunc_t LL, uint64_t side, double max_t, uint64_t n_points)
{
  Lfunc *L=(Lfunc *) LL;
  Lplot_t *Lp;
  Lp=(Lplot_t *)malloc(sizeof(Lplot_t));
  if(!Lp)
    return (Lplot_t *) NULL;

  if(max_t>512.0/(double) (L->degree*OUTPUT_RATIO))
    max_t=512.0/(double) (L->degree*OUTPUT_RATIO);
  double pts=max_t*L->A;
  uint64_t step=ceil(pts/(double)n_points);
  // we are going to return every <step>'th point upto max_t
  Lp->spacing=(double)step/L->A;
  Lp->n_points=ceil(max_t/Lp->spacing)+1;
  //printf("%f %lu %f %lu\n",pts,step,Lp->spacing,Lp->n_points);

  Lp->points=(double *)malloc(sizeof(double)*Lp->n_points);
  if(!Lp->points)
    {
      free(Lp);
      return (Lplot_t *) NULL;
    }

  for(uint64_t i=0;i<Lp->n_points;i++)
    Lp->points[i]=normalised(L,side,i*step,(double)i*Lp->spacing);
  return Lp;
}

arb_srcptr Lfunc_Taylor(Lfunc_t LL)
{
  Lfunc *L=(Lfunc *) LL;
  return (arb_srcptr) L->L_d;
}

void Lfunc_clear_plot(Lplot_t *Lp)
{
  if(Lp)
    {
      free(Lp->points);
      free(Lp);
    }
}

#ifdef __cplusplus
}
#endif
