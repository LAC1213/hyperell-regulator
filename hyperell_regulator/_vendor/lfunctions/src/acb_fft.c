#include <flint/acb.h>
#include "inttypes.h"

void acb_initfft(acb_t *w, uint64_t n, int64_t prec)
{
  arb_t I,IN;
  arb_init(I);
  arb_init(IN);
  acb_one(w[0]);
  for(uint64_t i=1;i<n/2;++i)
    {
      arb_set_ui(I,2*i);
      arb_div_ui(IN,I,n,prec);
      arb_sin_cos_pi(acb_imagref(w[i]),acb_realref(w[i]),IN,prec);
    }
  arb_clear(I);arb_clear(IN);
} /* acb_initfft */


// do inplace fft of x of length n a power of 2
// w[i]=e(i/n) i=0..n/2-1
void acb_fft(acb_t *x, uint64_t n, acb_t *w, int64_t prec) {
  acb_t tmp;
  acb_init(tmp);

  uint64_t i,j,k,l;
  acb_t *p,*xend=x+n;

  for(i=0,l=n>>1;i<l;++i)
    {
      for(k=1,j=0;k<n;k<<=1)
	{
	  j<<=1;
	  if(i&k)
	    j|=1;
	}
      if(i<j)
	acb_swap(x[i],x[j]);
      else if (i>j)
	acb_swap(x[n-1-i],x[n-1-j]);
      ++i;
      j|=l;
      acb_swap(x[i],x[j]);
    }

  for(k=1,l=n/2;k<n;k<<=1,l>>=1)
    for(p=x;p<xend;p+=k)
      for(j=0;j<n/2;j+=l,p++)
	{
	  acb_mul(tmp,p[k],w[j],prec);
	  acb_sub(p[k],p[0],tmp,prec);
	  acb_add(p[0],p[0],tmp,prec);
	}
  acb_clear(tmp);
} /* acb_fft */

// non normalised inverse dft
void acb_ifft(acb_t *x, uint64_t n, acb_t *w, uint64_t prec)
{
  acb_fft(x,n,w,prec);
  for(uint64_t i=1;i<n/2;++i)
    acb_swap(x[i],x[n-i]);
}

// x,y must be distinct
void acb_convolve(acb_t *res, acb_t *x, acb_t *y, uint64_t n, acb_t *w, uint64_t prec)
{
  acb_fft(x,n,w,prec);
  acb_fft(y,n,w,prec);
  for(uint64_t i=0;i<n;++i)
    acb_mul(res[i],x[i],y[i],prec);
  acb_ifft(res,n,w,prec);
  for(uint64_t i=0;i<n;++i)
    acb_div_ui(res[i],res[i],n,prec);
}

// x,y must be distinct. y has already been fft'd
void acb_convolve1(acb_t *res, acb_t *x, acb_t *y, uint64_t n, acb_t *w, uint64_t prec)
{
  acb_fft(x,n,w,prec);
  //fft(y,n,w,prec);
  for(uint64_t i=0;i<n;++i)
    acb_mul(res[i],x[i],y[i],prec);
  acb_ifft(res,n,w,prec);
  for(uint64_t i=0;i<n;++i)
    acb_div_ui(res[i],res[i],n,prec);
}

// x and y have already been fft'd
void acb_convolve2(acb_t *res, acb_t *x, acb_t *y, uint64_t n, acb_t *w, uint64_t prec) {
  //acb_fft(x,n,w,prec);
  //fft(y,n,w,prec);
  for(uint64_t i=0;i<n;++i)
    acb_mul(res[i],x[i],y[i],prec);
  acb_ifft(res,n,w,prec);
  for(uint64_t i=0;i<n;++i)
    acb_div_ui(res[i],res[i],n,prec);
}

