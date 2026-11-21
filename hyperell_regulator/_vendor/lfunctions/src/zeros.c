#include "glfunc.h"
#include "glfunc_internals.h"

#ifdef __cplusplus
extern "C"{
#endif


#define POS (1)
#define NEG (2)
#define UNK (3)
#define UP (1)
#define DOWN (2)
#define BOTH (3)
#define sign_t uint8_t
#define direction_t uint8_t

#define MAX_ZEROS (256)

sign_t sign(arb_t x) {
  if(arb_contains_zero(x))
    return UNK;
  if(arb_is_positive(x))
    return POS;
  return NEG;
}

// which way is the curve going?
direction_t direction(arb_t a, arb_t b, uint64_t prec) {
  static bool init = false;
  static arb_t tmp;
  if(!init)
  {
    init = true;
    arb_init(tmp);
  }
  arb_sub(tmp, a, b, prec);
  if(arb_contains_zero(tmp))
    return UNK;
  if(arb_is_negative(tmp)) // b>a
    return UP;
  return DOWN; // b<a
}

// binary chop - isolate as far as we can
Lerror_t isolate_zero(
    arb_t res,
    arb_t tt0, arb_t tt1,
    arb_t ff0, arb_t ff1,
    int8_t s0, // sign of ff0
    Lfunc *L, uint64_t side, uint64_t prec) {
  static bool init=false;
  static arb_t tmp1, tmp2, tmp3, t0, t1, f0, f1;
  if(!init) {
    init=true;
    arb_init(tmp1);
    arb_init(tmp2);
    arb_init(tmp3);
    arb_init(t0);
    arb_init(t1);
    arb_init(f0);
    arb_init(f1);
  }
  arb_set(t0, tt0);
  arb_set(t1, tt1);
  arb_set(f0, ff0);
  arb_set(f1, ff1);
  while(true) {
    arb_add(tmp1, t0, t1, prec);
    arb_mul_2exp_si(tmp1, tmp1, -1);
    if(!upsample_stride(tmp2, tmp1, L, side, prec))
      return ERR_UPSAMPLE;
    sign_t new_sign=sign(tmp2);
    if(new_sign==UNK) // can't do any better
    {
      arb_union(res, t0, t1, prec);
      return ERR_SUCCESS;
    }
    if(new_sign!=s0) // change is between t0 and tmp1
    {
      arb_set(t1, tmp1);
      arb_set(f1, tmp2);
    }
    else // change is between tmp1 and t1
    {
      arb_set(t0, tmp1);
      arb_set(f0, tmp2);
    }
  }
}


// called with suspected stationary point between m-1, m and m+1
// if !isolate_p, just confirm the two zeros to minimal precision
Lerror_t stat_point(arb_t z1, arb_t z2, uint64_t m, Lfunc *L, uint64_t side, uint64_t prec, bool isolate_p) {
  Lerror_t ecode=ERR_SUCCESS;
  static bool init=false;
  static arb_t t0, t1, f0, f1, t2, f2, t01, f01, t12, f12, tmp;
  if(!init)
  {
    init=true;
    arb_init(t0);
    arb_init(t1);
    arb_init(t2);
    arb_init(f0);
    arb_init(f1);
    arb_init(f2);
    arb_init(t01);
    arb_init(t12);
    arb_init(f01);
    arb_init(f12);
    arb_init(tmp);
  }

  arb_mul_ui(t0, L->one_over_A, m-1, prec);
  arb_mul_ui(t1, L->one_over_A, m, prec);
  arb_mul_ui(t2, L->one_over_A, m+1, prec);
  arb_mul_ui(tmp, t0, L->degree, prec);
  arb_mul(tmp, tmp, L->pi, prec);
  arb_mul_2exp_si(tmp, tmp, -2);
  arb_exp(tmp, tmp, prec);
  arb_mul(f0, L->u_values_off[side][m-1], tmp, prec);
  arb_mul_ui(tmp, t1, L->degree, prec);
  arb_mul(tmp, tmp, L->pi, prec);
  arb_mul_2exp_si(tmp, tmp, -2);
  arb_exp(tmp, tmp, prec);
  arb_mul(f1, L->u_values_off[side][m], tmp, prec);
  arb_mul_ui(tmp, t2, L->degree, prec);
  arb_mul(tmp, tmp, L->pi, prec);
  arb_mul_2exp_si(tmp, tmp, -2);
  arb_exp(tmp, tmp, prec);
  arb_mul(f2, L->u_values_off[side][m+1], tmp, prec);
  sign_t s=sign(f0);

  while(true) {
    if(verbose){printf("ti = ");arb_printd(t0, 20);printf(" ");arb_printd(t1, 20);printf(" ");arb_printd(t2, 20);printf("\n");}
    if(verbose){printf("fi = ");arb_printd(f0, 20);printf(" ");arb_printd(f1, 20);printf(" ");arb_printd(f2, 20);printf("\n");}
    arb_add(t01, t0, t1, prec);
    arb_mul_2exp_si(t01, t01, -1);
    if(!upsample_stride(f01, t01, L, side, prec))
      return ecode|ERR_STAT_POINT;
    if(verbose){printf("t01 = ");arb_printd(t01, 20);printf("\n");}
    if(verbose){printf("f01 = ");arb_printd(f01, 20);printf("\n");}
    sign_t s01=sign(f01);
    if(s01==UNK)
      return ecode|ERR_DBL_ZERO;
    if(s01!=s) {
      if(!isolate_p) {
        arb_union(z1, t0, t01, prec);
        arb_union(z2, t01, t1, prec);
        return ecode;
      }
      ecode|=isolate_zero(z1, t0, t01, f0, f01, s, L, side, prec);
      if(fatal_error(ecode))
        return ecode;
      ecode|=isolate_zero(z2, t01, t1, f01, f1, s01, L, side, prec);
      return ecode;
    }
    direction_t left=direction(f0, f01, prec);
    direction_t right=direction(f01, f1, prec);
    if(verbose){printf("f0 = ");arb_printd(f0, 20);printf("\n");}
    if(verbose){printf("f01 = ");arb_printd(f01, 20);printf("\n");}
    if(verbose){printf("f1 = ");arb_printd(f1, 20);printf("\n");}
    if((left==UNK)||(right==UNK))
      return ecode|ERR_DBL_ZERO;
    if(left!=right)
    {
      arb_set(t2, t1);
      arb_set(f2, f1);
      arb_set(t1, t01);
      arb_set(f1, f01);
      continue;
    }

    arb_add(t12, t1, t2, prec);
    arb_mul_2exp_si(t12, t12, -1);
    upsample_stride(f12, t12, L, side, prec);
    //printf("right middle = ");arb_printd(f12, 20);printf("\n");
    sign_t s12=sign(f12);
    if(s12==UNK)
      return ecode|ERR_DBL_ZERO;
    if(s12!=s)
    {
      if(!isolate_p)
      {
        arb_union(z1, t1, t12, prec);
        arb_union(z2, t12, t2, prec);
        return ecode;
      }
      ecode|=isolate_zero(z1, t1, t12, f1, f12, s, L, side, prec);
      if(fatal_error(ecode))
        return ecode;
      isolate_zero(z2, t12, t2, f12, f2, s12, L, side, prec);
      return ecode;
    }
    left=direction(f1, f12, prec);
    right=direction(f12, f2, prec);
    if((left==UNK)||(right==UNK))
      return ecode|ERR_DBL_ZERO;
    if(left!=right)
    {
      arb_set(t0, t1);
      arb_set(f0, f1);
      arb_set(t1, t12);
      arb_set(f1, f12);
      continue;
    }
    else
    {
      arb_set(t0, t01);
      arb_set(f0, f01);
      arb_set(t2, t12);
      arb_set(f2, f12);
      continue;
    }
  }
  return ecode;
}

// find some zeros
// errors:-
//   two UNK at start of data ERR_NO_DATA
//   direction unknown at start of data ERR_NO_DATA
//   UNK anywhere but at start ERR_SOME_DATA
Lerror_t find_zeros(Lfunc *L, uint64_t side)
{
  static bool init=false;
  static arb_t tmp, tmp1, tmp2, tmp3, t0, z1, z2;
  if(!init)
  {
    init=true;
    arb_init(tmp);
    arb_init(tmp1);
    arb_init(tmp2);
    arb_init(tmp3);
    arb_init(t0);
    arb_init(z1);
    arb_init(z2);
  }
  int64_t prec = L->wprec;
  bool stat_points = true;
  for(uint64_t z = 0; z < MAX_ZEROS; z++)
    arb_zero(L->zeros[side][z]);
  uint64_t count=0;

  uint64_t n=0;
  sign_t last_sign, this_sign=sign(L->u_values_off[side][n]);
  direction_t this_dir, last_dir;
  Lerror_t ecode=ERR_SUCCESS;

  if(this_sign == UNK) // central zero(s)
  {
    n++; // skip the central pt.
    this_sign=sign(L->u_values_off[side][n]);
    if(this_sign == UNK) // that was your last chance
      return ecode|ERR_NO_DATA;
  }

  this_dir=direction(L->u_values_off[side][n], L->u_values_off[side][n+1], prec);
  if(this_dir==UNK) {
    printf("Unknown direction at start of data.\n");
    return ecode|ERR_NO_DATA;
  }


  // now start searching for zeros and stat pts.
  while(true) {
    n++;
    if(n > L->fft_NN/OUTPUT_RATIO+L->fft_NN/TURING_RATIO)
      return ecode;
    last_sign=this_sign;
    this_sign=sign(L->u_values_off[side][n]);

    if(this_sign==UNK) // run out of precision
      return ecode|ERR_SOME_DATA;

    if(stat_points) {
      last_dir=this_dir;
      this_dir=direction(L->u_values_off[side][n-1], L->u_values_off[side][n], prec);
      if(this_dir==UNK) // time to stop doing stat points
        stat_points=false;
    }


    if(this_sign!=last_sign) // found a zero between n and n-1
    {
      arb_mul_ui(tmp1, L->one_over_A, n-1, prec);
      arb_mul_ui(tmp2, L->one_over_A, n, prec);
      if(verbose)
      {
        printf("zero found between ");arb_printd(tmp1, 20);
        printf(" and ");arb_printd(tmp2, 20);printf("\n");
      }
      ecode|=isolate_zero(L->zeros[side][count++], tmp1, tmp2, L->u_values_off[side][n-1], L->u_values_off[side][n], last_sign, L, side, prec);
      if(fatal_error(ecode)||(count==MAX_ZEROS))
        return ecode;
      continue;
    }

    // didn't find a sign change, so let's see about stationary points

    if(!stat_points) // not bothering any more
      continue;

    if(this_dir!=last_dir) { // change in direction
      if(((last_dir==UP)&&(this_sign==NEG))||((last_dir==DOWN)&&(this_sign==POS))) {
        if(verbose) printf("Stationary point detected.\n");
        ecode|=stat_point(z1, z2, n-1, L, side, prec, n<=L->fft_NN/OUTPUT_RATIO);
        if(fatal_error(ecode))
          return ecode;
        if(verbose)
        {
          printf("Stat zero found at ");arb_printd(z1, 20);
          printf("\nand                ");arb_printd(z2, 20);
          printf("\n");fflush(stdout);
        }
        if(verbose){printf("setting zeros %" PRIu64 " and %" PRIu64 " to ", count, count+1);
          arb_printd(z1, 20);printf(" ");arb_printd(z2, 20);printf("\n");}
        arb_set(L->zeros[side][count], z1);
        count++;
        if(count==MAX_ZEROS)
          return ecode;
        arb_set(L->zeros[side][count], z2);
        count++;
        if(count==MAX_ZEROS)
          return ecode;
      }
    }
  }
  return ecode;
}


#ifdef __cplusplus
}
#endif
