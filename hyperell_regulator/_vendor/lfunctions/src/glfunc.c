#include "glfunc.h"
#include "assert.h"
#include "glfunc_internals.h"
#include "math.h"
#include "stdlib.h"

#ifdef __cplusplus
extern "C" {
#endif

bool fatal_error(Lerror_t ecode) { return ecode & 0xFFFFFFFF; }

void fprint_errors(FILE *f, Lerror_t ecode) {
  // fatal errors
  if (ecode & ERR_OOM)
    fprintf(f, "Ran out of memory somewhere.\n");
  if (ecode & ERR_M_ERROR)
    fprintf(f, "Error computing M. Failed Lemma 2 of M_error.tex.\n");
  if (ecode & ERR_NO_DATA)
    fprintf(f, "Looks like we have no usable data.\n");
  if (ecode & ERR_ZERO_ERROR)
    fprintf(f, "Fatal error looking for zeros.\n");
#ifdef BUTHE
  if (ecode & ERR_BUT_ERROR)
    fprintf(
        f,
        "Error doing Buthe check. Estimate for Wf+Winf-Ws* must allow >=0.\n");
#endif
  if (ecode & ERR_UPSAMPLE)
    fprintf(f, "Error computing bounds for upsampling.\n");
  if (ecode & ERR_MU_HALF)
    fprintf(f, "We expect all mu's to be 1/2 non-negative integers.\n");
  if (ecode & ERR_STAT_POINT)
    fprintf(f, "Fatal error in stationary point routine.\n");
  if (ecode & ERR_SPEC_VALUE)
    fprintf(f, "Fatal error in special value routine.\n");
  // warnings
  if (ecode & ERR_INSUFF_EULER)
    fprintf(f, "Don't appear to have enough Euler factors.\n");

  if (ecode & ERR_SOME_DATA)
    fprintf(f, "Data became unusable after output region but before end of "
               "Turing zone.\n");
  if (ecode & ERR_ZERO_PREC)
    fprintf(f, "Couldn't isolate all zeros to requested precision.\n");
  if (ecode & ERR_NO_RANK)
    fprintf(f, "Could not determine rank of L.\n");
  if (ecode & ERR_CONFLICT_RANK)
    fprintf(f, "Computed rank did not agree with what we were told.\n");
  if (ecode & ERR_RH_ERROR)
    fprintf(f, "Failed to confirm RH for zeros in output region.\n");
  if (ecode & ERR_DBL_ZERO)
    fprintf(
        f,
        "Stationary point routine failed to converge. Possible double zero?\n");
  if (ecode & ERR_SPEC_PREC)
    fprintf(
        f, "Failed to achieve desired error bound in Special Value routine.\n");
  if (ecode & ERR_G_INFILE)
    fprintf(f, "Problem opening cached G data file.\n");
  if (ecode & ERR_G_OUTFILE)
    fprintf(f, "Problem opening file to cache G data.\n");
  if (ecode & ERR_G_EXTENT)
    fprintf(f,
            "G data does not extend low enough for this conductor (the cached "
            "or computed grid floor is too high).\n");
  if (ecode & ERR_BAD_DEGREE)
    fprintf(f, "The degree of the L-function must be between 1 and %d\n",
            MAX_DEGREE + 1);
  if (ecode & ERR_SPEC_NZ)
    fprintf(f, "Special value routine only works for Im s>=0.\n");
}

// what is the decay (in bits) in the gamma factors from 1/2 to 1/2+i(64/r)
uint64_t decay(Lfunc *L) {
  arb_t tmp1, tmp2, tmp3;
  acb_t s;
  arb_init(tmp1);
  arb_init(tmp2);
  arb_init(tmp3);
  acb_init(s);
  arb_set_d(acb_realref(s), 0.5);
  abs_gamma(tmp1, s, L, 100);
  arb_set_d(acb_imagref(s), 64.0 / (double)L->degree);
  abs_gamma(tmp2, s, L, 100);
  arb_div(tmp3, tmp1, tmp2, 100);
  if (verbose) {
    printf("Gamma factor decay is ");
    arb_printd(tmp3, 20);
    printf("\n");
  }
  arb_set_ui(tmp1, 1);
  for (uint64_t j = 0;; j++) {
    if (arb_lt(tmp3, tmp1)) {
      arb_clear(tmp1);
      arb_clear(tmp2);
      arb_clear(tmp3);
      acb_clear(s);
      return j;
    }
    arb_mul_2exp_si(tmp1, tmp1, 1);
  }
}

bool is_half_int(double x) {
  return (x >= 0.0) && ((2.0 * x) == ceil(2.0 * x)) &&
         ((2.0 * x) == floor(2.0 * x));
}

int double_comp(const void *a, const void *b) {
  double *x = (double *)a;
  double *y = (double *)b;
  if (*x < *y)
    return -1;
  else if (*x > *y)
    return 1;
  else
    return 0;
}

Lfunc_t Lfunc_init_advanced(Lparams_t *Lp, Lerror_t *ecode) {
  ecode[0] = ERR_SUCCESS;
  uint64_t i, j;
  arb_t tmp;

  if ((Lp->degree < 2) || (Lp->degree > MAX_DEGREE)) {
    ecode[0] |= ERR_BAD_DEGREE;
    return ((Lfunc_t)NULL);
  }

  Lfunc *L = (Lfunc *)malloc(sizeof(Lfunc));
  if (!L) {
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }
  L->degree = Lp->degree;
  L->normalisation = Lp->normalisation;
  L->conductor = Lp->conductor;
  L->mus = (double *)malloc(sizeof(double) * L->degree);
  if (!L->mus) {
    ecode[0] |= ERR_OOM;
    free(L);
    return (Lfunc_t)NULL;
  }
  for (i = 0; i < L->degree; i++) {
    L->mus[i] = Lp->mus[i] + Lp->normalisation; // alg->anal
    if (!is_half_int(L->mus[i])) {
      ecode[0] |= ERR_MU_HALF;
      // only L and L->mus are live here (no arb_t initialised yet)
      free(L->mus);
      free(L);
      return (Lfunc_t)NULL;
    }
  }

  // we sort the mus so we can name G cache files canonically
  qsort(L->mus, L->degree, sizeof(double), double_comp);

  L->target_prec = Lp->target_prec;
  arb_init(L->zero_prec);
  arb_set_ui(L->zero_prec, 1);
  arb_mul_2exp_si(L->zero_prec, L->zero_prec, -L->target_prec - 1);
  arb_init(L->zero_error);
  arb_add_error(L->zero_error, L->zero_prec);
  arb_init(L->pi);
  if (Lp->wprec > 0)
    L->wprec = Lp->wprec;
  else {
    arb_const_pi(L->pi, 100); // for now, needed by decay()
    // allow enough bits so we will get target_prec at height 1/2+i*64/degree
    L->wprec = L->target_prec + decay(L) + EXTRA_BITS;
    if (verbose)
      printf("working precision set to %" PRId64 "\n", L->wprec);
  }
  arb_const_pi(L->pi, L->wprec); // set it properly now we know what wprec is
  L->gprec = Lp->gprec;
  L->self_dual = Lp->self_dual;
  L->rank = Lp->rank;
  L->cache_dir = Lp->cache_dir;

  // See Lemma 2 of M_error1.pdf, Lemma 5 of g.pdf
  // r is always >=2
  L->nus = (arb_t *)malloc(sizeof(arb_t) * L->degree);
  if (!L->nus) {
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }

  for (j = 0; j < 2; j++) {
    arb_init(L->nus[j]);
    arb_set_d(L->nus[j], 0.5 * (L->mus[j] - 0.5));
  }
  for (j = 2; j < L->degree; j++) {
    arb_init(L->nus[j]);
    arb_set_d(L->nus[j], 0.5 * (L->mus[j] - 1.0));
  }

  // See Lemma 2/5 again
  arb_init(L->nu);
  arb_set(L->nu, L->nus[0]);
  for (j = 1; j < L->degree; j++)
    arb_add(L->nu, L->nu, L->nus[j], L->target_prec);
  arb_div_ui(L->nu, L->nu, L->degree, L->target_prec);
  arb_mul_2exp_si(L->nu, L->nu, 1);
  arb_init(tmp);
  arb_set_d(tmp, 0.5);
  arb_add(L->nu, L->nu, tmp, L->target_prec);

  // mu=-1/2+1/r(1+sum mu_j) See Lemma 5.2 of Artin
  arb_init(L->mu);
  arb_set_d(L->mu, -0.5);
  double smu = 1.0;
  for (i = 0; i < L->degree; i++)
    smu += L->mus[i];
  arb_set_d(tmp, smu);
  arb_div_ui(tmp, tmp, L->degree, L->target_prec);
  arb_add(L->mu, L->mu, tmp, L->target_prec);

  ecode[0] |= compute_g(L);
  if (fatal_error(ecode[0]))
    return (Lfunc_t)NULL;
  if (verbose) {
    printf("eq 5-9 error = ");
    arb_printd(L->eq59, 20);
    printf("\n");
  }

  arb_init(L->B);
  arb_set_d(tmp, L->one_over_B);
  arb_inv(L->B, tmp, L->wprec);

  arb_init(L->two_pi_by_B);
  arb_set_d(L->two_pi_by_B, L->one_over_B * 2.0);
  arb_mul(L->two_pi_by_B, L->two_pi_by_B, L->pi, L->wprec);

  L->fft_N = 1 << 11;  // length of DTF for convolutions
  L->fft_NN = 1 << 16; // final output length

  L->A = L->fft_NN * L->one_over_B;
  arb_init(L->arb_A);
  arb_set_d(L->arb_A, L->A);
  arf_init(L->arf_A);
  arf_set_d(L->arf_A, L->A);

  arb_init(L->one_over_A);
  arb_inv(L->one_over_A, L->arb_A, L->wprec);
  arf_init(L->arf_one_over_A);
  arf_ui_div(L->arf_one_over_A, 1, L->arf_A, L->wprec, ARF_RND_NEAR);

  L->G = (acb_t *)malloc(sizeof(acb_t) * L->fft_N);
  if (!L->G) {
    arb_clear(tmp);
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }
  for (i = 0; i < L->fft_N; i++)
    acb_init(L->G[i]);

  L->eta = 0.0;
  arb_init(L->delta);
  arb_mul_2exp_si(L->delta, L->pi, -1); // pi/2
  arb_set_d(tmp, 1.0 - L->eta);
  arb_mul(L->delta, L->delta, tmp, L->wprec); // (1-eta)pi/2
  arb_init(L->exp_delta);
  arb_neg(tmp, L->delta);
  arb_exp(L->exp_delta, tmp, L->wprec);

  L->w =
      (acb_t *)malloc(sizeof(acb_t) * L->fft_N / 2); // twiddles for little FFT
  if (!L->w) {
    arb_clear(tmp);
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }

  for (i = 0; i < L->fft_N / 2; i++)
    acb_init(L->w[i]);
  acb_initfft(L->w, L->fft_N, L->wprec); // set twiddles for little FFT

  L->ww =
      (acb_t *)malloc(sizeof(acb_t) * L->fft_NN / 2); // twiddles for big FFT
  if (!L->ww) {
    arb_clear(tmp);
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }
  for (i = 0; i < L->fft_NN / 2; i++)
    acb_init(L->ww[i]);
  acb_initfft(L->ww, L->fft_NN, L->wprec); // set twiddles for big FFT

  // space for the zeros once we isolate them
  L->zeros[0] = (arb_t *)malloc(sizeof(arb_t) * MAX_ZEROS);
  L->zeros[1] = (arb_t *)malloc(sizeof(arb_t) * MAX_ZEROS);
  if ((!L->zeros[0]) || (!L->zeros[1])) {
    arb_clear(tmp);
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }
  for (i = 0; i < MAX_ZEROS; i++) {
    arb_init(L->zeros[0][i]);
    arb_init(L->zeros[1][i]);
  }

  L->kres = (acb_t *)malloc(sizeof(acb_t) * L->fft_N);
  if (!L->kres) {
    arb_clear(tmp);
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }

  L->skm = (acb_t **)malloc(sizeof(acb_t *) * L->max_K);
  if (!L->skm) {
    arb_clear(tmp);
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }

  uint64_t k, n;
  for (k = 0; k < L->max_K; k++) {
    L->skm[k] = (acb_t *)malloc(sizeof(acb_t) * L->fft_N);
    if (!L->skm[k]) {
      arb_clear(tmp);
      ecode[0] |= ERR_OOM;
      return (Lfunc_t)NULL;
    }

    for (n = 0; n < L->fft_N; n++)
      acb_init(L->skm[k][n]);
  }

  L->res = (acb_t *)malloc(sizeof(acb_t) * L->fft_NN);
  if (!L->res) {
    arb_clear(tmp);
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }

  for (n = 0; n < L->fft_N; n++)
    acb_init(L->kres[n]);
  for (n = 0; n < L->fft_NN; n++)
    acb_init(L->res[n]);

  arb_init(L->pre_ftwiddle_error);
  arb_init(L->ftwiddle_error);
  init_ftwiddle_error(L, L->wprec);

  arb_init(L->one_over_root_N);
  arb_init(L->sum_ans);
  acb_init(L->sign);
  acb_init(L->sqrt_sign);
  L->allocated_M = 8192;
  L->ans = (acb_t *)malloc(sizeof(acb_t) * L->allocated_M);
  if (!L->ans) {
    arb_clear(tmp);
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }
  for (size_t i = 0; i < L->allocated_M; ++i)
    acb_init(L->ans[i]);

  arb_clear(tmp);

#ifdef BUTHE
  init_buthe(L, L->wprec); // setup stuff for Buthe zero check
#endif
#ifdef TURING
  arb_init(L->imint);
  arb_init(L->X);
#endif
  L->nmax_called = false; // noone has called nmax yet

  arb_init(L->Lam_d);
  arb_init(L->L_d);

  ecode[0] |= init_upsampling(L);

  return (Lfunc_t)L;
}

Lfunc_t Lfunc_init(uint64_t degree, uint64_t conductor, double normalisation,
                   const double *mus, Lerror_t *ecode) {
  Lparams_t Lp;
  Lp.degree = degree;
  Lp.conductor = conductor;
  Lp.normalisation = normalisation;
  Lp.mus = (double *)malloc(sizeof(double) * degree);
  if (!Lp.mus) {
    ecode[0] |= ERR_OOM;
    return (Lfunc_t)NULL;
  }
  for (size_t i = 0; i < degree; ++i)
    Lp.mus[i] = mus[i];
  Lp.target_prec = DEFAULT_TARGET_PREC;
  Lp.rank = DK;
  Lp.self_dual = DK;
  Lp.cache_dir = ".";
  Lp.gprec = 0; // We will try to do something sensible
  Lp.wprec = 0; // ditto

  // Lfunc_init_advanced copies Lp.mus into L->mus, so this scratch array is ours
  // to free; otherwise it leaks (Lp is a stack local that goes out of scope).
  Lfunc_t L = Lfunc_init_advanced(&Lp, ecode);
  free(Lp.mus);
  return L;
}

int64_t Lfunc_wprec(Lfunc_t Lf) {
  Lfunc *L;
  L = (Lfunc *)Lf;
  return L->wprec;
}

#ifdef __cplusplus
}
#endif
