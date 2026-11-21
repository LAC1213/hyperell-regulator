#ifndef _GLFUNC_INCLUDE
#define _GLFUNC_INCLUDE

#include <inttypes.h>
#include <flint/acb_poly.h>
#include <stdbool.h>

#define DK (-1) // tri-state "don't know" (for self_dual / rank)
#define YES (1) // tri-state "yes" (for self_dual / rank)
#define NO (0) // tri-state "no" (for self_dual / rank)

//#define BUTHE // if defined, verify RH via Buthe's method (off by default)
#define TURING // if defined, verify RH via Booker/Turing's method (default)

#define MAX_DEGREE (9) // if increasing, need more integrals in buthe.c
// and probably need to take a good look at g.c
#define MAX_R MAX_DEGREE // alias for MAX_DEGREE (max degree r); sizes the Buthe tables
#define MAX_ZEROS (256) // hard cap on the number of zeros found/stored per side
#define DEFAULT_TARGET_PREC (100) // default target precision in bits when none is given
// error codes
// those in lower 32 bits are fatal
// those in upper 32 bits are warnings
#define ERR_SUCCESS (0) // no error
// fatalities
#define ERR_NO_DATA (1) // the first two Lambda(t) values contained zero. Totally fatal
#define ERR_ZERO_ERROR (2) // some unexpected error isolating zeros
#ifdef BUTHE
#define ERR_BUT_ERROR (4) // Wf+Winf-Ws* must be negative (we cant find too many zeros)
#endif
#define ERR_OOM (8) // out of memory error
#define ERR_UPSAMPLE (16) // something bad happened trying to upsample
#define ERR_MU_HALF (32) // mus should be 1/2 integers
#define ERR_M_ERROR (64) // fatal: failed to compute the M (coeff-count) bound (M_error Lemma 2)
#define ERR_STAT_POINT (128) // fatal error in stat_point
#define ERR_SPEC_VALUE (256) // fatal error in Special Value
#define ERR_G_INFILE ((uint64_t) 512) // fatal error reading g_data from cache
#define ERR_BAD_DEGREE ((uint64_t) 1024) //fatal error when the degree is too low or too high
#define ERR_SPEC_NZ ((uint64_t) 2048) // special value routine requires Im s >= 0.
#define ERR_G_EXTENT ((uint64_t) 4096) // fatal: G grid does not extend low enough (conductor too large for the fixed grid floor, or a cached grid was reused)

// warnings
#define ERR_SOME_DATA ((uint64_t) 1<<32) // We had some sensible data, but not to end of Turing Zone
#define ERR_ZERO_PREC ((uint64_t) 1<<33) // couldn't isolate zeros to target_prec
#define ERR_RH_ERROR ((uint64_t) 1<<34) // RH check (using Buthe's or Turing's method) failed
#define ERR_INSUFF_EULER ((uint64_t) 1<<35) // ran out of Euler factors before we expected to
#define ERR_NO_RANK ((uint64_t) 1<<36) // could not determine rank of L
#define ERR_CONFLICT_RANK ((uint64_t) 1<<37) // rank we computed did not agree with what we were told
#define ERR_DBL_ZERO ((uint64_t) 1<<38) // stationary point failed to converge. Double zero?
#define ERR_SPEC_PREC ((uint64_t) 1<<39) // could not achieve target error bound in special value
#define ERR_G_OUTFILE ((uint64_t) 1<<40) // problem opening file to write g_data cache

#ifdef __cplusplus
extern "C"{
#endif

  typedef uint64_t Lerror_t;

  // keep details under wraps
  typedef void *Lfunc_t;

  typedef struct{
    uint64_t degree;
    uint64_t conductor;
    double normalisation;
    double *mus;
    int64_t target_prec;
    int64_t wprec;
    int64_t gprec;
    int self_dual; // -1 = DK, 0 = No, 1 = Yes
    int rank; // -1 = DK
    char *cache_dir;
  } Lparams_t;

  typedef struct{
    uint64_t n_points;
    double spacing;
    double *points;
  } Lplot_t;

  bool fatal_error(Lerror_t ecode); // are any errors in ecode considered fatal?
  void fprint_errors(FILE *f, Lerror_t ecode); // print all errors with newlines

  // return initialised Lfunc structure
  /* Input:
   *  - degree, the degree of the L-function
   *  - conductor, the conductor of the L-function. There is no separate
   *      "number of primes/coefficients" parameter: the coefficient range M is
   *      derived from the conductor, computed lazily on the first Lfunc_nmax
   *      call (triggered by Lfunc_use_all_lpolys / Lfunc_reduce_nmax) and then
   *      frozen. See Lfunc_nmax for the formula and the ERR_G_EXTENT caveat.
   *  - normalisation, the shift on s axis to go from the algebraic normalization to the analytic one, i.e., if an is the Dirichlet in the algebraic normalization, then an/n^{normalisation} is the Dirichlet coefficient of in the analytic normalization. Lambda(s)=eps Lambda(k-s) -> normalisation = (k-1)/2
   *  - mus, the shifts of Gamma_R, mu[i] + normalisation must be half integers
   *  - ecode, where we keep track of errors and warnings
   *
   *  Example:
   *    - EC/Q:
   *        mus = [0, 1] and normalisation = 0.5
   *      or
   *        mus = [0.5, 1.5] and normalisation = 0
   *
   *    - Classical modular form of weight 13 (motivic weight of Lfunc = 12)
   *        mus = [0, 1] and normalisation = 6
   *      or
   *        mus = [6, 7] and normalisation = 0
   */
  Lfunc_t Lfunc_init(uint64_t degree, uint64_t conductor, double normalisation, const double *mus, Lerror_t *ecode);
  // do the same but with more control. Lparams.conductor fixes the coefficient
  // range M implicitly, as in Lfunc_init (see Lfunc_nmax); no Lparams field sets
  // the prime/coefficient count directly.
  Lfunc_t Lfunc_init_advanced(Lparams_t *Lparams, Lerror_t *ecode);

  // Largest prime p (equivalently largest index M) for which an Euler factor /
  // coefficient a_n is expected. M is derived from the conductor, not passed in:
  //   dc = sqrt(conductor);  M = floor(dc * exp(2*pi*(hi_i + 0.5) / B))
  // (hi_i and B are fixed by degree/mus/precision at init). It is computed on
  // the first call and then frozen (the nmax_called flag); Lfunc_use_all_lpolys
  // and Lfunc_reduce_nmax trigger that first call, so the range is set once any
  // of them runs.
  //
  // ERR_G_EXTENT: the G data (gamma grid, built at init) is conductor-blind, its
  // floor set by the degree; the conductor only sets how far down the grid is
  // read at compute time. A large enough conductor pushes that read below the
  // floor, and Lfunc_compute returns the recoverable ERR_G_EXTENT.
  uint64_t Lfunc_nmax(Lfunc_t L);
  // If you can't get to nmax, declare how many Euler factors will be provided.
  // NB it takes your word for it and does not check. It calls Lfunc_nmax first
  // (so it triggers and freezes the conductor-derived M), then lowers M to nmax;
  // nmax must be strictly below the current M, else it returns false.
  bool Lfunc_reduce_nmax(Lfunc_t LL, uint64_t nmax);

  // lpoly_callback is called for each prime p <= M (the conductor-derived bound;
  // see Lfunc_nmax). The first call here triggers Lfunc_nmax, computing and
  // freezing M if not already done. Setting poly to zero stops the iteration,
  // lowers M to p-1, and flags ERR_INSUFF_EULER.
  Lerror_t Lfunc_use_all_lpolys(Lfunc_t L, void (*lpoly_callback) (acb_poly_t lpoly, uint64_t p, int d, int64_t prec, void *parm), void *param);

  // you provide one Euler polynomial at a time
  void Lfunc_use_lpoly(Lfunc_t L, uint64_t p, const acb_poly_t poly);

  // Once all polys have been provided, do the computation
  Lerror_t Lfunc_compute(Lfunc_t L);

  // Prepare certified samples and the root number for special-value evaluation,
  // without computing the analytic rank, Taylor coefficient, or zeros. This
  // retains the coefficient-tail, convolution, FFT and sampling error bounds.
  // Call exactly one of Lfunc_compute_values and Lfunc_compute on an object:
  // both destructively normalize its coefficient storage. Lfunc_Taylor and
  // Lfunc_zeros are not available after this values-only entry point.
  Lerror_t Lfunc_compute_values(Lfunc_t L);

  // what working precision did the computation use
  int64_t Lfunc_wprec(Lfunc_t L);

  // return the sign Lambda(s)=sign Lambda(k-s)
  acb_srcptr Lfunc_sign(Lfunc_t L);

  // return the sqrt of sign (chosen to make Lambda(delta)>0
  acb_srcptr Lfunc_sqrt_sign(Lfunc_t L);

  // return the zeros, side = 0,1 for L, conjugate L
  // if rank =0,1 this list is complete, providing
  // the error code did not have ERR_RH_ERROR set
  // otherwise zeros may be missing
  arb_srcptr Lfunc_zeros(Lfunc_t L, uint64_t side);

  // return rank
  // rank=0,1 is rigorous.
  // rank>1 isn't
  int64_t Lfunc_rank(Lfunc_t L);

  // return the first non-zero Taylor coefficient
  // L^(rank)((w + 1)/2)/rank!
  // where w/2 is the normalization (if algebraic, w = motivic weight)
  // same caveats as rank
  arb_srcptr Lfunc_Taylor(Lfunc_t L);

  // return roughly n_points of exp(i theta(t)) L(k/2+it)
  // covering t=[0,max_t]
  // for L or conjugate L
  // returned as doubles in an Lplot_t structure
  Lplot_t *Lfunc_plot_data(Lfunc_t L, uint64_t side, double max_t, uint64_t n_points);

  // reclaim memory from an Lplot_t structure
  void Lfunc_clear_plot(Lplot_t *Lp);

  // compute L(re+i*im)
  // warning, accuracy falls away rapidly as one moves away
  // from the critical line. Should return something sensible
  // for L(k) and L(0)
  // for re+i*im = (w + 1)/2, use Lfunc_Taylor
  Lerror_t Lfunc_special_value_choice(acb_t res, acb_t res_dash, Lfunc_t LL, double re, double im, bool lam_p, bool do_dash);
  Lerror_t Lfunc_special_value(acb_t res, Lfunc_t LL, double re, double im);

  // reclaim memory from an Lfunc_t structure
  void Lfunc_clear(Lfunc_t L);

#ifdef __cplusplus
}
#endif
#endif
