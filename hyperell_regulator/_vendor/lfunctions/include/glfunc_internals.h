#ifndef _GLFUNC_INTERNAL
#define _GLFUNC_INTERNAL

#include "inttypes.h"
#include <flint/acb.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include "glfunc.h"


#define OUTPUT_RATIO (8) // we will analyse 1/this of B
#define TURING_RATIO (16) // and use a further 1/this for Turing
#define EXTRA_BITS (35) // extra bits of precision for convolves etc.
#define verbose (false) // compile-time toggle for all diagnostic printf
#define BAD_64 (1LL<<62) // sentinel: an arb value did not pin to a unique integer

#ifdef BUTHE
// how many integrals have I precomputed for Buthe's method?
#define MAX_MUI (100)
#define MAX_MUI_2 (200) // = 2*MAX_MUI; a Buthe table row holds 2*MAX_MUI_2+1 entries
#define MAX_MU ((double) MAX_MUI) // largest mu Buthe can handle (MAX_MUI as a double)
#define MAX_MU_2 ((double) MAX_MUI_2) // MAX_MUI_2 as a double
#endif

#define MAX_L (10) // maximum differential allowed in upsampling


#define COMPUTE_ZEROS // compile-time switch: enable the zero-finding phase
#define COMPUTE_RANK // compile-time switch: enable the rank-determination phase


#ifdef __cplusplus
extern "C"{
#endif

  typedef struct{
    uint64_t degree; // degree r (number of Gamma_R factors), 2 <= r <= MAX_DEGREE
    uint64_t conductor; // arithmetic conductor N
    double normalisation; // algebraic->analytic s-shift = (k-1)/2; only used to form mus
    double *mus; // Gamma_R shifts (user mus + normalisation), sorted, non-neg half-ints
    int64_t target_prec; // target output precision in bits (default DEFAULT_TARGET_PREC=100)
    arb_t zero_prec; // intended zero half-width 2^(-target_prec-1); only feeds zero_error
    arb_t zero_error; // intended zero error ball of radius zero_prec; unused
    int64_t wprec; // internal working precision in bits = target_prec + decay + EXTRA_BITS
    int64_t gprec; // precision in bits used by the G/gamma-factor computation (>= wprec)
    char *cache_dir; // directory for cached G-data; used only in default-precision mode
    int self_dual; // self-dual? DK(-1)/NO(0)/YES(1); YES skips dual-side zero-finding (Turing then doubles the L-side)
    int rank; // analytic rank = order of vanishing at s=1/2; DK(-1) until computed
    arb_t mu; // scalar -1/2+(1+sum mus)/r (Artin Lemma 5.2); unused
    arb_t nu; // scalar 1/2+(2/r)*sum(nus); shift param for the error bounds (Lemma 2/5)
    arb_t *nus; // shifts nus[j]=(mus[j]-1/2)/2 for j<2 else (mus[j]-1)/2 (Lemma 2/5)
    arb_t C;           // coeff bound: |a_n| <= C*n^alpha
    arb_t alpha;       // alpha = 1 under Ramanujan (set in g.c)
    //arb_t k0; // no longer used
    //arb_t k2;
    double one_over_B; // 1/B = degree/512; the G-kernel u-grid step is 2*pi/B
    arb_t B; // scaling parameter B = 1/one_over_B; sets the t-spacing of Lambda samples
    arb_t two_pi_by_B; // 2*pi/B: u-spacing between consecutive G/Lambda sample indices
    arb_t pi; // the constant pi at working precision wprec
    int64_t low_i;     // bottom of G u-grid: floor(u_min/(2*pi/B)), u_min=-32*ln2; reaches the m=1 term
    int64_t hi_i;      // top of G u-grid: smallest i with |G(u_i)|<=2^-prec; sets M=sqrt(N)*exp(2*pi*(hi_i+.5)/B)
    uint64_t max_K;    // number of Taylor terms per grid point (first index of Gs)
    arb_t eq59;        // eq-(59) tail-truncation bound (~ 2^-prec)

    arb_t **Gs;        // Gs[k][i-low_i] = k-th Taylor coeff G^(k)(u_i)/k! at u_i=i*(2*pi/B); built in g.c

    // computation related
    uint64_t fft_N; // length of the short DFT for the Euler-product convolutions (2^11)
    uint64_t fft_NN; // length of the final output iFFT (2^16)
    double A; // output sampling rate = fft_NN/B; output sample i is Lambda at t = i/A
    arb_t arb_A; // A as a rigorous real ball (for the error analysis)
    arf_t arf_A; // A as an arf_t; only used to form arf_one_over_A
    arb_t one_over_A; // 1/A as a ball; the t-grid spacing (sample i has imag part i/A)
    arf_t arf_one_over_A; // 1/A as an arf_t; unused
    acb_t *G; // length-fft_N acb scratch staging one Gs[k] column per convolution
    acb_t *w; // twiddle factors for length fft_N
    acb_t *ww; // ditto for fft_NN
    arb_t *zeros[2]; // zero ordinates t on the critical line; [0]=L, [1]=dual L
    double eta; // contour-tilt knob in delta=(1-eta)*pi/2; =0, only feeds delta
    arb_t delta; // contour offset (1-eta)*pi/2 = pi/2; only feeds (unused) exp_delta
    arb_t exp_delta; // exp(-delta) = exp(-pi/2); unused
    acb_t *kres; // per-k convolution result, accumulated into res
    acb_t *res; // length-fft_NN buffer: convolution sum, then iFFT'd to Lambda samples
    acb_t **skm; // skm[k][n]: a_m/sqrt(m)*(log(m/sqrtN)-u_m)^k bucketed by freq index n
    arb_t pre_ftwiddle_error; // conductor-free bound on the truncated Dirichlet-series tail
    arb_t ftwiddle_error; // truncation bound = pre_ftwiddle_error*sqrt(N); added per sample

#ifdef BUTHE
    arb_t buthe_Wf; // prime/Euler-sum term of Buthe's inequality
    arb_t buthe_Winf; // archimedean (gamma) term of Buthe's inequality (uses buthe_ints)
    arb_t buthe_Ws; // sum-over-computed-zeros term of Buthe's inequality
    arb_t buthe_b; // height up to which RH is confirmed, b = B/OUTPUT_RATIO
    arb_t buthe_C; // Buthe constant = degree (Lemma 3.4); unused
    arb_t buthe_h; // Buthe test-function step h = BUTHE_H
    arb_t buthe_ints[(MAX_R-1)*(2*MAX_MUI_2+1)]; // precomputed gamma-integral table [deg-2][2*mu]
    uint64_t buthe_M; // prime-power cutoff for the Wf sum = floor(sqrt(#coeffs))
#endif

#ifdef TURING
    arb_t X; // Turing zero-counting verification height (see eq 4.10)
    arb_t imint; // integral of Im sum_j logGamma((1/2+it+mu_j)/2) over the Turing region
#endif
    
    arb_t one_over_root_N; // 1/sqrt(N), the N^{s/2} factor base in Lambda
    arb_t sum_ans; // running sum of |a_n|/sqrt(n), n<=M; scales the eq59 tail error
    acb_t sqrt_sign; // unit-modulus square root of the root number (sqrt_sign^2 = sign)
    acb_t sign; // root number eps in Lambda(s)=eps*Lambda(1-s); = sqrt_sign^2
    acb_t *ans; // Dirichlet coefficients: ans[n-1] = a_n (analytic norm), len allocated_M
    uint64_t M; // largest n with a_n used = floor(sqrt(N)*exp(2*pi*(hi_i+0.5)/B))
    uint64_t M0; // split point ceil(sqrt(N)/100): n<M0 summed directly, n>=M0 via FFT
    uint64_t allocated_M; // allocated capacity of ans (slots); M0, M <= allocated_M
    double dc; // sqrt(conductor); log-scale center mapping coeff index to FFT bin

    bool nmax_called; // true if user/system has called Lfunc_nmax

    int64_t offset; // FFT bin index of a_1 = calc_m(1); used only for a bounds check

    uint64_t u_N; // upsampling kernel half-width = ceil(A^2 h^2/2); kernel = 2*u_N+1 pts
    //arb_t *u_coshs;
    //arb_t *u_exps;
    arb_t u_H; // Gaussian window width h in the interpolation kernel
    arb_t u_pi_by_H2; // -pi/h^2: coefficient of (t-k/A)^2 in the Gaussian window
    arb_t u_A; // upsampling sample rate = A*u_stride (= A, since u_stride = 1)
    arb_t u_one_over_A; // 1/u_A: the upsampling grid spacing in t
    arb_t *u_values[2]; // per-side raw samples of Lambda, sign-normalized so Lambda(0)>=0 (full buffer incl. left guard)
    arb_t *u_values_off[2]; // view into u_values past the guard: index k = Lam(k/A)
    // Lam(t)=eps^1/2 N^it/2 prod gamma_r(1/2+it+mu_j) L(1/2+it)
    uint64_t u_no_values; // element count of each u_values[side]
    uint64_t u_no_values_off; // count of each u_values_off[side] = u_no_values - 2*u_N*u_stride
    uint64_t u_stride; // sample decimation; hardwired to 1 (use every FFT grid point)
    arb_t u_pi_A; // pi*u_A: sinc-argument scale and (pi*A)^d derivative factor
    arb_t upsampling_error; // rigorous upsampling error bound, added to every upsampled value

    arb_t Lam_d; // Lambda^(rank)(1/2)
    arb_t L_d; // L^(rank)(1/2)/rank!
  } Lfunc;

  // from glfunc_g.c
  Lerror_t compute_g(Lfunc *);

  // from acb_fft.c
  void acb_initfft(acb_t *w, uint64_t n, uint64_t prec);
  void acb_fft(acb_t *x, uint64_t n, acb_t *w, uint64_t prec);
  void acb_ifft(acb_t *x, uint64_t n, acb_t *w, uint64_t prec);
  void acb_convolve(acb_t *res, acb_t *x, acb_t *y, uint64_t n, acb_t *w, uint64_t prec);
  void acb_convolve1(acb_t *res, acb_t *x, acb_t *y, uint64_t n, acb_t *w, uint64_t prec);
  void acb_convolve2(acb_t *res, acb_t *x, acb_t *y, uint64_t n, acb_t *w, uint64_t prec);

  // from error.c
  void abs_gamma(arb_t res, acb_t s, Lfunc *L, int64_t prec);
  void init_ftwiddle_error(Lfunc *L, int64_t prec);
  void complete_ftwiddle_error(Lfunc *L, int64_t prec);
  Lerror_t do_pre_iFFT_errors(Lfunc *L);
  bool M_error(arb_t res, arb_t x, Lfunc *L, int64_t prec);
  bool M_error_at_cutoff(arb_t res, arb_t x, Lfunc *L, uint64_t cutoff, int64_t prec);
  int64_t calc_m(uint64_t i, double grid_step, double dc);

#ifdef BUTHE
  // from buthe.c
  void init_buthe(Lfunc *L, int64_t prec);
  void wf(Lfunc *L, uint64_t p, acb_poly_t fp1, acb_poly_t fp, int64_t prec);
  void buthe_Wf_error(Lfunc *L);
  Lerror_t buthe_check_RH(Lfunc *L);
#endif

#ifdef TURING
  // from turing.c
  Lerror_t turing_check_RH(Lfunc *L, int64_t);
#endif
  
  // from compute.c
  void lfunc_compute(Lfunc *L);

  //from upsample.c
  double upsample_error(long double M, long double H, long double h, long double A, double *mus, uint64_t r, uint64_t N, long double T, long double imz, uint64_t l);
  Lerror_t init_upsampling(Lfunc *L);
  bool newton(arb_ptr res, arb_ptr t0, Lfunc *L, uint64_t side, uint64_t prec);
  bool upsample_stride(arb_ptr res, arb_ptr t0, Lfunc *L, uint64_t side, uint64_t prec);
  Lerror_t arb_upsampling_error(arb_t res, double M,double H,double h,double A,double *mus,uint64_t r,uint64_t N,double T,arb_t imz, uint64_t l, arb_t pi, int64_t prec);

  //from zeros.c
  Lerror_t find_zeros(Lfunc *L, uint64_t side);

  // from rank.c
  //uint64_t guess_rank(Lfunc *L, uint64_t side, uint64_t prec);
  Lerror_t do_rank(Lfunc *L);

#ifdef __cplusplus
}
#endif
#endif
